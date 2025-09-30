# simple nimble model with just a spline on depth
# using the jagam object 

library(MCMCvis)
library(boot)
library(tidyverse)
library(mcmcplots)
library(ggplot2)
library(ggdist)
library(dplyr)
library(tidyr)
library(viridis)
library(patchwork)
library(nimble)

load("./ProcessedData/detect_data.RData")
load("ProcessedData/jagam_m1.0.RData")

# define the model
m1.0_nimble <- nimbleCode({
    eta[1:27540] <- X[1:27540, 1:10] %*% b[1:10] ## linear predictor (b is beta)
    for (i in 1:n) { mu[i] <-  ilogit(eta[i]) } ## expected response (not needed for full model)
    for (i in 1:n) { y[i] ~ dbin(mu[i],w[i]) } ## response (not needed for full model)
    ## Parametric effect priors CHECK tau=1/11^2 is appropriate!
    for (i in 1:1) { b[i] ~ dnorm(0,0.0087) }
    ## prior for s(depth)... 
    K1[1:9,1:9] <- S1[1:9,1:9] * lambda[1]  + S1[1:9,10:18] * lambda[2]
    b[2:10] ~ dmnorm(zero[2:10],K1[1:9,1:9]) 
    ## smoothing parameter priors CHECK...
    for (i in 1:2) {
      lambda[i] ~ dgamma(.05,.005)
      rho[i] <- log(lambda[i])
  }
})

# Data
data <- list(y = q1Model_m1.0$jags.data$y,
             X = q1Model_m1.0$jags.data$X,
             S1 = q1Model_m1.0$jags.data$S1)

# Constants
constants <- list(n = q1Model_m1.0$jags.data$n, # number of data points
                  w = q1Model_m1.0$jags.data$w, # not sure what this is
                  zero = q1Model_m1.0$jags.data$zero)

# Initial values
inits <- list(lambda = q1Model_m1.0$jags.ini$lambda, # 2 vec
              rho = log(q1Model_m1.0$jags.ini$lambda),
              b = q1Model_m1.0$jags.ini$b, # 10 vec
              eta = rep(0, 27540),
              mu = rep(0, 27540),
              K1 = matrix(rep(0, 9*9), nrow = 9))

# Run NIMBLE model
nimbleOut_m1.0 <- nimbleMCMC(code = m1.0_nimble, 
                             data = data, 
                             inits = inits,
                             constants = constants,
                             niter = 50000, 
                             nburnin = 10000, 
                             thin = 100, 
                             nchains = 4,
                             summary=TRUE,
                             samplesAsCodaMCMC = TRUE,
                             WAIC = TRUE)

# Gelman-Rubin diagnostic
MCMCsummary(nimbleOut_m1.0$samples)

# Visualize MCMC chains
mcmcplot(nimbleOut_m1.0$samples)

n.post <- 1600
post.samples <- rbind.data.frame(nimbleOut_m1.0$samples$chain1,
                                 nimbleOut_m1.0$samples$chain2,
                                 nimbleOut_m1.0$samples$chain3,
                                 nimbleOut_m1.0$samples$chain4)

# reconstruct the model expectation

mu.post <- matrix(rep(0, 27540*nrow(post.samples)), nrow = 27540)

# create a new lp matrix

predict.gam(object = m1.0, type = "lpmatrix")

X = q1Model_m1.0$jags.data$X

for (i in 1:nrow(post.samples)){
  eta.post <- X[1:27540, 1:10] %*% as.numeric(post.samples[i,1:10])
  mu.post[1:27540, i] <- as.numeric(ilogit(eta.post))}

mu.post.long <- as.data.frame(cbind(Depth = detect_species_meta$depth, mu.post)) %>%
  pivot_longer(cols = 2:(nrow(post.samples)+1), names_to = "Chain", values_to = "PDetect")

mu.post.med <- mu.post.long %>%
  group_by(Depth) %>%
  summarize(Med = median(PDetect),
            LCI = quantile(PDetect, 0.025),
            UCI = quantile(PDetect, 0.975))

p <- ggplot() +
  geom_ribbon(data = mu.post.med, 
            aes(x=Depth, ymin= LCI, ymax = UCI), fill = "lightgrey") +
  geom_line(data = mu.post.med, aes(x=Depth, y = Med)) +
  theme_bw()

ggsave(plot = p, file = "./Figures/m1.0_nimble.png", width = 4, height = 4, units = "in")

# this doesn't work -- not sure what m1.0_sePreds is/was?

names(m1.0_sePreds)[1] <- "Depth"
m1.0_compare <- left_join(m1.0_sePreds, mu.post.med, by = "Depth")

ggplot(m1.0_compare) +
  geom_line(aes(x=Depth, y = mu))+
  geom_line(aes(x=Depth, y = mu_jags), color = "blue")+
  geom_point(aes(x=Depth, y = Med), color = "green")+
  ylab("P(Detection)")+
  xlab("Depth")+
  theme_bw()
