# simple nimble model with just a spline on depth (for all cetaceans)
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
library(mgcv)

# Import the jagam object we created previously
load("ProcessedData/jagam_m1.0.RData")

m1.0 <- gam(Detected ~ s(depth, k = 5, bs = "bs"),  
            family = "binomial", data = detect_data, method="REML", 
            select = TRUE) # to create two lambdas


# Import the data
load("./ProcessedData/detect_data.RData")
mm.data <- detect_data


# because I've removed some of the unique biosample reference numbers with the above 
# methods removal, I need to renumber the unique biosamples so that they are consecutive
# and start with 1. I do this using the factor trick. 
mm.data$unique_biorep_numeric <- as.numeric(as.factor(mm.data$NWFSCsampleID))

# create a site variable that can hold multiple depths
mm.data$site.numeric <- as.numeric(factor(paste0(mm.data$utm.lat, mm.data$utm.lon)))

# create a primer variable
mm.data$primer.numeric <- as.numeric(factor(mm.data$primer))

# pull unique info for each bio sample (can add to these for enviro covariates)
biosamp_dat <- mm.data %>%
  group_by(unique_biorep_numeric) %>%
  slice(1) %>%
  ungroup()

# make data, index vectors and constants for nimble work
biosamp_station_index <- biosamp_dat$site.numeric
#biosamp_method_index <- biosamp_dat$Collection_method_numeric
Y_biosamp_index <- mm.data$unique_biorep_numeric
biosamp_Volume_filt_mL <- as.numeric(biosamp_dat$volume) - mean(as.numeric(biosamp_dat$volume)) #centered water volumes
biosamp_depth_Depth_m <- as.numeric(biosamp_dat$depth) - mean(as.numeric(biosamp_dat$depth)) #centered sample depths
Y_primer_index <- mm.data$primer.numeric

N <- dim(mm.data)[1]
n_sites <- length(unique(mm.data$site.numeric))
n_biosamples <- length(unique(mm.data$unique_biorep_numeric))
#n_methods <- length(unique(mm.data$Collection_method_numeric ))
n_primers <- length(unique(mm.data$primer.numeric))

Y <- mm.data$Detected 

#############
# SITE x DEPTH MODEL: Site-depth specific occupancy states
#############

# Create site-depth combinations for occupancy states
biosamp_dat$site_depth_combo <- paste(biosamp_dat$site.numeric, 
                                      round(biosamp_dat$depth, 1), 
                                      sep = "_")
unique_site_depth <- unique(biosamp_dat$site_depth_combo)
n_site_depth_states <- length(unique_site_depth)

# Create mapping from site-depth combo to numeric index
site_depth_lookup <- data.frame(
  combo = unique_site_depth,
  index = 1:n_site_depth_states
)

# Extract site and depth for each site-depth state
site_depth_lookup$site <- as.numeric(sapply(strsplit(site_depth_lookup$combo, 
                                                     "_"), `[`, 1))
site_depth_lookup$depth <- as.numeric(sapply(strsplit(site_depth_lookup$combo, 
                                                      "_"), `[`, 2))

# Create index vector for biosamples to site-depth states
biosamp_dat$site_depth_index <- match(biosamp_dat$site_depth_combo, 
                                      site_depth_lookup$combo)
biosamp_site_depth_index <- biosamp_dat$site_depth_index

# Center depths for site-depth states
site_depth_depths <- site_depth_lookup$depth - mean(site_depth_lookup$depth)

### MODEL

m1.0_nimble <- nimbleCode({
  # Priors
  ## Parametric effect priors CHECK tau=1/11^2 is appropriate!
  for (i in 1:1) { b_depth[i] ~ dnorm(0,0.0087) }
  ## prior for s(depth)... 
  K1[1:4,1:4] <- S1[1:4,1:4] * lambda[1]  + S1[1:4,5:8] * lambda[2]
  # S1 is a penalty matrix (the second bit 5:8 is S2)
  # K1 is the sum of the penalties each scaled by the smoothing parameters
  # effectively the same as select = TRUE, i.e., two penalties
  # bc if you didn't, it would be an improper prior
  b_depth[2:5] ~ dmnorm(zero[2:5], K1[1:4,1:4]) 
  ## smoothing parameter priors CHECK...
  for (i in 1:2) {
    lambda[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambda[i])
  }
  # Prior for effect of volume
  b_vol ~ dnorm(0, 1)
  # Capture probability parameters
  cap_prob_hat ~ dnorm(0, 1.7)
  cap_prob_SD ~ dexp(1)
  
  # primer detection probability priors
  prob_detection[1] ~ dbeta(1, 1)
  prob_detection[2] ~ dbeta(1, 1)
  prob_detection[3] ~ dbeta(1, 1)
  prob_detection[4] ~ dbeta(1, 1)
  
  
  # Linear predictor, effect of depth
  eta[1:N] <- X[1:N, 1:5] %*% b_depth[1:5] 
    
  # Depth-level occurrence
  for (i in 1:n_site_depth_states) { 
    prob_site_depth_occurrence[i] <- ilogit(eta[i]) # probability of occurrence at each depth
    site_depth_occurrence[i] ~ dbern(prob_site_depth_occurrence[i]) # OCCUPANCY (UNOBSERVED)
    } # end for n_site_depth_states

  # Site-level capture probability variation
  for (i in 1:n_sites) {
    cap_prob_logit_site[i] ~ dnorm(cap_prob_hat, cap_prob_SD)
    } # end for n_sites
  
  # Biosample-level capture 
  for (i in 1:n_biosamples) {
    logit(prob_capture[i]) <- ilogit(cap_prob_logit_site[biosamp_station_index[i]] +
    b_vol * biosamp_Volume_filt_mL[i]) 
    bio_capture[i] ~ dbern(site_depth_occurrence[biosamp_site_depth_index[i]] *
                               prob_capture[i])
    } # end for n_biosamples
    
  # Replicate-level likelihood of detecting in each of N lab samples
  for (i in 1:N) {
    Y[i] ~ dbern(bio_capture[Y_biosamp_index[i]] * prob_detection[Y_primer_index[i]]) 
    } # end for N
   
}) # end model definition

# Data
data <- list(X = q1Model_m1.0$jags.data$X,
             S1 = q1Model_m1.0$jags.data$S1,
             Y = Y,
             biosamp_Volume_filt_mL = biosamp_Volume_filt_mL)

# Constants
constants <- list(  n_sites = n_sites,
                    n_biosamples = n_biosamples,
                    n_site_depth_states = n_site_depth_states,
                    N = N,
                    biosamp_station_index = biosamp_station_index,
                    biosamp_site_depth_index = biosamp_site_depth_index,
                    Y_biosamp_index = Y_biosamp_index,
                    Y_primer_index = Y_primer_index,
                    zero = rep(0, 10))


# UPDATED: Using a function for initial values for the second model as well.
inits_fn_sitedepth <- function(){
  # Determine if a biosample had any positive detections
  max_obs_per_biosample <- tapply(Y, Y_biosamp_index, max)
  bio_capture_init <- rep(0, n_biosamples)
  bio_capture_init[as.numeric(names(max_obs_per_biosample))] <- max_obs_per_biosample
  
  # Determine if a site-depth combo had any positive biosamples
  site_depth_occurrence_init_raw <- tapply(bio_capture_init, 
                                           biosamp_site_depth_index, 
                                           max)
  site_depth_occurrence_init <- rep(0, n_site_depth_states)
  site_depth_occurrence_init[as.numeric(names(site_depth_occurrence_init_raw))] <- site_depth_occurrence_init_raw
  
  list(
    # Latent states (from before)
    bio_capture = bio_capture_init,
    site_depth_occurrence = site_depth_occurrence_init,
    
    # Model parameters 
    b_vol = 0,
    cap_prob_hat = 0,
    cap_prob_SD = 1,
    cap_prob_logit_site = rep(0, n_sites),
  #  prob_capture = rep(0.5, n_biosamples),
    prob_detection = rep(0.5, n_primers),
    lambda = q1Model_m1.0$jags.ini$lambda, # 2 vec
    rho = log(q1Model_m1.0$jags.ini$lambda),
    b_depth = m1.0$coefficients, # 5 vec
    eta = rep(0, 25088),
  K1 = matrix(rep(0, 4*4), nrow = 4)
  )
}

# Run NIMBLE model
nimbleOut_m1.0_Occ <- nimbleMCMC(code = m1.0_nimble, 
                             data = data, 
                             inits = inits_fn_sitedepth(),
                             constants = constants,
                             niter = 250000, 
                             nburnin = 225000, 
                             thin = 10, 
                             nchains = 4,
                             summary=TRUE,
                             samplesAsCodaMCMC = TRUE,
                             WAIC = TRUE)

save(nimbleOut_m1.0_Occ, file = "./Results/nimbleOut_m1.0_Occ.RData")

# Gelman-Rubin diagnostic
MCMCsummary(nimbleOut_m1.0_Occ$samples)

# Visualize MCMC chains
mcmcplot(nimbleOut_m1.0_Occ$samples)

n.post <- 10000
post.samples <- rbind.data.frame(nimbleOut_m1.0$samples$chain1,
                                 nimbleOut_m1.0$samples$chain2,
                                 nimbleOut_m1.0$samples$chain3,
                                 nimbleOut_m1.0$samples$chain4)

# reconstruct the model expectation

X = q1Model_m1.0$jags.data$X

mu.post <- matrix(rep(0, nrow(X)*nrow(post.samples)), nrow = nrow(X))

# create a new lp matrix

#predict.gam(object = m1.0, type = "lpmatrix")


for (i in 1:nrow(post.samples)){
  eta.post <- X[1:nrow(X), 1:5] %*% as.numeric(post.samples[i,1:5])
  mu.post[1:nrow(X), i] <- as.numeric(ilogit(eta.post))}

# note mu.post is [1:27540, 1:4000] 

mu.post.long <- as.data.frame(cbind(Depth = detect_data$depth, mu.post)) %>%
  pivot_longer(cols = 2:(nrow(post.samples)+1), names_to = "Chain", values_to = "PDetect")

mu.post.med <- mu.post.long %>%
  group_by(Depth) %>%
  summarize(Med = median(PDetect),
            LCI = quantile(PDetect, 0.025),
            UCI = quantile(PDetect, 0.975))

p <- ggplot() +
  geom_ribbon(data = mu.post.med, 
            aes(x=Depth, ymin= LCI, ymax = UCI), fill = "lightgrey") +
  geom_line(data = mu.post.med, aes(x=Depth, y = Med))+
  theme_bw()

ggsave(plot = p, file = "./Figures/m1.0_nimbleOccModel.png", width = 4, height = 4, units = "in")

names(m1.0_sePreds)[1] <- "Depth"
m1.0_compare <- left_join(m1.0_sePreds, mu.post.med, by = "Depth")

ggplot(m1.0_compare) +
  geom_line(aes(x=Depth, y = mu))+
  geom_line(aes(x=Depth, y = mu_jags), color = "blue")+
  geom_point(aes(x=Depth, y = Med), color = "green")+
  ylab("P(Detection)")+
  xlab("Depth")+
  theme_bw()
