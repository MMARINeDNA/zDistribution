#### 3D Distribution
#### Export best Q1 model to BUGS
#### February2025
#### EKJ&AVC

library(mgcv)
library(tidyverse)
library(rjags)

load("./ProcessedData/detect_species_meta.RData")
load("./ProcessedData/H1models.Rdata")

### Aggregate model AIC --------------------------------------------------------

modelAIC <- AIC(m1.0, m1.1, m1.2, m1.2a, m1.2b, m1.2c, m1.2d, m1.2e, m1.2f, m1.2g)
#m1.2f has lowest AIC

### Build jags model -----------------------------------------------------------

# simplest model to start with
q1Model_m1.0 <- jagam(Detected ~ s(depth), 
                      family = "binomial", data = detect_species_meta,
                      file = "./ProcessedData/m1.0.jag")

# depth by species
q1Model_m1.2 <- jagam(Detected ~ s(depth, by = as.factor(BestTaxon)), 
                           family = "binomial", data = detect_species_meta,
                           file = "./ProcessedData/m1.2.jag")

### Run jags model -------------------------------------------------------------

# m1.0
jm <-jags.model("./ProcessedData/m1.0.jag",
                data=q1Model_m1.0$jags.data,
                inits=q1Model_m1.0$jags.ini,
                n.chains=4) # changed to 4 so we can assess convergence

list.samplers(jm)
# simplest model has b and lambda
sam <- jags.samples(jm,c("b","lambda"),n.iter=10000,thin=10) # takes aaages
jam <- sim2jam(sam, q1Model_m1.0$pregam)
plot(jam)  # this is just the spline on depth

pd <- data.frame(depth = 0:500)
fv <- predict(jam,newdata=pd, scale = "response")
plot(pd$depth, exp(fv), type = "l")

# m1.2 (not run yet as it seriously takes ages)
jm <-jags.model("./ProcessedData/m1.2.jag",
                data=q1Model_m1.2$jags.data,
                inits=q1Model_m1.2$jags.ini,
                n.chains=4, n.adapt = 2500) # changed to 4 so we can assess convergence

list.samplers(jm)
# simplest model has b and lambda
sam <- jags.samples(jm,c("b","lambda"), n.iter=25000, thin=10) # takes aaages
jam <- sim2jam(sam, q1Model_m1.0$pregam)
plot(jam)  # this is just the spline on depth
