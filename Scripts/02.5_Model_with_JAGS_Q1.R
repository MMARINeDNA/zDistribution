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
q1Model <- jagam(Detected ~ s(time_per_m, by = as.factor(BestTaxon)) +
  s(depth, by = as.factor(BestTaxon)), 
family = "binomial", data = detect_species_divetime,
file = "./ProcessedData/Q1Model_bayes.jag")

### Run jags model -------------------------------------------------------------
# set random seed
q1Model$jags.ini$.RNG.name <- "base::Mersenne-Twister" ## setting RNG
q1Model$jags.ini$.RNG.seed <- 6 ## how to set RNG seed

jm <-jags.model("./ProcessedData/Q1Model_bayes.jag",
                data=q1Model$jags.data,
                inits=q1Model$jags.ini,
                n.chains=1)

### TODO:
### AVC: Running through test code, here I get a syntax error at line 2. Trying to
### learn more about JAGS syntax for priors...

list.samplers(jm)
sam <- jags.samples(jm,c("b","rho","scale","mu"),n.iter=10000,thin=10)
jam <- sim2jam(sam,jd$pregam)
plot(jam,pages=1)
jam
pd <- data.frame(x0=c(.5,.6),x1=c(.4,.2),x2=c(.8,.4),x3=c(.1,.1))
fv <- predict(jam,newdata=pd)

