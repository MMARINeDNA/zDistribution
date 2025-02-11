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

lags.data <- filter(detect_species_meta, BestTaxon == "Ziphius cavirostris")
#lags.data <- detect_species_meta
# each biological replicate needs a unique ID, IDK why it needs to be numeric/ordered
lags.data$Bio_UID <- as.numeric(factor(lags.data$NWFSCsampleID))

# create a site variable that can hold multiple depths
lags.data$Site <- as.numeric(factor(paste0(lags.data$utm.lat, lags.data$utm.lon)))

biosamp_data <- lags.data %>%
  group_by(Bio_UID) %>%
  slice(1) %>%
  ungroup()

# note this is each loc + depth as a unique biosamp station
biosamp_station_index <- biosamp_data$Site
Y_biosamp_index <- lags.data$Bio_UID
biosamp_Volume_filt_mL <- as.numeric(biosamp_data$volume) - mean(as.numeric(biosamp_data$volume))
biosamp_depth_Depth_m <- as.numeric(biosamp_data$depth) - mean(as.numeric(biosamp_data$depth))
Y_station_index <- lags.data$Site
Y_primer_index <- as.numeric(factor(lags.data$primer))

N <- nrow(lags.data)
n_sites <- length(unique(lags.data$Site))
n_biosamples <- length(unique(lags.data$Bio_UID))
n_primers <- length(unique(lags.data$primer))

Y <- lags.data$Detected


edna_code_vol_depth_meth_randCap <- nimbleCode({
  cap_prob_hat ~ dnorm(0,1.7) 
  cap_prob_SD~ dexp(1)
  
  for (i in 1:n_sites) {
    # Site-level occurrence probability
    site_occurrence[i] ~ dbern(prob_occurrence)
    cap_prob_logit[i] ~ dnorm(cap_prob_hat, cap_prob_SD)
  }
  
  for (i in 1:n_biosamples) {
    
    logit(prob_capture[i]) <- cap_prob_logit[biosamp_station_index[i]] +   # categorical differences in capture method # removed [biosamp_method_index[i]]
      b_depth * biosamp_depth_Depth_m[i] +    # continuous effect of depth
      b_vol * biosamp_Volume_filt_mL[i]      # continuous effect of volume
    
    # biosample-level occurrence probability
    bio_capture[i] ~ dbern(site_occurrence[biosamp_station_index[i]] *
                             prob_capture[i])
    # `biosamp_station_index` -- index vector that specifies which site each 
    # biosample belongs to
    # `biosamp_method_index` -- index vector that specifies which method each 
    # biosample belongs to
  }
  
  # Likelihood of detecting in each lab sample (tech replicate)
  for (i in 1:N) {
    Y[i] ~ dbern(bio_capture[Y_biosamp_index[i]] * prob_detection[Y_primer_index[i]]) 
    # `Y_biosamp_index[i]` -- index vector that specifies which biological sample each 
    # observation belongs to
  }
  
  # Priors for the parameters
  prob_occurrence <- 1
  b_depth ~ dnorm(0, 1)
  b_vol ~ dnorm(0, 1)
  prob_detection[1] ~ dbeta(1, 1) #detection prob for dloop, MiFish, MarVer
  prob_detection[2] ~ dbeta(1, 1) 
  prob_detection[3] ~ dbeta(1, 1)
})

# Define the constants, data, and initial values
constants <- list(
  n_sites = n_sites,
  n_biosamples = n_biosamples,
  n_primers = n_primers,
  N = N,
  biosamp_station_index = biosamp_station_index,
  Y_biosamp_index = Y_biosamp_index,
  Y_primer_index = Y_primer_index
)


# Define the data
data <- list(
  Y = Y, # detections by sample
  biosamp_Volume_filt_mL = biosamp_Volume_filt_mL, #centered water volumes
  biosamp_depth_Depth_m = biosamp_depth_Depth_m #centered sample depths
)

inits <- list(
  prob_occurrence = 0.5,  # Initial value for probability of site occurrence
  site_occurrence = rep(1, n_sites),  # Initial values for site occurrence, all set to 1 (can be set randomly between 0 and 1)
  #bio_capture = rep(1, n_biosamples),  # Initial values for biosample capture
  
  # Initial values for method-specific parameters
  #prob_detection = rep(0.5, n_methods),  # Initial values for detection probability
  
  # Initial values for coefficients
  b_depth = 0,  # Initial value for depth coefficient
  b_vol = 0     # Initial value for volume coefficient
)


# Run NIMBLE model
edna_code_vol_depth_meth_randCap.run <- nimbleMCMC(code = edna_code_vol_depth_meth_randCap, 
                                                   constants = constants, 
                                                   data = data, 
                                                   inits = inits,
                                                   niter = 200000, 
                                                   nburnin = 10000, 
                                                   thin = 100, 
                                                   nchains = 3,
                                                   summary=TRUE,
                                                   samplesAsCodaMCMC = TRUE,
                                                   WAIC = TRUE)

# Gelman-Rubin diagnostic
MCMCsummary(edna_code_vol_depth_meth_randCap.run$samples)

# Visualize MCMC chains
mcmcplot(edna_code_vol_depth_meth_randCap.run$samples)

n.post <- 1900
post.samples <- rbind.data.frame(edna_code_vol_depth_meth_randCap.run$samples$chain1,
                                 edna_code_vol_depth_meth_randCap.run$samples$chain2,
                                 edna_code_vol_depth_meth_randCap.run$samples$chain3)
b_depth <- post.samples$b_depth
b_vol <- post.samples$b_vol
prob_detection_1 <- post.samples$`prob_detection[1]`
prob_detection_2 <- post.samples$`prob_detection[2]`
prob_detection_3 <- post.samples$`prob_detection[3]`

# --- Define sequences for depth and volume ---
depth_min <- min(as.numeric(biosamp_data$depth))
depth_max <- max(as.numeric(biosamp_data$depth))
n_depth_steps <- 100
depth_seq <- seq(from = depth_min, to = depth_max, length.out = n_depth_steps)

vol_min <- min(as.numeric(biosamp_data$volume))
vol_max <- max(as.numeric(biosamp_data$volume))
n_vol_steps <- 100
vol_seq <- seq(from = vol_min, to = vol_max, length.out = n_vol_steps)

# --- Calculate predicted probabilities ---
plot.stor2 <- matrix(data = NA, n.post, length(depth_seq))
plot.stor <- matrix(data = NA, n.post, length(vol_seq))

for (i in 1:length(depth_seq)) {
  for (j in 1:n.post) {
    plot.stor2[j, i] <- inv.logit(b_depth[j] * (depth_seq[i]-mean(as.numeric(biosamp_data$depth))))
  }
}

for (i in 1:length(vol_seq)) {
  for (j in 1:n.post) {
    plot.stor[j, i] <- inv.logit(b_vol[j] * (vol_seq[i]-mean(as.numeric(biosamp_data$volume))))
  }
}

# --- Create plot data frames ---
plot_data_lineribbon_depth <- data.frame(
  x_value = rep(depth_seq, each = n.post),
  value = as.vector(plot.stor2)
)

plot_data_lineribbon_vol <- data.frame(
  x_value = rep(vol_seq, each = n.post),
  value = as.vector(plot.stor)
)

# --- Create plots ---
p2 <- ggplot(plot_data_lineribbon_depth, aes(x = x_value, y = value)) +
  stat_lineribbon(aes(y = value), alpha = 0.25, fill = "#808080", color = "#000000", .width = c(0.25, 0.5, 0.75)) +
  labs(x = "Depth Sampled", y = "Probability of Capture", title = "Probability of Capture Covariates") +
  theme_minimal()

p1 <- ggplot(plot_data_lineribbon_vol, aes(x = x_value, y = value)) +
  stat_lineribbon(aes(y = value), alpha = 0.25, fill = "#808080", color = "#000000", .width = c(0.25, 0.5, 0.75)) +
  labs(x = "Volume Sampled", y = "Probability of Capture") +
  theme_minimal()

# --- Combine plots ---
combined_data <- rbind(
  data.frame(variable = "Depth (m)", plot_data_lineribbon_depth),
  data.frame(variable = "Volume Filtered (mL)", plot_data_lineribbon_vol)
)

p <- ggplot(combined_data, aes(x = x_value, y = value)) +
  stat_lineribbon(aes(y = value), alpha = 0.25, fill = "#808080", color = "#000000", .width = c(0.25, 0.5, 0.75)) +
  facet_wrap(~ variable, scales = "free_x", nrow = 2) +
  labs(x = "", y = "Probability of Capture", title = "Prob. of Capture Covariates") +
  theme_minimal()

print(p)

# --- Plot 3: Probability of Detection (given occurrence) ---
# Prepare data frame
df_prob_detection <- data.frame(
  Probability = c(prob_detection_1, prob_detection_2, prob_detection_3),
  Method = rep(c("dloop", "Marver", "MiFish"), each = length(prob_detection_1))
)

# Calculate medians
medians_detection <- df_prob_detection %>%
  group_by(Method) %>%
  summarize(Median = median(Probability))

viridis_cols <- viridis(3, begin = 0.3, end = 0.7)
names(viridis_cols) <- c("Dloop", "MarVer", "MiFish")

p3 <- ggplot(df_prob_detection, aes(x = Probability, fill = Method, color = Method)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.3, position = "identity", bins = 30) +
  geom_density(size = 1.2) +
  geom_vline(data = medians_detection, aes(xintercept = Median, color = Method), linetype = "dashed", size = 1) +
  annotate("text", x = Inf, y = Inf, label = names(viridis_cols), color = viridis_cols, hjust = 1.1, vjust = c(3, 4.5, 6), size = 5) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.3, begin = 0.3, end = 0.7) +
  scale_color_viridis(discrete = TRUE, begin = 0.3, end = 0.7) +
  labs(x = "Detection Probability", y = "Density", title = "Primers") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")
p3 <- p3 + scale_x_continuous(limits = c(0, 0.25)) # Set x-axis limits

# Combine the plots using patchwork - Probability of Detection at the bottom
combined_plot <- p1 / p3 / p2

print(combined_plot)

