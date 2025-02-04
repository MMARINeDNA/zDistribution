# 3DDistribution
Modelling 3D distributions of marine mammal detections from eDNA samples, with the goal of answering the following questions:

1. Does the probability of a detecting cetaceans in eDNA samples vary with sample depth?
   - H0: Probability of detection does not vary with depth.
   - H1: Probability of detection varies across depth agnostic to species or functional group.
   - H2: Probability of detection varies across depth according to species or functional group.
2. Does the probability of a detecting cetaceans in eDNA samples change with the number of technical replicates?
   * NOTE should we wrap dilution into this question as well?
   - H0: Detection does not vary with # technical replicates.
   - H1: Detection varies with # technical replicated agnostic to species/fuctional group or depth.
   - H2: Detection varies with # technical replicates according to species/functional group, depth, or a combination of the two.
3. Does depth distribution of detections vary across xy spatial distribution?
   - H0: Depth distribution of detections does not vary across xy spatial distribution.
   - H1: Depth distribution of detection does vary across xy spatial distribution agnostic to oceanography (e.g. upwelling).
   - H2: Depth distribution of detection does vary across xy spatial distribution according to oceanography (e.g. upwelling).

# The Plan

Repo will contain code, data must be downloaded from [Google Drive](https://drive.google.com/drive/folders/1EZEfbxgRszwmN4RmaoQe7wh6S6zPgm5A?usp=drive_link) and placed in a folder named "Data".

Marine mammal metadata (taxonomy, species-specific time-at-depth and prey preference): contact AVC. Dive data come from [here](https://apps.dtic.mil/sti/tr/pdf/ADA560975.pdf) and [here](https://www.nepa.navy.mil/Portals/20/Documents/aftteis4/Dive%20Profile%20and%20Group%20Size_TR_2017_05_22.pdf).

EKJ TODO: clean up Analysis 3

# Analysis 1

- collapse all data across X and Y, ignore xy distribution
- model detection probability by depth ignoring species
- compare with models incorporating species, family, and prey category

# Analysis 2

- Using best model(s) from Analysis 1, retest varying number of replicates
- Potential alternative: run [Brice's replication model](https://github.com/BriceSemmens/eDNA_patch) without assuming species' presence and adding depth as a covariate

# Analysis 3

- develop 2D models of distribution for each species
- given that distribution, does sampling depth matter?
- does it vary by species or species type?
- develop explicit 3D models of distribution for each species
- assume that detections at depth reflect species z distribution
- incorporate spread and decay from Mod 1?
- does eDNA depth distribution interact with xy distribution due to differences in, e.g. oceanographic upwelling or downwelling?
- incorproate prey (Note for future consideration: this is one of the broad Mod 3 goals but not sure if it fits in this paper)



# Taking Brice's approach

We don't have observations of the presence or absence of marine mammals at each site, but we can still treat site occupancy as a latent variable (per species). 

We have: multiple sites, each of which contain biological samples at multiple depths (not replicates), then multiple primers x tech reps. Each has an associated volume and dilution. 

* need to add a dilution coefficient, or does this happen at the tech rep level?

### Hierarchical Model Structure

1. **Site-Level Occurrence**
   The presence of a marine mammal species at site $s$ is modeled as a Bernoulli random variable:

   $Z_s \sim \text{Bernoulli}(\psi)$

   where $Z_s$ is the site-level occurrence indicator, and $\psi$ is the overall occurrence probability, drawn from a Beta prior:

   $\psi \sim \text{Beta}(1,1)$

   --> note that in our case, particularly since we are working at the species level, psi might come from some field over X, Y

2. **Capture within a depth x primer x tech rep reaction **
   The logit-linear model for capture probability at depth $d$ is:

   $\text{logit}(p_{\text{capture},d}) = \beta_0 + \beta_{\text{vol}} \cdot X_{\text{vol},d} + \beta_{\text{depth}} \cdot X_{\text{depth},d} + \gamma_{s} + \delta_{\text{primer}[b]}$

   Where:
   - $p_{\text{capture},d}$ is the capture probability for depth $d$
   - $\beta_0$ is the intercept (site-level capture probability hyperparameter)
   - $\beta_{\text{vol}}$ is the volume coefficient
   - $\beta_{\text{depth}}$ is the depth coefficient
   - $X_{\text{vol},d}$ is the centered water volume
   - $X_{\text{depth},d}$ is the centered sampling depth
   - $\gamma_{s}$ is the site-specific random effect # EKJ note is this what we want?
   - $\delta_{\text{primer}[b]}$ is the primer-specific fixed effect, with method effects constrained such that:
   $\delta_{\text{MFU}} = 0$ and 
   $\delta_{\text{MV1}} \sim \text{Normal}(0, 1.7)$ # EKJ note may need to change these
   $\delta_{\text{DLP}} \sim \text{Normal}(0, 1.7)$
   
   The depth capture for a given site, volume, primer  is then modeled as:

   $Y_{\text{capture},d} \sim \text{Bernoulli}(Z_{s} \cdot p_{\text{capture},d})$

   Need to rewrite this bc it's only one Bernoulli trial in our case

   Conditional on capture in that volume x primer reaction, technical replicates are modeled as:

   $Y_{\text{detect},i} \sim \text{Bernoulli}(p_{\text{detect}} \cdot Y_{\text{capture},b[i]})$

   where $p_{\text{detect}}$ is the detection probability, drawn from a Beta prior:

   $p_{\text{detect}} \sim \text{Beta}(1,1)$
