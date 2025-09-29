#!/bin/bash

#SBATCH --account=coenv
#SBATCH --partition=cpu-g2 # Queue selection
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --mem=50G
#SBATCH --time=48:00:00 # Time limit hrs:min:sec\
#SBATCH --job-name=3Ddist_Q3_GAM # Job name
#SBATCH --mail-type=ALL # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=avancise@uw.edu # Where to send mail\
#SBATCH --chdir=/gscratch/coenv/avancise
#SBATCH --export=all
#SBATCH --output=3Ddist_Q3_GAM.out
#SBATCH --error=3Ddist_Q3_GAM.err

apptainer run --bind /gscratch /gscratch/coenv/shared/Apptainer/tidyverse_latest.sif Rscript /gscratch/coenv/avancise/3DDistribution/04_Model_with_GAMs_Q3.R