#!/bin/bash
#SBATCH --job-name=divergence
#SBATCH --time=0:30:00               # Time limit hrs:min:sec
#SBATCH --output=divergence.out      # Standard output and error log

module load r/3.6.3
cd ~/r
Rscript --vanilla divergence.R

## mv /gpfs/scratch/goldfk01/divergence.rda divergence.rda


