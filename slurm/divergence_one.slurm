#!/bin/bash
#SBATCH --job-name=divergence_one
#SBATCH --time=0:30:00                   # Time limit hrs:min:sec
#SBATCH --output=divergence_1.out      # Standard output and error log

module load r/4.0.0
cd ~/r

Rscript --vanilla divergence_one.R
