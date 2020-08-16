#!/bin/bash
#SBATCH --job-name=divergence_s
#SBATCH --partition=cpu_short
#SBATCH --time=0:30:00                   # Time limit hrs:min:sec
#SBATCH --output=binary_simple.out       # Standard output and error log

module load r/4.0.0
cd ~/r

Rscript --vanilla binary_simple.R
