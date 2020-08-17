#!/bin/bash
#SBATCH --job-name=binary_simple
#SBATCH --partition=cpu_short
#SBATCH --mem=12G
#SBATCH --time=0:30:00                   # Time limit hrs:min:sec
#SBATCH --output=binary_simple.out       # Standard output and error log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

module load r/3.6.3
cd /gpfs/home/goldfk01/r

Rscript --vanilla binary_simple.R
