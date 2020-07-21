#!/bin/bash
#SBATCH --job-name=bayes_interim
#SBATCH --time=1:00:00                    # Time limit hrs:min:sec
#SBATCH --output=bayes_interim.out        # Standard output and error log

module load r/3.6.3
cd ~/r
Rscript --vanilla bayes_interim.R


