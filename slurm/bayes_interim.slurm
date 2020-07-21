#!/bin/bash
#SBATCH --job-name=bayes_interim
#SBATCH --mail-type=END,FAIL              # send email if the job end or fail
#SBATCH --mail-user=dw2625@nyulangone.org
#SBATCH --partition=cpu_short
#SBATCH --time=1:00:00                    # Time limit hrs:min:sec
#SBATCH --output=bayes_interim.out        # Standard output and error log

module load r/3.6.3
cd ~/r
Rscript --vanilla bayes_interim.R


