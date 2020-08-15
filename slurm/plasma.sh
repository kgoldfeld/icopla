#!/bin/bash
#SBATCH --job-name=nc_10
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dw2625@nyulangone.org
#SBATCH --partition=cpu_short
#SBATCH --time=10:00:00
#SBATCH --mem=12GB
#SBATCH --output=nc_10.out

module load r/3.6.3
cd ~/R/x86_64-redhat-linux-gnu-library/file_danni
Rscript --vanilla interim_nc.R
