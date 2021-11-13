#!/bin/bash
#SBATCH --job-name=i_mult_1
#SBATCH --mail-type=END,FAIL                         # send email if the job end or fail
#SBATCH --mail-user=danni.wu@nyulangone.org
#SBATCH --partition=cpu_short
#SBATCH --time=12:00:00                               # Time limit hrs:min:sec
#SBATCH --output=multiple.out                        # Standard output and error log

module load r/4.0.3
cd /gpfs/data/troxellab/danniw/r

Rscript --vanilla multiple_v2.R