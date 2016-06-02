##--- to run on HPC ---

#!/usr/bin/env bash
#PBS -l walltime=00:20:00
#PBS -l select=1:ncpus=1:mem=500mb

##$-V                                      		# Using all environment Variables

module load R/3.2.0

Rscript /work/slevu/simulations/assort_matrix0.R $file
