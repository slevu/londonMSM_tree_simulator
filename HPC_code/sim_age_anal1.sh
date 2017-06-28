##--- to run on HPC ---
#!/usr/bin/env bash
#PBS -l walltime=00:20:00
#PBS -l select=1:ncpus=1:mem=500mb

cd $PBS_O_WORKDIR

module load R

Rscript code/assort_matrix1.R $file $PWD
