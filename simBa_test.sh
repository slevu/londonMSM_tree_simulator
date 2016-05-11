#!/usr/bin/env bash
#PBS -l walltime=00:02:00
#PBS -l select=1:ncpus=1:mem=50mb

module load R/3.2.0
print(getwd())
##Rscript model0-simulateBaseline0_slv.R

