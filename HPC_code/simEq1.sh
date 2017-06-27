#!/usr/bin/env bash
#PBS -l walltime=01:30:00
#PBS -l select=1:ncpus=1:mem=8GB

cd $PBS_O_WORKDIR

module load gcc/5.4.0 R/3.3.2
module unload gcc/4.9.1 # ?
unset R_HOME

Rscript code/model1-simEqualStage0.R $PWD

