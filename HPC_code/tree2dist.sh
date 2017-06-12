#!/usr/bin/env bash
#PBS -l walltime=00:10:00
#PBS -l select=1:ncpus=1:mem=4gb

cd $PBS_O_WORKDIR

## Compute cophenetic distance on year tree contained 
## in simulation results with mutation rate MU
## ex: SIM=data/simulations2/model0-simulateBaseline0/824.RData
BN=$(basename $SIM) # 824.RData
DN=$(dirname $SIM) # data/simulations2/model0-simulateBaseline0
DIROUT=${DN/%/-coph_distances} # data/simulations2/model0-simulateBaseline0-coph_distances
mkdir -p $DIROUT
OUT="$DIROUT/$BN" # data/simulations2/model0-simulateBaseline0-coph_distances/824.RData

module load gcc/5.4.0 R/3.3.2
module unload gcc/4.9.1 # ?
unset R_HOME

## parameters: mutation rate, distance threshold
MU=0.0043 #(Berry et al.)
DL=0.05 # don't report distances above that
SL=1000

Rscript code/tree2genDist.R $SIM $MU $DL $SL $OUT $PWD

