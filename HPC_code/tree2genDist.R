#!/usr/bin/Rscript --slave
args <- commandArgs(TRUE)
SIM <- args[1]
MU <- as.numeric(args[2])
DL <- as.numeric(args[3])
SL <- as.numeric(args[4])
OUT <- args[5]
wd <- args[6]

library(ape)

setwd(wd)

load(SIM)
source("code/fn_tree2CophDist.R")
D <- tree2CophDist(bdt, mu = MU, seqlength = SL, dlim = DL)
save(D, file = OUT)
