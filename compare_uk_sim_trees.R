## compare UK and sims dated tree for matching proportion into cluster (by adjusting rate)
# rm(list=ls())

library(ape)
library(phytools)
nh <- nodeHeights(uktree)
head(nh)
head(uktree$edge)
##---- load uk tree ----
if( any(grep("MacBook", Sys.info())) ){
  path.uk <- '../Box Sync/HPC/phylo-uk/BS_trees/MSM_B'
} else {
  path.uk <- '../Box/HPC/phylo-uk/BS_trees/MSM_B' # imac
}
l <- list.files(path.uk, full.names = TRUE)
uktree <- readRDS(l[1])[[2]]

##---- load sim ----
if( any(grep("MacBook", Sys.info())) ){
  path.sim <- '../Box Sync/HPC/simulations/model0-simulateBaseline0'
} else {
  path.sim <- '../Box/HPC/simulations/model0-simulateBaseline0' # imac
}
l2 <- list.files(path.sim, full.names = TRUE)
load(l2[1])
simtree <- bdt

## pick a clade