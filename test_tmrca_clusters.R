## get clusters by tmrca < 5 years
## scatter plot size of clusters for each tip by respective method
library(phytools)
source("load_sims_files.R")
load(list.sims[[1]][1])

# random node
set.seed(123)
i <- sample(1:bdt$Nnode, 1)
inode <- bdt$n + i

# a <- getDescendants(bdt, node =  inode)
# drop.tip(bdt, tip = bdt$tip.label[bdt$tip.label %in% a] )
sub_bdt <- extract.clade(bdt, node = inode, root.edge = 0) 
plot(sub_bdt)
axisPhylo()
nodelabels()
mrca(phy = sub_bdt)

library(phydynR)
st <- bdt$sampleTimes[sub_bdt$tip.label]
ss <- bdt$sampleStates[sub_bdt$tip.label, ]
sub <- DatedTree(sub_bdt, st)
sub$parentheights
