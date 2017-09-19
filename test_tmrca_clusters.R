## get clusters by tmrca < 5 years
## scatter plot size of clusters for each tip by respective method
# rm(list=ls())
library(phytools)
source("load_sims_files.R")
load(list.sims[[1]][1])

# random node
#set.seed(123)
set.seed(10)
i <- sample(1:bdt$Nnode, 1)
inode <- bdt$n + i

# a <- getDescendants(bdt, node =  inode)
# drop.tip(bdt, tip = bdt$tip.label[bdt$tip.label %in% a] )
sub_bdt <- extract.clade(bdt, node = inode, root.edge = 0) 

## ignoring DatedTree attributes
sub <- sub_bdt[c("edge","tip.label", "edge.length", "Nnode")]
class(sub) <- "phylo"
rm( list = setdiff(ls(), "sub") )

# mrca(phy = sub_bdt)
# library(phydynR)
# st <- bdt$sampleTimes[sub_bdt$tip.label]
# ss <- bdt$sampleStates[sub_bdt$tip.label, ]
# sub <- DatedTree(sub_bdt, st)
# sub$parentheights

plot(sub, label.offset = .5)
axisPhylo()
nodelabels()
tiplabels()

## loop tips
## aggregate tips which share a mrca within k years
tr <- sub
cbind(tr$edge, tr$edge.length)
mrcas <- mrca(sub)
clusters <- list()
k=5
i=1
j=1 #index cluster
ntips <- Ntip(tr)
count = 0 # number of tips explored
singleton <- c()

while (count < ntips) {
  # parent <- mrcas[i, i+1]
  ## get oldest ascendant possible
  oldest = 0
  current <- i # node evaluated
  print(paste("evaluating tip", current))
  repeat {
    # get mrca edge by edge
    #print(current)
    parent <- tr$edge[which(tr$edge[,2]== current), 1]
    #print(parent)
    d <- tr$edge.length[which(tr$edge[,2]== current)] # will break if root
    if( oldest + d > k ){
      parent <- current ## last valid mrca
      break
    } else {
      oldest <- oldest + d
      current <- parent
    }
  }
  print(paste("oldest mrca from tip", i, "is", parent))
  
 if (parent > ntips){ # is it a node ?
   #now parent is oldest ancestor within k
   clade <- extract.clade(tr, parent)
   depths <- node.depth.edgelength(clade)
   
   if ( sum(depths[1:Ntip(clade)] < k) < 2 ){ # singleton
     count <- count +1
     singleton <- c(singleton, tr$tip.label[parent])
     i <- which(tr$tip.label == setdiff(tr$tip.label, c(singleton, unlist(clusters)))[1])
     
   } else { # true cluster
     clusters[[j]] <- clade$tip.label[depths[1:Ntip(clade)] < k]
     names(clusters)[j] <- parent
     print(paste("from node", parent, "cluster =", clusters[[j]] ))
     count <- count + length( clusters[[j]] )
     i <- which(tr$tip.label == setdiff(tr$tip.label, c(singleton, unlist(clusters)) )[1]) # next non-evaluated (singleton or clustered) tip
     j=j+1
     if (length(i) == 0) break
   }
   
 } else { # if tip already farther than k, singleton
   count <- count +1
   singleton <- c(singleton, tr$tip.label[parent])
   i <- which(tr$tip.label == setdiff(tr$tip.label, c(singleton, unlist(clusters)))[1])
 }
}

colc <- lapply(clusters, function(x) which.edge(sub, x))
ecolor <- rep("darkgrey", dim(sub$edge)[1])
for (i in 1:length(clusters)){
  ecolor[ colc[[i]] ] <- palette()[i]
}
tcolor <- setNames(rep("darkgrey", length(sub$tip.label)), sub$tip.label)
for (i in 1:length(clusters)){
  tcolor[ clusters[[i]] ] <- palette()[i]
}
# tcolor <- rep(sample(palette(), length(clusters)), sapply(clusters, length))
# plot.phylo(sub, tip.color = color )
# axisPhylo()

plot.phylo(sub, edge.color = ecolor,tip.color = tcolor, edge.width = 2, label.offset = .5, cex = .8)
axisPhylo()
#edgelabels(round(sub$edge.length,1), bg="white", frame = "none", cex =.8)
?plot.phylo

# ## try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("ggtree")
# library(ggtree)
# tree <- groupClade(tr, node = 8)  
# ggtree(tree, aes(color=group, linetype=group)) + geom_tiplab()


###- tests
if (FALSE)
{
  which.edge(tr, c(1, 6))
  getMRCA(tr, c(1,4))
  findMRCA(tr, tips = c(1,2), type = "height") # height above the root of the MRCA
  fastHeight(tr, tr$tip.label[1], tr$tip.label[1])
  ?fastHeight
  ?dist.nodes()
  node.height(tr)
  node.depth.edgelength(tr)
}
