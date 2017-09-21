## get clusters by tmrca < 5 years
## scatter plot size of clusters for each tip by respective method
# rm(list=ls())
options(max.print = 50)
# ###- tests DatedTrees
library(phydynR)
library(RColorBrewer)
library(igraph)
#library(phytools)

if (FALSE) {
  rndDatedTree <- function(n){
  set.seed(123)
  rt <- rtree(n)
  st <- setNames(node.depth.edgelength(rt)[1:Ntip(rt)], rt$tip.label)
  drt <- DatedTree(rt, st)
  return(drt)
}
N <- 4
drt <- rndDatedTree(N)
} ##- random DatedTree

source("load_sims_files.R")
sim <- list.sims[[1]][1]
sim.name <- names(sim)
load( sim )

##- get subtree at random node
#set.seed(123)
set.seed(10)
tr <- bdt
random.subtree <- function(tr, min.subtree = NULL, max.subtree = Inf, iter = 10){
  i <- sample(Ntip(tr) + 1:tr$Nnode, 1)
  subtr <- extract.clade(tr, node = i, root.edge = 0)
  print(paste("subtree from node", i, "with", Ntip(subtr), "tips"))
  if (!is.null(min.subtree)) {
    n = 0
    repeat {
      subtr <- random.subtree(tr)
      n <- n + 1
      print(n)
      if ((Ntip(subtr) >= min.subtree & Ntip(subtr) <= max.subtree ) || n > iter) break
    }
  }
  return(subtr)
}
test <- random.subtree(tr, min.subtree = 20)

## ignoring DatedTree attributes
sub <- test[c("edge","tip.label", "edge.length", "Nnode")]
class(sub) <- "phylo"
rm( list = setdiff(ls(), "sub") )



##- real tree
drt <- DatedTree(sub, setNames(node.depth.edgelength(sub)[1:Ntip(sub)], sub$tip.label))
#N <- Ntip(drt)

if(Ntip(drt) < 50){
  plot(drt, label.offset = .1); axisPhylo(); nodelabels(); tiplabels()
  #edgelabels(round(drt$edge.length,2), bg="white",  cex =.8)
} else {
  plot.phylo(drt, show.tip.label = FALSE, type = "fan");  axisPhylo()
}

###--- calculate clusters
# tree = drt; threshold = 1; min.size = 0; plots = "all"
mrca.cluster <- function(tree, threshold = 5, min.size = 2, plots = c("none", "graph", "tree", "all")){
  ## create links between two tips if they share an ancestor within $k$ years for each
  require(ape, igraph)
  
  ##- matrix of mrca
  k <- threshold #5 #1.2
  print("calculating matrix of mrca")
  st1 <- system.time( m <- mrca(tree) )
  print(st1[3])
  ##- matrix of links
  print("calculating matrix of links")
  ml <- m; ml[] <-  0
  heights <- tree$heights
  st2 <- system.time(
  for (i in 1:dim(m)[1]) for (j in 1:dim(m)[2]){
    # m2[i,j] <- tree$heights[m[i,j]] - tree$heights[i] # time to mrca (for each daughter)
    ml[i,j] <- ifelse( i!=j && (heights[m[i,j]] - heights[i]) < k &&
                         (heights[m[j,i]] - heights[j]) < k, 1, 0)
  }
  )
  print(st2[3])
  ##- create graph
  print("create graph")
  st3 <- system.time(
    g <- graph_from_adjacency_matrix(ml, mode = "undirected")
  )
  print(st3[3])
  if (plots[1] %in% c("all", "graph")) plot.igraph(g, vertex.size = min(15, 15/log10(length(V(g)))), vertex.label.color="black", vertex.color="lightgrey", edge.color="black", vertex.label.cex = min(1,1/log10(length(V(g)))) )
  print("compute components")
  st4 <- system.time(
    clu <- components(g)
  )
  print(st4[3])
  clusters <- groups(clu)
  clusters <- clusters[clu$csize >= min.size] # remove singleton
  names(clusters) <- 1:length(clusters)
  
  ####- plot tree
  if (plots[1] %in% c("tree", "all")){
    
    plot.tree.cluster <- function(tr = tree, cluster.list = clusters){
      ##- palette
      require(RColorBrewer)
      darkcols <- brewer.pal(8, "Dark2")
      pal <- colorRampPalette(darkcols)
      colours <- pal(length(cluster.list))
      ##-  edge color
      colc <- lapply(clusters, function(x) {
        which.edge(tr, match(x, tr$tip.label ) )
      })
      ecolor <- rep("darkgrey", dim(tr$edge)[1])
      for (i in 1:length(clusters)){
        ecolor[ colc[[i]] ] <- colours[i]
      }
      ##- tip color
      tcolor <- setNames(rep("darkgrey", length(tr$tip.label)), tr$tip.label)
      for (i in 1:length(clusters)){
        tcolor[ clusters[[i]] ] <- colours[i]
      }
      ##- edge width
      ewidth <- ifelse(ecolor == "darkgrey", 1, 2)
      if (Ntip(tr) < 50){
        plot.phylo(tr, edge.color = ecolor,tip.color = tcolor, edge.width = 2, label.offset = .2, cex = .8)
        axisPhylo()
      } else {
        plot.phylo(tr, edge.color = ecolor, edge.width = ewidth, show.tip.label = FALSE, type = "fan")
      }
    }
    
    plot.tree.cluster()
  }
  return(clusters)
}


mclust <- mrca.cluster(tree = drt, threshold = 1, min.size = 0, plots = "all")

##- big
system.time(
  mclust <- lapply(c(5), function(k) mrca.cluster(bdt, threshold = k, min.size = 0, plots = "none"))
  )
saveRDS(mclust, file = paste0(path.results, '/', 'mrca.cluster.5years.', sim.name, '.rds'))

##- in a dataframe
mclust1 <- data.frame(SequenceID = as.integer(do.call(c, mclust)), ClusterID = as.integer(rep(names(mclust), sapply(mclust, length))) )


##- helper: add size to seq / cluster data.frame
get.cluster.size <- function(df){
  freqClust <- as.data.frame(table(df$ClusterID), stringsAsFactors = FALSE)
  newdf <- merge(x = df, y = freqClust, 
                 by.x = "ClusterID", by.y = "Var1", 
                 all.x = TRUE, sort = FALSE)
  return(newdf[order(newdf$SequenceID),])
}

mrca.clusters <- get.cluster.size(mclust1)

##- compare with hivclustering
x <- readRDS(paste0(path.results, "/", "list.hivclust.sim.Baseline0.rds"))
hivclust1 <- lapply(x, function(a) a[[sim.name]])
hivclustering.clusters <- lapply(hivclust1, get.cluster.size)

head(hivclustering.clusters[[1]][ order(hivclustering.clusters[[1]]$SequenceID), ])
sum(hivclustering.clusters[[1]]$SequenceID %in% bdt$tip.label)

d <- hivclustering.clusters[[1]]
##- add singletons
l <- lapply(hivclustering.clusters, function(d) {
  add <- data.frame(ClusterID = (nrow(d)+1):Ntip(bdt), 
             SequenceID = setdiff(as.integer(bdt$tip.label), d$SequenceID),
             Freq = 1L)
  newdf <- rbind(d, add)
  return(newdf[order(newdf$SequenceID),])
  })

setdiff(1:5, 4:8)
                 
###- SURPLUS -###
##- way to get heights from non DatedTree
heights <- max(node.depth.edgelength(drt)) - node.depth.edgelength(drt) # from most recent sample

if(FALSE){
  mrca2 <- function(tree){
    n <- Ntip(tree)
    m <- matrix(0, nrow = n, ncol = n, dimnames = list(tree$tip.label, tree$tip.label))
    for (i in tree$tip.label) for (j in tree$tip.label){
      m[i,j] <- fastMRCA(drt, i, j)
    }
    return(m)
  }
  library(microbenchmark)
  microbenchmark(mrca(drt), mrca2(drt))
} ##- test mrca


## an object that you want to recreate
m2 <- matrix(1:4,2,2)
## use capture.output to save structure as a string in a varible
xx <- capture.output(dput(m2))

## recreate the object 
m2_ <- eval(parse(text=xx))
image(z=m2_,col=rainbow(4))
dput(m)
