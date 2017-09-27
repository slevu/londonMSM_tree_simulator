## get clusters by tmrca < 5 years
## scatter plot size of clusters for each tip by respective method
# rm(list=ls())
options(max.print = 50)

##---- libs ----
library(phydynR)
library(RColorBrewer)
library(igraph)
library(Rcpp)
#library(phytools)

##---- rnd dated tree ----
rndDatedTree <- function(n){
  set.seed(123)
  rt <- rtree(n, br = rexp(n, rate = 1/5))
  st <- setNames(node.depth.edgelength(rt)[1:Ntip(rt)], rt$tip.label)
  drt <- DatedTree(rt, st)
  return(drt)
}
N <- 15
drt <- rndDatedTree(N)

if(Ntip(drt) < 50){
  plot(drt, label.offset = .6); axisPhylo(); nodelabels(); tiplabels()
  #edgelabels(round(drt$edge.length,2), bg="white",  cex =.8)
} else {
  plot.phylo(drt, show.tip.label = FALSE, type = "fan")
}

##---- function ----
##-- calculate clusters
mrca.cluster <- function(tree, threshold = 5, min.size = 2, plots = c("none", "graph", "tree", "all"), mrca_mat = NULL){
  ## create links between two tips if they share an ancestor within $k$ years for each
  require(ape, igraph, Rcpp)
  
  k <- threshold #5 ## limit of MRCA distance
  ##- matrix of mrca
  if (!is.null(mrca_mat)) {
    print("loading MRCA matrix")
    # saveRDS(m, file = MRCA.MAT)
    # mrca.mat = MRCA.MAT
    m <- readRDS(mrca_mat)
  } else {
    print("calculating matrix of mrca")
    st1 <- system.time( m <- mrca(tree) )
    print(st1[3])
  }
  
  ##- matrix of links
  print("calculating matrix of links")
  # heights <- max(node.depth.edgelength(tree)) - node.depth.edgelength(tree) # from most recent sample, without DatedTree object
  heights <- tree$heights
  if (FALSE) {
    ml <- m; ml[] <-  0
    st2 <- system.time(
      for (i in 1:dim(m)[1]) for (j in 1:dim(m)[2]){
        # m2[i,j] <- heights[m[i,j]] - heights[i] # time to mrca (for each daughter)
        ml[i,j] <- ifelse( i!=j && (heights[m[i,j]] - heights[i]) < k &&
                             (heights[m[j,i]] - heights[j]) < k, 1, 0)
      }
    )
    print(st2[3])
  } # R version
  ##- C++ loop ~500 times faster
  # sourceCpp("matrix_manip.cpp")
  cppFunction('NumericMatrix matrix_manip(NumericVector heights, IntegerMatrix m, double k){
    NumericMatrix ml(m.nrow(),m.ncol());
    for(int i = 0; i < m.nrow(); ++i){
      for(int j = 0; j < m.ncol(); ++j){
        if(i != j && (heights[m(i,j)-1] - heights[i]) < k && (heights[m(j,i)-1] - heights[j]) < k){
          ml(i,j) = 1;
        } else {
          ml(i,j) = 0;
        }
      }
    }
    return ml;
  }')
  st2 <- system.time(
    ml <- matrix_manip(heights, m, k)
  )
  rownames(ml) <- colnames(ml) <- rownames(m)
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

##---- on rnd tree ----
mclust0 <- mrca.cluster(tree = drt, threshold = 5, min.size = 2, plots = "all")

##---- load sim ----
source("load_sims_files.R")
sim <- list.sims[[1]][1]
sim.name <- names(sim)
load( sim )
MRCA_MAT <- paste0(path.results, '/', 'mrca.matrix.', sim.name, '.rds')

####---- subtree (if needed)
if (FALSE){
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
  
  drt <- DatedTree(sub, setNames(node.depth.edgelength(sub)[1:Ntip(sub)], sub$tip.label))
  #N <- Ntip(drt)
} ##- real subtree

##---- big ----
M_CLUST_LIST <- paste0(path.results, '/', 'mrca_clusters_list_', sim.name, '.rds')
if (FALSE){
  system.time(
    mclust <- lapply( list('1'=1,'2'=2,'5'=5,'10'=10), function(k) mrca.cluster(bdt, threshold = k, min.size = 0, plots = "none", mrca_mat = MRCA_MAT))
  ) # 87s
  # saveRDS(mclust, file = M_CLUST_LIST)
} else {
  mclust <- readRDS(file = M_CLUST_LIST)
}


##---- helpers
##- list of cluster in dataframes
list.clus.in.df <- function(lst) data.frame(SequenceID = as.integer(do.call(c, lst)), ClusterID = as.integer(rep(names(lst), sapply(lst, length))) )
##- add size to seq / cluster data.frame
get.cluster.size <- function(df){
  freqClust <- as.data.frame(table(df$ClusterID), stringsAsFactors = FALSE)
  newdf <- merge(x = df, y = freqClust, 
                 by.x = "ClusterID", by.y = "Var1", 
                 all.x = TRUE, sort = FALSE)
  return(newdf[order(newdf$SequenceID),])
}

mrca_clusters <- lapply(mclust, function(x) get.cluster.size(list.clus.in.df(x)))

##---- load hivclustering ----
x <- readRDS(paste0(path.results, "/", "list.hivclust.sim.Baseline0.rds"))
hivclust1 <- lapply(x, function(a) a[[sim.name]])
h1 <- lapply(hivclust1, get.cluster.size)
##- add singletons
hivclustering_clusters <- lapply(h1, function(d) {
  add <- data.frame(ClusterID = (nrow(d)+1):Ntip(bdt), 
             SequenceID = setdiff(as.integer(bdt$tip.label), d$SequenceID),
             Freq = 1L)
  newdf <- rbind(d, add)
  return(newdf[order(newdf$SequenceID),])
  })

##---- distribution cluster size ----
sapply(mrca_clusters, function(x) summary(x$Freq))
sapply(hivclustering_clusters, function(x) summary(x$Freq))

thr_mrca <-  "5" # "10" #
thr_ucsd <- "0.015" # "0.05"  # 

# par(mfrow=c(1,2))
# hist(mrca.clusters[[thr_mrca]]$Freq)
# hist(hivclustering.clusters[[thr_ucsd]]$Freq)
# dev.off()
# without singletons
par(mfrow=c(1,2))
hist(mrca_clusters[[thr_mrca]]$Freq[mrca_clusters[[thr_mrca]]$Freq > 1], xlab = "size", main = "mrca method, size > 1")
hist(hivclustering_clusters[[thr_ucsd]]$Freq[hivclustering_clusters[[thr_ucsd]]$Freq > 1], xlab = "size", main = "ucsd method, size > 1")
#dev.off()

##---- scatter ----
df <- merge(x = hivclustering_clusters[[thr_ucsd]][, c("SequenceID", "Freq")],
            y = mrca_clusters[[thr_mrca]][, c("SequenceID", "Freq")], by = "SequenceID" )
df_no_single <- df[df$Freq.x > 1 & df$Freq.y > 1, ]

x <- as.data.frame(table(df_no_single[, 2:3]), stringsAsFactors = FALSE)
par(mfrow=c(1,1))
radius <- sqrt(x$Freq/pi)
symbols(x$Freq.x, x$Freq.y, circles = radius, inches = 0.25, fg = "white",
        bg = "red",
        xlab = "Size (hivclustering)", # paste("hivclustering cluster size -", thr_ucsd),
        ylab = "Size (tMRCA)", # paste("tMRCA cluster size -", thr_mrca),
        main = "") # Sized by Freq

corcoef <- cor.test(df_no_single$Freq.x, df_no_single$Freq.y)$estimate

##---- geom_hex ----
library(ggplot2)
d <- ggplot(data= df_no_single, aes(Freq.x, Freq.y))
d + geom_hex()# stat_bin_hex() #geom_hex()

##---- stop ----

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

##- test C++ version
sourceCpp("matrix_manip.cpp")
system.time(ml <- matrix_manip(heights, m, k)) # 8s
