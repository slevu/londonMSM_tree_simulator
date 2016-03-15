##---- import ----
## from Rcolgem
rm(list=ls())
w0 <- read.csv("W0.sim.csv")
head(w0)
## add individual explanatory variates
df <- readRDS("demo.rds")
####---- Assort function ----
names(df)

##-- SA results
source("functions.r")
mix_age <- AssortMix(edgelist = w0, 
                     listID = df, dfvar = df, 
                     var = "age", 
                     pid = "patient", sid = "patient", 
                     graph = 1, nqt = NA, method = "SA")

mix_stage <- AssortMix(edgelist = w0, 
                       listID = df, dfvar = df, 
                       var = "stage", 
                       pid = "patient", sid = "patient", 
                       graph = 1, nqt = NA, method = "SA")

mix_risk <- AssortMix(edgelist = w0, 
                     listID = df, dfvar = df, 
                     var = "risk", 
                     pid = "patient", sid = "patient", 
                     graph = 1, nqt = NA, method = "SA")

mix_time <- AssortMix(edgelist = w0, 
                      listID = df, dfvar = df, 
                      var = "time", 
                      pid = "patient", sid = "patient", 
                      graph = 1, nqt = 10, method = "SA")

##-- Cluster results
listclus <- readRDS( file = "data/listclus_sim.rds")
head(listclus[[1]])
prepare.for.edgelist <- function(x){
  a <- as.vector(x[, "ClusterID"])
  names(a) <- x[, "id"]
  return(a)
}
ll <- lapply(listclus, prepare.for.edgelist)
head(ll[[2]])
lel <- lapply(ll, EdgeList)
names(lel)
head(lel[[2]])
## symetric transmission and add weigths (simply 1 here)
W <- lapply(lel, function(x) cbind(rbind(x, x[ , c(2,1)]), W = 1))
str(W)
head(W[[1]])
W_clust1 <- W[[2]]

##--- essai function AssortMix ---
##- load files
#load("../data/sub.RData") ## df AND s
#head(df)

age_mixing <- AssortMix (edgelist = W_clust1, 
                         listID = df, dfvar = df, 
                         var = "age", 
                         pid = "patient", sid = "patient",
                         nqt = NA, graph = 1, method = "clust")

stage_mixing <- AssortMix (edgelist = W_clust1, 
                           listID = df, dfvar = df, 
                           var = "stage", 
                           pid = "patient", sid = "patient",
                           nqt = NA, graph = 1, method = "clust")

risk_mixing <- AssortMix (edgelist = W_clust1, 
                           listID = df, dfvar = df, 
                           var = "risk", 
                           pid = "patient", sid = "patient",
                           nqt = NA, graph = 1, method = "clust") 

time_mixing <- AssortMix (edgelist = W_clust1, 
                          listID = df, dfvar = df, 
                          var = "time", 
                          pid = "patient", sid = "patient",
                          nqt = 10, graph = 1, method = "clust")

####--- calculate assortativity coefficient
##- igraph
library(igraph)

###--- Cluster
## edge list
wdf <- as.data.frame(W_clust1)
head(wdf)
#- vertices data with seqindex first and in the edge list
verts <- df[ df$patient %in% wdf$from | df$patient %in% wdf$to,
             c("patient", "age", "stage", "risk")]

testg <- graph.data.frame(wdf, directed = FALSE, vertices = verts)
testg

assortativity(testg, types1 = V(testg)$age, directed = FALSE)
assortativity(testg, types1 = V(testg)$stage, directed = FALSE)
assortativity(testg, types1 = V(testg)$risk, directed = FALSE)

###--- SA
## edge list
str(w0)

#- vertices data with seqindex first and in the edge list
verts <- df[ df$patient %in% w0$donor | df$patient %in% w0$recipient,
             c("patient", "age", "stage", "risk")]

testg <- graph.data.frame(w0, directed = FALSE, vertices = verts)

testg

assortativity(testg, types1 = V(testg)$age, directed = FALSE)
assortativity(testg, types1 = V(testg)$stage, directed = FALSE)
assortativity(testg, types1 = V(testg)$risk, directed = FALSE)


####---- surplus ----####
# from http://stackoverflow.com/questions/21727813/calculating-assortativity-in-igraph
set.seed(123)
A = data.frame(rnorm(10),rnorm(10),rnorm(10),rnorm(10))
inv<-cor(t(A))
inv[inv<0.5] <- 0
inv[inv==1] <- 0
g1 <- graph.adjacency(inv, mode = "undirected", diag=FALSE, weighted=TRUE)
as_edgelist(g1)

V(g1)$foo <- sample(1:3, replace=TRUE, vcount(g1))
assortativity.nominal(g1, types = V(g1)$foo)
###########
