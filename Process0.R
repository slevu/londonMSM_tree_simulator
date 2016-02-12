####---- run model0.R ----####
#### rm(list=ls())
####---- process o and tree ----####
## Load newest Rdata
l <- list.files(pattern="*.Rdata") # list.files(pattern="Rdata$") list.files(pattern="out")
load(l[length(l)])
str(o)
# require(deSolve)
# ?ode
##- ode(y, times, func, parms, method)


##- UPGMA
## get distances
##- matrix first
mtree <- cophenetic.phylo(tree)
mtree[1:6, 1:6]
str(tree)
##- into distances
dtree <- as.dist(mtree)
head(dtree)
## normalize
x <- dtree/ max(dtree)
head(x)


summary(dtree)
hist(dtree, breaks = 50, xlab = "distance", ylab = "frequency", main = "", col = "grey")
## same shape but very different values than UPGMA !

## cluster
hc <- hclust(x, method = "average") # UPGMA


breaks <- 10
up <- round(mean(x),1) # cut up to the mean
cluster <- cutree(hc,  h = seq(up / breaks, up, by = up / breaks))

# colnames(cluster) <- paste("c",colnames(cluster),sep='')
# colnames(cluster) <- paste("c", 1:8, sep='')
head(cluster)

##- Calculate size(=Freq) of each cluster across different threshold
freqClust <- apply(cluster, 2, function(x) as.data.frame(table(x))) # list

##- number of different clusters by threshold
sapply(freqClust, function(x) dim(x)[1])
##- cluster size
sapply(freqClust, function(x) summary(x$Freq))
##- percentiles
sapply(freqClust, function(x) round(quantile(x$Freq, probs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.95, 0.99, 1))))

head(freqClust[[1]])
hist(freqClust[[4]]$Freq)
##-- check assortativity 
## first use EdgeList at each level
## then function AssortMix


##- UCSD cluster
##
##
### ex cscaling
x <- matrix(1:10, ncol = 2)
(centered.x <- scale(x, scale = FALSE))
cov(centered.scaled.x <- scale(x)) # all 1