####---- run model0.R ----####
#### rm(list=ls())
####---- process o and tree ----####
## Load newest Rdata
l <- list.files(pattern="*.Rdata") # list.files(pattern="Rdata$") list.files(pattern="out")
load(l[length(l)])
# str(o)
# require(deSolve)
# ?ode
##- ode(y, times, func, parms, method)

####---- load ExaML tree ----####
t <- read.tree(file = "../phylo-uk/data/ExaML_result.subUKogC_noDRM.finaltree.000")
## drop OG
og <- c("Ref1", "Ref2", "Ref3", "Ref4", "Ref5", "Ref6", "HXB2")
t <- drop.tip(t, og ) 
# to drop all tips except wanted 
# t_og <- drop.tip(t, which(!t$tip.label %in% og))
# 
####---- compare----#####
####  cluster size to real data. Need to have same number of clusters...

## get distances
##- matrix first into distances
library(ape)
library(ggplot2)
tree
simtree <- as.dist(cophenetic.phylo(tree))
uktree <- as.dist(cophenetic.phylo(t))
head(simtree)
head(uktree)

## normalize
simx <- simtree / (max(simtree) - min(simtree))
ukx <- uktree / (max(uktree) - min(uktree))
# rm(simtree, uktree)

##- histogram distances
# summary(x)
hist(simx, breaks = 50, xlab = "distance", ylab = "frequency", main = "", col = "grey")
hist(ukx, breaks = 50, xlab = "distance", ylab = "frequency", main = "", col = "grey")

##- cluster UPGMA
simhc <- hclust(simx, method = "average") # UPGMA
ukhc <- hclust(ukx, method = "average")

##- cut
#- function of heights
nheights <- 10 # number of threshold
up <- round(mean(simx), 1) # cut up to the mean
# breaks <-  seq(1/nbreaks, 1-(1/nbreaks), by = 1/nbreaks)
simclus <- cutree(simhc, h = seq(up / nheights, up, by = up / nheights) ) # h = breaks
up <- round(mean(ukx),1) # cut up to the mean
ukclus <- cutree(ukhc,  h = seq(up / nheights, up, by = up /nheights) )
# rm(simx, ukx)

#- function of k groups
kgroups <- c(100, 500, 1000, 2000, 5000)
simclus <- cutree(simhc, k = kgroups ) # h = breaks
ukclus <- cutree(ukhc,  k = kgroups )

# colnames(cluster) <- paste("c",colnames(cluster),sep='')
# colnames(cluster) <- paste("c", 1:8, sep='')
head(simclus)
head(ukclus)

##- Calculate size(=Freq) of each cluster across different threshold
simfreqClust <- apply(simclus, 2, function(x) as.data.frame(table(x))) # list
ukfreqClust <- apply(ukclus, 2, function(x) as.data.frame(table(x)))
str(simfreqClust)
head(simfreqClust[[1]])
##- number of different clusters by threshold
sapply(simfreqClust, function(x) dim(x)[1])
sapply(ukfreqClust, function(x) dim(x)[1])

##- cluster size
sapply(simfreqClust, function(x) summary(x$Freq))
sapply(ukfreqClust, function(x) summary(x$Freq))
##- percentiles
# sapply(freqClust, function(x) round(quantile(x$Freq, probs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.95, 0.99, 1))))

##- plot: distributions of cluster sizes for different number of clusters
##- dataframe: 
df <- rbind( cbind(gr = "sim", simfreqClust[[3]]) ,
             cbind(gr = "uk", ukfreqClust[[3]]) )

g <- ggplot(data = df, aes(x = Freq, group = gr))
g + geom_histogram(binwidth = 10) + facet_grid( ~ gr)

bp <- ggplot(data = df, aes(gr, log(Freq)))
bp + geom_boxplot()  


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