####---- run model0.R ----####
#### rm(list=ls())
library(ape)
library(ggplot2)
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
if (F){
#- function of heights
nheights <- 10 # number of threshold
up <- round(mean(simx), 1) # cut up to the mean
# breaks <-  seq(1/nbreaks, 1-(1/nbreaks), by = 1/nbreaks)
simclus <- cutree(simhc, h = seq(up / nheights, up, by = up / nheights) ) # h = breaks
up <- round(mean(ukx),1) # cut up to the mean
ukclus <- cutree(ukhc,  h = seq(up / nheights, up, by = up /nheights) )
# rm(simx, ukx)
}

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

##- number of different clusters by threshold # if number varies !
sapply(simfreqClust, function(x) dim(x)[1])
sapply(ukfreqClust, function(x) dim(x)[1])

##- cluster size
sapply(simfreqClust, function(x) summary(x$Freq))
sapply(ukfreqClust, function(x) summary(x$Freq))
##- percentiles
# sapply(freqClust, function(x) round(quantile(x$Freq, probs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.95, 0.99, 1))))
# 

##- distr of cluster sizes: log(x) and Y untransformed
par(mfcol=c(2,5))
for (i in 1:length(kgroups)){
  h <- hist(log(ukfreqClust[[i]]$Freq), 
       main = paste("uk", names(ukfreqClust)[i]),
       xlab = "log(size)")
  hist(log(simfreqClust[[i]]$Freq), 
       main = paste("sim", names(simfreqClust)[i]),
       xlab = "log(size)")
}

##- distr of cluster sizes: log(x) and log(y)
par(mfcol=c(2,5))
for (i in 1:length(kgroups)){
  h <- hist(log(ukfreqClust[[i]]$Freq), plot = F)
  h$counts <- log1p(h$counts) # log(y)
  plot(h, ylab = "log(Freq)", 
       main = paste("uk", names(ukfreqClust)[i]),
       xlab = "log(size)")
  
  h <- hist(log(simfreqClust[[i]]$Freq), plot = F)
  h$counts <- log1p(h$counts) # log(y)
  plot(h, ylab = "log(Freq)",
          main = paste("sim", names(simfreqClust)[i]),
          xlab = "log(size)")
  
}

##- QQ plot
par(mfcol=c(2,5))
for (i in 1:length(kgroups)){
  qqplot(ukfreqClust[[i]]$Freq, 
         simfreqClust[[i]]$Freq,
         main = names(ukfreqClust)[i],
         xlab = "uk", ylab = "sim")
  qqline(y, col = 2)
  qqplot(log(ukfreqClust[[i]]$Freq), 
         log(simfreqClust[[i]]$Freq),
         main = names(ukfreqClust)[i],
         xlab = "log(uk)", ylab = "log(sim)")

}




##- ggplot: distributions of cluster sizes for different number of clusters
##- dataframe: 
df <- data.frame()
for (i in 2:length(kgroups)){
  df <- rbind(df,
              rbind( cbind(gr = "sim", 
                           k = names(simfreqClust)[i], 
                           simfreqClust[[i]]) ,
                     cbind(gr = "uk", 
                           k = names(ukfreqClust)[i],
                           ukfreqClust[[i]]) )
              )
}
head(df)
str(df)


##- UPGMA comparison of cluster sizes
g <- ggplot(data = df, aes(Freq))
g +  geom_histogram(origin = 0, binwidth = .5) + 
  scale_x_continuous(breaks=c(0,1,3,10,30,100,300,1000,5000), 
                     trans = "log1p",
                     name = "Cluster size") +
  scale_y_continuous(breaks=c(0, 30, 100, 300, 1000),
                     trans = "log1p",
                     name = "Number of clusters") +
  facet_grid(k  ~ gr, 
             scales = "free",
             space = "free") +
  theme_minimal()
 # theme_bw()


##-- quantile-quantile

##- data 
str(o)
dim(o)
str(tree)
o[990:1000, 120:122]
head(tree$tip.label)
##-- check assortativity 
## first use EdgeList at each level
## then function AssortMix

##- false discovery rate: trnamission rates are equal according to cluster or SA ?
## on 100 trees, find significant association
## 
## study variance in cluster size, outdegrees changes with respect to the risk levels

##- UCSD cluster
##
##
y <- rt(200, df = 5)
qqnorm(y); qqline(y, col = 2)
qqplot(y, rt(300, df = 5))

qqnorm(precip, ylab = "Precipitation [in/yr] for 70 US cities")