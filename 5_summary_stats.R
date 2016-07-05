##- summary stats cluster sims
# rm(list=ls())
##---- libs ----
library(ggplot2)
library(reshape2)
library(cowplot)

##---- load data ----
cw_Baseline0 <- readRDS(file = "data/sim_ucsd_results2/list.sim.clus-outdeg.Baseline0.rds" )
cw <- cw_Baseline0[-c(2)]
# names(cw)
##- rbind table
cw_bind <- lapply(cw, function(x) do.call(rbind, x[1:2]))
# str(cw_bind)

##---- plan ----
##- number of clusters (ucsd)
##- number of distinct neighborhood cluster >= 2 ?
##- size of clusters
##- size of neighborhood
##- outdegrees distribution (offspring)
##- size of cluster by risk, stage, age
##- nghbrhood size by risk, stage, age

####---- n clusters ----
## number of different clusters (counting size 1)
sapply(cw[-1], function(x){
  summary(sapply(x, function(x) {
    length(unique(x$ClusterID) )
  }))
})

####---- proportion membership ----
# sapply(cw[-1], function(x){
#   summary(sapply(x, function(x) sum(x$binclus) / length(x$binclus)))
# })
sapply(cw_bind[-1], function(x) mean(x$binclus))

####---- mean size ----
# sapply(cw[-1], function(x){
#  summary(sapply(x, function(x) mean(x$size)))
# })
sapply(cw_bind[-1], function(x) summary(x$size))

####---- mean nbhsize ----
# sapply(cw[-1], function(x){
#   summary(sapply(x, function(x) mean(x$nbhsize)))
# })
sapply(cw_bind[-1], function(x) summary(x$nbhsize))

####---- mean outdegree ----
x <- do.call(rbind, cw[[1]])
summary(x$outdegree)
# summary(cw_bind[[1]]$outdegree)

####---- median size 2 ----
## stats of median size
sapply(cw[-1], function(x){
  summary(sapply(x, function(x) median(x$size)))
})

## stats of max size
sapply(cw[-1], function(x){
  summary(sapply(x, function(x) max(x$size)))
})

sapply(cw[-1], function(x) lapply(x, function(df) aggregate(df$size, by = list("stage" = df$stage), mean )))

##- correlation
sapply(cw_bind[-1], function(x) cor(x$size, x$nbhsize))
str(cw_bind)

##- for plots
a <- melt(cw_bind[-1], id.vars = 'id', measure.vars = c('size','nbhsize'))

str(a)
head(a)

g1 <- ggplot(a, aes(L1, value, color = variable)) + geom_boxplot() + facet_wrap(~variable) 
g1
