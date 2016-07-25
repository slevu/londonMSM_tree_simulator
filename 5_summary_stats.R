##- summary stats cluster sims
# rm(list=ls())
##---- libs ----
library(ggplot2)
library(reshape2)
library(cowplot)

##---- load data ----
cw_Baseline0 <- readRDS(file = "data/sim_ucsd_results2/list.sim.clus-outdeg.Baseline0.rds" )
cw <- cw_Baseline0[c('SA', '0.001', '0.005', '0.015', '0.05')]
# names(cw)
##- rbind table
cw_bind <- lapply(cw, function(x) do.call(rbind, x[1:10]))
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
sapply(cw_bind[-1], function(x) sd(x$size))

####---- mean nbhsize ----
# sapply(cw[-1], function(x){
#   summary(sapply(x, function(x) mean(x$nbhsize)))
# })
sapply(cw_bind[-1], function(x) summary(x$nbhsize))
sapply(cw_bind[-1], function(x) sd(x$nbhsize))

####---- mean outdegree ----
x <- do.call(rbind, cw[[1]])
summary(x$outdegree)
sd(x$outdegree)
summary(x$indegree)
sd(x$indegree)
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

##---- plots ----
##- boxplot with base R
##- strip of 3 plots without outliers

bp_base <- function(ls , var , lbl, tran = identity){
  
  par(mfrow=c(1,3))    
  
  lapply(ls, function(x){
    if ( any(grepl('outdegree', names(x))) ) {
      
      boxplot(x[, 'outdegree'] ~ x[, var], outline = FALSE )
      title(main = 'Out-degree',  font.main = 1)
      
      } else {
        
      boxplot(x[, 'size'] ~ x[, var], outline = FALSE)
        title(main = 'Cluster size', xlab = lbl, cex.lab = 1.2,  font.main = 1)
        
      boxplot(x[, 'nbhsize'] ~ x[, var], outline = FALSE)
      title(main = 'Neighborhood size',  font.main = 1)
    }
  })
}
  
##---- bp base ----
p <- bp_base(ls = cw_bind[c('SA', '0.015')], var = 'risk', lbl = 'Risk level')
p <- bp_base(ls = cw_bind[c('SA', '0.015')], var = 'age', lbl = 'Age category')
p <- bp_base(ls = cw_bind[c('SA', '0.015')], var = 'stage', lbl = 'Infection stage')

