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
cw_bind <- lapply(cw, function(x) do.call(rbind, x[1:5]))
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
summary(x$indegree)
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
##- long table with value = od, size ...; method = SA, CL, NB; thr = ..
##- for plots
# source('functions.R')
molten <- function(cw_bind){
  .sa <- melt( cw_bind[[1]], id.vars = c('id', 'stage', 'age', 'risk'),  measure.vars = 'outdegree')
  # head(.sa)
  .thr <- melt(cw_bind[-1], id.vars = c('id', 'stage', 'age', 'risk'), measure.vars = c('size', 'nbhsize'))
  # head(.thr)
  tot <- rbind( cbind(.sa, "L1" = 'NA'), .thr )
  # head(tot)
  return(tot)
}

tot <- molten(cw_bind)
tot_2 <- molten(cw_bind[c('SA', '0.015')])



##- function plot
gg <- function(x , var , lbl, tran = identity){
  g1 <- ggplot(x, aes_string(var, 'tran(value)', color = 'variable')) + 
    geom_boxplot() + 
    facet_wrap( ~ variable + L1, scales = "free_y", ncol = 3, labeller= labeller(variable = c(outdegree = "Out-degree", size = "Cluster size", nbhsize = "Neighborhood size")))  + 
    xlab(lbl) + ylab("Value") + 
    theme(legend.position="none") +
    theme(strip.background = element_blank()) +
    background_grid()
  g1
}

##- function plot 2: start with list of binded df
##- and limit to boxplot gates
### http://stackoverflow.com/questions/5677885/ignore-outliers-in-ggplot2-boxplot
gg2 <- function(ls , var , lbl, tran = identity){
  
  ##- limit by boundaries of boxplot, by variable of interest
  ls_lim <- lapply(ls, function(x){
    if ( any(grepl('outdegree', names(x))) ) {
      ylim <- tapply(x$outdegree, x[,var], function(y) boxplot.stats(y)$stats[c(1, 5)] )
      for (i in 1:length(ylim)){
        x[(x[, 'outdegree'] < ylim[[i]][1] | x[, 'outdegree'] > ylim[[i]][2]) & x[, var] == i, 'outdegree'] <- NA
      }
      # boxplot.stats(x[x[, var]==i, ]$outdegree)$stats; boxplot(x$outdegree ~ x$risk)
      return(x)
      
    } else {
      
      ##- limits without 1: x$size > 1
      ylim_size <- tapply(x[,]$size, x[, var], function(y) boxplot.stats(y)$stats[c(1, 5)] )
      
      for (i in 1:length(ylim_size)){
        x[(x[, 'size'] < ylim_size[[i]][1] | x[, 'size'] > ylim_size[[i]][2]) & x[, var] == i, 'size'] <- NA
      }
      # hist(x$size); boxplot(x$size ~ x$risk)
      ##- limits without 0
      ylim_nbhsize <-  tapply(x[,]$nbhsize, x[, var], function(y) boxplot.stats(y)$stats[c(1, 5)] )
      #hist(x$nbhsize); boxplot(x$nbhsize ~ x$risk)
      
      for (i in 1:length(ylim_nbhsize)){
        x[(x[, 'nbhsize'] < ylim_nbhsize[[i]][1] | x[, 'nbhsize'] > ylim_nbhsize[[i]][2]) & x[, var] == i, 'nbhsize'] <- NA
      }
      return(x)
    }
  })
  
  ##- melt the limited list
  x <- molten(ls_lim)
  x[,var] <- as.factor(x[,var])
  
  ##- plot
  g1 <- ggplot(x, aes_string(var, 'tran(value)', color = 'variable')) + 
    geom_boxplot() + 
    facet_wrap( ~ variable + L1, scales = "free_y", ncol = 3, labeller= labeller(variable = c(outdegree = "Out-degree", size = "Cluster size", nbhsize = "Neighborhood size")))  + 
    xlab(lbl) + ylab("Value") + 
    theme(legend.position="none") +
    theme(strip.background = element_blank()) +
    background_grid()
  g1
}

##---- bp full ---
gg(tot, 'factor(risk)', 'Risk level')
gg(tot, 'factor(age)', 'Age category')
gg(tot, 'factor(stage)', 'Infection stage')

##--- bp 0.015 ---
gg(tot_2, 'factor(risk)', 'Risk level')
gg(tot_2, 'factor(age)', 'Age category')
gg(tot_2, 'factor(stage)', 'Infection stage')

##---- bp bounded ----
gg2(ls = cw_bind, var = 'risk', lbl = 'Risk level')
gg2(ls = cw_bind, var = 'age', lbl = 'Age category')
gg2(ls = cw_bind, var = 'stage', lbl = 'Infection stage')

