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

##--- plots ---
##- long table with value = od, size ...; method = SA, CL, NB; thr = ..
##- for plots
source('functions.R')
molten <- function(cw_bind){
  .sa <- melt( cw_bind[[1]], id.vars = c('id', 'stage', 'age', 'risk'),  measure.vars = 'outdegree')
  # head(.sa)
  .thr <- melt(cw_bind[-1], id.vars = c('id', 'stage', 'age', 'risk'), measure.vars = c('size', 'nbhsize'))
  # head(.thr)
  tot <- rbind( cbind(.sa, "L1" = 'NA'), .thr )
  # head(tot)
  return(tot)
}

##- winsorized (with function #1)
cw_bind_win1 <- lapply(cw_bind, function(x){
  if ( any(grepl('outdegree', names(x))) ) {
    x[, 'outdegree'] <- winsorize( x[, 'outdegree'], 0.02 )
    return(x)
  } else {
    x[, 'size'] <- winsorize( x[, 'size'], 0.02 )
    x[, 'nbhsize'] <- winsorize( x[, 'nbhsize'], 0.02 )
    return(x)
  }
})

##- winsorized (with function #2)
cw_bind_win2 <- lapply(cw_bind, function(x){
  if ( any(grepl('outdegree', names(x))) ) {
    x[, 'outdegree'] <- winsorize2( x[, 'outdegree'])
    return(x)
  } else {
    x[, 'size'] <- winsorize2( x[, 'size'])
    x[, 'nbhsize'] <- winsorize2( x[, 'nbhsize'])
    return(x)
  }
})

##- no minimum value
cw_bind_nozero <- lapply(cw_bind, function(x){
  if ( any(grepl('outdegree', names(x))) ) {
    x[x[, 'outdegree'] == 0, 'outdegree'] <- NA
    return(x)
  } else {
    x[x[, 'size'] == 1 , 'size'] <- NA
    x[x[, 'nbhsize'] == 0, 'nbhsize'] <- NA
    return(x)
  }
})

# str(cw_bind_nozero)
summary(cw_bind[[2]]$nbhsize)
summary(cw_bind_win[[2]]$nbhsize)

tot <- molten(cw_bind)
tot_win <- molten(cw_bind_win)
tot_win2 <- molten(cw_bind_win2)
tot_nozero <- molten(cw_bind_nozero)

# str(tot_win)

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

##- apply
gg(tot, 'factor(risk)', 'Risk level')
gg(tot_win, 'factor(risk)', 'Risk level')
gg(tot_win2, 'factor(risk)', 'Risk level')
gg(tot_nozero, 'factor(risk)', 'Risk level')

gg(tot, 'factor(age)', 'Age category')
gg(tot_win, 'factor(age)', 'Age category')
gg(tot_nozero, 'factor(age)', 'Age category')

gg(tot, 'factor(stage)', 'Infection stage')
gg(tot_win, 'factor(stage)', 'Infection stage')
gg(tot_nozero, 'factor(stage)', 'Infection stage')

###
# http://stackoverflow.com/questions/5677885/ignore-outliers-in-ggplot2-boxplot
# create a dummy data frame with outliers
df = data.frame(y = c(-100, rnorm(100), 100))

# create boxplot that includes outliers
p0 = ggplot(df, aes(y = y)) + geom_boxplot(aes(x = factor(1)))


# compute lower and upper whiskers
ylim1 = boxplot.stats(df$y)$stats[c(1, 5)]

# scale y limits based on ylim1
p1 = p0 + coord_cartesian(ylim = ylim1*1.05)
