# rm(list=ls())

##---- libs ----
require(phydynR)
library(reshape2)
library(ggplot2)
library(cowplot)
library(scales)

##---- source ----
source('model1.R') # contains parameter values for sims
source('Erik_code/stage_rate_estimator1.R')

##---- cohort effect control ----
THRESHOLD_YRS <- 2

##---- od by stage helper ----
## fn <- list.sims[['Baseline0']][1]
fns2od.by.stage <- function( fns )
{
  o <- list()
  for (fn in fns){
    load(fn) # daytree, bdt, W, cd4s, sampleDemes, plwhiv, newinf, MH
    samplesInCohort <- names(bdt$sampleTimes[bdt$sampleTimes > (max(bdt$sampleTimes)-THRESHOLD_YRS*365) ] )
    samplesInCohort_donors <- intersect( samplesInCohort, W$donor)
    
    i <- which( (W$donor %in% samplesInCohort) & (W$recip %in% samplesInCohort ) )
    WW <- list( donor = W$donor[i] , recip = W$recip[i], infectorProbability = W$infectorProbability[i] )
    
    system.time({
      od <- sapply( samplesInCohort , function(u){
        if (!(u%in% samplesInCohort_donors)) return(0)
        i <- which(WW$donor==u)
        sum(WW$infectorProbability[i] ) 
      })
    })
    
    rownames(bdt$sampleStates) <- bdt$tip.label
    stages <- setNames( sapply( 1:nrow(bdt$sampleStates), function(k){
      deme <- DEMES[ which.max(bdt$sampleStates[k,] )  ]
      stage <- strsplit( deme, '.' , fixed=T)[[1]][1]
      as.numeric( tail( strsplit(stage, '')[[1]], 1 ) )
    }), bdt$tip.label)
    
    stages <- stages[ names(stages) %in% samplesInCohort ]
    
    od_by_stage <- lapply( 1:5, function(stage){
      od[ names(stages)[which(stages==stage)] ]
    })
    names(od_by_stage) <- 1:5

    o[[fn]] <- od_by_stage
    # print(date())
  }
  o
}

##---- paths ----
if( any(grep("MacBook", Sys.info())) ){
  path.sims <- '../Box Sync/HPC/simulations/model1-sim' #'data/simulations2/model0-simulate'
  path.results <- '../Box Sync/HPC/simulations/model1-sim_ucsd'
} else {
  path.sims <- '../Box/HPC/simulations/model1-sim'
  path.results <- '../Box/HPC/simulations/model1-sim_ucsd' # imac
}

##- scenario
scenario <- c("Baseline0", "EqualStage0")
scenario <- setNames(scenario, scenario) # useful to name list in lapply

##- list of sims files and distances files
list.sims <- lapply(scenario, function(x){
  list.files('RData', full.names = TRUE, 
             path = paste(path.sims, x, sep = '') )
})
# str(list.sims)

##---- od by stage ----
## list of od vectors by stage, named by tips
ODSTAGE <- paste(path.results, 'od.by.stage.RData', sep = '/')
if(FALSE){
system.time(
obs_bl <- fns2od.by.stage( list.sims[["Baseline0"]] )
) # 106
system.time(
obs_er <- fns2od.by.stage( list.sims[["EqualStage0"]] )
) # 94

save(obs_bl, obs_er, file = ODSTAGE )
} else {
  load(ODSTAGE)
}

# str(obs_er[1])
# names(obs_bl)

##---- plot od by stage ----
##-- long dataframes
od_stage_bl <- melt(unname(obs_bl) )
od_stage_er <- melt(unname(obs_er) )
# head(od_stage_bl); tail(od_stage_er)
##-- one long table
od_stage <- rbind(
  cbind(od_stage_bl[,-3], scenario = 'BA'),
  cbind(od_stage_er[,-3], scenario = 'ER')
                  )
head(od_stage)

od_stage_boxplot <- function(df){
  p <- ggplot(df, aes(factor(L2), value, colour = scenario), outlier.colour = alpha(colours(), 0.001) ) + geom_boxplot() + # facet_wrap(~ scenario) + 
    theme(legend.position="top") + 
    background_grid() +
    theme(strip.background = element_blank()) + 
    xlab("Infection stage") + ylab("Out-degrees")
  p
}

bp <- od_stage_boxplot(od_stage)
bp

##- better and simpler without outliers ?
boxplot(value ~ L2*scenario, data = od_stage, outline = FALSE)

##---- U test ----
## U test for the equal rates simulations: 
wtest_er <- lapply( obs_er, function(obs) {
  wilcox.test( obs[[1]] / pstage[1], 
  # NOTE pstage[1] is propto average duration of EHI
  obs[[5]] / (1-pstage[1] ), 
  # c(obs[[2]],obs[[3]],obs[[4]], obs[[5]]) / (1-pstage[1] ),  
  # NOTE 1-pstage[1] is propto average duration of the rest of the infectious period
  alternative = 'greater' 
  # NOTE this is a one-tailed test. H1: EHI rate > !EHI rate
  )
})

## U test for the baseline simulations: 
wtest_bl <- lapply( obs_bl, function(obs) {
  wilcox.test( obs[[1]] / pstage[1] 
               , obs[[5]] / (1-pstage[1] ) 
               , alternative = 'greater'
  )
})

type1error <- mean( sapply( wtest_er, function(wt) wt$p.value) < .05 ) 
type2error <- mean( sapply(wtest_bl, function(wt) wt$p.value) > .05 ) 

type1error
type2error
# est.rd.batch
# rate.difference.ehi
# pstage

##---- apply rate difference ----
#- echo = FALSE
rd_er <- est.rd.batch( obs_er )
rd_br <- est.rd.batch( obs_bl )
# str(rd_br)

##---- plot errors SA ----
##- EHI/Late rate difference (Equal transmission rates)
boxplot( unname(rd_er), main = '' , xaxt="n")
axis(1, at = c(1,1:10*10))
title(ylab = 'Rate difference: EHI - non EHI', xlab = 'Simulation replicate', cex.lab = 1.2)
abline( h = 0, col = 'red' )

##- EHI/Late rate difference (Baseline)
boxplot( unname( rd_br ), main = '', xaxt="n")
axis(1, at = c(1,1:10*10))
title(ylab = 'Rate difference: EHI - non EHI', xlab = 'Simulation replicate', cex.lab = 1.2)
abline( h = 0, col = 'red' )

##---- apply rate ratio ----
rate_ratio_ehi
str(obs_bl) # nsim lists of 5 lists of od by stage 
od_by_stage <- obs_bl[[1]]
rrs <- rate_ratio_ehi(obs_bl[[1]] , nreps=1e3 )
## get rr = 0.5 and rd = -.35 ### problem, debug the functions
####-------------------------------------------------- lllllllaaaaaaa 29/6/17 -----!!!!!!

###################--- CLUSTERS ---###
##---- clusters ----
## questions:
## - how to restrict cluster sizes by THRESHOLD_YRS ?
## - explain why EHI rate is compared with stage 5 rate ?

##---- load cluster assignements
l_Baseline0 <- readRDS(file = paste(path.results, "list.sim.ucsd.Baseline0.rds", sep = '/') )
l_EqualStage0 <- readRDS(file = paste(path.results, "list.sim.ucsd.EqualStage0.rds", sep = '/') )
# names(l_Baseline0); str(l_Baseline0[[5]][[1]])

##---- pruning function ----
##- Pruning: cluster size and membership are recomputed as time restriction exclude patients
##- function to prune cluster according to a threshold of sampling time to control for cohort effect
##- recalculate size and cluster membership

prune.clus <- function(a, t , st){
  ## subset df by sampling times
  cohort <- names(st[st > (max(st) - t ) ] ) 
  aa <- a[a[,"id"] %in% cohort,]
  if(identical(aa,a)){
    print('do nothing')
  } else {
    ## for each clusterID, re-calculate size and binclus membership
    for (i in unique(aa[,"ClusterID"]) ){
      aa[ aa[,"ClusterID"] == i, "size" ] <- nrow(aa[ aa[,"ClusterID"] == i, ]) 
    }
    aa[,"binclus"] <- ifelse(aa[,"size"] < 2, 0, 1 ) 
  }
  return(aa) 
}

##---- pruning ---
##- sampling times (same across all sims)
load(list.sims[["Baseline0"]][[1]]) # for bdt
st_yrs <- days2years(bdt$sampleTimes)

system.time(
l_Baseline0.pruned <- lapply(l_Baseline0, function(x){
  lapply(x, function(df){
    prune.clus(a = df, t = THRESHOLD_YRS, st = st_yrs )
  })
})
) # 115
system.time(
  l_EqualStage0.pruned <- lapply(l_EqualStage0, function(x){
    lapply(x, function(df){
      prune.clus(a = df, t = THRESHOLD_YRS, st = st_yrs )
    })
  })
) # 117


##---- u test cluster rate ----
u.test.stage <- function(df){
  sizes_stage1 <-  df[df$stage == 1, ]$size 
  # sizes_otherstages <-  df[df$stage != 1, ]$size
  sizes_aids <-  df[df$stage == 5, ]$size
#   U <- wilcox.test(sizes_stage1 ,
#               sizes_otherstages ,
#               alternative = 'greater')
  U <- wilcox.test(sizes_stage1 / pstage[1],
                   sizes_aids / (1 - pstage[1]),
                   alternative = 'greater')
  return(U$p.value)
}

## test
# df = l_Baseline0.pruned[["0.015"]][[1]]
# summary(sizes_stage1)
# summary(sizes_aids)
# str(df)

p.cl_bl <- sapply(l_Baseline0.pruned[["0.015"]], u.test.stage)
p.cl_er <- sapply(l_EqualStage0.pruned[["0.015"]], u.test.stage)

##- type1err: 100% would find larger cluster size with EHI
mean(p.cl_er < 0.05)
##- type2err: 0% would not find difference when there is
mean(p.cl_bl > 0.05)

qqplot(p.cl_bl, p.cl_er)

##---- size by stage ----
# x <- tapply(df$size, df$stage, identity)

ll <- l_Baseline0[['0.015']]
names(ll)

xy <- lapply(ll, function(x){
  ( tapply(x$size, x$stage, identity) )
})

str(xy[1])

# xx <- melt( unlist(xy, recursive = F, use.names = F ) )
xx <- melt( xy )
str(xx)
tail(xx)
?unlist
?attributes
attributes(xy) <- NULL
