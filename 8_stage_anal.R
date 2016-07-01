# rm(list=ls())

##---- libs ----
require(phydynR)

##---- source ----
source('model0.R') # contains parameter values for sims
source('Erik_code/stage_rate_estimator1.R')

##---- cohort effect control ----
THRESHOLD_YRS <- 2

##---- od by stage helper ----
fns2od.by.stage <- function( fns )
{
  o <- list()
  for (fn in fns){
    load(fn) # daytree, bdt, W, cd4s, sampleDemes, plwhiv, newinf, MH
    samplesInCohort <- names(daytree$sampleTimes[daytree$sampleTimes > (max(daytree$sampleTimes)-THRESHOLD_YRS*365) ] )
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
    
    rownames(daytree$sampleStates) <- daytree$tip.label
    stages <- setNames( sapply( 1:nrow(daytree$sampleStates), function(k){
      deme <- DEMES[ which.max(daytree$sampleStates[k,] )  ]
      stage <- strsplit( deme, '.' , fixed=T)[[1]][1]
      as.numeric( tail( strsplit(stage, '')[[1]], 1 ) )
    }), daytree$tip.label)
    
    stages <- stages[ names(stages) %in% samplesInCohort ]
    
    od_by_stage <- lapply( 1:5, function(stage){
      od[ names(stages)[which(stages==stage)] ]
    })
    
    o[[fn]] <- od_by_stage
    # print(date())
  }
  o
}

####---- path sims ----
##- used for several rounds of simulations
path.sims <- 'data/simulations2/model0-simulate'
path.results <- 'data/sim_ucsd_results2'

####---- scenario ----
scenario <- c("Baseline0", "EqualStage0")
scenario <- setNames(scenario, scenario) # useful to name list in lapply

####--- list of sims files and distances files ---
list.sims <- lapply(scenario, function(x){
  list.files('RData', full.names = TRUE, 
             path = paste(path.sims, x, sep = '') )
})
# str(list.sims)

###---- od by stage ----
## list of od vectors by stage, named by tips
if(FALSE){
system.time(
obs_bl <- fns2od.by.stage( list.sims[["Baseline0"]] )
)
system.time(
obs_er <- fns2od.by.stage( list.sims[["EqualStage0"]] )
) # 84

# save(obs_bl, obs_er, file = 'data/od.by.stage.RData' )
} else {
  load('data/od.by.stage.RData')
}

# str(obs_er[[1]])

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

##---- plot errors SA ----
#- echo = FALSE
rd_er <- est.rd.batch( obs_er )
rd_br <- est.rd.batch( obs_bl )

##- EHI/Late rate difference (Equal transmission rates)
boxplot( unname(rd_er), main = '' )
abline( h = 0, col = 'red' )

##- EHI/Late rate difference (Baseline)
boxplot( unname( rd_br ), main = '')
abline( h = 0, col = 'red' )

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
st_yrs <- days2years(daytree$sampleTimes)

system.time(
l_Baseline0.pruned <- lapply(l_Baseline0, function(x){
  lapply(x, function(df){
    prune.clus(a = df, t = THRESHOLD_YRS, st = st_yrs )
  })
})
)
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
  sizes_otherstages <-  df[df$stage != 1, ]$size
#   U <- wilcox.test(sizes_stage1 ,
#               sizes_otherstages ,
#               alternative = 'greater')
  U <- wilcox.test(sizes_stage1 / pstage[1],
                   sizes_otherstages / (1 - pstage[1]),
                   alternative = 'greater')
  return(U$p.value)
}

## test
# df = l_Baseline0.pruned[["0.015"]][[1]]
# summary(sizes_stage1)
# summary(sizes_otherstages)

p.cl_bl <- sapply(l_Baseline0.pruned[["0.015"]], u.test.stage)
p.cl_er <- sapply(l_EqualStage0.pruned[["0.015"]], u.test.stage)

##- type1err: 100% would find larger cluster size with EHI
mean(p.cl_er < 0.05)
##- type2err: 0% would not find difference when there is
mean(p.cl_bl > 0.05)

qqplot(p.cl_bl, p.cl_er)


