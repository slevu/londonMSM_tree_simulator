require(phydynR)
source('model0.R') # contains parameter values for sims
source('Erik_code/stage_rate_estimator1.R')


THRESHOLD_YRS <- 2

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
    print(date())
  }
  o
}

####

####---- scenario ----
scenario <- c("EqualStage0", "Baseline0")

####--- list.trees ----
list.sims <- vector("list", length(scenario))
for (s in 1:length(scenario)){
  list.sims[[s]] <- list.files('RData', full.names = TRUE, 
                               path = paste('data/simulations/model0-simulate', 
                                            scenario[s], sep = '') )
  names(list.sims)[s] <- scenario[s]
}
str(list.sims)

###-
obs_bl <- fns2od.by.stage( list.sims[["Baseline0"]] )
obs_er <- fns2od.by.stage( list.sims[["EqualStage0"]] )
# save(obs_bl, obs_er, file = 'data/od.by.stage.RData' )
# load('data/od.by.stage.RData')
str(obs_er[[1]])

# U test for the equal rates simulations: 
wtest_er <- lapply( obs_er, function(obs) {
  wilcox.test( obs[[1]] / pstage[1] # NOTE pstage[1] is propto average duration of EHI
               #, obs[[5]] / (1-pstage[1] ) 
               , c(obs[[2]],obs[[3]],obs[[4]], obs[[5]]) / (1-pstage[1] )  #NOTE 1-pstage[1] is propto average duration of the rest of the infectious period
               , alternative = 'greater' #NOTE this is a one-tailed test. H1: EHI rate > !EHI rate
  )
})

# U test for the baseline simulations: 
wtest_bl <- lapply( obs_bl, function(obs) {
  wilcox.test( obs[[1]] / pstage[1] 
               , obs[[5]] / (1-pstage[1] ) 
               , alternative = 'greater'
  )
})

type1error <- mean( sapply( wtest_er, function(wt) wt$p.value) < .05 ) 
type2error <- mean( sapply(wtest_bl, function(wt) wt$p.value) > .05 ) 

est.rd.batch
rate.difference.ehi
pstage

###- do the same thing with clusters
l_Baseline0 <- readRDS(file = "data/sim_ucsd_results/list.sim.ucsd.Baseline0.rds" )
l_EqualStage0 <- readRDS(file = "data/sim_ucsd_results/list.sim.ucsd.EqualStage0.rds" )

u.test.stage <- function(df){
  sizes_stage1 <-  df[df$stage == 1, ]$size 
  sizes_otherstages <-  df[df$stage != 1, ]$size
  U <- wilcox.test(sizes_stage1 ,
              sizes_otherstages ,
              alternative = 'greater')
#   U <- wilcox.test(sizes_stage1 / pstage[1],
#                    sizes_otherstages / (1 - pstage[1]),
#                    alternative = 'greater')
  return(U$p.value)
}

## test
# df = l_Baseline0[["0.015"]][[1]]
# summary(sizes_stage1)
# summary(sizes_otherstages)

p.cl_bl <- sapply(l_Baseline0[["0.015"]], u.test.stage)
p.cl_er <- sapply(l_EqualStage0[["0.015"]], u.test.stage)

##- type1err
mean(p.cl_er < 0.05)
##- type2err
mean(p.cl_bl > 0.05)


###-
rd_er <- est.rd.batch( obs_er )
rd_br <- est.rd.batch( obs_bl )

##- EHI/Late rate difference (Equal transmission rates)
boxplot( unname(rd_er), main = '' )
abline( h = 0, col = 'red' )

##- EHI/Late rate difference (Baseline)
boxplot( unname( rd_br ), main = '')
abline( h = 0, col = 'red' )