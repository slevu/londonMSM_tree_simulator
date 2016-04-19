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