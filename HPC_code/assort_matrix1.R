###--- to run on HPC ---
args <- commandArgs(TRUE)
sim <- args[1] # sim <- "../Box/HPC/simulations/model1-simBaseline0/9131.RData"
pwd <- args[2]
print(pwd)

##- RData now also includes distance matrix
print(sim)
name.sim <- paste(regmatches(sim, regexpr("[0-9]{3,}", sim)), ".RData", sep = '')
FN <- paste('data/simulations2/model1_age/age-anal_', name.sim, sep = '') # filename results
print(FN)

if(TRUE){
  ####---- helper functions ----
  deme2age <- function(deme){ as.numeric(
    substr(regmatches( deme , regexpr( '\\.age[0-9]', deme ) ), 5, 5) 
  ) }
  deme2stage <- function(deme){as.numeric( 
    substr(regmatches( deme , regexpr( 'stage[0-9]', deme )), 6, 6)
  ) }
  
  
  THRESHOLD_YRS <- Inf #10 # only count donors sampled within this many years of the most recent sample (control for cohort effects) 
  # NOTE this doesnt appear to be important for age, whereas it was important for stage, so setting to Inf
  
  ##- function
  ## sim <- "data/simulations2/model1-simBaseline0/9131.RData"
  age_mat <- function(sim, thr)
  {
    load(sim)
    samplesInCohort <- names(bdt$sampleTimes[bdt$sampleTimes > (max(bdt$sampleTimes) - THRESHOLD_YRS*365) ] )
    samplesInCohort_donors <- intersect( samplesInCohort, W$donor)
    
    i <- which( (W$donor %in% samplesInCohort) & (W$recip %in% samplesInCohort ) )
    WW <- list( donor = W$donor[i] , recip = W$recip[i], infectorProbability = W$infectorProbability[i] )
    
    ##- outdegree with 0 if not donor
    # system.time({ 
    od <- sapply( samplesInCohort , function(u){
      if (!(u%in% samplesInCohort_donors)) return(0)
      i <- which(WW$donor==u)
      sum(WW$infectorProbability[i] ) 
    })
    # }) # 39
    
    age <- sapply( sampleDemes[ samplesInCohort ], deme2age )
    stage <- sapply( sampleDemes[ samplesInCohort ], deme2stage )
    
    amat <- matrix( 0, nrow=4, ncol=4)
    #system.time(
    for (k in 1:length(WW$donor)){
      u <- WW$donor[k]
      v <- WW$recip[k]
      au <- deme2age( sampleDemes[ u] )
      av <- deme2age( sampleDemes[v] )
      amat[ au, av ] <- amat[ au, av ] + WW$infectorProbability[k]
    }
    #) # 97s
    rownames(amat) = colnames(amat) <- 1:4
    
    ##- load distance mat
    #load(dfn) #D
    D <- distance[['D']]
    DD <- D
    DD <- DD[ (D[,1] %in% samplesInCohort) & (D[,2] %in% samplesInCohort), ]
    
    ##- age matrix based on neighborhood|threshold
    l <- list()
    # system.time(
    for (i in 1:length(thr)){
      CL_THRESHOLD <- as.numeric(thr[i])
      
      nbrhoodSize <- setNames( rep(0, length(samplesInCohort)), samplesInCohort )
      agemat2 <- matrix( 0, nrow = 4, ncol = 4 )
      
      ##- age matrix for neighborhood
      for (k in 1:nrow(DD)){
        if (DD[k,3] < CL_THRESHOLD ) {
          u <- as.character(DD[k,1])
          v <- as.character(DD[k,2])
          au <- deme2age( sampleDemes[u] )
          av <- deme2age( sampleDemes[v])
          agemat2[au, av] <- agemat2[ au,av] + 1
          agemat2[av, au] <- agemat2[ av,au] + 1
          nbrhoodSize[u] <- nbrhoodSize[u] + 1
          nbrhoodSize[v] <- nbrhoodSize[v] + 1
        }
      }
      l[[i]] <- list("agemat2" = agemat2, 
                     "nbrhoodSize" = nbrhoodSize)
      
      names(l)[i] <- thresholds[i]
      print(paste(name.sim, thresholds[i]))
    }
    # ) # 164s
    
    ll <- list("agemat_sa" = amat, 
               "tab_sa" = data.frame(od=od, age=age, stage=stage), 
               "agemat_cl" = l)
    
    return(ll)
  }
  
  ##---- apply function ----
  ## by thresholds
  thresholds <- c("0.005", "0.015", "0.02",  "0.05")
  
  system.time(
    ll <- age_mat(sim, thr = thresholds)
  )
  save(ll, file = FN)
  # str(ll, vec.len = 1, indent.str = "|", comp.str = "----")
  
}
