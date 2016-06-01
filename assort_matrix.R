# rm(list=ls())

####---- path sims ----
##- used for several rounds of simulations
path.sims <- 'data/simulations2/model0-simulate'

####---- scenario ----
scenario <- c("EqualStage0", "Baseline0")

####--- list.trees ----
list.sims <- vector("list", length(scenario))
for (.s in 1:length(scenario)){
  list.sims[[.s]] <- list.files('RData', full.names = TRUE, 
                                path = paste(path.sims, 
                                             scenario[.s], sep = '') )
  ## names of sims within scenario
  names(list.sims[[.s]]) <- lapply(list.sims[[.s]], function(x){
    regmatches(x, regexpr("[0-9]{2,}", x)) # numerical string of length >= 2
  })
  ## names of scenario
  names(list.sims)[.s] <- scenario[.s]
}

####---- list.dist ----
## load list of dist from sims
  distEqualStage0FNS <- list.files('RData', full.names=T, path = paste(path.sims, 'EqualStage0-distances', sep = '') )
  distBaseline0FNS <- list.files('RData', full.names=T, path = paste(path.sims, 'Baseline0-distances', sep = '') )

####---- helper functions ----

    deme2age <- function(deme){ as.numeric(
      substr(regmatches( deme , regexpr( '\\.age[0-9]', deme ) ), 5, 5) 
    ) }
    deme2stage <- function(deme){as.numeric( 
      substr(regmatches( deme , regexpr( 'stage[0-9]', deme )), 6, 6)
    ) }
    
  
THRESHOLD_YRS <- Inf #10 # only count donors sampled within this many years of the most recent sample (control for cohort effects) 
# NOTE this doesnt appear to be important for age, whereas it was important for stage, so setting to Inf

# Clustering threshold:
# CL_THRESHOLD <- .015

##- function
age_mat <- function(sim, dist, thr)
{
  name.sim <- regmatches(sim, regexpr("[0-9]{3,}", sim))
  dfn <- dist[ grep( name.sim, dist) ]
  
  load(sim)
  samplesInCohort <- names(daytree$sampleTimes[daytree$sampleTimes > (max(daytree$sampleTimes) - THRESHOLD_YRS*365) ] )
  samplesInCohort_donors <- intersect( samplesInCohort, W$donor)
  
  i <- which( (W$donor %in% samplesInCohort) & (W$recip %in% samplesInCohort ) )
  WW <- list( donor = W$donor[i] , recip = W$recip[i], infectorProbability = W$infectorProbability[i] )
  
  # system.time({
  od <- sapply( samplesInCohort , function(u){
    if (!(u%in% samplesInCohort_donors)) return(0)
    i <- which(WW$donor==u)
    sum(WW$infectorProbability[i] ) 
  })
 # }) # 39
  
  age <- sapply( sampleDemes[ samplesInCohort ], deme2age )
  stage <- sapply( sampleDemes[ samplesInCohort ], deme2stage )
  
  amat = agemat <- matrix( 0, nrow=4, ncol=4)
  for (k in 1:length(WW$donor)){
    u <- WW$donor[k]
    v <- WW$recip[k]
    au <- deme2age( sampleDemes[ u] )
    av <- deme2age( sampleDemes[v] )
    amat[ au, av ] <- amat[ au, av ] + WW$infectorProbability[k]
  }
  rownames(amat) = colnames(amat) <- 1:4
  
  #load distance mat
  load(dfn) #D
  DD <- D
  D1 <- daytree$tip.label[ D[1,] ]
  D2 <- daytree$tip.label[ D[2,] ]
  DD <- DD[ , (D1 %in% samplesInCohort) & (D2 %in% samplesInCohort) ]
  
  
  l <- list()
      for (i in 1:length(thr)){
      CL_THRESHOLD <- thr[i]
     
      nbrhoodSize <- setNames( rep(0, length(samplesInCohort)), samplesInCohort )
      agemat2 <- matrix( 0, nrow = 4, ncol = 4 )
      
    ##- age matrix for neighborhood
      for (k in 1:ncol(DD)){
        if (DD[3,k] < CL_THRESHOLD ) {
          u <- daytree$tip.label[ DD[1,k] ]
          v <- daytree$tip.label[ DD[2,k] ]
          au <- deme2age( sampleDemes[u] )
          av <- deme2age( sampleDemes[v])
          agemat2[au, av] <- agemat2[ au,av] + 1
          nbrhoodSize[u] <- nbrhoodSize[u] + 1
        }
      }
      l[[i]] <- list(agemat2 = agemat2, 
                     nbrhoodSize = nbrhoodSize)
      
      names(l)[i] <- thresholds[i]
      print(paste(name.sim, thresholds[i]))
    }
    
  ll <- list(agemat_sa = amat, tab_sa = data.frame(od=od, age=age, stage=stage),  agemat_cl = l)
  
  save(ll, file = paste('data/simulations2/age/age-anal_', name.sim, '.RData', sep = ''))
  return(ll)
}

##---- apply function ----
## by thresholds
thresholds <- c("0.001", "0.005", "0.015", "0.05")

system.time(
test <- lapply( list.sims[["Baseline0"]][1:2], function(sim){
  age_mat(sim, dist = distBaseline0FNS, thr = thresholds)
  })
) # 381s for 2 sims
str(test)
str(fns2agemat_tables)
tail(fns2agemat_tables[[2]][[3]])

mat2assortmat <- function(mat){
  rs <- rowSums(mat)
  cs <- colSums(mat)
  s <- sum(mat)
  M <- matrix(0, nrow=nrow(mat), ncol = ncol(mat))
  A <- matrix(0, nrow=nrow(mat), ncol = ncol(mat))
  k <- nrow(mat)
  for (i in 1:k) for (j in 1:k){
    M[i,j] <- rs[i] * cs[j] / s
    A[i,j] <- ( mat[i,j] - M[i,j] )/ M[i,j]
  }
  A
}

assor_mat <- mat2assortmat(fns2agemat_tables[[1]][["agemat"]])
assor_mat2 <- mat2assortmat(fns2agemat_tables[[1]][["agemat2"]])
sum(assor_mat)
sum(assor_mat2)
require(lattice)
# lattice.options(default.theme = standard.theme(color=F))
levelplot( assor_mat, col.regions = heat.colors)
levelplot( assor_mat2, col.regions = heat.colors)

