####---- compute true values from model0 ----
# rm(list=ls())

####---- helper functions ----
deme2age <- function(deme){ as.numeric(
  substr(regmatches( deme , regexpr( '\\.age[0-9]', deme ) ), 5, 5) 
) }
deme2stage <- function(deme){as.numeric( 
  substr(regmatches( deme , regexpr( 'stage[0-9]', deme )), 6, 6)
) }
deme2risk <- function(deme){as.numeric( 
  substr(regmatches( deme , regexpr( 'riskLevel[0-9]', deme )), 10, 10)
) }

##- run model
if (TRUE){
  require(phydynR)
  source('model0.R')
  o <- ode(y=y0, times=times_day, func=dydt, parms=list()  , method = 'adams')
  tfgy <- .tfgy( o )
  FF <- tfgy[[2]][[length(times_day)]]
  fmat <- matrix(0, nrow = 4, ncol = 4)
  for (k in 1:120) for (l in 1:120){
    fmat[ deme2age( DEMES[k] ), deme2age(DEMES[l]) ] <- fmat[ deme2age( DEMES[k] ), deme2age(DEMES[l]) ] + FF[k,l]
  }
}