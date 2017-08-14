# what is the true rate ratio / rate difference by risk ?
# a fixed proportion of patients are at high risk level, conferring a 10 time increase in transmission
# scenario with equal stage rates will be easier to look at
# 
# from input parameter OR 
# from transition matrix from the last year ?

source("HPC_code/model1.R")
## calculate transmission rate ratio and rate difference by risk level in diagnosed as recent
o <- ode(y=y0, times=times_day, func=dydt, parms=list()  , method = 'rk4')
tfgy <- .tfgy( o )

## util functions
deme2risk <- function(deme){as.numeric( 
  substr(regmatches( deme , regexpr( 'riskLevel[0-9]', deme )), 10, 10)
) }
deme2stage <- function(deme){as.numeric( 
  substr(regmatches( deme , regexpr( 'stage[0-9]', deme )), 6, 6)
) }

##
FF <- tfgy[[2]][[length(times_day)]]
fmat <- matrix(0, nrow = 2, ncol = 2) # matrix(0, nrow = 4, ncol = 4)
for (k in 1:120) for (l in 1:120){
  fmat[ deme2risk( DEMES[k] ), deme2risk(DEMES[l]) ] <- fmat[ deme2risk( DEMES[k] ), deme2risk(DEMES[l]) ] + FF[k,l]
}
fmat
#          [,1]      [,2]
# [1,] 1.203422 0.3008554
# [2,] 3.008555 0.7521388

## with stage
## lllllllllaaaaaaaaaaa
deme2stagerisk <- expand.grid(paste0('s',1:5), paste0('r',1:2),  stringsAsFactors = FALSE)
fmat_s <- matrix(0, nrow = 8, ncol = 8) # matrix(0, nrow = 4, ncol = 4)
for (k in 1:120) for (l in 1:120){
  fmat_s[ deme2risk( DEMES[k] ), deme2risk(DEMES[l]) ] <- fmat_s[ deme2risk( DEMES[k] ), deme2risk(DEMES[l]) ] + FF[k,l]
}
fmat
## 80% risk 1 with risk_wtransm = 1
## 20% risk 2 with risk_wtransm = 10