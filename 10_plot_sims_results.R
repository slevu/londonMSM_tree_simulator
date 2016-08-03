###--- First attempt to
###--- Explain simulations with plots of incidence, prevalence, diagnoses, treatment, by age, risk, stage 
###--- with bootstrap uncertainty
###--- and plot tree

source('model0.R')
tr <- sapply( times_day, function(t) diag.t( t, theta ) )
di <- sapply( times_day, function(t) tr.t( t ) )
inc <- sapply( times_day, function(t) inc.t( t , theta ) )

plot( times_year, inc )
plot( times_year, di )
plot( times_year, tr )

t <- 1e4 
incidence <- inc.t( t, theta )

o_a <- ode(y=y0, times=times_day, func=dydt, parms=list()  , method = 'adams')
o_e <- ode(y=y0, times=times_day, func=dydt, parms=list()  , method = 'euler')

str(o_a)
o_a[1:6,1:6]

##- data: global cum incidence
cum_incidence <- rowSums(o_a[,-1])
length(times_year)
plot(times_year, cum_incidence)
##- debug: agg = AGE_COORDS; n = agg[[1]]; desolve = o_a
desolve[, 1 + n]
names

## plots
snazzy.plot <- function( desolve, agg )
{
  oa <- ( sapply( names( agg ), function( n ) {
    rowSums( desolve[, 1 + agg[[n]] ] )  ## NOTE desolve mat has time as first col
  }) )
 # X11() 
  matplot( times_year, oa, type = 'l' )
  legend( x = 'topleft', legend = names(agg)
          , col = 1:length(agg) #SPCOLS[1:length(agg)] 
          , pch = 1
  )
}

snazzy.plot(  o_e, CARE_COORDS ) 
snazzy.plot(  o_a, AGE_COORDS ) 
snazzy.plot(  o_a, RISK_COORDS ) 
snazzy.plot(  o_a, NH_COORDS ) 

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
  # require(phydynR)
  source('model0.R')
  o <- ode(y=y0, times=times_day, func=dydt, parms=list()  , method = 'adams')
  tfgy <- .tfgy( o )
  FF <- tfgy[[2]][[length(times_day)]]
  fmat <- matrix(0, nrow = 4, ncol = 4)
  for (k in 1:120) for (l in 1:120){
    fmat[ deme2age( DEMES[k] ), deme2age(DEMES[l]) ] <- fmat[ deme2age( DEMES[k] ), deme2age(DEMES[l]) ] + FF[k,l]
  }
  
Y <-  tfgy[[4]]

length(Y)

##- from model0-simulate-Baseline0
#~ incidence and prevalence 
yfin <- tfgy[[4]][[length(times_day)]]
ffin <- tfgy[[2]][[length(times_day)]]
newinf <- sum(ffin[1:120, 1:120] ) * 365 
plwhiv <- sum( yfin[-length(yfin)] )
