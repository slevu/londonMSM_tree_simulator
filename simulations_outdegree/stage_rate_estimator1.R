# derive proportion of time spend in stage | not diagnosed:
{
	progparms<- c(	gamma0 = 1/(6/12) # EHI
	  , gamma1 = 0.157 
	  , gamma2 = 0.350
	  , gamma3 = 0.282
	  , gamma4 = .434 #aids 
	  , pstartstage1 = 0.58
	  , pstartstage2 = 0.23  #reg4 combine these 
	  , pstartstage3 = 0.16
	  , pstartstage4 = 0.03
	)
	pstartstage <- progparms[ c('pstartstage1','pstartstage2', 'pstartstage3', 'pstartstage4') ]
	DEMENAMES <- paste(sep='', 'stage', 0:4)
	
	# prior probability of each stage 
	# NOTE this would only be approx, since does not account for stage 0 going to stages > 1: 
	#pstage <- (1/progparms[c('gamma0', 'gamma1', 'gamma2', 'gamma3', 'gamma4')] ) / sum(1/progparms[c('gamma0', 'gamma1', 'gamma2', 'gamma3', 'gamma4')])
	prEverReach <- cumsum(  pstartstage )
	meanDurNRStages <- (1/progparms[c( 'gamma1', 'gamma2', 'gamma3', 'gamma4')] )
	pstage1 <- c( 1/progparms['gamma0'], prEverReach * meanDurNRStages)
	pstage <- pstage1 / sum(pstage1)
}



rate_ratio_ehi <- function(od_by_stage , nreps=1e3 )
{
#~ 	ir_stage <- isoreg( stages[names(od)] , od )
	mod_resamples <- lapply( 1:5, function(stage){
		n <- length(od_by_stage[[stage]])
		rnorm( nreps, mean(od_by_stage[[stage]]) , sd = sd(od_by_stage[[stage]]) / sqrt(n-1))
	})
	
	# rate ratio: 
	rr <- as.vector(sapply( 1:nreps, function(irep){
		mod_ehi <- mod_resamples[[1]][irep] # mean out degree ehi
		mod_non_ehi <- mod_resamples[[5]][irep] # mean out degree non ehi 
		rate_non_ehi <-  max(0,(mod_non_ehi - mod_ehi)) / (1 - pstage[1] ) 
		rate_ehi <- ( mod_ehi / pstage[1] )
		rate_ehi / rate_non_ehi
	}))
	rr
}

rate.difference.ehi <- function(od_by_stage , nreps = 1e3 )
{
	# od_by_stage list of vectors with out degree for each stage
	# use asymptotic distribution of sample mean to derive CI of difference between ehi and !ehi
	mod_resamples <- lapply( 1:5, function(stage){
		n <- length(od_by_stage[[stage]])
		rnorm( nreps, mean(od_by_stage[[stage]]) , sd = sd(od_by_stage[[stage]]) / sqrt(n-1))
	})
	
	rd <- as.vector(sapply( 1:nreps, function(irep){
		mod_ehi <- mod_resamples[[1]][irep] # mean out degree ehi
		mod_non_ehi <- mod_resamples[[5]][irep] # mean out degree non ehi 
		rate_non_ehi <-  max(0,(mod_non_ehi - mod_ehi)) / (1 - pstage[1] ) 
		rate_ehi <- ( mod_ehi / pstage[1] )
		rate_ehi - rate_non_ehi
	}))
	rd	
}

est.rd.batch <- function( od_by_stage ){
	# od_by_stage list of out degree by stage for multiple experiments
	# : list of lists with 5 components for each stage
	o <- list()
	for (iobs in 1:length(od_by_stage ))
	{
		obs <- od_by_stage[[iobs]] 
		rrs <- rate.difference.ehi(obs , nreps=1e3 )
		o[[iobs]] <- rrs 
		print(date())
		print( summary( rrs ))
		print(quantile( rrs, prob = c(.025, .5, .975)))
		cat('\n\n')
	}
	ix <- order ( sapply( o, median))
	o[ix] # sort in order of increasing median rates
}
