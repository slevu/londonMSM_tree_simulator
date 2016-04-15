# recomputing inf probs using RITA where available
require(phydynR)

MH <- 20 
newinf <- 2500 # c(1660, 4780)
plwhiv <- 43150 / 2 # c(43510 / 2, 43510 / 1.5)


cd4s.df <- read.csv('ukdrd_data/cd4s.csv' )
res.df <- read.csv('ukdrd_data/resistance.csv')
dem.df <- read.csv('ukdrd_data/demographics.csv')

# drop tips that have long cd4 sample time mismatch 
tre <- read.tree( 'LSD/t100.nwk.result.date.newick')
if (T){
	require(doMPI)
	cl <- startMPIcluster(count = 3)
	registerDoMPI(cl)
	todrop <- foreach( sid=tre$tip.label, .combine=cbind, .inorder=TRUE) %dopar% {
		siddt <- as.Date( paste(sep='/', '15', res.df$dbsample_my[sid==res.df$testindex]), format = '%d/%m/%Y')
		pid <- res.df$patientindex[res.df$testindex==sid]
		cd4s.df_ipid <- which( cd4s.df$patientindex==pid )
		cd4dts <- lapply( cd4s.df$cd4_date_my[cd4s.df_ipid], function(dtstring) {
			as.Date( paste(sep='/', '15', dtstring), format = '%d/%m/%Y')
		})
		dt_deltas <- as.vector( sapply(cd4dts, function(cd4dt) {
			abs( cd4dt - siddt )
		} ))
		
		xcd4 <- NA 
		if (length( dt_deltas) == 0) {
			dt_delta <- Inf 
		} else{
			dt_delta <- min(dt_deltas)
			xcd4 <- cd4s.df$cd4[ cd4s.df_ipid[ which.min(dt_deltas) ] ] 
		}
		#print(c( xcd4, dt_delta, dt_deltas) )
		
		#ifelse (dt_delta > 365 /2 , TRUE, FALSE) 
		if (dt_delta > 365 /2){
			 return( c( FALSE, xcd4 ))
		} else{ 
			return ( c(TRUE, xcd4) )
		}
	}
	closeCluster(cl)
	colnames( todrop) <- tre$tip.label
	saveRDS( todrop, file = 'tips2remove_cd4.rds' )
}

# check which RiTA valid
if (T)
{
	res.df2 <- res.df[res.df$testindex %in% tre$tip.label,]
	pids <- res.df2$patientindex[match(tre$tip.label, res.df2$testindex)] 
	siddts <- sapply(res.df2$dbsample_my[match(tre$tip.label, res.df2$testindex) ], function(sdt){
		as.Date( paste(sep='/', '15', sdt ), format='%d/%m/%Y')
	} )
	dem.df2 <- dem.df[ dem.df$patientindex %in% pids , ]
	diagdts <-  sapply( dem.df2$hivpos_my[ match( pids, dem.df2$patientindex ) ], function( sdt ){
		as.Date( paste(sep='/', '15', sdt ), format='%d/%m/%Y')
	})
	RITA_valid <- data.frame( sequenceID = tre$tip.label
	 , patientID = pids
	 , RITASampledCloseToSequence =  abs(siddts - diagdts) < 30 
	)
	saveRDS( RITA_valid, file = 'RITASampledCloseToSequence.rds' )
}
