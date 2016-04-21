CL_THRESHOLD <- .015

require(phydynR)
source('model0.R')
require(igraph) 

ratesBaselineFNS <- list.files('RData', full.names=T, path = 'model0-simulateBaseline0')

#~ fn <- ratesBaselineFNS[1]
fn2clust <- function(fn)
{
	load(fn)	
	#load distance mat
	samplesInCohort <- daytree$tip.label
	.dfn <-  strsplit( fn, '/' )[[1]]
	dfn <- paste( sep='/', paste(sep='-', .dfn[1], 'distances'), .dfn[2])
	load(dfn) #loads D
	
	i <- which( D[3,] < CL_THRESHOLD )
	.df <- data.frame( 
		   src = daytree$tip.labe[  floor(D[1,i]) ] 
		 , dst = daytree$tip.labe[  floor(D[2,i]) ] 
		)
	
	g <- simplify(make_graph( as.matrix(.df), directed=FALSE, isolates = setdiff( daytree$tip.label, union(.df$src, .df$dst) ) ) )
	
	comps <- components( g ) 
	
	csizes <- setNames(rep(1, daytree$n), daytree$tip.label)
	grps <- groups(comps)
	for (grp in grps){
		csizes[ grp] <- length(grp)
	}
	
	
	print(date())
	print(fn)
	
	list( 
	  df = data.frame( csizes = csizes, nbrhoodSize =  degree(g)[daytree$tip.label])
	, g =g
	)
}

#~ x <- fn2clust( ratesBaselineFNS[1] )
