CL_THRESHOLD <- .015
local_data <- "data/simulations/"
require(phydynR)
source('model0.R')
require(igraph) 

ratesBaselineFNS <- list.files('RData', full.names=T, path = paste(sep = '', local_data,'model0-simulateBaseline0'))

#~ fn <- ratesBaselineFNS[5]
fn2clust <- function(fn)
{
	load(fn)	
	#load distance mat
	samplesInCohort <- daytree$tip.label
	.dfn <-  strsplit( fn, '/' )[[1]]
	dfn <- paste( sep='/', .dfn[1], .dfn[2], paste(sep='-', .dfn[3] , 'distances'), .dfn[4])
	load(dfn) #loads D
	
	i <- which( D[3,] < CL_THRESHOLD )
	.df <- data.frame( 
		   src = daytree$tip.label[  floor(D[1,i]) ] 
		 , dst = daytree$tip.label[  floor(D[2,i]) ] 
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
head(x[[1]])
class(x[[2]])
hist(x[[1]]$csizes)
hist(x[[1]]$nbrhoodSize)

head(t(D))
hist(t(D)[,3])
## unrealistic distances < 0.07 ??