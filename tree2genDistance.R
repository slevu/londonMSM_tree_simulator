require(phydynR)

.tree2genDistance <- function(tree, mu = .0015, seqlength = 1e3, MH=20 ){
#~ tree : DatedTree, branch lenghts in years 
#~ mu : clock rate (subst/site/year)
#~ MH: dont bother computing distances when TMRCA beyond MH years
	tree$edge.length <- rpois(length(tree$edge.length), tree$edge.length * mu * seqlength ) / seqlength 
#~ 	NOTE the following does not work, because of insufficient memory:
#~ 	cophenetic.phylo( tree )# mem problems
#~ SO, we will do it the hard way: 
	n <- length(tree$tip.label)
	Nnode <- tree$Nnode 
	D <- matrix( NA, nrow = 3, ncol = n*(n-1)/2 ) #BIG 
	d2tip <- lapply( 1:(n+Nnode), function(i) c() ) # c( tip , distance to tip )
	for (u in 1:n){
		d2tip[[u]] <- matrix( c( u, 0. ), ncol = 1, nrow = 2)
	}
	daughterEdgeLengths <- cbind( tree$daughters
	 , t( sapply( 1:(n+Nnode), function(a){
		if (a <= n ) return( c(NA, NA ))
	    u <- tree$daughters[a,1]
	    v <- tree$daughters[a,2]
	    elu <- tree$edge.length[ which( tree$edge[,2]==u ) ]
	    elv <- tree$edge.length[ which( tree$edge[,2]==v ) ]
	    c( elu, elv)
	 }) )
	)
	
	# order to traverse nodes
	o <- rep(-1, n + Nnode )
	r <- which( is.na( tree$parent))
	o[r] <- 0
	ng <- tree$daughters[r,] 
	k <- 0
	while( length(ng) > 0 )
	{
		.ng <- c()
		for ( u in ng){
			o[u] <- k
			k <- k + 1
			.ng <- c( .ng, setdiff(tree$daughters[u,],1:n) )
		}
		.ng <- .ng[!is.na(.ng)]
		ng <- unique( .ng) 
	}
	o <- setdiff(order(o, decreasing=TRUE) , 1:n )
	
	k <- 1	
	for (a in o){
		if (tree$heights[a] < MH ){
			u <- tree$daughters[a,1]
			v <- tree$daughters[a,2]
			d2tip_u <- d2tip[[u]]
	if (length(d2tip_u)==0) browser()
	if (nrow(d2tip_u)<2) browser()
			d2tip_u[2,] <- d2tip_u[2,] + daughterEdgeLengths[a,3] 
			d2tip_v <- d2tip[[v]]
	if (length(d2tip_v)==0) browser()
	if (nrow(d2tip_v)<2) browser()
			d2tip_v[2,] <- d2tip_v[2,] + daughterEdgeLengths[a,4]
			d2tip[[a]] <- cbind( d2tip_u, d2tip_v )
			# update distances between tips
			for (iw in 1:ncol(d2tip_u) ) for (iz in 1:ncol(d2tip_v) ){
				w <- d2tip_u[1,iw]
				z <- d2tip_v[1,iz]
				D[1,k] <- w
				D[2,k] <- z
				D[3,k] <-  d2tip_u[2,iw] + d2tip_v[2,iz]
				k <- k + 1
			}
			print(paste(a, k, date()))
		}
	}
	D <- D[,1:(k-1)]
	D
}
