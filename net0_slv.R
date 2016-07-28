# compute mean inf probs for pairs in cluster 
# size : out degree
# color : age (heatmap) NOTE old = dark, young = light
# shape : high CD4 or not 
# Edge width : W 

THRESHOLD_W <- .01 
cd4s <- read.csv('../phylo-uk/source/cd4s.csv' )
dem <- read.csv('../phylo-uk/source/demographics.csv')

## load cluster sids 
require(ape)
tre <- read.tree( './Erik_code/clusterForBEAST0.interleaved_phyml_tree.txt')
cl_sids <- sapply(strsplit( tre$tip.label 
 , '_'), function(sp) sp[1] )
Wfns <- list.files("data/phydynR", pattern = "mh20", full.names = TRUE)[1:100] # could do regex to exclude 'rita' files

# combine
n <- length(cl_sids)
W <- matrix( 0, nrow=n, ncol =n)
rownames(W) = colnames(W) <- cl_sids
for ( wfn in Wfns ){
	w <- readRDS( wfn )
	i <- which ( w[[1]] %in% cl_sids  & w[[2]] %in% cl_sids)
	for (ii in i){
		W[ w[[1]][ii], w[[2]][ii]] <- W[ w[[1]][ii], w[[2]][ii]] + w[[3]][ii]
	}
}
W <- W / length( Wfns )


require(igraph)
# make edges
donor <- rep(NA, sum(W > THRESHOLD_W ) )
recip <- rep(NA, sum(W > THRESHOLD_W ) )
w <- rep(NA, sum( W> THRESHOLD_W ))

k <- 1;
for (u in cl_sids) for (v in cl_sids){
	if (W[u,v] > THRESHOLD_W) {
		donor[k] <- u
		recip[k] <- v
		w[k] <- W[u,v]
		k <- k + 1
	}
}
edges <- data.frame( donor = donor, recip = recip, w = w )
verts <- unique( c(donor, recip))
dobs <- as.numeric( dem$dob_y[ match(verts, dem$patientindex) ]  )
ages <-  setNames( 1980 + as.numeric( as.Date('1980-11-11') - as.Date('1980-01-01')) / 365 - dobs
 , verts) #approximate
cd40 <- cd4s$cd4[ match(verts, cd4s$patientindex ) ] > 500#700
cd40[is.na( cd40)]  <- FALSE
require(RColorBrewer)
pal <- brewer.pal( 9, "Blues")
age2col <- function(a){
	k <- min(9, max(1, floor( ecdf( ages )(a)*10 ) ) ) 
	if (is.na(k)) return( pal[4] )
	pal[ k ]
}
verts <- data.frame( vertices = verts
 , age = ages
 , cd40 = cd40
 , color = sapply( verts, function(w) age2col( ages[w] ) )
 , shape = sapply( cd40, function(x) c("circle", "square")[x+1] )
)
gdf <- graph.data.frame( edges, verts , directed = TRUE)

comps <- components( gdf )
compindex <- which.max( comps$csize )
comp <- names( comps$membership )[ comps$membership==compindex ]
subgdf <- induced_subgraph( gdf, match( comp, names( V(gdf))) )


## now plot subgdf

# size
od <- ( sapply( names(V(subgdf)), function(w){
	sum( edges$w[ edges$donor==w] )
}))
od <- od / max(od ) + .1


subgdf_el <- as_edgelist( subgdf )
ewidth <- sapply( 1:nrow(subgdf_el), function(k){
	w <- subgdf_el[k,1]
	z <- subgdf_el[k,2]
	edges$w[ edges$donor==w & edges$recip==z ]
})
ewidth <- sqrt( ( ewidth / median(ewidth) ) )


vsize <- sqrt(od)
vsize <- vsize / max(vsize)
vsize <- vsize * 15

#~ l <- layout.fruchterman.reingold(subgdf , weight=w)
#~ l <- layout.kamada.kawai(subgdf)
#~ l <- layout_with_fr( subgdf, weights = 1/ewidth^4 )
#~ l <- layout_with_fr( subgdf, weights = ewidth^4 * 1, niter = 1e4, dim =2 )
#~ l <- layout_with_kk( subgdf, weights = 1/ewidth )
#~ l <- layout_with_kk( subgdf, weights = ewidth^1.25 * 100) #*
#~ l <- layout_with_mds( subgdf)
l <- layout_with_dh( subgdf) #**


plot( subgdf
  , vertex.size= vsize / 1.5
  , edge.width = ewidth
  , edge.color='Black'
  , edge.arrow.size=.5
  , vertex.label=NA 
  , edge.curved=.4
  , layout=l
)



## now plot evo dist cluster using same layout
# get distances from phyml tree
THRESHOLD_D <- 0.015
vnames <- names(V(subgdf))
D <- cophenetic.phylo(tre)
rownames(D) <- sapply( strsplit( rownames(D), '_' ), function(spl) spl[1] )
colnames(D) <- sapply( strsplit( colnames(D), '_' ), function(spl) spl[1] )
D <- D[vnames,vnames]
diag(D) <- Inf
edgesD <- matrix(NA, nrow = sum(D < THRESHOLD_D), ncol = 2)
k <- 1
for (w in vnames) for (z in vnames){
	if (D[w,z] < THRESHOLD_D){
		edgesD[k, 1] <- w
		edgesD[k, 2] <- z
		k <- k + 1
	}
}
vertsD <- verts[verts$vertices %in% names(V(subgdf)),]
gD <- graph.data.frame( edgesD, vertsD , directed = FALSE)
plot( gD
  , vertex.size= 5
  , edge.width = .2
  , edge.color='Black'
  , vertex.label=NA 
  , edge.curved=.4
  , layout=l
)


## plot them together
pdf('net0_slv.0.pdf', 6, 3.8 )
par(mfrow = c(1,2)
  , mar = rep(0, 4)
)
plot( subgdf
  , vertex.size= vsize
  , edge.width = ewidth
  , edge.color='Black'
  , edge.arrow.size=.25
  , vertex.label=NA 
  , edge.curved=.4
  , layout=l
  , margin=rep(0, 4)
  , frame=T
#~   , asp=.5
)
plot( gD
  , vertex.size= 7
  , edge.width = .05
  , edge.color='Black'
  , vertex.label=NA 
  , edge.curved=.4
  , layout=l
  , margin=rep(0, 4)
  , frame=T
#~   , asp=.5
)
dev.off()


# save plots of each individually
pdf('net0.w.pdf', 4, 3.8 )
par(mar = rep(0, 4)
)
plot( subgdf
  , vertex.size= vsize
  , edge.width = ewidth
  , edge.color='Black'
  , edge.arrow.size=.25
  , vertex.label=NA 
  , edge.curved=.4
  , layout=l
  , margin=rep(0, 4)
  , frame=T
#~   , asp=.5
)
dev.off()


pdf('net0.cl.pdf', 4, 3.8 )
par(mar = rep(0, 4)
)
plot( gD
  , vertex.size= 7
  , edge.width = .05
  , edge.color='Black'
  , vertex.label=NA 
  , edge.curved=.4
  , layout=l
  , margin=rep(0, 4)
  , frame=T
#~   , asp=.5
)
dev.off()
