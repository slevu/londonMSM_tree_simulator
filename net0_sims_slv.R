##---- network figures with simulated data ----
# rm(list=ls())

##- for one simulation, get W and ucsd cluster at 5%
####---- path sims ----
##- used for several rounds of simulations
path.sims <- 'data/simulations2/model0-simulate'
path.results <- 'data/sim_ucsd_results2'

####---- scenario ----
scenario <- c("Baseline0", "EqualStage0")
scenario <- setNames(scenario, scenario) # useful to name list in lapply

####--- list of sims files and distances files ---
list.sims <- lapply(scenario, function(x){
  list.files('RData', full.names = TRUE, 
             path = paste(path.sims, x, sep = '') )
})
list.dist <- lapply(scenario, function(x){
  list.files('RData', full.names = TRUE, 
             path = paste(path.sims, x, '-distances', sep = '') )
})

##- load cluster
cw_Baseline0 <- readRDS(file = paste(path.results, 'list.sim.clus-outdeg.Baseline0.rds', sep = '/') )
names(cw_Baseline0)

##- at 5%, first replicate # second
##---- function to generate different plots ----
netplot <- function(r = 2, thr_D1 = 0.05, thr_D2 = 0.015, thr_W = 0.001){
##- debug
##- r = 5; thr_D1 = 0.05; thr_D2 = 0.015; thr_W = 0.001

  parms <- paste0(r, '_', thr_W)
  ##- select cluster at high threshold
  cl <- cw_Baseline0[[as.character(thr_D1)]][[r]] # at thr = 0.05, r th replicate
  # str(cl)
  ##- find size and largest cluster
  size <- as.data.frame(table(cl$ClusterID), stringsAsFactors = F)
  # head(size[order(size$Freq, decreasing = T),])
  largest_clusterID <- size[order(size$Freq, decreasing = T),]$Var1[1]
  ##- get names from largest cluster
  list_id <- cl[cl$ClusterID == largest_clusterID, 'id'] 

  ##- load corresponding sim and distances
  load(list.sims[['Baseline0']][r])
  load(list.dist[['Baseline0']][r])
  DD <- t(as.matrix(D))
  # head(DD)
  
  ##- get names instead of id in a edge list of distances
  # str(daytree)
  el <- data.frame(donor = daytree$tip.label[ DD[,1] ], 
                   recip = daytree$tip.label[ DD[,2] ],
                   dist = DD[,3])
  # head(el)
  i <- which ( el[,1] %in% list_id  & el[,2] %in% list_id)
  el_cl <- el[i,]
  
  ##- get inf probs
  # str(W)
  ##- change names to use Erik code
  cl_sids <- list_id

  ##--- start net0.R ---
  THRESHOLD_W <- thr_W 
  # cd4s <- read.csv('../phylo-uk/source/cd4s.csv' )
  # dem <- read.csv('../phylo-uk/source/demographics.csv')
  ##- load cluster sids 
  # require(ape)
  # tre <- read.tree( './Erik_code/clusterForBEAST0.interleaved_phyml_tree.txt')
  # cl_sids <- sapply(strsplit( tre$tip.label 
  #  , '_'), function(sp) sp[1] )
  # Wfns <- list.files("data/phydynR", pattern = "mh20", full.names = TRUE)[1:100] # could do regex to exclude 'rita' files
  Wfns <-  list.sims[['Baseline0']]
  ##- combine
  n <- length(cl_sids)
  W0 <- matrix( 0, nrow=n, ncol =n)
  rownames(W0) = colnames(W0) <- cl_sids
  for ( wfn in Wfns ){
	  load( wfn )
	  w <- W
	  i <- which ( w[[1]] %in% cl_sids  & w[[2]] %in% cl_sids)
	  for (ii in i){
		  W0[ w[[1]][ii], w[[2]][ii] ] <- W0[ w[[1]][ii], w[[2]][ii]] + w[[3]][ii]
	  }
  }
  # sum(W0); max(W0)
  W <- W0 / length( Wfns )

  require(igraph)
  ##- add triangle shape (from ?igraph.vertex.shapes )
  # triangle vertex shape
  mytriangle <- function(coords, v=NULL, params) {
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v]
    }
    vertex.size <- 1/200 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v]
    }
    
    symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
            stars=cbind(vertex.size, vertex.size, vertex.size),
            add=TRUE, inches=FALSE)
  }
  # clips as a circle
  add_shape("triangle", clip=shapes("circle")$clip,
            plot=mytriangle)
  
  ##- make edges
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
  dobs <- as.numeric( cl$age[ match(verts, cl$id) ]  )
  ages <-  setNames( dobs , verts) #approximate
  cd40 <- cl$stage[ match(verts, cl$id ) ] # keep stages for different shapes
  # cd40[is.na( cd40)]  <- FALSE
  require(RColorBrewer)
  pal <- brewer.pal( length(unique(ages)), "Blues")
  age2col <- function(a){
    	# k <- min(9, max(1, floor( ecdf( ages )(a)*10 ) ) ) 
  	k <- a
	  if (is.na(k)) return( pal[4] )
	  pal[ k ]
  }
  verts <- data.frame( vertices = verts
   , age = ages
   , cd40 = cd40
   , color = sapply( verts, function(w) age2col( ages[w] ) )
   , shape = sapply( cd40, function(x) c("triangle", rep("circle", 3), "square")[x] ) # assign shape
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

  # thickness of arrows 
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
  THRESHOLD_D <- thr_D2
    vnames <- names(V(subgdf))
# D <- cophenetic.phylo(tre)
# rownames(D) <- sapply( strsplit( rownames(D), '_' ), function(spl) spl[1] )
# colnames(D) <- sapply( strsplit( colnames(D), '_' ), function(spl) spl[1] )
# D <- D[vnames,vnames]
# diag(D) <- Inf
# edgesD <- matrix(NA, nrow = sum(D < THRESHOLD_D), ncol = 2)
# k <- 1
# for (w in vnames) for (z in vnames){
# 	if (D[w,z] < THRESHOLD_D){
# 		edgesD[k, 1] <- w
# 		edgesD[k, 2] <- z
# 		k <- k + 1
# 	}
# }

 ##- cluster corresponding to vertex
  el_cl_r <- el_cl[ el_cl$donor %in% vnames & el_cl$recip %in% vnames,]
  edgesD <- el_cl_r[el_cl_r$dist < THRESHOLD_D, 1:2]
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
  pdf(paste('./figure/net0_sims_slv.', parms,'.pdf',sep = ''), 6, 3.8 )
  
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
pdf(paste('./figure/net0_sims_slv.w.', parms,'.pdf',sep = ''), 4, 3.8 )
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


pdf(paste('./figure/net0_sims_slv.cl.', parms,'.pdf',sep = ''), 4, 3.8 )
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

}

netplot(r = 2)
netplot(r = 3)
netplot(r = 4)
netplot(r = 5)
netplot(r = 5, thr_W = 0.0005)
