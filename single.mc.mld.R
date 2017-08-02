##- toy example of W_ij
if(FALSE){
  d <- as.character(1:10)
  r <- as.character(1:3)
  W <- expand.grid('donor' = d, 'recip' = r, stringsAsFactors = FALSE)
  W$ip <- runif(dim(W)[1])/length(d)
  tapply(W$ip, W$recip, sum) # in-degree
}


##---- infector probs ----
## from a list, with split in df
single.mc.mld.ip <- function(x = W){ # x = edge list c(donor, recip, ip)
  nd <- length(unique(x$donor))
  ud <- setNames(rep(0, nd + 1), 
                 c(unique(x$donor), 'out'))# list of unique donors
  ur <- unique(x$recip) # unique recip
  ##- list of df named by recip
  df.r <- split(x, x$recip)
  ##- sample within donors' names + "outside source", with infector probs + (1 - indegree)
  mld <- sapply(df.r, function(r){
    sample(c(r$donor, 'out'), size = 1,
           prob = c(r$ip, max(0, 1 - sum(r$ip)) ) )
  })
  n.mld <- setNames(as.vector(table(mld)), names(table(mld))) # count
  ## add count of being mld to zero vector of all donors
  v <- c(ud, n.mld)
  new.ud <- tapply(v, names(v), sum)
  stopifnot(sum(new.ud) == length(ur)) # one donor found by recip
  return( setNames(as.vector(new.ud), names(new.ud)) )
}

##- toy example of distance edge list
if(FALSE){
  d <- as.character(1:10)
  r <- as.character(1:10)
  D <- expand.grid('ID1' = d, 'ID2' = r, stringsAsFactors = FALSE)
  D$distance <- runif(dim(D)[1])/length(d)/2
  D <- D[D$distance < 0.05,]
  rownames(D) <- NULL
  ## say 3 are incident
  index_incident <- sample(1:10, 3)
  print(paste(c(index_incident, 'are incident'), collapse = ' '))
  thr <- 0.015
}

##---- distance-based cluster ----
single.mc.mld.clust <- function(edgelist = D, # in form: ID1, ID2, distance
                             ids = index_incident,
                             thr = 0.015){
  u <- unique(c(edgelist[,1], edgelist[,2]))
  ud <- setNames(rep(0, length(u)), u)
  ids <- ids[ids %in% c(edgelist[,1], edgelist[,2])] # only incident recipient potentially in cluster
  #mld <- vector(mode = 'character', length = length(ids))
  for(i in 1:length(ids)){
    sub <- edgelist[ (edgelist[,1] == ids[i] | edgelist[,2] == ids[i]) & edgelist[,3] < thr, ] # only rows involving i and below threshold
    nbh <- unique(c(sub[sub[,'ID1']!=ids[i],'ID1'], sub[sub[,'ID2']!=ids[i],'ID2'])) # neighborhood
    if(length(nbh)!=0) {
      mld <- sample(nbh, 1) # select one nbh as donor
      ud[mld] <- ud[mld] + 1
    }
    print(paste(i, '/', length(ids)))
  }
  return(ud)
}

# mld <- single.mc.mld.clust(edgelist = D, ids = index_incident, thr = 0.015)
mlds <- sapply(1:2, function(x) {
  print(x)
  single.mc.mld.clust(edgelist = D, ids = index_incident, thr = 0.015)
})
# mlds <- replicate(10, single.mc.mld.clust(edgelist = D_inc, ids = index_incident))


