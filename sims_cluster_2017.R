# rm(list=ls())
####---- include ----
detail_knitr <- TRUE
source("functions.R")
startover <- TRUE

####---- lib ----
#library(ape)

##--- thrs
thresholds <- c("0.001", "0.005", "0.015", "0.05") # c("1e-05", "1e-04", "0.001", "0.005", "0.015", "0.05") 

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

##---- check
# str(list.sims); str(list.dist) # head(names(list.sims[[2]])) 
# load(list.sims[[2]][[1]]); head(sampleDemes); head(cd4s)

## plot distance
## check uk tree distances
duktree <- readRDS("data/uktree_dist.rds")
hist(duktree, breaks = 100, xlim = c(0,1))
dukTN93 <- read.csv("../phylo-uk/test/tn93_average/MSM_B_uk_tn93.csv") # dukTN93 <- readRDS(file = "../phylo-uk/source/subUKogC_noDRM_151202_ucsdTN93.rds" )
hist(dukTN93[,3], breaks = 100, xlim = c(0,1))

load(list.dist[[1]][[2]])
dim(D)
D[, 1:6]
hist(D[3,], breaks = 100)
load(list.sims[[1]][[2]])
library(ape)
#~ tree : DatedTree, branch lenghts in years 
#~ mu : clock rate (subst/site/year)
#~ MH: dont bother computing distances when TMRCA beyond MH years
tree <- bdt
mu = 4.3e-3 #.003
seqlength = 1e3;
# (year * susbt/site/year * n site) = number of substitutions 
tree$edge.length <- rpois(length(tree$edge.length), tree$edge.length * mu * seqlength ) / seqlength 
dsimtree <-  as.dist(cophenetic.phylo(tree))
summary(dsimtree)
hist(dsimtree, breaks = 100, xlim = c(0,1))

source("tree2genDistance.R")
dsimtree2 <- .tree2genDistance(bdt, mu = .0015, seqlength = 1e3, MH=20 )
summary(dsimtree2[3,])
hist(dsimtree2[3,], breaks = 100)
#plot(daytree)
summary(daytree$edge.length)
summary(tree$edge.length)
summary(bdt$edge.length)
summary(tree$heights)
hist(bdt$edge.length)


##-- test two functions for tree to distance
if (FALSE) {
  rt <- rtree(100)
  if (FALSE){
    ##- hack: exponential transform to get edge.length similar to bdt (year tree)
    s <- sample(bdt$edge.length, 198)
    y <- s[order(s)]
    x <- rt$edge.length[order(rt$edge.length)]
    expmod <- lm(log(y + 10e-6) ~ x)  # fit exponential relation
    parm <- expmod$coef
    summary(rt$edge.length)
    rt$edge.length <- exp(parm[1] + parm[2]*rt$edge.length)
    summary(rt$edge.length)
  }
  
  
  hist(rt$edge.length)
  tree2copheneticGenDist <- function(yeartree, mu = .0015, seqlength = 1e3, MH=20){
    yeartree$edge.length <- rpois(length(yeartree$edge.length), yeartree$edge.length * mu * seqlength ) / seqlength 
    m <- as.matrix( as.dist(cophenetic.phylo(yeartree)) )
    ###--- keep only the lower triangle by 
    ## filling upper with NA
    m[upper.tri(m, diag=TRUE)] <- NA
    ###--- create edge list if not there
    require(reshape2)
    el <- melt(m, na.rm = TRUE)
    colnames(el) <- c('ID1', 'ID2', 'distance')
    el
  }
  D1 <- tree2copheneticGenDist(rt)
  hist(D1[,3])
  rt$heights <- node.height(rt)
  D2 <- .tree2genDistance(rt, mu = .0015, seqlength = 1e3, MH=20 ) # doesn't work without tree$heights
} # not used

####---- ucsd clustering ----####
# ucsd_hivclust
if (startover == TRUE){
  thresholds  <-  as.numeric(thresholds) # c(0.00001, 0.0001)# c(1e-04, 5e-04) # c("0.001", "0.05") # c("0.015", "0.005") # c(0.005, 0.015, 0.02, 0.05) # c(0.005, 0.01, 0.02, 0.05, 0.1) 
  tmax <- max(thresholds) # limit of distance considered
  
## function: input list of dist filenames, output csv of clusters
#- debug: ldist = list.dist[["Baseline0"]]; i = 100
  ucsd <- function(ldist){
    
    for (i in 1:length(ldist)){
    ##- some processing
    load(ldist[i])
    name.sim <- regmatches(ldist[i], regexpr("[0-9]{3,}", ldist[i]))
    folder.sim <- substr(ldist[i], regexec("-simulate", ldist[i])[[1]][1] + 9, regexec("-distances", ldist[i])[[1]][1] -1)
    dd <- as.data.frame(t(D))
    names(dd) <- c('ID1', 'ID2', 'distance')
    ##- temp el
    temp.el.fn <- paste(tempdir(), "/", name.sim, "_el.rds", sep = '')
    saveRDS(dd, file = temp.el.fn )
    ##- ucsd hivclust
    ucsd_hivclust(path.el = temp.el.fn,
                          thr = thresholds, 
                          k = tmax, 
                          out = paste(path.results, folder.sim,'', sep = '/' ) )
    }
  } 
  
  lapply(list.dist[["Baseline0"]], ucsd)
  lapply(list.dist[["EqualStage0"]], ucsd)
}

####---- list.hivclust ----
if (startover == TRUE){
##-  get csv of clusters assignments in one list
    list.csv <- lapply(scenario, function(x){
    list.files("csv", full.names=T, 
               path = paste(path.results, x, sep = '/'))
  })
  

##- function n = 100; m = 1
list.hivclust <- function(csvs){
  ## Structure threshold > trees

  ## empty list of thresholds
  cl2 <- vector("list", length(thresholds))
  ## loop
    for (m in 1:length(thresholds) ){ # index of thr
     ## vector of csv at different tree for one thr
      filenames <- csvs[grep(thresholds[m], csvs)]
      ## empty list of different trees
      t <- vector("list", length(filenames))
      for (n in 1:length(filenames) ){ # index of trees
        t[[n]] <- read.csv( filenames[n] )
        names(t)[n] <- substr(filenames[n], 
                              regexec("/d", filenames[n])[[1]][1] + 2,
                              regexec("_ucsd_hivclust_output",
                                      filenames[n])[[1]][1] -1)
      }
      cl2[[m]] <- t
      names(cl2)[m] <- thresholds[m]
    }
return(cl2)
  }
  
  cl_Baseline0 <- list.hivclust(csvs = list.csv[["Baseline0"]] )
  cl_EqualStage0 <- list.hivclust(csvs = list.csv[["EqualStage0"]] )

# saveRDS(cl_Baseline0, file = paste(path.results, 'list.hivclust.sim.Baseline0.rds', sep = '/') )
# saveRDS(cl_EqualStage0, file = paste(path.results, 'list.hivclust.sim.EqualStage0.rds', sep = '/') )
} else {
  ##- load list.hivclust
  cl_Baseline0 <- readRDS( file = paste(path.results, 'list.hivclust.sim.Baseline0.rds', sep = '/') )
  cl_EqualStage0 <- readRDS( file = paste(path.results, 'list.hivclust.sim.EqualStage0.rds', sep = '/') )
}
# names(cl_EqualStage0[[1]]) 
# str(cl_Baseline0[[1]][1])

####---- nbhood size ----
  if(startover){
  ###-- start nbh ---
  ### ldist = list.dist[["Baseline0"]]; i = 1; j = 1 
  nbh <- function(ldist){
    ##- empty list of j thresholds * i sims
    thr <- as.numeric(thresholds)
    ll <- rep( list( vector("list", length(ldist)) ), length(thr) ) 
    names(ll) <- thresholds
    
    for (i in 1:length(ldist)){
      ##- load distances
      load(ldist[i])
      dd <- as.matrix(t(D))
      name.sim <- regmatches(ldist[i], regexpr("[0-9]{3,}", ldist[i]))
      
      ##- calculate neighborhood size
      ##- number of neighbour|threshold
      for (j in 1:length(thresholds)){
        print( paste(i, j, name.sim) )
        
        .t1  <- tapply(dd[,3], dd[,1], function(x) sum(x < thr[j]))
        .t2  <- tapply(dd[,3], dd[,2], function(x) sum(x < thr[j]))
        ##- add from and to neighbors
        .t <- tapply(c(.t1, .t2), names(c(.t1, .t2)), sum)
        
        ll[[j]][[i]] <- cbind("id" = as.numeric(names(.t)), "nbhsize" = unname(.t))
        names(ll[[j]])[i] <- name.sim
      }
    }
    return(ll)
  }
  ##--- end nbh ---
  
   system.time(
    nbh_Baseline0 <- nbh(list.dist[["Baseline0"]] )
  ) # 55s
  system.time(
    nbh_EqualStage0 <- nbh(list.dist[["EqualStage0"]] )
  ) # 40s

  # saveRDS(nbh_Baseline0, file = paste(path.results, 'list.nbhsize.sim.Baseline0.rds', sep = '/') )
  # saveRDS(nbh_EqualStage0, file = paste(path.results, 'list.nbhsize.sim.EqualStage0.rds', sep = '/') )
  } else {
    ##- load list.hivclust
    nbh_Baseline0 <- readRDS( file = paste(path.results, 'list.nbhsize.sim.Baseline0.rds', sep = '/') )
    nbh_EqualStage0 <- readRDS( file = paste(path.results, 'list.nbhsize.sim.EqualStage0.rds', sep = '/') )
  }
  ##- check
  #   head(nbh_Baseline0[[4]][[1]])
  #   a <- lapply(nbh_Baseline0, function(x) do.call(rbind, x))
  #   str(a)
  #   sapply(a, function(x) mean(x[,2]))
  
  ##- Todo
  ##- convert tip ID in tip label
  ##- incorporate in clust.stats
  
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

## At this stage, list of cluster assignements have an index of tips not tip labels
####---- clus.stat ----
if (startover == TRUE){
  ###- check if sim.name = clus.name
  # sim = list.sims[["Baseline0"]] ; clus = cl_Baseline0; nbh = nbh_Baseline0; i = 1 ; thr = "0.015"
clus.stat <- function(clus, nbh, sim){
  
  ##- Start loop
  ##- empty list of n thresholds
  ll <- vector("list", length(clus))
  names(ll) <- names(clus)
  ##- loop threshold
  for (thr in names(clus)){
    ##- loop sims
    for (i in 1:length(clus[[thr]]) ){
      ##- check same sim and clus
        print( paste( 
          names(clus[[thr]])[i], # name cluster
          thr, # threshold
          i # num sim
          ) )
      ## load cluster assignement and nbhood size
      cl <- clus[[thr]][[i]]
      nb <- as.data.frame(nbh[[thr]][[i]])
     
      ##- load sim
      .m <- grep(paste0('/', names(clus[[thr]])[i], '.RData'), sim )
      load( sim[.m] )
      
      ## get ALL tips names and demes
      tip.states <- data.frame(
        "id" = daytree$tip.label,
        "stage" = sapply( sampleDemes, deme2stage ),
        "age" = sapply( sampleDemes, deme2age ),
        "risk" = sapply( sampleDemes, deme2risk ),
        stringsAsFactors = FALSE)
      rownames(tip.states) <- NULL
      
        ## get size of clusters
        freqClust <- as.data.frame(table(cl$ClusterID), stringsAsFactors = FALSE)
        
        ## get name instead of index in clus
        cl$id <- daytree$tip.label[ cl$SequenceID ]
        nb$id <- daytree$tip.label[ nb$id ]
        
        ## merge cluster number (without SequenceID)
        a <- merge(x = tip.states, y = cl[, 2:3],
                   by.x = "id", by.y = "id",
                   all.x = TRUE, sort = FALSE)
        
        ## merge cluster size (with NA)
        aa <- merge(x = a, y = freqClust, 
                   by.x = "ClusterID", by.y = "Var1", 
                   all.x = TRUE, sort = FALSE)
        
        ## merge neighborhood size (with NA)
        b <- merge(x = aa, y = nb, by.x = "id", by.y = "id", 
                    all.x = TRUE, sort = FALSE)
        
        #- size 1 if not into a cluster
        b$Freq[is.na(b$Freq)] <- 1
        
        #- size 0 if NA
        b$nbhsize[is.na(b$nbhsize)] <- 0
        
        #- binary clustering variable
        b$binclus <- ifelse(b$Freq > 1 & !is.na(b$Freq), 1, 0)
        # b[sample(nrow(b),10),]
        
        ##- pseudo ClusterID when size = 1
        ## starting from max(ClusterID+1)
        ## important for down-sampling
        .start <- max(b$ClusterID, na.rm = T)
        .n <- dim(b[is.na(b$ClusterID),] )[1]
        b$ClusterID[is.na(b$ClusterID)]  <- .start + 1:.n 
        
        #- colnames
        colnames(b)[which(colnames(b) =="Freq")] <- "size"
        ll[[thr]][[i]] <- b
        names(ll[[thr]])[i] <- names(clus[[thr]][i])
    } 
  }
  return(ll)
}
####---- fin clus.stat ----####
}

if (startover == TRUE){
  system.time(
    l_Baseline0 <- clus.stat(clus = cl_Baseline0, 
                             nbh = nbh_Baseline0,
                             sim = list.sims[["Baseline0"]])
  ) # 766
  system.time(
    l_EqualStage0 <- clus.stat(clus = cl_EqualStage0,
                               nbh = nbh_EqualStage0,
                               sim = list.sims[["EqualStage0"]])
  ) # 750
  
}

  ####---- saved listUKclus ----
 # saveRDS(l_Baseline0, file = paste(path.results, 'list.sim.ucsd.Baseline0.rds', sep = '/') )
 # saveRDS(l_EqualStage0, file = paste(path.results, 'list.sim.ucsd.EqualStage0.rds', sep = '/') )

  l_Baseline0 <- readRDS(file = paste(path.results, 'list.sim.ucsd.Baseline0.rds', sep = '/') )
  l_EqualStage0 <- readRDS(file = paste(path.results, 'list.sim.ucsd.EqualStage0.rds', sep = '/') )

###- check
  #   names(l_Baseline0)
  #   str( l_Baseline0[[3]][[1]] )
  #   length(l_Baseline0[[1]])
# 
#   system.time(
#   z <- lapply(l_Baseline0, function(x){
#     sapply(x, function(m) {
#     # faster than `merge(df, s, all.x = TRUE)`
#     #- add covariates
#      # m <-  cbind(df, s[ match(df$id, s$id), -1 ])
#     #- mean y by x
#       nb <-  tapply(m$nbhsize, m$stage, mean )
#       return(nb)
#       }) 
#     })
#   )
#   sapply(z, function(x) apply(x,1, mean))

####---- add degrees ----

  ##- loop to merge all W to all clusters
#   clus = l_Baseline0; sim = list.sims[["Baseline0"]];  i <- 2

  ####---- add.w ----
  add.w <- function(clus, sim){
    ##- empty list of n sims (of 1st threshold level)
    cl <- clus[[1]]
    ll <- vector("list", length(cl))
      ##- loop sims
      for (i in 1:length(cl) ){

        ##- monitor loop
        print( paste(names(cl)[i], # name sim
                    i ) ) # num sim
        
        ##- load sim with W matrix of infect probs
        .m <- grep(paste0('/', names(cl)[i], '.RData'), sim )
        load( sim[.m] )
        
        ##- outdegree
        outd <- tapply(W$infectorProbability, W$donor, sum)
        ##- indegree
        ind <- tapply(W$infectorProbability, W$recip, sum)
        ##- merge out and in degrees
        b <-  data.frame('patient' = names(outd), 
                         'outdegree' = as.vector(outd), 
                         'indegree' = as.vector(ind), 
                         row.names = NULL, stringsAsFactors = FALSE)
      
        ##- load covariates by tip
        a <- cl[[i]][, c(1,3:5)]
        
        ##- merge all
        c <- merge(a, b, by.x = "id", by.y = "patient", all.x = TRUE )
        
        ##- NA degree to 0
        c$outdegree[is.na(c$outdegree)] <- 0
        c$indegree[is.na(c$indegree)] <- 0
        
        ##- names
        ll[[i]] <- c
        names(ll)[i] <- names(cl)[i]
      }
    ##- concatenate lists SA + thresholds
    return( c("SA" = list(ll), clus) )
    }

  ####---- fin add.w ----####
  
  
  if (startover = TRUE){
    system.time(
      cw_Baseline0 <- add.w(clus = l_Baseline0, 
                            sim = list.sims[["Baseline0"]] )
    ) # 58
    system.time(
      cw_EqualStage0 <- add.w(clus = l_EqualStage0,
                              sim = list.sims[["EqualStage0"]] )
    ) # 52
  }
  
  
  ####---- saved listUKclus ----
  # saveRDS(cw_Baseline0, file = paste(path.results, 'list.sim.clus-outdeg.Baseline0.rds', sep = '/') )
  # saveRDS(cw_EqualStage0, file = paste(path.results, 'list.sim.clus-outdeg.EqualStage0.rds', sep = '/') )
  
  cw_Baseline0 <- readRDS(file = paste(path.results, 'list.sim.clus-outdeg.Baseline0.rds', sep = '/') )
  cw_EqualStage0 <- readRDS(file = paste(path.results, 'list.sim.clus-outdeg.EqualStage0.rds', sep = '/') )
  
  ##---- end ----
  