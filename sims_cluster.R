# rm(list=ls())
####---- include ----
detail_knitr <- TRUE
source("functions.R")

####---- lib ----
library(ape)

####---- scenario ----
scenario <- c("EqualStage0", "Baseline0")

####--- list.trees ----
list.sims <- vector("list", length(scenario))
for (s in 1:length(scenario)){
  list.sims[[s]] <- list.files('RData', full.names = TRUE, 
             path = paste('data/simulations/model0-simulate', 
                          scenario[s], sep = '') )
  names(list.sims)[s] <- scenario[s]
}
str(list.sims)

####---- list.dist ----
#### load list of dist from sims
if(TRUE){
distEqualStage0FNS <- list.files('RData', full.names=T, 
                           path = 'data/simulations/model0-simulateEqualStage0-distances')
distBaseline0FNS <- list.files('RData', full.names=T,
                              path = 'data/simulations/model0-simulateBaseline0-distances')
}

####---- ucsd clustering ----####
# ucsd_hivclust
if(FALSE){
  thresholds  <-  c(0.005, 0.015, 0.02, 0.05) # c(0.005, 0.01, 0.02, 0.05, 0.1) 
  tmax <- max(thresholds) # limit of distance considered
## function: input list of dist filenames
  ucsd <- function(ldist){
    
    for (i in 1:length(ldist)){
    ##- some processing
    load(ldist[i])
    name.sim <- substr(ldist[i], 
                       regexec(".RData",ldist[i])[[1]][1] - 5, 
                       regexec(".RData", ldist[i])[[1]][1] -1)
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
                          out = paste("data/sim_ucsd_results/", folder.sim, '/', sep = '' ) )
    }
  } 
  
  lapply(distBaseline0FNS, ucsd)
  lapply(distEqualStage0FNS, ucsd)
}

    
####---- list.hivclust ----
if(0){
###  get csv of clusters assignments in one list ###
  # getwd()
csvEqualStage0FNS <- list.files("csv", full.names=T, 
                          path = 'data/sim_ucsd_results/EqualStage0')
csvBaseline0FNS <- list.files("csv", full.names=T, 
                          path = 'data/sim_ucsd_results/Baseline0')
##- function n = 1; m = 1
list.hivclust <- function(list.csv){
  ## Structure threshold > trees
  thresholds <- c("0.005", "0.015", "0.02", "0.05") # c("0.01", "0.02", "0.05")
  ## empty list of thresholds
  cl2 <- vector("list", length(thresholds))
  ## loop
    for (m in 1:length(thresholds) ){ # index of thr
     ## vector of csv at different tree for one thr
      filenames <- list.csv[grep(thresholds[m], list.csv)]
      ## empty list of different trees
      t <- vector("list", length(filenames))
      for (n in 1:length(filenames) ){ # index of trees
        t[[n]] <- read.csv( filenames[n] )
        names(t)[n] <- substr(filenames[n], 
                              regexec("/d", filenames[n])[[1]][1] + 2,
                              regexec("_ucsd_hivclust_output", filenames[n])[[1]][1] -1)
      }
      cl2[[m]] <- t
      names(cl2)[m] <- thresholds[m]
    }
return(cl2)
  }
  
  cl_Baseline0 <- list.hivclust(list.csv = csvBaseline0FNS )
  cl_EqualStage0 <- list.hivclust(list.csv = csvEqualStage0FNS )
 # names(cl_EqualStage0[[1]]) 
  # str(cl_Baseline0[[1]][1])
  
# saveRDS(cl_Baseline0, file = "data/sim_ucsd_results/list.hivclust.sim.Baseline0.rds")
# saveRDS(cl_EqualStage0, file = "data/sim_ucsd_results/list.hivclust.sim.EqualStage0.rds")
}

####---- load list.hivclust ----
  cl_Baseline0 <- readRDS( file = "data/sim_ucsd_results/list.hivclust.sim.Baseline0.rds" )
  cl_EqualStage0 <- readRDS( file = "data/sim_ucsd_results/list.hivclust.sim.EqualStage0.rds" )

  ## not the same number of simulation at threshold 4:
  identical(names(cl_Baseline0[[1]]), names(cl_Baseline0[[4]]))
  cl_Baseline0[[4]] <- cl_Baseline0[[4]][ names(cl_Baseline0[[1]]) ]
  
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
  deme2care <- function(deme){as.numeric( 
    substr(regmatches( deme , regexpr( 'care[0-9]', deme )), 5, 5)
  ) }
  
## At this stage, list of cluster assignements have an index of tips not tip labels

####---- clust.stats ----
 
  ### --- start function clust.stats --- ###
  ##- calculate both numclus and sizeclus for each seqindex into a LIST
  ##- with same variable names
  ##- comprises demes states
  ##- For UCSD files, no ID if no cluster (size < 2) !
# clus = cl_Baseline0[[1]]; i=2; sim = list.sims[["Baseline0"]]; str(clus[[1]])

  ## need to laod sim one by one
  clust.stats <- function(clus, sim){
    
    ## empty list
    ll <- list()
    
    ## loop over sims
    for (i in 1:length(clus) ) {
      load(sim[i])
      
      ## get ALL tips names and demes
      tip.names <- data.frame(
        "id" = daytree$tip.label,
        "stage" = sapply( sampleDemes, deme2stage ),
        "age" = sapply( sampleDemes, deme2age ),
        "risk" = sapply( sampleDemes, deme2risk ),
        stringsAsFactors = FALSE)
      rownames(tip.names) <- NULL
      
      ## get size of clusters
      freqClust <- #lapply(clus, function(x){
        as.data.frame(table(clus[[i]]$ClusterID),
                      stringsAsFactors = FALSE)
      #})
      ## get name instead of index in clus
      clus[[i]]$seq.label <- daytree$tip.label[ clus[[i]]$SequenceID ]
      # head(clus[[i]]); i = 1
      
      ## merge cluster number (without SequenceID)
      a <- merge(x = tip.names, y = clus[[i]][, 2:3],
                 by.x = "id", by.y = "seq.label",
                 all.x = TRUE, sort = FALSE)
      
      ## merge cluster size (with NA)
      b <- merge(x = a, y = freqClust, 
                 by.x = "ClusterID", by.y = "Var1", 
                 all.x = TRUE, sort = FALSE)
      
      #- size 1 if not into a cluster
      b$Freq[is.na(b$Freq)] <- 1
      
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
      ll[[i]] <- b
      names(ll)[i] <- names(clus[i])
      print(paste(names(clus[i]), i))
    }
    return(ll)
  }
  ### --- end function clust.stats --- ###
 
if(0){
  
#  clus = cl_Baseline0[[1]]; i=2; sim = list.sims[["Baseline0"]]; str(clus[[1]])
  
 test <-  clust.stats(clus = cl_Baseline0[[1]], sim = list.sims[["Baseline0"]] )
 
  system.file(
  l_Baseline0 <- lapply(cl_Baseline0, 
                        function(x) clust.stats(clus = x, sim = list.sims[["Baseline0"]] ))
  )
  names(cl_Baseline0)
  
  l_EqualStage0 <- lapply(cl_EqualStage0, 
                        function(x) clust.stats(clus = x, sim = list.sims[["EqualStage0"]] ))
  
#   names(l_Baseline0)
#   names( l_Baseline0[[1]][[1]] )
#   length(l_Baseline0[[1]])
#   test <- l_Baseline0[[1]][[2]]
##   test <- b
##   lapply(test, function(x) aggregate(x$size, by = list("stage" = x$stage), mean ))
#   aggregate(test$size, by = list("risk" = test$risk), mean )
#   aggregate(test$size, by = list("age" = test$age), mean )
 
  ####---- saved listUKclus ----
  # saveRDS(l_Baseline0, file = "data/sim_ucsd_results/list.sim.ucsd.Baseline0.rds" )
  # saveRDS(l_EqualStage0, file = "data/sim_ucsd_results/list.sim.ucsd.EqualStage0.rds" )

  l_Baseline0 <- readRDS(file = "data/sim_ucsd_results/list.sim.ucsd.Baseline0.rds" )
  l_EqualStage0 <- readRDS(file = "data/sim_ucsd_results/list.sim.ucsd.EqualStage0.rds" )
}  



