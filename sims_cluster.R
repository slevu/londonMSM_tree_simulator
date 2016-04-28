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
for (.s in 1:length(scenario)){
  list.sims[[.s]] <- list.files('RData', full.names = TRUE, 
             path = paste('data/simulations/model0-simulate', 
                          scenario[.s], sep = '') )
  ## names of sims within scenario
  names(list.sims[[.s]]) <- lapply(list.sims[[.s]], function(x){
    substr(x, regexec(".RData", x)[[1]][1] - 5, regexec(".RData", x)[[1]][1] -1)
  })
  ## names of scenario
  names(list.sims)[.s] <- scenario[.s]
}
str(list.sims) # head(names(list.sims[[2]]))

####---- list.dist ----
#### load list of dist from sims
if(FALSE){
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
                              regexec("_ucsd_hivclust_output",
                                      filenames[n])[[1]][1] -1)
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
####---- clus.stat ----
if(1){
  ###- check if sim.name = clus.name
  # sim = list.sims[["Baseline0"]] ; clus = cl_Baseline0; # i = 1
clus.stat <- function(clus, sim){
  ##- empty list of n thresholds
  ll <- vector("list", length(clus))
  ##- loop threshold
  for (thr in 1:length(clus)){
    ##- loop sims
    for (i in 1:length(clus[[thr]]) ){
      ##- check same sim and clus
      if (!identical(names(sim[i]), names(clus[[thr]][i]) )) {
          stop(paste("not the same sims", names(sim[i]), names(clus[[thr]][i]) ))
      } else {
        print(paste( names(clus)[thr], i))
        
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
        freqClust <- as.data.frame(table(clus[[thr]][[i]]$ClusterID),
                        stringsAsFactors = FALSE)
        
        ## get name instead of index in clus
        clus[[thr]][[i]]$seq.label <- daytree$tip.label[ clus[[thr]][[i]]$SequenceID ]
        
        ## merge cluster number (without SequenceID)
        a <- merge(x = tip.names, y = clus[[thr]][[i]][, 2:3],
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
        ll[[thr]][[i]] <- b
        names(ll[[thr]])[i] <- names(clus[[thr]][i])
      }
    } 
    names(ll)[thr] <- names(clus)[[thr]]
  }
  return(ll)
}
####---- fin clus.stat ----####
}

system.time(
l_Baseline0 <- clus.stat(clus = cl_Baseline0, 
                         sim = list.sims[["Baseline0"]])
) # 459
system.time(
l_EqualStage0 <- clus.stat(clus = cl_EqualStage0,
                           sim = list.sims[["EqualStage0"]])
) 
#  [1] "0.05 85"
#  Error in clus.stat(clus = cl_EqualStage0, sim = list.sims[["EqualStage0"]]) : 
#  not the same sims 24308 24339 


# lapply(l_Baseline0, function(x) lapply(x, function(df) aggregate(df$size, by = list("stage" = df$stage), mean )))
  
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




