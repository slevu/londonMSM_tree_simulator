# rm(list=ls())
####---- include ----
detail_knitr <- TRUE
#source("functions.R")
startover <- FALSE

####---- lib ----
#library(ape)

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

#  for (.s in 1:length(scenario)){
#   ## names of sims within scenario
#   names(list.sims[[.s]]) <- lapply(list.sims[[.s]], function(x){
#     regmatches(x, regexpr("[0-9]{2,}", x)) # numerical string of length >= 2
#           })
# }
# str(list.sims); str(list.dist) # head(names(list.sims[[2]])) 
# load(list.sims[[2]][[1]]); head(sampleDemes); head(cd4s)

####---- ucsd clustering ----####
# ucsd_hivclust
if (startover == TRUE){
  thresholds  <-  c("0.001", "0.05") # c("0.015", "0.005") # c(0.005, 0.015, 0.02, 0.05) # c(0.005, 0.01, 0.02, 0.05, 0.1) 
  tmax <- max(thresholds) # limit of distance considered
  
## function: input list of dist filenames, output csv of clusters
#- debug: ldist = distBaseline0FNS; i = 2
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
                          out = paste(path.results, folder.sim, sep = '/' ) )
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
  thresholds <- c("0.001", "0.005", "0.015", "0.05") 
  
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
 # names(cl_EqualStage0[[1]]) 
  # str(cl_Baseline0[[1]][1])
  
# saveRDS(cl_Baseline0, file = paste(path.results, 'list.hivclust.sim.Baseline0.rds', sep = '/') )
# saveRDS(cl_EqualStage0, file = paste(path.results, 'list.hivclust.sim.EqualStage0.rds', sep = '/') )
}

####---- load list.hivclust ----

  cl_Baseline0 <- readRDS( file = paste(path.results, 'list.hivclust.sim.Baseline0.rds', sep = '/') )
  cl_EqualStage0 <- readRDS( file = paste(path.results, 'list.hivclust.sim.EqualStage0.rds', sep = '/') )

if (startover == TRUE){
  ##- Keep first 100 cluster results corresponding to a simulation
  ##- reference for names of sims
  ref_e <- names(list.sims[["EqualStage0"]])
  ref_b <- names(list.sims[["Baseline0"]])
  if(length(ref_e) > 100){
     cl_EqualStage0 <- lapply(cl_EqualStage0, 
                  function(a){
                    a[names(a) %in% ref_e][1:100]
                         })
  }
  if (length(ref_b) > 100){
    cl_Baseline0 <- lapply(cl_Baseline0, 
                  function(a){
                    a[names(a) %in% ref_b][1:100]
                         })
  }
  rm(ref_b, ref_e)
  # identical(names(cl_EqualStage0[[1]]), names(cl_EqualStage0[[2]]) )
  # sims at threshold 4 EqualStage0 are different
}

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
  # sim = list.sims[["Baseline0"]] ; clus = cl_Baseline0; i = 100 ; thr = 4
clus.stat <- function(clus, sim){
  
  ##- Start loop
  ##- empty list of n thresholds
  ll <- vector("list", length(clus))
  ##- loop threshold
  for (thr in 1:length(clus)){
    ##- loop sims
    for (i in 1:length(clus[[thr]]) ){
      ##- check same sim and clus
        print(paste( 
          names(clus[[thr]])[i], # name cluster
          names(clus)[thr], # threshold
          i # num sim
          ))
      cl <- clus[[thr]][[i]]
      ##- find sim that match first clus
      ##- watch for embedded number
      ##- all deme states are the same within one model 
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
        
        ## merge cluster number (without SequenceID nor covariates)
        a <- merge(x = tip.states, y = cl[, 2:3],
                   by.x = "id", by.y = "id",
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
    names(ll)[thr] <- names(clus)[[thr]]
  }
  return(ll)
}
####---- fin clus.stat ----####
}

if (startover == TRUE){
  system.time(
    l_Baseline0 <- clus.stat(clus = cl_Baseline0, 
                             sim = list.sims[["Baseline0"]])
  ) # 469
  system.time(
    l_EqualStage0 <- clus.stat(clus = cl_EqualStage0,
                               sim = list.sims[["EqualStage0"]])
  ) # 463
  
}

  ####---- saved listUKclus ----
 # saveRDS(l_Baseline0, file = "data/sim_ucsd_results2/list.sim.ucsd.Baseline0.rds" )
 # saveRDS(l_EqualStage0, file = "data/sim_ucsd_results2/list.sim.ucsd.EqualStage0.rds" )

  l_Baseline0 <- readRDS(file = "data/sim_ucsd_results2/list.sim.ucsd.Baseline0.rds" )
  l_EqualStage0 <- readRDS(file = "data/sim_ucsd_results2/list.sim.ucsd.EqualStage0.rds" )

###- check
  #   names(l_Baseline0)
  #   head( l_Baseline0[[1]][[1]] )
  #   length(l_Baseline0[[1]])
# 
#   system.time(
#   z <- lapply(l_Baseline0, function(x){
#     sapply(x, function(m) {
#     # faster than `merge(df, s, all.x = TRUE)`
#     #- add covariates
#      # m <-  cbind(df, s[ match(df$id, s$id), -1 ])
#     #- mean y by x
#       agg <-  tapply(m$size, m$stage, mean )
#       return(agg)
#       }) 
#     })
#   )
  

####---- add degrees and neighborhood size ----

  ##- loop to merge all W to all clusters
#   clus <- l_Baseline0
#   sim <- list.sims[["Baseline0"]]
#   thr <- 1
#   i <- 15
  ####---- add.w ----
  add.w <- function(clus, sim, dist){
    ##- empty list of n thresholds
    ll <- vector("list", length(clus))
    ##- loop threshold
    for (thr in 1:length(clus)){
      ##- loop sims
      for (i in 1:length(clus[[thr]]) ){
        ##- check same sim and clus
        #           if (!identical(names(clus[[thr]])[i], # name cluster
        #                          names(sim[ names(clus[[thr]])[i] ]))  # name sim
        #               ) stop("not the same sim results")
        
        ##- monitor loop
        print(paste(names(clus)[thr], # threshold
                    i )) # num sim
        
        ##- load sim with W matrix of infect probs
        load(sim[ names(clus[[thr]])[i] ])
        
        ##- outdegree
        out0 <- aggregate(x = list(outdegree = W$infectorProbability), 
                          by = list(patient = W$donor), FUN = sum)
        ##- indegree
        in0 <- aggregate(x = list(indegree = W$infectorProbability),
                         by = list(patient = W$recip), FUN = sum)
        
        ##- merge out and in degrees
        b <- merge(out0,in0)
        
        ###- Add neighborhood size
        load(dist[[i]])
        dd <- as.data.frame(t(D))
        
        ##- calculate neighborhood size by tip label
          ##- number of neighbour|threshold
        .t  <- tapply(dd$V3, dd$V1, function(x) sum(x < thr))
          ##- convert tip ID in tip label
        names(.t) <- lapply(names(.t), function(x) daytree$tip.label[ as.numeric(x) ])
        nbrhood <- data.frame("id" = names(test), "nbhsize" = test, stringsAsFactors = FALSE)
        
        ##- load cluster assignement
        a <- clus[[thr]][[i]]
        
        ##- merge all
        c <- merge(a, b, by.x = "id", by.y = "patient", all.x = TRUE )
        d <- merge(c, nbrhood, by = "id")
        
        ##- names
        ll[[thr]][[i]] <- d
        names(ll[[thr]])[i] <- names(clus[[thr]][i])
      } 
      names(ll)[thr] <- names(clus)[[thr]]
    }
    return(ll)
  }
  ####---- fin add.w ----####
  
  
  if (startover = TRUE){
    system.time(
      cw_Baseline0 <- add.w(clus = l_Baseline0, 
                            sim = list.sims[["Baseline0"]],
                            dist = list.dist[["Baseline0"]] )
    ) # 464
    system.time(
      cw_EqualStage0 <- add.w(clus = l_EqualStage0,
                              sim = list.sims[["EqualStage0"]],
                              dist = list.dist[["EqualStage0"]])
    ) # 404
  }
  
  
  ####---- saved listUKclus ----
  # saveRDS(cw_Baseline0, file = "data/sim_ucsd_results2/list.sim.clus-outdeg.Baseline0.rds" )
  # saveRDS(cw_EqualStage0, file = "data/sim_ucsd_results2/list.sim.clus-outdeg.EqualStage0.rds" )
  
  cw_Baseline0 <- readRDS(file = "data/sim_ucsd_results2/list.sim.clus-outdeg.Baseline0.rds" )
  cw_EqualStage0 <- readRDS(file = "data/sim_ucsd_results2/list.sim.clus-outdeg.EqualStage0.rds" )
 
  
  