###--- not used, now in TreeToEdgeList adn EdgeListToCluster
### rm(list=ls())
##----------------------------##
####---- function
####---- HIV clustering ----
####---- for simulated tree 
##----------------------------##

getwd()

###--- could scale the time based distance by consensus mutation rate to obtain 'standard' substitution per site distances
###--- first, construct edge list with distance as input for ucsd software
###--- second, loop through distance thresholds 

if(FALSE){
  
  dsimtree <- readRDS("data/simtree_dist.rds")
  ## normalize
  dsimtree <- round(dsimtree / 
                      (max(dsimtree) - min(dsimtree)), 4)
  str(dsimtree)

  ## as matrix
  m <- as.matrix(dsimtree)
  dim(m)
  # m <- m[1:6,1:6]

  ## keep only the lower triangle by 
  ## filling upper with NA
  m[upper.tri(m, diag=TRUE)] <- NA

  library(reshape2)
  system.time(
    el <- melt(m, na.rm = TRUE)
  ) # the na.rm removal takes most of time
  # head(el)
  # str(el)
  colnames(el) <- c('ID1', 'ID2', 'distance')
  ## omit the NA values
  # el <- na.omit(el)

  ## get rid of distance > k = mean
  k <- round(mean(el$distance),1)
  subel <- el[el$distance < k,]
  head(subel)

  ## write csv without rownames !
  inputCSV <- paste(tempdir(), "/input.csv", sep = "")
  write.csv(subel, file = inputCSV, row.names = FALSE )
  
  rm(m, dsimtree, subel, el)
}

####---- loop threshold
for (thr in (2:6)/10){
  inputCSV <- paste(tempdir(), "/input.csv", sep = "")
  ## output
  outputCSV <- paste("ucsd_hivclustering_output_", thr,".csv", sep = '')
  ## parms
  parms <- paste("-t", thr, "-f plain")
  ## full path needed 
  exec <- "~/Documents/softwares/hivclustering/scripts/hivnetworkcsv"
  ## command
  cmd_hivclustering <- paste(exec, "-i", inputCSV, "-c", outputCSV, parms )
  cmd_hivclustering
  
  system(cmd_hivclustering)
}


###--- Analyse cluster
###t <- 3
    cl <- data.frame()
    for(t in 2:6/10){
      cl <- rbind(cl,
              cbind(t,
               read.csv(
                 paste("ucsd_hivclustering_output_", t, ".csv", sep = '')
    )))
    }
  str(cl)
  
  ##- Calculate size(=Freq) of each cluster across different threshold
  freqClust <- apply(cluster, 2, function(x) as.data.frame(table(x))) # list
  
  ##- number of different clusters by threshold
  sapply(freqClust, function(x) dim(x)[1])
  ##- cluster size
  sapply(freqClust, function(x) summary(x$Freq))
  ##- percentiles
  sapply(freqClust, function(x) round(quantile(x$Freq, probs = c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.95, 0.99, 1))))
  
  
  
  simfreqClust <- tapply(cl$ClusterID, cl$t, 
                         function(x) as.data.frame(table(x),
                                                   stringsAsFactors = FALSE))
  # str(freqClust)
  # head(freqClust[[5]])
  
  ##- number of different clusters by threshold
  sapply(simfreqClust, function(x) dim(x)[1])
  ##- cluster size
  sapply(simfreqClust, function(x) summary(x$Freq))

