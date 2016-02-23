rm(list=ls())
##----------------------------##
####---- function
####---- HIV clustering ----
####---- for simulated tree 
##----------------------------##
# cd /Documents/softwares:
# hivnetworkcsv -i dist_pour_hivclust.csv -c outhivclust.csv -t 0.04 -f plain

getwd()
# setwd("../")
# getwd()

# test <- read.csv(file = "~/Documents/phylo-uk/source/subUKogC_noDRM_151202_ucsdTN93.csv")
# dim(test)/12164 *12164 /2
# head(test)
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
  )
  # head(el)
  str(el)
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

for (thr in (3:6)/10){
  inputCSV <- paste(tempdir(), "/input.csv", sep = "")
  ## output
  outputCSV <- paste("ucsd_hivclustering_output_", thr,".csv", sep = '')
  ## parms
  parms <- paste("-t", thr, "-f plain")
  ## full path needed 
  exec <- "~/Documents/softwares/hivclustering/scripts/hivnetworkcsv"
  ## command
  cmd_hivclustering <- paste(exec, " -i", inputCSV, "-c", outputCSV, parms )
  cmd_hivclustering
  
  system(cmd_hivclustering)
}
    
  cl <- read.csv(outputCSV)

head(cl)

##### surplus #####
# thr <- 0.5
# for (thr in c(0.5)){
#   
#   inputCSV <- paste(getwd(), "/source/","subUKogC_noDRM_151202_ucsdTN93.csv", sep = '')
#   outputCSV <- paste(getwd(), "/data/","subUKogC_noDRM_151202_UCSDclust_", thr,".csv", sep = '')
#   
#   parms <- paste("-t", thr/100, "-f plain")
#   
#   ## full path needed 
#   exec <- "~/Documents/softwares/hivclustering/scripts/hivnetworkcsv"
#   ## command
#   cmd_hivclustering <- paste(exec, " -i", inputCSV, "-c", outputCSV, parms )
#   cmd_hivclustering
#   
#   system(cmd_hivclustering)
#   
# }