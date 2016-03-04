##- todo :
## use system2
## set stderr = TRUE to capture warnings

##----------------------------##
####---- function
####---- HIV clustering ----
####---- for a tree 
##----------------------------##

###--- could scale the time based distance by consensus mutation rate to obtain 'standard' substitution per site distances
###--- first, construct edge list with distance as input for ucsd software (d = distance matrix)
###--- second, loop through distance thresholds ( = vector of threshold based on quantiles)

ucsd_hivclust <- function(d, quant = c(5e-4, 1e-3, 1e-2, 1e-1, 0.25, 0.5) ){
  ## get var.name for output path
  var.name <- substitute(d)
  
  ## normalize (if not normalized)
  if(max(d) > 1) {
    d <- round(d / max(d), 4)
  }
  
  ## as matrix
  m <- as.matrix(d)

  ## keep only the lower triangle by 
  ## filling upper with NA
  m[upper.tri(m, diag=TRUE)] <- NA
  
  ## create edge list if not there
  require(reshape2)
  out_edge_list <- paste("data/", var.name, "_el.rds", sep = '')
  
  if (file.exists(out_edge_list)){
    el <- readRDS(out_edge_list)
    } else {
  # system.time(
      el <- melt(m, na.rm = TRUE)
  # ) # the na.rm removal takes most of time
      colnames(el) <- c('ID1', 'ID2', 'distance')
      saveRDS(el, file = out_edge_list)
  }

  # choose threshold based on quantiles
    qt <- quantile(el$distance, 
                   probs = quant )
    
  ## get rid of distance > larger quantile (= median)
  k <- round(qt[length(qt)], 2)
  subel <- el[el$distance < k,]
  # rm(el, m, subel, d)

  ## write csv without rownames and get input path
  inputCSV <- paste(tempdir(), "/input.csv", sep = "")
  write.csv(subel, file = inputCSV, row.names = FALSE )

  ## full path needed 
  exec <- '~/Documents/softwares/hivclustering/scripts/hivnetworkcsv'
  
####---- loop threshold (first 3 qt = 0.05, 0.1, 1, 10)

  thr <- round(qt[1:3], 2)
  ## empty results
 cmd <- vector( mode= "character" )
 warn <- list()
 
  for ( t in thr ){
    print(paste("threshold =", t))
  ## output
  outputCSV <- paste(var.name, "_ucsd_hivclust_output_", 
                     t, ".csv", sep = '')
  
  ## parms
  parms <- paste("-t", t, "-f plain")
 
  ## command
  cmd_hivclustering <- paste(exec, "-i", inputCSV, "-c", outputCSV, parms )
  
  print(cmd_hivclustering)
  
  ## too long, the fitting of degree
  ## distribution takes most time
  ## => just issue command
  ## to run on terminal
  if ( !file.exists(outputCSV) ){
   stderr <-  system(
     paste(cmd_hivclustering, "2>&1"),
           intern = TRUE)
  }
  # save commands and 'stderr' warnings
  cmd <- c(cmd, cmd_hivclustering)
  warn <- c(warn, c(t, stderr))
  }
  
###--- bin table in one list
  cl <- list()
  for(i in 1:length(thr)){
    ## add table i
    CSV <- paste(var.name, 
                 "_ucsd_hivclust_output_", thr[i],
                 ".csv", sep = '')
    if(file.exists(CSV)){
       cl[[i]] <- read.csv(CSV)
       # name threshold
       names(cl)[i] <- thr[i]
       }
  }
 
 return(list(qt, cmd, cl))
}
