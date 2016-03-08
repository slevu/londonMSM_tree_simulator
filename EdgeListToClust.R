##- todo :
## Kill external program: http://goo.gl/AQVYRP

##----------------------------##
####--- Apply UCSD
####--- HIV clustering
####--- to a dataframe of 
####--- pairwise distances 
##----------------------------##

###--- loops through distance thresholds ( = vector of threshold based on quantiles) 

ucsd_hivclust <- function(el, quant = c(1e-4, 5e-4, 1e-3, 1e-2, 1e-1) ){
  ## get var.name for output path
  var.name <- "dsimtree"
  
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
  
####---- loop threshold (first 3 qt = 0.05, 0.1, 1, 10) t <- thr[2]

  thr <- round(qt[1:4], 2)
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
