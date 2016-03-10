##- todo :
## Kill external program: http://goo.gl/AQVYRP

##----------------------------##
####--- Apply UCSD
####--- HIV clustering
####--- to a dataframe of 
####--- pairwise distances 
##----------------------------##

###- input: path to RDS file of edge list
###--- loops through distance thresholds ( = vector of threshold based on quantiles) 

ucsd_hivclust <- function(path.el, thr = NA, quant = c(5e-5, 1e-4, 5e-4, 1e-3, 1e-2, 1e-1) ){
  
  ## read RDS file
  el <- readRDS(file = path.el)
  
  ## get var.name for output path (Regex !)
  var.name <- paste("d", substr(path.el,
                    regexpr("\\/[^\\/]*$", path.el)[[1]][1] +1,
                    regexpr("\\_el.rds", path.el)[[1]][1] - 1 ),
                    sep = '')
  
  ## Either thr given or based on quantiles
  if(is.na(thr[1])){
    
      # choose threshold based on quantiles (first 4 for now)
      qt <- quantile(el$distance, 
                   probs = quant )
      thr <- round(qt[1:4], 2)
    
     ## get rid of distance > larger quantile
#      k <- round(qt[length(qt)], 2)
#      subel <- el[el$distance < k,]
     
  } else { 
    thr <-  thr
    qt <- "Not used"
    }
  
  ## write csv without rownames and get input path
  inputCSV <- paste(tempdir(), "/input.csv", sep = "")
  write.csv(el, file = inputCSV, row.names = FALSE )

  ## full path needed 
  exec <- '~/Documents/softwares/hivclustering/scripts/hivnetworkcsv'
  
####---- loop threshold (first 4 qt = 0.05, 0.1, 1, 10) 

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
 
 return(list(qt = qt, cmd = cmd, warn = warn, cl = cl))
}
