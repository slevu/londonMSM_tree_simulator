rm(list=ls())
##----------------------------##
####---- HIV clustering ----####
##----------------------------##
# cd /Documents/softwares:
# hivnetworkcsv -i dist_pour_hivclust.csv -c outhivclust.csv -t 0.04 -f plain

getwd()
# setwd("../")
# getwd()

list.files("../softwares/")  
  inputCSV <- "../softwares/dist_pour_hivclust.csv"
  outputCSV <- "../softwares/outhivclust_system2.csv"
  
  parms <- paste("-t0.04 -f plain")
  
  ## full path needed 
  exec <- "~/Documents/softwares/hivclustering/scripts/hivnetworkcsv"
#   exec <- shQuote("hivnetworkcsv")
#   exec <- "ls"

  ## command
  cmd_hivclustering <- paste(exec, " -i", inputCSV, "-c", outputCSV, parms )
  cmd_hivclustering
  
  ## args
  args <- paste(" -i", inputCSV, "-c", outputCSV, parms )
  
  system(cmd_hivclustering)

  ## capture stderr
  stderr <- system(paste( cmd_hivclustering, "2>&1"), intern = TRUE)
  
## does not work 
#  a <-  system2(command = "hivnetworkcsv", args = args, stderr = T, stdout = T)
#   ?system2
