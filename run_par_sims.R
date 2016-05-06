
library(doParallel)
registerDoParallel(cores = 4)

system.time(
foreach(i = 1:4) %dopar% {
  source('Erik_code/model0-simulateBaseline0.R')
}
)
# getDoParWorkers()
# getDoParName()
?detectCores()

##- test nested
registerDoParallel(cores = 5)
foreach(k = 1:3, .combine = cbind) %:% 
foreach(i = 1:5, .combine = c) %dopar% {
  Sys.getpid() # Sys.info()[['nodename']]
}

?"foreach"
