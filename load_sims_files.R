##---- paths ----
if( any(grep("stephane", Sys.info())) ){
  path.sims <- '../Box Sync/HPC/simulations/model1-sim' #'data/simulations2/model0-simulate'
  path.results <- '../Box Sync/HPC/simulations/model1-sim_ucsd'
} else {
  path.sims <- '../Box/HPC/simulations/model1-sim'
  path.results <- '../Box/HPC/simulations/model1-sim_ucsd' # imac
}

##- scenario
scenario <- c("Baseline0", "EqualStage0")
scenario <- setNames(scenario, scenario) # useful to name list in lapply

##- list of sims files and distances files and name the vector by sim
list.sims <- lapply(scenario, function(x){
  v <- list.files('RData', full.names = TRUE, 
             path = paste(path.sims, x, sep = '') )
  names(v) <- sapply(v, function(i) sub(".RData", "", basename(i)))
  return(v)
})

