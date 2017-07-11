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

##- list of sims files and distances files
list.sims <- lapply(scenario, function(x){
  list.files('RData', full.names = TRUE, 
             path = paste(path.sims, x, sep = '') )
})
# str(list.sims)