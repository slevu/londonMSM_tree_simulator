##- load data simulations
if( any(grep("MacBook", Sys.info())) ){
  path.results <- '../Box Sync/HPC/simulations'
} else {
  path.results <- '../Box/HPC/simulations' # imac
}

{
  cw_Baseline0 <- readRDS(file = paste(path.results, 'model1-sim_ucsd/list.sim.clus-outdeg.Baseline0.rds', sep = '/') )
}