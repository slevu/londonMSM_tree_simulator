##- calculate genetic distances from simulated year-tree
##- keep name of simulation
setwd(getwd())

####---- scenario ----
scenario <- c("EqualStage0", "Baseline0")

####--- list.trees ----
list.sims <- vector("list", length(scenario))
for (.s in 1:length(scenario)){
  list.sims[[.s]] <- list.files(
    'RData', full.names = TRUE, 
    path = paste(getwd(), '/data/simulations2/model0-simulate', 
    scenario[.s], sep = '') )
  
  ## names of sims within scenario
  names(list.sims[[.s]]) <- lapply(list.sims[[.s]], function(x){
    regmatches(x, regexpr("[0-9]{2,}", x)) # numerical string of length >= 2
  })
  ## names of scenario
  names(list.sims)[.s] <- scenario[.s]
}
str(list.sims) # head(names(list.sims[[2]]))

source("tree2genDistance.R")

##- loop over 2 scenarii * 100 sims
for (i in 1:length(scenario) ){
  for (j in 1:length(list.sims[[i]]) ){
    load(list.sims[[i]][j])
    name <- names(list.sims[[i]][j])
    path.dist <- paste(getwd(), '/data/simulations2/model0-simulate', 
                       scenario[i], '-distances/', name, '.RData', sep = '')
    D <-  .tree2genDistance(bdt, mu = .0015, seqlength = 1e3, MH=20 )
    save(D, file = path.dist)
    print(paste(i,j))
  }
}
