##- MLD regression on simulations
##- find estimates of W_ij
##- define incident infection
##- find MLDs
##- correlates of being MLS vs not
getwd()
PATHW <- "./data/simulations2/model0-simulateBaseline0"
l <- list.files(PATHW, full.names = TRUE)

##- test on first
load(l[1])

####---- helper functions ----
deme2age <- function(deme){ as.numeric(
  substr(regmatches( deme , regexpr( '\\.age[0-9]', deme ) ), 5, 5) 
) }
deme2stage <- function(deme){as.numeric( 
  substr(regmatches( deme , regexpr( 'stage[0-9]', deme )), 6, 6)
) }
deme2risk <- function(deme){as.numeric( 
  substr(regmatches( deme , regexpr( 'riskLevel[0-9]', deme )), 10, 10)
) }

## get ALL tips names and demes
tip.states <- data.frame(
  "id" = daytree$tip.label,
  "stage" = sapply( sampleDemes, deme2stage ),
  "age" = sapply( sampleDemes, deme2age ),
  "risk" = sapply( sampleDemes, deme2risk ),
  stringsAsFactors = FALSE)
rownames(tip.states) <- NULL

head(tip.states)

##- select ID of stage 1
##- apply MLD