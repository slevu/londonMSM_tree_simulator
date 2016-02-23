### source attribution on simulated trees ###
### ------------------------------------- ###
### First assume simulated tree is dated 
### as would the real tree after LSD (later try BEAST) 
### Second input dates and tree on Rcolgem
### Third analyse out-degrees
rm(list=ls())
####---- lib ----
library(ape)
library(ggplot2)

####---- load sim ----
## Load newest Rdata
l <- list.files(pattern="*.Rdata") # list.files(pattern="Rdata$") list.files(pattern="out")
load(l[length(l)])
ls()

## look tree
class(tree)
str(tree)

is.rooted(tree)

head(tree$sampleTimes)
tail(tree$sampleTimes)

##- see W0_rerunSLV.R for re-run of LSD on real data
STFN <- '../phylo-uk/data/ExaML_result.subUKogC_noDRM.reroot_dropOG.dates'
dates <-( read.table(STFN, skip=1, colClasses=c('character', 'numeric') ) )
head(dates)
sampleTimes <- scan( file = 'sampleTimes' )
head(sampleTimes)
