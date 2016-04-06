### pipeline for UK bootstrapped trees ###
##1. import
##2. LSD
##3. clustering
##4. SA

rm(list=ls())
###--- LSD ---###

library(ape)
getwd()

##-- names
 ## tree path
 path.trees <- "../HPC/phylo-msm-uk/in_HPC_20160404/"
 list.files(path.trees)
 ## generic part of tree name
 name.tree <- "ExaML_result.subUKogC_noDRM_151202.finaltree."
 ## dates
 inputDAT <-  "data/LSD/t000.dates"


##- get dates in calendar year
res <- read.csv( '../phylo-uk/source/resistance.csv'  )
origin <- as.Date('0000-01-01')
x <- setNames( paste(sep='', '15/', res$dbsample_my ), res$testindex )
head(x)
x <- setNames( as.numeric(as.Date( x, format = '%d/%m/%Y' ) - origin)/365 , res$testindex )
dates <- x[ t$tip.label ]
## "The first line must indicate the number of leaves"
write.table( cbind( names(dates), dates), file = inputDAT, col.names=c(as.character(length(dates)), '') , quote = F, row.names=F) 

##-- loop LSD cmd
## filenames
a <- seq(1, 100)
num <- sprintf("%03d", a)

## outgroup
og <- c( paste("Ref", 1:6, sep = ''), "HXB2" )

## parms
parms <- paste("-c -s 1200 -v -b 10")

## empty results

lsd_stderr <- list()

for (i in 1:length(num)){
  filename <- paste(path.trees, name.tree, num[i], sep = '')
  t <- read.tree(filename)
  t <- drop.tip(t, og)
  ## rooted newick tree without OG
  inputNWK <-  paste("data/LSD/t", num[i], ".nwk", sep = "")
  ## create newick trees
  write.tree(t, file = inputNWK)
  
  ## command
  cmd_lsd <- paste("lsd", " -i", inputNWK, "-d", inputDAT, parms )
  print(cmd_lsd)
  
      stderr <-  system(
      paste(cmd_lsd, "2>&1"),
      intern = TRUE)
  
  lsd_stderr[[i]] <- stderr
  # name threshold
  names(lsd_stderr)[i] <- num[i]
}
##- save stderr
# saveRDS(lsd_stderr, file = "data/LSD/lsd_stderr.rds")
# lsd_stderr <- readRDS(file = "data/LSD/lsd_stderr.rds")

##- extract rate and tMRCA
b <- lsd_stderr
rate <- tMRCA <- vector()
for (i in 1:100){
  string <- unlist( strsplit(b[[i]][3], "\t") )
  rate <- c(rate, as.numeric(string[2]) )
  tMRCA <- c(tMRCA, as.numeric(string[4]) )
}
hist(rate)
hist(tMRCA)


# lsd -c  -i ukdrdb/ExaML_trees/ExaML_result.subUKogC_noDRM.reroot_dropOG.nwk -d ukdrdb/ExaML_trees/ExaML_result.subUKogC_noDRM.reroot_dropOG.dates   -s 1200   -v   -b  10
# options: -c temporal constraints recommanded; -i tree file newick format; -d dates file; -s length of seq (default 1000); -v variances; -b (default 10); -r to estimate root position

  