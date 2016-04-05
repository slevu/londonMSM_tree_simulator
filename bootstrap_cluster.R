# rm(list=ls())

## tree path
path.trees <- "../HPC/phylo-msm-uk/in_HPC_20160404/"
a <- list.files(path.trees)
l <- a[grep("result", a)] # list of tree results
## get rid of tree 000
l <- l[-1]

## filenames number
num <- sprintf("%03d", seq(1, 100))

## outgroup
og <- c( paste("Ref", 1:6, sep = ''), "HXB2" )

###--- Clustering ---###
source("functions.R")

##- list of trees (minus outgroup tips)
first <- 5
last <- 6
## empty list 
list.of.trees <- vector("list", last - first + 1)
## index j of list elements 
j <- 1
for (i in first:last){
  t <- read.tree(file = paste(path.trees, l[i], sep = ''))
  list.of.trees[[j]] <- drop.tip(t, og)
  names(list.of.trees)[j] <- paste("bs_tree_", num[i], sep="")
  j <- j + 1
}
str(list.of.trees)

##- Save edge lists
test <-  function(i){
          TreeToEdgeList(t = list.of.trees[[i]],
                         name.output = names(list.of.trees)[[i]],
                             rate = 1,
                             output = "data/bootstrap/", 
                             stats = FALSE,
                             plot = FALSE)
}

system.time(
 lapply(seq_along(list.of.trees), test) 
)


# lapply(seq_along(list.of.trees), function(i) names(list.of.trees)[[i]])

#######

for (i in first:last) {
  filename <- paste(path.trees, name.tree, num[i], sep = '')
t <- read.tree(filename)

## name tree
ntree <- paste("bs_tree_", num[i], sep="")
assign( ntree, drop.tip(t, og) )
print(paste("Tree object is named", ntree ))
a <- get(ntree)
class(deparse(substitute(ntree)) )
ls()

system.time(
  path001 <-  TreeToEdgeList(t = a,
                             rate = 1,
                             output = "data/bootstrap/", 
                             stats = FALSE,
                             plot = FALSE)
)
ucsd_hivclust
?assign
?parse
parse(ntree)
?get
}
library(parallel)
detectCores()
?mclapply
