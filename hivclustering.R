rm(list=ls())
##----------------------------##
####---- function
####---- HIV clustering ----
####---- for simulated tree 
##----------------------------##
# cd /Documents/softwares:
# hivnetworkcsv -i dist_pour_hivclust.csv -c outhivclust.csv -t 0.04 -f plain

getwd()
# setwd("../")
# getwd()

test <- read.csv(file = "~/Documents/phylo-uk/source/subUKogC_noDRM_151202_ucsdTN93.csv")
dim(test)/12164 *12164 /2
head(test)

dsimtree <- readRDS("data/simtree_dist.rds")
## normalize
d <- round(dsimtree / (max(dsimtree) - min(dsimtree)), 4)
rm(dsimtree)

m <- as.matrix(d)
dim(m)

lt <- m[lower.tri(m, diag = FALSE)]
dim(lt)
library(reshape2)
system.time(
  el <- melt(lt)
)
head(el)
## get rid of distance > k = mean
k <- round(mean(d),1)
dd <- as.dist(d[d < k])
str(dd)
vect <- d < k
sum(vect)
as.dist(yoursymmetricmatrix[k, k])

rm(d)
str(dd)
m <- as.matrix(d)
dim(m)
mm <- m[12160:12164, 12160:12164]
dd <- as.dist(mm)

 require(igraph)
 g <- graph.adjacency(mm, mode = "lower", diag = FALSE)
 plot(g)
 get.edgelist(g)
#     ?graph.adjacency
melt(mm)

mm[lower.tri(mm, diag = FALSE)]
a <- m[ mm<0.3 ]
str(a)
ucsd_hivclust <- function(d){
  if ( class(d) == "dist" ){
    ## make a csv edge list (ID1, ID2, distance)
     ## convert matrix to edge list
     ## essai igraph
#     require(igraph)
#     g <- graph(d)
#     edge <- get.edgelist(g)
#     ?graph.adjacency
#     
      m <- as.matrix(d)
      edges <- matrix(0, nrow = nrow(m)*nrow(m), ncol = 3)
      system.time(
       for (i in 1:nrow(m)) {
        for (j in 1:ncol(m)) {
         edges[(i-1)*nrow(m)+j , ] <- c(rownames(m)[i], 
                         rownames(m)[j], 
                         m[i,j])
      }
         # progress
         print(i)
         # update GUI console
         flush.console() }
      )
      warnings()
      ## 37s for 3 rows = 6s by row, 6*12164 = 72984 / 60 /60 = 20h !!!
      ## 
    colnames(edges) <- c('ID1', 'ID2', 'distance')
    df <- as.data.frame(edges, stringsAsFactors = F)
    df$distance <- as.numeric(df$distance)
    
    } else stop("d is not a distance matrix")
  
#     head(edges)
#     str(df)
#     
    ## write csv without rownames !
    inputCSV <- paste(tempdir(), "/input.csv", sep = "")
    write.csv(df, file = inputCSV, row.names = FALSE )

  ## threshold
  thr <- 0.9 # or for (thr in c(0.5)){}
  
  ## output
  outputCSV <- paste(tempdir(), "/input_", thr,".csv", sep = '')
  ## parms
  parms <- paste("-t", thr, "-f plain")
  ## full path needed 
  exec <- "~/Documents/softwares/hivclustering/scripts/hivnetworkcsv"
  ## command
  cmd_hivclustering <- paste(exec, " -i", inputCSV, "-c", outputCSV, parms )
  cmd_hivclustering
  
  system(cmd_hivclustering)
  
  read.csv(outputCSV)
}

# thr <- 0.5
for (thr in c(0.5)){
  
  inputCSV <- paste(getwd(), "/source/","subUKogC_noDRM_151202_ucsdTN93.csv", sep = '')
  outputCSV <- paste(getwd(), "/data/","subUKogC_noDRM_151202_UCSDclust_", thr,".csv", sep = '')
  
  parms <- paste("-t", thr/100, "-f plain")
  
  ## full path needed 
  exec <- "~/Documents/softwares/hivclustering/scripts/hivnetworkcsv"
  ## command
  cmd_hivclustering <- paste(exec, " -i", inputCSV, "-c", outputCSV, parms )
  cmd_hivclustering
  
  system(cmd_hivclustering)
  
}