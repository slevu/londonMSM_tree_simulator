##-----------------------------------------##
####---- function converting timed tree
####---- in edge list of pairwise distances 
##-----------------------------------------##

###--- input: timed tree, substitution rate
###--- output: returns path for dataframe (ID1, ID2, distance)

###--- Scales the time based distance by consensus mutation rate to obtain 'standard' substitution per site distances

TreeToEdgeList <- function(t, rate = 4.3e-3/365, 
                           output = "data/", fig = "figure/" , 
                           stats = TRUE, plot = TRUE ){
  
  ###--- object names for output path
  tree.name <- substitute(t)
  dist.mat.name <- paste(tree.name, "_dist", sep = '')
  dist.mat.path <- paste(output, dist.mat.name,".rds", sep = '')
  plot.name <- paste(fig, dist.mat.name, ".png", sep ='')
                     
  ###--- get distances
  if (file.exists(dist.mat.path)){
    d <- readRDS( dist.mat.path )
  } else {
    d <- as.dist(cophenetic.phylo(t))
    print("Matrix of cophenetic distances was saved")
    saveRDS(d, file = dist.mat.path )
  }
  
  ###--- time to rate
  d <-  d * rate
 
  ###--- stats
  if(stats == TRUE){ print(summary(d)) }
  if(plot == TRUE){ 
    png(filename = plot.name, type="cairo",
        units="in", width=5, height=4, 
        pointsize=12, res=96)
    hist(d, breaks = 50, 
         xlab = "distance", ylab = "frequency",
         main = '', # Tree's distances in subst/site
         col = "grey")
    dev.off()
   }
 
  ###--- normalize (if not normalized)
  if(max(d) > 1) {
    d <- round(d / max(d), 4)
  }
  
  ###--- as matrix
  m <- as.matrix(d)

  ###--- keep only the lower triangle by 
  ## filling upper with NA
  m[upper.tri(m, diag=TRUE)] <- NA
  
  ###--- create edge list if not there
  require(reshape2)
  out_edge_list <- paste(output, tree.name, "_el.rds", sep = '')
  if (file.exists(out_edge_list)){
    el <- readRDS(out_edge_list)
    } else {
  # system.time(
      el <- melt(m, na.rm = TRUE)
  # ) # the na.rm removal takes most of time
      colnames(el) <- c('ID1', 'ID2', 'distance')
      saveRDS(el, file = out_edge_list)
      print("Edge list was saved")
    }
  return(out_edge_list)
}
  
