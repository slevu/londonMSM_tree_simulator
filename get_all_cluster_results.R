source("load_sims_files.R")
cw_Baseline0 <- readRDS(file = paste(path.results, 'list.sim.clus-outdeg.Baseline0.rds', sep = '/') )
sim.name <- "10222" ## do equal rate sim also
hclust <- lapply(cw_Baseline0, function(x) x[sim.name])

##- load mrca clusters
M_CLUST_LIST <- paste0(path.results, '/', 'mrca_clusters_list_', sim.name, '.rds')
mclust <- readRDS(file = M_CLUST_LIST)

list.clus.in.df <- function(lst) data.frame(SequenceID = as.integer(do.call(c, lst)), ClusterID = as.integer(rep(names(lst), sapply(lst, length))) )
##- add size to seq / cluster data.frame
get.cluster.size <- function(df){
  freqClust <- as.data.frame(table(df$ClusterID), stringsAsFactors = FALSE)
  newdf <- merge(x = df, y = freqClust, 
                 by.x = "ClusterID", by.y = "Var1", 
                 all.x = TRUE, sort = FALSE)
  return(newdf[order(newdf$SequenceID),])
}

mrca_clusters <- lapply(mclust, function(x) get.cluster.size(list.clus.in.df(x)))

# id <- mrca_clusters[[1]]$SequenceID
# head( as.character(id[order(id)]) )
# tail( as.character(id[order(id)]) )
# 
# id2 <- as.integer( test[[1]][[1]]$id )
# head( as.character(id2[order(id2)]) )
# tail( as.character(id2[order(id2)]) )

ordered_hclust <- lapply(hclust, function(x) x[[1]][order(as.integer(x[[1]]$id)), ] )


o2 <- lapply(names(ordered_hclust[-1]), function(n) {
  df <- ordered_hclust[[n]][, c("ClusterID", "size", "nbhsize", "binclus")]
  names(df) <- c(paste0('clus', n), 
                 paste0('size', n), 
                 paste0('nbhsize', n), 
                 paste0('binclus', n))
  return(df)
  }) # without SA

o3 <- c(list(ordered_hclust[[1]]), o2)

o4 <- c(list('id' = mrca_clusters[[1]]$SequenceID),
        lapply(names(mrca_clusters), function(n) {
  df <- mrca_clusters[[n]][, c("ClusterID", "Freq")]
  df$binclus <- ifelse(df$Freq < 2, 0, 1)
  names(df) <- c(paste0('clus', n), 
                 paste0('size', n),
                 paste0('binclus', n))
  return(df)
  })
)

df_hclust <- do.call(cbind, o3)
df_mclust <- do.call(cbind, o4)

str(df_hclust)
str(df_mclust)

df <- cbind(df_hclust, df_mclust[,-1])

library(ggplot2)
d <- ggplot(data= df, aes(size0.015, size5))
d + geom_hex()# stat_bin_hex() #geom_hex()

d <- ggplot(data= df[df$size0.015 > 1 & df$size5 > 1 ,], aes(size0.015, size5)) # without singleton
d + geom_hex()# stat_bin_hex() #geom_hex()

cor(df$size0.015, df$size5)

x <- as.data.frame( table(df[df$size0.015 > 1 & df$size5 > 1 , c("size0.015", "size5")]), stringsAsFactors = FALSE)

radius <- sqrt(x$Freq/pi)
symbols(x$size0.015, x$size5, circles = radius, inches = 0.25, fg = "white",
        bg = "red",
        xlab = paste("hivclustering cluster size -"),
        ylab = paste("mrca cluster size -"),
        main = "Sized by Freq")
##- make sure id are ordered
##- keep id in first list on both 

names(df)
##---- plots ----
##- boxplot with base R
##- strip of 3 plots without outliers

bp_base <- function(df, Y , X , lbl, tran = identity, cx = 1.5){
  
  par(mfrow=c(1, length(Y)),
      oma = c(2.5, .5, .5, 1), # b,r,t,l
      mar = c(2, 4.5, 1, 1),
      cex.axis = cx, # x
      cex.lab = cx) # y
  
  lapply(1:length(Y), function(y){
    if ( any(grepl('outdegree', Y[y])) ) {
      
      boxplot(df[, 'outdegree'] ~ df[, X], outline = FALSE )
      title(main='', ylab = 'Out-degree', xlab = '')
      
    } else {
      
      boxplot(df[, Y[y]] ~ df[, X], outline = FALSE)
      title(main='', ylab = Y[y], xlab = '')
      
      #boxplot(x[, 'nbhsize'] ~ x[, var], outline = FALSE)
      #title(main = 'Neighborhood size',  font.main = 1)
      mtext(lbl, side = 1, outer = TRUE, line = 1, cex = cx)
    }
  })
}

##---- stop ----
p <- bp_base(df = df, Y = c('outdegree', 'size0.015', 'size5'), X = 'risk', lbl = 'Risk level')
p <- bp_base(df = df, Y = c('outdegree', 'size0.015', 'size5'), X =  'age', lbl = 'Age category')
p <- bp_base(df = df, Y = c('outdegree', 'size0.015', 'size5'), X = 'stage', lbl = 'Infection stage')

##---- bp_base ----
p <- bp_base(df = df, Y = c('outdegree', 'size0.05', 'size5'), X = 'risk', lbl = 'Risk level')
p <- bp_base(df = df, Y = c('outdegree', 'size0.05', 'size5'), X =  'age', lbl = 'Age category')
p <- bp_base(df = df, Y = c('outdegree', 'size0.05', 'size5'), X = 'stage', lbl = 'Infection stage')



