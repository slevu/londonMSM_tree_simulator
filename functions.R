###--- Functions ---###
### used for LondonMSM_tree_simulator

######################################
####---- Start TreeToEdgeList ----####
######################################
##-----------------------------------------##
####---- function converting timed tree
####---- in edge list of pairwise distances 
##-----------------------------------------##

###--- input: timed tree, substitution rate
###--- output: returns path for dataframe (ID1, ID2, distance)

###--- Scales the time based distance by consensus mutation rate to obtain 'standard' substitution per site distances
# tree.name <- "simtree" tree.name <- "uktree"
TreeToEdgeList <- function(t, rate = 1, 
                           output = "data/", fig = "figure/" , 
                           stats = TRUE, plot = TRUE ){
  
  ###--- object names for output path
  tree.name <- substitute(t)
  dist.mat.name <- paste(tree.name, "_dist", sep = '')
  dist.mat.path <- paste(output, dist.mat.name,".rds", sep = '')
  plot.name <- paste(fig, dist.mat.name, ".png", sep ='')
  out_edge_list <- paste(output, tree.name, "_el.rds", sep = '')
  
  ## if file exists, don't run and return the path 
  if (file.exists(out_edge_list) ){
    print(paste('Edge list already exists here:', out_edge_list))
    return(out_edge_list)
  } else {
    
    ###--- get distances
    if (file.exists(dist.mat.path)){
      d <- readRDS( dist.mat.path )
    } else {
      d <- as.dist(cophenetic.phylo(t))
      print("Matrix of cophenetic distances was saved")
      saveRDS(d, file = dist.mat.path )
    }
    
    ###--- time to rate
    if (rate != 1){
      print(paste('distances scaled by rate =', rate, 'subst/site/day'))
      d <-  d * rate
    }
    
    ###--- stats
    if(stats == TRUE){ print(summary(d)) }
    if(plot == TRUE){ 
      print(paste("png plot saved:", plot.name))
      png(filename = plot.name, type="cairo",
          units="in", width=5, height=4, 
          pointsize=12, res=96)
      hist(d, breaks = 50, xlim = c(0, 1),
           xlab = "distance", ylab = "frequency",
           main = '', #tree.name # Tree's distances in subst/site
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
    if (file.exists(out_edge_list)){
      el <- readRDS(out_edge_list)
    } else {
      # system.time(
      el <- melt(m, na.rm = TRUE)
      # ) # the na.rm removal takes most of time
      colnames(el) <- c('ID1', 'ID2', 'distance')
      saveRDS(el, file = out_edge_list)
      print(paste("Edge list was saved:", out_edge_list))
    }
    return(out_edge_list)
  }
} 
####################################
####---- End TreeToEdgeList ----####
####################################


#####################################
####---- Start ucsd_hivclust ----####
#####################################
##----------------------------##
####--- Apply UCSD
####--- HIV clustering
####--- to a dataframe of 
####--- pairwise distances 
##----------------------------##

###- input: path to RDS file of edge list
###- output: cluster assignements
###--- loops through distance thresholds ( = vector of threshold based on quantiles) 

ucsd_hivclust <- function(path.el, thr = NA, quant = c(5e-5, 1e-4, 5e-4, 1e-3, 1e-2, 1e-1) ){
  
  ## read RDS file
  el <- readRDS(file = path.el)
  
  ## get var.name for output path (Regex !)
  var.name <- paste("d", substr(path.el,
                                regexpr("\\/[^\\/]*$", path.el)[[1]][1] +1,
                                regexpr("\\_el.rds", path.el)[[1]][1] - 1 ),
                    sep = '')
  
  ## Either thr given or based on quantiles
  if(is.na(thr[1])){
    
    # choose threshold based on quantiles (first 4 for now)
    qt <- quantile(el$distance, 
                   probs = quant )
    thr <- round(qt[1:4], 2)
    
    ## get rid of distance > larger quantile
    #      k <- round(qt[length(qt)], 2)
    #      subel <- el[el$distance < k,]
    
  } else { 
    thr <-  thr
    qt <- "Not used"
  }
  
  ## write csv without rownames and get input path
  inputCSV <- paste(tempdir(), "/input.csv", sep = "")
  write.csv(el, file = inputCSV, row.names = FALSE )
  
  ## full path needed 
  exec <- '~/Documents/softwares/hivclustering/scripts/hivnetworkcsv'
  
  ####---- loop threshold (first 4 qt = 0.05, 0.1, 1, 10) 
  
  ## empty results
  cmd <- vector( mode= "character" )
  warn <- list()
  
  for ( t in thr ){
    print(paste("threshold =", t))
    ## output
    outputCSV <- paste(var.name, "_ucsd_hivclust_output_", 
                       t, ".csv", sep = '')
    
    ## parms
    parms <- paste("-t", t, "-f plain")
    
    ## command
    cmd_hivclustering <- paste(exec, "-i", inputCSV, "-c", outputCSV, parms )
    
    print(cmd_hivclustering)
    
    ## too long, the fitting of degree
    ## distribution takes most time
    ## => just issue command
    ## to run on terminal
    if ( !file.exists(outputCSV) ){
      stderr <-  system(
        paste(cmd_hivclustering, "2>&1"),
        intern = TRUE)
    }
    # save commands and 'stderr' warnings
    cmd <- c(cmd, cmd_hivclustering)
    warn <- c(warn, c(t, stderr))
  }
  
  ###--- bin table in one list
  cl <- list()
  for(i in 1:length(thr)){
    ## add table i
    CSV <- paste(var.name, 
                 "_ucsd_hivclust_output_", thr[i],
                 ".csv", sep = '')
    if(file.exists(CSV)){
      cl[[i]] <- read.csv(CSV)
      # name threshold
      names(cl)[i] <- thr[i]
    }
  }
  
  return(list(qt = qt, cmd = cmd, warn = warn, cl = cl))
}
###################################
####---- End ucsd_hivclust ----####
###################################


###############################
####---- Start reg.sum ----####
###############################
###-- Summary results for both
### linear and logistic regression
### return only parameters, p-values and R squared
### input: list of df, regression function, model formula,
### ellipsis for extre arguments, like link function

reg.sum <- function(ls, reg, model, ...){
  ## pvalue by threshold
  pvalue <- sapply(ls , function(x){
    coef(summary(reg(formula = model, data = x, ...)))[,4]
  })
  ## coded by significance
  pvalue.code <- as.data.frame(
    apply(pvalue,2, function(x){
      cut(x,
          breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
          labels=c('***','**', '*', '.', ''))
    }))
  row.names(pvalue.code) <- row.names(pvalue)
  ## parameter by threshold
  param <- signif(sapply(ls , function(x){
    coef(summary(reg(formula = model, data = x, ...)))[,1]
  }), 2)
  ## R square, only for lm()
  if(identical(reg, lm)){
    r2 <- signif(sapply(ls , function(x){
      summary(reg(model, data = x))$r.squared
    }), 3)
    
    return(list("model" = model, "parameter" = param, "pvalue" = pvalue.code, "r.squared" = r2)) 
  } else {
    
    return(list("model" = model, "parameter" = param, "pvalue" = pvalue.code))
  }
}
#############################
####---- End reg.sum ----####
#############################


#####################################
####---- Start size.vs.covar ----####
#####################################

####---- plot correlation ----
#### input: list of
#### df containing indepvar and depvar
size.vs.covar <- function(l, depvar = "size",
                          indepvar = c("stage", "risk", "age", "time")){
  
  par(mfcol=c(length(indepvar), length(l) ), 
      mar = c(4,3,3,2)+0.1, oma = c(0, 0, 2, 0), bty = 'n') # b,l,t,r
  for(c in 1:length(l)){
    for(r in 1:length(indepvar)){
      plot( x = l[[c]][, indepvar[r]], y = l[[c]][, depvar], 
            ylab = '', xlab = '', font.main = 1,
            col="#00000050",
            main = paste(indepvar[r], names(l)[c]))
    }
  }
  mtext(paste(depvar, "(y) vs co-variates (x) "), outer = TRUE, cex = 1)
}

# size.vs.covar(down)
# dev.off()
###################################
####---- End size.vs.covar ----####
###################################