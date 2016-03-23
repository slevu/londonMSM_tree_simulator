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
TreeToEdgeList <- function(t, rate = 1, seqlength = NA,
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
    
    ###--- for time-based tree,
    ##- estimate subst along branches
    if (rate != 1){
      print(paste('distances scaled by rate =', rate, 'subst/site/day'))
      t$edge.length <- rpois(length(t$edge.length), 
                             t$edge.length * rate * seqlength ) / seqlength ## (random number of subst)/ number sites
    }
    
    ###--- get distances between tips
    if (file.exists(dist.mat.path)){
      d <- readRDS( dist.mat.path )
    } else {
      d <- as.dist(cophenetic.phylo(t))
      print("Matrix of cophenetic distances was saved")
      saveRDS(d, file = dist.mat.path )
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
### ellipsis for extra arguments, like link function

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

### From phylo-uk
#---- functions ----#

#-------------------#
####---- EdgeList----
#-------------------#
## calculate edge list matrix (2 columns)
## from cluster assignement: vector of cluster number with ID as names

EdgeList <- function(x) {
  # empty edge list
  el <- matrix(nrow=0, ncol = 2)
  # trnasform matrix in df
  x <- as.data.frame(cbind(as.integer(names(x)), 
                           as.integer(x)))
  colnames(x) <- c("id","cluster")
  # for each unique cluster number
  for(clust in unique(x[,"cluster"])){
    # all id belonging to that cluster
    tmp <- x[x[,"cluster"] == clust, "id"]
    # combine elements of tmp by 2, in 2 columns
    # unless cluster size = 1 !
    if(length(tmp) > 1) {
      co <- t(combn(tmp , 2))
      # add rows
      el <-  rbind(el, co)
    }
  }
  colnames(el) <- c("from","to")
  el <- el[order( c(el[, "from"])), ]
  return(el)
}

#-------------------#
####---- AssortMix ----
#-------------------#
## function to 
## match seqindex with patientindex from list of both
## determine from patientID adn to patientID
## add covariate
## if numeric, calculate quantiles
## calculate quantiles-quantiles matrix
## claculate assortativity coefficient

##- eventually program regression analyses

## require:
## edgelist: 3 columns matrix from, to (test or sequence index) and weight W
## listID: dataframe with match between testindex and patientindex
## dfvar: dataframe containing independant variable VAR by patientindex
## var: variable
## option: name of patientindex, name of seqindex or testindex, number of quantiles, 
## png graph to be issued

AssortMix <- function(edgelist = W, listID = df, dfvar = df, var = "dob_y", 
                      pid = "patientindex", sid = "testindex",
                      nqt = 10, graph = 0, trans = 0, method = "") {
  ##- force dataframe
  test <- as.data.frame(edgelist)
  names(test) <- c("from", "to", "W")
  
  ## patientindex at position of listID where edgelist$from = listID$testindex
  test$pid_from <- listID[, pid][match(test[,"from"], listID[, sid])]
  test$pid_to <- listID[, pid][match(test[,"to"], listID[, sid])] 
  
  ##- values of VAR for from_patient and to_patient
  test$var_from <- dfvar[, var][match(test$pid_from, dfvar[, pid])]
  test$var_to <- dfvar[, var][match(test$pid_to, dfvar[, pid])]
  
  ##- if numeric (and if nqt given or quantile option activated ?), 
  ## do quantiles
  if(is.na(nqt)){
    qt_from <- qt_to <- sort(unique(test$var_from))
  } else {
    qt_from <- qt_to <-  quantile(test$var_from, 
                      prob = seq(1/nqt, 1, length.out = nqt),
                      na.rm = TRUE)
  }
  
  #     qt_to <- quantile(test$var_to, 
  #                       prob = seq(1/nqt, 1, length.out = nqt),
  #                       na.rm = TRUE) 
  
  ##- cut in nqt quantiles of continuous VAR
  #   test$var_qt_from <-  cut(test$var_from, 
  #                            qt_from,
  #                            labels=FALSE)
  #   test$var_qt_to <-  cut(test$var_to, 
  #                          qt_to,
  #                          labels=FALSE)
  
  var2q <- function(var) {
    k <- 1
    for (i in qt_from ){
      if (!is.na(var) & var <= i) return (k)
      k <- k + 1
    }
    return (min( k, nqt ))
  } 
  
  test$var_qt_from <- sapply( test$var_from, var2q )
  test$var_qt_to <- sapply( test$var_to, var2q )
  
  ##- empty matrix
  mixing <- matrix(0, nrow = length(qt_from), ncol = length(qt_from))
  
  ##- loop
  #- for each pair
  for(i in 1:nrow(test)){
    #- if both values of VAR exist
    if (!is.na(test$var_qt_from[i]) & !is.na(test$var_qt_to[i])) {
      #- sum weights for each pair of quantiles encountered (row=from, listIDl=to)
      mixing[ test$var_qt_from[i], test$var_qt_to[i]] <- 
        mixing[ test$var_qt_from[i], test$var_qt_to[i]] + test$W[i]
    }
  }
  
  ## optional graph
  if(graph == 1){
    date <- format(Sys.time(),"%Y%m%d_%H%M")
    file <- paste( "figure/", var,"_mixing_", method, date, ".png", sep = "")
    # png(paste( "figure/", var,"_mixing_", date, ".png", sep = ""))
    
    ##- trnasformation sqrt
    if (trans == 0) {
      image(qt_from, qt_to, mixing)
    } else {
      image(sqrt(qt_from), sqrt(qt_to), mixing)
    }
    
    dev.copy(device = png, filename = file, width = 500, height = 500)
    dev.off()
  }
  
  return(mixing)
}

##- verif
# sum(mixing)

#-------------------#
##---- AgeOfInf ----
#-------------------#
## age of infection in years
## loop over iter and patient id to calculate 
## time from intercept to cd4 according to 
## HPA. Longitudinal ... 2011

## sqrt(cd4) = beta0 + beta1_ethn + b2 * time + beta3_ethn * time + beta4 * age * time
## time = [ sqrt(cd4) - beta0 - beta1_ethn] / [beta2 + beta3_ethn + beta4 * age]

AgeOfInf <- function(df, age = "age", ethn = "ethn", 
                     cd4 = "cd4", id = "patientindex",
                     iter = 1){
  
  ##- parameter values from HPA. Longitudinal ... 2011
  
  ethnicity <- c( "black", "other","white")
  
  mu_beta0 <- 24.24
  sd_beta0 <- (24.93 - mu_beta0)/1.96
  mu_beta1 <- c(-1.95, -0.92, 0)
  names(mu_beta1) <- ethnicity # e.g. beta1["black"]
  sd_beta1 <- (c( -1.38, -0.27, 0) - mu_beta1)/1.96 
  names(sd_beta1) <- ethnicity
  
  mu_beta2 <- -0.002
  sd_beta2 <- (-0.0013 - mu_beta2)/1.96
  
  mu_beta3 <- c( 0.0015, 0.0005, 0)
  names(mu_beta3) <- ethnicity
  sd_beta3 <- (c( 0.0020, 0.0012, 0) - mu_beta3)/1.96
  names(sd_beta3) <- ethnicity
  
  mu_beta4 <- -0.00006
  sd_beta4 <- ( -0.00004 - mu_beta4)/1.96
  
  ##- initialize
  # time <- vector()
  # intercept <- vector()
  rr <- matrix(0,  ncol = 4, 
               nrow = length(df[, id]) * iter)
  
  ##- default value for age
  mean_age <- mean(df[ , age], na.rm = TRUE)
  
  ##- loop1
  for (k in 1:iter){
    ##- loop2
    for (i in 1:length(df[, id])){
      
      ##- assign default value for age and ethnicity
      a <- if(is.na(df[ , age][i])) mean_age
      else df[ , age][i]
      e <- if(is.na(df[ , ethn][i])) "white"
      else df[ , ethn][i]
      
      ##- draw parameters
      beta0 <- rnorm(1, mu_beta0, sd_beta0)
      beta1 <- rnorm(1, mu_beta1[e], sd_beta1[e])
      beta2 <- rnorm(1, mu_beta2, sd_beta2)
      beta3 <- rnorm(1, mu_beta3[e], sd_beta3[e])
      beta4 <- rnorm(1, mu_beta4, sd_beta4)
      
      time <-  max( (sqrt(df[ , cd4][i]) - beta0 - beta1) / 
                      (beta2 + beta3 + beta4 * a), 0) # positive values in days
      ##- fill matrix by row number         
      rr[(k - 1)*length(df[, id]) + i, ] <- c(k, df[, id][i], beta0^2, round(time / 365, 2))
    }
  }
  # names
  colnames(rr) <- c("iter", "id", "intercept", "aoi")
  
  ##- statistics by patient
  p <- as.data.frame(
    as.list(
      aggregate(aoi ~ id, data = rr, 
                FUN = function(x) c(m = round(mean(x), 2),
                                    se = round(sd(x)/sqrt(length(x)),2),
                                    lo = round(max(mean(x) - 1.96*sd(x)/sqrt(length(x)), 0), 2),
                                    up = round(mean(x) + 1.96*sd(x)/sqrt(length(x)), 2))) 
    )) # need as.data.frame(as.list()) to get df
  
  ##- return matrix of individual values over iter (ind) and table of summary statistics (pop)
  return(list(ind = rr, pop = p))
}
