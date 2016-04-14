###--- start function
###- summarize regression on bootstrap
compare.reg.bs <- function(ls, reg, model, alpha = 0.05, ...){
  
  ## coef by threshold and by tree
  coef <- lapply(ls, function(x){
    lapply(x , function(x){
      ## choose model
      if ("size" %in% names(x)){
        full.model <- sub("y", "scale(size)", model)
      } else if("outdegree" %in% names(x)){
        full.model <- sub("y", "scale(outdegree)", model)
      } else stop("cannot find y")
      
      coef(summary(reg(formula = full.model, data = x, ...)))
    })
  })
  
  # names(coef[[1]])
  
  ## pvalue by threshold and by tree
  pvalue <- lapply(coef, function(x){
    sapply(x , function(x){
      identity(x[,4])
    })
  })
  
  ##- number of p-value < 0.05
  sum.signif <- sapply(pvalue, function(x){
    apply(x, 1, function(x) sum(x < alpha) / length(x))
  }
  )
  
  ## parameter by threshold
  param <-  lapply(coef, function(x){
    sapply(x , function(x){
      identity(x[,1])
    })
  })
  
  ## mean of parameter
  mean.parms <- signif(sapply(param, function(x){
    apply(x, 1, mean)
  }), 2)
  
  ## R square, only for lm()
  if(identical(reg, lm)){
    r2 <- lapply(ls, function(x){
      sapply(x , function(x){
        ## choose model
        if ("size" %in% names(x)){
          full.model <- sub("y", "scale(size)", model)
        } else if("outdegree" %in% names(x)){
          full.model <- sub("y", "scale(outdegree)", model)
        } else stop("cannot find y")
        
        summary(reg(full.model, data = x))$r.squared
      })
    })
    ## mean R2
    mean.r2 <- signif(sapply(r2, function(x){
      mean(x)
    }), 3)
    
    return(list("model" = model, "mean parameter" = mean.parms, "signif pvalue" = sum.signif, "mean r.squared" = mean.r2)) 
  } else {
    
    return(list("model" = model, "mean parameter" = mean.parms, "signif pvalue" = sum.signif))
  }
}
###--- end function 