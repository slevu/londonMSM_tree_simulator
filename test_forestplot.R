#### TODO
#### Multiple confidence bands
#### for naive vs down-sampling
#### at the same threshold

### 1. same model at diffreent threshold
library(forestplot)

  a <- synth_hsh[synth_hsh$cage2 == 'Total'
                 & synth_hsh$reg == 'IDF',]
head(a)

  forestplot(paste(a$study,a$year), a$p_hsh12, a$lo_p_hsh12, a$up_p_hsh12,
             cex  = .5, lineheight = "auto",
             xlab = "",
             col=fpColors(box = gray(0:1/2)), 
             #                         line="darkblue", summary="royalblue", 
             #                         hrz_lines = "#444444"),
             vertices = TRUE,
             new_page = TRUE)
  
  dev.off()
  
fit <- lapply(listclus, function(x) lm(model, data = x))


co <- lapply(fit, coef)
ic <- lapply(fit, confint)


mean <- t(matrix(unlist(co), ncol = length(co[[1]]), byrow = TRUE))

low_up <- t(matrix(unlist(ic), ncol = length(ic[[1]]), byrow = TRUE))
low <- low_up[1:(length(ic[[1]])/2) , ]
up <- low_up[((length(ic[[1]])/2)+1):length(ic[[1]]) , ]


# summary(fit)
# str(confint(fit))
# cbind( coef(fit), confint(fit))

tabletext <- rownames(ic[[1]])
#num <- cbind( coef(fit), confint(fit))

forestplot(tabletext, mean, low, up, 
           legend = names(co),
           xticks = seq(-0.5, 0.5, by = 0.1),
           cex  = 1, lineheight = "auto",
           xlab = "",
           col=fpColors(box = rainbow(4)), 
#                         line="darkblue", summary="royalblue", 
#                         hrz_lines = "#444444"),
           vertices = TRUE,
           new_page = TRUE)

###---
### naive and downsampling
###---
tail(mean.down_uk[[3]])
tail(listUKclus[[3]])
lm_model_uk  <-  "scale(size) ~ scale(agediag) + scale(sqrt(cd4)) +  scale(ydiag)"
## anonymous function
fit1 <- (function(x) lm(lm_model_uk, data = x))(listUKclus[[2]])
fit2 <- (function(x) lm(lm_model_uk, data = x))(mean.down_uk[[2]])
summary(fit1)
summary(fit2)

mean <- cbind(coef(fit1), coef(fit2))
low <- cbind(confint(fit1)[,1] , confint(fit2)[,1])
up <- cbind(confint(fit1)[,2] , confint(fit2)[,2] )

tabletext <- names(coef(fit1))
#num <- cbind( coef(fit), confint(fit))

## at threshold level 2 (0.02)
forestplot(tabletext, mean, low, up, 
           legend = c("naive", "down-sampled"),
           #xticks = seq(-0.2, 0.4, by = 0.1),
           #boxsize = .05,
           cex  = .5, lineheight = "auto",
           xlab = "",
           col=fpColors(box = gray(0:1/2)), 
           #                         line="darkblue", summary="royalblue", 
           #                         hrz_lines = "#444444"),
           vertices = TRUE,
           new_page = TRUE)
dev.off()
?"forestplot"
