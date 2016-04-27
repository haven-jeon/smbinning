#' PSI(Population Stability Index)
#' 
#' @param base_vec base_vec used for reference PSI
#' @param eval_vec eval_vec used for target PSI
#' @param base_ym base_ym reference ymd 
#' @param eval_ym eval_ym target ymd
#' @import dplyr
psi_f <- function(base_vec, eval_vec, base_ym, eval_ym) {
  
  psi_raw <- cbind(
    prop.table(table(base_vec, useNA = 'always')),
    prop.table(table(eval_vec, useNA = 'always'))
  )
  
  psi_raw <- data.frame(psi_raw, row.names = NULL)
  
  f1 <- base_ym
  f2 <- eval_ym
  names(psi_raw)[1:2] <- c(f1, f2)
  #setnames(psi_raw, 2:3, c(f1, f2))
  
  psi_raw <- psi_raw %>% mutate(psi=(get(f2,envir = as.environment(psi_raw))-get(f1, envir = as.environment(psi_raw)))*log(get(f2, envir = as.environment(psi_raw))/get(f1, envir = as.environment(psi_raw)))) 
  
  return(list(psi_detail=data.frame(psi_raw), psi=sum(psi_raw$psi,na.rm = T)))
}

#' kolmogorov-smirnov plot
#' 
#' @param bad_score bad_score scores explain BAD
#' @param n_bad_score n_bad_score scores explain NORMAL
#' @param bad_txt bad_txt text for BAD
#' @param n_bad_txt  n_bad_txt text for NORMAL
#' @return maximum difference in cumulative fraction
plot_ks <- function(bad_score, n_bad_score, bad_txt='BAD', n_bad_txt='NORM' ){
  
  y <- bad_score
  n  <- n_bad_score
  
  f.a <- ecdf(y)
  f.b <- ecdf(n)
  
  x <- seq(min(y, n), max(y, n), length.out=length(n))
  x0 <- x[which( abs(f.a(x) - f.b(x)) == max(abs(f.a(x) - f.b(x))) )]
  y0 <- f.a(x0)
  y1 <- f.b(x0)
  
  plot(f.a, verticals=TRUE, do.points=FALSE, col="blue",main='K-S(Kolmogorov-Smirnov)', xlab='score', ylab='Accumulative Fraction')
  plot(f.b, verticals=TRUE, do.points=FALSE, col="green", add=TRUE)
  
  points(c(x0, x0), c(y0, y1), pch=16, col="red")
  segments(x0, y0, x0, y1, col="red", lty="dotted")
  text(max(x0), min(y0) - 0.1,labels = round(max(y1- y0), 3))
  legend("topleft", legend=c(bad_txt, n_bad_txt), col=c('blue', 'green'), lty = c(1,1))
  return(max(y1- y0))
}

#' Divergence plot 
#' 
#' @param truth factor like TRUE/FALSE 
#' @param pred predicted value 
#' @return Divergence statistics 
#' @import ggplot2 
#' @import dplyr
plot_div <- function(truth, pred){
  restbl <- data.frame(truth=factor(truth), pred=pred)
  restbl_ <- restbl %>% group_by(truth) %>% summarize(sds=sd(pred), avg=mean(pred)) %>% ungroup
  #restbl_ <- restbl[,list(sds=sd(pred), avg=mean(pred)),truth]
  print(ggplot(restbl, aes(pred)) + geom_density(aes(fill=truth), alpha=0.7))
  return((diff(restbl_$avg))^2/(sum(restbl_$sds)/2))
}







