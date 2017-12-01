#' iMK compute all the MKT extensions 
#' Date = 30/11/2016
#' Author = Sergi Hervás, Marta Coronado
#'
#' EXPLANATION
#'
#'
#' @param x dad file
#' @param y divergence file
#' @param xlow fit curve
#' @param ylow fit curv
#'
#' @return None
#'
#' @examples
#' asymptoticMK
#' mkt_DGPR
#' mkt_fww
#' mkt_standard
#' @export
#' 


## iMK
## compute d (strongly deleterious sites), b (weakly del) and f (neutral)
## includes asymptoticMK output

### IMPORTANT MESSAGE: ###
## If we want to compare datasets, they must contain a very similar number of unrelated samples, as the
## whole methodology depends on the ratio: segregating sites / analyzed sites.
## If one data file contains many more samples and therefore, more segregating sites (assuming that bigger
## the sample size, higher the number of rare SNPs detected), results may be strongly biased.

iMK <- function(x, y, xlow, xhigh) {
  
  check<-check_input(x, y, xlow, xhigh)
  if(check$data==FALSE)
    stop(check$print_errors)
  
  
  m <- as.numeric(y$m0f)
  m0 <- as.numeric(y$m4f)
  d <- as.numeric(y$D0f)
  d0 <- as.numeric(y$D4f)
  
  ####Not defined###
  #names(x) <- c("daf","pN","pS")
   
  # t1 <- asymptoticMK(x, y, xlow, xhigh)
  # a_asym <- as.numeric(t1$alpha_asymptotic)
  # a_original <- as.numeric(t1$alpha_original)
  # a_low <- t1$CI_low ##CIs
  
  x$N <- x$pN/sum(x$pN)           #relative proportion of i
  x$S <- x$pS/sum(x$pS)           #relative prop of 0
  x$pN <- as.numeric(x$pN)
  x$pS <- as.numeric(x$pS)
  
  ## alpha for each bin
  x$alpha <- 1 - ((d0*x$pN)/(d*x$pS))
  
  ## neutral sites 
  syn <- sum(x$pS)/m0     #syn ratio
  nonsyn <- sum(x$pN)/m   #nonsyn ratio
  f <- nonsyn/syn
  
  ## strongly del (letal) sites
  d <- 1 - f #sites not segregating  
  
  x <- na.omit(x); x <- droplevels(x)
  
  ## weakly deleterious sites
  wd <- 0
  for (i in 1:nrow(x)) {
    row <- x[i, ]
    
    if (row$alpha < a_low) {
      wd <- wd + ((a_asym - row$alpha) * row$N)
                   # row$N = row$pN / sum(x$pN)
    } else { break }
  }
  
  wd <- wd / (a_asym - min(x$alpha,na.rm=T))
  b <- wd * f
  
  ## real neural sites!
  f <- f - b
  
  ## include DGRP alphas
  dgrp <- DGRP(x, y)
  a_05 <- dgrp[1,2]
  a_10 <- dgrp[2,2] 
  
  ## iMKT information
  out1 <- rbind(cbind("d",d),cbind("b",b),cbind("f",f),cbind("alpha_original",a_original),cbind("alpha_asymptotic",a_asym),cbind("alpha_DGRP05",a_05),cbind("alpha_DGRP10",a_10))
  out1 <- as.data.frame(out1)
  names(out1) <- c("Parameter","Value")
  out1$Value <- as.numeric(as.character(out1$Value))
  out1$Value <- round(out1$Value, 4)
  rownames(out1) <- NULL
  
  ## alpha asymptotic information
  out2 <- t1
  out2[2:8] <- round(out2[2:8],4)
  
  ## alpha DGRP information
  out3 <- dgrp
  out3[2:3] <- round(out3[2:3],4)
  
  out <- list(out1,out2,out3)
  names(out) <- c("iMK","asymptoticMK","DGRP")
  return(out) 
}
