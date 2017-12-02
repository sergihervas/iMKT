#' iMK compute all the MKT extensions 
#' Date = 30/11/2016
#' Author = Sergi Hervás, Marta Coronado
#'
#' EXPLANATION
#'
#'
#' @param daf dad file
#' @param divergence div file
#' @param xlow fit curve
#' @param xhigh curv
#'
#' @return None
#'
#' @examples
#' @import knitr 
#' @import utils
#' @import stats
#' @import ggplot2 
#' @importFrom ggplot2 ggsave
#' @import ggthemes
#' @importFrom ggthemes theme_map
#' @import cowplot
#' @import grid 
#' @import gridExtra
#' @import scales
#' @import reshape2 
#' @export
#' 

iMK <- function(daf, divergence, xlow, xhigh) {
  
  check<-check_input(daf, divergence, xlow, xhigh)
  if(check$data==FALSE)
    stop(check$print_errors)
  
  
  m <- as.numeric(divergence$m0f)
  m0 <- as.numeric(divergence$m4f)
  d <- as.numeric(divergence$D0f)
  d0 <- as.numeric(divergence$D4f)
  
  ####Not defined###
  names(daf) <- c("daf","pN","pS")
   
  t1 <- asymptoticMK(daf, divergence, xlow, xhigh)
  a_asym <- as.numeric(t1$alpha_asymptotic)
  a_original <- as.numeric(t1$alpha_original)
  a_low <- t1$CI_low ##CIs
  
  daf$N <- daf$pN/sum(daf$pN)           #relative proportion of i
  daf$S <- daf$pS/sum(daf$pS)           #relative prop of 0
  daf$pN <- as.numeric(daf$pN)
  daf$pS <- as.numeric(daf$pS)
  
  ## alpha for each bin
  daf$alpha <- 1 - ((d0*daf$pN)/(d*daf$pS))
  
  ## neutral sites 
  syn <- sum(daf$pS)/m0     #syn ratio
  nonsyn <- sum(daf$pN)/m   #nonsyn ratio
  f <- nonsyn/syn
  
  ## strongly del (letal) sites
  d <- 1 - f #sites not segregating  
  
  daf <- na.omit(daf); daf <- droplevels(daf)
  
  ## weakly deleterious sites
  wd <- 0
  for (i in 1:nrow(daf)) {
    row <- daf[i, ]
    
    if (row$alpha < a_low) {
      wd <- wd + ((a_asym - row$alpha) * row$N)
                   # row$N = row$pN / sum(daf$pN)
    } else { break }
  }
  
  wd <- wd / (a_asym - min(daf$alpha,na.rm=T))
  b <- wd * f
  
  ## real neural sites!
  f <- f - b
  
  ## include DGRP alphas
  dgrp <- mkt_DGRP(daf, divergence)
  a_05 <- dgrp$Results[2,2]
  a_10 <- dgrp$Results[3,2]
  
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
  #out3[2:3] <- round(out3[2:3],4)
  
  out <- list(out1,out2,out3)
  names(out) <- c("iMK","asymptoticMK","DGRP")
  return(out) 
}
