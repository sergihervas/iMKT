## iMK
## compute d (strongly deleterious sites), b (weakly del) and f (neutral)
## includes asymptoticMK output

iMK <- function(x, y, xlow, xhigh) {
  
  m <- as.numeric(y$m0f)
  m0 <- as.numeric(y$m4f)
  d <- as.numeric(y$D0f)
  d0 <- as.numeric(y$D4f)
  
  names(x) <- c("daf","pN","pS")
   
  t1 <- asymptoticMK(x, y, xlow, xhigh)
  a_asym <- as.numeric(t1$alpha_asymptotic)
  a_low <- t1$CI_low ##CIs
  
  x$N <- x$pN/sum(x$pN)           #relative proportion of i
  x$S <- x$pS/sum(x$pS)           #relative prop of 0
  x$pN <- as.numeric(x$pN)
  x$pS <- as.numeric(x$pS)
  x$m <- (x$pN+x$pS)/(sum(x$pN)+sum(x$pS))
  
  ## alpha for each bin
  x$alpha <- 1 - ((d0*x$pN)/(d*x$pS))
  
  ## neutral sites 
  syn <- sum(x$pS)/m0     #syn ratio
  nonsyn <- sum(x$pN)/m   #nonsyn ratio
  f <- nonsyn/syn
  
  ## strongly del (letal) sites
  d <- 1 - f #sites not segregating  
  
  ## weakly deleterious sites
  wd <- 0
  x <- na.omit(x); x <- droplevels(x)
  for (i in 1:nrow(x)) {
    row <- x[i, ]
    
    if (row$alpha < a_low) {
      wd <- wd + ((a_asym - row$alpha) * row$m)
    } else { break }
  }
  
  wd <- wd / (a_asym - min(x$alpha,na.rm=T))
  b <- wd * f
  
  ## real neural sites!
  f <- f - b
  
  ## iMKT information
  out1 <- rbind(cbind("d",d),cbind("b",b),cbind("f",f),cbind("alpha_asymptotic",a_asym))
  out1 <- as.data.frame(out1)
  names(out1) <- c("Parameter","Value")
  out1$Value <- as.numeric(as.character(out1$Value))
  out1$Value <- round(out1$Value, 4)
  rownames(out1) <- NULL
  
  ## alpha asymptotic information
  out2 <- t1
  out2[2:8] <- round(out2[2:8],4)
  
  ## double output: table from iMK and table from asymptoticMK
  out <- list(out1,out2)
  names(out) <- c("iMK","asymptoticMK") 
  return(out) 
}
