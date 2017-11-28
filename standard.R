## Standard DGRP

standard <- function(x, y) {
  
  out <- NULL
  ##### Not defined####
  #names(x) <- c("daf","pN","pS")
  mi <- as.numeric(y$m0f)
  m0 <- as.numeric(y$m4f)
  di <- as.numeric(y$D0f)
  d0 <- as.numeric(y$D4f)
  

    #create table with raw values and < and > remove% values
    t <- rbind(cbind(sum(x$pS),d0),cbind(sum(x$pN),di))
  
    P0 <- t[1,1]; Pi <- t[2,1]
    D0 <- t[1,2]; Di <- t[2,2]
    
    alpha <- 1 - ((Pi/P0) * (D0/Di))
    
    m <- matrix(c(P0, Pi,D0,Di), ncol=2)
    m1 <- fisher.test(m); m1 <- m1$p.value
    
    g1 <- cbind(alpha,m1)
    names(g1) <- c("alpha","pval")
    out <- rbind(out,g1)
    names(out) <- c("alpha","pval")
  
  
  out <- as.data.frame(out)
  names(out) <- c("alpha","pval")
  return(out)
}
