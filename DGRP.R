## DGRP alpha removing 5 and 10%
##Remove % DGRP
### MUST OPTIMIZE! ###

DGRP <- function(x, y) {
  
  out <- NULL
  names(x) <- c("daf","pN","pS")
  mi <- as.numeric(y$m0f)
  m0 <- as.numeric(y$m4f)
  di <- as.numeric(y$D0f)
  d0 <- as.numeric(y$D4f)
  
  ns <- c(0.05, 0.1)
  for (n in ns) {
    
    x1 <- x[x$daf < n, ] #below DAF
    x2 <- x[x$daf > n, ] #over DAF
    
    #create table with raw values and < and > remove% values
    t <- rbind(cbind(sum(x$pS),d0),cbind(sum(x$pN),di))
    tab1 <- rbind(cbind(sum(x1$pS),sum(x2$pS)),cbind(sum(x1$pN),(sum(x2$pN))))
    t <- cbind(t,tab1); t <- as.data.frame(t)
    
    P0 <- t[1,1]; Pi <- t[2,1]
    D0 <- t[1,2]; Di <- t[2,2]
    P0b <- t[1,3]; Pib <- t[2,3]
    P0o <- t[1,4]; Pio <- t[2,4]
    
    fneutralb <- P0b / P0
    Pineutralb <- Pi * fneutralb
    Pi_wd <- Pib - Pineutralb
    Pineutral <- round(Pineutralb + Pio)
    if(P0 == "NaN") {P0 <- 0}
    if(Pineutral == "NaN") {Pineutral <- 0}
    alpha <- 1 - ((Pineutral/P0) * (D0/Di))
    
    m <- matrix(c(P0,Pineutral,D0,Di), ncol=2)
    m1 <- fisher.test(m); m1 <- m1$p.value
    
    g1 <- cbind(n,alpha,m1)
    names(g1) <- c("cutoff","alpha","pval")
    out <- rbind(out,g1)
    names(out) <- c("cutoff","alpha","pval")
  } 
  
  out <- as.data.frame(out)
  names(out) <- c("cutoff","alpha","pval")
  return(out)
}
