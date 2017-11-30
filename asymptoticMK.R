## only fits exponential model, depends on fitMKmodel and predictNLS
## includes watchdog
## returns only a table with results, no plot

asymptoticMK <- function(x, y, xlow=0, xhigh=1) {
  
  ## loading required packages
  require(MASS, quietly=TRUE)
  require(nls2, quietly=TRUE)
  
  ## check data: if there is an error, watchdog stops computation
  check<-check_input(x, y, xlow, xhigh)
  if(check$data==FALSE)
    stop(check$print_errors)
  
  tryCatch({
    if(any(x[3]<=0)){
      result_df <- data.frame(model="exponential", a=NaN, b=NaN, c=NaN, alpha_asymptotic=NaN, CI_low=NaN, CI_high=NaN, alpha_original=NaN, row.names=NULL)
      return(result_df)
      }}, error=cat("Your input files x has the following(s) errors: \n x$P0 has one o more 0 values, cannot compute asymptotic. Returning results with NaN \n \n")
    )
  
  ## assign proper names to the columns of x and y 
  #####Not defined####
  #names(x) <- c("daf", "pN", "pS")
  #names(y) <- c("Chr",  "Pop",  "m0f",  "D0f",  "m4f",  "D4f")
  
  ## parse the data from argument x
  f <- x$daf #derived alelle frequencies
  p <- x$pN #non-synonymous polymorphism 
  p0 <- x$pS #synonymous polymorphism
  
  ## parse the data from argument y
  m <- as.numeric(y$m0f) #number of non-synonymous analyzed positions   
  m0 <- as.numeric(y$m4f) ##number of synonymous analyzed positions
  d <- as.numeric(y$D0f) #non-synonymous divergence
  d0 <- as.numeric(y$D4f) #synonymous divergence
  # pop <- y$Pop #population name, used for plotting...
  
  ## compute alpha values and trim
  alpha <- 1 - (d0/d) * (p/p0)
  cutoff_f1 <- xlow
  cutoff_f2 <- xhigh
  trim <- ((f >= cutoff_f1) & (f <= cutoff_f2))
  f_trimmed <- f[trim]
  alpha_trimmed <- alpha[trim]
  
  ## compute the original MK alpha
  alpha_nonasymp <- 1 - (d0/d) * (sum(p[trim])/sum(p0[trim])) #using trimmed values
      
  ## two-step nls2() model fit at a given level of precision (res)
  mod1 <- fitMKmodel(alpha_trimmed, f_trimmed, 10)
  
  ## if mod1 did not work, try a deeper scan for a decent fit (res=20)
  if (length(mod1) == 0) {
    mod1 <- fitMKmodel(alpha_trimmed, f_trimmed, 20)
  } 
  
  tryCatch({
    mod2 <- lm(alpha_trimmed ~ f_trimmed)
  }, error=function(cond) {})
  
  ## compute confidence intervals of alpha using predictNLS 
  ci_pred <- predictNLS(mod1, newdata=data.frame(f_trimmed=1.0))
  alpha_1_low <- ci_pred[6]
  alpha_1_high <- ci_pred[7]
  
  ## preparation of ouput (alpha asym, a, b, c)
  alpha_1_est <- predict(mod1, newdata=data.frame(f_trimmed=1.0))
  const_a <- coef(mod1)["const_a"]
  const_b <- coef(mod1)["const_b"]
  const_c <- coef(mod1)["const_c"]
  
  ## output table
  result_df <- data.frame(model="exponential", a=const_a, b=const_b, c=const_c, alpha_asymptotic=alpha_1_est, CI_low=ci_pred[6], CI_high=ci_pred[7], alpha_original=alpha_nonasymp, row.names=NULL)
  return(result_df)
}

