#' @title Asymptotic MKT method
#' 
#' @description MKT calculation using asymptoticMK method (Messer and Petrov 2012 PNAS; Haller and Messer 2017 G3)
#'
#' @details In the standard McDonald and Kreitman test, the estimate of adaptive evolution (alpha) can be easily biased by the segregation of slightly deleterious non-synonymous substitutions. Specifically, slightly deleterious mutations contribute more to polymorphism than they do to divergence, and thus, lead to an underestimation of alpha. Messer and Petrov proposed a simple asymptotic extension of the MK test that yields accurate estimates of alpha. Briefly, this method first estimates alpha for each DAF category using its specific Pi and P0 values and then fits an exponential function to this values, of the form: alpha Fit(x) = a + b exp(-cx). Although the exponential function is generally expected to provide the best fit, a linear function is also fit to the data, of the form: alpha Fit(x) = a + bx. Finally, the asymptotic alpha estimate is obtained by extrapolating the value of this function to x = 1: alpha Asymptotic = alpha Fit (x=1). The exponential fit is always reported, except if the exponential fit fails to converge or if the linear fit is superior according to AIC. The code of this function is adapted from Haller and Messer 2017 G3 (http://github.com/MesserLab/asymptoticMK).
#'
#' @param daf data frame containing DAF, Pi and P0 values
#' @param divergence data frame containing divergent and analyzed sites for selected (i) and neutral (0) classes
#' @param xlow lower limit for asymptotic alpha fit
#' @param xhigh higher limit for asymptotic alpha fit
#' @param seed seed value (optional). No seed by default
#'
#' @return Estimation of asymptotic alpha and details about the model fit (function parameters, confidence intervals, etc.)
#'
#' @examples
#' asymptoticMK(myDafData, myDivergenceData, xlow=0, xhigh=0.9)
#'
#' @import utils
#' @import stats
#' @importFrom MASS mvrnorm
#' @importFrom nls2 nls2
#'
#' @keywords MKT
#' @export

asymptoticMK <- function(daf, divergence, xlow, xhigh, seed) {
  
  ## Check data
  check <- checkInput(daf, divergence, xlow, xhigh)
  if(check$data == FALSE) {
    stop(check$print_errors) }

  if (any(daf$P0 == 0)){ ## Warning P0
    warning("Input daf file contains P0 values = 0.\nThis can bias the function fitting and the estimation of alpha.")}

  ## Check seed
  if(missing(seed)) {
    seed <- NULL
  } else {
   set.seed(seed)
  }
  
  ## Parse the data from argument x
  f <- daf$daf #derived alelle frequencies
  p <- daf$Pi #non-synonymous polymorphism 
  p0 <- daf$P0 #synonymous polymorphism
  
  ## Parse the data from argument y
  m <- divergence$mi #number of non-synonymous analyzed positions   
  m0 <- divergence$m0 ##number of synonymous analyzed positions
  d <- divergence$Di #non-synonymous divergence
  d0 <- divergence$D0 #synonymous divergence
  
  ## Compute alpha values and trim
  alpha <- 1 - (d0/d) * (p/p0)
  cutoff_f1 <- xlow
  cutoff_f2 <- xhigh
  trim <- ((f >= cutoff_f1) & (f <= cutoff_f2))
  f_trimmed <- f[trim]
  alpha_trimmed <- alpha[trim]
  
  ## Compute the original MK alpha
  alpha_nonasymp <- 1 - (d0/d) * (sum(p[trim])/sum(p0[trim])) #using trimmed values
      
  ## Two-step nls2() model fit at a given level of precision (res)
  fitMKmodel <- function(alpha_trimmed, f_trimmed, res) {
    
    ## First fitting using starting values (st)
    mod <- tryCatch({
      
      ## Starting values to fit the model  
      st <- expand.grid(const_a=seq(-1,1,length.out=res + 1), const_b=seq(-1,1,length.out=res), const_c=seq(1,10,length.out=res + 1))
      
      ## Fitting
      nls2(alpha_trimmed ~ const_a + const_b * exp(-const_c* f_trimmed), start=st, algorithm="brute-force", control=nls.control(maxiter=NROW(st)))
      
    }, error=function(cond) {}) ## Return condition of error when unable to fit
    
    ## If mod fails...
    if (length(mod) == 0) { return(NULL) }
    
    ## Second fitting, starting from previous fit (mod)
    mod2 <- tryCatch({
      nls2(alpha_trimmed ~ const_a + const_b * exp(-const_c* f_trimmed), start = mod, control=nls.control(maxiter=200))
      
    }, error=function(cond) {}) ## Same error handling than the previous step
    
    ## If mod2 fails...
    if (length(mod2) == 0) { return(NULL) }
    
    ## Return mod2 if fitted
    return(mod2)
  }
  
  mod1 <- fitMKmodel(alpha_trimmed, f_trimmed, 10)
  
  ## If mod1 did not work, try a deeper scan for a decent fit (res=20)
  if (length(mod1) == 0) {
    mod1 <- fitMKmodel(alpha_trimmed, f_trimmed, 20)
  } 

  tryCatch({
    mod2 <- lm(alpha_trimmed ~ f_trimmed)
  }, error=function(cond) {})
  
  ## Compute confidence intervals of alpha using predictNLS 
  ## Get a CI using Monte Carlo simulation based upon a fitted model.  
  ## Thanks to Andrej-Nikolai Spiess (http://www.dr-spiess.de) for this code.
  predictNLS <- function(object, newdata, level = 0.95, nsim = 10000) {
    
    ## get right-hand side of formula
    RHS <- as.list(object$call$formula)[[3]]
    EXPR <- as.expression(RHS)
    
    ## all variables in model
    VARS <- all.vars(EXPR)
    
    ## coefficients
    COEF <- coef(object)
    
    ## extract predictor variable    
    predNAME <- setdiff(VARS, names(COEF))  
    
    ## take fitted values, if 'newdata' is missing
    if (missing(newdata)) {
      newdata <- eval(object$data)[predNAME]
      colnames(newdata) <- predNAME
    }
    
    ## check that 'newdata' has same name as predVAR
    if (names(newdata)[1] != predNAME) stop("newdata should have name '", predNAME, "'!")
    
    ## get parameter coefficients
    COEF <- coef(object)
    
    ## get variance-covariance matrix
    VCOV <- vcov(object)
    
    ## augment variance-covariance matrix for 'mvrnorm' 
    ## by adding a column/row for 'error in x'
    NCOL <- ncol(VCOV)
    ADD1 <- c(rep(0, NCOL))
    ADD1 <- matrix(ADD1, ncol = 1)
    colnames(ADD1) <- predNAME
    VCOV <- cbind(VCOV, ADD1)
    ADD2 <- c(rep(0, NCOL + 1))
    ADD2 <- matrix(ADD2, nrow = 1)
    rownames(ADD2) <- predNAME
    VCOV <- rbind(VCOV, ADD2) 
    
    ## iterate over all entries in 'newdata' as in usual 'predict.' functions
    NR <- nrow(newdata)
    respVEC <- numeric(NR)
    seVEC <- numeric(NR)
    varPLACE <- ncol(VCOV)   
    
    ## define counter function
    counter <- function(i) {
      if (i%%10 == 0) { cat(i) 
      } else { cat(".") }
      if (i%%50 == 0) { cat("\n") }
      flush.console()
    }
    
    ## create output matrix (df)
    outMAT <- NULL 
    
    for (i in 1:NR) {
      
      ## get predictor values and optional errors
      predVAL <- newdata[i, 1]
      if (ncol(newdata) == 2) predERROR <- newdata[i, 2] else predERROR <- 0
      names(predVAL) <- predNAME  
      names(predERROR) <- predNAME  
      
      ## create mean vector for 'mvrnorm'
      MU <- c(COEF, predVAL)
      
      ## create variance-covariance matrix for 'mvrnorm'
      ## by putting error^2 in lower-right position of VCOV
      newVCOV <- VCOV
      newVCOV[varPLACE, varPLACE] <- predERROR^2
      
      ## create MC simulation matrix
      simMAT <- mvrnorm(n = nsim, mu = MU, Sigma = newVCOV, empirical = TRUE)
      
      ## evaluate expression on rows of simMAT
      EVAL <- try(eval(EXPR, envir = as.data.frame(simMAT)), silent = TRUE)
      if (inherits(EVAL, "try-error")) stop("There was an error evaluating the simulations!")
      
      ## collect statistics
      PRED <- data.frame(predVAL)
      colnames(PRED) <- predNAME   
      FITTED <- predict(object, newdata = data.frame(PRED))
      MEAN.sim <- mean(EVAL, na.rm = TRUE)
      SD.sim <- sd(EVAL, na.rm = TRUE)
      MEDIAN.sim <- median(EVAL, na.rm = TRUE)
      MAD.sim <- mad(EVAL, na.rm = TRUE)
      QUANT <- quantile(EVAL, c((1 - level)/2, level + (1 - level)/2))
      RES <- c(FITTED, MEAN.sim, SD.sim, MEDIAN.sim, MAD.sim, QUANT[1], QUANT[2])
      outMAT <- rbind(outMAT, RES)
    }
    
    colnames(outMAT) <- c("fit", "mean", "sd", "median", "mad", names(QUANT[1]), names(QUANT[2]))
    rownames(outMAT) <- NULL   
    return(outMAT)      
  }
  
  ## Compare linear and exponential fit
  linear_better <- FALSE

  if ((length(mod1) == 0) || (AIC(mod2) < AIC(mod1))) {
    linear_better <- TRUE }

  ## If linear is not better, check wide of confidence intervals of exp fit
  if (!linear_better) {
    
    tryCatch({
      ci_pred <- predictNLS(mod1, newdata=data.frame(f_trimmed=1.0))
      alpha_1_low <- ci_pred[6]
      alpha_1_high <- ci_pred[7]
  
      if ((alpha_1_low < -100) || (alpha_1_high > 100)) {
        linear_better <- TRUE }
  
    }, error=function(cond) {cat("Could not compute CI for the exponential alpha fit.\n")})
  }

  ## If linear fit better than exp fit
  if (linear_better) {

    ## Predict linear model confidence, not prediction
    ci_pred <- predict.lm(mod2, newdata=data.frame(f_trimmed=1.0), interval="confidence")
    alpha_1_low <- ci_pred[2]
    alpha_1_high <- ci_pred[3]
    
    alpha_1_est <- predict(mod2, newdata=data.frame(f_trimmed=1.0))
    const_a <- coef(mod2)["(Intercept)"]
    const_b <- coef(mod2)["f_trimmed"]
    const_c <- NA
  
  ## If exp is the best fit
  } else {
    ## Preparation of ouput (alpha asym, a, b, c)
    alpha_1_est <- predict(mod1, newdata=data.frame(f_trimmed=1.0))
    const_a <- coef(mod1)["const_a"]
    const_b <- coef(mod1)["const_b"]
    const_c <- coef(mod1)["const_c"]
    alpha_1_low <- ci_pred[6]
    alpha_1_high <- ci_pred[7]
  }

  ## Output table
  result_df <- data.frame(model=(if ((length(mod1) == 0) || linear_better) "linear" else "exponential"), a=const_a, b=const_b, c=const_c, alpha_asymptotic=alpha_1_est, CI_low=alpha_1_low, CI_high=alpha_1_high, alpha_original=alpha_nonasymp, row.names=NULL)
  return(result_df)
}
