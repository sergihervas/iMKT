#' Get a CI using Monte Carlo simulation based upon a fitted model.  This is necessary because
# getting confidence intervals for non-linear models is a complicated business, apparently. Authors thanks to Andrej-Nikolai Spiess (http://www.dr-spiess.de) for this code.
#' Directly retrieved from: https://github.com/MesserLab/asymptoticMK
#' Date = 30/11/2016
#' Author = REFER AUTHROS? CONTACT WITH PETROV LAB
#'
#'
#' @param alpha_trimmed correspond to the vectors with alpha and frequency (DAF) values
#' @param f_trimmed  correspond to the vectors with alpha and frequency (DAF) values
#' @param res ASK SERGI
#'
#' @return None
#'
#' @examples
#' 
#' @export
#' 


# Get a CI using Monte Carlo simulation based upon a fitted model.  This is necessary because
# getting confidence intervals for non-linear models is a complicated business, apparently.
# Thanks to Andrej-Nikolai Spiess (http://www.dr-spiess.de) for this code.
# See: https://www.r-bloggers.com/predictnls-part-1-monte-carlo-simulation-confidence-intervals-for-nls-models/
# Or, if that link goes stale: http://stats.stackexchange.com/a/251501/141766


predictNLS <- function(object, newdata, level = 0.95, nsim = 10000) {
  
  ## define required packages
  require(MASS, quietly = TRUE)
   
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
    #counter(i)		# show a counter for lengthy fits; commented out to reduce noise here...
     
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
   
  #cat("\n")	# commented out along with the call to counter() above
   
  return(outMAT)  
}
