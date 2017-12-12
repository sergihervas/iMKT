#' @title fitMKmodel
#' 
#' @description \code{fitMKmodel} Core code for two-step nls2() model fit at a given level of precision (res, usually set at 10)
#' Adapted from: https://github.com/MesserLab/asymptoticMK
#'
#' @details nls2() model fit used to estimate asymptotic alpha value. This function is adapted from the code developed in "Haller BC, Messer PW. asymptoticMK: A Web-Based Tool for the Asymptotic McDonald-Kreitman Test. G3 (Bethesda). 2017 May 5;7(5):1569-1575".
#' 
#' @param alpha_trimmed correspond to the vectors with alpha and frequency (DAF) values
#' @param f_trimmed  correspond to the vectors with alpha and frequency (DAF) values
#' @param res refers to the level of precision, residues to use (10 or 20)
#' @return None
#'
#' @importFrom nls2 nls2
#'
#' @export
 
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
