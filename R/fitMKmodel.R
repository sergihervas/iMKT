#' @title fitMKmodel
#' 
#' @description Core code for two-step nls2() model fit at a given level of precision (res, usually set at 10)
#' Directly retrieved from: https://github.com/MesserLab/asymptoticMK
#'
#' @details Put details here
#' 
#' @param alpha_trimmed correspond to the vectors with alpha and frequency (DAF) values
#' @param f_trimmed  correspond to the vectors with alpha and frequency (DAF) values
#' @param res ASK SERGI
#'
#' @return None
#'
#' @export
 
fitMKmodel <- function(alpha_trimmed, f_trimmed, res) {

  
  ## first fitting using starting values (st)
  mod <- tryCatch({
  
    ## starting values to fit the model  
    st <- expand.grid(const_a=seq(-1,1,length.out=res + 1), const_b=seq(-1,1,length.out=res), const_c=seq(1,10,length.out=res + 1))
    
    ## fitting
    nls2(alpha_trimmed ~ const_a + const_b * exp(-const_c* f_trimmed), start=st, algorithm="brute-force", control=nls.control(maxiter=NROW(st)))
  
  }, error=function(cond) {}) ## return condition of error when unable to fit
  
  ## if mod fails...
  if (length(mod) == 0) { return(NULL) }
  
  ## second fitting, starting from previous fit (mod)
  mod2 <- tryCatch({
    nls2(alpha_trimmed ~ const_a + const_b * exp(-const_c* f_trimmed), start = mod, control=nls.control(maxiter=200))
  
  }, error=function(cond) {}) ## same error handling than the previous step
  
  ## if mod2 fails...
  if (length(mod2) == 0) { return(NULL) }
  
  ## return mod2 if fitted
  return(mod2)
}
