#' @title asymptoticMK
#' 
#' @description \code{asymptoticMK} developed in "Haller BC, Messer PW. asymptoticMK: A Web-Based Tool for the Asymptotic McDonald-Kreitman Test. G3 (Bethesda). 2017 May 5;7(5):1569-1575".
#' Adapted from: http://github.com/MesserLab/asymptoticMK
#'
#' @details The standard McDonald and Kreitman test (MKT) is used to detect the signature of selection at the molecular level. The MKT compares the amount of variation within a species (polymorphism, P) to the divergence (D) between species at two types of sites, one of which is putatively netral and used as the reference to detect selection at the other type of site. In the standard MKT, these sites are synonymous (putatively neutral, 0) and non-synonymous sites (selected sites, i) in a coding region. Under strict neutrality, the ratio of the number of selected and neutral polymorphic sites (Pi/P0) is equal to the ratio of the number of selected and neutral divergence sites (Di/D0).
# The null hypothesis of neutrality is rejected in a MKT when Di/D0 > Pi/P0. The excess of divergence relative to polymorphism for class i, is interpreted as adaptive selection for a subset of sites i. The fraction of adaptive fixations, α, is estimated from 1-(Pi/P0)(Ds/Dn). The significance of the test can be assesed with a Fisher exact test.
#'
#' @param daf data frame containing DAF, Pi and P0 values
#' @param divergence data frame containing divergent and analyzed sites for selected (i) and neutral (0) classes
#' @param xlow trimming values below this daf threshold
#' @param xhigh trimming values above this daf threshold
#' @param seed seed value (optional). No seed by default
#'
#' @return Estimation of alpha asymptotic value and details of the model fit
#'
#' @examples
#' # asymptoticMK(mydafdata, mydivergencedata, 0, 0.9)
#'
#' @import utils
#' @import stats
#' @importFrom MASS mvrnorm
#' @importFrom nls2 nls2
#'
#' @export

asymptoticMK <- function(daf, divergence, xlow, xhigh, seed) {
  
  ## Check data
  check <- check_input(daf, divergence, xlow, xhigh)
  if(check$data == FALSE) {
    stop(check$print_errors) }

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
  mod1 <- fitMKmodel(alpha_trimmed, f_trimmed, 10)
  
  ## If mod1 did not work, try a deeper scan for a decent fit (res=20)
  if (length(mod1) == 0) {
    mod1 <- fitMKmodel(alpha_trimmed, f_trimmed, 20)
  } 
  
  tryCatch({
    mod2 <- lm(alpha_trimmed ~ f_trimmed)
  }, error=function(cond) {})
  
  ## Compute confidence intervals of alpha using predictNLS 
  ci_pred <- predictNLS(mod1, newdata=data.frame(f_trimmed=1.0))
  alpha_1_low <- ci_pred[6]
  alpha_1_high <- ci_pred[7]
  
  ## Preparation of ouput (alpha asym, a, b, c)
  alpha_1_est <- predict(mod1, newdata=data.frame(f_trimmed=1.0))
  const_a <- coef(mod1)["const_a"]
  const_b <- coef(mod1)["const_b"]
  const_c <- coef(mod1)["const_c"]
  
  ## Output table
  result_df <- data.frame(model="exponential", a=const_a, b=const_b, c=const_c, alpha_asymptotic=alpha_1_est, CI_low=ci_pred[6], CI_high=ci_pred[7], alpha_original=alpha_nonasymp, row.names=NULL)
  return(result_df)
}
