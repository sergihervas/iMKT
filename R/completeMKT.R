#' @title completeMKT
#'
#' @description \code{completeMKT()} put details here
#'
#' @details put description here
#'
#' @param daf data frame containing DAF, Pi and P0 values
#' @param divergence data frame containing divergent and analyzed sites for selected (i) and neutral (0) classes
#' @param xlow lower limit for asymptotic alpha fit
#' @param xhigh higher limit for asymptotic alpha fit
#' @param seed seed value (optional). No seed by default
#'
#' @return Execute all the MKT extensions
#'
#' @examples 
#' completeMKT(mydafdata, mydivergencedata, 0, 0.9)
#' 
#' @import utils
#' @import stats
#'
#' @export

completeMKT <- function(daf, divergence, xlow, xhigh, seed) {
  
  ## Check data
  check <- checkInput(daf, divergence, xlow, xhigh)
  if(check$data == FALSE) {
     stop(check$print_errors) }
  
  ## Check seed
  if(missing(seed)) {
    seed <- NULL
  } else {
    set.seed(seed)
  }

  ## Create output list
  fullResults<-list()
  
  ## Standard MKT
  fullResults[["StandardMK"]] <- standard(daf,divergence)
  
  ## FWW MKT
  fullResults[["FWW"]] <- FWW(daf,divergence)
  
  ## DGRP MKT
  fullResults[["DGRP"]] <- DGRP(daf,divergence)
  
  ## Asymptotic MKT
  fullResults[["Asymptotic"]] <- asymptoticMK(daf,divergence,xlow,xhigh)
  
  ## iMK
  fullResults[["iMK"]] <- iMK(daf,divergence,xlow,xhigh)
  
  ## Output
  return(fullResults)
}

