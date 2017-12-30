#' @title Complete MK methodologies
#'
#' @description MKT calculation using all methodologies included in the package: standardMK, FWW, DGRP, asymptoticMK, iMK.
#' 
#' @details Perform all MKT derived methodologies at the same time using the same input parameters.
#'
#' @param daf data frame containing DAF, Pi and P0 values
#' @param divergence data frame containing divergent and analyzed sites for selected (i) and neutral (0) classes
#' @param xlow lower limit for asymptotic alpha fit
#' @param xhigh higher limit for asymptotic alpha fit
#' @param seed seed value (optional). No seed by default
#'
#' @return List with all MKT results: standardMK, FWW, DGRP, asymptoticMK, iMK
#'
#' @examples 
#' completeMK(myDafData, myDivergenceData, 0, 0.9)
#' 
#' @import utils
#' @import stats
#'
#' @export

completeMK <- function(daf, divergence, xlow, xhigh, seed) {
  
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
  fullResults <- list()
  
  ## Standard MKT
  fullResults[["StandardMK"]] <- standardMK(daf,divergence)
  
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

