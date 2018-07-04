#' @title Complete MK methodologies
#'
#' @description MKT calculation using all methodologies included in the package: standardMKT, FWW, DGRP, asymptoticMKT, iMKT.
#' 
#' @details Perform all MKT derived methodologies at once using the same input data and parameters.
#'
#' @param daf data frame containing DAF, Pi and P0 values
#' @param divergence data frame containing divergent and analyzed sites for selected (i) and neutral (0) classes
#' @param xlow lower limit for asymptotic alpha fit
#' @param xhigh higher limit for asymptotic alpha fit
#' @param seed seed value (optional). No seed by default
#'
#' @return List with the diverse MKT results: standardMKT, FWW, DGRP, asymptoticMKT, iMKT
#'
#' @examples 
#' completeMKT(myDafData, myDivergenceData, xlow=0, xhigh=0.9)
#' 
#' @import utils
#' @import stats
#'
#' @keywords MKT
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
  fullResults <- list()
  
  ## Standard MKT
  fullResults[["StandardMKT"]] <- standardMKT(daf,divergence)
  
  ## FWW MKT
  fullResults[["FWW"]] <- FWW(daf,divergence)
  
  ## DGRP MKT
  fullResults[["DGRP"]] <- DGRP(daf,divergence)
  
  ## Asymptotic MKT
  fullResults[["Asymptotic"]] <- asymptoticMKT(daf,divergence,xlow,xhigh)
  
  ## iMKT
  fullResults[["iMKT"]] <- iMKT(daf,divergence,xlow,xhigh)
  
  ## Output
  return(fullResults)
}

