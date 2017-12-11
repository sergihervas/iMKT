#' @title completeMKT
#'
#' @description \code{completeMKT()} put details here
#'
#' @details put description here
#'
#' @param daf daf file
#' @param divergence div file
#' @param xlow trimming values below this daf threshold
#' @param xhigh trimming values above this daf threshold
#' @param seed seed value (optional). No seed by default
#'
#' @return Execute all the MKT extensions
#'
#' @examples 
#' ## Load your Derived Allele Frequency and Divergence files
#' daf <- mydafdata
#' div <- mydivergencedata
#' ## Run the function
#' # completeMKT(daf, div, 0, 0.9)
#'
#' @import knitr 
#' @import utils
#' @import stats
#' @import grid 
#' @import gridExtra
#' @import scales
#' @import ggplot2
#' @importFrom ggthemes theme_foundation
#' @importFrom cowplot plot_grid
#'
#' @export

completeMKT <- function(daf, divergence, xlow, xhigh, seed) {
  
  ## Check data
  check <- check_input(daf, divergence, xlow, xhigh)
  if(check$data == FALSE) {
     stop(check$print_errors) }

  fullResutls<-list()
  
  ## Standard MKT
  fullResutls[["standardMKT"]] <- standard(daf,divergence)
  
  ## FWW MKT
  fullResutls[["FWW"]] <- FWW(daf,divergence)
  
  ## DGRP MKT
  fullResutls[["DGRP"]] <- DGRP(daf,divergence)
  
  ## Asymptotic MKT
  fullResutls[["Asymptotic"]] <- asymptoticMK(daf,divergence,xlow,xhigh)
  
  ## iMK
  fullResutls[["iMK"]] <- iMK(daf,divergence,xlow,xhigh)
  
  ## Output
  return(fullResutls)
}

