#' @title integrativeMKT
#'
#' @description \code{itnegrativeMKT()} put details here
#'
#' @details put description here
#'
#' @param daf daf file
#' @param divergence div file
#' @param xlow fit curve
#' @param xhigh fit curv
#' @param seed seed value (optional). No seed by default
#'
#' @return Execute all the MKT extensions
#'
#' @examples 
#' #Load your Derived Allele Frequency file and Divergence file
#' daf<-read.table("/home/jmurga/MKT/Test/data.daf.txt",header=TRUE)
#' div<-read.table("/home/jmurga/MKT/Test/data.divergence.txt",header=TRUE)
#' #Run the function!
#' integrativeMKT(daf,div,0,1)
#'
#' @import knitr 
#' @import utils
#' @import stats
#' @import grid 
#' @import gridExtra
#' @import scales
#' @import reshape2 
#' @import ggplot2
#' @importFrom ggthemes theme_foundation
#' @importFrom cowplot plot_grid
#'
#' @export


integrativeMKT <-  function(daf, divergence, xlow, xhigh, seed) {
  # Print errors if data not correct
  check<-check_input(daf, divergence, 0, 1)
  if(check$data==FALSE)
     stop(check$print_errors)

  fullResutls<-list()
  ## standard MKT
  fullResutls[["standardMKT"]] <- standard(daf,divergence)
  # FWW MKT
  fullResutls[["FWW"]] <- FWW(daf,divergence)
  ## DGRP MKT
  fullResutls[["DGRP"]] <- DGRP(daf,divergence)
  ## Asymptotic MKT
  fullResutls[["Asymptotic"]] <- asymptoticMK(daf,divergence,xlow,xhigh)
  ## iMK
  fullResutls[["iMK"]] <- iMK(daf,divergence,xlow,xhigh)
  
  return(fullResutls)
  
}

