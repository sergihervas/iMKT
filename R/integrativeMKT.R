#' @title mkt_standard
#'
#' @details \code{mkt_standard()} MKT standard formula
#'Load your Derived Allele Frequency file (remember you can use 10 or 20 categories) and Divergence file (it contains de divergents and analysed sites in the synonimous and non synonimous categories)
#'
#' @description The standard McDonald and Kreitman test (MKT) is used to detect the signature of selection at the molecular level. The MKT compares the amount of variation within a species (polymorphism, P) to the divergence (D) between species at two types of sites, one of which is putatively netral and used as the reference to detect selection at the other type of site. In the standard MKT, these sites are synonymous (putatively neutral, 0) and non-synonymous sites (selected sites, i) in a coding region. Under strict neutrality, the ratio of the number of selected and neutral polymorphic sites (Pi/P0) is equal to the ratio of the number of selected and neutral divergence sites (Di/D0).
# The null hypothesis of neutrality is rejected in a MKT when Di/D0 > Pi/P0. The excess of divergence relative to polymorphism for class i, is interpreted as adaptive selection for a subset of sites i. The fraction of adaptive fixations, alpha.symbol, is estimated from 1-(PN/PS)(Ds/Dn). The significance of the test can be assesed with a Fisher exact test.
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
#' integrativeMKT(daf,div)
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
#' @export
#' 


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

