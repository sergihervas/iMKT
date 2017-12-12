#' @title mkt_standard
#'
#' @description \code{mkt_standard()} MKT standard formula
#'
#' @details The standard McDonald and Kreitman test (MKT) is used to detect the signature of selection at the molecular level. The MKT compares the amount of variation within a species (polymorphism, P) to the divergence (D) between species at two types of sites, one of which is putatively netral and used as the reference to detect selection at the other type of site. In the standard MKT, these sites are synonymous (putatively neutral, 0) and non-synonymous sites (selected sites, i) in a coding region. Under strict neutrality, the ratio of the number of selected and neutral polymorphic sites (Pi/P0) is equal to the ratio of the number of selected and neutral divergence sites (Di/D0).
# The null hypothesis of neutrality is rejected in a MKT when Di/D0 > Pi/P0. The excess of divergence relative to polymorphism for class i, is interpreted as adaptive selection for a subset of sites i. The fraction of adaptive fixations, alpha.symbol, is estimated from 1-(Pi/P0)(Ds/Dn). The significance of the test can be assesed with a Fisher exact test.
#'
#' @param daf data frame containing DAF, Pi and P0 values
#' @param divergence data frame containing divergent and analyzed sites for selected (i) and neutral (0) classes
#' 
#' @return Standard McDonald and Kreitman Test
#'
#' @examples 
#' # standard(mydafdata, mydivergencedata)
#' 
#' @import utils
#' @import stats
#' @importFrom knitr kable 
#'
#' @export

standard <- function(daf, divergence) {
  
  ## Check data
  check <-check_input(daf, divergence, 0, 1)
    if(check$data == FALSE) {
     stop(check$print_errors) }

  ## Declare output data frame
  output <- data.frame(alpha = numeric(0), pvalue = integer(0))
  
  ## Create MKT table 
  mkt_table <- data.frame(Polymorphism = c(sum(daf$P0), sum(daf$Pi)), Divergence=c(divergence$D0,divergence$Di), row.names = c("Neutral class","Selected class"))
  
  ## Estimation of alpha and fisher exact test p-value
  alpha <- 1-(mkt_table[2,1]/mkt_table[1,1])*(mkt_table[1,2]/mkt_table[2,2])
  pvalue <- fisher.test(mkt_table)$p.value
  
  ## Ka, Ks, omega, omegaA, omegaD
  Ka <- divergence$Di/divergence$mi
  Ks <- divergence$D0/divergence$m0
  omega <- Ka/Ks
  omegaA <- omega*alpha
  omegaD <- omega-omegaA
  divergence_metrics <- data.frame(Ka, Ks, omega, omegaA, omegaD)
  
  ## Output  
  output <- list(`alpha.symbol`=alpha, 
                 `Fishers exact test P-value`= pvalue, 
                 `MKT table`= kable(mkt_table),
                 `Divergence metrics` = kable(divergence_metrics))
  return(output)
}

