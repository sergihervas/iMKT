#' MKT standard formula
#'
#' Date = 30/11/2016
#' Author = Sergi Hervás, Marta Coronado
#'
#' The standard McDonald and Kreitman test (MKT) is used to detect the signature of selection at the molecular level. The MKT compares the amount of variation within a species (polymorphism, P) to the divergence (D) between species at two types of sites, one of which is putatively netral and used as the reference to detect selection at the other type of site. In the standard MKT, these sites are synonymous (putatively neutral, 0) and non-synonymous sites (selected sites, i) in a coding region. Under strict neutrality, the ratio of the number of selected and neutral polymorphic sites (Pi/P0) is equal to the ratio of the number of selected and neutral divergence sites (Di/D0).
# The null hypothesis of neutrality is rejected in a MKT when Di/D0 > Pi/P0. The excess of divergence relative to polymorphism for class i, is interpreted as adaptive selection for a subset of sites i. The fraction of adaptive fixations, α, is estimated from 1-(PN/PS)(Ds/Dn). The significance of the test can be assesed with a Fisher exact test.
#'
#'
#' @param x dad file
#' @param y divergence file
#' 
#' @return None
#'
#' @examples
#' mkt_DGRP(x,y)
#'
#' Standar MKT
#'$α
#'[1] 0.07371226
#'
#'$`Fisher's exact test P-value`
#'[1] 2.17533e-14
#'
#'$`MKT table`
#'
#'
#'|               | Polymorphism| Divergence|
#'|:--------------|------------:|----------:|
#'|Neutral class  |        32308|      52537|
#'|Selected class |        31125|      54641|
#'#'
#' @export
#' 


mkt_standard <- function(daf = "Data frame containing the DAF, Pn and Ps", 
                     divergence = "Data frame that contains the number of non-synonymous sites analyzed (m0f), number of synonymous sites analyzed (m4f), total number of non-synonymous divergent sites (D0f) and number of synonymous divergent sites (D4f)") {
  # Print errors if data not correct
 # check<-check_input(daf, divergence, 0, 1)
  #  if(check$data==FALSE)
  # stop(check$print_errors)
  # Shows a message when using the function
  # packageStartupMessage("MKT standard function")
  check<-check_input(x, y, 0, 1)
  
  # Declare output data frame
  output <- data.frame(alpha = numeric(0), pvalue = integer(0))
  #out <- NULL
  
  #Create MKT table 
  mkt_table <- data.frame(Polymorphism = c(sum(daf$pS), sum(daf$pN)), Divergence=c(divergence$D4f,divergence$D0f),row.names = c("Neutral class","Selected class"))
  
  # Estimation of alpha
  alpha <- 1-(mkt_table[2,1]/mkt_table[1,1])*(mkt_table[1,2]/mkt_table[2,2])
    
  # Fisher's exact test p-value from the MKT
  pvalue <- fisher.test(mkt_table)$p.value
  
  # Store output  
  output <- list(`α`=alpha, 
                 `Fishers exact test P-value`= pvalue, 
                 `MKT table` = knitr::kable(mkt_table))
  # Return output in list format
  return(output)
}

