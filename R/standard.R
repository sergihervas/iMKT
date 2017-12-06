#' @title mkt_standard
#'
#' @details \code{mkt_standard()} MKT standard formula
#'Load your Derived Allele Frequency file (remember you can use 10 or 20 categories) and Divergence file (it contains de divergents and analysed sites in the synonimous and non synonimous categories)
#'
#' @description The standard McDonald and Kreitman test (MKT) is used to detect the signature of selection at the molecular level. The MKT compares the amount of variation within a species (polymorphism, P) to the divergence (D) between species at two types of sites, one of which is putatively netral and used as the reference to detect selection at the other type of site. In the standard MKT, these sites are synonymous (putatively neutral, 0) and non-synonymous sites (selected sites, i) in a coding region. Under strict neutrality, the ratio of the number of selected and neutral polymorphic sites (Pi/P0) is equal to the ratio of the number of selected and neutral divergence sites (Di/D0).
# The null hypothesis of neutrality is rejected in a MKT when Di/D0 > Pi/P0. The excess of divergence relative to polymorphism for class i, is interpreted as adaptive selection for a subset of sites i. The fraction of adaptive fixations, alpha.symbol, is estimated from 1-(Pi/P0)(Ds/Dn). The significance of the test can be assesed with a Fisher exact test.
#'
#' @param daf daf file
#' @param divergence divergence file
#' 
#' @return Standard McDonald and Kreitman Test
#'
#' @examples 
#' #Load your Derived Allele Frequency file and Divergence file
#' daf<-read.table("/home/jmurga/MKT/Test/data.daf.txt",header=TRUE)
#' div<-read.table("/home/jmurga/MKT/Test/data.divergence.txt",header=TRUE)
#' #Run the function!
#' standard(daf,div)
#' @import knitr 
#' @import utils
#' @import stats
#'
#' @export 


standard <- function(daf = "Data frame containing the DAF, Pi and P0", 
                     divergence = "Data frame that contains sites analyzed and divergencen 0fold and 4fold") {
  # Print errors if data not correct
 # check<-check_input(daf, divergence, 0, 1)
  #  if(check$data==FALSE)
  # stop(check$print_errors)
  # Shows a message when using the function
  # packageStartupMessage("MKT standard function")
  check<-check_input(daf, divergence, 0, 1)
  
  # Declare output data frame
  output <- data.frame(alpha = numeric(0), pvalue = integer(0))
  #out <- NULL
  
  #Create MKT table 
  mkt_table <- data.frame(Polymorphism = c(sum(daf$P0), sum(daf$Pi)), Divergence=c(divergence$D0,divergence$Di),row.names = c("Neutral class","Selected class"))
  
  # Estimation of alpha
  alpha <- 1-(mkt_table[2,1]/mkt_table[1,1])*(mkt_table[1,2]/mkt_table[2,2])
    
  # Fisher's exact test p-value from the MKT
  pvalue <- fisher.test(mkt_table)$p.value
  
  # Store output  
  output <- list(`alpha.symbol`=alpha, 
                 `Fishers exact test P-value`= pvalue, 
                 `MKT table` = knitr::kable(mkt_table))
  # Return output in list format
  return(output)
}

