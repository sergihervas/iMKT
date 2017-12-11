#' @title KaKs
#'
#' @description KaKs
#'
#' @details KaKs
#'
#' @param daf daf file
#' @param divergence divergence file
#' 
#' @return Standard McDonald and Kreitman Test
#'
#' @examples 
#' ## Load your Derived Allele Frequency and Divergence files
#' daf <- mydafdata
#' div <- mydivergencedata
#' ## Run the function
#' # divergenceParameters(daf, div)
#' 
#' @import knitr 
#' @import utils
#' @import stats
#'
#' @export 


divergenceParameters <- function(daf = "Data frame containing the DAF, Pi and P0",
                     divergence = "Data frame that contains sites analyzed and divergencen 0fold and 4fold") {
  # Print errors if data not correct
 # check<-check_input(daf, divergence, 0, 1)
  #  if(check$data==FALSE)
  # stop(check$print_errors)
  # Shows a message when using the function
  # packageStartupMessage("MKT standard function")
  check<-check_input(daf, divergence, 0, 1)
  
  # Declare output data frame
  output <- data.frame(KaKs = numeric(0))
  #out <- NULL
  
  #Create MKT table 
  kaks <- (sum(daf$Pi)/divergence$mi)/(sum(daf$P0)/divergence$m0)
  # Estimation of alpha
  output<-kaks
# 
  return(output)
}

