#' @title FWW
#'
#' @description \code{FWW()} Perform MKT corrected with FWW method  
#'
#' @details The standard McDonald and Kreitman test (MKT) is used to detect the signature of selection at the molecular level. The MKT compares the amount of variation within a species (polymorphism, P) to the divergence (D) between species at two types of sites, one of which is putatively netral and used as the reference to detect selection at the other type of site. In the standard MKT, these sites are synonymous (putatively neutral, 0) and non-synonymous sites (selected sites, i) in a coding region. Under strict neutrality, the ratio of the number of selected and neutral polymorphic sites (Pi/P0) is equal to the ratio of the number of selected and neutral divergence sites (Di/D0). The null hypothesis of neutrality is rejected in a MKT when Di/D0 > Pi/P0. The excess of divergence relative to polymorphism for class i, is interpreted as adaptive selection for a subset of sites i. The fraction of adaptive fixations, alpha.symbol, is estimated from 1-(Pi/P0)(Ds/Dn). The significance of the test can be assesed with a Fisher exact test.The estimate of alpha.symbol can be easily biased by the segregation of slightly deleterious non-synonymous substitutions. Specifically, slightly deleterious mutations tend to contribute more to polymorphism than to divergence, and thus, lead to an underestimation of alpha. Bevause they tend to segregate at lower frequencies than do neutral mutations, they can be apartially controled for by removing low frequency polymorphisms from the analysis (Fay et al. 2001). This is known as the FWW method.
#'
#' @param daf data frame containing DAF, Pi and P0 values
#' @param divergence data frame containing divergent and analyzed sites for selected (i) and neutral (0) classes
#' @param list_cutoffs list of cutoffs to use (optional). Default cutoffs are: 0, 0.05, 0.2
#' 
#' @return MKT corrected by the FWW method
#'
#' @examples
#' ## Using default cutoffs
#' # FWW(mydafdata, mydivergencedata)
#' ## Using custom cutoffs
#' # FWW(mydafdata, mydivergencedata, c(0.05, 0.1, 0.15))
#' 
#' @import utils
#' @import stats
#' @import scales
#' @import ggplot2
#' @importFrom ggthemes theme_foundation
#' @importFrom cowplot plot_grid
#'
#' @export

FWW <- function(daf, divergence, list_cutoffs=c(0, 0.05, 0.1), plot=FALSE) {
  
  ## Check data
  check <- check_input(daf, divergence, 0, 1)
  if(check$data == FALSE) {
     stop(check$print_errors) }

  ## Declare outputs
  output <- list()
  mkt_tables <-  list()
  div_metrics <- list()
  div_cutoff <- list()
  
  ## Divergence metrics
  Ka <- divergence$Di/divergence$mi
  Ks <- divergence$D0/divergence$m0
  omega <- Ka/Ks
  div_table <- data.frame(Ka, Ks, omega)
  names(div_table) <- c("Ka", "Ks", "omega")
  
  ## Iterate along cutoffs
  for (cutoff in list_cutoffs) {
    
    daf_remove <-daf[daf$daf > cutoff,] 
  
    ## Create MKT table 
    mkt_table <- data.frame(Polymorphism=c(sum(daf_remove$P0), sum(daf_remove$Pi)), Divergence=c(divergence$D0,divergence$Di), row.names=c("Neutral class","Selected class"))
  
    ## Estimate of alpha
    alpha <- 1-(mkt_table[2,1]/mkt_table[1,1])*(mkt_table[1,2]/mkt_table[2,2])
    
    ## Fisher exact test p-value from the MKT
    pvalue <- fisher.test(mkt_table)$p.value
    
    ## Omega A and Omega D
    omegaA <- omega*alpha
    omegaD <- omega-omegaA
    
    ## Store output  
    output[[paste("Cutoff = ",cutoff)]] <- c(cutoff, alpha,pvalue)
    div_cutoff[[paste("Cutoff = ",cutoff)]] <- c(cutoff, omegaA, omegaD)
    mkt_tables[[paste("Cutoff = ",cutoff)]]  <- kable(mkt_table,caption = "cutoff")
  }
  
  ## Output format
  output <- as.data.frame(do.call("rbind",output))
  colnames(output) <- c("cutoff", "alpha", "pvalue")
 
  ## Divergence metrics
  div_cutoff <- as.data.frame(do.call("rbind",div_cutoff))
  names(div_cutoff) <- c("cutoff", "omegaA", "omegaD")
  
  if(plot == TRUE) {
  ## Cut-offs graphs
  plot <- ggplot(output, aes(x=as.factor(cutoff), y=alpha, group=1)) +
    geom_line(color="#386cb0") + 
    geom_point(size=2.5, color="#386cb0")+
    theme_Publication() +
    xlab("Cut-off") + ylab(expression(bold(paste("Adaptation (",alpha,")"))))
  
  ## Re-format outputs
  output <- output[,c(2,3)]
  names(output) <- c("alpha.symbol","Fishers exact test P-value")
  div_cutoff <- div_cutoff[,c(2,3)]
  colnames(div_cutoff) <- c("omegaA.symbol", "omegaD.symbol")
  div_metrics <- list(div_table, div_cutoff)
  names(div_metrics) <- c("Global metrics", "Estimates by cutoff")
  
  ## Return list
  list_output <-list(output, plot, div_metrics, mkt_tables)
  names(list_output) <- c("Results", "Graph", "Divergence metrics", "MKT tables")
  } else if (plot==FALSE)
  {
    ## Re-format outputs
    output <- output[,c(2,3)]
    names(output) <- c("alpha.symbol","Fishers exact test P-value")
    div_cutoff <- div_cutoff[,c(2,3)]
    colnames(div_cutoff) <- c("omegaA.symbol", "omegaD.symbol")
    div_metrics <- list(div_table, div_cutoff)
    names(div_metrics) <- c("Global metrics", "Estimates by cutoff")
    
    ## Return list
    list_output <-list(output, div_metrics, mkt_tables)
    names(list_output) <- c("Results", "Divergence metrics", "MKT tables")
  }
  return(list_output)
}
