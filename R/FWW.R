#' @title FWW correction method
#'
#' @description MKT calculation corrected using FWW method (Fay et al. 2001 Genetics).
#'
#' @details In the standard McDonald and Kreitman test, the estimate of adaptive evolution (alpha) can be easily biased by the segregation of slightly deleterious non-synonymous substitutions. Specifically, slightly deleterious mutations contribute more to polymorphism than they do to divergence, and thus, lead to an underestimation of alpha. Because they tend to segregate at lower frequencies than do neutral mutations, they can be partially controled by removing low frequency polymorphisms from the analysis. This is known as the FWW method.
#'
#' @param daf data frame containing DAF, Pi and P0 values
#' @param divergence data frame containing divergent and analyzed sites for selected (i) and neutral (0) classes
#' @param listCutoffs list of cutoffs to use (optional). Default cutoffs are: 0, 0.05, 0.1
#' @param plot report plot (optional). Default is FALSE
#' 
#' @return MKT corrected by the FWW method. List with alpha results, graph (optional), divergence metrics, MKT tables and negative selection fractions
#'
#' @examples
#' ## Using default cutoffs
#' FWW(myDafData, myDivergenceData)
#' ## Using custom cutoffs and rendering plot
#' FWW(myDafData, myDivergenceData, c(0.05, 0.1, 0.15), plot=TRUE)
#' 
#' @import utils
#' @import stats
#' @import ggplot2
#' @importFrom ggthemes theme_foundation
#' @importFrom cowplot plot_grid
#'
#' @keywords MKT
#' @export

FWW <- function(daf, divergence, listCutoffs=c(0,0.05,0.1), plot=FALSE) {
  
  ## Check data
  check <- checkInput(daf, divergence, 0, 1)
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
  for (cutoff in listCutoffs) {
    
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
  
  ## Perform plot
  if(plot == TRUE) {
    ## Cut-offs graphs
    plot <- ggplot(output, aes(x=as.factor(cutoff), y=alpha, group=1)) +
      geom_line(color="#386cb0") + 
      geom_point(size=2.5, color="#386cb0")+
      themePublication() +
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
  
  ## If no plot to perform  
  } else if (plot == FALSE) {
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
  
  ## Return output
  return(list_output)
}
