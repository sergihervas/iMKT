#' @title FWW
#'
#' @description \code{FWW()} Perform MKT corrected with FWW method  
#'
#'Load your Derived Allele Frequency file (remember you can use 10 or 20 categories) and Divergence file (it contains de divergents and analysed sites in the synonimous and non synonimous categories)
#'
#'@details The standard McDonald and Kreitman test (MKT) is used to detect the signature of selection at the molecular level. The MKT compares the amount of variation within a species (polymorphism, P) to the divergence (D) between species at two types of sites, one of which is putatively netral and used as the reference to detect selection at the other type of site. In the standard MKT, these sites are synonymous (putatively neutral, 0) and non-synonymous sites (selected sites, i) in a coding region. Under strict neutrality, the ratio of the number of selected and neutral polymorphic sites (Pi/P0) is equal to the ratio of the number of selected and neutral divergence sites (Di/D0). The null hypothesis of neutrality is rejected in a MKT when Di/D0 > Pi/P0. The excess of divergence relative to polymorphism for class i, is interpreted as adaptive selection for a subset of sites i. The fraction of adaptive fixations, alpha.symbol, is estimated from 1-(Pi/P0)(Ds/Dn). The significance of the test can be assesed with a Fisher exact test.The estimate of alpha.symbol can be easily biased by the segregation of slightly deleterious non-synonymous substitutions. Specifically, slightly deleterious mutations tend to contribute more to polymorphism than to divergence, and thus, lead to an underestimation of alpha. Bevause they tend to segregate at lower frequencies than do neutral mutations, they can be apartially controled for by removing low frequency polymorphisms from the analysis (Fay et al. 2001). This is known as the FWW method.
#'
#' @param daf daf file
#' @param divergence divergence file
#' @param list_cutoffs Default cutoff c(0, 0.05, 0.2)
#' 
#' @return MKT corrected by the FWW method
#'
#' @examples
#' #Load your Derived Allele Frequency file and Divergence file
#' #daf<-read.table("/home/jmurga/MKT/Test/data.daf.txt",header=TRUE)
#' #div<-read.table("/home/jmurga/MKT/Test/data.divergence.txt",header=TRUE)
#' #Run the function!
#' #FWW(daf,div)
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



####################################################
################# MKT-FWW function #################
####################################################

FWW<- function(daf = "Data frame containing the DAF, Pi and P0", 
                     divergence = "Data frame that contains sites analyzed and divergencen 0fold and 4fold",list_cutoffs=c(0, 0.05, 0.1)) {
  # Shows a message when using the function
  packageStartupMessage("MKT corrected by the FWW method")
  
  # Declare output data frame
  output <- list()
  
  mkt_tables <-  list()
  # list_cutoffs <- c(0, 0.05, 0.1)
  
  for (cutoff in list_cutoffs) {
  
    daf_remove <-daf[daf$daf > cutoff,] 
  
    #Create MKT table 
    mkt_table <- data.frame(Polymorphism = c(sum(daf_remove$P0), sum(daf_remove$Pi)), Divergence=c(divergence$D0,divergence$Di),row.names = c("Neutral class","Selected class"))
  
    # Estimation of alpha
    alpha <- 1-(mkt_table[2,1]/mkt_table[1,1])*(mkt_table[1,2]/mkt_table[2,2])
    
    # Fisher's exact test p-value from the MKT
    pvalue <- fisher.test(mkt_table)$p.value
    
    # Store output  
    output[[paste("Cutoff = ",cutoff)]] <- c(cutoff, alpha,pvalue)
    
    
    mkt_tables[[paste("Cutoff = ",cutoff)]]  <- knitr::kable(mkt_table,caption = "cutoff")
  }
  
  output <- as.data.frame(do.call("rbind",output))
  colnames(output) <- c("cutoff", "alpha", "pvalue")
  plot <- ggplot(output, aes(x=as.factor(cutoff), y=alpha, group=1)) +
    geom_line(color="#386cb0") + 
    geom_point(size=2.5, color="#386cb0")+
    theme_Publication() +
    xlab("Cut-off") + ylab(expression(bold(paste("Adaptation (",alpha,")"))))
  plot
  
  names(output) <- c("alpha.symbol","Fisher's exact test P-value")
  
  list_output <-list(output,plot,mkt_tables)
  names(list_output) <- c("Results","Graph", "MKT tables")

    return(list_output)
}

