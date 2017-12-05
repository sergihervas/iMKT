#' @title mkt_DGRP
#'
#' @description \code{mktDGRP()} Perform MKT corrected with DGRP method 
#'
#' Load your Derived Allele Frequency file (remember you can use 10 or 20 categories) and Divergence file (it contains de divergents and analysed sites in the synonimous and non synonimous categories)
#'
#'Date = 30/11/2016 Author = Sergi Hervas, Marta Coronado
#'
#' @details The standard McDonald and Kreitman test (MKT) is used to detect the signature of selection at the molecular level. The MKT compares the amount of variation within a species (polymorphism, P) to the divergence (D) between species at two types of sites, one of which is putatively netral and used as the reference to detect selection at the other type of site. In the standard MKT, these sites are synonymous (putatively neutral, 0) and non-synonymous sites (selected sites, i) in a coding region. Under strict neutrality, the ratio of the number of selected and neutral polymorphic sites (Pi/P0) is equal to the ratio of the number of selected and neutral divergence sites (Di/D0).The null hypothesis of neutrality is rejected in a MKT when Di/D0 > Pi/P0. The excess of divergence relative to polymorphism for class i, is interpreted as adaptive selection for a subset of sites i. The fraction of adaptive fixations, alpha.symbol, is estimated from 1-(Pi/P0)(D0/Di). The significance of the test can be assesed with a Fisher exact test. The estimate of alpha.symbol can be easily biased by the segregation of slightly deleterious non-synonymous substitutions. Specifically, slightly deleterious mutations tend to contribute more to polymorphism than to divergence, and thus, lead to an underestimation of alpha. Because adaptive mutations and weakly deleterious selection sct in opposite directions on the MKT, alpha and the fraction of substitutions that are sligholty deleterious, b, will be both underestimated when the two selection regimes occur. To take adaptive and slighlty deleterious mutations mutually into account, Pi, the count off segregatning sites in class i, should be seaprated into the number of neutral variants and the number of weakly deleterious variants, Pi = Pineutral + Pi weak del. Alpha is then estimated as 1-(Pineutral/P0)(D0/Di)
#'
#' @param daf Derived Alelle Frequency File
#' @param divergence Divergence and analyzed sites
#' @return MKT corrected by the DGRP method
#'
#' @examples
#' #Load your Derived Allele Frequency file and Divergence file
#' daf<-read.table("/home/jmurga/MKT/data.daf.txt",header=TRUE)
#' div<-read.table("/home/jmurga/MKT/data.divergence.txt",header=TRUE)
#' #Run the function!
#' mkt_DGRP(daf,div)
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

mkt_DGRP <- function(daf = "Data frame containing the DAF, Pn and Ps", 
                 divergence = "Data frame that contains sites analyzed and divergencen 0fold and 4fold") {
  
  # Shows a message when using the function
  packageStartupMessage("MKT corrected by the DGRP method")
  
  # Declare output data frame
  output <- list()
  
  mkt_tables <-  list()
  fractions <- data.frame(row.names = c("d","f","b"))
  list_cutoffs <- c(0, 0.05, 0.2)
  
  for (cutoff in list_cutoffs) {
    
    daf_below_cutoff <- daf[daf$daf <= cutoff, ] #below DAF
    daf_above_cutoff <- daf[daf$daf > cutoff, ] #over DAF
    
    P0 <- sum(daf$pS) 
    Pi <- sum(daf$pN) 
    
    #Create MKT table 
    mkt_table_standard <- data.frame(Polymorphism = c(sum(daf$pS), sum(daf$pN)), Divergence=c(divergence$D4f,divergence$D0f),row.names = c("Neutral class","Selected class"))
    
    mkt_table <- data.frame(`DAF below cutoff` = c(sum(daf_below_cutoff$pS), sum(daf_below_cutoff$pN)), `DAF above cutoff`=c(sum(daf_above_cutoff$pS), sum(daf_above_cutoff$pN)),row.names = c("Neutral class","Selected class"))
    
    f_neutral <- mkt_table[1,1]/sum(daf$pS)
    Pi_neutral_below_cutoff <- Pi * f_neutral
    
    Pi_wd <- mkt_table[2,1] - Pi_neutral_below_cutoff
    
    Pi_neutral <- round(Pi_neutral_below_cutoff + mkt_table[2,2])
    
    # Estimation of alpha
    alpha <- 1-((Pi_neutral/P0)*(mkt_table_standard[1,2]/mkt_table_standard[2,2]))
    
    # Estimation of b: weakly deleterious
    b <- (Pi_wd/P0)*(divergence$m4f/divergence$m0f)
    
    # ??????? Estimation of y: recently neutral sites
    y <- ((Pi_neutral/P0)-(mkt_table_standard[2,2]/mkt_table_standard[1,2]))*(divergence$m4f/divergence$m0f)
    
    # Estimation of f: neutral sites
    f <- (divergence$m4f*Pi_neutral)/(as.numeric(divergence$m0f)*as.numeric(P0))
    
    # d, strongly deleterious sites
    d <- 1-(f+b)
    
    # Fishers exact test p-value from the MKT
    m <- matrix(c(P0,Pi_neutral,divergence$D4f,divergence$D0f), ncol=2)
    pvalue <- fisher.test(m)$p.value

    # Fractions graph
    fraction <-data.frame(c(d,f,b))
    names(fraction) <- cutoff
    fractions <- cbind(fractions,fraction)
    
    # Store output  
    output[[paste("Cutoff = ",cutoff)]] <- c(cutoff, alpha, pvalue)
    
    
    mkt_table <- knitr::kable(mkt_table,caption = "cutoff")
    
    mkt_tables[[paste("Number of segregating sites by MAF category - Cutoff = ",cutoff)]]  <- mkt_table
    
  } 
  mkt_tables[[paste("MKT standard table")]]  <- mkt_table_standard
  output <- as.data.frame(do.call("rbind",output))
  
  colnames(output) <- c("cutoff", "alpha", "pvalue")
  plot <- ggplot(output, aes(x=as.factor(cutoff), y=alpha, group=1)) +
    geom_line(color="#386cb0") + 
    geom_point(size=2.5, color="#386cb0")+
    theme_Publication() +
    xlab("Cut-off") + ylab(expression(bold(paste("Adaptation (",alpha,")")))) 
  plot
  
  names(output) <- c("alpha.symbol","Fishers exact test P-value")
  
  fractions_melt <- melt(fractions, id.vars=NULL) 
  fractions_melt$Fraction <-  rep(c("d", "f", "b"),3)
  
  plotfraction <- ggplot(fractions_melt) + geom_bar(stat = "identity",aes_string(x="variable", y = "value", fill = "Fraction"),color="black") + coord_flip() + theme_Publication() + ylab(label = "Fraction") + xlab(label = "Cut-off") + scale_fill_manual(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33"), breaks=c("d","f","b"), labels=c(expression(italic("d")), expression(italic("f")),expression(italic("b"))))  + theme(axis.line = element_blank())  + scale_y_discrete(limit=seq(0,1,0.25), expand = c(0, 0))
  plotfraction
  
  plot<-plot_grid(plot,plotfraction, nrow = 2,  labels = c("A", "B"), rel_heights = c(2, 1))
  
  list_output <-list(output,plot,mkt_tables, fractions)
  names(list_output) <- c("Results","Graph", "MKT tables","Fractions")
  

 
  return(list_output)
}

