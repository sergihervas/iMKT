#' @title DGRP
#'
#' @description \code{DGRP()} Perform MKT corrected with DGRP method 
#'
#' @details The standard McDonald and Kreitman test (MKT) is used to detect the signature of selection at the molecular level. The MKT compares the amount of variation within a species (polymorphism, P) to the divergence (D) between species at two types of sites, one of which is putatively netral and used as the reference to detect selection at the other type of site. In the standard MKT, these sites are synonymous (putatively neutral, 0) and non-synonymous sites (selected sites, i) in a coding region. Under strict neutrality, the ratio of the number of selected and neutral polymorphic sites (Pi/P0) is equal to the ratio of the number of selected and neutral divergence sites (Di/D0).The null hypothesis of neutrality is rejected in a MKT when Di/D0 > Pi/P0. The excess of divergence relative to polymorphism for class i, is interpreted as adaptive selection for a subset of sites i. The fraction of adaptive fixations, alpha.symbol, is estimated from 1-(Pi/P0)(D0/Di). The significance of the test can be assesed with a Fisher exact test. The estimate of alpha.symbol can be easily biased by the segregation of slightly deleterious non-synonymous substitutions. Specifically, slightly deleterious mutations tend to contribute more to polymorphism than to divergence, and thus, lead to an underestimation of alpha. Because adaptive mutations and weakly deleterious selection sct in opposite directions on the MKT, alpha and the fraction of substitutions that are sligholty deleterious, b, will be both underestimated when the two selection regimes occur. To take adaptive and slighlty deleterious mutations mutually into account, Pi, the count off segregatning sites in class i, should be seaprated into the number of neutral variants and the number of weakly deleterious variants, Pi = Pineutral + Pi weak del. Alpha is then estimated as 1-(Pineutral/P0)(D0/Di)
#'
#' @param daf data frame containing DAF, Pi and P0 values
#' @param divergence data frame containing divergent and analyzed sites for selected (i) and neutral (0) classes
#' @param list_cutoffs list of cutoffs to use (optional). Default cutoffs are: 0, 0.05, 0.1
#' @param plot report plot (optional). Default is FALSE
#' 
#' @return MKT corrected by the DGRP method. List with alpha results, graph, divergence metrics, MKT tables and negative selection fractions
#'
#' @examples
#' ## Using default cutoffs
#'  DGRP(mydafdata, mydivergencedata)
#' ## Using custom cutoffs
#' DGRP(mydafdata, mydivergencedata, c(0.05, 0.1, 0.15))
#'
#' @import utils
#' @import stats
#' @import ggplot2
#' @importFrom ggthemes theme_foundation
#' @importFrom cowplot plot_grid
#' @importFrom reshape2 melt 
#' @importFrom knitr kable 
#'
#' @export


####################################################
################# MKT-FWW function #################
####################################################

DGRP <- function(daf, divergence, list_cutoffs=c(0, 0.05, 0.1), plot=FALSE) {
  
  ## Check data
  check <- checkInput(daf, divergence, 0, 1)
  if(check$data == FALSE) {
     stop(check$print_errors) }

  ## Declare output lists and data frames
  output <- list()
  mkt_tables <- list()
  fractions <- data.frame(row.names = c("d","f","b"))
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
    
    daf_below_cutoff <- daf[daf$daf <= cutoff, ] ## below DAF
    daf_above_cutoff <- daf[daf$daf > cutoff, ] ## over DAF 
    P0 <- sum(daf$P0) 
    Pi <- sum(daf$Pi) 

    ## Create MKT table 
    mkt_table_standard <- data.frame(Polymorphism = c(sum(daf$P0), sum(daf$Pi)), Divergence=c(divergence$D0,divergence$Di),row.names = c("Neutral class","Selected class"))
    mkt_table <- data.frame(`DAF below cutoff` = c(sum(daf_below_cutoff$P0), sum(daf_below_cutoff$Pi)), `DAF above cutoff`=c(sum(daf_above_cutoff$P0), sum(daf_above_cutoff$Pi)),row.names = c("Neutral class","Selected class"))
    
    ## Estimate fractions
    f_neutral <- mkt_table[1,1]/sum(daf$P0)
    Pi_neutral_below_cutoff <- Pi * f_neutral
    Pi_wd <- mkt_table[2,1] - Pi_neutral_below_cutoff
    Pi_neutral <- round(Pi_neutral_below_cutoff + mkt_table[2,2])
    
    ## Estimation of alpha
    alpha <- 1-((Pi_neutral/P0)*(mkt_table_standard[1,2]/mkt_table_standard[2,2]))
    
    ## Estimation of b: weakly deleterious
    b <- (Pi_wd/P0)*(divergence$m0/divergence$mi)
    
    ## Estimation of f: neutral sites
    f <- (divergence$m0*Pi_neutral)/(as.numeric(divergence$mi)*as.numeric(P0))
    
    ## Estimation of d, strongly deleterious sites
    d <- 1-(f+b)
    
    ## Fisher exact test p-value from the MKT
    m <- matrix(c(P0,Pi_neutral,divergence$D0,divergence$Di), ncol=2)
    pvalue <- fisher.test(m)$p.value
    
    ## Omega A and Omega D
    omegaA <- omega*alpha
    omegaD <- omega-omegaA
    
    ## Fractions
    fraction <-data.frame(c(d,f,b))
    names(fraction) <- cutoff
    fractions <- cbind(fractions,fraction)
    
    ## Store output  
    output[[paste("Cutoff = ",cutoff)]] <- c(cutoff, alpha, pvalue)
    div_cutoff[[paste("Cutoff = ",cutoff)]] <- c(cutoff, omegaA, omegaD)
    mkt_table <- kable(mkt_table, caption = "cutoff")
    mkt_tables[[paste("Number of segregating sites by DAF category - Cutoff = ",cutoff)]]  <- mkt_table
  } 

  ## MKT tables
  mkt_tables[[paste("MKT standard table")]]  <- kable(mkt_table_standard)
  
  ## Results table
  output <- as.data.frame(do.call("rbind",output))
  colnames(output) <- c("cutoff", "alpha", "pvalue")
  
  ## Divergence metrics
  div_cutoff <- as.data.frame(do.call("rbind",div_cutoff))
  names(div_cutoff) <- c("cutoff", "omegaA", "omegaD")
  
  ## Render plot
  if (plot == TRUE) {
    
    ## Cut-offs graph
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
  
    ## Melt fractions data
    fractions_melt <- melt(fractions, id.vars=NULL) 
    fractions_melt$Fraction <-  rep(c("d", "f", "b"),length(fractions_melt$variable)/3)
  
    ## Fractions graph
    plotfraction <- ggplot(fractions_melt) + geom_bar(stat="identity", aes_string(x="variable", y="value", fill="Fraction"), color="black") +
      coord_flip() + themePublication() + ylab(label="Fraction") + xlab(label="Cut-off") +
      scale_fill_manual(values=c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33"), breaks=c("f","d","b"), labels=c(expression(italic("f")),expression(italic("d")),expression(italic("b")))) +
      theme(axis.line=element_blank())  + scale_y_discrete(limit=seq(0,1,0.25), expand=c(0,0))
  
    plot <- plot_grid(plot, plotfraction, nrow=2, labels=c("A","B"), rel_heights=c(2,1))
  
    ## Return list output  
    list_output <- list(output, plot, div_metrics, mkt_tables, fractions)
    names(list_output) <- c("Results","Graph", "Divergence metrics", "MKT tables","Fractions")
  
  ## If no plot to render
  } else if (plot==FALSE) {
      ## Re-format outputs
      output <- output[,c(2,3)]
      names(output) <- c("alpha.symbol","Fishers exact test P-value")
      div_cutoff <- div_cutoff[,c(2,3)]
      colnames(div_cutoff) <- c("omegaA.symbol", "omegaD.symbol")
      div_metrics <- list(div_table, div_cutoff)
      names(div_metrics) <- c("Global metrics", "Estimates by cutoff")
    
      ## Melt fractions data
      fractions_melt <- melt(fractions, id.vars=NULL) 
      fractions_melt$Fraction <-  rep(c("d", "f", "b"),length(fractions_melt$variable)/3)
      
      ## Return list output  
      list_output <- list(output, div_metrics, mkt_tables, fractions)
      names(list_output) <- c("Results", "Divergence metrics", "MKT tables","Fractions")
  }
  
  ## Return output
  return(list_output)
}

