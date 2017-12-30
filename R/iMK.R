#' @title integrative MKT method
#' 
#' @description iMK: MKT using asymptoticMK method and estimation of negative selection fractions (d, b, f)
#'
#' @details The integrative MKT (iMKT) allows the estimation of the rate of adaptive evolution (alpha) and the diverse negative selection regimens. iMKT uses asymptotic MK method (Messer and Petrov 2012 PNAS; Haller and Messer 2017 G3) to estimate alpha and the diverse negative selection fractions (d: strongly deleterious, b: weakly deleterious, f: neutral), based on the assumption that weakly deleterious mutations usually do not reach high allele frequencies and therefore, produce the underestimation of alpha at low DAF categories. The fraction of strongly deleterious mutations is estimated as the difference between neutral (0) and selected (i) polymorphic sites relative to the number of analyzed sites: d = 1 - (P0/m0 / Pi/mi). The fraction of weakly deleterious sites (b) corresponds to the relative proportion of selected polymorphic sites that cause the underestimation of alpha at low DAF categories. Finally, the fraction of neutral sites (f) is estimated as: f = 1 - d - b.
#' 
#' @param daf data frame containing DAF, Pi and P0 values
#' @param divergence data frame containing divergent and analyzed sites for selected (i) and neutral (0) classes
#' @param xlow lower limit for asymptotic alpha fit
#' @param xhigh higher limit for asymptotic alpha fit
#' @param seed seed value (optional). No seed by default
#' @param plot report plots of daf, alpha and negative selection fractions (optional). Default is FALSE
#'
#' @return iMKT method. List with asymptotic MK table and values, fractions of sites and graphs of DAF, asymptotic alpha model and negative selection fractions (optional).
#'
#' @examples
#' ## Without plot
#' iMK(myDafData, myDivergenceData, xlow=0, xhigh=0.9)
#' ## With plot
#' iMK(myDafData, myDivergenceData, xlow=0, xhigh=0.9, plot=TRUE)
#' 
#' @import utils
#' @import stats
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom ggthemes theme_foundation
#' @importFrom cowplot plot_grid
#' 
#' @keywords MKT
#' @export

iMK <- function(daf, divergence, xlow, xhigh, seed, plot=FALSE) {
  
  ## Check data
  check <- checkInput(daf, divergence, xlow, xhigh)
  if (check$data == FALSE) {
    stop (check$print_errors) }
    
  ## Check seed
  if(missing(seed)) {
    seed <- NULL
  } else {
    set.seed(seed)
  }
  
  ## Create MKT table standard
  mkt_table_standard <- data.frame(Polymorphism = c(sum(daf$P0), sum(daf$Pi)), 
                                   Divergence=c(divergence$D0,divergence$Di),
                                   row.names = c("Neutral class","Selected class"))
  
  ## Total number of sites analyzed 
  mi <- divergence$mi
  m0 <- divergence$m0
  
  ## Run asymptotic MK and retrieve alphas 
  asymptoticMK_table <- asymptoticMK(daf, divergence, xlow, xhigh)
  alpha_asymptotic <- as.numeric(asymptoticMK_table$alpha_asymptotic)
  alpha_standard <- as.numeric(asymptoticMK_table$alpha_original)
  alpha_CI_low <- asymptoticMK_table$CI_low 
  
  ## Estimate the relative proportion of non-synonymous and synonymous substitutions
  daf$N <- daf$Pi/sum(daf$Pi)   
  daf$S <- daf$P0/sum(daf$P0)
   
  ## Estimate alpha for each DAF category
  daf$alpha <- 1-((mkt_table_standard[1,2]*daf$Pi)/(mkt_table_standard[2,2]*daf$P0))
  
  ## Estimate the synonymous and non-synonymous ratio
  synonymous_ratio <- mkt_table_standard[1,1]/m0    
  nonsynonymous_ratio <- mkt_table_standard[2,1]/mi  
  
  ## Estimate the fraction of neutral sites (f)
  f <- nonsynonymous_ratio/synonymous_ratio
  
  ## Estimate the fraction of strongly deleleterious sites (d)
  d <- 1-f   
  
  ## Estimate the fraction of weakly deleterious sites (b)
  daf <- na.omit(daf); daf <- droplevels(daf)
  wd <- 0
  for (i in 1:nrow(daf)) {
    row <- daf[i,]
    if (row$alpha < alpha_CI_low) {
      wd <- wd + ((alpha_asymptotic-row$alpha)*row$N)
    } else { break }
  }
  wd <- wd/(alpha_asymptotic-min(daf$alpha,na.rm=T))
  b <- wd*f
  
  ## Re-estimate the truly number of neutral sites, removing the slightly deleterious 
  f <- f-b
  
  ## Fraction of f, b and d sites
  fraction <- data.frame(Fraction=c(d,f,b), Type=c("d","f","b"), MKT=rep("Asymptotic MK", 3))
  
  ## Perform plots 
  if(plot == TRUE) {
    
    ## plot fractions

   plotfraction <- ggplot(fraction) + geom_bar(stat="identity", aes_string(x="MKT", y="Fraction", fill="Type"), color="black") +
      coord_flip() + themePublication() + ylab(label="Fraction") + xlab(label="Cut-off") +
      scale_fill_manual(values=c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33"), breaks=c("f","d","b"), labels=c(expression(italic("f")),expression(italic("d")),expression(italic("b")))) +
      theme(axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.line.x =element_blank())  + scale_y_discrete(limit=seq(0,1,0.25), expand=c(0,0))
    
    ## DAF graph
    daf_graph <- daf[c("daf","Pi","P0")]
    daf_graph<-melt(daf_graph,id.vars = "daf")
  
    plotdaf <-  ggplot(daf_graph) +
      geom_point(aes_string(x="daf", y="value", color="variable"), size=3) +
      themePublication() + scale_color_manual(values=c("#386cb0","#fdb462"), name="Type", breaks=c("Pi","P0"), labels=c("Non-synonymous","Synonymous")) +
      xlab("Derived Allele Frequency") + ylab("Number of Sites")
  
    ## Alpha graph
    y1 <- function(daf) {
      asymptoticMK_table$a+asymptoticMK_table$b*exp(-asymptoticMK_table$c*daf) }
    xs <- seq(xlow, xhigh, length.out=nrow(daf)*5)
    ysmax <- rep(asymptoticMK_table$alpha_asymptotic, length(xs))
    ysmin <- y1(xs)
    shader_df <- data.frame(xs, ysmin, ysmax)
 
    plot_alpha <- ggplot(daf, aes_string(x="daf", y="alpha"))  +
      ## Confidence intervals
      geom_rect(data=data.frame(xmin=-Inf, xmax=Inf, ymin=asymptoticMK_table$CI_low, ymax=asymptoticMK_table$CI_high), aes_string(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax"), fill="gray30", alpha=0.4, inherit.aes=F) +
      ## Function curve
      stat_function(fun=y1, color="#ef3b2c", size=2) +
      ## Asymptotic alpha
      geom_hline(yintercept=asymptoticMK_table$alpha_asymptotic, color="#662506", linetype="dashed") +  
      ## Alpha derived via classic MKT
      geom_hline(yintercept=asymptoticMK_table$alpha_original, color="#386cb0", linetype="dashed") +
      ## Cut-offs
      geom_vline(xintercept=c(xlow, xhigh), color="gray10", linetype="dotted") +    
      ## Points
      geom_point(size=3, color="gray15") +
      ## Shade the fraction of WDMs
      geom_ribbon(data=shader_df, aes_string(x="xs", ymin="ysmin", ymax="ysmax"), fill="gray30", alpha=0.2, inherit.aes=F) + 
      ## Customization
      themePublication() + 
      xlab("Derived allele frequency") + ylab(expression(bold(paste("Adaptation (",alpha,")")))) +
      ## Alphas labels
      annotate("text", x=xhigh-0.2, y=asymptoticMK_table$alpha_asymptotic-0.2, label=paste0('alpha [asymptotic] == ', round(asymptoticMK_table$alpha_asymptotic, digits = 3)), parse=T, color="#662506", size=4) +
      annotate("text",x=xhigh-0.2, y=asymptoticMK_table$alpha_original-0.1, label=paste0('alpha [standard] == ', round(asymptoticMK_table$alpha_original,digits = 3)), parse=T, color="#386cb0", size=4) 
  
    ## Render plots with labels
    plots_iMKT <- plot_grid(plotdaf, plot_alpha, plotfraction, nrow=3,  labels=c("A","B","C"), rel_heights=c(2,2,1))
    
    ## Store output
    ## iMKT output
    fraction <- fraction[c("Type","Fraction")]
    asymptoticMK_table[2:8] <- round(asymptoticMK_table[2:8],4)
    output <- list(asymptoticMK_table, fraction, plots_iMKT)
    names(output) <- c("Asymptotic MK table", "Fractions of sites", "Graphs")
    
  ## If no plot to perform  
  } else if (plot == FALSE) {
    ## iMKT output
    fraction <- fraction[c("Type","Fraction")]
    asymptoticMK_table[2:8] <- round(asymptoticMK_table[2:8],4)
    output <- list(asymptoticMK_table, fraction)
    names(output) <- c("Asymptotic MK table", "Fractions of sites")
  }
  
  ## Return output
  return(output) 
}
