#' @title iMK
#'
#' @description \code{iMK()} iMK method
#'
#' @details put some details here
#'
#' @param daf daf file
#' @param divergence div file
#' @param xlow trimming values below this daf threshold
#' @param xhigh trimming values above this daf threshold
#' @param seed seed value (optional). No seed by default
#'
#' @return None
#'
#' @examples
#' ## Load your Derived Allele Frequency and Divergence files
#' daf <- mydafdata
#' div <- mydivergencedata
#' ## Run the function
#' # iMK(daf, div, 0, 0.9)
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

iMK <- function(daf, divergence, xlow, xhigh, seed) {
  
  check <- check_input(daf, divergence, xlow, xhigh)
  if (check$data == FALSE) {
    stop (check$print_errors) }
  
  
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
  mi <- as.numeric(divergence$mi)
  m0 <- as.numeric(divergence$m0)
  
  ## Run asymptotic MK and retrieve alphas 
  asymptoticMK_table <- asymptoticMK(daf, divergence, xlow, xhigh)
  alpha_asymptotic <- as.numeric(asymptoticMK_table$alpha_asymptotic)
  alpha_standard <- as.numeric(asymptoticMK_table$alpha_original)
  alpha_CI_low <- asymptoticMK_table$CI_low 
  
  ## Estimate the relative proportion of non-synonymous and synonymous substitutions
  daf$Pi <- as.numeric(daf$Pi) # Esto va aqui? En principio si el whatchdog (check data) lo comprueba ya ni hace falta
  daf$P0 <- as.numeric(daf$P0)
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
  
  daf <- na.omit(daf); daf <- droplevels(daf)
  
  ## Estimate the fraction of weakly deleterious sites (b)
  wd <- 0
  for (i in 1:nrow(daf)) {
    row <- daf[i,]
    if (row$alpha < alpha_CI_low) {
      wd <- wd+((alpha_asymptotic-row$alpha)*row$N)
    } else { break }
    }
  
  wd <- wd/(alpha_asymptotic-min(daf$alpha,na.rm=T))
  b <- wd*f
  
  ## Re-estimate the truly number of neutral sites, removing the slightly deleterious 
  f <- f-b
  
  ## Fraction of f, b and d sites
  fraction <-data.frame(Fraction=c(d,f,b), Type = c("d","f","b"), MKT = rep("Asymptotic MK", 3))
  
  plotfraction <- ggplot(fraction) + geom_bar(stat = "identity",aes_string(x = "MKT", y = "Fraction", fill = "Type")) + coord_flip() + theme_Publication() + ylab(label = "Fraction") + scale_fill_manual(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33"), breaks=c("d","f","b"), labels=c(expression(italic("d")), expression(italic("f")),expression(italic("b")))) + theme(axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(),   axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1))  +  scale_y_discrete(limit=seq(0,1,0.25), expand = c(0, 0))
  plotfraction
  
  ## DAF graph
  daf_graph <- daf[c("daf","Pi","P0")]
  daf_graph <- melt(daf_graph,id.vars = "daf")
  
  plotdaf <- ggplot(daf_graph) +
    geom_point(aes_string(x="daf", y="value", color="variable"), size=3) +
    theme_Publication() + scale_color_manual(values =  c("#386cb0","#fdb462"), name ="Type", breaks=c("Pi", "P0"), labels=c("Non-synonymous", "Synonymous")) +
    xlab("Derived Allele Frequency") + ylab("Number of Sites")
  plotdaf
 
  ## Alpha graph
  y1 <- function(daf) {
    asymptoticMK_table$a+asymptoticMK_table$b*exp(-asymptoticMK_table$c*daf)
  }
 
  xs <- seq(xlow, xhigh, length.out=nrow(daf)*5)
  ysmax <- rep(asymptoticMK_table$alpha_asymptotic, length(xs))
  ysmin <- y1(xs)
  shader_df <- data.frame(xs, ysmin, ysmax)
 
  plot_alpha <- ggplot(daf, aes_string(x="daf", y="alpha"))  +
    ## Confidence intervals
    geom_rect(data=data.frame(xmin=-Inf, xmax=Inf, ymin=asymptoticMK_table$CI_low ,ymax=asymptoticMK_table$CI_high), aes_string(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax"), fill="gray30", alpha=0.4, inherit.aes=F) +
    ## Function curve
    stat_function(fun=y1, color="#ef3b2c", size = 2) +
    ## Asymptotic alpha
    geom_hline(yintercept=asymptoticMK_table$alpha_asymptotic, color="#662506", linetype="dashed") +  
    ## Alpha derived via classic MKT
    geom_hline(yintercept=asymptoticMK_table$alpha_original, color="#386cb0", linetype="dashed") +
    ## Cut-offs
    geom_vline(xintercept=c(xlow, xhigh), color="gray10", linetype="dotted") +    
    ## Points
    geom_point(size=3, color = "gray15") +
    ## Shade the fraction of WDMs
    geom_ribbon(data=shader_df, aes_string(x="xs", ymin="ysmin", ymax="ysmax"), fill="gray30", alpha=0.2,inherit.aes=F) + 
    ## Customization
    theme_Publication() + xlab("Derived allele frequency") + ylab(expression(bold(paste("Adaptation (",alpha,")")))) +
    ## Alphas labels
    annotate("text", x=xhigh-0.2, y=asymptoticMK_table$alpha_asymptotic-0.2, label=paste0('alpha [asymptotic] == ', round(asymptoticMK_table$alpha_asymptotic, digits = 3)), parse=T, color="#662506", size=4) +
    annotate("text",x=xhigh-0.2, y=asymptoticMK_table$alpha_original-0.1, label=paste0('alpha [standard] == ', round(asymptoticMK_table$alpha_original,digits = 3)), parse=T, color="#386cb0", size=4) 
  plot_alpha
  
  ## Fraction Graph
  plots_iMKT <- plot_grid(plot_alpha,plotdaf, plotfraction, nrow = 3,  labels = c("A", "B", "C"), rel_heights = c(2, 2, 1))
  
  ## iMKT output
  fraction <- fraction[c("Type","Fraction")]
  asymptoticMK_table[2:8] <- round(asymptoticMK_table[2:8],4)
  
  output <- list(asymptoticMK_table, fraction, plots_iMKT)
  names(output) <- c("Asymptotic MK table", "Fractions of sites", "Graphs")

  invisible(output)
  # return(output)
}

