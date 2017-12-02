#' iMK compute all the MKT extensions 
#' Date = 30/11/2016
#' Author = Sergi Hervás, Marta Coronado
#'
#' EXPLANATION
#'
#'
#' @param daf dad file
#' @param divergence div file
#' @param xlow fit curve
#' @param xhigh curv
#'
#' @return None
#'
#' @examples
#' @import knitr 
#' @import utils
#' @import stats
#' @import grid 
#' @import gridExtra
#' @import scales
#' @import reshape2
#' @importFrom ggplot2
#' @importFrom ggthemes theme_foundation theme
#' @import cowplot plot_grid
#' @export
#' 

theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

iMK <- function(daf, divergence, xlow, xhigh) {
  
  check <- check_input(daf, divergence, xlow, xhigh)
  if (check$data==FALSE)
    stop (check$print_errors)
  
  # Create MKT table standard
  mkt_table_standard <- data.frame(Polymorphism = c(sum(daf$pS), sum(daf$pN)), 
                                   Divergence=c(divergence$D4f,divergence$D0f),
                                   row.names = c("Neutral class","Selected class"))
  
  # Total number of sites analyzed 
  mi <- as.numeric(divergence$m0f)
  m0 <- as.numeric(divergence$m4f)
  
  # names(daf) <- c("daf","pN","pS") Borrar? Se supone que se chequea que asi sea el header no?
  
  # Run asymptotic MK and retrieve alphas 
  asymptoticMK_table <- asymptoticMK(daf, divergence, xlow, xhigh)
  alpha_asymptotic <- as.numeric(asymptoticMK_table$alpha_asymptotic)
  alpha_standard <- as.numeric(asymptoticMK_table$alpha_original)
  alpha_CI_low <- asymptoticMK_table$CI_low 
  
  # Estimate the relative proportion of non-synonymous and synonymous substitutions
  daf$N <- daf$pN/sum(daf$pN)   
  daf$S <- daf$pS/sum(daf$pS)
  daf$pN <- as.numeric(daf$pN) # Esto va aqui?
  daf$pS <- as.numeric(daf$pS)
  
  ## Estimate alpha for each DAF category
  daf$alpha <- 1-((mkt_table_standard[1,2]*daf$pN)/(mkt_table_standard[2,2]*daf$pS))
  
  ## Estimate the synonymous and non-synonymous ratio
  synonymous_ratio <- mkt_table_standard[1,1]/m0    
  nonsynonymous_ratio <- mkt_table_standard[2,1]/mi  
  
  ## Estimate the fraction of neutral sites (f)
  f <- nonsynonymous_ratio/synonymous_ratio
  
  ## Estimate the fraction of strongly deleleterious sites (d)
  d <- 1-f   
  
  daf <- na.omit(daf); daf <- droplevels(daf) # UY ??? esto para que es? es necesario?
  
  ## Estimate the fraction of weakly deleterious sites (b)
  wd <- 0
  for (i in 1:nrow(daf)) {
    row <- daf[i,]
    if (row$alpha < alpha_CI_low) {
      wd <- wd+((alpha_asymptotic-row$alpha)*row$N)
      } 
    else { 
      break
      }
    }
  
  wd <- wd/(alpha_asymptotic-min(daf$alpha,na.rm=T)) # na.rm por que?
  b <- wd*f
  
  # Re-estimate the truly number of neutral sites, removing the slightly deleterious 
  f <- f-b
  
  # Fraction of f, b and d sites
  fraction <-data.frame(Fraction=c(d,f,b), Type = c("d","f","b"), MKT = rep("Asymptotic MKT", 3))
  
  plotfraction <- ggplot(fraction) + geom_bar(stat = "identity",aes_string(x = "MKT", y = "Fraction", fill = "Type")) + coord_flip() + theme_Publication() + ylab(label = "Fraction") + scale_fill_manual(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33"), breaks=c("d","f","b"), labels=c(expression(italic("d")), expression(italic("f")),expression(italic("b")))) + theme(axis.title.y=element_blank())
  plotfraction
  
  # DAF graph
  daf_graph <- daf[c("daf","pN","pS")]
  daf_graph<-melt(daf_graph,id.vars = "daf")
  
 plotdaf <-  ggplot(daf_graph) +
    geom_point(aes_string(x="daf", y="value", color="variable"), size=3) +
    theme_Publication() + scale_color_manual(values =  c("#386cb0","#fdb462"), name ="Type", breaks=c("pN", "pS"), labels=c("Non-synonymous", "Synonymous")) +
    xlab("Derived Allele Frequency") + ylab("Number of Sites")
 plotdaf
 
 # Alpha graph
 
 y1 <- function(daf) {
   asymptoticMK_table$a+asymptoticMK_table$b*exp(-asymptoticMK_table$c*daf)
   }
 
 xs <- seq(xlow, xhigh, length.out=nrow(daf)*5)
 ysmax <- rep(asymptoticMK_table$alpha_asymptotic, length(xs))
 ysmin <- y1(xs)
 shader_df <- data.frame(xs, ysmin, ysmax)
 
 plot_alpha <- ggplot(daf, aes_string(x="daf", y="alpha"))  +
   #confidence intervals
   geom_rect(data=data.frame(xmin=-Inf, xmax=Inf, ymin=asymptoticMK_table$CI_low ,ymax=asymptoticMK_table$CI_high), aes_string(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax"), fill="gray30", alpha=0.4, inherit.aes=F) +
   #function curve
   stat_function(fun=y1, color="#ef3b2c", size = 2) +
   #asymptotic alpha
   geom_hline(yintercept=asymptoticMK_table$alpha_asymptotic, color="#662506", linetype="dashed") +  
   #alpha derived via classic MKT
   geom_hline(yintercept=asymptoticMK_table$alpha_original, color="#386cb0", linetype="dashed") +
   #cut-offs
   geom_vline(xintercept=c(xlow, xhigh), color="gray10", linetype="dotted") +    
   # points
   geom_point(size=3, color = "gray15") +
 
   #shade the fraction of WDMs
   geom_ribbon(data=shader_df, aes_string(x="xs", ymin="ysmin", ymax="ysmax"), fill="gray30", alpha=0.2,inherit.aes=F) + 
   #customization
   theme_Publication() + #theme(panel.grid = element_blank()) +
   xlab("Derived allele frequency") + ylab("Alpha") +
   #alphas labels
   annotate("text", x=xhigh-0.2, y=asymptoticMK_table$alpha_asymptotic-0.2, label=paste0("alpha asymptotic = ", round(asymptoticMK_table$alpha_asymptotic, digits = 3)), color="#662506", size=4) +
   annotate("text",x=xhigh-0.2, y=asymptoticMK_table$alpha_original-0.1, label=paste0("alpha original = ",round(asymptoticMK_table$alpha_original, digits = 3)), color="#386cb0", size=4)
  plot_alpha
  
  plots_iMKT <- plot_grid(plot_alpha,plotdaf, plotfraction, nrow = 3,  labels = c("A", "B", "C"), rel_heights = c(2, 2, 1))
  
  ## iMKT output
 fraction <- fraction[c("Type","Fraction")]
  
  asymptoticMK_table
  asymptoticMK_table[2:8] <- round(asymptoticMK_table[2:8],4)
  
  output <- list(asymptoticMK_table,fraction, plots_iMKT)
  names(output) <- c("Asymptotic MK table", "Fractions of sites", "Graphs")
  
  return(output) 
}
