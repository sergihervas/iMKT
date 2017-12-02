#' MKT corrected with DGRP method
#'
#' Date = 30/11/2016
#' Author = Sergi Hervás, Marta Coronado
#'
#' The standard McDonald and Kreitman test (MKT) is used to detect the signature of selection at the molecular level. The MKT compares the amount of variation within a species (polymorphism, P) to the divergence (D) between species at two types of sites, one of which is putatively netral and used as the reference to detect selection at the other type of site. In the standard MKT, these sites are synonymous (putatively neutral, 0) and non-synonymous sites (selected sites, i) in a coding region. Under strict neutrality, the ratio of the number of selected and neutral polymorphic sites (Pi/P0) is equal to the ratio of the number of selected and neutral divergence sites (Di/D0).
# The null hypothesis of neutrality is rejected in a MKT when Di/D0 > Pi/P0. The excess of divergence relative to polymorphism for class i, is interpreted as adaptive selection for a subset of sites i. The fraction of adaptive fixations, α, is estimated from 1-(Pi/P0)(D0/Di). The significance of the test can be assesed with a Fisher exact test.
#' The estimate of α can be easily biased by the segregation of slightly deleterious non-synonymous substitutions. Specifically, slightly deleterious mutations tend to contribute more to polymorphism than to divergence, and thus, lead to an underestimation of alpha. Because adaptive mutations and weakly deleterious selection sct in opposite directions on the MKT, alpha and the fraction of substitutions that are sligholty deleterious, b, will be both underestimated when the two selection regimes occur. To take adaptive and slighlty deleterious mutations mutually into account, Pi, the count off segregatning sites in class i, should be seaprated into the number of neutral variants and the number of weakly deleterious variants, Pi = Pineutral + Pi weak del. Alpha is then estimated as 1-(Pineutral/P0)(D0/Di)
#'
#'
#' @param x dad file
#' @param y divergence file
#' 
#' @return None
#'
#' @examples
#' mkt_DGRP(x,y)
#'
#'MKT corrected by the DGRP method
#'$Results
#' Cut-off          α Fishers exact test P-value
#'1    0.00 0.07371226                2.175330e-14
#'2    0.05 0.25019053               1.122252e-169
#'3    0.20 0.20888331               4.789022e-115
#'
#'$Graph
#'
#'$`MKT tables`
#'$`MKT tables`$`Number of segregating sites by MAF category - Cutoff =  0`
#'
#'
#'|               | DAF.below.cutoff| DAF.above.cutoff|
#'|:--------------|----------------:|----------------:|
#'|Neutral class  |                0|            32308|
#'|Selected class |                0|            31125|
#'
#'$`MKT tables`$`Number of segregating sites by MAF category - Cutoff =  0.05`
#'
#'
#'|               | DAF.below.cutoff| DAF.above.cutoff|
#'|:--------------|----------------:|----------------:|
#'|Neutral class  |            17189|            15119|
#'|Selected class |            22490|             8635|
#'
#'$`MKT tables`$`Number of segregating sites by MAF category - Cutoff =  0.2`
#'
#'
#'|               | DAF.below.cutoff| DAF.above.cutoff|
#'|:--------------|----------------:|----------------:|
#'|Neutral class  |            21969|            10339|
#'|Selected class |            25707|             5418|
#'
#'$`MKT tables`$`MKT standard table`
#'               Polymorphism Divergence
#'Neutral class         32308      52537
#'Selected class        31125      54641
#'
#' @export
#' 



#### Libraries ####
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

library(ggplot2)
library(gridExtra)
library(cowplot)
library(reshape2)

####################################################
################# MKT-FWW function #################
####################################################

mkt_DGRP <- function(daf = "Data frame containing the DAF, Pn and Ps", 
                 divergence = "Data frame that contains the number of non-synonymous sites analyzed (m0f), number of synonymous sites analyzed (m4f), total number of non-synonymous divergent sites (D0f) and number of synonymous divergent sites (D4f)") {
  
  # Shows a message when using the function
  packageStartupMessage("MKT corrected by the DGRP method")
  
  # Declare output data frame
  output <- data.frame(cutoff = numeric(0), alpha = numeric(0), slighlty_deleterious = numeric(0), neutral = numeric(0), pvalue = integer(0))
  
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
    
    # Fisher's exact test p-value from the MKT
    m <- matrix(c(P0,Pi_neutral,divergence$D4f,divergence$D0f), ncol=2)
    pvalue <- fisher.test(m)$p.value

    # Fractions graph
    fraction <-data.frame(c(d,f,b))
    names(fraction) <- cutoff
    fractions <- cbind(fractions,fraction)
    
    # Store output  
    output_cutoff <- cbind(cutoff,alpha,pvalue)
    output <- rbind(output,output_cutoff)
    
    mkt_table <- knitr::kable(mkt_table,caption = "cutoff")
    
    mkt_tables[[paste("Number of segregating sites by MAF category - Cutoff = ",cutoff)]]  <- mkt_table
    
  } 
  mkt_tables[[paste("MKT standard table")]]  <- mkt_table_standard
  output <- as.data.frame(output)
  
  plot <- ggplot(output, aes(x=as.factor(cutoff), y=alpha, group=1)) +
    geom_line(color="#386cb0") + 
    geom_point(size=2.5, color="#386cb0")+
    theme_Publication() +
    xlab("Cut-off") + ylab("α") 
  plot
  
  names(output) <- c("Cut-off","α","Fisher's exact test P-value")
  
  fractions_melt <- melt(fractions, id.vars=NULL) 
  fractions_melt$Fraction <-  rep(c("d", "f", "b"),3)
  
  plotfraction <- ggplot(fractions_melt) + geom_bar(stat = "identity",aes(x=variable, y = value, fill = Fraction)) + coord_flip() + theme_Publication() + ylab(label = "Fraction") + xlab(label = "Cut-off") + scale_fill_manual(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33"), breaks=c("d","f","b"), labels=c(expression(italic("d")), expression(italic("f")),expression(italic("b"))))
  plotfraction
  
  plot<-plot_grid(plot,plotfraction, nrow = 2,  labels = c("A", "B"), rel_heights = c(2, 1))
  
  list_output <-list(output,plot,mkt_tables, fractions)
  names(list_output) <- c("Results","Graph", "MKT tables","Fractions")
  

  
  return(list_output)
}

#mkt_DGRP(daf,divergence)
