
## iMK horizontal plot 
## input: table from iMK (1 observation) = iMK$iMK

plotIMK <- function(df) {
  
  #load packages
  require(ggplot2, quietly=TRUE)
  
  names(df) <- c("param","value"); rownames(df) <- NULL
  
  df <- df[df$param != "alpha_original", ]
  df <- df[df$param != "alpha_asymptotic", ]
  df <- df[df$param != "alpha_DGRP05", ]
  df <- df[df$param != "alpha_DGRP10", ]
  df <- droplevels(df)
  
  #replace names
  df$param <- as.character(df$param)
  df$param[df$param == "d"] <- "d: strongly deleterious"
  df$param[df$param == "b"] <- "b: weakly deleterious"
  df$param[df$param == "f"] <- "f: neutral sites"
  df$param <- as.factor(df$param)
  df$param <- factor(df$param, levels=c("d: strongly deleterious","b: weakly deleterious","f: neutral sites"))
  df$value <- as.numeric(as.character(df$value))
  
  #label position
  pos <- cumsum(df$value) - (df$value * 0.5) 
  df$pos <- pos
  
  p <- ggplot(df, aes(x="", y=value, fill=param)) + geom_bar(stat="identity", color="black") +
    geom_text(aes(label=round(value,3), y=pos)) +
    scale_fill_manual(values=c("darkgreen","limegreen","greenyellow")) +
    theme_classic() + xlab("") + ylab("") +
    theme(axis.ticks.y=element_blank()) + theme(axis.text.y=element_blank()) +
    theme(legend.title=element_blank()) + theme(legend.position="bottom") +
    theme(legend.text=element_text(size=rel(1.5))) + 
    theme(axis.text.x=element_text(size=rel(1.5))) +  
    theme(line = element_blank()) +
    #ggtitle("Integrative MKT") + theme(plot.title=element_text(size=rel(2))) +
    coord_flip()
    
  return(p)
}
