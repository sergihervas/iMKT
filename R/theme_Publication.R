#' MKT corrected with DGRP method
#'
#' Date = 30/11/2016
#' Author = Sergi Hervás, Marta Coronado
#'
#' Default theme use for plot images. It has enough quality to publish anything
#'
#'
#' @param base_size base size required from theme_Publication
#' @param base_family font to load in theme_Publication
#' 
#' @return MKT corrected by the DGRP method
#'
#' @examples
#' @import cowplot
#' @import grid 
#' @import gridExtra
#' @import ggplot2 
#' @import ggthemes
#' @importFrom ggplot2 ggsave
#' @importFrom ggthemes theme_map
#' @export
theme_Publication <- function(base_size=14, base_family="helvetica") {
  # library('grid')
  # library('ggthemes')
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
