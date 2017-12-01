#' Daf plot 
#' Date = 30/11/2016
#' Author = Marta Coronado, Sergi Hervàs
#'
#'
#' @param x daf file: daf, pN, pS
#'
#' @return None
#'
#' @examples
#' @export
#' 
## input 

plotDAF <- function(x) {
  
  #load package
  require(ggplot2, quietly=T)
  
  ###Column names not associated
  ###names(x)[1] <- "daf"
  x$pN <- as.numeric(x$pN)
  x$N <- x$pN/sum(x$pN)
  x$pS <- as.numeric(x$pS)  
  x$S <- x$pS/sum(x$pS)

  p1 <- ggplot() +
    geom_point(data=x, aes(x=daf, y=pS, color="p0"), size=3) + 
    geom_point(data=x, aes(x=daf, y=pN, color="pi"), size=3) +
    theme_classic() +
    theme(axis.text=element_text(size=rel(1.5))) +
    theme(axis.title=element_text(size=rel(1.5))) +
    xlab("Derived Allele Frequency") + ylab("Number of Sites") +
    ggtitle("DAF distribution\nNumber of Sites") + 
    theme(plot.title=element_text(size=rel(1.5))) +
    scale_colour_manual(name="", values=c(p0="blue", pi="red")) +
    theme(legend.text = element_text(size=rel(1.5)))
                        
  p2 <- ggplot() +
    geom_point(data=x, aes(x=daf, y=S, color="p0"), size=3) + 
    geom_point(data=x, aes(x=daf, y=N, color="pi"), size=3) +
    theme_classic() +
    theme(axis.text=element_text(size=rel(1.5))) +
    theme(axis.title=element_text(size=rel(1.5))) +
    xlab("Derived Allele Frequency") + ylab("Proportion of Sites") +
    ggtitle("DAF distribution\nProportion of sites") + 
    theme(plot.title=element_text(size=rel(1.5))) +
    scale_colour_manual(name="", values=c(p0="blue", pi="red")) +
    theme(legend.text = element_text(size=rel(1.5)))
  
  p <- list(p1, p2)
  names(p) <- c("NumberSites", "ProportionSites")
  return(p)
}
