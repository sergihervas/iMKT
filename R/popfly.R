#' Extract a list of genes from PopFly
#' Date = 30/11/2016
#' Author = Jesús Murga, Marta Coronado
#'
#'
#' @param genes Drosophila gene list
#' @param pops Drosophila populations
#' @param recomb group genes according to recombination values (must specify bin)
#' @param bin number of recombination bins to compute
#' @return None
#' @examples
#' @import utils
#' @export
#' 

################################NEED A FOLDER WITH FILES TO CHECK HOW THE OUTPUT IS RETURNED################################

#data should be pre-loaded, so the function takes as input (genes, pops, recomb, bins)
subsetPopFly <- function(genes=c("gene1","gene2","..."), pops=c("pop1","pop2","..."), recomb=TRUE/FALSE, bins) { 
  
  if (length(genes) == 0 || length(pops) == 0) #check this...
      stop("You must specify at least one gene and one population.")
  
  #if data is loaded in a sth called "popflydata"
  subsetGenes <- popflydata[ ,(genes %in% popflydata$Name & POP &in% popflydata$Pop)] 
  
  if (recomb == FALSE) {
    #do not keep recomb cols
    subsetGenes <- subsetGenes[,-17] #dnt remember col name
 
    #subsetGene contains the desired info. now, transform format and execute iMK
    subsetGenes <- droplevels(subsetGenes)
    for (i in levels(subsetGenes$Pop)) {
      x <- subsetGenes[subsetGenes$Pop == i, ]

      #x file
      x$DAF0f <- as.character(x$DAF0f)
      daf0f <- unlist(strsplit(x$DAF0f, split=";"))
      x$DAF4f <- as.character(x$DAF4f)
      daf4f <- unlist(strsplit(x$DAF4f, split=";"))
      x1 <- cbind(daf0f,daf4f)
      x1 <- as.data.frame(x1)
      x1$daf <- seq(0.05,0.95,0.1)

      #y file (must change input to add m values)
      y <- rbind(subsetGenes$m0f, subsetGenes$D0f, subsetGenes$m4f, subsetGenes$D4f) #check order
      y <- as.data.frame(y)
      
      #perform iMK and think how to return...
      iMK(x, y, 0, 0.9)
    }
  }
  
  else if (recomb == TRUE) { #if recomb=T user must specify bins!
    if (bins == 0 || is.na(bins) || !is.integer(bins))
      stop("If recomb=T, you must specify a number of bins (between 1 and 100)") #100? xd
    
    #first bins loop. assign bin id to each gene 
    
    #group bins and sum pi, p0, di, d0, mi, m0
    
    #perform iMK for each bin and think how to return...
    iMK(x, y, 0, 0.9) 
  }
}
  
