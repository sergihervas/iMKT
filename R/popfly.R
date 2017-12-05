#' Extract a list of genes from PopFly
#' Date = 30/11/2016
#' Author = Sergi Hervas, Jesús Murga, Marta Coronado
#'
#'
#' @param genes Drosophila gene list
#' @param pops Drosophila populations
#' @param recomb group genes according to recombination values (must specify number of bins)
#' @param bin number of recombination bins to compute (mandatory if recomb = TRUE)
#' @return None
#' @examples
#' @import utils
#' @export
#' 

#######################
### IN CONSTRUCTION ###
#######################

### DO NOT REMOVE: ###
  
  ## data SHOULD BE pre-loaded in: genesPopFly
  # genesPopFly <- read.table("/home/sergi/Desktop/GenesData_recomb.tab", header=T)
  ## test data 
  # subsetPopFly("FBgn0000008",c("AM","AUS","SA"), F, 3)

### END OF DO NOT REMOVE ###

subsetPopFly <- function(genes=c("gene1","gene2","..."), pops=c("pop1","pop2","..."), recomb=TRUE/FALSE, bins){ 
  
  ## CHECK INPUT VARIABLES, ERROR HANDLING ##
  ## numer of arguments
  if (nargs() != 3 && nargs() != 4) {
    stop("You must specify 3 arguments at least: genes, pops, and recomb (T/F)") }
  
  ## argument genes
  if (length(genes) == 0 || genes == "" || !is.character(genes)) {
    stop("You must specify at least one gene.") }
  if (!all(genes %in% genesPopFly$Name) == TRUE) { ## --> report which genes cause problem <-- ##
    stop("MKT data is not available for the requested gene(s).\nRemember to use FlyBase IDs (FBgn...).") }
  
  ## argument pops  
  if (length(pops) == 0 || pops == "" || !is.character(pops)) { 
    stop("You must specify at least one population.") }
  if (!all(pops %in% genesPopFly$Pop) == TRUE) {
    stop("Select at least one of the following populations:\nAM, AUS, CHB, EA, EF, EG, ENA, EQA, FR, RAL, SA, SD, SP, USI, USW, ZI") }
  
  ## argument recomb
  if (recomb != TRUE && recomb != FALSE) {
    stop("Parameter recomb must be TRUE or FALSE.")
  }
  
  ## argument bins
  if (recomb == TRUE) {
    if (nargs() != 4 || !is.numeric(bins) || bins == 0  || bins == 1) {
      stop("If recomb = TRUE, you must specify the number of bins to use (>1).") }
  }
  if (recomb == FALSE && !is.null(bins)) {
    cat("[Warning]: parameter bins not used! (recomb=F selected)\n\n")
  }
  
  ## PERFORM SUBSET ##
  subsetGenes <- genesPopFly[(genesPopFly$Name %in% genes & genesPopFly$Pop %in% pops), ]
  
  ## IF RECOMB == TRUE ##
  if (recomb == TRUE) {
    ## create bins
  }
  
  ## IF RECOMB == FALSE ##
  else if (recomb == FALSE) {
    subsetGenes$Name <- as.factor(subsetGenes$Name)
    subsetGenes <- droplevels(subsetGenes) 
    
    #Gen1: pop1, pop2 || Pop1: gen1, gen2?. Let's do Gen1:pop1, pop2
    for (i in levels(subsetGenes$Name)) {
      print(i)
      x <- subsetGenes[subsetGenes$Name == i, ]
      
      for (j in levels(x$Pop)) {
        print(j)
        x1 <- x[x$Pop == j, ]

        #daf file
        x1$DAF0f <- as.character(x1$DAF0f)
        x1$DAF4f <- as.character(x1$DAF4f)
        daf <-cbind(unlist(strsplit(x1$DAF0f, split=";")), unlist(strsplit(x1$DAF4f, split=";")))
        daf <- as.data.frame(daf)
        daf$daf <- seq(0.05,0.95,0.1)
        daf <- daf[,c(3,1,2)]
        names(daf) <- c("daf","pi","p0")
        
        #div file
        div <- cbind(x1$mi, x1$di, x1$m0, x1$d0) # check order
        div <- as.data.frame(div)
        names(div) <- c("mi","di","m0","d0")
        
        ## PERFORM iMK ##
        # iMK(daf, div, 0, 0.9)
        # think how to report!
      }
    }
  }
}

