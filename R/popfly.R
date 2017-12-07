#' @title subsetPopfly
#'
#' @description Extract a list of genes from PopFly
#'
#' @details put details here
#'
#' @param data input PopFly or PopHuman data
#' @param genes Drosophila gene list
#' @param pops Drosophila populations
#' @param recomb group genes according to recombination values (must specify number of bins)
#' @param bins number of recombination bins to compute (mandatory if recomb = TRUE)
#'
#' @return None
#'
#' @examples
#'
#' @import utils
#'
#' @export
#' 

#######################
### IN CONSTRUCTION ###
#######################


### DO NOT REMOVE: ###
  ## test data 
  ## subsetPopFly("popfly", c("FBgn0053196", "FBgn0000008"),c("RAL","ZI"), F, 3)
  ## To DO:
    ## if predictNLS fails, error handling
    ## if no data for that population / gene, notify
    ## think how to report / return --> list by pop?
    ## recomb = T
### END OF DO NOT REMOVE ###

subsetPopFly <- function(data=c("popfly","pophuman"), genes=c("gene1","gene2","..."), pops=c("pop1","pop2","..."), recomb=TRUE/FALSE, bins){ 
  
  if (data == "popfly"){
    if (exists("PopFlyData") == TRUE) {
      data <- PopFlyData
    } else {
      stop("Load PopFly data, with the command loadPopFly(), before using this function.")
    }
  }

  ## CHECK INPUT VARIABLES, ERROR HANDLING ##
  ## numer of arguments
  if (nargs() != 4 && nargs() != 5) {
    stop("You must specify 3 arguments at least: genes, pops, and recomb (T/F)") }

  ## argument genes
  if (length(genes) == 0 || genes == "" || !is.character(genes)) {
    stop("You must specify at least one gene.") }
  if (!all(genes %in% data$Name) == TRUE) { ## --> report which genes cause problem <-- ##
    stop("MKT data is not available for the requested gene(s).\nRemember to use FlyBase IDs (FBgn...).") }

  ## argument pops
  if (length(pops) == 0 || pops == "" || !is.character(pops)) {
    stop("You must specify at least one population.") }
  if (!all(pops %in% data$Pop) == TRUE) {
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
    warning("Parameter bins not used! (recomb=F selected)")
  }

  ## PERFORM SUBSET ##
  subsetGenes <- data[(data$Name %in% genes & data$Pop %in% pops), ]

  ## IF RECOMB == TRUE ##
  if (recomb == TRUE) {
    ## create bins
  }

  ## IF RECOMB == FALSE ##
  else if (recomb == FALSE) {
    subsetGenes$Name <- as.factor(subsetGenes$Name)
    subsetGenes <- droplevels(subsetGenes)

    for (i in levels(subsetGenes$Pop)) {
      print(i)
      x <- subsetGenes[subsetGenes$Pop == i, ]
      
      ## set counters to 0
      Pi <- c(0,0,0,0,0,0,0,0,0,0)
      P0 <- c(0,0,0,0,0,0,0,0,0,0)
      f <- seq(0.05,0.95,0.1)
      mi <- 0; m0 <- 0
      Di <- 0; D0 <- 0
      
      ## group genes
      for (j in levels(x$Name)) {
        print(j)
        x1 <- x[x$Name == j, ]

        #daf file
        x1$DAF0f <- as.character(x1$DAF0f)
        x1$DAF4f <- as.character(x1$DAF4f)
        daf0f <- unlist(strsplit(x1$DAF0f, split=";"))
        daf4f <- unlist(strsplit(x1$DAF4f, split=";"))
        daf0f <- as.numeric(daf0f)
        daf4f <- as.numeric(daf4f)
        Pi <- Pi + daf0f
        P0 <- P0 + daf4f

        #div file
        mi <- mi + x1$mi
        m0 <- m0 + x1$m0
        Di <- Di + x1$di
        D0 <- D0 + x1$d0
      }
      
      daf <- cbind(f, Pi, P0)
      daf <- as.data.frame(daf)
      names(daf) <- c("daf","Pi","P0")
      div <- cbind(mi, Di, m0, D0)
      div <- as.data.frame(div)
      names(div) <- c("mi","Di","m0","D0")
      
      
      ## PERFORM iMK ##
      out <- iMK(daf, div, 0, 0.9)
      print(out)[[1]]
      
    }
  }
}

