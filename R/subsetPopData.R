#' @title subsetPopData
#'
#' @description Perform iMK using a subset of PopFly or PopHuman data defined by custom genes and populations lists 
#'
#' @details put details here
#'
#' @param data input PopFly or PopHuman data
#' @param genes list of genes
#' @param pops populations to use
#' @param recomb group genes according to recombination values (must specify number of bins)
#' @param bins number of recombination bins to compute (mandatory if recomb = TRUE)
#'
#' @return None
#'
#' @examples
#' ## Load PopFly data
#' # loadPopFly()
#' ## Perform analysis
#' # subsetPopData("PopFly", c("FBgn0053196", "FBgn0000008"),c("RAL","ZI"), F)
#'
#' @import utils
#'
#' @export
#' 

#######################
### IN CONSTRUCTION ###
#######################

## test data 
## gen <- c("FBgn0000008","FBgn0000014","FBgn0000015","FBgn0000017","FBgn0000018","FBgn0000022","FBgn0000024","FBgn0000028","FBgn0000032","FBgn0000036","FBgn0000037","FBgn0000038","FBgn0000039","FBgn0000042","FBgn0000043","FBgn0000044","FBgn0000045","FBgn0000046","FBgn0000047","FBgn0000052","FBgn0000053","FBgn0000054","FBgn0000055","FBgn0000056","FBgn0000057","FBgn0000061","FBgn0000063","FBgn0000064","FBgn0000071","FBgn0000075","FBgn0000077","FBgn0000078","FBgn0000079","FBgn0000083","FBgn0000084","FBgn0000092","FBgn0000094","FBgn0000097","FBgn0000099","FBgn0000100","FBgn0000108","FBgn0000109","FBgn0000114","FBgn0000115","FBgn0000116","FBgn0000117","FBgn0000119","FBgn0000120","FBgn0000121","FBgn0000137","FBgn0000139","FBgn0000140","FBgn0000146","FBgn0000147","FBgn0000150","FBgn0000152","FBgn0000153","FBgn0000157","FBgn0000158","FBgn0000163","FBgn0000166","FBgn0000173","FBgn0000179","FBgn0000180","FBgn0000181","FBgn0000182","FBgn0000183","FBgn0000206","FBgn0000210","FBgn0000212","FBgn0000216","FBgn0000221","FBgn0000227","FBgn0000228","FBgn0000229","FBgn0000233","FBgn0000239","FBgn0000241","FBgn0000244","FBgn0000246","FBgn0000247","FBgn0000250","FBgn0000251","FBgn0000253","FBgn0000256","FBgn0000257","FBgn0000259","FBgn0000261","FBgn0000273","FBgn0000274","FBgn0000276","FBgn0000277","FBgn0000278","FBgn0000279","FBgn0000283","FBgn0000286","FBgn0000287","FBgn0000289","FBgn0000299","FBgn0000303")
## subsetPopData("PopFly", gen ,c("RAL","ZI"), T, 3)
## To Do:
  ## if predictNLS fails, error handling
  ## think how to report / return --> list by pop?

subsetPopData <- function(data=c("PopFly","PopHuman"), genes=c("gene1","gene2","..."), pops=c("pop1","pop2","..."), recomb=TRUE/FALSE, bins=0){ 
  
  ## GET POPFLY or POPHUMAN PRELOADED DATA ##
  if (data == "PopFly"){
    if (exists("PopFlyData") == TRUE) {
      data <- get("PopFlyData")
    } else {
      stop("Load PopFly data, with the command loadPopFly(), before using this function.") }}
  
  else if (data == "PopHuman"){
    if (exists("PopHumanData") == TRUE) {
      data <- get("PopHumanData")
    } else {
      stop("Load PopHuman data, with the command loadPopHuman(), before using this function.") }}
  
  
  ## CHECK INPUT VARIABLES, ERROR HANDLING ##
  ## numer of arguments
  if (nargs() != 4 && nargs() != 5) {
    stop("You must specify 4 arguments at least: data, genes, pops, and recomb (T/F)") }

  ## argument genes
  if (length(genes) == 0 || genes == "" || !is.character(genes)) {
    stop("You must specify at least one gene.") }
  if (!all(genes %in% data$Name) == TRUE) {
    difGenes <- setdiff(genes, data$Name)
    difGenes <- paste(difGenes, collapse=", ")
    stopMssg <- paste0("MKT data is not available for the requested gene(s).\nRemember to use FlyBase IDs (FBgn...)\nThe genes that caused the error are: ", difGenes, ".")
    stop(stopMssg) }

  ## argument pops
  if (length(pops) == 0 || pops == "" || !is.character(pops)) {
    stop("You must specify at least one population.") }
  if (!all(pops %in% data$Pop) == TRUE) {
    correctPops <- c("AM","AUS","CHB","EA","EF","EG","ENA","EQA","FR","RAL","SA","SD","SP","USI","USW","ZI")
    difPops <- setdiff(pops, correctPops)
    difPops <- paste(difPops, collapse=", ")
    stopMssg <- paste0("MKT data is not available for the sequested populations(s).\nSelect among the following populations:\nAM, AUS, CHB, EA, EF, EG, ENA, EQA, FR, RAL, SA, SD, SP, USI, USW, ZI.\nThe populations that caused the error are: ", difPops, ".")
    stop(stopMssg) }

  ## argument recomb
  if (recomb != TRUE && recomb != FALSE) {
    stop("Parameter recomb must be TRUE or FALSE.")
  }

  ## argument bins
  if (recomb == TRUE) {
    if (nargs() != 5 || !is.numeric(bins) || bins == 0  || bins == 1) {
      stop("If recomb = TRUE, you must specify the number of bins to use (> 1).") }
    if (bins > round(length(genes)/2)) {
      stop("Parameter bins > (genes/2). At least 2 genes for each bin are required.") }
  }
  if (recomb == FALSE && bins != 0) {
    warning("Parameter bins not used! (recomb=F selected)")
  }

  ## PERFORM SUBSET ##
  subsetGenes <- data[(data$Name %in% genes & data$Pop %in% pops), ]
  subsetGenes$Name <- as.factor(subsetGenes$Name)
  subsetGenes <- droplevels(subsetGenes)

  ## IF RECOMB == TRUE ##
  if (recomb == TRUE) {
    
    for (k in levels(subsetGenes$Pop)) {
      print(k)
      x <- subsetGenes[subsetGenes$Pop == k, ]
      x <- x[order(x$p_bp), ]
      
      ## create bins
      binsize <- round(nrow(x)/bins) ## number of genes for each bin
      count <- 1; dat <- NULL

      for (i in 0:nrow(x)) {
        if (i%%binsize == 0) { ## only if reminder of division = 0 (equally sized bins)
          i1 <- i + binsize - 1
          if (i == 0) {
            g1 <- x[i:binsize,]
            group <- count
            g1$Group <- group
            dat <- rbind(dat,g1)
            count <- count+1
          } else if (i1 < nrow(x)) {
            g1 <- x[i:i1,]
            group <- count
            g1$Group <- group
            dat <- rbind(dat,g1)
            count <- count+1
          }
        }
      }
      dat$Group <- as.factor(dat$Group)
      dat <- unique(dat)
      
      ## iterate through each recomb bin
      for (j in levels(dat$Group)) {
        print(j)
        x1 <- dat[dat$Group == j, ]
        ## set counters to 0
        Pi <- c(0,0,0,0,0,0,0,0,0,0)
        P0 <- c(0,0,0,0,0,0,0,0,0,0)
        f <- seq(0.05,0.95,0.1)
        mi <- 0; m0 <- 0
        Di <- 0; D0 <- 0

        ## group genes
        x1 <- droplevels(x1)
        for (l in levels(x1$Name)) {
          x2 <- x1[x1$Name == l, ]

          ## DAF
          x2$DAF0f <- as.character(x2$DAF0f); x2$DAF4f <- as.character(x2$DAF4f)
          daf0f <- unlist(strsplit(x2$DAF0f, split=";"))
          daf4f <- unlist(strsplit(x2$DAF4f, split=";"))
          daf0f <- as.numeric(daf0f); daf4f <- as.numeric(daf4f)
          Pi <- Pi + daf0f; P0 <- P0 + daf4f

          ## Divergence
          mi <- mi + x2$mi; m0 <- m0 + x2$m0
          Di <- Di + x2$di; D0 <- D0 + x2$d0
       }

       ## Proper formats
       daf <- cbind(f, Pi, P0); daf <- as.data.frame(daf)
       names(daf) <- c("daf","Pi","P0")
       div <- cbind(mi, Di, m0, D0); div <- as.data.frame(div)
       names(div) <- c("mi","Di","m0","D0")

       ## PERFORM iMK ##
       out <- iMK(daf, div, 0, 0.9)
      }
    }
    
    ## Warning if some genes are lost. Bins must be equally sized.
    if (nrow(dat) != length(genes)) {
      missingGenes <- round(length(genes) - nrow(dat))
      genesNames <- as.vector(tail(x, missingGenes)$Name)
      genesNames <- paste(genesNames, collapse=", ")
      warningMssg <- paste0("The ",missingGenes," gene(s) with highest recombination rate estimates (", genesNames, ") was/were excluded from the analysis in order to get equally sized bins.\n")
      warning(warningMssg)
    }
  }

  ## IF RECOMB == FALSE ##
  else if (recomb == FALSE) {

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
        x1 <- x[x$Name == j, ]

        ## DAF
        x1$DAF0f <- as.character(x1$DAF0f); x1$DAF4f <- as.character(x1$DAF4f)
        daf0f <- unlist(strsplit(x1$DAF0f, split=";"))
        daf4f <- unlist(strsplit(x1$DAF4f, split=";"))
        daf0f <- as.numeric(daf0f); daf4f <- as.numeric(daf4f)
        Pi <- Pi + daf0f; P0 <- P0 + daf4f

        ## Divergence
        mi <- mi + x1$mi; m0 <- m0 + x1$m0
        Di <- Di + x1$di; D0 <- D0 + x1$d0
      }
      
      ## Proper formats
      daf <- cbind(f, Pi, P0); daf <- as.data.frame(daf)
      names(daf) <- c("daf","Pi","P0")
      div <- cbind(mi, Di, m0, D0); div <- as.data.frame(div)
      names(div) <- c("mi","Di","m0","D0")
      
      ## PERFORM iMK ##
      out <- iMK(daf, div, 0, 0.9)
    }
  }
}

