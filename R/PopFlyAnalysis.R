#' @title PopFlyAnalysis
#'
#' @description Perform any MK test using a subset of PopFly data defined by custom genes and populations lists 
#'
#' @details put details here
#'
#' @param genes list of genes
#' @param pops list of populations
#' @param recomb group genes according to recombination values (must specify number of bins). TRUE/FALSE
#' @param bins number of recombination bins to compute (mandatory if recomb = TRUE)
#' @param test which test to perform. Options include: standard (default), DGRP, FWW, asymptotic, iMK
#' @param xlow lower limit for asymptotic alpha fit (default=0)
#' @param xhigh higher limit for asymptotic alpha fit (default=1)
#' 
#' @return None
#'
#' @examples
#' ## Load PopFly data into your workspace
#' # loadPopFly()
#' ## Perform analysis
#' # mygenes <- c("FBgn0053196", "FBgn0000008")
#' # PopFlyAnalysis(mygenes , c("RAL","ZI"), recomb=F)
#' # mygenes <- c("FBgn0000008","FBgn0000014","FBgn0000015","FBgn0000017","FBgn0000018","FBgn0000022",
#' #              "FBgn0000024","FBgn0000028","FBgn0000032","FBgn0000036","FBgn0000037","FBgn0000038",
#' #              "FBgn0000039","FBgn0000042","FBgn0000043","FBgn0000044","FBgn0000045","FBgn0000046")
#' # PopFlyAnalysis(mygenes , c("RAL","ZI"), recomb=T, bins=3, test="DGRP", xlow=0, xhigh=0.9)
#'
#' @import utils
#' @import stats
#'
#' @export

PopFlyAnalysis <- function(genes=c("gene1","gene2","..."), pops=c("pop1","pop2","..."), recomb=TRUE/FALSE, bins=0, test=c("standard","DGRP","FWW","asymptotic","iMK"), xlow=0, xhigh=1) { 
  
  ## Get PopFly data
  if (exists("PopFlyData") == TRUE) {
    data <- get("PopFlyData")
  } else {
    loadPopFly()
    data <- get("PopFlyData") }

  ## Check input variables
  ## Numer of arguments
  if (nargs() < 3 && nargs()) {
    stop("You must specify 3 arguments at least: genes, pops, recomb (T/F).\nIf test = asymptotic or test = iMK, you must specify xlow and xhigh values.") }
  
  ## Argument genes
  if (length(genes) == 0 || genes == "" || !is.character(genes)) {
    stop("You must specify at least one gene.") }
  if (!all(genes %in% data$Name) == TRUE) {
    difGenes <- setdiff(genes, data$Name)
    difGenes <- paste(difGenes, collapse=", ")
    stopMssg <- paste0("MKT data is not available for the requested gene(s).\nRemember to use FlyBase IDs (FBgn...)\nThe genes that caused the error are: ", difGenes, ".")
    stop(stopMssg) }
  
  ## Argument pops
  if (length(pops) == 0 || pops == "" || !is.character(pops)) {
    stop("You must specify at least one population.") }
  if (!all(pops %in% data$Pop) == TRUE) {
    correctPops <- c("AM","AUS","CHB","EA","EF","EG","ENA","EQA","FR","RAL","SA","SD","SP","USI","USW","ZI")
    difPops <- setdiff(pops, correctPops)
    difPops <- paste(difPops, collapse=", ")
    stopMssg <- paste0("MKT data is not available for the sequested populations(s).\nSelect among the following populations:\nAM, AUS, CHB, EA, EF, EG, ENA, EQA, FR, RAL, SA, SD, SP, USI, USW, ZI.\nThe populations that caused the error are: ", difPops, ".")
    stop(stopMssg) }
  
  ## Argument recomb
  if (recomb != TRUE && recomb != FALSE) {
    stop("Parameter recomb must be TRUE or FALSE.") }
  
  ## Argument bins
  if (recomb == TRUE) {
    if (!is.numeric(bins) || bins == 0  || bins == 1) {
      stop("If recomb = TRUE, you must specify the number of bins to use (> 1).") }
    if (bins > round(length(genes)/2)) {
      stop("Parameter bins > (genes/2). At least 2 genes for each bin are required.") }
  }
  if (recomb == FALSE && bins != 0) {
    warning("Parameter bins not used! (recomb=F selected)") }
  
  ## Argument test and xlow + xhigh (when necessary)
  if(missing(test)) {
    test <- "standard"
  }
  else if (test != "standard" && test != "DGRP" && test != "FWW" && test != "asymptotic" && test != "iMK") {
    stop("Parameter test must be one of the following: standard, DGRP, FWW, asymptotic, iMK")
  }
  if (length(test) > 1) {
    stop("Select only one of the following tests to perform: standard, DGRP, FWW, asymptotic, iMK") }
  if ((test == "standard" || test == "DGRP" || test == "FWW") && (xlow != 0 || xhigh != 1)) {
    warningMssgTest <- paste0("Parameters xlow and xhigh not used! (test = ",test," selected)")
    warning(warningMssgTest) }
  
  ## Arguments xlow, xhigh features (numeric, bounds...) checked in check_input()
  
  ## Perform subset
  subsetGenes <- data[(data$Name %in% genes & data$Pop %in% pops), ]
  subsetGenes$Name <- as.factor(subsetGenes$Name)
  subsetGenes <- droplevels(subsetGenes)
  
  ## If recomb analysis is selected
  if (recomb == TRUE) {
    
    ## Declare output list (each element 1 pop)
    outputList <- list()
    
    for (k in levels(subsetGenes$Pop)) {
      
      print(paste0("Population = ", k))
      ## Declare bins output list (each element 1 bin)
      outputListBins <- list()
      
      x <- subsetGenes[subsetGenes$Pop == k, ]
      x <- x[order(x$p_bp), ]
      
      ## create bins
      binsize <- round(nrow(x)/bins) ## Number of genes for each bin
      count <- 1; dat <- NULL
      
      for (i in 0:nrow(x)) {
        if (i%%binsize == 0) { ## Only if reminder of division = 0 (equally sized bins)
          i1 <- i + binsize - 1
          if (i == 0) {
            g1 <- x[i:binsize,]
            group <- count
            g1$Group <- group
            dat <- rbind(dat,g1)
            count <- count+1 }
          else if (i1 < nrow(x)) {
            g1 <- x[i:i1,]
            group <- count
            g1$Group <- group
            dat <- rbind(dat,g1)
            count <- count+1 }
        }
      }
      dat$Group <- as.factor(dat$Group)
      dat <- unique(dat)
      
      ## Iterate through each recomb bin
      for (j in levels(dat$Group)) {
        
        print(paste0("Recombination bin = ", j))
        x1 <- dat[dat$Group == j, ]
        
        ## Set counters to 0
        Pi <- c(0,0,0,0,0,0,0,0,0,0)
        P0 <- c(0,0,0,0,0,0,0,0,0,0)
        f <- seq(0.05,0.95,0.1)
        mi <- 0; m0 <- 0
        Di <- 0; D0 <- 0
        
        ## Group genes
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
        
        ## Check data inside each test!
        ## Perform test
        if(test == "standard") {
          output <- standard(daf, div) }
        else if(test == "DGRP") {
          output <- DGRP(daf, div) }
        else if(test == "FWW") {
          output <- FWW(daf, div) }
        else if(test == "asymptotic") {
          output <- asymptoticMK(daf, div, xlow, xhigh) }
        else if(test == "iMK") {
          output <- iMK(daf, div, xlow, xhigh) }
        
        ## Fill list with each bin
        outputListBins[[paste("Recombination bin = ",j)]] <- output
      }
      
      ## Fill list with each pop
      outputList[[paste("Population = ",k)]] <- outputListBins
    }
    
    ## Warning if some genes are lost. Bins must be equally sized.
    if (nrow(dat) != length(genes)) {
      missingGenes <- round(length(genes) - nrow(dat))
      genesNames <- as.vector(tail(x, missingGenes)$Name)
      genesNames <- paste(genesNames, collapse=", ")
      warningMssg <- paste0("The ",missingGenes," gene(s) with highest recombination rate estimates (", genesNames, ") was/were excluded from the analysis in order to get equally sized bins.\n")
      warning(warningMssg) }
    
    ## Return output
    cat("\n")
    return(outputList)
  }
  
  ## If NO recombination analysis selected
  else if (recomb == FALSE) {
    
    ## Declare output list (each element 1 pop)
    outputList <- list()
    
    for (i in levels(subsetGenes$Pop)) {
      
      print(paste0("Population = ", i))
      x <- subsetGenes[subsetGenes$Pop == i, ]
      
      ## Set counters to 0
      Pi <- c(0,0,0,0,0,0,0,0,0,0)
      P0 <- c(0,0,0,0,0,0,0,0,0,0)
      f <- seq(0.05,0.95,0.1)
      mi <- 0; m0 <- 0
      Di <- 0; D0 <- 0
      
      ## Group genes
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
      
      ## Check data inside each test!
      ## Perform test
      if(test == "standard") {
        output <- standard(daf, div) }
      else if(test == "DGRP") {
        output <- DGRP(daf, div) }
      else if(test == "FWW") {
        output <- FWW(daf, div) }
      else if(test == "asymptotic") {
        output <- asymptoticMK(daf, div, xlow, xhigh) }
      else if(test == "iMK") {
        output <- iMK(daf, div, xlow, xhigh) }
      
      ## Fill list with each pop
      outputList[[paste("Population = ",i)]] <- output
    }
    
    ## Return output
    cat("\n")
    return(outputList)
  }
}

