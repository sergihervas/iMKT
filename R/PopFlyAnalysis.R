#' @title iMKT using PopFly data
#'
#' @description Perform any MKT method using a subset of PopFly data defined by custom genes and populations lists
#'
#' @details Execute any MKT method (standardMKT, FWW, DGRP, asymptoticMKT, iMKT) using a subset of PopFly data defined by custom genes and populations lists. It uses the dataframe PopFlyData, which can be already loaded in the workspace (using loadPopFly()) or is directly loaded when executing this function. It also allows deciding whether to analyze genes groupped by recombination bins or not, using recombination rate estimates from Comeron et al. 2012 Plos Genetics. 
#' 
#' @param genes list of genes to analyze
#' @param pops list of populations to analyze
#' @param recomb group genes according to recombination values (TRUE/FALSE)
#' @param bins number of recombination bins to compute (mandatory if recomb=TRUE)
#' @param test which test to perform. Options include: standardMKT (default), DGRP, FWW, asymptoticMKT, iMKT
#' @param xlow lower limit for asymptotic alpha fit (default=0)
#' @param xhigh higher limit for asymptotic alpha fit (default=1)
#' @param plot report plot (optional). Default is FALSE
#' 
#' @return List of lists with the default test output for each selected population (and recombination bin when defined)
#'
#' @examples
#' ## List of genes
#' mygenes <- c("FBgn0053196", "FBgn0086906", "FBgn0261836", "FBgn0031617", 
#'              "FBgn0260965", "FBgn0028899", "FBgn0052580", "FBgn0036181",
#'              "FBgn0263077", "FBgn0013733", "FBgn0031857", "FBgn0037836")
#' ## Perform analyses
#' PopFlyAnalysis(genes=mygenes, pops="RAL", recomb=FALSE, test="iMKT", xlow=0, xhigh=0.9, plot=TRUE)
#' PopFlyAnalysis(genes=mygenes, pops=c("RAL","ZI"), recomb=TRUE, bins=3, test="DGRP", plot=FALSE)
#' 
#' @import utils
#' @import stats
#'
#' @keywords PopData
#' @export

PopFlyAnalysis <- function(genes=c("gene1","gene2","..."), pops=c("pop1","pop2","..."), recomb=TRUE/FALSE, bins=0, test=c("standardMKT","DGRP","FWW","asymptoticMKT","iMKT"), xlow=0, xhigh=1, plot=FALSE) { 
  
  ## Get PopFly data
  if (exists("PopFlyData") == TRUE) {
    data <- get("PopFlyData")
  } else {
    loadPopFly()
    data <- get("PopFlyData") }

  ## Check input variables
  ## Numer of arguments
  if (nargs() < 3 && nargs()) {
    stop("You must specify 3 arguments at least: genes, pops, recomb (T/F).\nIf test = asymptoticMKT or test = iMKT, you must specify xlow and xhigh values.") }
  
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
    test <- "standardMKT"
  }
  else if (test != "standardMKT" && test != "DGRP" && test != "FWW" && test != "asymptoticMKT" && test != "iMKT") {
    stop("Parameter test must be one of the following: standardMKT, DGRP, FWW, asymptoticMKT, iMKT")
  }
  if (length(test) > 1) {
    stop("Select only one of the following tests to perform: standardMKT, DGRP, FWW, asymptoticMKT, iMKT") }
  if ((test == "standardMKT" || test == "DGRP" || test == "FWW") && (xlow != 0 || xhigh != 1)) {
    warningMssgTest <- paste0("Parameters xlow and xhigh not used! (test = ",test," selected)")
    warning(warningMssgTest) }
  
  ## Arguments xlow, xhigh features (numeric, bounds...) checked in checkInput()
  
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
      x <- x[order(x$cM_Mb), ]
      
      ## create bins
      binsize <- round(nrow(x)/bins) ## Number of genes for each bin
      count <- 1
      x$Group <- ""
      dat <- x[FALSE, ] ## Create df with colnames
      
      for (i in 0:nrow(x)) {
        if (i%%binsize == 0) { ## Only if reminder of division = 0 (equally sized bins)
          i1 <- i + binsize
          if (i == 0) {
            g1 <- x[i:binsize,]
            group <- count
            g1$Group <- group
            dat[i:binsize,] <- g1
            count <- count+1 }
          else if (i1 <= nrow(x)) {
            ii <- i+1
            g1 <- x[ii:i1,]
            group <- count
            g1$Group <- group
            dat[ii:i1,] <- g1
            count <- count+1 }
        }
      }
      dat$Group <- as.factor(dat$Group)

      ## Iterate through each recomb bin
      for (j in levels(dat$Group)) {
        print(paste0("Recombination bin = ", j))
        x1 <- dat[dat$Group == j, ]
        
        ## Recomb stats from bin j
        numGenes <- nrow(x1)
        minRecomb <- min(x1$cM_Mb, na.rm=T)
        medianRecomb <- median(x1$cM_Mb, na.rm=T)
        meanRecomb <- mean(x1$cM_Mb, na.rm=T)
        maxRecomb <- max(x1$cM_Mb, na.rm=T)
        recStats <- cbind(numGenes,minRecomb,medianRecomb,meanRecomb,maxRecomb)
        recStats <- as.data.frame(recStats)
        recStats <- list("Recombination bin Summary"=recStats)
        
        ## Set counters to 0
        Pi <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
        P0 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
        f <- seq(0.025,0.975,0.05)
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
        
        ## Transform daf20 to daf10 (faster fitting) for asymptoticMKT and iMKT
        if (nrow(daf) == 20) {
          daf1 <- daf
          daf1$daf10 <- sort(rep(seq(0.05,0.95,0.1),2)) ## Add column with the daf10 frequencies
          daf1 <- daf1[c("daf10","Pi","P0")] ## Keep new frequencies, Pi and P0
          daf1 <- aggregate(. ~ daf10, data=daf1, FUN=sum)  ## Sum Pi and P0 two by two based on daf
          colnames(daf1)<-c("daf","Pi","P0") ## Final daf columns name
        }
        ## Perform test
        if(test == "standardMKT") {
          output <- standardMKT(daf, div) 
          output <- c(output, recStats) } ## Add recomb summary for bin j
        else if(test == "DGRP" && plot == FALSE) {
          output <- DGRP(daf, div) 
          output <- c(output, recStats) }
        else if(test == "DGRP" && plot == TRUE) {
          output <- DGRP(daf, div, plot=TRUE) 
          output <- c(output, recStats) }
        else if(test == "FWW" && plot == FALSE) {
          output <- FWW(daf, div)           
          output <- c(output, recStats) }
        else if(test == "FWW" && plot == TRUE) {
          output <- FWW(daf, div, plot=TRUE)           
          output <- c(output, recStats) }
        else if(test == "asymptoticMKT") {
          output <- asymptoticMKT(daf1, div, xlow, xhigh) 
          output <- c(output, recStats) }
        else if(test == "iMKT" && plot == FALSE) {
          output <- iMKT(daf1, div, xlow, xhigh)
          output <- c(output, recStats) }
        else if(test == "iMKT" && plot == TRUE) {
          output <- iMKT(daf1, div, xlow, xhigh, plot=TRUE)
          output <- c(output, recStats) }
        
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
      Pi <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
      P0 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
      f <- seq(0.025,0.975,0.05)
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
      
      ## Transform daf20 to daf10 (faster fitting) for asymptoticMKT and iMKT
      if (nrow(daf) == 20) {
        daf1 <- daf
        daf1$daf10 <- sort(rep(seq(0.05,0.95,0.1),2)) ## Add column with the daf10 frequencies
        daf1 <- daf1[c("daf10","Pi","P0")] ## Keep new frequencies, Pi and P0
        daf1 <- aggregate(. ~ daf10, data=daf1, FUN=sum)  ## Sum Pi and P0 two by two based on daf
        colnames(daf1)<-c("daf","Pi","P0") ## Final daf columns name
      }
      
      ## Perform test
      if(test == "standardMKT") {
        output <- standardMKT(daf, div) }
      else if(test == "DGRP" && plot == FALSE) {
        output <- DGRP(daf, div) }
      else if(test == "DGRP" && plot == TRUE) {
        output <- DGRP(daf, div, plot=TRUE) }
      else if(test == "FWW" && plot == FALSE) {
        output <- FWW(daf, div) }
      else if(test == "FWW" && plot == TRUE) {
        output <- FWW(daf, div, plot=TRUE) }
      else if(test == "asymptoticMKT") {
        output <- asymptoticMKT(daf1, div, xlow, xhigh) }
      else if(test == "iMKT" && plot == FALSE) {
        output <- iMKT(daf1, div, xlow, xhigh) }
      else if(test == "iMKT" && plot == TRUE) {
        output <- iMKT(daf1, div, xlow, xhigh, plot=TRUE) }
      
      ## Fill list with each pop
      outputList[[paste("Population = ",i)]] <- output
    }
    
    ## Return output
    cat("\n")
    return(outputList)
  }
}

