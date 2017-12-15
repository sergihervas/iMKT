#' @title multipleDatasets
#' 
#' @description  Perform any MK test using all files (or a subset of them) in a given directory.
#'
#' @details Files in directory must be named: file1*daf*, file1*divergence*, file2*daf*, file2*divergence*, ...
#'
#' @param directory directory (path/to/files/) where daf and divergence files are stored in your local machine
#' @param test which test to perform. Options include: standard (default), DGRP, FWW, asymptotic, iMK
#' @param fullAnalysis decide whether to analyze all files in directory or not (default=TRUE) 
#' @param idList used when fullAnalysis = F, list of IDs to analyze
#' @param xlow lower limit for asymptotic alpha fit (default=0)
#' @param xhigh higher limit for asymptotic alpha fit (default=1)
#'
#' @return None
#'
#' @examples 
#' ##example here
#' 
#' @import utils
#' @import stats
#'
#' @export

#multipleDatasets(directory="/home/sergi/testiMK/", test="standard", fullAnalysis=F, idList=c("id1","id2"))


multipleDatasets <- function(directory=directory, test=c("standard","DGRP","FWW","asymptotic","iMK"), xlow=0, xhigh=1, fullAnalysis=TRUE/FALSE, idList='NA') {
  
  ## Set working directory and list files (daf and divergence)
  wd <- directory; setwd(wd) 
  files <- list.files(wd) 
  subset_daf <- list.files(wd, pattern = "daf.*$")
  subset_divergence<-list.files(wd,pattern = "divergence.*$")
  
  ## Create output empty list
  result <- list()

  ## If fullAnalysis == F, check there is idList!
  if (fullAnalysis == FALSE && idList == 'NA') {
  	stop("You must specify a list of IDs to analyze (fullAnalysis = F selected).")
  }

  ## If idList, keep only matching names from subset_daf and subset_divergence
  if (fullAnalysis == FALSE && idList != 'NA') {
  	
    ## Analyze each ID
  	for (i in idList) {
      daf <- (grep(paste0(i,".daf"), subset_daf, value=T))
      divergence <- (grep(paste0(i,".divergence"), subset_divergence, value=T))
      daf <- read.table(daf, header=T)
      divergence <- read.table(divergence, header=T)
      
      ## Perform test
      if(test == "standard") {
        temp <- list("ID"=i, "StandardMK"=standard(daf, divergence))
        result[[i]] <- temp
      }
      else if(test == "DGRP") {
        temp <- list("ID"=i, "DGRP"=DGRP(daf, divergence))
        result[[i]] <- temp
      }
      else if(test == "FWW") {
        temp <- list("ID"=i, "FWW"=FWW(daf, divergence))
        result[[i]] <- temp
      }
      else if(test == "asymptotic") {
        temp <- list("ID"=i, "asymptotic"=asymptoticMK(daf, divergence, xlow, xhigh))
        result[[i]] <- temp
      }
      else if(test == "iMK") {
        temp <- list("ID"=i, "iMK"=iMK(daf, divergence, xlow, xhigh))
        result[[i]] <- temp
      }
  	}
    
    ## Return output
    return(result)
  }
  
  else if (fullAnalysis == TRUE) {
    names <- sub("*.daf.*", "", subset_daf)
    
    ## Analyze each ID
    for (i in 1:length(names)){
      daf <- (grep(paste0(names[i],".daf"), subset_daf, value=T))
      divergence<-(grep(paste0(names[i],".divergence"), subset_divergence, value=T))
      daf <- read.table(daf, header=T)
      divergence <- read.table(divergence, header=T)
      
      ## Perform test
      if(test == "standard") {
        temp <- list("ID"=i, "StandardMK"=standard(daf, divergence))
        result[[i]] <- temp
      }
      else if(test == "DGRP") {
        temp <- list("ID"=i, "DGRP"=DGRP(daf, divergence))
        result[[i]] <- temp
      }
      else if(test == "FWW") {
        temp <- list("ID"=i, "FWW"=FWW(daf, divergence))
        result[[i]] <- temp
      }
      else if(test == "asymptotic") {
        temp <- list("ID"=i, "asymptotic"=asymptoticMK(daf, divergence, xlow, xhigh))
        result[[i]] <- temp
      }
      else if(test == "iMK") {
        temp <- list("ID"=i, "iMK"=iMK(daf, divergence, xlow, xhigh))
        result[[i]] <- temp
      }
    }
    
    ## Return output
    return(result)
  }
}