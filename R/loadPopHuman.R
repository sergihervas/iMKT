#' @title Load PopHuman dataset
#'
#' @description Load PopHuman dataset with information regarding protein coding gene annotations
#'
#' @details This function loads PopHuman data (Mulet et al. 2017 NAR, http://pophuman.uab.cat/) into the current workspace. Data is stored in a dataframe named PopHumanData, which includes population genetics estimates (nucleotide diversity, divergence, basic tests of neutrality, recombination rates, etc.) regarding each protein coding gene for 26 worldwide Homo sapiens populations from the 1000 Genomes Project (The 1000 Genomes Project Consortium 2012 Nature, The 1000 Genomes Project Consortium 2015 Nature).
#' 
#' @return PopHumanData object loaded in the workspace
#'
#' @examples
#' ## Load PopHuman data if necessary. This process may take several seconds to complete.
#' # loadPopHuman()
#' 
#' @import utils
#'
#' @keywords PopData
#' @export

loadPopHuman <- function() {
  PopHumanData <- ""
  cat("Loading PopHuman data into your workspace.\nThis process may take several seconds to complete, please be patient.\n")
  x <- read.table("http://pophuman.uab.cat/files/genes/GenesData_ALL_iMKT.tab", header=T,sep='\t')
  PopHumanData <<- x
}
