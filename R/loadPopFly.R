#' @title Load PopFly dataset
#'
#' @description Load PopFly dataset with information regarding protein coding gene annotations
#'
#' @details This function loads PopFly data (Hervas et al. 2017 Bioinformatics, http://popfly.uab.cat/) into the current workspace. Data is stored in a dataframe named PopFlyData, which includes population genetics estimates (nucleotide diversity, divergence, basic tests of neutrality, recombination rates, etc.) regarding each protein coding gene for 16 worldwide wild-derived Drosophila melanogaster populations from the Drosophila Genome Nexus project (Lack et al. 2015 Genetics, Lack et al. 2016 MBE).
#'
#' @return PopFlyData object loaded in the workspace
#'
#' @examples
#' ## Load PopFly data if necessary. This process may take several seconds to complete.
#' # loadPopFly()
#' 
#' @import utils
#'
#' @export

loadPopFly <- function() {
  PopFlyData <- ""
  cat("Loading PopFly data into your workspace.\nThis process may take several seconds to complete, please be patient.\n")
  x <- read.table("http://popfly.uab.cat/files/genes/GenesData_recomb_comeron.tab", header=T)
  PopFlyData <<- x
}
