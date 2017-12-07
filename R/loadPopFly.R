#' @title loadPopFly
#'
#' @description Load PopFly dataset
#'
#' @details This function loads PopFly data into the current workspace. Data is stored in a dataframe named PopFlyData.
#'
#' @return None
#'
#' @examples
#' loadPopFly()
#' 
#' @import utils
#'
#' @export
#' 

## load popfly data
loadPopFly <- function() {
  PopFlyData <- ""
  cat("Loading PopFly data into your workspace.\nThis process may take some seconds to complete, please be patient.")
  x <- read.table("http://popfly.uab.cat/files/genes/GenesData_recomb.tab", header=T)
  assign("PopFlyData", x, envir = .GlobalEnv)
}
