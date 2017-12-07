#' @title loadPopHuman
#'
#' @description Load PopHuman dataset
#'
#' @details This function loads PopHuman data (http://pophuman.uab.cat/) into the current workspace. Data is stored in a dataframe named PopFlyData.
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

loadPopHuman <- function() {
  PopHumanData <- ""
  cat("Loading PopHuman data into your workspace.\nThis process may take some seconds to complete, please be patient.")
  x <- read.table("http://pophuman.uab.cat/files/genes/GenesData_recomb.tab", header=T)
  assign("PopHumanData", x, envir = .GlobalEnv)
}
