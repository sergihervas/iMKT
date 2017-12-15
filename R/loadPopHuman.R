#' @title loadPopHuman
#'
#' @description Load PopHuman dataset
#'
#' @details This function loads PopHuman data (http://pophuman.uab.cat/) into the current workspace. Data is stored in a dataframe named PopFlyData.
#'
#' @return None
#'
#' @examples
#' # loadPopHuman()
#' 
#' @import utils
#'
#' @export

loadPopHuman <- function() {
  PopHumanData <- ""
  cat("Loading PopHuman data into your workspace.\nThis process may take some seconds to complete, please be patient.\n")
  x <- read.table("http://pophuman.uab.cat/files/genes/GenesData_ALL.tab", header=T)
  PopHumanData <<- x
}
