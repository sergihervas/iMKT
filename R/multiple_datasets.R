#' Execute any test  with all the files in a directory. The files need to be called id.daf.txt and id.divergence.txt
#' Date = 30/11/2016
#' Author = Jesús Murga, Marta Coronado
#'
#'
#' @param directory dad file
#' @param test divergence file
#'
#' @return None
#'
#' @examples
#' @import knitr 
#' @import utils
#' @import stats
#' @import ggplot2 
#' @import ggthemes
#' @import cowplot
#' @import grid 
#' @import gridExtra
#' @import scales
#' @import reshape2 
#' @export
#' 

multiple_datasets<-function(directory="Directory",test="MKT type"){
  wd<-directory
  setwd(wd)
  files<-list.files(wd)
  subset_daf<-list.files(wd,pattern = "daf.*$")
  subset_divergence<-list.files(wd,pattern = "divergence.*$")
  names<-sub("*.daf.*", "", subset_daf)

  if(test=="standard"){
    for (i in names){
      daf<-(grep(paste0(i,".daf"),subset_daf,value=T))
      divergence<-grep(paste0(i,".divergence"),subset_divergence,value = T)
      z<-read.table(daf,header = T)
      c<-read.table(divergence,header = T)
      mkt_standard(z,c)}
  }
  else if(test=="FWW"){
    for (i in names){
      daf<-(grep(paste0(i,".daf"),subset_daf,value=T))
      divergence<-grep(paste0(i,".divergence"),subset_divergence,value = T)
      z<-read.table(daf,header = T)
      c<-read.table(divergence,header = T)
      mkt_standard(z,c)}
  }
  else if(test=="DGRP"){
    for (i in names){
      daf<-(grep(paste0(i,".daf"),subset_daf,value=T))
      divergence<-grep(paste0(i,".divergence"),subset_divergence,value = T)
      z<-read.table(daf,header = T)
      c<-read.table(divergence,header = T)
      mkt_standard(z,c)}
  }
}





