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

multiple_datasets<-function(directory="Directory",test="mkt type",fullanalyis=TRUE/FALSE,idlist="listfile"){
  
  wd<-directory;setwd(wd) #SET AND CREATE A WORKING DIRECTORY
  files<-list.files(wd) #LIST ALL FILES IN DIRECTORY
  subset_daf<-list.files(wd,pattern = "daf.*$")  #SEPARATE FILES IN DAF AND DIVERGENCE
  subset_divergence<-list.files(wd,pattern = "divergence.*$")  #SEPARATE FILES IN DAF AND DIVERGENCE

  result<-list() #EMPTY RESULT
  
  if (fullanalyis==FALSE & idlist!='NA'){
    list<-as.matrix(read.table(idlist))
    if(test=="standard"){
      for (i in list){
        daf<-(grep(paste0(i,".daf"),subset_daf,value=T))
        divergence<-(grep(paste0(i,".divergence"),subset_divergence,value = T))
        z<-read.table(daf,header = T)
        c<-read.table(divergence,header = T)
        # result<-append(result,)
        temp<-list("ID"=i,"StandarMKT"=mkt_standard(z,c))
        result<-append(result,temp)
      }
    }
    return(result)
  }
else{
  names<-sub("*.daf.*", "", subset_daf)
  if(test=="standard"){
    for (i in 1:length(names)){
      daf<-(grep(paste0(names[i],".daf"),subset_daf,value=T))
      divergence<-(grep(paste0(names[i],".divergence"),subset_divergence,value = T))
      z<-read.table(daf,header = T)
      c<-read.table(divergence,header = T)
      # result<-append(result,)
      temp<-list("ID"=names[i],"StandarMKT"=mkt_standard(z,c))
      result<-append(result,temp)
    }
  }
  
  else if(test=="FWW"){
    for (i in 1:length(names)){
      daf<-(grep(paste0(names[i],".daf"),subset_daf,value=T))
      divergence<-(grep(paste0(names[i],".divergence"),subset_divergence,value = T))
      z<-read.table(daf,header = T)
      c<-read.table(divergence,header = T)
      # result<-append(result,)
      temp<-list("ID"=i,"StandarMKT"=mkt_fww(z,c))
      result<-append(result,temp)
    }
  }
  else if(test=="DGRP"){
    for (i in 1:length(names)){
      daf<-(grep(paste0(names[i],".daf"),subset_daf,value=T))
      divergence<-(grep(paste0(names[i],".divergence"),subset_divergence,value = T))
      z<-read.table(daf,header = T)
      c<-read.table(divergence,header = T)
      # result<-append(result,)
      temp<-list("ID"=names[i],"StandarMKT"=mkt_DGRP(z,c))
      result<-append(result,temp)
    }
  }
  
}
  return(result)
 }

  
c<-multiple_datasets(directory = "~/Test",test = "standard",fullanalyis = TRUE, idlist = "NA")
c
