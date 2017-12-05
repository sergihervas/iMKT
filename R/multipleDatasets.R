#' Execute any test  with all the files in a directory. The files need to be called id.daf.txt and id.divergence.txt
#' Date = 30/11/2016
#' Author = Jesús Murga, Marta Coronado
#'
#'
#' @param directory dad file
#' @param test divergence file
#' @param fullanalyis TRUE indicates execute the analysis for all the files in a folder 
#' @param idlist List of genes or datasets to analyze
#' @return None
#'
#' @examples
#' @import knitr 
#' @import utils
#' @import stats
#' @import grid 
#' @import gridExtra
#' @import scales
#' @import reshape2 
#' @import ggplot2
#' @importFrom ggthemes theme_foundation
#' @importFrom cowplot plot_grid
#' @export
#' 

multipleDatasets<-function(directory="Directory",test=c("DGRP","FWW"),fullanalyis=TRUE/FALSE,idlist="listfile"){
  
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
          temp<-list("ID"=i,"StandarMKT"=standard(z,c))
          result[[i]]<-temp}
          }
  }
# 
else{
  names<-sub("*.daf.*", "", subset_daf)
  if(test=="standard" & idlist == 'NA'){
    for (i in 1:length(names)){
      daf<-(grep(paste0(names[i],".daf"),subset_daf,value=T))
      divergence<-(grep(paste0(names[i],".divergence"),subset_divergence,value = T))
      z<-read.table(daf,header = T)
      c<-read.table(divergence,header = T)
      # result<-append(result,)
      temp<-list("ID"=names[i],"StandarMKT"=standard(z,c))
      result[[names[i]]]<-temp
    }
  }

  else if(test=="FWW"){
    for (i in 1:length(names)){
      daf<-(grep(paste0(names[i],".daf"),subset_daf,value=T))
      divergence<-(grep(paste0(names[i],".divergence"),subset_divergence,value = T))
      z<-read.table(daf,header = T)
      c<-read.table(divergence,header = T)
      # result<-append(result,)
      temp<-list("ID"=i,"FWW"=FWW(z,c))
      result[[names[i]]]<-temp
    }
  }
  else if(test=="DGRP"){
    for (i in 1:length(names)){
      daf<-(grep(paste0(names[i],".daf"),subset_daf,value=T))
      divergence<-(grep(paste0(names[i],".divergence"),subset_divergence,value = T))
      z<-read.table(daf,header = T)
      c<-read.table(divergence,header = T)
      # result<-append(result,)
      temp<-list("ID"=names[i],"DGRP"=DGRP(z,c))
      result[[names[i]]]<-temp
    }
  }
  else if(test=="iMK"){
    for (i in 1:length(names)){
      daf<-(grep(paste0(names[i],".daf"),subset_daf,value=T))
      divergence<-(grep(paste0(names[i],".divergence"),subset_divergence,value = T))
      z<-read.table(daf,header = T)
      c<-read.table(divergence,header = T)
      # result<-append(result,)
      temp<-list("ID"=names[i],"iMKT"=iMK(z,c))
      result[[names[i]]]<-temp
    }
  }
}
  return(result)
}

  
# multiple_datasets(directory = "~/MKT/Test",test = "standard",fullanalyis = TRUE, idlist = "NA")
# multiple_datasets(directory = "~/MKT/Test",test = "standard",fullanalyis = FALSE, idlist = "idlist")
