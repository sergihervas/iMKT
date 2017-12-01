x <- read.table("~/GitHub/iMKT/example/RAL_Chr2L.daf10.txt", header=T) #polymorphism file (DAF)
y <- read.table("~/GitHub/iMKT/example/RAL_Chr2L_div.txt", header=T) #divergence and m file
# 
w<-mkt_standard(x, y);w

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

f<-cmpfun(multiple_datasets)

for(i in 1:3){
  print(system.time(for(i in 1:1000){multiple_datasets(directory = "~/Test",test = "standard")}))
}

for(i in 1:3){
  print(system.time(for(i in 1:1000){f(directory = "~/Test",test = "standar")}))
}




