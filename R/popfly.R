#' Extract a list of genes from PopFly
#' Date = 30/11/2016
#' Author = Jesús Murga, Marta Coronado
#'
#'
#' @param datos PopFly dataset preload in package
#' @param genes Drosophila gene list
#' @param POP Drosophila gene population
#' @param recomb Retrieve genes by recombinations values
#' @param bin recombination value
#' @return None
#' @examples
#' @import utils
#' @export
#' 

################################NEED A FOLDER WITH FILES TO CHECK HOW THE OUTPUT IT'S RETURNED################################

subsetPopfly<-function(datos='data input',genes=c("gene1"),POP=c("population"), recomb=TRUE/FALSE,bin){
  listofdaf<-list()
  listofdivergence<-list()
  if (recomb == FALSE){
    if (length(genes)==0 || POP == "")
      stop("You need to specify at least one gene or one population!")
    
    else if(length(genes) == 1){
      for (j in 1:length(POP)){
        gene<-datos[datos[2]==genes & datos[1]==POP[j],]
        subsetGenes<-rbind(subsetGenes,gene)}
    }
    
    else if(length(genes)!=1){
      subsetGenes<-data.frame()
      for (i in 1:length(genes)){
        for (j in 1:length(POP)){
          gene<-datos[datos[2]==genes[i] & datos[1]==POP[j],]
          # print(gene)
          subsetGenes<-rbind(subsetGenes,gene)}
      }
    }
    
    else if(POP == 'ALL'){
      subsetGenes<-data.frame()
      for (i in 1:length(genes)){
        gene<-datos[datos[2]==genes[i],]
        subsetGenes<-rbind(subsetGenes,gene)}
    }
  }
  else if (recomb == TRUE){
    bin
    for (j in 1:length(POP)){}
    
  }
  return(subsetGenes)
}

# system.time(a<-subsetPopfly(popfly,genes=genes,POP=c("RAL","ZI","AM"),recomb = F))

