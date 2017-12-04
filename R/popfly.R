#' Extract a list of genes from PopFly
#' Date = 30/11/2016
#' Author = Jesús Murga, Marta Coronado
#'
#'
#' @param Genes Drosophila gene list
#' @return None
#'
#' @examples
#' @import utils
#' @export
#' 

################################NEED A FOLDER WITH FILES TO CHECK HOW THE OUTPUT IT'S RETURNED################################

subsetPopfly<-function(datos='data input',genes=c("gene1"),POP="population"){
  listofdaf<-list()
  listofdivergence<-list()
  if (length(genes)==0 || POP == "")
    stop("You need to specify at least one gene or one population!")
  else if(length(genes)==1){
    gene<-datos[datos[2] == genes & datos[1]==POP,]
    return(gene)
  }
  else if(length(genes)!=1){
      subsetGenes<-data.frame()
      for (i in 1:length(genes)){
        gene<-datos[datos[2]==genes[i] & datos[1]==POP,]
        print(gene)
        subsetGenes<-rbind(subsetGenes,gene)}
      return(subsetGenes)
    }
}

subsetPopfly(popfly,genes=c("FBgn0000055","FBgn0000055"),POP="RAL")

