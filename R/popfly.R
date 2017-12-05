#' Extract a list of genes from PopFly
#' Date = 30/11/2016
#' Author = Jesús Murga, Marta Coronado, Sergi Hervas
#'
#'
#' @param genes Drosophila gene list
#' @param pops Drosophila populations
#' @param recomb group genes according to recombination values (must specify bin)
#' @param bin number of recombination bins to compute
#' @return None
#' @examples
#' @import utils
#' @export
#' 

################################NEED A FOLDER WITH FILES TO CHECK HOW THE OUTPUT IS RETURNED################################

#data should be pre-loaded, so the function takes as input (genes, pops, recomb, bins). Altough data is preload in the package you need to pass it in the function
subsetPopFly <- function(data="data",genes=c("gene1","gene2","..."), pops=c("pop1","pop2","..."), recomb=TRUE/FALSE, bins){ 
  daflist<-list()
  divlist<-list()
  if (length(genes) == 0 || length(pops) == 0){ #check this...
      stop("You must specify at least one gene and one population.")
    }
  
  else{
    subsetGenes <- data[(data$Name %in% genes & data$Pop %in% pops),]
    # print(subsetGenes)
    if(recomb == FALSE) {
      #do not keep recomb cols
      subsetGenes <- subsetGenes[,-17] #dnt remember col name
      # print(subsetGenes)
      #subsetGene contains the desired info. now, transform format and execute iMK
      # subsetGenes <- droplevels(subsetGenes)
      for (i in subsetGenes$Name){
        print(i)
      
        x <- subsetGenes[subsetGenes$Name == i]
        #x file
        print(x$DAF4f[2])
        print('+++++++++++++++')
        x$DAF0f <- as.character(x$DAF0f)
        daf0f <- unlist(strsplit(x$DAF0f, split=";"))
        print(daf0f)
        x$DAF4f <- as.character(x$DAF4f)
        daf4f <- unlist(strsplit(x$DAF4f, split=";"))
        x1 <- cbind(daf0f,daf4f)
        print(x1)
        x1 <- as.data.frame(x1)
        x1$daf <- seq(0.1,1,0.1)
        x1<-x1[,c("daf","daf0f","daf4f")]
        # print(x1)
        # daflist[[i]]<-x1
        # #y file (must change input to add m values)
        # y <- rbind(subsetGenes$m0f, subsetGenes$D0f, subsetGenes$m4f, subsetGenes$D4f) #check order
        # y <- as.data.frame(y)
        }
    }
  }
  # return(x1)
}
  
