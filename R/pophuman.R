#' Extract a list of genes from PopHuman
#' Date = 30/11/2016
#' Author = Jesús Murga, Marta Coronado
#'
#'
#' @param Genes Human gene list
#' @return None
#'
#' @examples
#' @import utils
#' @export
#' 

################################NEED A FOLDER WITH FILES TO CHECK HOW THE OUTPUT IT'S RETURNED################################
# popfly<-function(Genes=c("gene1","gene2")){
#   listofdaf<-list()
#   listofdivergence<-list()
#   URL<-""
#   if (length(Genes==0))
#     stop("You need to specify at least one gene!")
#   else if(length(Genes)==1){
#     temporal <- tempfile()

#     daf<-downloadfile(paste0(URL,gene,"daf"),temporal)
#     divergence<-downloadfile(paste0(URL,gene,"divergence"),temporal) #save file in temporal localization and save a variable with the content

#     return(daf)
#     return(divergence)
#   }
#     else if(length(Genes)!=1){
#       temporal <- tempfile()

#       for (i in Genes){
        
#         filesdaf<-downloadfile((paste0(URL,i,"daf")) #save file in temporal localization and save a variable with the content
#         filesdivergence<-downloadfile((paste0(URL,i,"daf"))
                                       
#         listofdaf<-append(listofdataframes,filesdaf)
#         listofdivergence<-append(listofdataframes,filesdivergence)
#       }
#     }
#   return(listofdataframes)
# }

