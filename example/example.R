##' MKT corrected with DGRP method
#'
#' Date = 30/11/2016
#' Author = Sergi Hervás, Marta Coronado
#'
#' The standard McDonald and Kreitman test (MKT) is used to detect the signature of selection at the molecular level. The MKT compares the amount of variation within a species (polymorphism, P) to the divergence (D) between species at two types of sites, one of which is putatively netral and used as the reference to detect selection at the other type of site. In the standard MKT, these sites are synonymous (putatively neutral, 0) and non-synonymous sites (selected sites, i) in a coding region. Under strict neutrality, the ratio of the number of selected and neutral polymorphic sites (Pi/P0) is equal to the ratio of the number of selected and neutral divergence sites (Di/D0).
# The null hypothesis of neutrality is rejected in a MKT when Di/D0 > Pi/P0. The excess of divergence relative to polymorphism for class i, is interpreted as adaptive selection for a subset of sites i. The fraction of adaptive fixations, α, is estimated from 1-(Pi/P0)(D0/Di). The significance of the test can be assesed with a Fisher exact test.
#' The estimate of α can be easily biased by the segregation of slightly deleterious non-synonymous substitutions. Specifically, slightly deleterious mutations tend to contribute more to polymorphism than to divergence, and thus, lead to an underestimation of alpha. Because adaptive mutations and weakly deleterious selection sct in opposite directions on the MKT, alpha and the fraction of substitutions that are sligholty deleterious, b, will be both underestimated when the two selection regimes occur. To take adaptive and slighlty deleterious mutations mutually into account, Pi, the count off segregatning sites in class i, should be seaprated into the number of neutral variants and the number of weakly deleterious variants, Pi = Pineutral + Pi weak del. Alpha is then estimated as 1-(Pineutral/P0)(D0/Di)
#'
#'
#' @param x dad file
#' @param y divergence file
#' 
#' @return NONE
#'
#' @examples
#' mkt_DGRP(x,y)

x <- read.table("~/MKT/iMKT/example/RAL_Chr2L.daf10.txt", header=T) #polymorphism file (DAF)
y <- read.table("~/MKT/iMKT/example/RAL_Chr2L_div.txt", header=T) #divergence and m file
w<-asymptoticMK(x, y, 0, 1)
w
## 3rd. perform analysis ##
# print("Time iMKT")
# for(i in 1:3){
#   print(system.time(iMK(x, y, 0, 1)))
#   w #check out results
# }
# print("Time PlotDAF")
# for (i in 1:3){
#   print(system.time(plotDAF(x))) #plots
# }

# print("Time plotIMKT")
# for (i in 1:3){
#   print(system.time(plotIMK(w$iMK)))
# }
# print("Time plotALPHA")
# for (i in 1:3){
#   print(system.time(plotALPHA(x, y, w$asymptoticMK, 0, 1)))
# }
# #check watchdog:
# # w <- iMK(x, y, 0, 1.4)


# xlow <- 0
# xhigh <- 1

# result <- iMK(x, y, xlow, xhigh)
# result




plotDAF(x)


# plotIMK(result$iMK)


plotALPHA(x, y, result$asymptoticMK, xlow, xhigh)

