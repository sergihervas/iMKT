## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
	echo = TRUE,
	fig.align = "center",
	fig.height = 7,
	fig.width = 7,
	collapse = TRUE,
	comment = "#>"
)

## ----Loading the package, echo=TRUE--------------------------------------
# install.packages("iMKT") ##If CRAN  
# library(devtools) ##If github
# install_github("sergihervas/iMKT")  
  
library(iMKT)

## ----Checking the example data, echo=TRUE--------------------------------
mydafdata
mydivergencedata

exampleDaf<-mydafdata
exampleDiverge<-mydivergencedata

## ----Standard MKT, echo=TRUE---------------------------------------------
standard<-standard(daf = mydafdata,divergence = mydivergencedata)

## ----Result, echo=TRUE---------------------------------------------------
standard$alpha.symbol
standard$`Fishers exact test P-value`
standard$`MKT table`

## ----FWW, echo=TRUE------------------------------------------------------
FWW(daf = mydafdata,divergence = mydivergencedata)

## ----Result FWW, echo=TRUE-----------------------------------------------
methodFWW<-FWW(daf = mydafdata,divergence = mydivergencedata)
methodFWW$Results
methodFWW$`MKT tables`
methodFWW$`Divergence metrics`

## ----FWW plot, echo=TRUE-------------------------------------------------
methodFWW<-FWW(daf = mydafdata, divergence = mydivergencedata,list_cutoff=c(0, 0.05, 0.1,0.15,0.2),
               plot=TRUE)
savePlotInVariable<-methodFWW$Graph
savePlotInVariable

## ----DGRP, echo=TRUE-----------------------------------------------------
DGRP(daf = mydafdata, divergence = mydivergencedata)

## ----Result DGRP, echo=TRUE----------------------------------------------
methodDGRP<-DGRP(daf = mydafdata, divergence = mydivergencedata)
methodDGRP$Results
methodDGRP$Fractions
methodDGRP$`MKT tables`

## ---- echo=TRUE----------------------------------------------------------
methodDGRP<-DGRP(daf = mydafdata, divergence = mydivergencedata,list_cutoff=c(0, 0.05,0.2),plot=TRUE)
savePlotInVariable<-methodDGRP$Graph
savePlotInVariable

## ----Asymptotic MKT, echo=TRUE-------------------------------------------
asymptoticMK(daf = mydafdata, divergence = mydivergencedata, xlow = 0,xhigh = 0.9)

## ----iMK, echo=TRUE------------------------------------------------------
iMK(daf = mydafdata, divergence = mydivergencedata, xlow = 0, xhigh = 0.9)

## ---- echo=TRUE----------------------------------------------------------
methodiMK<-iMK(daf = mydafdata,divergence = mydivergencedata, xlow = 0, xhigh = 0.9)
methodiMK$Results
methodiMK$`Divergence metrics`
methodiMK$`MKT tables`

## ----fig.width=7,fig.height=9,echo=TRUE----------------------------------
methodiMK<-iMK(daf = mydafdata, divergence = mydivergencedata, xlow = 0, xhigh = 0.9 ,plot=TRUE)
savePlotInVariable<-methodiMK$Graph
savePlotInVariable

## ----multipleDatasets, eval=FALSE, include=FALSE-------------------------
#  # multipleDatasets(directory = directory, test = c("standard", "DGRP", "FWW","asymptotic", "iMK"), xlow = 0, xhigh = 1, fullAnalysis = TRUE/FALSE,idList = "NA")

