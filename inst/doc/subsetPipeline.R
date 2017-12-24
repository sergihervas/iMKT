## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
	echo = TRUE,
	fig.align = "center",
	fig.height = 7,
	fig.width = 7,
	collapse = TRUE,
	comment = "#>"
)

## ----popfly data, echo=TRUE, fig.width=7---------------------------------
library(iMKT)
loadPopFly()
ls() ## new object created
names(PopFlyData)

## ----PopFly data manual retrieve, echo=FALSE, eval=FALSE-----------------
#  ## Preparing RAL Adh
#  adhRAL<-dataPopfly[dataPopfly$Name=='FBgn0000055' & dataPopfly$Pop=='RAL',]
#  
#  adhRAL$DAF0f <- as.character(adhRAL$DAF0f); adhRAL$DAF4f <- as.character(adhRAL$DAF4f)
#  adhRAL0f<-unlist(strsplit(adhRAL$DAF0f, split=";"))
#  adhRAL4f <- unlist(strsplit(adhRAL$DAF4f, split=":"))
#  adhRAL0f <- as.numeric(adhRAL0f); adhRAL4f <- as.numeric(adhRAL0f)
#  
#  f <- seq(0.05,0.95,0.1)
#  mi <- adhRAL$mi; m0 <- adhRAL$m0
#  Di <- adhRAL$di; D0 <- adhRAL$d0
#  
#  dafAdhRAL <- cbind(f, adhRAL0f, adhRAL4f); dafAdhRAL <- as.data.frame(dafAdhRAL)
#  names(dafAdhRAL) <- c("daf","Pi","P0")
#  divAdhRAL <- cbind(mi, Di, m0, D0); divAdhRAL <- as.data.frame(divAdhRAL)
#  names(divAdhRAL) <- c("mi","Di","m0","D0")
#  
#  ## Preparing ZI Adh
#  adhZI<-dataPopfly[dataPopfly$Name=='FBgn0000055' & dataPopfly$Pop=='ZI',]
#  
#  adhZI$DAF0f <- as.character(adhZI$DAF0f); adhZI$DAF4f <- as.character(adhZI$DAF4f)
#  adhZI0f<-unlist(strsplit(adhZI$DAF0f, split=";"))
#  adhZI4f <- unlist(strsplit(adhZI$DAF4f, split=":"))
#  adhZI0f <- as.numeric(adhZI0f); adhZI4f <- as.numeric(adhZI0f)
#  
#  f <- seq(0.05,0.95,0.1)
#  mi <- adhZI$mi; m0 <- adhZI$m0
#  Di <- adhZI$di; D0 <- adhZI$d0
#  
#  dafAdhZI <- cbind(f, adhZI0f, adhZI4f); dafAdhZI <- as.data.frame(dafAdhZI)
#  names(dafAdhZI) <- c("daf","Pi","P0")
#  divAdhZI <- cbind(mi, Di, m0, D0); divAdhZI <- as.data.frame(divAdhZI)
#  names(divAdhZI) <- c("mi","Di","m0","D0")

## ----PopFly data retrieve Adh RAL, echo=FALSE, eval=F--------------------
#  standard(daf = dafAdhRAL, divergence = divAdhRAL)
#  DGRP(daf = dafAdhRAL, divergence = divAdhRAL,plot = TRUE)

## ----PopFly data retrieve Adh ZI, echo=F, eval=F-------------------------
#  standard(daf = dafAdhZI, divergence = divAdhZI)
#  DGRP(daf = dafAdhZI, divergence = divAdhZI,plot = TRUE)

## ----PopFly data retrieve automated no recomb, echo=TRUE-----------------
PopFlyAnalysis(genes=c("FBgn0000055","FBgn0003016"), pops=c("RAL","ZI"), recomb=F, test="DGRP")

## ----PopFly data retrieve automated, echo=TRUE---------------------------
geneList <- c("FBgn0053196","FBgn0086906","FBgn0261836","FBgn0031617","FBgn0260965",
              "FBgn0028899","FBgn0052580","FBgn0036181","FBgn0263077","FBgn0013733",
              "FBgn0031857","FBgn0037836")

PopFlyAnalysis(genes=geneList , pops=c("RAL"), recomb=T, bins=2, test="iMK", xlow=0, xhigh=0.9)

