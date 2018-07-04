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
## Load the iMKT library
library(iMKT)

## Load PopFlyData
loadPopFly()

## Check data object
ls()
names(PopFlyData)

## ----PopFly data retrieve automated no recomb, echo=TRUE-----------------
PopFlyAnalysis(genes=c("FBgn0000055","FBgn0003016"), pops=c("RAL","ZI"), recomb=F, test="DGRP", plot=TRUE)

## ----PopFly data retrieve automated, echo=TRUE---------------------------
geneList <- c("FBgn0053196","FBgn0086906","FBgn0261836","FBgn0031617","FBgn0260965",
              "FBgn0028899","FBgn0052580","FBgn0036181","FBgn0263077","FBgn0013733",
              "FBgn0031857","FBgn0037836")
PopFlyAnalysis(genes=geneList , pops="RAL", recomb=T, bins=2, test="iMKT", xlow=0, xhigh=0.9, plot=TRUE)
rm(geneList)

