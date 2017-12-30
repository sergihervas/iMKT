## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
	echo = TRUE,
	fig.align = "center",
	fig.height = 7,
	fig.width = 7,
	collapse = TRUE,
	comment = "#>"
)

## ----load package and see sample data, echo=T----------------------------
## Load package
# install.packages("devtools")
# devtools::install_github("sergihervas/iMKT")
library(iMKT)

## Sample daf data
head(myDafData)

## Sample divergence data
myDivergenceData

## ----Standard MKT, echo=TRUE---------------------------------------------
standardMK(daf=myDafData, divergence=myDivergenceData)

## ----FWW, echo=TRUE------------------------------------------------------
FWW(daf=myDafData, divergence=myDivergenceData)

## ----FWW plot, echo=TRUE, fig.width=6, fig.height=4----------------------
FWW(daf=myDafData, divergence=myDivergenceData, listCutoffs=c(0.05, 0.15,0.25,0.35), plot=TRUE)

## ----DGRP, echo=TRUE-----------------------------------------------------
DGRP(daf=myDafData, divergence=myDivergenceData)

## ---- echo=TRUE, fig.width=6, fig.height=6-------------------------------
DGRP(daf=myDafData, divergence=myDivergenceData, listCutoffs=c(0.05, 0.15,0.25,0.35), plot=TRUE)

## ----Asymptotic MKT, echo=TRUE-------------------------------------------
asymptoticMK(daf=myDafData, divergence=myDivergenceData, xlow=0, xhigh=0.9)

## ----iMK, echo=TRUE, fig.width=6, fig.height=9---------------------------
iMK(daf=myDafData, divergence=myDivergenceData, xlow=0, xhigh=0.9, plot=TRUE)

## ----summary, echo=FALSE-------------------------------------------------
results <- data.frame("Standard"=0.2365, "FWW_0.05"=0.5409, "FWW_0.1"=0.5798, "DGRP_0.05"=0.4249, "DGRP_0.1"=0.4249, "asymptotic_iMK"=0.6259)
knitr::kable(results, align="c")
rm(results)

