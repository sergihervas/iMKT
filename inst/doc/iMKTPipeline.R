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
head(mydafdata)

## Sample divergence data
mydivergencedata

## ----Standard MKT, echo=TRUE---------------------------------------------
standard(daf=mydafdata, divergence=mydivergencedata)

## ----FWW, echo=TRUE------------------------------------------------------
FWW(daf=mydafdata, divergence=mydivergencedata)

## ----FWW plot, echo=TRUE, fig.width=6, fig.height=4----------------------
FWW(daf=mydafdata, divergence=mydivergencedata, list_cutoff=c(0.05, 0.15,0.25,0.35), plot=TRUE)

## ----DGRP, echo=TRUE-----------------------------------------------------
DGRP(daf=mydafdata, divergence=mydivergencedata)

## ---- echo=TRUE, fig.width=6, fig.height=6-------------------------------
DGRP(daf=mydafdata, divergence=mydivergencedata, list_cutoff=c(0.05, 0.15,0.25,0.35), plot=TRUE)

## ----Asymptotic MKT, echo=TRUE-------------------------------------------
asymptoticMK(daf=mydafdata, divergence=mydivergencedata, xlow=0, xhigh=0.9)

## ----iMK, echo=TRUE, fig.width=6, fig.height=9---------------------------
iMK(daf=mydafdata, divergence=mydivergencedata, xlow=0, xhigh=0.9, plot=TRUE)

