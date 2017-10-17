### example procedure to execute the package.

## 1st. load all functions ##
#(...)

## 2nd. read data ##
#setwd("...")
x <- read.table("RAL_Chr2L.txt", header=T) #polymorphism file (DAF)
y <- read.table("RAL_Chr2L_div.txt", header=T) #divergence and m file

## 3rd. perform analysis ##
w <- iMK(x, y, 0, 1)
w #check out results
plotDAF(x) #plots
plotIMK(w$iMK)
plotALPHA(x, y, w$asymptoticMK, 0, 1)

#check watchdog:
# w <- iMK(x, y, 0, 1.4)


xlow <- 0
xhigh <- 1

result <- iMK(x, y, xlow, xhigh)
result




plotDAF(x)


plotIMK(result$iMK)


plotALPHA(x, y, result$asymptoticMK, xlow, xhigh)

