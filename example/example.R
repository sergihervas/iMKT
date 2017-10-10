### example procedure to execute the package.

## 1st. load all functions 
#(...)

## 2nd. read data
#setwd("")
x <- read.table("RAL_Chr2L.txt", header=T) #polymorphism file (DAF)
y <- read.table("RAL_Chr2L_div.txt", header=T) #divergence and m file

## 3rd. perform analysis
w <- iMK(x, y, 0, 1)
#check watchdog:
# w <- iMK(x, y, 0, 1.4)
plot1(w$iMK)
plot2(x, y, w$asymptoticMK, 0, 1)
