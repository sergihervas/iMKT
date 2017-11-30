### example procedure to execute the package.

## 1st. load all functions ##
#(...)

## 2nd. read data ##
#setwd("...")
x <- read.table("~/MKT/iMKT/example/RAL_Chr2L.daf10.txt", header=T) #polymorphism file (DAF)
y <- read.table("~/MKT/iMKT/example/RAL_Chr2L_div.txt", header=T) #divergence and m file
w<-asymptoticMK(x, y, 0, 1)
w
## 3rd. perform analysis ##
print("Time iMKT")
for(i in 1:3){
  print(system.time(iMK(x, y, 0, 1)))
  w #check out results
}
print("Time PlotDAF")
for (i in 1:3){
  print(system.time(plotDAF(x))) #plots
}

print("Time plotIMKT")
for (i in 1:3){
  print(system.time(plotIMK(w$iMK)))
}
print("Time plotALPHA")
for (i in 1:3){
  print(system.time(plotALPHA(x, y, w$asymptoticMK, 0, 1)))
}
#check watchdog:
# w <- iMK(x, y, 0, 1.4)


xlow <- 0
xhigh <- 1

result <- iMK(x, y, xlow, xhigh)
result




plotDAF(x)


plotIMK(result$iMK)


plotALPHA(x, y, result$asymptoticMK, xlow, xhigh)

