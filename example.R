### example procedure to execute the package.

## 1st. load all functions (...)

## 2nd. read data (must change)

## x
#       x    pN    pS
#1  0.025 22490 17189
#2  0.075  3217  4780
#3  0.125  1616  2874
#4  0.175   999  2088
#5  0.225   754  1685
#6  0.275   679  1443
#7  0.325   575  1264
#8  0.375   484  1232
#9  0.425   427  1148
#10 0.475   437  1068
#11 0.525   378   986
#12 0.575   341   928
#13 0.625   310   893
#14 0.675   335   928
#15 0.725   315   945
#16 0.775   297   822
#17 0.825   326   885
#18 0.875   369   953
#19 0.925   448  1086
#20 0.975  1019  1904

## y
#   Chr Pop     m0f   D0f    m4f   D4f
#1 Chr2L RAL 2598805 54641 620019 52537


## 3rd. perform analysis
w <- iMK(x, y, 0, 1)
#check watchdog:
# w <- iMK(x, y, 0, 1.4)
#plot
plot1(w$iMK)
plot2(x, y, w$asymptoticMK, 0, 1)
