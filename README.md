[![Build Status](https://travis-ci.org/sergihervas/iMKT.svg?branch=master)](https://travis-ci.org/sergihervas/iMKT)

# iMKT: integrative McDonald and Kreitman Test


Overview
--------
iMKT is an R package to compute the McDonald and Kreitman test (McDonald and Kreitman, 1991, Nature) on polymorphism and divergence genomic data provided by the user or automatically downloaded from PopFly (Hervas et al. 2017 Bioinformatics) or PopHuman (Casillas et al., 2017 Nucleic Acids Res.). It includes five MK derived methodologies which allow inferring the rate of adaptive evolution (α) as well as the fraction of strongly deleterious (d), weakly deleterious (b), and neutral (f) sites.


Installation
------------
The package is deposited in GitHub and can be installed using the devtools library.
``` r
## Install devtools package if necessary
install.packages("devtools")

## Install iMKT package from GitHub
devtools::install_github("sergihervas/iMKT")

## Load iMKT library
library(iMKT)
```


Usage
-----
In summary, iMKT allows performing diverse MK-derived tests using the number of polymorphic (P, classified in Derived Allel Frequency (DAF) categories), divergent (D) and analyzed (m) sites for neutral (0) and selected (i) classes. Briefly, most functions require two input parameters: ```daf``` (data frame containing DAF, Pi and P0) and ```divergence``` (data frame containing mi, Di, m0, D0) and return the estimation of α together with specific details of the methodology.

The package includes two sample data frames (```myDafData```, ```myDivergenceData```). The vignettes and manual documentation contain detailed descriptions and examples of each function and type of analysis, and instructions on how to use PopFly and PopHuman genomic data.

The following example, which shows how to perform a standard MK test using sample data, is aimet at illustrating how the iMKT package works.
``` r
## Sample daf data included in the package
head(myDafData)
#>     daf    Pi    P0
#> 1 0.025 22490 17189
#> 2 0.075  3217  4780
#> 3 0.125  1616  2874
#> 4 0.175   999  2088
#> 5 0.225   754  1685
#> 6 0.275   679  1443

## Sample divergence data included in the package
myDivergenceData
#>        mi    Di     m0    D0
#> 1 2598805 54641 620019 52537

## Perform standard MKT
standardMK(myDafData, myDivergenceData)
#> $alpha.symbol
#> [1] 0.2364499

#> $`Fishers exact test P-value`
#> [1] 1.480943e-183

#> $`MKT table`
#> |               | Polymorphism| Divergence|
#> |:--------------|------------:|----------:|
#> |Neutral class  |        45101|      52537|
#> |Selected class |        35816|      54641|

#> $`Divergence metrics`
#> |        Ka|        Ks|     omega|   omegaA|   omegaD|
#> |---------:|---------:|---------:|--------:|--------:|
#> | 0.0210254| 0.0847345| 0.2481331| 0.058671| 0.189462|
```


Citation
--------
Citation to paper


Licence
-------
Licence of package


Development & Contact
---------------------
iMKT has been developed by Sergi Hervas (sergi.hervas@uab.cat), Marta Coronado (marta.coronado@uab.cat) and Jesús Murga (jesus.murga@uab.cat), from the Bioinformatics of Genome Diversity group from the Universitat Autònoma de Barcelona (UAB) and the Institut de Biotecnologia i Biomedicina (IBB).

If you have any feedback or feature requests regarding iMKT, please contact antonio.barbadilla@uab.cat or jesus.murga@uab.cat.