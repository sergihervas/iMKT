[![Build Status](https://travis-ci.com/sergihervas/iMKT.svg?token=zpMDJ1ixtfzEon7BwxwR&branch=master)](https://travis-ci.com/sergihervas/iMKT)

# iMKT: integrative McDonald and Kreitman Test


Overview
--------
iMKT is an R package to compute the integrative McDonald and Kreitman test. It includes several MK derived methodologies which allow inferring the rate of adaptive evolution (α) and two new approximations to quantify the different negative selection regimes (d: strongly deleterious, b: weakly deleterious, f: neutral) from polymorphism and divergence genomic data.


Installation
------------
The package is deposited in GitHub and must be installed using the devtools library.
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
Shortly, iMKT allows performing diverse MK-derived tests using the number of polymorphic (P, classified in DAF categories), divergent (D) and analyzed (m) sites for neutral (0) and selected (i) classes. In brief, most functions require two input parameters: ```daf``` (data frame containing daf, Pi and P0) and ```divergence``` (data frame containing mi, Di, m0, D0) and return the estimation of α together with specific details of the methodology.

The package includes two sample data frames (```myDafData```, ```myDivergenceData```). The vignettes and manual documentation contains detailed descriptions and examples regarding each function, types of analyses, how to use PopFly and PopHuman genomic data, etc.

The following example shows how to perform standard MKT using sample data:
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



---------------------------------------------------------------------
---------------------------------------------------------------------

## List of things to do:

### R Package
- Update Readme.md
- Update manual.pdf: "R CMD Rd2pdf /home/sergi/git/iMKT/"
- Build package:
	- Extra stuff about Roxygen2 and the .Rd format
	- Document correctly functions
	- Vignettes 
	- Check functions independently / Tests
	- Rename variables (wordWord2). **Do it for all variables inside all functions or only for functionNames and dataFrames?**
	- Check consistency of variables and style along functions
- Functions:
	- Put in functions comparision scripts (Jesus) **What is this?**
	- PopHuman functions: **Sergi**
		- **S** Check times of download functions: loadPopFly() & loadPopHuman. Around 20 seconds the first one	
		- Update pophuman data file!		
- Reference Messer & Haller code (Question: Shall we write to them to let them know we're implementing their code into another package?). **S** We reference their code in the asymptoticMK function. We can write them later, as Antonio suggested.


### Server
- Implement GUI through web-server (with Django) (Ask Esteve for help? **M** the web looks awesome Jesus!!!)(Jesús **in progress**)  
- Rewrite funcionts to run from terminal (receiving input/inputs) and generate html (Jesús **in progress**)
- daf.pl: Modify perl functions to extract the correct files.


### Manuscript
- **Follow the manuscript progress at: https://github.com/marta-coronado/statistics_mkt/tree/master/Report**
- Prepare list of genes analyzable with the iMKT (not the same as asymptotic because the exponential model is forced)
- **S** It would also be great to know the minimum number of sites required to run the different tests. This makes sense for asymptotic alpha basically (and the iMK function). Should be done with simulated data, but is a tough work. **M** This was also commented on one of the meetengs and that I should definitely take a look on this, also is super important and nobody tested it yet. And not only the minimum number, but also the minimum number of concatenated genes for the asymptotic and DFE-alpha.

- Manuscript: **M** suggests: 

	- first, start writing the **Methods** (Sergi/Jesús: explain the data pipeline and such; Marta: I'm defining the different statistical test performed)
		- Drosophila genome data: input seqs (DGN), reference annotations and outgroup species (refer PopFly).
		- Human genome data: same as before (refer PopHuman). Maybe join 1 and 2 in Data description or sth like this.
		- Main pipeline / workflow: Fasta/VCF > Recodification > DAF / div 
		- Statistical tests (iMK negative selection is completely new and unpublished before, although it is based on DGRP idea. In fact what is new is the way to estimate b -weakly deleterious fraction. Also the DGRP method was never mentioned, only in the suplementary material where has been forgotten. So this time we should really promote them!)
		- Simulations.
		- R package
		- Web server
		- Things about candidate genes. Statistical tests used, GO...
	
	- **Results**: **M** I was thinking in organizing the following parts:
		- Comprison of methodologies: using Drosophila data, apply the different MKT tests (standard, DGRP, FWW, integrative). **S** And human data? I say so because of point 4. Or here we just assess the best method and apply this one in point 4? **M** Okay, at the moment is done with Drosophila because we still don't have the human data
		- DFE-Based Extensions of the MK: add a part devoted to DFE-alpha (**S** all methods are in here, right? DGRP, FWW and asymptotic are based on DFE assumptions which allow spliting mutations according to their fitness. **M** No, DFE-Based extensions is only DFE-alpha (and some more which I don't know), the rest no, they are extensions of MK test. The difference is that DFE-alpha try to estimate how many non-adaptive substitutions will become fix given an inferred DFE. The others, assume that the non-adaptive substitutions can be removed removing low-frequency variants (FWW) or a more sophisticated solution (DGRP) or extrapolating alpha using a function that fits data at different frequencyes (asymptotic), but nothing to do with inferring the DFE of mutations. Hope its clear :) **S** Okay! I thought they also assumed an specific DFE to remove non-adaptive substitutions, my bad.)
		- Simulations: comparison of the different methodologies against simulated data
		- Adaptation in the human and D. melanogaster genome: genes that have alfa positive and significative, and study of them (eg: GO, networks...)
		- The package and the webpage. The pipeline for obtaining the DAF? **S** This pipeline is a result or a method? I always doubt on this kind of things. However, it is something referees can criticize a lot because it is not perfect and I think the work is very complete and long enough, so I would not talk about it in the Results section, just in the Methods. **M** Okay, then we only comment it on the methods!
		- Something else?
