## iMKT  
### List of things to do:

- Update Readme.md
- ~~Discuss input, maybe create parser~~
- Discuss output, update plots and include table -> html/pdf for non expert R users?  (Marta **in progress**)
    - Discuss invisible(output)
- ~~Plot functions (plotalpha, plotdaf, plotimkr) are not necessary, the plots are done inside the iMK function: remove them (Marta)~~
- ~~Plot c: y axis text should be bold (Sergi, Marta -- in reality the "Asymptotic MKT" text shouldn't appear xD)~~
- Do you think it's useful to compute ω, Ka and Ks (internally, then the user can acces it)? It would be easier to compute ωα and ωd. [**S** yes we should include them. maybe inside the standard mkt function? **J** I agree with **S**]
- ~~When yo run multiple times asymptotic_MK, the CI are different. We should add a set.seed parameter in case the user wants to make it reproducible. If the argument is empty, then set.seed is NULL. (Marta)~~
- ~~Update watchdog: adapt to diverse functions~~
- ~~Delete watchdog and add stops and messages inside the functions (Jesus) (selected: i, neutral: 0)~~
- ~~Error handling when function fails (Jesus):~~
	- ~~Currently checked: data.frame size(NCOL(x)==3;NROW(x)==9,19[daf10/20 categories] & NROW(x)!=0;NCOL(y)==6; NROW(y)==3 & NCOL(y)!=0), global columns names (daf,Pi,P0/Chr,pop,mi,Di,m0,D0), columns values (boundaries,types,NaN,0 [need to check possible values in Pi and P0];done)~~
	- ~~Watchdog only in assymtoticMKT(It is Needed in all functions? If just in assymptotic, include it inside the function)~~
	- ~~Check if some dataset not execute the asymototic (p0<=0). Include inside asysmtotic( check_input(...)) parameter to stop this functions but continue on loop (require for multiple_dataset(...) and iMK (...) (Jesus)~~

- System.time(Re-run iMKT):  
	- ~~First approach, three datasets from concatenate.genes in 1000 iteration loop 100seconds aprox (Jesus)~~ 
	- Check times between genes and concatenate.genes dataset. Same expected time ~ 
	- Functions timing in concatenate dataset (daf10): 
		- standard =  0.082 seconds (best of three)
		- FWW (with plot) =   0.191 seconds (best of three)
		- DGRP (with plot) = 0.373 seconds (best of three)
		- asymptotic = 0.228 seconds (best of three)
		- iMK = 0.457 seconds (best of three)
		- integrativeMKT = 1.408 seconds(best of three)
	- Functions timing in genes dataset (daf10): 
		- standard =  0.003 seconds (best of three)
		- FWW (with plot) =   0.013 seconds (best of three)
		- DGRP (with plot) = 0.139 seconds (best of three)
		- asymptotic = 0.228 seconds (best of three)
		- iMK = 0.481  seconds (best of three)
		- integrativeMKT = 0.885 seconds (best of three)
- Build package:
	- ~~Main structure (Jesus)~~
	- ~~Review license (Need to talk, probably GPL-3)~~
	- Extra stuff about Roxygen2 and the .Rd format
	- ~~See section 'Suggested packages' in the 'Writing R Extensions' manual.  requireNamespace(). Packages: knitr, MASS, ggplot2, ggthemes, grid, nls2, extragrid, cowplot (Jesus)~~
	- Github, CRAN or bioconductor. Marta: I'd say the package should be kept in Github, at least privately, so we can still improve it, easy to handle etc. Between CRAN/Bioconductor I'd choose Bioconductor for the following reasons: (1) The package will be announced on Twitter by core Bioconductor members and listed on Bioconductor website --> a lot of impact in the community. (2) Professional suppport. (3) It's easy to publish it if it's already on Bioconductor. (4) Packages can be submitted in Bioconductor through Github. (5) We'll have a doi and it's citable. Disadvantatges: much tougher review process than CRAN, must have a vignette, the code should be 80 characters wide, we can't use S3 class (???), the package must be re-updated each 6 months. See: https://bioconductor.org/developers/how-to/efficient-code/, https://bioconductor.org/developers/how-to/coding-style/   
	- ~~Hidden datasets (Jesus)~~
	- Tests (Jesus)
	- ~~Document correctly functions (Marta, Jesus)~~
		- Just need to rewrite correctly.
		- ~~Change RD files to S4 class (bioconductor) (Jesús). FINALLY NOT BIOCONDUCTOR~~
	- User manual (Marta, Jesus)
	- Vignettes (Marta, Jesus)

- ~~Check consistency of variables and style along functions (Marta)~~
- Functions:
	- ~~Standard: done~~
	- ~~FWW:~~
		- ~~done~~  
		- ~~add loop~~  
		- ~~graph cutoffs simple (Marta)~~
	- ~~DGRP: to do, " " (Marta)~~
		- ~~done~~  
		- ~~add loop and graph cuttoff (Marta)~~
		- Asked Sònia about b, y - when they are negative is because you haven't completely removed the slightly deleterious variants. She sets the values to 0 then. Ask Antonio for a better estimation?
		- ~~graph cutoffs simple (Marta)~~
	- Modify perl functions to extract the correct files. Review the categories. (Sergi **in progress**)
	- Put in functions comparision scripts (Jesus)
	- ~~Multiple_datasets(...): (Jesus)~~
		- ~~Add a variable to check file list (Jesus)~~
		- Check Pre-allocate and fill’ (Jesus: Marta check the function pohfavó)
	- ~~Assymptotic: done (check variable and tables names)~~
	- iMK (Marta **in progress**)
	- ~~Retrieve information by gene or genelist from PopFly and PopHuman (Jesus)~~
		- ~~Just need a folder to check the URLs and the format output (Jesus)~~
			- http://popfly.uab.cat/files/genes/GenesData_ALL.tab
			- http://pophuman.uab.cat/files/genes/GenesData_ALL.tab
		- Need talk about methods: ~~(1)download and process whole dataset in R~~, (2)preload files in package (function to process them), ~~(3)daf+div files in folder, call each one from R function.~~
		- ~~Best of three times to load http://popfly.uab.cat/files/genes/GenesData_ALL.tab in a data.frame using read.table(URL) = 3.391 seconds~~
		- ~~Best of three times to load http://pophuman.uab.cat/files/genes/GenesData_ALL.tab in a data.frame using read.table(URL) = 94.418 seconds~~
		- ~~**S** check times using directories structure with andromeda IP.~~
		- Add recomb else if to extract genes by recombination bin and population (Sergi)
- ~~Reference Messer & Haller code (Question: Shall we write to them to let them know we're implementing their code into another package?)~~


- Update sample data (Marta, Jesus)
- Your input files x has the following(s) errors: x$P0 has one o more 0 values, cannot compute asymptotic. Returning results with NaN
- Example tutorial with sample data (Marta, Jesus)

- Implement GUI through web-server (with Django) (Ask Esteve for help?)  

### Beta Tests
- Check inputs works with check_input(...) (Marta, Jesus)
- Check functions independently. ~~Error in DGRP (Marta, Jesus)~~
- Review tryCatch({...}) as proper error (Jesus)
