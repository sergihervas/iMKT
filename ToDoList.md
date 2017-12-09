## iMKT  
### List of things to do:

- Update Readme.md
- ~~Discuss input, maybe create parser~~
- Discuss output, update plots and include table -> html/pdf for non expert R users?  (Marta **in progress**)
    - Discuss invisible(output)
- ~~Plot functions (plotalpha, plotdaf, plotimkr) are not necessary, the plots are done inside the iMK function: remove them (Marta)~~
- ~~Plot c: y axis text should be bold (Sergi, Marta -- in reality the "Asymptotic MKT" text shouldn't appear xD)~~
- Do you think it's useful to compute ω, Ka and Ks (internally, then the user can acces it)? It would be easier to compute ωα and ωd. **S** yes we should include them. maybe inside the standard mkt function? **J** I agree with **S**. **M** Yes, but those values should be available for all tests performed: how about: function "divergence_parameters.R" which run inside all tests and the output is invisible? **S** I do not understand why these parameters should be in all funcions. Aren't they like other statistics to estimate the fraction of adaptive evolution, such as alpha? We can discuss this on Tuesday, I guess I am missing something here. I would suggest one single independent function to compute w, kaks, wa and wd.
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
	- **S** I do not understand this comparison. Differences in time depend on the values of Pi, P0, Di and D0, but not on the type of data used (1 gene or multiple genes) if the input is the same (daf and div files). Of course the multiple genes test takes more time, but because P and D values are higher. Maybe it would make more sense to use "simulated" (I mean inventados xd) datasets controlling Pi, P0, Di, D0 or different number of daf categories to perform the speed tests.
	- **S** It would also be great to know the minimum number of sites required to run the different tests. This makes sense for asymptotic alpha basically (and the iMK function). Should be done with simulated data, but is a tough work.
	- Functions timing in concatenate dataset (daf10): 
		- standard =  0.082 seconds (best of three)
		- FWW (with plot) =   0.191 seconds (best of three)
		- DGRP (with plot) = 0.373 seconds (best of three)
		- asymptotic = 0.228 seconds (best of three)
		- iMK = 0.457 seconds (best of three)
		- completeMKT = 1.408 seconds(best of three) **M**: why it takes so long if it should be like iMK + DGRP??? **S** I think this function computes all the above tests calling each function one by one (0.1+0.2+0.4+0.2+0.5=1.4)
	- Functions timing in genes dataset (daf10): 
		- standard =  0.003 seconds (best of three)
		- FWW (with plot) =   0.013 seconds (best of three)
		- DGRP (with plot) = 0.139 seconds (best of three)
		- asymptotic = 0.228 seconds (best of three)
		- iMK = 0.481  seconds (best of three)
		- completeMKT = 0.885 seconds (best of three)
- Build package:
	- ~~Main structure (Jesus)~~
	- ~~Review license (Need to talk, probably GPL-3)~~
	- Extra stuff about Roxygen2 and the .Rd format
	- ~~See section 'Suggested packages' in the 'Writing R Extensions' manual.  requireNamespace(). Packages: knitr, MASS, ggplot2, ggthemes, grid, nls2, extragrid, cowplot (Jesus)~~
	- ~~Github, CRAN or bioconductor. Marta: I'd say the package should be kept in Github, at least privately, so we can still improve it, easy to handle etc. Between CRAN/Bioconductor I'd choose Bioconductor for the following reasons: (1) The package will be announced on Twitter by core Bioconductor members and listed on Bioconductor website --> a lot of impact in the community. (2) Professional suppport. (3) It's easy to publish it if it's already on Bioconductor. (4) Packages can be submitted in Bioconductor through Github. (5) We'll have a doi and it's citable. Disadvantatges: much tougher review process than CRAN, must have a vignette, the code should be 80 characters wide, we can't use S3 class (???), the package must be re-updated each 6 months. See: https://bioconductor.org/developers/how-to/efficient-code/, https://bioconductor.org/developers/how-to/coding-style/~~   
	- ~~Hidden datasets (Jesus)~~
	- Tests (Jesus, Sergi)
	- ~~Document correctly functions (Marta, Jesus)~~
		- Just need to rewrite correctly the documentation.
		- ~~Change RD files to S4 class (bioconductor) (Jesús). FINALLY NOT BIOCONDUCTOR~~
	- User manual (Marta, Jesus)
	- Vignettes (Marta, Jesus)

- Check consistency of variables and style along functions (Marta, Sergi, Jesus) **Should be also done at the end of coding**
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
		- **S** Test DGRP estimates of f,b,d with real and simulated data. I do not like them, in fact I think the formulae are wrong, and that's why I did not include them in the original DGRP function, just the alpha value.
		- ~~graph cutoffs simple (Marta)~~
	- Modify perl functions to extract the correct files. Review the categories. (Sergi **in progress**) **S** What do u mean with: Review the categories?
	- Put in functions comparision scripts (Jesus) **S** What is this?
	- ~~Multiple_datasets(...): (Jesus)~~
		- ~~Add a variable to check file list (Jesus)~~
		- Check Pre-allocate and fill’ (Jesus: Marta check the function pohfavó. **M**: checking in progress) **S** If this has something to do with for loops let me know cause I use lots of them. As far as I understood it is about emtpy lists which are then filled and returned in the function, rigth?
	- ~~Assymptotic: done (check variable and tables names)~~
	- ~~iMK (Marta)~~
	- subsetPopData() **Sergi**
		- ~~Need talk about methods:(1)download and process whole dataset in R, (2)preload files in package (function to process them), (3)daf+div files in folder, call each one from R function.~~
		- **S** Check times of download functions: loadPopFly() & loadPopHuman. Around 20 seconds the first one. 
		- Discuss output of the function. **S** I would suggest (and I will implement it soon if you agree) to allow deciding which test to perform with the subseted data from PopFly or PopHuman, so the function would call compareMK function and give its output in lists for populations.

- Reference Messer & Haller code (Question: Shall we write to them to let them know we're implementing their code into another package?)		

- Update sample data (Marta, Jesus)
- Example tutorial with sample data (Marta, Jesus)

- Implement GUI through web-server (with Django) (Ask Esteve for help? **M** the web looks awesome Jesus!!!)(Jesús **in progress**)  
	- Rewrite funcionts to run from termninal (receiving input/inptus) and generate html (Jesús **in progress**) **S** U mean using curl? Like the asymptoticMK webpage? For the server we can use FastR which allows running R online and generating reports in markdown-like style. That's what we use in PopFly / PopHuman for MKT gene report.

- Manuscript: **M** suggests: 
	-- first, start writing the **Methods** (Sergi/Jesús: explain the data pipeline and such; Marta: I'm defining the different statistical test performed)
		1. Drosophila genome data: input seqs (DGN), reference annotations and outgroup species (refer PopFly).
		2. Human genome data: same as before (refer PopHuman). Maybe join 1 and 2 in Data description or sth like this.
		3. Main pipeline / workflow
			- Fasta/VCF > Recodification > DAF / div 
		4. Statistical tests (iMK negative selection is completely new and unpublished before, although it is based on DGRP idea. In fact what is new is the way to estimate b -weakly deleterious fraction.)
		5. Simulations
		6. R package
		7. Web server
		8. Things about Results 4. Statistical tests used
	
	-- **Results**: **M** I was thinking in organizing the following parts:
		1. Comprison of methodologies: using Drosophila data, apply the different MKT tests (standard, DGRP, FWW, integrative). **S** And human data? I say so because of point 4. Or here we just assess the best method and apply this one in point 4?
		2. DFE-Based Extensions of the MK: add a part devoted to DFE-alpha (**S** all methods are in here, right? DGRP, FWW and asymptotic are based on DFE assumptions which allow spliting mutations according to their fitness)
		3. Simulations: comparison of the different methodologies against simulated data
		4. Adaptation in the human and D. melanogaster genome: genes that have alfa positive and significative, and study of them (eg: GO, networks...)
		5. The package and the webpage. The pipeline for obtaining the DAF? **S** This pipeline is a result or a method? I always doubt on this kind of things. However, it is something referees can criticize a lot because it is not perfect and I think the work is very complete and long enough, so I would not talk about it in the Results section, just in the Methods.
		6. Something else?
		**S** I would suggest permuting points 4 and 5. Hence, we first present the package and server and then the adaptation results which we obtained using the previously described software. This way we demonstrate it is useful.
		
### Beta Tests
- Check inputs works with check_input(...) (Marta, Jesus)
- Check functions independently. ~~Error in DGRP (Marta, Jesus)~~
- ~~Review tryCatch({...}) as proper error (Jesus)~~
- Error handling predictNLS (Sergi)
- check()
	- NOTE: import / export. too many packages...
