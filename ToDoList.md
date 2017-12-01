## iMKT  
### List of things to do:

- ~~Discuss input, maybe create parser~~
- Discuss output, update plots and include table -> html/pdf for non expert R users?  
- ~~Update watchdog: adapt to diverse functions~~
- ~~Delete watchdog and add stops and messages inside the functions (Jesus) (selected: i, neutral: 0)~~
- ~~Error handling when function fails (Jesus):~~
	- ~~Currently checked: data.frame size(NCOL(x)==3;NROW(x)==9,19[daf10/20 categories] & NROW(x)!=0;NCOL(y)==6; NROW(y)==3 & NCOL(y)!=0), global columns names (daf,Pi,P0/Chr,pop,mi,Di,m0,D0), columns values (boundaries,types,NaN,0 [need to check possible values in Pi and P0];done)~~
	- ~~Watchdog only in assymtoticMKT(It is Needed in all functions? If just in assymptotic, include it inside the function)~~
	- ~~Check if some dataset not execute the asymototic (p0<=0). Include inside asysmtotic( check_input(...)) parameter to stop this functions but continue on loop (require for multiple_dataset(...) and iMK (...) (Jesus)~~

- System.time(Re-run iMKT):  
	- ~~First	approach, three datasets from concatenate.genes in 1000 iteration loop ~ 100" (Jesus) ~~ 
	- Check times between genes and concatenate.genes dataset. Same expected time ~ 
- Build package:
	- ~~Main structure (Jesus)~~
	- ~~Review license (Need to talk, probably GPL-3)~~
	- Extra stuff about Roxygen2 and the .Rd format 
	- See section 'Suggested packages' in the 'Writing R Extensions' manual.  requireNamespace(). Packages: knitr, MASS, ggplot2, ggthemes, grid, nls2, extragrid, cowplot (Jesus **in progress**)
	- Github, CRAN or bioconductor 
	- Tests and hidden datasets (Jesus **to do**)

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
	- Modify perl functions to extract the correct files. Review the categories.(Sergi)
	- Put in functions comparision scripts (Jesus)
	- ~~Multiple_datasets(...): (Jesus)~~
		- Just need to add a variable to check file list (Jesus) (Jesus **to do**)
	- ~~Assymptotic: done (check variable and tables names)~~
	- iMK (Marta **in progress**)

- Reference Messer & Haller code (Question: Shall we write to them to let them know we're implementing their code into another package?)

- ~~Documentation of all functions (Marta, Jesus)~~
	- Main structure done, need to rewrite correctly.
- User manual (Marta, Jesus)
- Vignette (Marta, Jesus)
- Update sample data (Marta, Jesus)
- Example tutorial with sample data (Marta, Jesus)

- Implement GUI through web-server (with Django) (Ask Esteve for help?)  

## Beta Tests
- Check inputs works with check_input(...) (Marta, Jesus)
- Check functions independently. ~~Error in DGRP (Marta, Jesus)~~
- Review tryCatch({...}) as proper error (Jesus)
