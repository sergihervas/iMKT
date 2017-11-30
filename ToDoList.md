## iMKT  
### List of things to do:

- ~~Discuss input, maybe create parser~~
- Discuss output, update plots and include table -> html/pdf for non expert R users?  
- ~~Update watchdog: adapt to diverse functions~~
- ~~Delete watchdog and add stops and messages inside the functions (Jesus) (selected: i, neutral: 0)~~
- ~~Error handling when function fails (Jesus):~~
	- ~~Currently checked: data.frame size(NCOL(x)==3;NROW(x)==9,19[daf10/20 categories] & NROW(x)!=0;NCOL(y)==6; NROW(y)==3 & NCOL(y)!=0), global columns names (daf,Pi,P0/Chr,pop,mi,Di,m0,D0), columns values (boundaries,types,NaN,0 [need to check possible values in Pi and P0];done)~~
	- ~~Watchdog only in assymtoticMKT(It is Needed in all functions? If just in assymptotic, include it inside the function)~~
	- Check if some dataset not execute the asymototic (p0<=0). Include inside asysmtotic( check_input(...)) parameter to stop this functions but continue on loop (require for multiple_dataset(...) and iMK (...) (Jesus **in progress**)

- System.time(Re-run iMKT) and compare (gene.input/concatgenes.input) (Jesus **in progress**;)  
- Build package (connect to other packages (require knitr, ggplot2, nls2), CRAN or bioconductor, tests and hidden datasets) (Jesus)

- ~~Rename certain variables in asymptoticMK and iMK (0f, 4f) (Marta)~~
- ~~Check consistency of variables and style along functions (Marta)~~
- Functions:
	- ~~Standard: done~~
	- ~~FWW:~~
		- ~~done~~  
		- ~~add loop~~  
		- ~~graph cutoffs simple (Marta)~~
	- ~~DGRP: to do, " " (Marta)~~
	- Put in functions comparision scripts (Jesus)
	- Multiple_datasets(...): (Jesus)
	- Assymptotic: done (check variable and tables names)

- Reference Messer & Haller code

- Documentation of all functions
- User manual
- Vignette 
- Update sample data

- Implement GUI through web-server (Django)  

## Beta Tests
- Check inputs works with check_input(...) (Marta, Jesus)
- Check functions independently (Marta, Jesus)
