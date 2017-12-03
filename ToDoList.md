## iMKT  
### List of things to do:

- ~~Discuss input, maybe create parser~~
- Discuss output, update plots and include table -> html/pdf for non expert R users?  (Marta **on progres**)
- Plot functions (plotalpha, plotdaf, plotimkr) are not necessary, the plots are done inside the iMK function: remove them (Marta)
- Plot c: y axis text should be bold (Sergi, Marta -- in reality the "Asymptotic MKT" text shouldn't appear xD)
- Do you think it's useful to compute ω, Ka and Ks (internally, then the user can acces it)? It's would be easier to compute ωα and ωd.
- ~~When yo run multiple times asymptotic_MK, the CI are different. We should add a set.seed parameter in case the user wants to make it reproducible. If the argument is empty, then set.seed is NULL. (Marta) ~~
- ~~Update watchdog: adapt to diverse functions~~
- ~~Delete watchdog and add stops and messages inside the functions (Jesus) (selected: i, neutral: 0)~~
- ~~Error handling when function fails (Jesus):~~
	- ~~Currently checked: data.frame size(NCOL(x)==3;NROW(x)==9,19[daf10/20 categories] & NROW(x)!=0;NCOL(y)==6; NROW(y)==3 & NCOL(y)!=0), global columns names (daf,Pi,P0/Chr,pop,mi,Di,m0,D0), columns values (boundaries,types,NaN,0 [need to check possible values in Pi and P0];done)~~
	- ~~Watchdog only in assymtoticMKT(It is Needed in all functions? If just in assymptotic, include it inside the function)~~
	- ~~Check if some dataset not execute the asymototic (p0<=0). Include inside asysmtotic( check_input(...)) parameter to stop this functions but continue on loop (require for multiple_dataset(...) and iMK (...) (Jesus)~~

- System.time(Re-run iMKT):  
	- ~~First	approach, three datasets from concatenate.genes in 1000 iteration loop 100seconds aprox (Jesus)~~ 
	- Check times between genes and concatenate.genes dataset. Same expected time ~ 
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
		- Change RD files to S4 class (bioconductor) (Jesús)
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
		- Just need a folder to check the URLs and the format output (Jesus **in progess**) 
			- http://popfly.uab.cat/files/genes/GenesData_ALL.tab
			- http://pophuman.uab.cat/files/genes/GenesData_ALL.tab
		- Need talk about methods: (1)download and process whole dataset in R, (2)preload files in package (function to process them), (3)daf+div files in folder, call each one from R function.
		- Best of three times to load http://popfly.uab.cat/files/genes/GenesData_ALL.tab in a data.frame using read.table(URL) = 3.391 seconds  
		- Best of three times to load http://pophuman.uab.cat/files/genes/GenesData_ALL.tab in a data.frame using read.table(URL) = 94.418 seconds
- Reference Messer & Haller code (Question: Shall we write to them to let them know we're implementing their code into another package?)


- Update sample data (Marta, Jesus)
- Your input files x has the following(s) errors: x$P0 has one o more 0 values, cannot compute asymptotic. Returning results with NaN
- Example tutorial with sample data (Marta, Jesus)

- Implement GUI through web-server (with Django) (Ask Esteve for help?)  

### Bioconductor code requeriments (Jesus). 
Package is already functionally but we need to solve some bugs, finish a couple of extra functions and rewrite the functions in order to pass the bioconductor requiriments
#### Build aspects
- Must pass R CMD build (or R CMD INSTALL --build) and pass R CMD check with no errors and no warnings using a recent R-devel (**We have warning due to the high number of packages used**)
- Packages must also pass R CMD BiocCheck with no errors and no warnings. The BiocCheck package is a set of tests that encompass Bioconductor Best Practices (**not tested yet**)
- ~~Do not use filenames that differ only in case, as not all file systems are case sensitive~~
- ~~The source package resulting from running R CMD build should occupy less than 4MB on disk~~
- ~~The package should require less than 5 minutes to run R CMD check --no-build-vignettes~~
- Using the --no-build-vignettes option ensures that the vignette is built only once
- Vignette and man page examples should not use more than 3GB of memory since R cannot allocate more than this on 32-bit Windows
- ~~Choose a descriptive name. An easy way to check whether your name is already in use is to check that the following command fails (checked: no mkt or imkt)~~
- Core Bioconductor packages are typically licensed under Artistic-2.0. To specify a non-standard license, include a file named LICENSE in your package (containing the full terms of your license) and use the string “file LICENSE” (without the double quotes) in the “License:” field of your DESCRIPTION file.
- Package must:
	- Contain a vignette that demonstrates how to use the package to accomplish a task (more on this below)
	- Include examples in all man pages
	- Specify one or more biocViews categories ¿?
	- ~~Contain a NAMESPACE file to define the functions, classes, and methods that are imported into the name space, and exported for users~~
	- Document data structures used and, if different from data structures used by similar packages, explain why a different data structure was used
	- ~~Contain only code that can be redistributed according to the package license. Be aware of the licensing agreements for packages you are depending on in your package. Messer & Petrov GLP-3~~
	- ~~Not contain unnecessary files such as .DS_Store, .project, .svn, cache file, log files, etc.~~.
- ~~Packages you depend on must be available via Bioconductor or CRAN~~
- ~~Reuse, rather than re-implement or duplicate, well-tested functionality from other packages. Specify package dependencies in the DESCRIPTION file, listed as follows~~

#### Coding style
- S4 Classes and Methods (don't understand the bioconductor description)
- Vectorize and Pre-allocate and fill’ if iterations are necessary
- Check 1:n or 1:length(x) in loops. Use seq_len(n) or seq_along(x)
- Use 4 spaces for indenting. No tabs (SERIOUSLY BIOCONDUCTOR?!)
- No lines longer than 80 characters (check ggplot lines)
- Use camelCaps: initial lowercase, then alternate case between words
- Do not use ‘.’
- Filename extension for R code should be ‘.R’. Use the prefix ‘methods-‘ for S4 class methods, e.g., ‘methods-coverage.R’. Generic definitions can be listed in a single file, ‘AllGenerics.R’, and class definitions in ‘AllClasses.R’.
- Filename extension for man pages should be ‘.Rd’.
- Always use space after a comma. This: a, b, c.
- No space around “=” when using named arguments to functions. This: somefunc(a=1, b=2)
- Space around all binary operators: a == b.
- Use “##” to start full-line comments.
- Indent at the same level as surrounding code.
- Use <- not = for assignment.
- Import all symbols used from packages other than “base”. Except for default packages (base, graphics, stats, etc.) or when overly tedious, fully enumerate imports.
- Export all symbols useful to end users. Fully enumerate exports.
- Use dev.new() to start a graphics device if necessary. Avoid using x11() or X11() for it can only be called on machines that have access to an X server ¿?

####Package Author and Maintainer Responsibilities
Acceptance of packages into Bioconductor brings with it ongoing responsibility for package maintenance. These responsibilities include:
- Subscription to the bioc-devel mailing list.
- Registration on the support site.
- Response to bug reports and questions from users regarding your package, as posted on the Bioconductor support site or directly to developers. Add a BugReports: field to the DESCRIPTION file if reports should be directed to a particular web page rather than the package maintainer. You should register on the support site and edit your profile, changing the “Watched Tags” field to include all packages you maintain, so you will be notified when anyone posts a question about your package.
- Package maintenance through software release cycles, including prompt updates to software and documentation necessitated by underlying changes in R.
- All authors mentioned in the package DESCRIPTION file are entitled to modify package source code. Changes to package authorship require consent of all authors.

### Beta Tests
- Check inputs works with check_input(...) (Marta, Jesus)
- Check functions independently. ~~Error in DGRP (Marta, Jesus)~~
- Review tryCatch({...}) as proper error (Jesus)
