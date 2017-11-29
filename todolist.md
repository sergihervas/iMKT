##iMKT ToDoList


 ̶-̶ ̶D̶i̶s̶c̶u̶s̶s̶ ̶i̶n̶p̶u̶t̶,̶ ̶m̶a̶y̶b̶e̶ ̶c̶r̶e̶a̶t̶e̶ ̶p̶a̶r̶s̶e̶r̶
 
- Discuss output, update plots and include table -> html/pdf for non expert R users?

- Update watchdog: adapt to diverse functions
-- Delete watchdog and add stops and messages inside the functions (Jesus) (selected: i, neutral: 0)
- Error handling when function fails (in progress Jesus).
-- Currently checked: data.frame size(NCOL(x)==3;NROW(x)==9,19[daf10/20 categories] & NROW(x)!=0;NCOL(y)==6; NROW(y)==3 & NCOL(y)!=0), global columns names (daf,Pi,P0/Chr,pop,mi,Di,m0,D0), columns values (boundaries,types,NaN,0[need to check possible values in Pi and P0],)
-- Watchdog only in assymtoticMKT(It is Needed in all functions? If just in assymptotic, include it inside the function)
- System.time(Re-run iMKT) and compare(gene.input/concatgenes.input) (Jesus)
- Build package (connect to other packages, CRAN, tests and hidden datasets) (Jesus)

- Rename certain variables in asymptoticMK and iMK (0f, 4f) (Marta)
- Check consistency of variables and style along functions (Marta)
- Functions:
-- [x] Standard: done
-- FWW: done, but add loop, graph (in progress Marta)
-- DGRP: to do, " " (Marta)
-- Multiple datasets (Jesus)
-- Assymptotic: done (revisar)

- Reference Messer & Haller code

- Documentation of all functions
- User manual
- Vignette 
- Update sample data

- SERGI: Li function (+ include in iMK); StandardMKT function (+ include in iMK)
-- standardMKT done, Li in progress (Marta)

- Implement GUI through web-server (Django)