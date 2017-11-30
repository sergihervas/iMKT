## ERROR HANDLING
## function which checks input data for asymptoticMK
## should check and discuss many things. Use as an example for now! NEED TO DISSCUSS predictNLS.R: Error in as.list(object$call$formula)[[3]] : subscript out of bounds 
## mix between Messer & Haller
# Author = Sergi Hervás, Jesús Murga

check_input <- function(x, y, xlow, xhigh){
  
  data_is_good <- TRUE
  main_errors<-'Your input files has the following(s) errors: '
  ######### Parsing NROW and NCOL ############
  #error handling: check data formatting
  if (NCOL(x) != 3 || NCOL(x) <3){
    data_is_good<-FALSE
    error<-"Argument x does not contain three tab-separated columns.Argument x, at least three data rows are required to constrain the fit"
    main_errors<-append(main_errors,error)}
  if (NROW(x) < 0){
    data_is_good<-FALSE
    error<-"Argument x contains no data rows or x categories are not correct define. Check the header!"
    main_errors<-append(main_errors,error)}
  if(NROW(x)!=20 & NROW(x)!=10){
    data_is_good<-FALSE
    error<-"Daf categories in x are not correct define. Check the header!"
    main_errors<-append(main_errors,error)}
  if (NCOL(y) != 6){
    data_is_good<-FALSE
    error<-"Argument y does not contain six tab-separated columns"
    main_errors<-append(main_errors,error)}
  if (NROW(y) <= 0 & (NROW(y)<3)){
    data_is_good<-FALSE
    error<-"Argument y contains no data rows"
    main_errors<-append(main_errors,error)}
  
  ##Assign proper names to the columns of x and y
  names(x) <- c("daf", "Pi", "P0")
  names(y) <- c("Chr",  "Pop",	"mi",	"Di",	"m0",	"D0")
  
  #error handling: input colnames 
  #for arg x
  
  # NO LO ENTIENDO,
  # x.colnames <- names(x)
  # suppressWarnings(#the goal is to generate NAs here, so we dont need no warnings
  #   if (!is.na(as.numeric(x.colnames[1])) ||
  #         !is.na(as.numeric(x.colnames[2])) ||
  #         !is.na(as.numeric(x.colnames[3])))
  #     error<-"argument x has a numeric column name;
  #          probably the required header row is missing")
  # 
  # #for arg y
  # y.colnames <- names(y)
  # suppressWarnings(
  #   if (!is.na(as.numeric(y.colnames[1])) ||
  #         !is.na(as.numeric(y.colnames[2])) ||
  #         !is.na(as.numeric(y.colnames[3])) ||
  #         !is.na(as.numeric(y.colnames[4])) ||
  #         !is.na(as.numeric(y.colnames[5])) ||
  #         !is.na(as.numeric(y.colnames[6])))
  #     error<-"argument y has a numeric column name;
  #          probably the required header row is missing"
  #     )
  
  
  ##parse the data from argument x
  ####### dependent on the column name
  f <- x$daf #derived alelle frequencies
  Pi <- x$Pi #non-synonymous polymorphism 
  P0 <- x$P0 #synonymous polymorphism
  
  #error handling: check if variables are good
  if (any(is.na(f))){
    data_is_good<-FALSE
    error<-"f contains NA values (not allowed)"
    main_errors<-append(main_errors,error)}
  if (any(is.na(Pi))){
    data_is_good<-FALSE
    error<-"pi contains NA values (not allowed)"
    main_errors<-append(main_errors,error)}
  if (any(is.na(P0))){
    data_is_good<-FALSE
    error<-"p0 contains NA values (not allowed)"
    main_errors<-append(main_errors,error)}
  #error handling: check if variables are numeric
  if (!is.numeric(f)){
    data_is_good<-FALSE
    error<-"f is not numeric"
    main_errors<-append(main_errors,error)}
  if (!is.numeric(Pi)){
    data_is_good<-FALSE
    error<-"pi is not numeric"
    main_errors<-append(main_errors,error)}
  if (!is.numeric(P0)){
    data_is_good<-FALSE
    error<-"p0 is not numeric"
    main_errors<-append(main_errors,error)}
  ##JUST CHECK IF THE VARIABLE IS NUMERIC, IF NOT, THEN stop. NOT NEED NULL. JUST stop IT IF AREN'T NUMERIC
  # if (is.null(f))
  #   error<-"f malformed (must be numeric)"
  # if (is.null(pi))
  #   error<-"p malformed (must be numeric)"
  # if (is.null(p0))
  #   error<-"p0 malformed (must be numeric)"
  # 
  # #error handling: check if variables are numeric
  # if (!is.numeric(f))
  #   error<-"f is not numeric"
  # if (!is.numeric(Pi))
  #   error<-"p is not numeric"
  # if (!is.numeric(P0))
  #   error<-"p0 is not numeric"
  
  #error handling: check if variables are not out of bounds.
  if (any(f < 0.0) || any(f > 1.0)){
    data_is_good<-FALSE
    error<-"f contains values out of the required range [0,1]"
    main_errors<-append(main_errors,error)}
  if (all(f == 0)){
    data_is_good<-FALSE
    error<-"f contains all values == 0 (not allowed)"
    main_errors<-append(main_errors,error)}
  if (any(Pi <= 0)||all(Pi == 0)){   # note that zero is allowed, although not recommended
    data_is_good<-FALSE
    error<-"Pi contains values < 0 (not allowed) Pi contains all values == 0 (not allowed)"
    main_errors<-append(main_errors,error)}
  # if (all(Pi == 0))    # not all can be zero, however
  #   error<-"p contains all values == 0 (not allowed)"
  if (any(P0 < 0)||all(P0 == 0)){
    data_is_good<-FALSE
    error<-"P0 contains values <= 0 (not allowed) or P0 contains all values == 0 (not allowed)"
    main_errors<-append(main_errors,error)}
  # if (all(p0 == 0))
  #   error<-"p0 contains all values == 0 (not allowed)"
  
  
  #Incorpore at begging when check the input length
  # error handling: check if argument x has enough data points
  if (NROW(x) < 3)
  error<-"Argument x: at least three data rows are required to constrain the fit"

  ##parse the data from argument y and force it to be numeric... NOT SENSE CHECK THEN IS.NULL OR NUMERIC
  # NOT FORCE TO BE NUMERIC, CHECK IF VALUES ARE NOT NUMERIC [not force+not is.null]
  mi <- as.numeric(y$mi) #number of non-synonymous sites
  m0 <- as.numeric(y$m0) ##number of synonymous sites
  Di <- as.numeric(y$Di) #non-synonymous divergence
  D0 <- as.numeric(y$D0) #synonymous divergence
  
  #error handling: check if variables are good
  if (is.na(mi) || !is.numeric(mi)){
    data_is_good<-FALSE
    error<-"malformed m (must be numeric)"
    data_is_good<-FALSE
    main_errors<-append(main_errors,error)}
  if (is.na(m0) || !is.numeric(m0)){
    error<-"malformed m0 (must be numeric)"
    data_is_good<-FALSE
    main_errors<-append(main_errors,error)}
  if (is.na(D0) || !is.numeric(D0)){
    data_is_good<-FALSE
    error<-"malformed D0 (must be numeric)"
    main_errors<-append(main_errors,error)}
  if (is.na(Di) || !is.numeric(Di)){
    data_is_good<-FALSE
    error<-"malformed d (must be numeric)"
    main_errors<-append(main_errors,error)}
  if (is.na(xlow) || is.null(xlow)){
    data_is_good<-FALSE
    error<-"malformed xlow (must be numeric)"
    main_errors<-append(main_errors,error)}
  if (is.na(xhigh) || is.null(xhigh)){
    data_is_good<-FALSE
    error<-"malformed xhigh (must be numeric)"
    main_errors<-append(main_errors,error)}
  # Checks if number of sites (m, m0) is not higher than divergenge (d, D0), the sum of the polimorphisms (p|p0) or the divergenge + the sum of the polimorphsms
  if (Di>mi || sum(Pi)>mi || sum(Pi)+Di>mi){
    data_is_good<-FALSE
    error<-"mi must be higher than Pi, Di and sum(Pi+Di)"
    main_errors<-append(main_errors,error)}
  if (D0>m0 || sum(P0)>m0 || sum(P0)+D0>m0){
    data_is_good<-FALSE
    error<-"m0 must be higher than P0, D0 and sum(P0+D0)"
    main_errors<-append(main_errors,error)}
  # #error handling: check if variables are numeric
  # #it makes no sence as we forced them to be numeric... should check that!
  # if (!is.numeric(m))
  #   error<-"m is not numeric"
  # if (!is.numeric(m0))
  #   error<-"m0 is not numeric"
  # if (!is.numeric(d))
  #   error<-"d is not numeric"
  # if (!is.numeric(D0))
  #   error<-"D0 is not numeric"
  
  #error handling: check if variables are not out of bounds
  if (mi <= 0){
    data_is_good<-FALSE
    error<-"m must be greater than zero"
    main_errors<-append(main_errors,error)}
  if (m0 <= 0){
    data_is_good<-FALSE
    error<-"m0 must be greater than zero"
    main_errors<-append(main_errors,error)}
  if (D0 <= 0){
    data_is_good<-FALSE
    error<-"D0 must be greater than zero"
    main_errors<-append(main_errors,error)}
  if (Di<= 0){
    data_is_good<-FALSE
    error<-"Di must be greater than zero"
    main_errors<-append(main_errors,error)}
  #error handling: check if cutoff values not null and numeric
  if (is.na(xlow) || is.null(xlow) || !is.numeric(xlow)){
    data_is_good<-FALSE
    error<-"malformed xlow (must be numeric)"
    main_errors<-append(main_errors,error)}
  if (is.na(xhigh) || is.null(xhigh) || !is.numeric(xhigh)){
    data_is_good<-FALSE
    error<-"malformed xhigh (must be numeric)"
    main_errors<-append(main_errors,error)}
  #error handling: check if cutoff values are numeric
  # if (!is.numeric(xlow))
  #   error<-"xlow is not numeric"
  # if (!is.numeric(xhigh))
  #   error<-"xhigh is not numeric"
  
  #error handling: check if cutoff values are not out of bounds
  if ((xlow < 0.0) || (xlow > 1.0)){
    data_is_good<-FALSE
    error<-"xlow must be in the interval [0,1]"
    main_errors<-append(main_errors,error)}
  if ((xhigh < 0.0) || (xhigh > 1.0)){
    data_is_good<-FALSE
    error<-"xhigh must be in the interval [0,1]"
    main_errors<-append(main_errors,error)
  }
  if (xhigh <= xlow){
    data_is_good<-FALSE
    error<-"xhigh must be greater than xlow"
    main_errors<-append(main_errors,error)}
  cutoff_f1 <- xlow
  cutoff_f2 <- xhigh
  trim <- ((f >= cutoff_f1) & (f <= cutoff_f2))
  
  #error handling: check if trimmed f has enough data points
  if (sum(trim) < 3){
    data_is_good<-FALSE
    error<-"Argument x: at least 3 data rows are required to constrain the fit;
         after trimming the frequency range there are less than 3.
         Consider changing cutoff values or not trimming your data."
    main_errors<-append(main_errors,error)}
  # data_is_good<-TRUE
  #writeLines("data is good: TRUE\n"
  
  return(list(data=data_is_good,print_errors=main_errors))
}

