## ERROR HANDLING
## function which checks input data for asymptoticMK
## should check and discuss many things. Use as an example for now!
## mix between Messer & Haller and Vasily work

check_input <- function(x, y, xlow, xhigh){
  
  data_is_good <- FALSE
  
  ######### Detecting if the colnames are the correct ones ############
  # if (!("daf" %in% colnames(x) & "pS" %in% colnames(x) & "pN" %in% colnames(x)))
  #   stop("required header doesn't have the correct names: daf, pN, pS")
  
  #error handling: check data formatting
  if (NCOL(x) != 3 || NCOL(x) <3)
    stop("Argument x does not contain three tab-separated columns.Argument subsitute(x), at least three data rows are required to constrain the fit")
  if (NROW(x) < 0) 
    stop("Argument x contains no data rows or x categories are not correct define. Check the header!")
  if(NROW(x)!=20 & NROW(x)!=10)
    stop("Daf categories in x are not correct define. Check the header!")
  if (NCOL(y) != 6)
    
    stop("Argument y does not contain six tab-separated columns")
  if (NROW(y) <= 0 & (NROW(y)<3))
    
    stop("Argument y contains no data rows")
  
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
  #     stop("argument x has a numeric column name;
  #          probably the required header row is missing"))
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
  #     stop("argument y has a numeric column name;
  #          probably the required header row is missing")
  #     )
  
  
  ##parse the data from argument x
  ####### dependent on the column name
  f <- x$daf #derived alelle frequencies
  Pi <- x$Pi #non-synonymous polymorphism 
  P0 <- x$P0 #synonymous polymorphism
  
  #error handling: check if variables are good
  if (any(is.na(f)))
    
    stop("f contains NA values (not allowed)")
  if (any(is.na(Pi)))
    
    stop("pi contains NA values (not allowed)")
  if (any(is.na(P0)))
    
    stop("p0 contains NA values (not allowed)")
  #error handling: check if variables are numeric
  if (!is.numeric(f))
    
    stop("f is not numeric")
  if (!is.numeric(Pi))
    
    stop("pi is not numeric")
  if (!is.numeric(P0))
    
    stop("p0 is not numeric")
  
  ##JUST CHECK IF THE VARIABLE IS NUMERIC, IF NOT, THEN stop. NOT NEED NULL. JUST stop IT IF AREN'T NUMERIC
  # if (is.null(f))
  #   stop("f malformed (must be numeric)")
  # if (is.null(pi))
  #   stop("p malformed (must be numeric)")
  # if (is.null(p0))
  #   stop("p0 malformed (must be numeric)")
  # 
  # #error handling: check if variables are numeric
  # if (!is.numeric(f))
  #   stop("f is not numeric")
  # if (!is.numeric(Pi))
  #   stop("p is not numeric")
  # if (!is.numeric(P0))
  #   stop("p0 is not numeric")
  
  #error handling: check if variables are not out of bounds.
  if (any(f < 0.0) || any(f > 1.0))
    
    stop("f contains values out of the required range [0,1]")
  if (all(f == 0))
    
    stop("f contains all values == 0 (not allowed)")
  if (any(Pi < 0)||all(Pi == 0))   # note that zero is allowed, although not recommended
    
    stop("Pi contains values < 0 (not allowed) Pi contains all values == 0 (not allowed)")
  # if (all(Pi == 0))    # not all can be zero, however
  #   stop("p contains all values == 0 (not allowed)")
  if (any(P0 <= 0)||all(P0 == 0))
    
    stop("P0 contains values <= 0 (not allowed) or P0 contains all values == 0 (not allowed)")
  # if (all(p0 == 0))
  #   stop("p0 contains all values == 0 (not allowed)")
  
  
  #Incorpore at begging when check the input length
  # error handling: check if argument x has enough data points
  # if (NROW(x) < 3)
  # stop("Argument x: at least three data rows are required to constrain the fit")
  # 
  ##parse the data from argument y and force it to be numeric... NOT SENSE CHECK THEN IS.NULL OR NUMERIC
  # NOT FORCE TO BE NUMERIC, CHECK IF VALUES ARE NOT NUMERIC [not force+not is.null]
  mi <- as.numeric(y$mi) #number of non-synonymous sites
  m0 <- as.numeric(y$m0) ##number of synonymous sites
  Di <- as.numeric(y$Di) #non-synonymous divergence
  D0 <- as.numeric(y$D0) #synonymous divergence
  
  #error handling: check if variables are good
  if (is.na(mi) || !is.numeric(mi))
    
    stop("malformed m (must be numeric)")
  if (is.na(m0) || !is.numeric(m0))
    
    stop("malformed m0 (must be numeric)")
  if (is.na(D0) || !is.numeric(D0))
    
    stop("malformed D0 (must be numeric)")
  if (is.na(Di) || !is.numeric(Di))
    
    stop("malformed d (must be numeric)")
  if (is.na(xlow) || is.null(xlow))
    
    stop("malformed xlow (must be numeric)")
    
  if (is.na(xhigh) || is.null(xhigh))
    stop("malformed xhigh (must be numeric)")
  
  # Checks if number of sites (m, m0) is not higher than divergenge (d, D0), the sum of the polimorphisms (p|p0) or the divergenge + the sum of the polimorphsms
  if (Di>mi || sum(Pi)>mi || sum(Pi)+Di>mi)
    
    stop("mi must be higher than Pi, Di and sum(Pi+Di)")
  if (D0>m0 || sum(P0)>m0 || sum(P0)+D0>m0)
    
    stop("m0 must be higher than P0, D0 and sum(P0+D0)")
  
  # #error handling: check if variables are numeric
  # #it makes no sence as we forced them to be numeric... should check that!
  # if (!is.numeric(m))
  #   stop("m is not numeric")
  # if (!is.numeric(m0))
  #   stop("m0 is not numeric")
  # if (!is.numeric(d))
  #   stop("d is not numeric")
  # if (!is.numeric(D0))
  #   stop("D0 is not numeric")
  
  #error handling: check if variables are not out of bounds
  if (mi <= 0)
    
    stop("m must be greater than zero")
  if (m0 <= 0)
    
    stop("m0 must be greater than zero")
  if (D0 <= 0)
    
    stop("D0 must be greater than zero")
  if (Di<= 0)
    
    stop("Di must be greater than zero")
  
  #error handling: check if cutoff values not null and numeric
  if (is.na(xlow) || is.null(xlow) || !is.numeric(xlow))
    
    stop("malformed xlow (must be numeric)")
  if (is.na(xhigh) || is.null(xhigh) || !is.numeric(xhigh))
    
    stop("malformed xhigh (must be numeric)")
  
  #error handling: check if cutoff values are numeric
  # if (!is.numeric(xlow))
  #   stop("xlow is not numeric")
  # if (!is.numeric(xhigh))
  #   stop("xhigh is not numeric")
  
  #error handling: check if cutoff values are not out of bounds
  if ((xlow < 0.0) || (xlow > 1.0))
    
    stop("xlow must be in the interval [0,1]")
  if ((xhigh < 0.0) || (xhigh > 1.0))
    
    stop("xhigh must be in the interval [0,1]")
  if (xhigh <= xlow)
    
    stop("xhigh must be greater than xlow")
  
  cutoff_f1 <- xlow
  cutoff_f2 <- xhigh
  trim <- ((f >= cutoff_f1) & (f <= cutoff_f2))
  
  #error handling: check if trimmed f has enough data points
  if (sum(trim) < 3)
    
    stop("Argument x: at least 3 data rows are required to constrain the fit;
         after trimming the frequency range there are less than 3.
         Consider changing cutoff values or not trimming your data.")
  
  data_is_good<-TRUE
  #writeLines("data is good: TRUE\n")
  
  return(data_is_good)
}

