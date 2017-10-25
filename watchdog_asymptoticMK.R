## ERROR HANDLING
## function which checks input data for asymptoticMK
## should check and discuss many things. Use as an example for now!
## mix between Messer & Haller and Vasily work

watchdog <- function(x, y, xlow, xhigh){
  
  data_is_good <- FALSE
  
  ######### Detecting if the colnames are the correct ones ############
  if (!("daf" %in% colnames(x) & "pS" %in% colnames(x) & "pN" %in% colnames(x)))
    stop("required header doesn't have the correct names: daf, pN, pS")
  
  #error handling: check data formatting
  if (NCOL(x) != 3)
    stop("argument x does not contain exactly three tab-separated columns")
  if (NROW(x) <= 0)
    stop("argument x contains no data rows")
  if (NCOL(y) != 6)
    stop("argument y does not contain exactly six tab-separated columns")
  if (NROW(y) <= 0)
    stop("argument y contains no data rows")
  
  ##assign proper names to the columns of x and y. MUST CHECK THAT!
  ####Not defined###### 
  #names(x) <- c("daf", "pN", "pS")
  #names(y) <- c("Chr",  "Pop",	"m0f",	"D0f",	"m4f",	"D4f")
  
  #error handling: input colnames 
  #for arg x
  x.cols <- names(x)
  suppressWarnings(#the goal is to generate NAs here, so we dont need no warnings
    if (!is.na(as.numeric(x.cols[1])) ||
          !is.na(as.numeric(x.cols[2])) ||
          !is.na(as.numeric(x.cols[3])))
      stop("argument x has a numeric column name;
           probably the required header row is missing"))
  #for arg y
  y.cols <- names(y)
  suppressWarnings(
    if (!is.na(as.numeric(y.cols[1])) ||
          !is.na(as.numeric(y.cols[2])) ||
          !is.na(as.numeric(y.cols[3])) ||
          !is.na(as.numeric(y.cols[4])) ||
          !is.na(as.numeric(y.cols[5])) ||
          !is.na(as.numeric(y.cols[6])))
      stop("argument y has a numeric column name;
           probably the required header row is missing")
      )
  
  ##parse the data from argument x
  ####### dependent on the column name
  f <- x$daf #derived alelle frequencies
  p <- x$pN #non-synonymous polymorphism 
  p0 <- x$pS #synonymous polymorphism
  
  #error handling: check if variables are good
  if (any(is.na(f)))
    stop("f contains NA values (not allowed)")
  if (any(is.na(p)))
    stop("p contains NA values (not allowed)")
  if (any(is.na(p0)))
    stop("p0 contains NA values (not allowed)")
  if (is.null(f))
    stop("f malformed (must be numeric)")
  if (is.null(p))
    stop("p malformed (must be numeric)")
  if (is.null(p0))
    stop("p0 malformed (must be numeric)")
  
  #error handling: check if variables are numeric
  if (!is.numeric(f))
    stop("f is not numeric")
  if (!is.numeric(p))
    stop("p is not numeric")
  if (!is.numeric(p0))
    stop("p0 is not numeric")
  
  #error handling: check if variables are not out of bounds
  if (any(f < 0.0) || any(f > 1.0))
    stop("f contains values out of the required range [0,1]")
  if (all(f == 0))
    stop("f contains all values == 0 (not allowed)")
  if (any(p < 0))   # note that zero is allowed, although not recommended
    stop("p contains values < 0 (not allowed)")
  if (all(p == 0))    # not all can be zero, however
    stop("p contains all values == 0 (not allowed)")
  if (any(p0 <= 0))
    stop("p0 contains values <= 0 (not allowed)")
  if (all(p0 == 0))
    stop("p0 contains all values == 0 (not allowed)")
  
  #error handling: check if argument x has enough data points
  if (NROW(x) < 3)
    stop("argument x: at least three data rows are required to constrain the fit")
  
  ##parse the data from argument y and force it to be numeric...
  m <- as.numeric(y$m0f) #number of non-synonymous sites   
  m0 <- as.numeric(y$m4f) ##number of synonymous sites
  d <- as.numeric(y$D0f) #non-synonymous divergence
  d0 <- as.numeric(y$D4f) #synonymous divergence
  
  #error handling: check if variables are good
  if (is.na(m) || is.null(m))
    stop("malformed m (must be numeric)")
  if (is.na(m0) || is.null(m0))
    stop("malformed m0 (must be numeric)")
  if (is.na(d0) || is.null(d0))
    stop("malformed d0 (must be numeric)")
  if (is.na(d) || is.null(d))
    stop("malformed d (must be numeric)")
  if (is.na(xlow) || is.null(xlow))
    stop("malformed xlow (must be numeric)")
  if (is.na(xhigh) || is.null(xhigh))
    stop("malformed xhigh (must be numeric)")
    
  # Checks if number of sites (m, m0) is not higher than divergenge (d, d0), the sum of the polimorphisms (p|p0) or the divergenge + the sum of the polimorphsms
  if (d>m || sum(p)>m || sum(p)+d>m)
    stop("m must be higher than p, d and p+d")
  if (d0>m0 || sum(p0)>m0 || sum(p)+d0>m0)
    stop("m0 must be higher than p0, d0 and p0+d0")
  
  #error handling: check if variables are numeric
  #it makes no sence as we forced them to be numeric... should check that!
  if (!is.numeric(m))
    stop("m is not numeric")
  if (!is.numeric(m0))
    stop("m0 is not numeric")
  if (!is.numeric(d))
    stop("d is not numeric")
  if (!is.numeric(d0))
    stop("d0 is not numeric")
  
  #error handling: check if variables are not out of bounds
  if (m <= 0)
    stop("m must be greater than zero")
  if (m0 <= 0)
    stop("m0 must be greater than zero")
  if (d0 <= 0)
    stop("d0 must be greater than zero")
  if (d <= 0)
    stop("d must be greater than zero")
  
  #error handling: check if cutoff values not null
  if (is.na(xlow) || is.null(xlow))
    stop("malformed xlow (must be numeric)")
  if (is.na(xhigh) || is.null(xhigh))
    stop("malformed xhigh (must be numeric)")
  
  #error handling: check if cutoff values are numeric
  if (!is.numeric(xlow))
    stop("xlow is not numeric")
  if (!is.numeric(xhigh))
    stop("xhigh is not numeric")
  
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
  
  data_is_good <- TRUE
  #writeLines("data is good: TRUE\n")
  
  return(data_is_good)
}

