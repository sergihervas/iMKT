#' @title Check input data
#' 
#' @description Check input data and return detailed errors when it is malformed
#' 
#' @details Check input data used in most package's functions (arguments daf, divergence, xlow and xhigh) and return a brief description of the error(s) found. This function is called within each analysis function (standardMK, FWW, DGRP, asymptoticMK, iMK) and if data does not pass checkInput() without errors, the requested analysis is not performed. 
#'
#' @param daf data frame containing DAF, Pi and P0 values
#' @param divergence data frame containing divergent and analyzed sites for selected (i) and neutral (0) classes
#' @param xlow lower limit for asymptotic alpha fit
#' @param xhigh higher limit for asymptotic alpha fit
#'
#' @import utils
#' @import stats
#'
#' @export

checkInput <- function(daf, divergence, xlow, xhigh){
    
    dataIsGood <- TRUE
    mainErrors <- "Your input files have the following errors: "

    ## Check data formatting
    if (NCOL(daf) != 3){
        dataIsGood <- FALSE
        error <- "Argument daf does not contain three tab-separated columns (daf, Pi, P0)."
        mainErrors <- append(mainErrors,error) }
    if(NROW(daf) < 10){ 
        dataIsGood <- FALSE
        error <- "Argument daf contains < 10 categories. At least 10 categories are required. Check the header or change you input!"
        mainErrors <- append(mainErrors,error) }
    if (NCOL(divergence) != 4){
        dataIsGood <- FALSE
        error <- "Argument divergence does not contain four tab-separated columns (mi, Di, m0, D0)."
        mainErrors <- append(mainErrors,error) }
    if (NROW(divergence) != 1){
        dataIsGood <- FALSE
        error <- "Argument divergence contains no data rows or it is malformed. Check the header!"
        mainErrors <- append(mainErrors,error) }

    ## Check proper names to the columns of daf and divergence. As we collect values based on name, data must have names well defined!
    if (!("daf" %in% colnames(daf) & "Pi" %in% colnames(daf) & "P0" %in% colnames(daf))){
        stop("Daf file doesn't have the correct column names: daf, Pi, P0.") }
    if (!("mi" %in% colnames(divergence) & "Di" %in% colnames(divergence) & "m0" %in% colnames(divergence) & "D0" %in% colnames(divergence))){
        stop("Divergence file doesn't have the correcto column names: mi, Di, m0, D0.") }

    ## Parse the data from argument daf
    f <- daf$daf ## Derived alelle frequencies
    Pi <- daf$Pi ## Non-synonymous polymorphism 
    P0 <- daf$P0 ## Synonymous polymorphism

    ## Check NAs and numeric for daf
    if (any(is.na(f))){
        dataIsGood <- FALSE
        error <- "Daf contains NA values (not allowed)."
        mainErrors <- append(mainErrors,error) }
    if (any(is.na(Pi))){
        dataIsGood <- FALSE
        error <- "Pi contains NA values (not allowed)."
        mainErrors <- append(mainErrors,error) }
    if (any(is.na(P0))){
        dataIsGood <- FALSE
        error <- "P0 contains NA values (not allowed)."
        mainErrors <- append(mainErrors,error) }
    if (!is.numeric(f)){
        dataIsGood <- FALSE
        error <- "Daf is not numeric"
        mainErrors <- append(mainErrors,error) }
    if (!is.numeric(Pi)){
        dataIsGood <- FALSE
        error <- "Pi is not numeric"
        mainErrors<-append(mainErrors,error) }
    if (!is.numeric(P0)){
        dataIsGood <- FALSE
        error <- "P0 is not numeric"
        mainErrors<-append(mainErrors,error) }

    ## Check if daf variables are not out of bounds.
    if (any(f < 0.0) || any(f > 1.0)){
        dataIsGood <- FALSE
        error <- "Daf contains values out of the required range [0,1]."
        mainErrors <- append(mainErrors,error) }
    if (all(f == 0)){
        dataIsGood <- FALSE
        error <- "Daf contains all values == 0 (not allowed)."
        mainErrors <- append(mainErrors,error) }
    if (any(Pi < 0) ){ 
        dataIsGood <- FALSE
        error <- "Pi contains values < 0 (not allowed)."
        mainErrors <- append(mainErrors,error) }
    if (all(Pi == 0)){
        warning("Pi contains all values == 0.") }
    if (any(P0 < 0) || all(P0 == 0)){ 
        dataIsGood <- FALSE
        error <- "P0 contains values < 0 (not allowed) or P0 contains all values == 0 (not allowed)."
        mainErrors <- append(mainErrors,error) } 
    
    # if (any(P0 == 0)){  ## Only for asymptoticMK and iMK
    #    warning("Input daf file contains P0 values = 0.\nThis can bias the function fitting and the estimation of alpha.")}

    ## Check if daf argument has enough data points
    if (NROW(daf) < 3){
        dataIsGood <- FALSE
        error <- "Argument daf: at least three data rows are required to constrain the fit."
        mainErrors <- append(mainErrors,error) }

    ## Parse the data from argument divergence
    mi <- divergence$mi ## Number of non-synonymous sites
    m0 <- divergence$m0 ## Number of synonymous sites
    Di <- divergence$Di ## Non-synonymous divergence
    D0 <- divergence$D0 ## Synonymous divergence

    ## Check NAs and numeric for divergence
    if (is.na(mi) || !is.numeric(mi)){
        dataIsGood <- FALSE
        error <- "Malformed mi (must be numeric)."
        mainErrors <- append(mainErrors,error) }
    if (is.na(m0) || !is.numeric(m0)){
        dataIsGood <- FALSE
        error <- "Malformed m0 (must be numeric)."
        mainErrors <- append(mainErrors,error) }
    if (is.na(D0) || !is.numeric(D0)){
        dataIsGood <- FALSE
        error <- "Malformed D0 (must be numeric)."
        mainErrors <- append(mainErrors,error) }
    if (is.na(Di) || !is.numeric(Di)){
        dataIsGood <- FALSE
        error <- "Malformed Di (must be numeric)."
        mainErrors <- append(mainErrors,error) }
    if (is.na(xlow) || is.null(xlow)){
        dataIsGood <- FALSE
        error <- "Malformed xlow (must be numeric)."
        mainErrors <- append(mainErrors,error) }
    if (is.na(xhigh) || is.null(xhigh)){
        dataIsGood <- FALSE
        error <- "Malformed xhigh (must be numeric)."
        mainErrors <- append(mainErrors,error) }

    ## Check if divergence variables are not out of bounds
    if (mi <= 0){
        dataIsGood <- FALSE
        error <- "mi must be greater than zero."
        mainErrors <- append(mainErrors,error) }
    if (m0 <= 0){
        dataIsGood <- FALSE
        error <- "m0 must be greater than zero."
        mainErrors <- append(mainErrors,error) }
    if (D0 < 0){
        dataIsGood <- FALSE
        error <- "D0 must be greater than zero."
        mainErrors <- append(mainErrors,error) }
    if (D0 == 0) {
        warning("D0 == 0.") }
    if (Di <= 0){
        dataIsGood <- FALSE
        error <- "Di must be greater than zero."
        mainErrors <- append(mainErrors,error) }

    ## Check if number of sites (m, m0) is not higher than divergenge (d, D0) + polimorphisms (p|p0)
    if (Di > mi || sum(Pi) > mi || sum(Pi) + Di > mi){
        dataIsGood <- FALSE
        error <- "mi must be higher than Pi, Di and sum(Pi+Di)."
        mainErrors <- append(mainErrors,error) }
    if (D0 > m0 || sum(P0) > m0 || sum(P0) + D0 > m0){
        dataIsGood <- FALSE
        error <- "m0 must be higher than P0, D0 and sum(P0+D0)."
        mainErrors <- append(mainErrors,error) }

    ## Check if cutoff values not null and numeric
    if (is.na(xlow) || is.null(xlow) || !is.numeric(xlow)){
        dataIsGood <- FALSE
        error <- "Malformed xlow (must be numeric)."
        mainErrors <- append(mainErrors,error) }
    if (is.na(xhigh) || is.null(xhigh) || !is.numeric(xhigh)){
        dataIsGood <- FALSE
        error <- "Malformed xhigh (must be numeric)."
        mainErrors <- append(mainErrors,error) }

    ## Check if cutoff values are not out of bounds
    if ((xlow < 0.0) || (xlow > 1.0)){
        dataIsGood <- FALSE
        error <- "xlow must be in the interval [0,1]."
        mainErrors <- append(mainErrors,error) }
    if ((xhigh < 0.0) || (xhigh > 1.0)){
        dataIsGood <- FALSE
        error <- "xhigh must be in the interval [0,1]."
        mainErrors <- append(mainErrors,error) }
    if (xhigh <= xlow){
        dataIsGood <- FALSE
        error <- "xhigh must be greater than xlow."
        mainErrors <- append(mainErrors,error) }

    cutoff_f1 <- xlow
    cutoff_f2 <- xhigh
    trim <- ((f >= cutoff_f1) & (f <= cutoff_f2))

    ## Check if trimmed f has enough data points
    if (sum(trim) < 3){
        dataIsGood <- FALSE
        error <- "Argument x: at least 3 data rows are required to constrain the fit;
         after trimming the frequency range there are less than 3.
         Consider changing cutoff values or not trimming your data."
        mainErrors <- append(mainErrors,error) }
    
    ## Return T/F + errors
    return(list(data=dataIsGood, print_errors=mainErrors))
}

