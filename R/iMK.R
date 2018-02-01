#' @title integrative MKT method
#' 
#' @description iMK: MKT using asymptoticMK method and estimation of negative selection fractions (d, b, f)
#'
#' @details The integrative MKT (iMKT) allows the estimation of the rate of adaptive evolution (alpha) and the diverse negative selection regimens. iMKT uses asymptotic MK method (Messer and Petrov 2012 PNAS; Haller and Messer 2017 G3) to estimate alpha and the diverse negative selection fractions (d: strongly deleterious, b: weakly deleterious, f: neutral), based on the assumption that weakly deleterious mutations usually do not reach high allele frequencies and therefore, produce the underestimation of alpha at low DAF categories. The fraction of strongly deleterious mutations is estimated as the difference between neutral (0) and selected (i) polymorphic sites relative to the number of analyzed sites: d = 1 - (P0/m0 / Pi/mi). The fraction of weakly deleterious sites (b) corresponds to the relative proportion of selected polymorphic sites that cause the underestimation of alpha at low DAF categories. Finally, the fraction of neutral sites (f) is estimated as: f = 1 - d - b. iMK() only fits an exponential model for the computation of alpha.
#' 
#' @param daf data frame containing DAF, Pi and P0 values
#' @param divergence data frame containing divergent and analyzed sites for selected (i) and neutral (0) classes
#' @param xlow lower limit for asymptotic alpha fit
#' @param xhigh higher limit for asymptotic alpha fit
#' @param seed seed value (optional). No seed by default
#' @param plot report plots of daf, alpha and negative selection fractions (optional). Default is FALSE
#'
#' @return iMKT method. List with asymptotic MK table and values, fractions of sites and graphs of DAF, asymptotic alpha model and negative selection fractions (optional).
#'
#' @examples
#' ## Without plot
#' iMK(myDafData, myDivergenceData, xlow=0, xhigh=0.9)
#' ## With plot
#' iMK(myDafData, myDivergenceData, xlow=0, xhigh=0.9, plot=TRUE)
#' 
#' @import utils
#' @import stats
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom ggthemes theme_foundation
#' @importFrom cowplot plot_grid
#' 
#' @keywords MKT
#' @export

iMK <- function(daf, divergence, xlow, xhigh, seed, plot=FALSE) {
  
  ## Check data
  check <- checkInput(daf, divergence, xlow, xhigh)
  if (check$data == FALSE) {
    stop (check$print_errors) }
  
  ## Warning P0
  if (any(daf$P0 == 0)){
    warning("Input daf file contains P0 values = 0.\nThis can bias the function fitting and the estimation of alpha.") }
  
  ## Check seed
  if(missing(seed)) {
    seed <- NULL
  } else {
    set.seed(seed)
  }

  ## Create MKT table standard
  mkt_table_standard <- data.frame(Polymorphism = c(sum(daf$P0), sum(daf$Pi)), 
                                   Divergence=c(divergence$D0,divergence$Di),
                                   row.names = c("Neutral class","Selected class"))
  
  ## Total number of sites analyzed 
  mi <- divergence$mi
  m0 <- divergence$m0
  
  ## Run asymptotic MK and retrieve alphas 
  ## iMK only fits exponential model in asymptotic alpha. it uses asymptoticMKExp()
  
  asymptoticMKExp <- function(daf, divergence, xlow, xhigh, seed) {
    
    ## Check data
    check <- checkInput(daf, divergence, xlow, xhigh)
    if(check$data == FALSE) {
      stop(check$print_errors) }
    
    if (any(daf$P0 == 0)){ ## Warning P0
      warning("Input daf file contains P0 values = 0.\nThis can bias the function fitting and the estimation of alpha.")}
    
    ## Check seed
    if(missing(seed)) {
      seed <- NULL
    } else {
      set.seed(seed)
    }
    
    ## Parse the data from argument x
    f <- daf$daf #derived alelle frequencies
    p <- daf$Pi #non-synonymous polymorphism 
    p0 <- daf$P0 #synonymous polymorphism
    
    ## Parse the data from argument y
    m <- divergence$mi #number of non-synonymous analyzed positions   
    m0 <- divergence$m0 ##number of synonymous analyzed positions
    d <- divergence$Di #non-synonymous divergence
    d0 <- divergence$D0 #synonymous divergence
    
    ## Compute alpha values and trim
    alpha <- 1 - (d0/d) * (p/p0)
    cutoff_f1 <- xlow
    cutoff_f2 <- xhigh
    trim <- ((f >= cutoff_f1) & (f <= cutoff_f2))
    f_trimmed <- f[trim]
    alpha_trimmed <- alpha[trim]
    
    ## Compute the original MK alpha
    alpha_nonasymp <- 1 - (d0/d) * (sum(p[trim])/sum(p0[trim])) #using trimmed values
    
    ## Two-step nls2() model fit at a given level of precision (res)
    fitMKmodel <- function(alpha_trimmed, f_trimmed, res) {
      
      ## First fitting using starting values (st)
      mod <- tryCatch({
        
        ## Starting values to fit the model  
        st <- expand.grid(const_a=seq(-1,1,length.out=res + 1), const_b=seq(-1,1,length.out=res), const_c=seq(1,10,length.out=res + 1))
        
        ## Fitting
        nls2(alpha_trimmed ~ const_a + const_b * exp(-const_c* f_trimmed), start=st, algorithm="brute-force", control=nls.control(maxiter=NROW(st)))
        
      }, error=function(cond) {}) ## Return condition of error when unable to fit
      
      ## If mod fails...
      if (length(mod) == 0) { return(NULL) }
      
      ## Second fitting, starting from previous fit (mod)
      mod2 <- tryCatch({
        nls2(alpha_trimmed ~ const_a + const_b * exp(-const_c* f_trimmed), start = mod, control=nls.control(maxiter=200))
        
      }, error=function(cond) {}) ## Same error handling than the previous step
      
      ## If mod2 fails...
      if (length(mod2) == 0) { return(NULL) }
      
      ## Return mod2 if fitted
      return(mod2)
    }
    
    mod1 <- fitMKmodel(alpha_trimmed, f_trimmed, 10)
    
    ## If mod1 did not work, try a deeper scan for a decent fit (res=20)
    if (length(mod1) == 0) {
      mod1 <- fitMKmodel(alpha_trimmed, f_trimmed, 20)
    } 
    
    tryCatch({
      mod2 <- lm(alpha_trimmed ~ f_trimmed)
    }, error=function(cond) {})
    
    ## Compute confidence intervals of alpha using predictNLS 
    ## Get a CI using Monte Carlo simulation based upon a fitted model.  
    ## Thanks to Andrej-Nikolai Spiess (http://www.dr-spiess.de) for this code.
    predictNLS <- function(object, newdata, level = 0.95, nsim = 10000) {
      
      ## get right-hand side of formula
      RHS <- as.list(object$call$formula)[[3]]
      EXPR <- as.expression(RHS)
      
      ## all variables in model
      VARS <- all.vars(EXPR)
      
      ## coefficients
      COEF <- coef(object)
      
      ## extract predictor variable    
      predNAME <- setdiff(VARS, names(COEF))  
      
      ## take fitted values, if 'newdata' is missing
      if (missing(newdata)) {
        newdata <- eval(object$data)[predNAME]
        colnames(newdata) <- predNAME
      }
      
      ## check that 'newdata' has same name as predVAR
      if (names(newdata)[1] != predNAME) stop("newdata should have name '", predNAME, "'!")
      
      ## get parameter coefficients
      COEF <- coef(object)
      
      ## get variance-covariance matrix
      VCOV <- vcov(object)
      
      ## augment variance-covariance matrix for 'mvrnorm' 
      ## by adding a column/row for 'error in x'
      NCOL <- ncol(VCOV)
      ADD1 <- c(rep(0, NCOL))
      ADD1 <- matrix(ADD1, ncol = 1)
      colnames(ADD1) <- predNAME
      VCOV <- cbind(VCOV, ADD1)
      ADD2 <- c(rep(0, NCOL + 1))
      ADD2 <- matrix(ADD2, nrow = 1)
      rownames(ADD2) <- predNAME
      VCOV <- rbind(VCOV, ADD2) 
      
      ## iterate over all entries in 'newdata' as in usual 'predict.' functions
      NR <- nrow(newdata)
      respVEC <- numeric(NR)
      seVEC <- numeric(NR)
      varPLACE <- ncol(VCOV)   
      
      ## define counter function
      counter <- function(i) {
        if (i%%10 == 0) { cat(i) 
        } else { cat(".") }
        if (i%%50 == 0) { cat("\n") }
        flush.console()
      }
      
      ## create output matrix (df)
      outMAT <- NULL 
      
      for (i in 1:NR) {
        
        ## get predictor values and optional errors
        predVAL <- newdata[i, 1]
        if (ncol(newdata) == 2) predERROR <- newdata[i, 2] else predERROR <- 0
        names(predVAL) <- predNAME  
        names(predERROR) <- predNAME  
        
        ## create mean vector for 'mvrnorm'
        MU <- c(COEF, predVAL)
        
        ## create variance-covariance matrix for 'mvrnorm'
        ## by putting error^2 in lower-right position of VCOV
        newVCOV <- VCOV
        newVCOV[varPLACE, varPLACE] <- predERROR^2
        
        ## create MC simulation matrix
        simMAT <- mvrnorm(n = nsim, mu = MU, Sigma = newVCOV, empirical = TRUE)
        
        ## evaluate expression on rows of simMAT
        EVAL <- try(eval(EXPR, envir = as.data.frame(simMAT)), silent = TRUE)
        if (inherits(EVAL, "try-error")) stop("There was an error evaluating the simulations!")
        
        ## collect statistics
        PRED <- data.frame(predVAL)
        colnames(PRED) <- predNAME   
        FITTED <- predict(object, newdata = data.frame(PRED))
        MEAN.sim <- mean(EVAL, na.rm = TRUE)
        SD.sim <- sd(EVAL, na.rm = TRUE)
        MEDIAN.sim <- median(EVAL, na.rm = TRUE)
        MAD.sim <- mad(EVAL, na.rm = TRUE)
        QUANT <- quantile(EVAL, c((1 - level)/2, level + (1 - level)/2))
        RES <- c(FITTED, MEAN.sim, SD.sim, MEDIAN.sim, MAD.sim, QUANT[1], QUANT[2])
        outMAT <- rbind(outMAT, RES)
      }
      
      colnames(outMAT) <- c("fit", "mean", "sd", "median", "mad", names(QUANT[1]), names(QUANT[2]))
      rownames(outMAT) <- NULL   
      return(outMAT)      
    }
    
    tryCatch({
      ci_pred <- predictNLS(mod1, newdata=data.frame(f_trimmed=1.0))
      alpha_1_low <- ci_pred[6]
      alpha_1_high <- ci_pred[7]

      ## Preparation of ouput (alpha asym, a, b, c)
      alpha_1_est <- predict(mod1, newdata=data.frame(f_trimmed=1.0))
      const_a <- coef(mod1)["const_a"]
      const_b <- coef(mod1)["const_b"]
      const_c <- coef(mod1)["const_c"]
      
      ## Output table
      result_df <- data.frame(model="exponential", a=const_a, b=const_b, c=const_c, alpha_asymptotic=alpha_1_est, CI_low=alpha_1_low, CI_high=alpha_1_high, alpha_original=alpha_nonasymp, row.names=NULL)
      return(result_df)
    }, error=function(cond) {cat("Could not fit exponential model for the computation of asymptotic alpha.\n")})
  }
  
  asymptoticMK_table <- asymptoticMKExp(daf, divergence, xlow, xhigh)
  alpha_asymptotic <- as.numeric(asymptoticMK_table$alpha_asymptotic)
  alpha_standard <- as.numeric(asymptoticMK_table$alpha_original)
  alpha_CI_low <- asymptoticMK_table$CI_low 
  
  ## Estimate the relative proportion of non-synonymous and synonymous substitutions
  daf$N <- daf$Pi/sum(daf$Pi)   
  daf$S <- daf$P0/sum(daf$P0)
   
  ## Estimate alpha for each DAF category
  daf$alpha <- 1-((mkt_table_standard[1,2]*daf$Pi)/(mkt_table_standard[2,2]*daf$P0))
  
  ## Estimate the synonymous and non-synonymous ratio
  synonymous_ratio <- mkt_table_standard[1,1]/m0    
  nonsynonymous_ratio <- mkt_table_standard[2,1]/mi  
  
  ## Estimate the fraction of neutral sites (f)
  f <- nonsynonymous_ratio/synonymous_ratio
  
  ## Estimate the fraction of strongly deleleterious sites (d)
  d <- 1-f   
  
  ## Estimate the fraction of weakly deleterious sites (b)
  daf <- na.omit(daf); daf <- droplevels(daf)
  wd <- 0
  for (i in 1:nrow(daf)) {
    row <- daf[i,]
    if (row$alpha < alpha_CI_low) {
      wd <- wd + ((alpha_asymptotic-row$alpha)*row$N)
    } else { break }
  }
  wd <- wd/(alpha_asymptotic-min(daf$alpha,na.rm=T))
  b <- wd*f
  
  ## Re-estimate the truly number of neutral sites, removing the slightly deleterious 
  f <- f-b
  
  ## Fraction of f, b and d sites
  fraction <- data.frame(Fraction=c(d,f,b), Type=c("d","f","b"), MKT=rep("Asymptotic MK", 3))
  
  ## Perform plots 
  if(plot == TRUE) {
    
    ## plot fractions

   plotfraction <- ggplot(fraction) + geom_bar(stat="identity", aes_string(x="MKT", y="Fraction", fill="Type"), color="black") +
      coord_flip() + themePublication() + ylab(label="Fraction") + xlab(label="Cut-off") +
      scale_fill_manual(values=c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33"), breaks=c("f","d","b"), labels=c(expression(italic("f")),expression(italic("d")),expression(italic("b")))) +
      theme(axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.line.x =element_blank())  + scale_y_discrete(limit=seq(0,1,0.25), expand=c(0,0))
    
    ## DAF graph
    daf_graph <- daf[c("daf","Pi","P0")]
    daf_graph<-melt(daf_graph,id.vars = "daf")
  
    plotdaf <-  ggplot(daf_graph) +
      geom_point(aes_string(x="daf", y="value", color="variable"), size=3) +
      themePublication() + scale_color_manual(values=c("#386cb0","#fdb462"), name="Type", breaks=c("Pi","P0"), labels=c("Non-synonymous","Synonymous")) +
      xlab("Derived Allele Frequency") + ylab("Number of Sites")
  
    ## Alpha graph
    y1 <- function(daf) {
      asymptoticMK_table$a+asymptoticMK_table$b*exp(-asymptoticMK_table$c*daf) }
    xs <- seq(xlow, xhigh, length.out=nrow(daf)*5)
    ysmax <- rep(asymptoticMK_table$alpha_asymptotic, length(xs))
    ysmin <- y1(xs)
    shader_df <- data.frame(xs, ysmin, ysmax)
 
    plot_alpha <- ggplot(daf, aes_string(x="daf", y="alpha"))  +
      ## Confidence intervals
      geom_rect(data=data.frame(xmin=-Inf, xmax=Inf, ymin=asymptoticMK_table$CI_low, ymax=asymptoticMK_table$CI_high), aes_string(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax"), fill="gray30", alpha=0.4, inherit.aes=F) +
      ## Function curve
      stat_function(fun=y1, color="#ef3b2c", size=2) +
      ## Asymptotic alpha
      geom_hline(yintercept=asymptoticMK_table$alpha_asymptotic, color="#662506", linetype="dashed") +  
      ## Alpha derived via classic MKT
      geom_hline(yintercept=asymptoticMK_table$alpha_original, color="#386cb0", linetype="dashed") +
      ## Cut-offs
      geom_vline(xintercept=c(xlow, xhigh), color="gray10", linetype="dotted") +    
      ## Points
      geom_point(size=3, color="gray15") +
      ## Shade the fraction of WDMs
      geom_ribbon(data=shader_df, aes_string(x="xs", ymin="ysmin", ymax="ysmax"), fill="gray30", alpha=0.2, inherit.aes=F) + 
      ## Customization
      themePublication() + 
      xlab("Derived allele frequency") + ylab(expression(bold(paste("Adaptation (",alpha,")")))) +
      ## Alphas labels
      annotate("text", x=xhigh-0.2, y=asymptoticMK_table$alpha_asymptotic-0.2, label=paste0('alpha [asymptotic] == ', round(asymptoticMK_table$alpha_asymptotic, digits = 3)), parse=T, color="#662506", size=4) +
      annotate("text",x=xhigh-0.2, y=asymptoticMK_table$alpha_original-0.1, label=paste0('alpha [standard] == ', round(asymptoticMK_table$alpha_original,digits = 3)), parse=T, color="#386cb0", size=4)
  
    ## Render plots with labels
    plots_iMKT <- plot_grid(plotdaf, plot_alpha, plotfraction, nrow=3,  labels=c("A","B","C"), rel_heights=c(2,2,1))
    
    ## Store output
    ## iMKT output
    fraction <- fraction[c("Type","Fraction")]
    asymptoticMK_table[2:8] <- round(asymptoticMK_table[2:8],4)
    output <- list(asymptoticMK_table, fraction, plots_iMKT)
    names(output) <- c("Asymptotic MK table", "Fractions of sites", "Graphs")
    
  ## If no plot to perform  
  } else if (plot == FALSE) {
    ## iMKT output
    fraction <- fraction[c("Type","Fraction")]
    asymptoticMK_table[2:8] <- round(asymptoticMK_table[2:8],4)
    output <- list(asymptoticMK_table, fraction)
    names(output) <- c("Asymptotic MK table", "Fractions of sites")
  }
  
  ## Return output
  return(output) 
}
