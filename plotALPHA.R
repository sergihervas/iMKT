## alpha plot 
## input: x, y, alpha asymptotic output (iMK$asymptoticMK), xlow and xhigh for asymptoticMK (or iMK)

plotALPHA <- function(x, y, df, xlow, xhigh) {
  
  ## load packages
  require(ggplot2, quietly=TRUE)
  
  ## parse data: should include watchdog adapted
  ## Not defined colnames
  ##names(x) <- c("daf","pN","pS")
  d <- as.numeric(y$D0f)
  d0 <- as.numeric(y$D4f)
  
  #ci values
  CI_low <- df$CI_low
  CI_high <- df$CI_high
  #alpha asymp
  alpha_asymptotic <- as.numeric(df$alpha_asymptotic)
  alpha_original <- as.numeric(df$alpha_original)
  
  ## alpha for each bin
  x$alpha <- 1 - ((d0*x$pN)/(d*x$pS))
  
  ## plot
  #y = a + b * exp(-c*x) function written in a way ggplot can parse it
  y1 <- function(x) { df$a + df$b * exp(-df$c * x) }
  # prepare shader dataframe
  xs <- seq(xlow, xhigh, length.out=nrow(x)*5)
  ysmax <- rep(alpha_asymptotic, length(xs))
  ysmin <- y1(xs)
  shader_df <- data.frame(xs, ysmin, ysmax)

  #plotting
  p <- ggplot(x, aes(x=daf, y=alpha)) + geom_point(size=3) +
    #function curve
    stat_function(fun=y1, color="red") +
    #asymptotic alpha
    geom_hline(yintercept=alpha_asymptotic, color="red", linetype="dashed") +  
    #alpha derived via classic MKT
    geom_hline(yintercept=alpha_original, color="blue", linetype="dashed") +
    #cut-offs
    geom_vline(xintercept=c(xlow, xhigh), color="black", linetype="dashed") +    
    #confidence intervals
    geom_rect(data=data.frame(xmin=-Inf, xmax=Inf, ymin=CI_low ,ymax=CI_high), aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="gray30", alpha=0.4, inherit.aes=F) +
    #shade the fraction of WDMs
    geom_ribbon(data=shader_df, aes(x=xs, ymin=ysmin, ymax=ysmax), fill="gray30", alpha=0.2,inherit.aes=F) + 
    #customization
    theme_classic() +
    theme(axis.text=element_text(size=rel(1.5))) +
    theme(axis.title=element_text(size=rel(1.5))) +
    xlab("derived allele frequency") + ylab("Adaptation (alpha)") +
    #alphas labels
    annotate("text", x=xhigh-0.3, y=alpha_asymptotic-0.2, label=paste0("alpha asymptotic = ", alpha_asymptotic), family="", color="red", size=6) +
    annotate("text",x=xhigh-0.3, y=alpha_original-0.1, label=paste0("alpha original = ",alpha_original), family="", color="blue", size=6)
  
  return(p)
}
