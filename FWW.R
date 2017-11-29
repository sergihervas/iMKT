#################################################################
################# MKT corrected with FWW method #################
#################################################################

# Date = 28/11/2016
# Author = Sergi Hervás, Marta Coronado

# The standard McDonald and Kreitman test (MKT) is used to detect the signature of selection at the molecular level. The MKT compares the amount of variation within a species (polymorphism, P) to the divergence (D) between species at two types of sites, one of which is putatively netral and used as the reference to detect selection at the other type of site. In the standard MKT, these sites are synonymous (putatively neutral, 0) and non-synonymous sites (selected sites, i) in a coding region. Under strict neutrality, the ratio of the number of selected and neutral polymorphic sites (Pi/P0) is equal to the ratio of the number of selected and neutral divergence sites (Di/D0).
# The null hypothesis of neutrality is rejected in a MKT when Di/D0 > Pi/P0. The excess of divergence relative to polymorphism for class i, is interpreted as adaptive selection for a subset of sites i. The fraction of adaptive fixations, α, is estimated from 1-(PN/PS)(Ds/Dn). The significance of the test can be assesed with a Fisher exact test.
# The estimate of α can be easily biased by the segregation of slightly deleterious non-synonymous substitutions. Specifically, slightly deleterious mutations tend to contribute more to polymorphism than to divergence, and thus, lead to an underestimation of alpha. Bevause they tend to segregate at lower frequencies than do neutral mutations, they can be apartially controled for by removing low frequency polymorphisms from the analysis (Fay et al. 2001). This is known as the FWW method.

####################################################
################# MKT-FWW function #################
####################################################

mkt_fww <- function(daf = "Data frame containing the DAF, Pn and Ps", 
                     divergence = "Data frame that contains the number of non-synonymous sites analyzed (m0f), number of synonymous sites analyzed (m4f), total number of non-synonymous divergent sites (D0f) and number of synonymous divergent sites (D4f)") {
  # Shows a message when using the function
  packageStartupMessage("MKT corrected by the FWW method")
  
  # Declare output data frame
  output <- data.frame(cutoff = numeric(0), alpha = numeric(0), pvalue = integer(0))
  #out <- NULL
  
  list_cutoffs <- c(0, 0.05, 0.1)
  
  for (cutoff in list_cutoffs) {
  # Remove daf below a frequency of 0.05
  daf_remove <-daf[daf$daf > cutoff,] 
  
  #Create MKT table 
  mkt_table <- data.frame(Polymorphism = c(sum(daf_remove$pS), sum(daf_remove$pN)), Divergence=c(divergence$D4f,divergence$D0f),row.names = c("Neutral class","Selected class"))

  # Estimation of alpha
  alpha <- 1-(mkt_table[2,1]/mkt_table[1,1])*(mkt_table[1,2]/mkt_table[2,2])
  
  # Fisher's exact test p-value from the MKT
  pvalue <- fisher.test(mkt_table)$p.value
  
  # Store output  
  
  output_n1 <- cbind(cutoff,alpha,pvalue)
  output <- rbind(output,output_n1)
  

  }
  output <- as.data.frame(output)
  names(output) <- c("Cut-off","α","Fisher's exact test P-value")
   return(output)
}
