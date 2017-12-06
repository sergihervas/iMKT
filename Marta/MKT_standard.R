# Load genes
genes <- read.table(file="../data_RAL_sim_20DAF/genes",header = F) # List of genes RAL, DAF20, D. Simulans as outgroup
genes <- as.character(genes$V1)

# Create output dataframe
output_1 <- data.frame(gene = character(), Pi = numeric(0), P0 = numeric(0), D0f = numeric(0), D4f = numeric(0),alpha = numeric(0), pvalue = numeric(0), stringsAsFactors=FALSE)

  for (gene in genes) {
    # Load information for each gene
    daf <- read.table(paste0("../data_RAL_sim_20DAF/FORMAT/",gene,"/daf"), header=T)
    divergence <- read.table(paste0("../data_RAL_sim_20DAF/FORMAT/",gene,"/divergence"), header=T)
    
    Pi <- sum(daf$pN)
    P0 <- sum(daf$pS)
    
    if (P0 > 0 & divergence$D0f > 0) {
      #Create MKT table 
      mkt_table <- data.frame(Polymorphism = c(P0, Pi), Divergence=c(divergence$D4f,divergence$D0f),row.names = c("Neutral class","Selected class"))
      
      # Estimation of alpha
      alpha <- 1-(mkt_table[2,1]/mkt_table[1,1])*(mkt_table[1,2]/mkt_table[2,2])
      
      # Fisher's exact test p-value from the MKT
      pvalue <- fisher.test(mkt_table)$p.value
      
      gene_result_cutoff <- data.frame(as.character(gene),as.numeric(Pi),as.numeric(P0),as.numeric(divergence$D0f),as.numeric(divergence$D4f),as.numeric(alpha), as.numeric(pvalue),stringsAsFactors=FALSE)
      
      output_1 <- rbind(output_1,gene_result_cutoff)
    }
}
colnames(output_1) <- c("gene", "Pi", "P0", "D0f", "D4f", "alpha", "pvalue")

## Keep only significant & positive genes ##
output_2 <- subset(output_1, output_1$pvalue < 0.05 & output_1$alpha > 0)

output_1$test <- "MKT"
output_2$test <- "MKT"
mean <- mean(output_2$alpha)
ggplot(output_2, aes(y=alpha,x=test)) +
  geom_boxplot(fill="#386cb0", alpha=0.7) +
  theme_Publication() +
  xlab("McDonald & Kreitman standard test") + ylab("α") +
  annotate("text", x=Inf, y = Inf, label=paste0('alpha [standard] == ', round(mean,digits = 3)), parse=T, size=4, vjust=1, hjust=1) + theme(axis.text.x = element_blank())
