# Load genes
genes <- read.table(file="genes_asym",header = F) # List of genes RAL, DAF20, D. Simulans as outgroup
genes <- as.character(genes$V1)
genes<-genes[1]
# Load genes
genes <- read.table(file="../data_RAL_sim_20DAF/genes",header = F) # List of genes RAL, DAF20, D. Simulans as outgroup
genes <- as.character(genes$V1)



# Create output dataframe
output_1 <- data.frame(gene = character(), Pi = numeric(0), P0 = numeric(0), D0f = numeric(0), D4f = numeric(0),alpha = numeric(0), pvalue = numeric(0), stringsAsFactors=FALSE)

list_cutoffs <- c(0, 0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325)

for (cutoff in list_cutoffs) {
    for (gene in genes) {
      # Load information for each gene
      daf <- read.table(paste0("../data_RAL_sim_20DAF/FORMAT/",gene,"/daf"), header=T)
      divergence <- read.table(paste0("../data_RAL_sim_20DAF/FORMAT/",gene,"/divergence"), header=T)
      
      daf_below_cutoff <- daf[daf$daf <= cutoff, ] #below DAF
      daf_above_cutoff <- daf[daf$daf > cutoff, ] #over DAF
      
      P0 <- sum(daf$pS) 
      Pi <- sum(daf$pN) 
      
      if (P0 > 0 & divergence$D0f > 0) {
        #Create MKT table 
        mkt_table_standard <- data.frame(Polymorphism = c(sum(daf$pS), sum(daf$pN)), Divergence=c(divergence$D4f,divergence$D0f),row.names = c("Neutral class","Selected class"))
        
        mkt_table <- data.frame(`DAF below cutoff` = c(sum(daf_below_cutoff$pS), sum(daf_below_cutoff$pN)), `DAF above cutoff`=c(sum(daf_above_cutoff$pS), sum(daf_above_cutoff$pN)),row.names = c("Neutral class","Selected class"))
        
        f_neutral <- mkt_table[1,1]/sum(daf$pS)
        Pi_neutral_below_cutoff <- Pi * f_neutral
        
        Pi_neutral <- round(Pi_neutral_below_cutoff + mkt_table[2,2])
        
        # Estimation of alpha
        alpha <- 1-((Pi_neutral/P0)*(mkt_table_standard[1,2]/mkt_table_standard[2,2]))
        
        # Fisher's exact test p-value from the MKT
        m <- matrix(c(P0,Pi_neutral,divergence$D4f,divergence$D0f), ncol=2)
        pvalue <- fisher.test(m)$p.value
        
        gene_result_cutoff <- data.frame(as.character(gene),as.numeric(Pi),as.numeric(P0),as.numeric(divergence$D0f),as.numeric(divergence$D4f),as.numeric(alpha), as.numeric(pvalue),as.numeric(cutoff),stringsAsFactors=FALSE)
        
        output_1 <- rbind(output_1,gene_result_cutoff)
      }
      
    }
}
colnames(output_1) <- c("gene", "Pi", "P0", "D0f", "D4f", "alpha", "pvalue", "cutoff")

## Keep only significant & positive genes ##
output_1 <- subset(output_1, output_1$pvalue < 0.05 & output_1$alpha > 0)

## Aggregate data: akpha mean and number of genes analyzed in each cutoff
data <- aggregate(alpha ~ cutoff, data=output_1, function(x) c(mean = mean(x), length = length(x), sd = sd(x)))
data_cutoff <- data$cutoff
data<-as.data.frame(data$alpha)
data$cutoff <- data_cutoff

ggplot(data, aes(x=as.factor(cutoff), y = mean, group=1, label=length)) +
  geom_line(color="#386cb0") + 
  geom_point(size=2.5, color="#386cb0")+
  theme_Publication() +
  xlab("Cut-off") + ylab("α") +
  geom_text_repel(aes(label=length))

######## Check genes analyzed with a cutoff with a lower cut-off
output_2 <- data.frame(gene = character(), Pi = numeric(0), P0 = numeric(0), D0f = numeric(0), D4f = numeric(0),alpha = numeric(0), pvalue = numeric(0), cutoff = numeric(0), cutoff_original = numeric(0), stringsAsFactors=FALSE)


for (cutoff_original in head(rev(list_cutoffs),-1)) {
  gene_list <- subset(output_1$gene,output_1$cutoff==cutoff_original)
  
  for (cutoff in list_cutoffs[list_cutoffs<cutoff_original])
  {
    gene_cutoff<- output_1[ which(output_1$cutoff==cutoff),]
    
    data <- gene_cutoff[gene_cutoff$gene %in% gene_list,] 
    data$cutoff_original <- cutoff_original
    output_2 <- rbind(output_2,data)
  }
}
## Keep only significant genes
output_2 <- subset(output_2, output_2$pvalue < 0.05 & output_2$alpha > 0)

output_1$cutoff_original <- output_1$cutof

# Merge with previpus result
data_plot <- rbind(output_1,output_2)

## Aggregate data: akpha mean and number of genes analyzed in each cutoff
data <- aggregate(alpha ~ cutoff + cutoff_original, data=data_plot, function(x) c(mean = mean(x), length = length(x),sd=sd(x)))

data_cutoff <- data$cutoff
data_cutoff_original <- data$cutoff_original
data<-as.data.frame(data$alpha)
data$cutoff <- data_cutoff
data$cutoff_original <- data_cutoff_original

data_labels <- data[data$cutoff==data$cutoff_original,]

ggplot(data, aes(x=as.factor(cutoff), y=mean,color=as.factor(cutoff_original),label=length)) +
  geom_text_repel(aes(label=length),data = data_labels,nudge_x = c(0.2,0.2),  segment.color = 'white') +
  geom_point()+
  geom_line(aes(group=as.factor(cutoff_original))) +
  scale_colour_Publication(name = "Original cut-off") +
  theme_Publication() +
  xlab("Cut-off") + ylab("α") 
