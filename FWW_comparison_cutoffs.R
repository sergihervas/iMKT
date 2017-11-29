# Genes
genes <- read.table(file="../data_RAL_sim_20DAF/genes",header = F)
genes <- as.character(genes$V1)

out <- data.frame(gene = character(), Pi = numeric(0), P0 = numeric(0), D0f = numeric(0), D4f = numeric(0),alpha = numeric(0), pvalue = numeric(0), stringsAsFactors=FALSE)

list_cutoffs <- c(0, 0.025, 0.05, 0.075, 0.125, 0.175)

for (cutoff in list_cutoffs) {
  
    for (gene in genes) {
    
      daf <- read.table(paste0("../data_RAL_sim_20DAF/FORMAT/",gene,"/daf"), header=T)
      divergence <- read.table(paste0("../data_RAL_sim_20DAF/FORMAT/",gene,"/divergence"), header=T)
      
      Pi <- sum(daf[daf$daf>cutoff,]$pN)
      P0 <- sum(daf[daf$daf>cutoff,]$pS)
      
      if (P0 > 0 & divergence$D0f > 0) {
        
        #Create MKT table 
        mkt_table <- data.frame(Polymorphism = c(P0, Pi), Divergence=c(divergence$D4f,divergence$D0f),row.names = c("Neutral class","Selected class"))
        
        # Estimation of alpha
        alpha <- 1-(mkt_table[2,1]/mkt_table[1,1])*(mkt_table[1,2]/mkt_table[2,2])
        
        # Fisher's exact test p-value from the MKT
        pvalue <- fisher.test(mkt_table)$p.value
        
        g1 <- data.frame(as.character(gene),as.numeric(Pi),as.numeric(P0),as.numeric(divergence$D0f),as.numeric(divergence$D4f),as.numeric(alpha), as.numeric(pvalue),as.numeric(cutoff),stringsAsFactors=FALSE)
        
        out <- rbind(out,g1)
      }
      
    }
}
colnames(out) <- c("gene", "Pi", "P0", "D0f", "D4f", "alpha", "pvalue", "cutoff")

#data <- aggregate(alpha ~ cutoff, data=out, mean)

data <- aggregate(alpha ~ cutoff, data=out, function(x) c(mean = mean(x), length = length(x)))
data_cutoff <- data$cutoff
data<-as.data.frame(data$alpha)
data$cutoff <- data_cutoff

ggplot(data, aes(x=as.factor(cutoff), y = mean, group=1, label=length)) +
  geom_line(color="#386cb0") + 
  geom_point(size=2.5, color="#386cb0")+
  theme_Publication() +
  xlab("Cut-off") + ylab("α") +
  geom_text_repel(aes(label=length))

out2 <- data.frame(gene = character(), Pi = numeric(0), P0 = numeric(0), D0f = numeric(0), D4f = numeric(0),alpha = numeric(0), pvalue = numeric(0), cutoff = numeric(0), cutoff_original = numeric(0), stringsAsFactors=FALSE)


for (cutoff_original in head(rev(list_cutoffs),-1)) {
  gene_list <- subset(out$gene,out$cutoff==cutoff_original)
  
  for (cutoff in list_cutoffs[list_cutoffs<cutoff_original])
  {
    gene_cutoff<- out[ which(out$cutoff==cutoff),]
    
    data <- gene_cutoff[gene_cutoff$gene %in% gene_list,] 
    data$cutoff_original <- cutoff_original
    out2 <- rbind(out2,data)
  }
}

out$cutoff_original <- out$cutof

data_plot <- rbind(out,out2)

#data <- aggregate(alpha ~ cutoff + cutoff_original, data=data_plot, mean)

data <- aggregate(alpha ~ cutoff + cutoff_original, data=data_plot, function(x) c(mean = mean(x), length = length(x)))
data_cutoff <- data$cutoff
data_cutoff_original <- data$cutoff_original
data<-as.data.frame(data$alpha)
data$cutoff <- data_cutoff
data$cutoff_original <- data_cutoff_original
#data <- aggregate(alpha ~ cutoff + cutoff_original, data=data_plot, function(x) c(mean = mean(x), length = length(x)))

ggplot(data, aes(x=as.factor(cutoff), y=mean,color=as.factor(cutoff_original),label=length)) +
  geom_point()+
  geom_line(aes(group=as.factor(cutoff_original))) +
  scale_colour_Publication() +
  theme_Publication() +
  xlab("Cut-off") + ylab("α") + 
  geom_text_repel(aes(label=length))




