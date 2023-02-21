DESeqplot_annogene<-function(dataframe,log2FC_cutoff,padj_cutoff,...)
{
res1<-dataframe %>%
mutate(gene_type = case_when(log2FoldChange >=2  & padj <= 0.05 ~ "up",
                               log2FoldChange <= -2 & padj <= 0.05 ~ "down",
                               TRUE ~ "ns")) 
cols <- c("down" = '#3E4A89FF', "up" = '#E69F00', "ns" = "grey") 
sizes <- c("up" = 3, "down" = 3, "ns" = 1) 
alphas <- c("up" = 1, "down" = 1, "ns" = 0.3)
sig_genes <-  res1[which(res1$gene_type=="down"|res1$gene_type=="up"),] 
up_genes <-  res1[which(res1$gene_type=="up"),]
down_genes <-  res1[which(res1$gene_type=="down"),]
png(paste0("/projects/Deseq","/DEseq.png"),width = 6,height = 4, units = "in", res = 400)  
ggplot(res1,aes(x=log2FoldChange,
             y = -log(padj),
             fill = gene_type,    
             size = gene_type,
             alpha = gene_type)) + 
   geom_point(shape = 21,    
             colour = "black") + 
   ggtitle('Valcano')+           
   geom_hline(yintercept = -log(0.05),
             linetype = "dashed") + 
   geom_vline(xintercept = c(-2, 2),
             linetype = "dashed") +
   scale_fill_manual(values = cols) + # Modify point colour
   scale_size_manual(values = sizes) + # Modify point size
   scale_alpha_manual(values = alphas)+
   theme_prism(base_size = 13)+
   theme(legend.position = "none")+
   geom_label_repel(data = sig_genes, # Add labels last to appear as the top layer  
                   aes(label = rownames(sig_genes)),
                   force = 2,
                   cex=5,
                   color="white",
                   label.size = NA,
                   nudge_y = 1)
dev.off() 
}
