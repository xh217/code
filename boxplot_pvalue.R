boxplot_clear<-function(dataframe,y_p_pos,title_text,ylab_text...)
{
result <- wilcox.test(dataframe[,2]~ dataframe[,1], data = dataframe)$p.value
p_value <- signif(result, digits = 3)
title_text<-title_text
ylab_text<-ylab_text
p_vals <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position,
  unique(dataframe[,1])[1],   unique(dataframe[,1])[2],     p_vals, 1) 
ggplot(RCB_hrd,aes(x = dataframe[,1], y = dataframe[,2])) + 
  geom_violin(aes(colour = dataframe[,1], fill = dataframe[,1]), trim = FALSE) + 
  geom_boxplot(aes(fill = dataframe[,1]), width = 0.2, colour = "black") +
  scale_color_manual(values=c('#E69F00',"#3E4A89FF"))+
  scale_fill_manual(values=c('#E69F00',"#3E4A89FF"))+
  scale_y_continuous(limits = c(min(dataframe[,2]), max(dataframe[,2])))  +                         
  ggtitle(title)+
  xlab('') +
  ylab(ylab_text)+
  guides(y = "prism_minor") + 
  theme_prism(base_size = 13) + 
  theme(legend.position = "none") + 
  add_pvalue(p_vals, label = "p = {p.adj}", 
                     remove.bracket = TRUE, x = 1.5,label.size=5,y_p_pos)
 }                   
