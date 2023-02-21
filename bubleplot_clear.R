bubleplot_clear<-function(buble_df)
{
#set the SBS order for visualization           
buble_df$variable=factor(buble_df$variable,levels=rev(c("SBS1","SBS3","SBS4","SBS5","SBS6","SBS7","SBS8","SBS9","SBS10","SBS11","SBS12","SBS14","SBS15","SBS16","SBS17","SBS18","SBS19","SBS20","SBS21","SBS22","SBS23","SBS24","SBS25","SBS26","SBS27","SBS28","SBS29","SBS30")))
ggplot(buble_df, aes(x=Tumor_Sample_Barcode, y=variable)) + 
       geom_point(aes(size=value), shape=21, color="white", fill="#4B0082") + 
       scale_size_area(max_size=5, guide=FALSE) + 
       theme(axis.text.x = element_blank(),axis.text.y = element_text(hjust=1,size=9),axis.ticks.x=element_blank()) +
       theme(axis.title.x = element_blank(),axis.title.y = element_blank())
}            
