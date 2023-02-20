##Read multiple SNV VCF files from a directory and consolidate them into a single file, and plot mutational signature burden for each sample
vcf_dir<-c() 
vcf_suffix<-c() # for example "_S_somatic.snvs.vcf.gz$"

library(VariantAnnotation)
library(dplyr)
library(tidyr)
filenames <- list.files(vcf_dir, full.names = T, pattern = vcf_suffix, recursive = TRUE)
my_read_csv <- function(x) 
   {
      out<-readVcf(x)
        if(dim(out)==0)
          {
  	        df=data.frame(
  	        chr=as.character(0),
  	        start=as.numeric(0),
  	        ref=as.character(0),
  	        alt=as.character(0),
  	        sample=as.character(0)
  	        )
          }
        else
         {
          df<-data.frame(info(out))
          df<-data.frame(rownames(df))
          df<-df %>% separate(colnames(df), c("chr","start","ref","alt"), "[:_/]")
          sample<-gsub(paste0(vcf_dir,"/"),"",x)
          sample<-gsub(vcf_suffix,"",sample)
          df$sample<-sample
         }
     return(df)
   }
tbl_vcf <- lapply(filenames, my_read_csv)
Matrix_vcf = do.call(rbind, tbl_vcf)  
Matrix_vcf<-as.data.frame(Matrix_vcf) #chromosome position should be as numeric vectors
n2.input <- mut.to.sigs.input(mut.ref = Matrix_vcf, sample.id = "sample", chr = "Chr", pos = "Start", ref = "Ref", alt = "Alt", bsg=BSgenome.Hsapiens.UCSC.hg38)
#set the output sample order in your project, case and control visulization 
#order_event<-c("S4","S8","S10","S14","S16","S21","S25","S29","S50","S52","S56","S63","S65","S19","S30","S34","S53","S61","S3","S44","S17","S42","S48","S28","S32","S18","S49","S7","S45","S55","S39","S66","S58","S13","S9","S11","S43","S22","S31","S36")
n2.input<-n2.input[rownames(n2.input)%in%order_event,]
n2.input1<-n2.input[order(match(rownames(n2.input),order_event)),]
fit_res <- fit_to_signatures(t(n2.input1), signatures=t(signatures.cosmic))
fit_res_retain<-fit_res$contribution[rowSums(fit_res$contribution)>5000,] # Set a cutoff to exclude mutational signatures and make the plot less crowded
Others<-colSums(fit_res$contribution)-colSums(fit_res_retain)
fit_res_all<-rbind(fit_res_retain,Others)
#assign specific color to highlight particular mutational signatures
col1<-c("orange","yellow")
col2<-"red"
colfunc <- colorRampPalette(c("lightblue","blue"))
col3<-colfunc(3)
col4<-"purple"
colfunc <- colorRampPalette(c("lightgreen","white"))
col5<-colfunc(7)
plot_contribution(fit_res_all,
  coord_flip = T,
  mode = "relative",
  palette=c(col1,col2,col3,col4,col5))
  theme_prism(base_size = 5) # optimize ggplot to improve visulization
  guides(fill=guide_legend(ncol=2),size =5)


   
   
