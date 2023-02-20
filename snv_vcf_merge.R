##Read multiple SNV VCF files from a directory and consolidate them into a single file
vcf_dir<-c() 
vcf_suffix<-c() # for example "_S_somatic.snvs.vcf.gz$"

extract_sub_vcf<-function(dir,suffix){
vcf_dir<-dir
suffix<-vcf_suffix
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
