##filter SNV based on alllele frequency and supporting alternative allele
#usage: Rscript strelka_snv_filter.R filelist

library(VariantAnnotation)
library(tidyverse)
library(tidyr)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
snv_bcftools_path=/projects/SNV/snv_bcftools
AF_threshold=0.05                        #set cutoff for allele frequency
alt_count=4                              #set cutoff for supporting alternative allele 
filelist<-args[1]
x<-paste0(snv_bcftools_path,"/",filelist)
out<-readVcf(x)
Sample<-gsub(paste0(snv_bcftools_path,"/"),"",x)
Sample<-gsub(".gz","",Sample)
T_DP<-as.vector(geno(out)$DP[,"TUMOR"])
a<-geno(out)$AU[,,1]
b<-geno(out)$CU[,,1]
c<-geno(out)$GU[,,1]
d<-geno(out)$TU[,,1]
all<-cbind(T_DP,a,b,c,d)
all<-all[,-c(2,4,6,8)]
all<-as.data.frame(all) 
all<-rownames_to_column(all, var = "name")
colnames(all)<-c("name","T_DP","TUMOR_A","TUMOR_C","TUMOR_G","TUMOR_T") 
lat<-substr(all$name, nchar(all$name)-2,nchar(all$name)) 
lat1<-substr(all$name, nchar(all$name)-2,nchar(all$name)-2) 
lat2<-substr(all$name, nchar(all$name),nchar(all$name)) 
library(stringi)
lat1 <- stri_replace_all_regex(lat1,
                                  pattern=c('A','C','G','T'),
                                  replacement=c(3,4,5,6),
                                  vectorize=FALSE)
lat2 <- stri_replace_all_regex(lat2,
                                  pattern=c('A','C','G','T'),
                                  replacement=c(3,4,5,6),
                                  vectorize=FALSE)                                  
all$lat1<-lat1
all$lat2<-lat2
all$lat1<-as.numeric(all$lat1)
all$lat2<-as.numeric(all$lat2)
AF<-vector(mode = "list", length = nrow(all))
 for (j in 1:nrow(all))
 {
 	tem<-as.vector(all[j,])
 	ref<-as.numeric(tem[7])
 	alt<-as.numeric(tem[8])
 	if(as.numeric(tem[alt])>=alt_count)       #remove alternative allele with less than alt_count reads
    {
 	     AF[[j]]<-as.numeric(tem[alt])/(as.numeric(tem[ref])+as.numeric(tem[alt]))
 	  }
  else
    {
    AF[[j]]<-0
    }
 }
all$AF<-unlist(AF)
all1<-subset(all,AF>=AF_threshold)
all1$sample<-Sample
all1<-all1 %>% separate(name,sep=":",c("chr","pos")) %>% separate(pos,sep="_",c("pos","alt")) 
all1<-all1[,c("sample","chr","pos","alt","AF")]
keepchr<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")
all1<-all1[which(all1$chr%in%keepchr),]
write.table(all1,paste0(snv_bcftools_path,"/",Sample,".snv"),quote=F,row.names=F)

