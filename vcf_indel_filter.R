##filter Indel based on alllele frequency and supporting alternative allele
#usage: Rscript strelka_indel_filter.R filelist

library(VariantAnnotation)
library(tidyverse)
library(tidyr)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
AF_threshold=0.05                       #set cutoff for allele frequency
alt_count=4                             #set cutoff for supporting alternative allele 
indel_bcftools_path=/projects/Indel/Indel_bcftools
filelist<-args[1]
x<-paste0(indel_bcftools_path,"/",filelist)
out<-readVcf(x)
Sample<-gsub(paste0(indel_bcftools_path,"/"),"",x)
Sample<-gsub(".gz","",Sample)
T_DP<-as.vector(geno(out)$DP[,"TUMOR"])
a<-geno(out)$TAR[,,1]                   #extract ref allele
b<-geno(out)$TIR[,,1]                   #extract alternative allele
all<-cbind(T_DP,a,b)
all<-all[,-c(2,4)]
all<-as.data.frame(all) 
all<-rownames_to_column(all, var = "name")
colnames(all)<-c("name","T_DP","TUMOR_REF","TUMOR_ALT") 
all$AF<-all$TUMOR_ALT/(all$TUMOR_REF+all$TUMOR_ALT)
all<-subset(all,TUMOR_ALT>=alt_count)
all1<-subset(all,AF>=AF_threshold)
all1$sample<-Sample
all1<-all1 %>% separate(name,sep=":",c("chr","pos")) %>% separate(pos,sep="_",c("pos","alt")) 
all1<-all1[,c("sample","chr","pos","alt","AF")]
keepchr<-paste0("chr",c(seq(1:22),"X"))
all1<-all1[which(all1$chr%in%keepchr),]
write.table(all1,paste0(indel_bcftools_path,"/",Sample,".indel"),quote=F,row.names=F)
