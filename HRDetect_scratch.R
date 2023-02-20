library(VariantAnnotation)
library(signature.tools.lib)

##set the file names and file directory 
sample_names <- c()
SNV_tab_files<- list.files()
SV_bedpe_files <- list.files()
CNV_tab_files <- list.files()
Indels_vcf_files <- list.files()
SNV_output_dir<-c()
Indel_output_dir<-c()

##name the vectors with the sample names
names(SNV_tab_files)<-sample_names
names(SV_bedpe_files) <- sample_names
names(Indels_vcf_files) <- sample_names
names(CNV_tab_files) <- sample_names

##prep snv input
#calculate snv catalogy
SNVcat_list <- list()
for (i in 1:length(SNV_tab_files))
{
 tmp<-read.table(SNV_tab_files[i],header=T,sep='\t',stringsAsFactors=F)
 colnames(tmp)<-c("chr","position","position1","REF","ALT")
 tmp$chr<-gsub("chr","",tmp$chr)
 res <- tabToSNVcatalogue(subs = tmp,genome.v = "hg38")
 colnames(res$catalogue) <- sample_names[i]
 SNVcat_list[[i]] <- res$catalogue
}
#bind all catalogues
SNV_catalogues <- do.call(cbind,SNVcat_list)
#set snv input matrix
col_hrdetect <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
input_matrix <- matrix(NA,nrow = length(sample_names),
                       ncol = length(col_hrdetect),
                       dimnames = list(sample_names,col_hrdetect))

#fit the 12 breast cancer signatures using the bootstrap signature fit approach
sigsToUse <- c(1,2,3,5,6,8,13,17,18,20,26,30)
subs_fit_res <- SignatureFit_withBootstrap_Analysis(outdir = SNV_output_dir,
                                    cat = SNV_catalogues,
                                    signature_data_matrix = COSMIC30_subs_signatures[,sigsToUse],
                                    type_of_mutations = "subs",
                                    nboot = 1000,nparallel = 4)
#The signature exposures can be found here and correspond to the median of the boostrapped runs followed by false positive filters
snv_exp <- subs_fit_res$E_median_filtered
input_matrix[colnames(snv_exp),"SNV3"] <- snv_exp["Signature.3",]
input_matrix[colnames(snv_exp),"SNV8"] <- snv_exp["Signature.8",]

##prep indel input
expected_chroms <- paste0("chr",c(seq(1:22),"X"))
Indelscat_list <- list()
for (i in 1:length(Indels_vcf_files))
{
  tmp<-read.table(Indels_vcf_files[i],header=T,sep='\t',stringsAsFactors=F)
  colnames(tmp)<-c("chr","position","position1","REF","ALT")
  tmp$chr<-gsub("chr","",tmp$chr)
  res <- tabToIndelsClassification(tmp,sample_names[i],genome.v = "hg38")
  Indelscat_list[[i]] <- res$count_proportion
}
#bind the catalogues in one table
Indel_catalogues <- do.call(rbind,Indelscat_list)
rownames(Indel_catalogues) <- sample_names
input_matrix[rownames(Indel_catalogues),"del.mh.prop"] <- Indel_catalogues[,"del.mh.prop"]


#run the HRDetect pipeline
res <- HRDetect_pipeline(input_matrix,
                         genome.v = "hg38",
                         SV_bedpe_files = SV_bedpe_files,
                         Indels_vcf_files = Indels_vcf_files,
                         CNV_tab_files = CNV_tab_files,
                         nparallel = 2)
