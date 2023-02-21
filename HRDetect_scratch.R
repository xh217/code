library(VariantAnnotation)
library(signature.tools.lib)

##set the file names and file directory 
sample_names <- c()             #sample order
SNV_tab_files<- list.files()    #SNV files directory with file suffix
SV_bedpe_files <- list.files()  #SV files directory with file suffix
CNV_files <- c()                #cnv files directory
Indels_vcf_files <- c()         #indel files directory


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
expected_chroms <- paste0("chr",c(seq(1:22),"X","Y"))
Indelscat_list <- list()
for (i in 1:length(Indels_vcf_files)){
tmp<-read.table(Indels_vcf_files[i],header=F,sep='\t',stringsAsFactors=F)
 colnames(tmp)<-c("chr","position","position1","REF","ALT")
 tmp$chr<-gsub("chr","",tmp$chr)
  res <- tabToIndelsClassification(tmp,sample_names[i],genome.v = "hg38")
  Indelscat_list[[i]] <- res$count_proportion
}
#bind the catalogues in one table
Indel_catalogues <- do.call(rbind,Indelscat_list)
rownames(Indel_catalogues) <- sample_names
input_matrix[colnames(Indel_catalogues),"del.mh.prop"] <- Indel_catalogues[,"del.mh.prop"]
#bind all catalogues
SNV_catalogues <- do.call(cbind,SNVcat_list)
#set snv input matrix
col_hrdetect <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
input_matrix <- matrix(NA,nrow = length(sample_names),
                       ncol = length(col_hrdetect),
                       dimnames = list(sample_names,col_hrdetect))

#fit the 12 breast cancer signatures using the bootstrap signature fit approach
sigsToUse <- c(1,2,3,5,6,8,13,17,18,20,26,30)      #used if for breast cancer
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

##prep indel input
Indels_vcf_files <- list.files(Indels_vcf_files, full.names = T, pattern = "\\.indel$", recursive = F)
expected_chroms <- paste0("chr",c(seq(1:22),"X"))
Indelscat_list <- list()
for (i in 1:length(Indels_vcf_files))
{
tmp<-read.table(Indels_vcf_files[i],header=F,sep='\t',stringsAsFactors=F)
  colnames(tmp)<-c("chr","position","position1","REF","ALT")
  tmp$chr<-gsub("chr","",tmp$chr)
  res <- tabToIndelsClassification(tmp,sample_names[i],genome.v = "hg38")
  Indelscat_list[[i]] <- res$count_proportion
}
#bind the catalogues into SNV input_matrix
Indel_catalogues <- do.call(rbind,Indelscat_list)
rownames(Indel_catalogues) <- sample_names
input_matrix[colnames(Indel_catalogues),"del.mh.prop"] <- Indel_catalogues[,"del.mh.prop"]

#prep cnv input
output<-c()
cnv_filenames <- list.files(CNV_files, full.names = F, pattern = "\\.cnv", recursive = TRUE)
for(i in 1:(length(cnv_filenames)))
{
cnv_df <- cnv_filenames[i]
file_a<-read.table(cnv_df,header=F, sep='\t',stringsAsFactors=F)
file_a_extract<-file_a[,c("V1","V2")]
file_a_extract$V1<-gsub("chr","",file_a_extract$V1)
cnv_vcf <- readVcf(cnv_df, "hg38")
cnv_vcf_extract<-data.frame(info(vcf1))[,c("END","TCN_EM","LCN_EM")]
combined<-cbind(cnv_vcf,cnv_vcf_extract)
combined<-combined[!(combined$TCN_EM == 2 & combined$LCN_EM == 1) & !is.na(combined$TCN_EM) & !is.na(combined$LCN_EM),]
combined$total.copy.number.inNormal<-2
combined$minor.copy.number.inNormal<-1
combined$seg_no<-1:nrow(combined) 
#combined$V1<-gsub("chr","",combined$V1)
colnames(combined)[1:5]<-c("Chromosome","chromStart","chromEnd","total.copy.number.inTumour","minor.copy.number.inTumour")
combined_extract<-combined[,c('seg_no', 'Chromosome', 'chromStart', 'chromEnd', 'total.copy.number.inNormal', 'minor.copy.number.inNormal', 'total.copy.number.inTumour', 'minor.copy.number.inTumour')]
write.table(combined_extract, paste0(output,cnv_df),sep='\t',quote=F)
}

#prep SV input
SV_dir=/projects/SV_dir
for i in `ls *_convert_bedpe`;
do
cat $i | grep -v "#"| awk '{print $1,$2,$3,$4,$5,$6,$9,$10,$23,$7}' |awk  -F' |:|,'  'BEGIN {OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$10+$11,$14}' > $SV_dir/$i;
done
for i in `ls *_bedpe`; do 
sed -i "s/$/\t$i/" $i;
sed -i "s/chr//g" $i; 
awk 'BEGIN{print "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tstrand1\tstrand2\tpe_support\tsvclass\tsample"}1' $i > ${i}1;
done
for i in `ls *_bedpe1`; do 
sed -i 's/Manta//g' $i;
sed -i "s/DEL/deletion/g" $i;
sed -i "s/BND/translocation/g" $i;
sed -i "s/INV/inversion/g" $i;
sed -i "s/DUP/tandem-duplication/g" $i;
done

#run the HRDetect pipeline
res <- HRDetect_pipeline(input_matrix,
                         genome.v = "hg38",
                         SV_bedpe_files = SV_bedpe_files,
                         #Indels_vcf_files = Indels_vcf_files, # indel matrix has been incorporated into SNV input_matrix
                         CNV_tab_files = CNV_tab_files,
                         nparallel = 2)
