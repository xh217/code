input_dir<-"/scratch/xh217/KMT2C_work"
snv_path<-file.path(input_dir,'TCGA_MAF_file_extract')

#create output path
out_dir<-"/scratch/xh217/KMT2C_work/TCGA_anno_case_control"
output<-file.path(out_dir,'counts')
#TCGA maf
snv_file_list=c( 'TCGA.BLCA.somatic.maf',
				         'TCGA.BRCA.somatic.maf',
				         'TCGA.CESC.somatic.maf',
                 'TCGA.COAD.somatic.maf',
                 'TCGA.ESCA.somatic.maf',
                 'TCGA.GBM.somatic.maf',
                 'TCGA.HNSC.somatic.maf',        
                 'TCGA.KIRC.somatic.maf',
                 'TCGA.KIRP.somatic.maf',
                 'TCGA.LIHC.somatic.maf',
                 'TCGA.LUAD.somatic.maf',
                 'TCGA.LUSC.somatic.maf',
                 'TCGA.OV.somatic.maf',
                 'TCGA.PRAD.somatic.maf',
                 'TCGA.SKCM.somatic.maf',
                 'TCGA.STAD.somatic.maf',
                 'TCGA.UCEC.somatic.maf') 

total=c()
snv_df <- file.path(snv_path,snv_file_list[m])
test<-read.table(text = gsub(" ", "\t", readLines(snv_df)),header=T,stringsAsFactors=FALSE)
func_loci=c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","Nonstop_Mutation","Translation_Start_Site","In_Frame_Del","In_Frame_Ins")
sample_name<-test[which(test$Hugo_Symbol=="KMT2C" 
                                        & (test$Variant_Classification==func_loci[1]
                                         | test$Variant_Classification==func_loci[2]
                                         | test$Variant_Classification==func_loci[3]
                                         | test$Variant_Classification==func_loci[4]
                                         | test$Variant_Classification==func_loci[5]
                                         | test$Variant_Classification==func_loci[6]
                                         | test$Variant_Classification==func_loci[7]
                                         | test$Variant_Classification==func_loci[8]
                                          )),4]

uniq_name<-unique(sample_name)
uniq_name<-cbind.data.frame(uniq_name,anno=1)
colnames(uniq_name)<-c("Tumor_Sample_Barcode", "KMT2C")
output_name<-strsplit(snv_file_list[m], '.somatic', fixed= TRUE)[[1]][1]
# add one column description of MLL3 status - MLL3 mutant(1), MLL3 wild type(0))
data1<-full_join(test, uniq_name, by = "Tumor_Sample_Barcode")
data1$KMT2C[is.na(data1$KMT2C)] <-0
