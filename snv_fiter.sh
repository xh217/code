file_name=xx
pbs_path=/projects/Strelka_filter 
strelka_path=/projects/strelka  
snv_bcftools_path=/projects/SNV/snv_bcftools 
annovar_path=/projects/SNV/annovar
ref_genome_version=hg38
annovar_ref_path=/home/annovar/humandb
common_snv_path=/projects/common_snv

##keep PASS,QSS>=40, normal and tumor read depth >=20   
bcftools filter  -i 'FILTER=="PASS" && QSS >= 40 && FORMAT/DP[0] >= 10 && FORMAT/DP[1] >= 20' $strelka_path/$file_name/results/variants/somatic.snvs.vcf.gz --output $snv_bcftools_path/$file_name.gz --output-type z

##filter SNV AF<=0.05 and alternative allele with less than 4 reads 
Rscript --vanilla $pbs_path/strelka_snv_filter.R $file_name.gz 

##annovar filteration,fitering >1% in dbSNPv150, 1000 genomes (The 1000 Genomes Project Consortium 2015), segmental duplications (genomicSuperDups) and centromeres 

#annotate_variation.pl -buildver $ref_genome_version -downdb -webfrom annovar refGene $annovar_ref_path/ 

cat $snv_bcftools_path/$file_name.snv|sed '1d' |awk -F' |/' '{print $2,$3,$3,$4,$5,"unknown",".",100,50}' |awk -F' |/' 'BEGIN {FS=" ";OFS="t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' > $annovar_path/$file_name.snv.annovar  

table_annovar.pl $annovar_path/$file_name.snv.annovar $annovar_ref_path  -remove -buildver $ref_genome_version -out $annovar_path/$file_name -protocol avsnp150,1000g2015aug_all,1000g2015aug_amr,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,refGene -operation f,f,f,f,f,f,f,g -nastring .

awk 'NR==FNR{a[$5]=$0} NR!=FNR{if(!($6 in a)) {print $0} }'  $common_snv_path/snp150Common.txt    $file_name.hg38_multianno.txt  > $file_name.snp.snv

cat $file_name.snp.snv |awk '{if($7<0.01 && $8<0.01 && $9<0.01 && $10<0.01 && $11<0.01 && $12<0.01)print $0}' > $file_name.snv.final
