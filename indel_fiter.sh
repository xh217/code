file_name=xx
pbs_path=/projects/Strelka_filter 
strelka_path=/projects/strelka  
indel_bcftools_path=/projects/Indel/indel_bcftools 
annovar_path=/projects/Indel/annovar
ref_genome_version=hg38
annovar_ref_path=/home/annovar/humandb
common_indel_path=/projects/common_indel

##keep PASS,QSI>=40, normal and tumor read depth >=20   
bcftools filter  -i 'FILTER=="PASS" && QSI >= 40 && FORMAT/DP[0] >= 10 && FORMAT/DP[1] >= 20' $strelka_path/$file_name/results/variants/somatic.indels.vcf.gz --output $indel_bcftools_path/$file_name.gz --output-type z

##filter Indel AF<=0.05 and alternative allele with less than 4 reads 
Rscript --vanilla $pbs_path/strelka_indel_filter.R $file_name.gz 

##annovar filteration,fitering >1% in dbSNPv150, 1000 genomes (The 1000 Genomes Project Consortium 2015), segmental duplications (genomicSuperDups) and centromeres 
#annotate_variation.pl -buildver $ref_genome_version -downdb -webfrom annovar refGene $annovar_ref_path/ 
cat $indel_bcftools_path/$file_name.indel|sed '1d' |awk -F' |/' '{print $2,$3,$3,$4,$5,"unknown",".",100,50}' |awk -F' |/' 'BEGIN {FS=" ";OFS="t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' > $annovar_path/$file_name.indel.annovar  
table_annovar.pl $annovar_path/$file_name.indel.annovar $annovar_ref_path  -remove -buildver $ref_genome_version -out $annovar_path/$file_name -protocol avsnp150,1000g2015aug_all,1000g2015aug_amr,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,refGene -operation f,f,f,f,f,f,f,g -nastring .
awk 'NR==FNR{a[$5]=$0} NR!=FNR{if(!($6 in a)) {print $0} }'  $common_indel_path/snp150Common.txt    $file_name.hg38_multianno.txt  > $file_name.indel
cat $file_name.indel |awk '{if($7<0.01 && $8<0.01 && $9<0.01 && $10<0.01 && $11<0.01 && $12<0.01)print $0}' > $file_name.indel.final
