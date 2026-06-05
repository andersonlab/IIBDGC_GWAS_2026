# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

#  files downloaded here:
# /path/to/project


# All requested files are downloaded to
# /path/to/project You can find the files in two
# folders: GCST90301001-GCST90302000 and GCST90302001-GCST90303000.

# list of all studies in GWAS Catalog explore and select those appropiate

# singularity exec iibdgc_postprocess_10_singularity.sif
# singularity exec polyfun_2_singularity.sif

path_gwas="/path/to/ibdgwas/IIBDGC/"

# reference file downloaded from:
# https://www.ebi.ac.uk/gwas/docs/file-downloads 

# gzip /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/gwas-catalog-v1.0.3.1-studies-r2025-05-13.tsv


###########################################
# SELECT FILES TO DOWNLOAD AMONG THOSE WITH RELATED PHENOTYPES, CYTOKINES, 

MEM=1000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R 

library(data.table)
library(ggplot2)
library(ggpubr)
library(colorspace)
library(rcartocolor)
library(qvalue)
library(stringr)

path_gwas="/path/to/ibdgwas/IIBDGC/"

dat<-fread(paste(path_gwas,"resources/gwas_summary_statistics/gwas-catalog-v1.0.3.1-studies-r2025-05-13.tsv.gz",sep=""),quote="")
# dat[which(dat$'DISEASE/TRAIT'=="Psoriasis"  & dat$'FULL SUMMARY STATISTICS'=="yes"),]

# european ancestry studies, larger number of hits and sample size, harmonised data in b38 available

# metabolic biomarkers from adam's paper - downloaded in /path/to/project
cyto<-dat[which(dat$'PUBMED ID'==38448586),]
cyto$map_trait<-gsub("http://www.ebi.ac.uk/efo/","",cyto$MAPPED_TRAIT_URI)


# blood cell counts from Chen MH - only EUR for ldsc
cc<-dat[which(dat$'PUBMED ID'==32888493),]
cc<-cc[grep("European ancestry",cc$'INITIAL SAMPLE SIZE'),]
cc<-cc[!grepl("Afro-Caribbean",cc$'INITIAL SAMPLE SIZE'),]

# blood and urine biomarkers 33462484
# bub<-dat[which(dat$'PUBMED ID'==33462484),]

# cancer studies - breast prostate endometrial
cancer<-dat[which(dat$'STUDY ACCESSION' %in% c("GCST004988","GCST90274714","GCST006464",
"GCST004748","GCST90018921")),]

# immune mediated studies
imm<-dat[which(dat$'STUDY ACCESSION' %in% c("GCST90027161","GCST005523",
"GCST90044763","GCST90018847","GCST005531","GCST90061440","GCST004030","GCST002318","GCST003156",
"GCST90014023","GCST90018926","GCST005529","GCST005527",
"GCST90132223","GCST011096","GCST90044158","GCST90472771","GCST90480502","GCST003566","GCST90013445","GCST90029017")),]


# Other - Alzheimer, height, bmi, clonal hematopoyesis
cm<-dat[which(dat$'STUDY ACCESSION' %in% c("GCST90027158","GCST90018959","GCST90018947")),]

# cardiovascular traits:
cv<-dat[which(dat$'STUDY ACCESSION'  %in% c("GCST90476007","GCST90475936","GCST90475967","GCST90475929","GCST90435254","GCST90475214")),]

all<-rbind(cc,cancer,imm,cm,cv)
all$map_trait<-gsub("http://www.ebi.ac.uk/efo/","",all$MAPPED_TRAIT_URI)
all$map_trait<-gsub("http://purl.obolibrary.org/obo/","",all$map_trait)

dim(all)
# [1] 50 26

# GCST003566,GCST90476007 processed separately - issue with harmonization, see separated script
all<-all[which(!all$'STUDY ACCESSION' %in% c("GCST003566","GCST90476007")),]

dim(all)
# [1] 48 26

names(table(all$'DISEASE/TRAIT'))

# two formats of IDs
new<-c("GCST90132223","GCST011096","GCST90472771","GCST90480502","GCST90435254","GCST90476007","GCST90475929","GCST90475936","GCST90475967","GCST90475214","GCST90013445","GCST90274714")
rm(list_tmp)
for (i in 1:nrow(all)) {

    if (!all$'STUDY ACCESSION'[i] %in% new) {
        tmp<-paste("wget ",all$'SUMMARY STATS LOCATION'[i],"/harmonised/",all$'PUBMED ID'[i],"-",all$'STUDY ACCESSION'[i],"-",all$map_trait[i],".h.tsv.gz -O ",path_gwas,"/resources/gwas_summary_statistics/",all$'STUDY ACCESSION'[i],".h.tsv.gz",sep="")
    } else {
        tmp<-paste("wget ",all$'SUMMARY STATS LOCATION'[i],"/harmonised/",all$'STUDY ACCESSION'[i],".h.tsv.gz -O ",path_gwas,"/resources/gwas_summary_statistics/",all$'STUDY ACCESSION'[i],".h.tsv.gz",sep="")
    }
        
    if(!exists("list_tmp")) {
        list_tmp<-tmp
    } else {
        list_tmp<-c(list_tmp,tmp)
    }
}


write.table(list_tmp,paste(path_gwas,"resources/gwas_summary_statistics/list_files_to_download.sh",sep=""),col.names=F,row.names=F,quote=F,sep="")


# create a file with the files downloaded by hgi and final location of the files downloaded locally:
dim(cyto)
# [1] 233  24

cyto$summary_stats<-gsub("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/","",cyto$'SUMMARY STATS LOCATION')
cyto$summary_stats<-paste("/path/to/project",cyto$summary_stats,"/harmonised/",cyto$'STUDY ACCESSION',".h.tsv.gz",sep="")
cyto$'STUDY ACCESSION'


all$summary_stats<-paste("/path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/",all$'STUDY ACCESSION',".h.tsv.gz",sep="")

all<-rbind(cyto,all)

all_list<-c(all$summary_stats)
write.table(all_list,paste(path_gwas,"resources/gwas_summary_statistics/list_path_summary_stats_gwas_catalog",sep=""),col.names=F,row.names=F,quote=F,sep="")

all_list<-c(all$'STUDY ACCESSION')
write.table(all_list,paste(path_gwas,"resources/gwas_summary_statistics/list_study_ids_summary_stats_gwas_catalog",sep=""),col.names=F,row.names=F,quote=F,sep="")

# save sample size:
all$N_tmp<-gsub("[A-z]*","",all$'INITIAL SAMPLE SIZE')
all$N_tmp<-gsub(",","",all$N_tmp)

tmp<-str_split(all$N_tmp," ")

all$N<-NA
for(i in c(1:nrow(all))) {
    all$N[i]<-sum(as.numeric(tmp[[i]]),na.rm=T)
}
write.table(all$N,paste(path_gwas,"resources/gwas_summary_statistics/list_study_sample_size_gwas_catalog",sep=""),col.names=F,row.names=F,quote=F,sep="")


# # subset of cyto - colesterol:
# all1<-all[grepl("cholesterol|Cholesterol",all$'DISEASE/TRAIT'),]

# all_list<-c(all1$summary_stats)
# write.table(all_list,paste(path_gwas,"resources/gwas_summary_statistics/list_path_summary_stats_cholesterol",sep=""),col.names=F,row.names=F,quote=F,sep="")


# for(i in c(1:nrow(imm))) {
#     # print(imm[,c('STUDY ACCESSION',"MAPPED_TRAIT"),])
#     system(paste0("cp ",path_gwas,"resources/gwas_summary_statistics/",imm$'STUDY ACCESSION'[i],"_edited_2.tsv.gz /path/to/project"))
# }

##########################################

q("no")

###############################################################################################
# DOWNLOAD THE FILES

MEM=250
path_gwas="/path/to/ibdgwas/IIBDGC/"

bsub -J"downloads" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G humgen-priority \
-e ${path_gwas}post_imputation/log/download_gwas_catalog_stderr \
-o ${path_gwas}post_imputation/log/download_gwas_catalog_stdout \
"sh ${path_gwas}resources/gwas_summary_statistics/list_files_to_download.sh"

# bsub -J"downloads" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
# -e ${path_gwas}post_imputation/log/download_gwas_catalog_may2025_stderr \
# -o ${path_gwas}post_imputation/log/download_gwas_catalog_may2025_stdout \
# "sh ${path_gwas}resources/gwas_summary_statistics/list_files_to_download_may2025.sh"

cat ${path_gwas}resources/gwas_summary_statistics/list_files_to_download.sh | grep GCST90018959



# cd /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/
# # wget https://portals.broadinstitute.org/collaboration/giant/images/f/f7/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz
# # b37, only 1K variants

# # Sakaue S	GCST90018959	2021-09-30	Nat Genet	A cross-population atlas of genetic associations…	Height	body height	-	• 360388 European
# # • 165056 East Asian		723	FTP Download
# wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90018001-GCST90019000/GCST90018959/harmonised/34594039-GCST90018959-EFO_0004339.h.tsv.gz

# zcat 34594039-GCST90018959-EFO_0004339.h.tsv.gz | wc -l
# # 19172469

# zcat ${path_gwas}resources/gwas_summary_statistics/34594039-GCST90018959-EFO_0004339.h.tsv.gz | cut -f1 | sort | uniq -d 

# # subset to keep variants represented in at least 50% of the samples, and maf>1%
# # reformat variant name - to match chr15:19787290:C:A
# # exclude duplicated IDs

# cat <(zcat 34594039-GCST90018959-EFO_0004339.h.tsv.gz | head -1 | cut -f1,3-7,17,19,20 | sed 's/hm_chrom/CHR/g' | sed 's/hm_pos/BP/g' | sed 's/EFFECT_ALLELE/A1/g' | sed 's/hm_variant_id/SNP/g' | sed 's/hm_other_allele/A0/g' | sed 's/hm_effect_allele/A1/g' | sed 's/hm_beta/BETA/g' | sed 's/effect_allele_frequency/A1FREQ/g' | sed 's/standard_error/SE/g' | sed 's/p_value/PVALUE/g' | awk -F'\t' '{$10="N"; print $0}') \
# <(zcat 34594039-GCST90018959-EFO_0004339.h.tsv.gz | cut -f1,3-7,17,19,20 | awk '{ if ($7 >= 0.001 && $7 <= 0.999) print $0 }' \
# | awk -v OFS=$'\t' '{ $1="chr" $1; print}' | sed 's/_/:/g' | sed 's/X/23/g' | grep -Ev 'chr10:17924539:C:G|chr10:50214925:C:T|chr10:50215018:C:T|chr10:50240776:T:C|chr15:82372034:G:A|chr19:20396180:C:T|chr19:327683:G:A|chr19:39890893:C:T|chr4:157285557:T:TTG|chr4:74446969:A:G|chr4:74450339:C:T|chr4:74456775:A:T|chr5:161312018:T:TTG|chr7:495710:A:G|chr8:139977662:CTCAT:C|chr8:141738107:A:C|chr8:144235357:C:T|chr8:2623362:GGA:G|chr8:530450:C:G' | awk -F'\t' '{$10=525444; print $0}') | \
# gzip > ${path_gwas}resources/gwas_summary_statistics/sakaue_2021_height_list_variants_maf_0.001.tsv.gz

# zcat ${path_gwas}resources/gwas_summary_statistics/sakaue_2021_height_list_variants_maf_0.001.tsv.gz | cut -f7 | sort -n | sed -n '1s/^/min=/p; $s/^/max=/p' 
# # min=EFFECT_ALLELE_FREQ
# # max=0.999

# zcat ${path_gwas}resources/gwas_summary_statistics/sakaue_2021_height_list_variants_maf_0.001.tsv.gz | wc -l
# # 17257606

# zcat ${path_gwas}resources/gwas_summary_statistics/sakaue_2021_height_list_variants_maf_0.001.tsv.gz | head

# # reformat so that it preserves a common name id:
# MEM=9000
# bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
# -e ${path_gwas}post_imputation/log/heigh_sumstats_nohla_common_0.01_mungen_format_stderr \
# -o ${path_gwas}post_imputation/log/heigh_sumstats_nohla_common_0.01_mungen_format_stdout \
# "munge_polyfun_sumstats.py \
# --sumstats ${path_gwas}resources/gwas_summary_statistics/sakaue_2021_height_list_variants_maf_0.001.tsv.gz \
# --min-info 0 \
# --min-maf 0 \
# --out ${path_gwas}resources/gwas_summary_statistics/sakaue_2021_height_list_variants_maf_0.001_munged.parquet"
# # bjobs -l 101972

#########################################################################################################################################################
### reformat all summary statistics as prioritised in IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_45_compute_genetic_correlation.R

cd /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/


files=($(cat /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/list_path_summary_stats_gwas_catalog )) 
echo ${#files[@]}
# 281

# # new files May 2025
# files=($(cat /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/list_path_summary_stats_gwas_catalog_may2025))
# echo ${#files[@]}
# # 15

gwas_id=($(cat /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/list_path_summary_stats_gwas_catalog | \
sed 's/\/lustre\/scratch124\/humgen\/projects_v2\/ibdgwas\/IIBDGC\/resources\/gwas_summary_statistics\///g' | \
sed 's/\/lustre\/scratch125\/humgen\/resources_v2\/GWAScatalog\/summary_statistics\/.*\///g' | sed  's/.h.tsv.gz//g')) 
echo ${#gwas_id[@]}
# 281


gwas_n=($(cat /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/list_study_sample_size_gwas_catalog))
echo ${#gwas_n[@]}
# 281

for i in {0..280}
do
echo ${i} && echo ${files[i]} && echo ${gwas_id[i]} && echo ${gwas_n[i]} 
done


for i in {0..280}
do
ls -la ${files[i]}
done


## find format of files:

# beta (no hm_beta)
for i in {0..280}
do
echo ${i} && echo ${gwas_id[i]} && zcat ${files[i]} | head -1 | grep -w beta | grep -v hm_beta
done 
# for i in {0..232} 249 258 267 270 276

# odds_ratio (no hm_odds_ratio)
for i in {0..280}
do
echo ${i} && echo ${gwas_id[i]} && zcat ${files[i]} | head -1 | grep -w odds_ratio | grep -v hm_odds_ratio
done 
# for i in 254 255 272 {277..280}

# hm_beta
for i in {0..280}
do
echo ${i} && echo ${gwas_id[i]} && zcat ${files[i]} | head -1 | grep -w hm_beta
done 
# for i in {233..248} {250..254} 256 257 {259..266} 268 271 {273..275}


#################################################################################
### FILE FORMATS


############################################
# 1.- EFFECT SIZE IN COLUMN BETA OR ODDS_RATIO  
# no hm_beta flags in the harmonised files, files partially harmonised - manually harmonize beta and frequency

MEM=15000  
for i in {0..232} 249 254 255 258 267 269 270 272 {277..280}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${gwas_id[i]}_may2025_reformat_stderr \
-o ${path_gwas}post_imputation/log/${gwas_id[i]}_may2025_reformat_stdout \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/reformat_partially_harmonized_gwas_catalog_beta_files.R ${gwas_id[i]} ${gwas_n[i]} ${files[i]} > \
${path_gwas}post_imputation/2022/log/eformat_partially_harmonized_gwas_catalog_beta_files_${gwas_id[i]}.Rout"
done

for i in {0..232} 249 254 255 258 267 270 272 {277..280}
do
echo ${i} && tail -30 ${path_gwas}post_imputation/log/${gwas_id[i]}_may2025_reformat_stdout | grep -E "Successfully|Exited"
done


for i in {0..232} 249 255 254 258 267 270 272 {277..280}
do
echo ${i} && cat ${path_gwas}post_imputation/2022/log/eformat_partially_harmonized_gwas_catalog_beta_files_${gwas_id[i]}.Rout && echo "############" && echo "############"
done

############################################
# 2.- EFFECT SIZE IN COLUMN HM_BETA 
# # default harmonised format from GWAS catalog - with beta as effect estimator, and hm flags

MEM=500
for i in {233..248} {250..253} 256 257 {259..266} 268 269 271 {273..275}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${gwas_id[i]}_may2025_reformat_stderr \
-o ${path_gwas}post_imputation/log/${gwas_id[i]}_may2025_reformat_stdout \
"cat <(zcat ${files[i]} | head -1 | csvcut -t -d $'\t' -c hm_variant_id,hm_chrom,hm_pos,hm_effect_allele,hm_other_allele,effect_allele_frequency,hm_beta,standard_error,p_value | \
sed 's/,/\\t/g' | \
sed 's/chromosome/CHR/g;s/hm_chrom/CHR/g' | sed 's/hm_pos/BP/g;s/base_pair_location/BP/g' | \
sed 's/EFFECT_ALLELE/A1/g;s/hm_effect_allele/A1/g;s/effect_allele/A1/g;s/effect_allele/A1/g' | \
sed 's/hm_variant_id/SNP/g;s/variant_id/SNP/g' | \
sed 's/hm_other_allele/A0/g;s/other_allele/A0/g' | sed 's/hm_beta/BETA/g;s/beta/BETA/g' | \
sed 's/effect_allele_frequency/A1FREQ/g;s/A1_frequency/A1FREQ/g' | sed 's/standard_error/SE/g' | sed 's/p_value/PVALUE/g' | \
awk -F'\t' '{\$10=\"N\"; print \$0}') \
<(zcat ${files[i]} | csvcut -t -d $'\t' -c hm_variant_id,hm_chrom,hm_pos,hm_effect_allele,hm_other_allele,hm_effect_allele_frequency,hm_beta,standard_error,p_value | sed 's/,/\t/g' | \
awk -v OFS=$'\t' '{\$1=\"chr\"  \$1; print}' | sed 's/_/:/g' | sed 's/X/23/g' | awk -F'\t' '{\$10=${gwas_n[i]}; print \$0}' | tail -n +2) | \
gzip > ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_edited.tsv.gz"
done

for i in {233..248} {250..253} 256 257 {259..266} 268 269 271 {273..275}
do
echo ${i} && tail -30 ${path_gwas}post_imputation/log/${gwas_id[i]}_may2025_reformat_stdout | grep -E "Successfully|Exited"
done

rm ~/tmp.txt
for i in  {0..280}
do
echo ${i} && ls -la ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_edited.tsv.gz
done


# TO DO - DOUBLE CHECK FREQUENCY BOUNDARY OF VARIANTS
for i in  {0..280}
do
echo ${gwas_id[i]} && zcat ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_edited.tsv.gz | head -2
done

########################################

# manually add additional files - add frequency later

# files that fail:
# i= 257 259 260 261 262 255 254 256 258 253

# GCST005531_edited.tsv.gz - 253 Multiple sclerosis no effect allele frequency
# GCST005523_edited.tsv.gz - 254 Celiac disease - no effect allele frequency
# GCST005527_edited.tsv.gz - 256 Psoriais no effect allele frequency
# GCST003156_edited.tsv.gz - 257 lupus, no effect allele frequency
# GCST005529_edited.tsv.gz - 258 Ankylosing spondylitis no effect allele frequency


##### ONLY IIBDGC SNPS

# list variants present in IBD CD or UC summary statistics with similar frequency, exclude duplicates, and keep only variants selected in 
# IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_45_compute_genetic_correlation.R
# path_gwas,"post_imputation/2022/log/subset_gwas_variants_for_genetic_correlation.snplist


# DEFINE THE LIST OF GOOD QUALITY SNPS TO USE ARE REFERNCE:
# ~/git/IIBDGC_GWAS/scripts/other/create_reference_file_allchr_ibd_cd_uc_eur_tier2_list_variants_rate_0.5_info.R
# ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_ibd_cd_uc_uc_eur_tier2_list_variants_rate_05_info_0.9.tsvgz

MEM=30000

for i in {0..280}
do
bsub -J"subset_gwas" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-o ${path_gwas}post_imputation/2022/log/subset_gwas_variants_may2025_${gwas_id[i]}_stdout \
-e ${path_gwas}post_imputation/2022/log/subset_gwas_variants_may2025_${gwas_id[i]}_stderr \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/subset_gwas_variants_for_genetic_correlation.R ${gwas_id[i]} > \
${path_gwas}post_imputation/2022/log/subset_gwas_variants_for_genetic_correlation_${gwas_id[i]}.Rout"
done

# CONTINUE HERE


for i in {0..280}
do
echo ${i} && tail -50 ${path_gwas}post_imputation/2022/log/subset_gwas_variants_may2025_${gwas_id[i]}_stdout  | grep -E "Successfully|Exit"
# echo ${i} && tail -30 ${path_gwas}post_imputation/2022/log/subset_gwas_variants_may2025_${gwas_id[i]}_stdout | grep -E "Max Memory"
done

for i in {0..280}
do
echo ${i} && ls -la ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_edited_2.tsv.gz
done

rm ~/tmp.txt
for i in {0..280}
do
echo ${i} >> ~/tmp.txt && cat ${path_gwas}post_imputation/2022/log/subset_gwas_variants_for_genetic_correlation_${gwas_id[i]}.Rout >> ~/tmp.txt && \
echo "###########" >> ~/tmp.txt && echo "###########" >> ~/tmp.txt
done

#######################################################################################################
# reformat to munged file

path_gwas="/path/to/ibdgwas/IIBDGC/"

MEM=7000

for i in {0..280}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${gwas_id[i]}_ldsc_sumstats_only_SNPs_polyfun_mungen_format_stderr \
-o ${path_gwas}post_imputation/log/${gwas_id[i]}_ldsc_sumstats_only_SNPs_polyfun_mungen_format_stdout \
"munge_polyfun_sumstats.py \
--sumstats ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_edited_2.tsv.gz \
--min-info 0 \
--min-maf 0 \
--remove-strand-ambig \
--out ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_only_SNPs_sumstats_munged.parquet"
done

for i in {0..280}
do
echo ${i} && tail -50 ${path_gwas}post_imputation/log/${gwas_id[i]}_ldsc_sumstats_only_SNPs_polyfun_mungen_format_stdout  | grep -E "Successfully|Exit"
done

for i in {0..280}
do
echo ${i} && ls -la ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_only_SNPs_sumstats_munged.parquet
done


# SAVE PARQUET FILE AS FLAT TEXT FILE AS INPUT FOR LDSC/LDSC.PY

module unload HGI/softpack/groups/team152/iibdgc_postprocess/10
module unload HGI/softpack/groups/team152/polyfun-2/1

MEM=3000
for i in {0..280}
do
bsub -J"subset_gwas" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-o ${path_gwas}post_imputation/2022/log/munge_parquet_to_sumstats_${gwas_id[i]}_stdout \
-e ${path_gwas}post_imputation/2022/log/munge_parquet_to_sumstats_${gwas_id[i]}_stderr \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/reformat_munge_parquet_to_sumtats_for_genetic_correlation.R ${gwas_id[i]} > \
${path_gwas}post_imputation/2022/log/reformat_munge_parquet_to_sumtats_for_genetic_correlation_${gwas_id[i]}.Rout"
done

for i in {0..280}
do
echo ${i} && tail -50 ${path_gwas}post_imputation/2022/log/munge_parquet_to_sumstats_${gwas_id[i]}_stdout | grep -E "Successfully|Exit"
done

for i in {0..280}
do
echo ${i} && ls -la ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_only_SNPs_sumstats_munged.sumstats
done

# COMPLETED JUNE 25


# interactive submission:
# cat <(zcat ${files[i]} | head -1 | cut -f1,3-7,17,19,20 | \
# sed 's/chromosome/CHR/g;s/hm_chrom/CHR/g' | sed 's/hm_pos/BP/g;s/base_pair_location/BP/g' | \
# sed 's/EFFECT_ALLELE/A1/g;s/hm_effect_allele/A1/g;s/effect_allele/A1/g;s/effect_allele/A1/g' | \
# sed 's/hm_variant_id/SNP/g;s/variant_id/SNP/g' | \
# sed 's/hm_other_allele/A0/g;s/other_allele/A0/g' | sed 's/hm_beta/BETA/g;s/beta/BETA/g' | sed 's/effect_allele_frequency/A1FREQ/g;s/A1_frequency/A1FREQ/g' | sed 's/standard_error/SE/g' | sed 's/p_value/PVALUE/g' \
# | awk -F'\t' '{$10="N"; print $0}') \
# <(zcat ${files[i]} | cut -f1,3-7,17,19,20 | awk '{ if ($7 >= 0.001 && $7 <= 0.999) print $0 }' \
# | awk -v OFS=$'\t' '{ $1="chr" $1; print}' | sed 's/_/:/g' | sed 's/X/23/g' | awk -F'\t' '{$10="'"${gwas_n[i]}"'"; print $0}') | \
# gzip > ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_edited.tsv.gz

# zcat ${files[i]} | head -10 | csvcut -t -d $'\t' -c variant_id,chromosome,base_pair_location,effect_allele,other_allele,effect_allele_frequency,beta,standard_error,p_value | sed 's/,/\t/g'
# zcat ${files[i]} | head -10 | csvcut -t -d $'\t' -c hm_variant_id,hm_chrom,hm_pos,hm_effect_allele,hm_other_allele,effect_allele_frequency,hm_beta,standard_error,p_value | sed 's/,/\t/g'

# find and replace duplicates
# zcat ${files[i]} | cut -f1 | sort | uniq -d | sed -e 'H;${x;s/\n/|/g;s/^|//;p;};d' | sed '1s/^/"/;1s/$/"/'

# zcat ${files[i]} | grep -Ev -f <(zcat ${files[i]} | csvcut -t -d $'\t' -c variant_id | sort | uniq -d | sed -e 'H;${x;s/\n/|/g;s/^|//;p;};d' | sed '1s/^/"/;1s/$/"/') | wc -l
# zcat ${files[i]} | grep -Ev -f <(zcat ${files[i]} | csvcut -t -d $'\t' -c hm_variant_id | sort | uniq -d | sed -e 'H;${x;s/\n/|/g;s/^|//;p;};d' | sed '1s/^/"/;1s/$/"/') | wc -l


# slightly different format in variant ID, so that equals chr_posb38_alt_ref from new files from Adam's study, 
# # after being harmonised by GWAS catalog, modify to be able to compare later

# # to filter by AF:
# # awk '{ if (\$6 >= 0.001 && \$6 <= 0.999) print \$0 }' | \



# ### OLD SUBMISSIONS:
# ############################################
# # 1.- EFFECT SIZE IN COLUMN BETA 
# # no hm_beta flags in the harmonised files, files partially harmonised

# MEM=250 
# for i in {0..232}
# do
# bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
# -e ${path_gwas}post_imputation/log/${gwas_id[i]}_reformat_stderr \
# -o ${path_gwas}post_imputation/log/${gwas_id[i]}_reformat_stdout \
# "cat <(zcat ${files[i]} | head -1 | csvcut -t -d $'\t' -c variant_id,chromosome,base_pair_location,other_allele,effect_allele,effect_allele_frequency,beta,standard_error,p_value | \
# sed 's/,/\\t/g' | \
# sed 's/chromosome/CHR/g;s/hm_chrom/CHR/g' | sed 's/hm_pos/BP/g;s/base_pair_location/BP/g' | \
# sed 's/EFFECT_ALLELE/A1/g;s/hm_effect_allele/A1/g;s/effect_allele/A1/g;s/effect_allele/A1/g' | \
# sed 's/hm_variant_id/SNP/g;s/variant_id/SNP/g' | \
# sed 's/hm_other_allele/A0/g;s/other_allele/A0/g' | sed 's/hm_beta/BETA/g;s/beta/BETA/g' | \
# sed 's/effect_allele_frequency/A1FREQ/g;s/A1_frequency/A1FREQ/g' | sed 's/standard_error/SE/g' | sed 's/p_value/PVALUE/g' | \
# awk -F'\t' '{\$10=\"N\"; print \$0}') \
# <(zcat ${files[i]} | csvcut -t -d $'\t' -c variant_id,chromosome,base_pair_location,other_allele,effect_allele,effect_allele_frequency,beta,standard_error,p_value | sed 's/,/\t/g' | \
# awk '{print \$2 \":\" \$3 \":\" \$4 \":\" \$5,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' | \
# awk -v OFS=$'\t' '{\$1=\"chr\" \$1; print}' | sed 's/_/:/g' | sed 's/X/23/g' | awk -F'\t' '{\$10=${gwas_n[i]}; print \$0}' | tail -n +2) | \
#  gzip > ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_edited.tsv.gz"
# done

# ############################################
# # 2.- EFFECT SIZE IN COLUMN HM_BETA 
# # # default harmonised format from GWAS catalog - with beta as effect estimator, and hm flags

# # for i in {233..251} {253..258} {260..268}
# for i in 3 6
# do
# bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
# -e ${path_gwas}post_imputation/log/${gwas_id[i]}_reformat_stderr \
# -o ${path_gwas}post_imputation/log/${gwas_id[i]}_reformat_stdout \
# "cat <(zcat ${files[i]} | head -1 | csvcut -t -d $'\t' -c hm_variant_id,hm_chrom,hm_pos,hm_effect_allele,hm_other_allele,effect_allele_frequency,hm_beta,standard_error,p_value | \
# sed 's/,/\\t/g' | \
# sed 's/chromosome/CHR/g;s/hm_chrom/CHR/g' | sed 's/hm_pos/BP/g;s/base_pair_location/BP/g' | \
# sed 's/EFFECT_ALLELE/A1/g;s/hm_effect_allele/A1/g;s/effect_allele/A1/g;s/effect_allele/A1/g' | \
# sed 's/hm_variant_id/SNP/g;s/variant_id/SNP/g' | \
# sed 's/hm_other_allele/A0/g;s/other_allele/A0/g' | sed 's/hm_beta/BETA/g;s/beta/BETA/g' | \
# sed 's/effect_allele_frequency/A1FREQ/g;s/A1_frequency/A1FREQ/g' | sed 's/standard_error/SE/g' | sed 's/p_value/PVALUE/g' | \
# awk -F'\t' '{\$10=\"N\"; print \$0}') \
# <(zcat ${files[i]} | csvcut -t -d $'\t' -c hm_variant_id,hm_chrom,hm_pos,hm_effect_allele,hm_other_allele,hm_effect_allele_frequency,hm_beta,standard_error,p_value | sed 's/,/\t/g' | \
# awk -v OFS=$'\t' '{\$1=\"chr\"  \$1; print}' | sed 's/_/:/g' | sed 's/X/23/g' | awk -F'\t' '{\$10=${gwas_n[i]}; print \$0}' | tail -n +2) | \
# gzip > ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_edited.tsv.gz"
# done


# ############################################
# # 3.- EFFECT SIZE IN COLUMN HM_ODDS_RATIO 
# # # # GCST004030_edited.tsv.gz - 259 psc - no beta, only OR, RESUBMIT - estimate beta from OR
# # # default harmonised format from GWAS catalog - with only OR as effect estimator

# MEM=500
# # i=259
# for i in 0
# do
# bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
# -e ${path_gwas}post_imputation/log/${gwas_id[i]}_reformat_stderr \
# -o ${path_gwas}post_imputation/log/${gwas_id[i]}_reformat_stdout \
# "cat <(zcat ${files[i]} | head -1 | csvcut -t -d $'\t' -c hm_variant_id,hm_chrom,hm_pos,hm_effect_allele,hm_other_allele,effect_allele_frequency,hm_beta,standard_error,p_value | \
# sed 's/,/\\t/g' | \
# sed 's/chromosome/CHR/g;s/hm_chrom/CHR/g' | sed 's/hm_pos/BP/g;s/base_pair_location/BP/g' | \
# sed 's/EFFECT_ALLELE/A1/g;s/hm_effect_allele/A1/g;s/effect_allele/A1/g;s/effect_allele/A1/g' | \
# sed 's/hm_variant_id/SNP/g;s/variant_id/SNP/g' | \
# sed 's/hm_other_allele/A0/g;s/other_allele/A0/g' | sed 's/hm_beta/BETA/g;s/beta/BETA/g;s/hm_odds_ratio/BETA/g' | \
# sed 's/effect_allele_frequency/A1FREQ/g;s/A1_frequency/A1FREQ/g' | sed 's/standard_error/SE/g' | sed 's/p_value/PVALUE/g' | \
# awk -F'\t' '{\$10=\"N\"; print \$0}') \
# <(zcat ${files[i]} | csvcut -t -d $'\t' -c hm_variant_id,hm_chrom,hm_pos,hm_effect_allele,hm_other_allele,hm_effect_allele_frequency,hm_odds_ratio,standard_error,p_value | sed 's/,/\t/g' | \
# awk '{\$10=log(\$7); print \$1,\$2,\$3,\$4,\$5,\$6,\$10,\$8,\$9}' | \
# awk -v OFS=$'\t' '{\$1=\"chr\"  \$1; print}' | sed 's/_/:/g' | sed 's/X/23/g' | awk -F'\t' '{\$10=${gwas_n[i]}; print \$0}' | tail -n +2) | \
# gzip > ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_edited.tsv.gz"
# done

# zcat /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/GCST002318.h.tsv.gz | head -10 | csvcut -t -d $'\t' -c hm_variant_id,hm_chrom,hm_pos,hm_effect_allele,hm_other_allele,hm_effect_allele_frequency,hm_odds_ratio,standard_error,p_value | sed 's/,/\t/g' | awk '{$10=log($7); print $1,$2,$3,$4,$5,$6,$10,$8,$9}' | awk -v OFS=$'\t' '{$1="chr"  $1; print}' | sed 's/_/:/g' | sed 's/X/23/g' | awk -F'\t' '{$10=79799; print $0}' | tail -n +2

# ############################################
# # 4.- EFFECT SIZE IN COLUMN ODDS_RATIO  
# #  no hm_ flags and only OR:

# MEM=500
# # i=259
# for i in 8 {10..14}
# do
# bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
# -e ${path_gwas}post_imputation/log/${gwas_id[i]}_reformat_stderr \
# -o ${path_gwas}post_imputation/log/${gwas_id[i]}_reformat_stdout \
# "cat <(zcat ${files[i]} | head -1 | csvcut -t -d $'\t' -c variant_id,chromosome,base_pair_location,effect_allele,other_allele,effect_allele_frequency,odds_ratio,standard_error,p_value | \
# sed 's/,/\\t/g' | \
# sed 's/chromosome/CHR/g;s/hm_chrom/CHR/g' | sed 's/hm_pos/BP/g;s/base_pair_location/BP/g' | \
# sed 's/EFFECT_ALLELE/A1/g;s/hm_effect_allele/A1/g;s/effect_allele/A1/g;s/effect_allele/A1/g' | \
# sed 's/hm_variant_id/SNP/g;s/variant_id/SNP/g' | \
# sed 's/hm_other_allele/A0/g;s/other_allele/A0/g' | sed 's/hm_beta/BETA/g;s/beta/BETA/g;s/odds_ratio/BETA/g' | \
# sed 's/effect_allele_frequency/A1FREQ/g;s/A1_frequency/A1FREQ/g' | sed 's/standard_error/SE/g' | sed 's/p_value/PVALUE/g' | \
# awk -F'\t' '{\$10=\"N\"; print \$0}') \
# <(zcat ${files[i]} | csvcut -t -d $'\t' -c variant_id,chromosome,base_pair_location,effect_allele,other_allele,effect_allele_frequency,odds_ratio,standard_error,p_value | sed 's/,/\t/g' | \
# awk '{\$10=log(\$7); print \$2 \":\" \$3 \":\" \$4 \":\" \$5,\$2,\$3,\$4,\$5,\$6,\$10,\$8,\$9}' | \
# awk -v OFS=$'\t' '{\$1=\"chr\"  \$1; print}' | sed 's/_/:/g' | sed 's/X/23/g' | awk -F'\t' '{\$10=${gwas_n[i]}; print \$0}' | tail -n +2) | \
# gzip > ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_edited.tsv.gz"
# done
