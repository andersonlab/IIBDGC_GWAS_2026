# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# double check we are not leaving out any region by using COJO to determine GWAS significant hits:

# give bfile as input

# --bfile test 
# Note: --bgen, --pfile, --bpfile, --mbfile, --mbgen, --mpfile and --mbpfile are currently only supported in the GRM and fastGWA analyses. 
# More functions will be available after rewriting some of the legacy codes. All the QC flags (e.g. --keep, --extract, --maf) in GCTA 1.92.4 are 
# currently applicable to these two formats.


# need summary statistics in the format they request

# The choice of reference sample for GCTA-COJO analysis
# 1) If the summary data are from a single cohort based GWAS, the best reference sample is the GWAS sample itself.
# 2) For a meta-analysis where individual-level genotype data are not available, you could use one of the large participating cohorts. For example, 
# Wood et al. 2014 Nat Genet used the ARIC cohort (data available from dbGaP).
# 3) We suggest you use a reference sample with a sample size > 4000 (see Supplementary Figure 4 of Yang et al. 2012 Nat Genet).
# 4) We do NOT suggest you use HapMap or 1000G panels as the reference sample. The sample sizes of HapMap and 1000G are not large enough.


# --sample ${path}imputed_data/release_interval_ibdbioresource_dec22/combined_hrc_uk10k/${chr}_combined_hrc_uk10k_${i}_with_consensus_ids.sample \
# --phenoFile ${path}/username/phenotype_data/freeze_20230428/pheno_${ph}_${j}_interval_ibdbioresource \
# --covarFile ${path}/username/phenotype_data/freeze_20230428/covariates_${j}_interval_ibdbioresource \



##########################################
# 1.- Define 

# 1.1- subset by list samples included in each analysis (CD|UC|IBD) - largest dataset is ibdbioresource:

path=/path/to/ibdgwas_bioresource/
path_gwas=/path/to/ibdgwas/IIBDGC/
  
pheno=(cd uc ibd)
  
for ph in ${pheno[@]}
do
awk -F '\t' '$3 != "NA" { print }' ${path}/username/phenotype_data/freeze_20230428/pheno_${ph}_sas_interval_ibdbioresource \
| cut -f1-2 > ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_list_samples_included_in_${ph}_analysis
done

for ph in ${pheno[@]}
do 
echo ${ph} && wc -l ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_list_samples_included_in_${ph}_analysis
done

# cd
# 1129 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_list_samples_included_in_cd_analysis
# uc
# 1332 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_list_samples_included_in_uc_analysis
# ibd
# 1882 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_list_samples_included_in_ibd_analysis


# 1.2- subset by variants with MAF 0.001 and info 0.4 as in regenie output for gsa

for chr in {1..22} X
do
for ph in ${pheno[@]}
do
zcat ${path}username/analysis/regenie/${ph}/interval_ibdbioresource_step2_${ph}_noneur_chip_sex_PCs_chr${chr}_${ph}_sas.regenie.gz | \
awk -F ' ' '$6 >= 0.001 && $6 <= 0.999 && $9 >= 0.4 { print $3 }' \
> ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr${chr}_list_variants_included_in_${ph}_analysis
done
done


for ph in ${pheno[@]}
do 
echo ${ph} && for chr in {1..22} X
do
echo ${chr} && wc -l  ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr${chr}_list_variants_included_in_${ph}_analysis
done
done


##################
# cd
# 1
# 1354390 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr1_list_variants_included_in_cd_analysis
# 2
# 1477156 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr2_list_variants_included_in_cd_analysis
# 3
# 1244207 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr3_list_variants_included_in_cd_analysis
# 4
# 1256013 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr4_list_variants_included_in_cd_analysis
# 5
# 1123436 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr5_list_variants_included_in_cd_analysis
# 6
# 1150526 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr6_list_variants_included_in_cd_analysis
# 7
# 1015015 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr7_list_variants_included_in_cd_analysis
# 8
# 968010 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr8_list_variants_included_in_cd_analysis
# 9
# 747697 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr9_list_variants_included_in_cd_analysis
# 10
# 877444 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr10_list_variants_included_in_cd_analysis
# 11
# 869243 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr11_list_variants_included_in_cd_analysis
# 12
# 842735 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr12_list_variants_included_in_cd_analysis
# 13
# 629108 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr13_list_variants_included_in_cd_analysis
# 14
# 569043 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr14_list_variants_included_in_cd_analysis
# 15
# 507148 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr15_list_variants_included_in_cd_analysis
# 16
# 547123 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr16_list_variants_included_in_cd_analysis
# 17
# 475172 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr17_list_variants_included_in_cd_analysis
# 18
# 483842 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr18_list_variants_included_in_cd_analysis
# 19
# 393822 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr19_list_variants_included_in_cd_analysis
# 20
# 373961 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr20_list_variants_included_in_cd_analysis
# 21
# 226286 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr21_list_variants_included_in_cd_analysis
# 22
# 234250 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr22_list_variants_included_in_cd_analysis
# X
# 641914 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chrX_list_variants_included_in_cd_analysis
# uc
# 1
# 1352790 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr1_list_variants_included_in_uc_analysis
# 2
# 1477175 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr2_list_variants_included_in_uc_analysis
# 3
# 1241277 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr3_list_variants_included_in_uc_analysis
# 4
# 1253789 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr4_list_variants_included_in_uc_analysis
# 5
# 1120218 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr5_list_variants_included_in_uc_analysis
# 6
# 1149593 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr6_list_variants_included_in_uc_analysis
# 7
# 1014707 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr7_list_variants_included_in_uc_analysis
# 8
# 967417 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr8_list_variants_included_in_uc_analysis
# 9
# 745275 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr9_list_variants_included_in_uc_analysis
# 10
# 874929 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr10_list_variants_included_in_uc_analysis
# 11
# 869938 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr11_list_variants_included_in_uc_analysis
# 12
# 842269 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr12_list_variants_included_in_uc_analysis
# 13
# 627531 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr13_list_variants_included_in_uc_analysis
# 14
# 568445 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr14_list_variants_included_in_uc_analysis
# 15
# 505937 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr15_list_variants_included_in_uc_analysis
# 16
# 547773 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr16_list_variants_included_in_uc_analysis
# 17
# 471825 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr17_list_variants_included_in_uc_analysis
# 18
# 482746 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr18_list_variants_included_in_uc_analysis
# 19
# 393423 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr19_list_variants_included_in_uc_analysis
# 20
# 373497 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr20_list_variants_included_in_uc_analysis
# 21
# 226344 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr21_list_variants_included_in_uc_analysis
# 22
# 233389 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr22_list_variants_included_in_uc_analysis
# X
# 640678 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chrX_list_variants_included_in_uc_analysis
# ibd
# 1
# 1350899 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr1_list_variants_included_in_ibd_analysis
# 2
# 1475883 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr2_list_variants_included_in_ibd_analysis
# 3
# 1240732 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr3_list_variants_included_in_ibd_analysis
# 4
# 1253852 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr4_list_variants_included_in_ibd_analysis
# 5
# 1121025 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr5_list_variants_included_in_ibd_analysis
# 6
# 1149047 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr6_list_variants_included_in_ibd_analysis
# 7
# 1013829 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr7_list_variants_included_in_ibd_analysis
# 8
# 966448 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr8_list_variants_included_in_ibd_analysis
# 9
# 746123 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr9_list_variants_included_in_ibd_analysis
# 10
# 874524 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr10_list_variants_included_in_ibd_analysis
# 11
# 869035 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr11_list_variants_included_in_ibd_analysis
# 12
# 841809 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr12_list_variants_included_in_ibd_analysis
# 13
# 626886 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr13_list_variants_included_in_ibd_analysis
# 14
# 567644 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr14_list_variants_included_in_ibd_analysis
# 15
# 505860 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr15_list_variants_included_in_ibd_analysis
# 16
# 546145 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr16_list_variants_included_in_ibd_analysis
# 17
# 472033 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr17_list_variants_included_in_ibd_analysis
# 18
# 482419 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr18_list_variants_included_in_ibd_analysis
# 19
# 392647 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr19_list_variants_included_in_ibd_analysis
# 20
# 372990 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr20_list_variants_included_in_ibd_analysis
# 21
# 226097 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr21_list_variants_included_in_ibd_analysis
# 22
# 233060 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr22_list_variants_included_in_ibd_analysis
# X
# 635338 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chrX_list_variants_included_in_ibd_analysis##################

##########################################
# 2.- Convert bgen to bfile to run cojo

# --hard-call-threshold 0.2 (~80% PP calls)
# --hard-call-threshold 0.4 (~60% PP calls)

path_gwas=/path/to/ibdgwas/IIBDGC/
pheno=(ibd cd uc)

MEM=4000

for ph in ${pheno[@]}
do 
for chr in {1..22} X
do
bsub -J"bgen_bed" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 4 \
-o ${path_gwas}post_imputation/2022/log/interval_ibdbioresource_${ph}_${chr}_bgen_to_bed_stdout \
-e ${path_gwas}post_imputation/2022/log/interval_ibdbioresource_${ph}_${chr}_bgen_to_bed_stderr \
"/path/to/software/./plink2 \
--bgen ${path}imputed_data/release_interval_ibdbioresource_dec22/combined_hrc_uk10k/${chr}_combined_hrc_uk10k_noneur_with_consensus_ids.bgen ref-first \
--sample ${path}imputed_data/release_interval_ibdbioresource_dec22/combined_hrc_uk10k/${chr}_combined_hrc_uk10k_noneur_with_consensus_ids.sample \
--hard-call-threshold 0.4 \
--keep ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_list_samples_included_in_${ph}_analysis \
--threads 4 \
--keep-allele-order --allow-no-sex \
--extract ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr${chr}_list_variants_included_in_${ph}_analysis \
--make-bed --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr${chr}_subset_included_in_${ph}_analysis"
done
done

for ph in ${pheno[@]}
do 
echo ${ph} && for chr in {1..22} X
do
echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/interval_ibdbioresource_${ph}_${chr}_bgen_to_bed_stdout | grep "Successfully"
done
done



##########################################
# 2.1 liftover to b38 

#### 
# MEM=10000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1


path_gwas="/path/to/ibdgwas/IIBDGC/"
path="/path/to/ibdgwas_bioresource/"

library(data.table)


pheno<-c("ibd","cd","uc")

for (ph in pheno) {
  
  print(ph)
  
  for (chr in c(1:22,"X")) {
    
    print(chr)
    
    bim<-fread(paste(path_gwas,"post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr",chr,"_subset_included_in_",ph,"_analysis.bim",sep=""),head=F)
    print(paste("N variants ",nrow(bim)))
    
    info<-fread(paste(path,"imputed_data/release_interval_ibdbioresource_dec22/combined_hrc_uk10k/",chr,"_combined_hrc_uk10k_eur_with_consensus_ids_imputation_info_hwe_ctr_marker_names_b37_b38.tsv.gz",sep=""),head=T)
    info<-as.data.frame(info)
    
    bim<-merge(bim,info[,c("id_b37","id_b38","position_b38","Ref_b38","Alt_b38")],by.x="V2",by.y="id_b37",all.x=T,sort=T)
    
    list_exclude_1<-bim[which( !(bim$V5==bim$Alt_b38 & bim$V6==bim$Ref_b38) & !(bim$V6==bim$Alt_b38 & bim$V5==bim$Ref_b38) ),]
    list_exclude_2<-bim[which(is.na(bim$id_b38)),]
    
    list_exclude<-rbind(list_exclude_1,list_exclude_2)
    
    print(paste("N variants to exclude",nrow(list_exclude)))
    write.table(list_exclude$V2,paste(path_gwas,"post_imputation/2022/analysis/conditional_analysis/sas/list_variants_chr",chr,"_",ph,"_to_exclude",sep=""),
                col.names=F,row.names=F,quote=F,sep="\t")
    
    
    list_to_update<-bim[which(!bim$V2 %in% list_exclude$V2),c("V2","Ref_b38")]
    print(paste("N variants to update",nrow(list_to_update)))
    write.table(list_to_update,paste(path_gwas,"post_imputation/2022/analysis/conditional_analysis/sas/list_alleles_chr",chr,"_",ph,"_to_update",sep=""),
                col.names=F,row.names=F,quote=F,sep="\t")
    
    
    rm(info,bim,list_to_update,list_exclude,list_exclude_1,list_exclude_2)
  }
  
}

q("no")


####################
####################

path_gwas=/path/to/ibdgwas/IIBDGC/
pheno=(ibd cd uc)

MEM=4000

for ph in ${pheno[@]}
do 
for chr in {1..22}
do
bsub -J"bgen_bed" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 4 \
-o ${path_gwas}post_imputation/2022/log/interval_ibdbioresource_${ph}_${chr}_ref_all_stdout \
-e ${path_gwas}post_imputation/2022/log/interval_ibdbioresource_${ph}_${chr}_ref_all_stderr \
"/path/to/software/./plink2 \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr${chr}_subset_included_in_${ph}_analysis \
--exclude ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/list_variants_chr${chr}_${ph}_to_exclude \
--ref-allele 'force' ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/list_alleles_chr${chr}_${ph}_to_update 2 1 \
--make-bed --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr${chr}_subset_included_in_${ph}_analysis_tmp"
done
done



# for ph in ${pheno[@]}
# do 
# for chr in {1..22} X
# do
# bsub -J"bgen_bed" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 4 \
# -o ${path_gwas}post_imputation/2022/log/interval_ibdbioresource_${ph}_${chr}_freq_stdout \
# -e ${path_gwas}post_imputation/2022/log/interval_ibdbioresource_${ph}_${chr}_frea_stderr \
# /path/to/software/./plink2 \
# --bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr${chr}_subset_included_in_${ph}_analysis_tmp \
# --freq \
# --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr${chr}_subset_included_in_${ph}_analysis_tmp
# 

####################

# MEM=10000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1


path_gwas="/path/to/ibdgwas/IIBDGC/"
path="/path/to/ibdgwas_bioresource/"

library(data.table)

pheno<-c("cd","uc")

for (ph in pheno) {
  
  print(ph)
  
  for (chr in c(1:22)) {
    
    print(chr)
    
    bim<-fread(paste(path_gwas,"post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr",chr,"_subset_included_in_",ph,"_analysis_tmp.bim",sep=""),head=F)
    
    info<-fread(paste(path,"imputed_data/release_interval_ibdbioresource_dec22/combined_hrc_uk10k/",chr,"_combined_hrc_uk10k_eur_with_consensus_ids_imputation_info_hwe_ctr_marker_names_b37_b38.tsv.gz",sep=""),head=T)
    info<-as.data.frame(info)
    
    bim<-merge(bim,info[,c("id_b37","id_b38","position_b38","Ref_b38","Alt_b38","Alt_allele_freq_ctr_b37")],by.x="V2",by.y="id_b37",all.x=T,sort=F)

    print(table(bim$V5==bim$Alt_b38))
    print(table(bim$V6==bim$Ref_b38))
    
    bim<-bim[,c("V1","id_b38","V3","position_b38","V5","V6","V2")]
    
    bim_test<-fread(paste(path_gwas,"post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr",chr,"_subset_included_in_",ph,"_analysis_tmp.bim",sep=""),head=F)
    table(bim_test$V2==bim$V2)
    
    bim<-bim[,1:6]
    
  
    write.table(bim,paste(path_gwas,"post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr",chr,"_subset_included_in_",ph,"_analysis_tmp_edited.bim",sep=""),
                col.names=F,row.names=F,quote=F,sep="\t")
    
    rm(bim,bim_test,info)
  }
  
}

q("no")

################################


path_gwas="/path/to/ibdgwas/IIBDGC/"
path="/path/to/ibdgwas_bioresource/"
pheno=(ibd cd uc)

for ph in ${pheno[@]}
do 
for chr in {1..22}
do
bsub -J"bgen_bed" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 4 \
-o ${path_gwas}post_imputation/2022/log/interval_ibdbioresource_${ph}_${chr}_ref_all_1_stdout \
-e ${path_gwas}post_imputation/2022/log/interval_ibdbioresource_${ph}_${chr}_ref_all_1_stderr \
"/path/to/software/./plink2 \
--bed ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr${chr}_subset_included_in_${ph}_analysis_tmp.bed \
--bim ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr${chr}_subset_included_in_${ph}_analysis_tmp_edited.bim \
--fam ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr${chr}_subset_included_in_${ph}_analysis_tmp.fam \
--make-pgen --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr${chr}_subset_included_in_${ph}_analysis_b38"
done
done


for ph in ${pheno[@]}
do 
echo ${ph} && for chr in {1..22}
do
echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/interval_ibdbioresource_${ph}_${chr}_ref_all_1_stdout | grep "Successfully"
done
done



for ph in ${pheno[@]}
do 
for chr in {1..22}
do
bsub -J"bgen_bed" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 4 \
-o ${path_gwas}post_imputation/2022/log/interval_ibdbioresource_${ph}_${chr}_ref_all_2_stdout \
-e ${path_gwas}post_imputation/2022/log/interval_ibdbioresource_${ph}_${chr}_ref_all_2_stderr \
"/path/to/software/./plink2 \
--pfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr${chr}_subset_included_in_${ph}_analysis_b38 \
--sort-vars --make-pgen --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr${chr}_subset_included_in_${ph}_analysis_b38"
done
done


for ph in ${pheno[@]}
do 
echo ${ph} && for chr in {1..22}
do
echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/interval_ibdbioresource_${ph}_${chr}_ref_all_2_stdout | grep "Successfully"
done
done


for ph in ${pheno[@]}
do 
for chr in {1..22} X
do
bsub -J"bgen_bed" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 4 \
-o ${path_gwas}post_imputation/2022/log/interval_ibdbioresource_${ph}_${chr}_ref_all_3_stdout \
-e ${path_gwas}post_imputation/2022/log/interval_ibdbioresource_${ph}_${chr}_ref_all_3_stderr \
"/path/to/software/./plink2 \
--pfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr${chr}_subset_included_in_${ph}_analysis_b38 \
--make-bed --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr${chr}_subset_included_in_${ph}_analysis_b38"
done
done

for ph in ${pheno[@]}
do 
echo ${ph} && for chr in {1..22}
do
echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/interval_ibdbioresource_${ph}_${chr}_ref_all_3_stdout | grep "Successfully"
done
done

for ph in ${pheno[@]}
do 
echo ${ph} && ls -la ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr*_subset_included_in_${ph}_analysis_b38.bed | wc -l
done




# merge all files into one

rm ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/bed_list_*.txt

for ph in ${pheno[@]}
do 
for chr in $(seq 2 22)
do
echo "${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr${chr}_subset_included_in_${ph}_analysis_b38.bed ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr${chr}_subset_included_in_${ph}_analysis_b38.bim ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr${chr}_subset_included_in_${ph}_analysis_b38.fam" \
>> ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/bed_list_${ph}.txt
done
done


MEM=30000 
# - needs more - running interactively on hpc-server with constraint
# Allocated 4686 MB successfully, after larger attempt(s) failed.

pheno=(ibd cd uc)

for ph in ${pheno[@]}
do
bsub -J"bgen_bed" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-o ${path_gwas}post_imputation/2022/log/interval_ibdbioresource_${ph}_bed_merge_stdout \
-e ${path_gwas}post_imputation/2022/log/interval_ibdbioresource_${ph}_bed_merge_stderr \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_chr1_subset_included_in_${ph}_analysis_b38 \
--merge-list ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/bed_list_${ph}.txt \
--memory ${MEM} \
--make-bed --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_allchr_subset_included_in_${ph}_analysis"
done

# ibd
# 17309347 variants and 1881 people pass filters and QC.

# cd
# 17345840 variants and 1128 people pass filters and QC.

# uc
# 17320086 variants and 1331 people pass filters and QC.

# MEM=20000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1


path_gwas="/path/to/ibdgwas/IIBDGC/"
path="/path/to/ibdgwas_bioresource/"

library(data.table)

pheno<-c("cd","uc")

for (ph in pheno) {
  
  print(ph)
  bim<-fread(paste(path_gwas,"post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_allchr_subset_included_in_",ph,"_analysis.bim",sep=""),head=F)
  
  bim$V2<-gsub("_",":",bim$V2)
  bim$V2<-paste("chr",bim$V2,sep="")
  
  write.table(bim,paste(path_gwas,"post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_allchr_subset_included_in_",ph,"_analysis_edited.bim",sep=""),
              col.names=F,row.names=F,quote=F,sep="\t")
  
}

q("no")

##############

MEM=30000 
# - needs more - running interactively on hpc-server with constraint
# Allocated 4686 MB successfully, after larger attempt(s) failed.

pheno=(ibd cd uc)

for ph in ${pheno[@]}
do
bsub -J"bgen_bed" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-o ${path_gwas}post_imputation/2022/log/interval_ibdbioresource_${ph}_bed_merge_stdout \
-e ${path_gwas}post_imputation/2022/log/interval_ibdbioresource_${ph}_bed_merge_stderr \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_allchr_subset_included_in_${ph}_analysis.bed \
--bim ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_allchr_subset_included_in_${ph}_analysis_edited.bim \
--fam ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_allchr_subset_included_in_${ph}_analysis.fam \
--make-bed --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_allchr_subset_included_in_${ph}_analysis_2"
done


##########################################
# 2.- Rerformat summary statistics

# --cojo-file test.ma
# Input the summary-level statistics from a meta-analysis GWAS (or a single GWAS).
# 
# Input file format
# test.ma
# SNP A1 A2 freq b se p N 
# rs1001 A G 0.8493 0.0024 0.0055 0.6653 129850 
# rs1002 C G 0.0306 0.0034 0.0115 0.7659 129799 
# rs1003 A C 0.5128 0.0045 0.0038 0.2319 129830
# ...
# Columns are SNP, the effect allele, the other allele, frequency of the effect allele, effect size, standard error, p-value and sample size. 
# The headers are not keywords and will be omitted by the program. Important: "A1" needs to be the effect allele with "A2" being the other allele and 
# "freq" should be the frequency of "A1".
# 
# Note: 1) For a case-control study, the effect size should be log(odds ratio) with its corresponding standard error. 2) Please always input the summary 
# statistics of all SNPs even if your analysis only focuses on a subset of SNPs because the program needs the summary data of all SNPs to calculate the 
# phenotypic variance. You can use one of the --extract options (Data management) to limit the COJO analysis in a certain genomic region.


echo "SNP A1 A2 freq b se p N" \
> ${path_gwas}post_imputation/2022/analysis/conditional_analysis/header.ma

path_gwas=/path/to/ibdgwas/IIBDGC/
pheno=(ibd cd uc)

for chr in {1..22} X
do
for ph in ${pheno[@]}
do 
zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/${chr}_${ph}_meta_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz \
| awk -F '\t' '{ print $1,$4,$3,$16,$5,$6,$7,$15}' | awk -F ' ' '$8 >= 1000 { print }' > \
${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/${chr}_${ph}_meta_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format.ma
done
done

for ph in ${pheno[@]}
do cat ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/*_${ph}_meta_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format.ma \
| grep -v "MarkerName" \
> ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/allchr_${ph}_meta_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format_tmp.ma
done

for ph in ${pheno[@]}
do 
cat ${path_gwas}post_imputation/2022/analysis/conditional_analysis/header.ma  ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/allchr_${ph}_meta_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format_tmp.ma \
> ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/allchr_${ph}_meta_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format.ma
done

for ph in ${pheno[@]}
do 
rm ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/allchr_${ph}_meta_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format_tmp.ma
done

for ph in ${pheno[@]}
do
sed -i 's/chrX:/chr23:/g' ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/allchr_${ph}_meta_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format.ma
done

for ph in ${pheno[@]}
do echo ${ph} && wc -l ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/allchr_${ph}_meta_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format.ma
done

# ibd
# 20443378 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/allchr_ibd_meta_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format.ma
# cd
# 16767157 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/allchr_cd_meta_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format.ma
# uc
# 20469725 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/sas/allchr_uc_meta_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format.ma

# --cojo-slct
# Perform a stepwise model selection procedure to select independently associated SNPs. Results will be saved in a *.jma file with additional file 
# *.jma.ldr showing the LD correlations between the SNPs.
# 
# --cojo-p 5e-8
# Threshold p-value to declare a genome-wide significant hit. The default value is 5e-8 if not specified. This option is only valid in conjunction 
# with the option --cojo-slct.

# keep only individuals that were included in the analysis in plink files


##############################
# INFORMED APPROACH:

# create file per chr to conditioning on known signals
# original data saved in spreadsheed: list_known_index_gwas_fm_2023 in qc_steps_summary_2023

# /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/forward_regression/list_known_index_gwas_fm_2023_with_iibdgc_tier1_eas_summary_stats.tsv
# list of variants in ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/cond_chr${chr}.snplist

MEM=20000
pheno=(ibd cd uc)

for ph in ${pheno[@]}
do
for chr in {1..22} X
do
bsub -J"cojo_i" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-o ${path_gwas}post_imputation/2022/log/cojo_sas_${ph}_${chr}_conditioning_known_ibd_stdout \
-e ${path_gwas}post_imputation/2022/log/cojo_sas_${ph}_${chr}_conditioning_known_ibd_stderr \
"/path/to/software/username/gcta-1.94.1-linux-kernel-3-x86_64/./gcta-1.94.1 \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_allchr_subset_included_in_${ph}_analysis_2 \
--chr ${chr} \
--cojo-file ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/allchr_${ph}_meta_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format.ma \
--cojo-cond ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/cond_chr${chr}.snplist \
--out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/${chr}_${ph}_independent_index_variants_conditioning_on_known_ibd_index"
done
done


for ph in ${pheno[@]}
do
echo ${ph} && for chr in {1..22}
do
echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/cojo_sas_${ph}_${chr}_conditioning_known_ibd_stdout | grep -E "Successfully|Exited"
done
done


### see how many variants are being used:

for ph in ${pheno[@]}
do
echo ${ph} && for chr in {1..22}
do
echo ${chr} && wc -l ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/cond_chr${chr}.snplist && \
tail -100 ${path_gwas}post_imputation/2022/log/cojo_sas_${ph}_${chr}_conditioning_known_ibd_stdout | grep -E "of them are matched"
done
done


##### RUN COJO NON-INFORMED - FIND ANY CONDITIONALLY SIGNIF VARIANT AT 5e-8

path_gwas=/path/to/ibdgwas/IIBDGC/
pheno=(ibd cd uc)
MEM=20000

for ph in ${pheno[@]}
do
for chr in {1..22} X
do
bsub -J"cojo_i" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-o ${path_gwas}post_imputation/2022/log/cojo_sas_${ph}_${chr}_no_conditioning_stdout \
-e ${path_gwas}post_imputation/2022/log/cojo_sas_${ph}_${chr}_no_conditioning__stderr \
"/path/to/software/username/gcta-1.94.1-linux-kernel-3-x86_64/./gcta-1.94.1 \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/interval_ibdbioresource_allchr_subset_included_in_${ph}_analysis_2 \
--chr ${chr} \
--cojo-file ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/allchr_${ph}_meta_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format.ma \
--cojo-slct \
--cojo-p 5e-8 \
--out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/sas/${chr}_${ph}_independent_index_variants_no_conditioning"
done
done

for ph in ${pheno[@]}
do
echo ${ph} && for chr in {1..22} X
do
echo ${chr} && ls -la ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${chr}_${ph}_independent_index_variants_no_conditioning.cma.cojo
done
done
done

