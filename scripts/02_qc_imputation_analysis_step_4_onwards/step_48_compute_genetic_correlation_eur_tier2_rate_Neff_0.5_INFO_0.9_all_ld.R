# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# singularity exec iibdgc_postprocess_10_singularity.sif
# singularity exec polyfun_2_singularity.sif

path_gwas="/path/to/ibdgwas/IIBDGC/"

######################################################################################
# 1.- create ld regression files with only SNPs - rg function cannot handle non-SNPs

MEM=12000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q yesterday R

library(data.table)
library(R.utils)
library(ggplot2)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

ibd<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_ibd_eur_tier2_list_variants_rate_0.5_info_0.9.tsv.gz",sep=""),head=T)
cd<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_cd_eur_tier2_list_variants_rate_0.5_info_0.9.tsv.gz",sep=""),head=T)
cd<-cd[which(!cd$SNP %in% ibd$SNP),]
ibd<-rbind(ibd,cd)
rm(cd)
uc<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_uc_eur_tier2_list_variants_rate_0.5_info_0.9.tsv.gz",sep=""),head=T)
uc<-uc[which(!uc$SNP %in% ibd$SNP),]
ibd<-rbind(ibd,uc)
rm(uc)

dim(ibd)
# [1] 12979489       11

ibd<-ibd[which(nchar(ibd$ALLELE2)==1 & nchar(ibd$ALLELE1)==1), ]
dim(ibd)
# [1] 12132640       12

table(ibd$ALLELE2,ibd$ALLELE1)
#           A       C       G       T
#   A       0  445147 1757072  411276
#   C  536262       0  519102 2397748
#   G 2401060  517284       0  537009
#   T  410020 1755025  445635       0

fwrite(ibd[,c("SNP")],paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_ibd_cd_uc_eur_tier2_list_variants_rate_0.5_info_0.9_only_SNPs",sep=""),
col.names=T,row.names=F,quote=F,sep="\t")

fwrite(ibd[,c("SNP","ALLELE1","ALLELE2","INFO","A1FREQ")],paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_ibd_cd_uc_eur_tier2_list_variants_rate_0.5_info_0.9_only_SNPs.tsv.gz",sep=""),
col.names=T,row.names=F,quote=F,sep="\t")


ibd<-ibd[,c("SNP","ALLELE1","ALLELE2")]
colnames(ibd)<-c("SNP","A1","A2")

fwrite(ibd,paste(path_gwas,"resources/gwas_summary_statistics/subset_allchr_ibd_cd_uc_eur_tier2_variants_for_genetic_correlation.snplist",sep=""),
col.names=T,row.names=F,quote=F,sep="\t")

q("no")


##################################################################################################################################################################################
## CREATE .BGEN FILES FROM GSA USING SAME FILES (and N samples) USED FOR REGENIE STEP2, BUT CONVERTED AND SUBSET BY SNPs INCLUDED IN FINAL ANALYSIS:

MEM=22000

for chr in {1..22}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 4 \
-e ${path_gwas}post_imputation/log/${chr}_ldsc_ldsc_plink_stderr \
-o ${path_gwas}post_imputation/log/${chr}_ldsc_ldsc_plink_stdout \
"plink2 \
--bfile ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/allarrays_chr${chr}_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_rate_0.5_info_0.9 \
--extract ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_ibd_cd_uc_eur_tier2_list_variants_rate_0.5_info_0.9_only_SNPs \
--threads 4 \
--make-bed --out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/allarrays_chr${chr}_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_rate_0.5_info_0.9_only_SNPs"
done

for chr in {1..22}
do
echo ${chr} && tail -30 ${path_gwas}post_imputation/log/${chr}_ldsc_ldsc_plink_stdout  | grep -E "Successfully|Exit"
done


for chr in {1..22}
do
ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/allarrays_chr${chr}_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_rate_0.5_info_0.9_only_SNPs.bim
done

wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_ibd_cd_uc_eur_tier2_list_variants_rate_0.5_info_0.9_only_SNPs
# 12132641 - with header

wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/allarrays_chr*_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_rate_0.5_info_0.9_only_SNPs.bim
# 12121043

### update reference files

MEM=12000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q yesterday R

library(data.table)
library(R.utils)
library(ggplot2)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

for (chr in 1:22) {
    tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/allarrays_chr",chr,"_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_rate_0.5_info_0.9_only_SNPs.bim",sep=""),head=F)
    if(chr==1) {
        bim<-tmp
    } else {
        bim<-rbind(bim,tmp)
    }
    rm(tmp)
}

dim(bim)
# 12121043

ibd<-fread(paste(path_gwas,"resources/gwas_summary_statistics/subset_allchr_ibd_cd_uc_eur_tier2_variants_for_genetic_correlation.snplist",sep=""))
dim(ibd)
# 12132640

ibd<-ibd[which(ibd$SNP %in% bim$V2),]
dim(ibd)
# [1] 12121043        3

fwrite(ibd[,c("SNP")],paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_ibd_cd_uc_eur_tier2_list_variants_rate_0.5_info_0.9_only_SNPs",sep=""),
col.names=T,row.names=F,quote=F,sep="\t")


fwrite(ibd,paste(path_gwas,"resources/gwas_summary_statistics/subset_allchr_ibd_cd_uc_eur_tier2_variants_for_genetic_correlation.snplist",sep=""),
col.names=T,row.names=F,quote=F,sep="\t")

q("no")

### estimate regression weights, no partitioned:

# MEM=60000
# for chr in 2
# MEM=55000
# for chr in 1 {3..6} 
# MEM=45000
# for chr in {7..8}
# MEM=40000
# for chr in {9..12}
# MEM=30000
# for chr in {13..16}
# MEM=25000
# for chr in {17..19}
# MEM=20000
# for chr in {20..22}
do
bsub -J"ldsc_w_SNP" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q week \
-e ${path_gwas}post_imputation/log/${chr}_ldsc_w_stderr \
-o ${path_gwas}post_imputation/log/${chr}_ldsc_w_stdout \
"ldsc.py \
--bfile ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_april2025/allarrays_chr${chr}_union_variants_cd_uc_ibd_eur_tier_2_metaanalysis_rate_0.5_info_0.9_only_SNPs \
--l2 \
--ld-wind-kb 1000 \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/regression_weight_files/weights_only_SNPs.april2025.${chr}"
done



for chr in {1..22}
do
echo ${chr} && tail -30 ${path_gwas}post_imputation/log/${chr}_ldsc_w_stdout  | grep -E "Successfully|Exit"
done


##################################################################################################################################
# subset the summary statistics so they only include SNPs:

# The ldsc .sumstats format requires six pieces of information for each SNP:
# A unique identifier (e.g., the rs number)
# Allele 1 (effect allele)
# Allele 2 (non-effect allele)
# Sample size (which often varies from SNP to SNP)
# A P-value
# A signed summary statistic (beta, OR, log odds, Z-score, etc)

MEM=12000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q yesterday R

library(data.table)
library(R.utils)
library(ggplot2)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"
pheno<-c("ibd","cd","uc")

list_snps<-fread(paste(path_gwas,"resources/gwas_summary_statistics/subset_allchr_ibd_cd_uc_eur_tier2_variants_for_genetic_correlation.snplist",sep=""),head=T)
dim(list_snps)
# [1] 12121043        3
for (ph in pheno) {

    print(ph)

    tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_",ph,"_eur_tier2_list_variants_rate_0.5_info_0.9.tsv.gz",sep=""),head=T)
    print(dim(tmp))
    tmp<-tmp[which(tmp$SNP %in% list_snps$SNP),]
    print(dim(tmp))

    tmp<-tmp[,c("SNP","BP","ALLELE2","ALLELE1","BETA","SE","P-value","INFO","A1FREQ","N","CHR")]

    fwrite(tmp,paste(path_gwas,"resources/gwas_summary_statistics/allchr_",ph,"_eur_tier2_list_variants_rate_0.5_info_0.9_only_SNPs.tsv.gz",sep=""),
    col.names=T,row.names=F,quote=F,sep="\t")
    rm(tmp)

}

q("no")



###################################################
# REFORMAT IBD SUMMARY STATISTICS - rg from polyfun/ldsc does not work use a combination of munge_polyfun_sumstats -> parquet to sumstats -> ldsc -rg (from ldsc)

# path_gwas="/path/to/ibdgwas/IIBDGC/"
path_gwas="/path/to/ibdgwas/IIBDGC/"

# version for polyfun/munge
# singularity exec polyfun_2_singularity.sif


pheno=(ibd cd uc)

for ph in ${pheno[@]}
do 
echo ${ph} &&  zcat ${path_gwas}resources/gwas_summary_statistics/allchr_${ph}_eur_tier2_list_variants_rate_0.5_info_0.9_only_SNPs.tsv.gz | wc -l
done 
# ibd
# 11845840
# cd
# 11991582
# uc
# 11989744

MEM=12000

for ph in ${pheno[@]}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G humgen-priority \
-e ${path_gwas}post_imputation/log/${ph}_allchr_ldsc_sumstats_only_SNPs_nohla_polyfun_mungen_format_stderr \
-o ${path_gwas}post_imputation/log/${ph}_allchr_ldsc_sumstats_only_SNPs_nohla_polyfun_mungen_format_stdout \
"munge_polyfun_sumstats.py \
--sumstats ${path_gwas}resources/gwas_summary_statistics/allchr_${ph}_eur_tier2_list_variants_rate_0.5_info_0.9_only_SNPs.tsv.gz \
--min-info 0 \
--min-maf 0 \
--remove-strand-ambig \
--out ${path_gwas}resources/gwas_summary_statistics/${ph}_eur_tier2_list_variants_rate_0.5_info_0.9_only_SNPs_nohla_sumstats_munged.parquet"
done



less ${path_gwas}post_imputation/log/ibd_allchr_ldsc_sumstats_only_SNPs_nohla_polyfun_mungen_format_stdout 
# [INFO]  Reading sumstats file...
# [INFO]  Done in 49.23 seconds
# [INFO]  11845839 SNPs are in the sumstats file
# [INFO]  Removing 1813877 SNPs with strand ambiguity
# [INFO]  Removing 55100 HLA SNPs
# [INFO]  9985018 SNPs with sumstats remained after all filtering stages
# [INFO]  Saving munged sumstats of 9985018 SNPs to /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/ibd_eur_ti
# er2_list_variants_rate_0.5_info_0.9_only_SNPs_nohla_sumstats_munged.parquet
# [INFO]  Done

less ${path_gwas}post_imputation/log/cd_allchr_ldsc_sumstats_only_SNPs_nohla_polyfun_mungen_format_stdout 
# [INFO]  Reading sumstats file...
# [INFO]  Done in 38.91 seconds
# [INFO]  11991581 SNPs are in the sumstats file
# [INFO]  Removing 1836361 SNPs with strand ambiguity
# [INFO]  Removing 55584 HLA SNPs
# [INFO]  10107872 SNPs with sumstats remained after all filtering stages
# [INFO]  Saving munged sumstats of 10107872 SNPs to /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/cd_eur_ti
# er2_list_variants_rate_0.5_info_0.9_only_SNPs_nohla_sumstats_munged.parquet
# [INFO]  Done

less ${path_gwas}post_imputation/log/uc_allchr_ldsc_sumstats_only_SNPs_nohla_polyfun_mungen_format_stdout 
# [INFO]  Reading sumstats file...
# [INFO]  Done in 38.68 seconds
# [INFO]  11989743 SNPs are in the sumstats file
# [INFO]  Removing 1835833 SNPs with strand ambiguity
# [INFO]  Removing 55163 HLA SNPs
# [INFO]  10106917 SNPs with sumstats remained after all filtering stages
# [INFO]  Saving munged sumstats of 10106917 SNPs to /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/uc_eur_ti
# er2_list_variants_rate_0.5_info_0.9_only_SNPs_nohla_sumstats_munged.parquet
# [INFO]  Done



# SAVE PARQUET FILE AS FLAT TEXT FIEL AS INPUT FOR LDSC/LDSC.PY


MEM=6000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q yesterday R

library(arrow)
library(data.table)

pheno<-c("ibd","cd","uc")

for (ph in pheno) {

    print(ph)

    df<-read_parquet(paste("/path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/",ph,"_eur_tier2_list_variants_rate_0.5_info_0.9_only_SNPs_nohla_sumstats_munged.parquet",sep=""))
    print(head(df))

    print(dim(df))

    fwrite(df[,2:9],paste("/path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/",ph,"_eur_tier2_list_variants_rate_0.5_info_0.9_only_SNPs_nohla_sumstats_munged.sumstats",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    rm(df)

}

# [1] "ibd"
#   index                SNP        BP A2 A1        MAF      N CHR          Z
# 1     0 chr1:153913869:G:A 153913869  G  A 0.01925376 240386   1  0.7905983
# 2     1  chr1:40462393:G:A  40462393  G  A 0.02262502 216515   1 -1.8523207
# 3     4 chr1:109049812:C:T 109049812  C  T 0.04416600 240386   1 -0.1543210
# 4     6 chr1:212561652:A:G 212561652  A  G 0.26622786 240386   1 -0.5945946
# 5     7 chr1:117432712:C:T 117432712  C  T 0.14571169 240386   1  1.3191489
# 6     8  chr1:22652153:G:A  22652153  G  A 0.27606190 240386   1  1.7222222
# [1] 9985018       9
# [1] "cd"                                                                                                                             
#   index                SNP        BP A2 A1        MAF      N CHR           Z
# 1     0 chr1:153913869:G:A 153913869  G  A 0.01920181 118922   1  1.31306991
# 2     1  chr1:40462393:G:A  40462393  G  A 0.02285365 115451   1 -1.24117647
# 3     4 chr1:109049812:C:T 109049812  C  T 0.04395874 118922   1  0.03097345
# 4     6 chr1:212561652:A:G 212561652  A  G 0.26625780 118922   1 -0.60576923
# 5     7 chr1:117432712:C:T 117432712  C  T 0.14599057 118922   1  0.45112782
# 6     8  chr1:22652153:G:A  22652153  G  A 0.27639155 118922   1  0.96039604
# [1] 10107872        9
# [1] "uc"                                                                                                                             
#   index                SNP        BP A2 A1        MAF      N CHR          Z
# 1     0 chr1:153913869:G:A 153913869  G  A 0.01907813 145099   1  0.1280277
# 2     1  chr1:40462393:G:A  40462393  G  A 0.02295090 134977   1 -0.9413793
# 3     4 chr1:109049812:C:T 109049812  C  T 0.04405045 145099   1 -0.4572864
# 4     6 chr1:212561652:A:G 212561652  A  G 0.26656346 145099   1 -0.1208791
# 5     7 chr1:117432712:C:T 117432712  C  T 0.14638562 145099   1  1.1565217
# 6     8  chr1:22652153:G:A  22652153  G  A 0.27612573 145099   1  1.1931818
# [1] 10106917        9

q("no")


##################################################################################################################################

# CONTINUE HERE - part 2



# works with sumstats file + lsdsc from ldsc:

module unload HGI/softpack/groups/team152/polyfun-2/1
# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=7500

pheno=(ibd cd uc)

for j in {0..2}
do 
for i in {0..2}
do 
bsub -J"ldsc_pf" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${pheno[i]}_${pheno[j]}_ldsc_only_SNPs_pf_mungen_format_sumstats_ldsc_script_stderr \
-o ${path_gwas}post_imputation/log/${pheno[i]}_${pheno[j]}_ldsc_only_SNPs_pf_mungen_format_sumstats_ldsc_script_stdout \
"ldsc.py \
--rg ${path_gwas}resources/gwas_summary_statistics/${pheno[i]}_eur_tier2_list_variants_rate_0.5_info_0.9_only_SNPs_nohla_sumstats_munged.sumstats,${path_gwas}resources/gwas_summary_statistics/${pheno[j]}_eur_tier2_list_variants_rate_0.5_info_0.9_only_SNPs_nohla_sumstats_munged.sumstats \
--no-check-alleles \
--w-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/regression_weight_files/weights_only_SNPs.april2025. \
--ref-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/regression_weight_files/weights_only_SNPs.april2025. \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/genetic_correlation/${pheno[i]}_${pheno[j]}_gc"
done
done


##########################################################################################################################################################
##########################################################################################################################################################
# REFORMAT OTHER SUMMARY STATISTICS FROM GWAS CATALOG - reformat data to munge.sumstats format:

~/git/IIBDGC_GWAS/scripts/other/reformat_other_gwas_summary_statistics_to_munge_snpstats.R


##########################################################################################################################################################
##########################################################################################################################################################
# run the analysis with the other gwas summary stats:

# list the studies:

files=($(cat /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/list_path_summary_stats_gwas_catalog)) 
echo ${#files[@]}
# 281

# GCST90129505 - Fernandez rozadilla 2023 to be added
# GCST90029070 CRP to be added
# GCST90165267 clonal haematopoyesis added
# GCST003566 Multiple sclerosis
# GCST90476007 venous embolism and thrombosis

gwas_id=($(cat /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/list_path_summary_stats_gwas_catalog <(echo GCST90129505) <(echo GCST90029070) <(echo GCST90165267) <(echo GCST003566) <(echo GCST90476007)| \
sed 's/\/lustre\/scratch124\/humgen\/projects_v2\/ibdgwas\/IIBDGC\/resources\/gwas_summary_statistics\///g' | \
sed 's/\/lustre\/scratch125\/humgen\/resources_v2\/GWAScatalog\/.*\///g' | sed  's/.h.tsv.gz//g' && ls /path/to/project | \
sed 's/-/_/g' | sed 's/(/_/g' | sed 's/)/_/g' | grep ".tsv.gz$" | sed 's/.tsv.gz//g')) 
echo ${#gwas_id[@]}
# 349

for i in {0..348}
do 
ls -la ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_only_SNPs_sumstats_munged.sumstats
done

module unload HGI/softpack/groups/team152/polyfun-2/1
# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=7500

pheno=(ibd cd uc)

for j in {0..2}
do 
for i in {0..348}
do 
echo ${pheno[j]}_${gwas_id[i]}
done
done


for j in {0..2}
do 
for i in  {0..348}
do 
bsub -J"ldsc_pf" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G humgen-priority \
-e ${path_gwas}post_imputation/log/${pheno[j]}_${gwas_id[i]}_ldsc_only_SNPs_pf_mungen_format_sumstats_ldsc_script_stderr \
-o ${path_gwas}post_imputation/log/${pheno[j]}_${gwas_id[i]}_ldsc_only_SNPs_pf_mungen_format_sumstats_ldsc_script_stdout \
"ldsc.py \
--rg ${path_gwas}resources/gwas_summary_statistics/${pheno[j]}_eur_tier2_list_variants_rate_0.5_info_0.9_only_SNPs_nohla_sumstats_munged.sumstats,${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_only_SNPs_sumstats_munged.sumstats \
--no-check-alleles \
--w-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/regression_weight_files/weights_only_SNPs.april2025. \
--ref-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/regression_weight_files/weights_only_SNPs.april2025. \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/genetic_correlation/${pheno[j]}_eur_tier2_list_variants_rate_0.5_info_0.9_${gwas_id[i]}_gc"
done
done


rm ~/tmp.txt

for j in {0..2}
do 
for i in {0..348}
do echo ${pheno[j]}_${gwas_id[i]} >> ~/tmp.txt && tail -50 ${path_gwas}post_imputation/log/${pheno[j]}_${gwas_id[i]}_ldsc_only_SNPs_pf_mungen_format_sumstats_ldsc_script_stdout | grep -E "Successfully|Exited" >> ~/tmp.txt
done
done

rm ~/tmp.txt
for i in {0..348}
do echo ${i} >> ~/tmp.txt && echo ${gwas_id[i]} >> ~/tmp.txt
done



# exclude:

# files that fail:

# GCST002318_edited.tsv.gz - 253 Rheumatoid arthritis no beta, no SE, no allele frequency - EXCLUDE
# GCST90044763_edited.tsv.gz - 269 Eczema no beta or se
# GCST90480502 - h2's was out of bounds - duplicated  AS, exclude
# GCST003566 - multiple sclerosis - 284 duplicated trait 
# GCST003156 - lupus - 259, duplicated trait
# GCST90014023 - type 1 - 256 diabetes, duplicated trait
GCST90044158 - celiac disease - 268 duplicated
GCST005527 - psoriasis - 261 duplicated trait

rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/genetic_correlation/*_${gwas_id[i]}_gc.log


# plot the results
IIBDGC_GWAS/scripts/other/plot_genetic_correlation_results.R


