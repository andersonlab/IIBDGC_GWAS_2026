# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
############################################################################################################################################################
# RUN THE SAME ANALYSIS JUST INCLUDING HAPMAP3 VARIANTS:

# create the reference file to subset and compare:

# singularity exec iibdgc_postprocess_10_singularity.sif
# singularity exec polyfun_2_singularity.sif

path_gwas="/path/to/ibdgwas/IIBDGC/"

MEM=8000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q yesterday R

library(data.table)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

list_hapmap_variants<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/w_hm3.snplist",sep=""))

for (chr in c(1:22)) {

    tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/GRCh38/plink_files/1000G.EUR.hg38.",chr,".bim",sep=""))
    tmp<-tmp[which(tmp$V2 %in% list_hapmap_variants$SNP),]

    if(chr==1) {
        all<-tmp
    } else {
        all<-rbind(all,tmp)
    }
    rm(tmp)
}

# double check that Alele notation in bim does not reflect frequency
dim(all)
# [1] 1189841       6

all$SNP<-paste("chr",all$V1,":",all$V4,":",all$V6,":",all$V5,sep="")
all$X<-paste("chr",all$V1,":",all$V4,sep="")

list_snps<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier2_metaanalysis_rate_0.5_info_0.9",sep=""),head=F)
list_snps$X<-gsub(":[A-Z]{1}:[A-Z]{1}","",list_snps$V1)
colnames(list_snps)[1]<-"SNP"

dim(all[which(all$X %in% list_snps$X),])
# [1] 1181702       8

dim(all[which(all$SNP %in% list_snps$SNP),])
# [1] 832734      7

# subset the original list to the 1185395 variants present in both
dim(list_snps)
# [1] 12979489        4
list_snps<-list_snps[which(list_snps$X %in% all$X),]
dim(list_snps)
# [1] 1181711       4

# rename:
list_snps<-merge(list_snps,all[,c("X","V2")],by="X")
head(list_snps)
dim(list_snps)
[1] 1181711       4

# subset the summary statistics so they only those SNPs

path_gwas<-"/path/to/ibdgwas/IIBDGC/"
pheno<-c("ibd","cd","uc")


for (ph in pheno) {

    print(ph)
    tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_",ph,"_eur_tier2_list_variants_rate_0.5_info_0.9.tsv.gz",sep=""),head=T)
    tmp<-tmp[which(tmp$SNP %in% list_snps$SNP),]
    tmp<-merge(tmp,list_snps[,c("SNP","V2")],by="SNP")

    tmp<-tmp[,c("V2","CHR","BP","ALLELE2","ALLELE1","BETA","SE","P-value","INFO","N","A1FREQ")]
    colnames(tmp)<-c("SNP","CHR","BP","ALLELE2","ALLELE1","BETA","SE","P-value","INFO","N","A1FREQ")

    fwrite(tmp,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_",ph,"_eur_tier2_list_variants_rate_0.5_info_0.9_hm3.tsv.gz",sep=""),
    col.names=T,row.names=F,quote=F,sep="\t")
    rm(tmp)

}

q("no")


###################################################
# REFORMAT IBD SUMMARY STATISTICS - rg from polyfun/ldsc does not work use a combination of munge_polyfun_sumstats -> parquet to sumstats -> ldsc -rg (from ldsc)




pheno=(ibd cd uc)

for ph in ${pheno[@]}
do
echo ${ph} && zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_eur_tier2_list_variants_rate_0.5_info_0.9_hm3.tsv.gz | wc -l
done
# ibd
# 1180744
# cd
# 1181057
# uc
# 1181602

MEM=8000

for ph in ${pheno[@]}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${ph}_allchr_ldsc_sumstats_only_hm3_nohla_polyfun_mungen_format_stderr \
-o ${path_gwas}post_imputation/log/${ph}_allchr_ldsc_sumstats_only_hm3_nohla_polyfun_mungen_format_stdout \
"munge_polyfun_sumstats.py \
--sumstats ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${ph}_eur_tier2_list_variants_rate_0.5_info_0.9_hm3.tsv.gz \
--min-info 0 \
--min-maf 0 \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/${ph}_eur_tier2_list_variants_rate_0.5_info_0.9_hm3_nohla_sumstats_munged.parquet"
done

for ph in ${pheno[@]}
do
echo ${ph} && tail -30  ${path_gwas}post_imputation/log/${ph}_allchr_ldsc_sumstats_only_hm3_nohla_polyfun_mungen_format_stdout | grep "Successfully"
done

# CONTINUE HERE


MEM=10000

pheno=(ibd cd uc)
ppreval=(0.007 0.003 0.004)
spreval=(0.5 0.5 0.5)
# using now Neff, so the sample prevalence = 0.5

for i in {0..2}  
do echo ${pheno[${i}]} && echo ${ppreval[${i}]} && echo ${spreval[${i}]}
done


MEM=8000
for i in {0..2}  
do
bsub -J"h2_py2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/log/${pheno[${i}]}_eur_tier2_rate_0.5_info_0.9_allarrays_h2_ldsc_py2_stderr \
-o ${path_gwas}post_imputation/log/${pheno[${i}]}_eur_tier2_rate_0.5_info_0.9_allarrays_h2_ldsc_py2_stdout \
"ldsc.py \
--h2 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/${pheno[${i}]}_eur_tier2_list_variants_rate_0.5_info_0.9_hm3_nohla_sumstats_munged.parquet \
--w-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/GRCh38/weights/weights.hm3_noMHC. \
--ref-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/GRCh38/weights/weights.hm3_noMHC. \
--not-M-5-50 \
--samp-prev ${spreval[${i}]} \
--pop-prev ${ppreval[${i}]} \
--no-check-alleles \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/heritability/${pheno[${i}]}_eur_tier2_list_variants_rate_0.5_info_0.9_hm3_allarrays_heritability_estimation_liability_scale"
done

for i in {0..2}  
do echo ${pheno[${i}]} && cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/heritability/${pheno[${i}]}_eur_tier2_list_variants_rate_0.5_info_0.9_hm3_allarrays_heritability_estimation_liability_scale.log | \
grep -E "Total Liability scale h2|Lambda GC|Intercept|Ratio|SNPs remain"
done

# ibd
# After merging with reference panel LD, 1177805 SNPs remain.
# After merging with regression SNP LD, 1177805 SNPs remain.
# Total Liability scale h2: 0.1257 (0.0091)
# Lambda GC: 1.5617
# Intercept: 1.2118 (0.0167)
# Ratio: 0.1728 (0.0136)
# cd
# After merging with reference panel LD, 1178118 SNPs remain.
# After merging with regression SNP LD, 1178118 SNPs remain.
# Total Liability scale h2: 0.1685 (0.0153)
# Lambda GC: 1.4297
# Intercept: 1.143 (0.0139)
# Ratio: 0.1488 (0.0145)
# uc
# After merging with reference panel LD, 1178664 SNPs remain.
# After merging with regression SNP LD, 1178664 SNPs remain.
# Total Liability scale h2: 0.1332 (0.0095)
# Lambda GC: 1.4297
# Intercept: 1.1542 (0.0139)
# Ratio: 0.179 (0.0162)

#########################################################################################################################################

# rationale to keep lower freq variants, ie NealeLab using 1% in:
# https://nealelab.github.io/UKBB_ldsc/methods.html

# however to test this:

# keeping only variants with at different MAF thresholdss, still excluding HLA:


### CREATE MUNGE.PARQUET FILES

MEM=8000

for i in {0..2}  
do
echo ${pheno[${i}]} 
for maf in 0.005 0.01 0.05
do 
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${pheno[${i}]}_allchr_ldsc_sumstats_eur_tier2_rate_0.5_info_0.9_hm3_${maf}_nohla_mungen_format_stderr \
-o ${path_gwas}post_imputation/log/${pheno[${i}]}_allchr_ldsc_sumstats_eur_tier2_rate_0.5_info_0.9_hm3_${maf}_nohla_mungen_format_stdout \
"munge_polyfun_sumstats.py \
--sumstats ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_${pheno[${i}]}_eur_tier2_list_variants_rate_0.5_info_0.9_hm3.tsv.gz \
--min-info 0 \
--min-maf ${maf} \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/${pheno[${i}]}_eur_tier2_list_variants_rate_0.5_info_0.9_hm3_maf_${maf}_nohla_sumstats_munged.parquet"
done
done



for i in {0..2}  
do
echo ${pheno[${i}]} 
for maf in 0.005 0.01 0.05
do echo ${maf} && tail -50 ${path_gwas}post_imputation/log/${pheno[${i}]}_allchr_ldsc_sumstats_eur_tier2_rate_0.5_info_0.9_hm3_${maf}_nohla_mungen_format_stdout | grep "Successfully"
done
done

for i in {0..2}  
do
echo ${pheno[${i}]} 
for maf in 0.005 0.01 0.05
do echo ${maf} && ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/${pheno[${i}]}_eur_tier2_list_variants_rate_0.5_info_0.9_maf_${maf}_nohla_sumstats_munged.parquet
done
done

### ESTIMATE HERITABILITY

MEM=8000
for i in {0..2}  
do
echo ${pheno[${i}]} 
for maf in 0.005 0.01 0.05
do 
bsub -J"h2_py2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/log/${pheno[${i}]}_h2_ldsc_py2_common_0.005_${maf}_stderr \
-o ${path_gwas}post_imputation/log/${pheno[${i}]}_h2_ldsc_py2_common_0.005_${maf}_stdout \
"ldsc.py \
--h2 ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/${pheno[${i}]}_eur_tier2_list_variants_rate_0.5_info_0.9_hm3_maf_${maf}_nohla_sumstats_munged.parquet \
--w-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/GRCh38/weights/weights.hm3_noMHC. \
--ref-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/GRCh38/weights/weights.hm3_noMHC. \
--not-M-5-50 \
--samp-prev ${spreval[${i}]} \
--pop-prev ${ppreval[${i}]} \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/heritability/${pheno[${i}]}_eur_tier2_list_variants_rate_0.5_info_0.9_hm3_allarrays_maf_${maf}_heritability_estimation_liability_scale"
done
done

for i in {0..2}  
do 
echo ${pheno[${i}]} &&
for maf in 0.005 0.01 0.05
do echo ${maf} && cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/heritability/${pheno[${i}]}_eur_tier2_list_variants_rate_0.5_info_0.9_hm3_allarrays_maf_${maf}_heritability_estimation_liability_scale.log | \
grep -E "Total Liability scale h2|Lambda GC|Intercept|Ratio|SNPs remain"
done
done

# this output also feeds:
# ~/git/IIBDGC_GWAS/scripts/other/plot_heritability_liability_scale_different_thresholds.R