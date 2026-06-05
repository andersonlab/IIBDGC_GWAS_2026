# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# reformat partially harmonised files:

# # singularity exec iibdgc_postprocess_10_singularity.sif


# MEM=10000
# bsub -J"gsmr_pqtl" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 2 \
# -o ${path_gwas}post_imputation/2022/log/subset_gwas_variants_${gwas_id[i]}_stdout \
# -e ${path_gwas}post_imputation/2022/log/subset_gwas_variants_${gwas_id[i]}_stderr \
# "Rscript ~/git/IIBDGC_GWAS/scripts/other/reformat_partially_harmonized_gwas_catalog_beta_files.R ${gwas_id[i]} ${gwas_n[i]} > \
# ${path_gwas}post_imputation/2022/log/eformat_partially_harmonized_gwas_catalog_beta_files_${gwas_id[i]}.Rout"

# MEM=25000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
library(R.utils)
library(ggplot2)

rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

args = commandArgs(trailingOnly=TRUE)
study_id<-args[1]
N<-args[2]
path_file<-args[3]

# study_id<-"GCST90044763"
# N<-"400449"
# path_file<-"/path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/GCST90044763.h.tsv.gz"

print(study_id)
print(N)
print(path_file)

df<-fread(path_file)    
dim(df)

head(df[which(df$hm_code==13),])

# subset the files if allele frequency field is not missing:
if (nrow(df[which(is.na(df$effect_allele_frequency)),])!=nrow(df)) {
    df<-df[which(df$effect_allele_frequency>=0.001 & df$effect_allele_frequency<=0.999),]
    print(dim(df))
}


df$SNP<-paste("chr",df$chromosome,":",df$base_pair_location,":",df$other_allele,":",df$effect_allele,sep="")
df<-df[which(df$effect_allele!=df$other_allele),]
dim(df)
#

table(df$hm_code)

# hm_code	Number	Percentage	Explanation
# 10	Forward strand; Correct orientation; Already harmonised
# 11	Forward strand; Flipped orientation; Requires harmonisation
# 12	Reverse strand; Correct orientation; Already harmonised
# 13	Reverse strand; Flipped orientation; Requires harmonisation
# 5	Palindromic; Assume forward strand; Correct orientation; Already harmonised
# 6	Palindromic; Assume forward strand; Flipped orientation; Requires harmonisation

# if only odds_ratio collected:
if (nrow(df[is.na(df$beta)])==nrow(df)) {
    df$beta<-log(df$odds_ratio)
} else if (is.null(df$beta)) {
    df$beta<-log(df$odds_ratio)
}


# if no SE collected, but instead ci
if (nrow(df[is.na(df$standard_error)])==nrow(df)) {
    if (nrow(df[is.na(df$ci_upper)])!=nrow(df)) {
        df$SE<-(log(df$ci_upper)-log(df$ci_lower))/3.919928
    }
} else {
    df$SE<-df$standard_error
    if (all(df$SE==0)) {
        df$SE<-NA
    }
}

# update beta
df$BETA<-NA
df$BETA[which(df$hm_code %in% c(10,12,11,13))]<-df$beta[which(df$hm_code %in% c(10,12,11,13))]
# df$BETA[which(df$hm_code %in% c(11,13))]<-df$beta[which(df$hm_code %in% c(11,13))]*(-1) - no needed, all harmonised

# exclude the variants that were not harmonised, assumption on plus strand
df<-df[which(!df$hm_code %in% c(5,6)),]
dim(df)

# exclude variants with no beta
df<-df[which(!is.na(df$BETA)),]

# if no allelle frequency - keep NA and later add IIBDGC frequency:
if (nrow(df[which(is.na(df$effect_allele_frequency)),])==nrow(df)) {
    df$A1FREQ<-df$effect_allele_frequency
} else {
    df$A1FREQ<-NA
    df$A1FREQ[which(df$hm_code %in% c(10,12,11,13))]<-df$effect_allele_frequency[which(df$hm_code %in% c(10,12,11,13))]
    # df$A1FREQ[which(df$hm_code %in% c(11,13))]<-1-(df$effect_allele_frequency[which(df$hm_code %in% c(11,13))]) - no needed, all harmonised
}

df$N<-N
df$INFO<-1

summary(df$BETA)
summary(df$effect_allele_frequency)
summary(df$A1FREQ)


df<-df[,c("SNP","chromosome","base_pair_location","other_allele","effect_allele","BETA","SE","p_value","INFO","N","A1FREQ")]

colnames(df)<-c("SNP","CHR","BP","A0","A1","BETA","SE","PVALUE","INFO","N","A1FREQ")
dim(df)


# remove duplicated ids:
dups<-df$SNP[duplicated(df$SNP)]
length(dups)

df<-df[!df$SNP %in% dups,]
dim(df)



fwrite(df,paste(path_gwas,"resources/gwas_summary_statistics/",study_id,"_edited.tsv.gz",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

q("no")