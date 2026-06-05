# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# double check overlap of sample ids:

# # singularity exec iibdgc_postprocess_10_singularity.sif

# MEM=15000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
library(R.utils)
library(ggplot2)

rm(list=ls())

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

ibd<-ibd[,c("SNP","CHR","BP","ALLELE2","ALLELE1","INFO","A1FREQ","N","rate_N")]

fwrite(ibd,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_ibd_cd_uc_uc_eur_tier2_list_variants_rate_0.5_info_0.9.tsv.gz",sep=""),
col.names=T,row.names=F,quote=F,sep="\t")

q("no")
