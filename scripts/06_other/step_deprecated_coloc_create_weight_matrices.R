# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# # singularity exec iibdgc_postprocess_10_singularity.sif
# MEM=6000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q cpu-interactive R

library(data.table)
library(rtracklayer)
library(stringr)
library(rtracklayer)
library(qvalue)

rm(list=ls())


path_gwas<-"/path/to/ibdgwas/IIBDGC/"

coloc_results<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_filtered_by_gwas_r2_0.6_allpheno.tsv",sep=""),head=T)
coloc_results<-as.data.frame(coloc_results)

all<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_with_fm_signals_no_overlap.tsv.gz",sep=""))
all<-as.data.frame(all)
dim(all)
# [1] 554 218

dim(all[which(all$MarkerName %in% coloc_results$IIBDGC_GWAS_index_variant),])
# [1] 138 211


for (i in 1:nrow(all)) {

    tmp<-co

}