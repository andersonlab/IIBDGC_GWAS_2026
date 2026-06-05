# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

# # singularity exec iibdgc_postprocess_10_singularity.sif
# MEM=4000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q cpu-interactive R

library(data.table)
library(plyr)
library(dplyr)

rm(list=ls())

# provide PHENOTYPE
args = commandArgs(trailingOnly=TRUE)

pheno<-args[1]
# pheno<-"cd"
# pheno<-"ibd"

print(pheno)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"
path_eqtl<-"/path/to/project"

# LOAD THE COLOCALIZATION RESULTS - SEE ~/git/IIBDGC_GWAS/scripts/other/coloc_subset_results_by_ld_leads.R
coloc_results<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_filtered_by_gwas_r2_0.6_allpheno_noFM_with_beta.tsv",sep=""),head=T)
coloc_results<-as.data.frame(coloc_results)
dim(coloc_results)
# 10420

# LOAD THE MAP FOR COLOCALIZATION RESULTS - SEE ~/git/snakemake_colocalisation/scripts/create_maps_for_IIBDGC.R
map<-fread("/path/to/project",head=T)
map<-as.data.frame(map)
map<-map[which(map$cohort!="edQTL_GTEX"),]
map<-map[which(map$quant_method %in% c("aptamer","ge","microarray")),]

# LOAD THE ENRICHMENT RESULTS - SEE ~/git/IIBDGC_GWAS/scripts/other/plot_enrichment_analyses.R
enrich<-fread("~/git/IIBDGC_GWAS/plots/paper_tables/Source_data_enrichment_analyses_eqtl_pqtl.tsv",head=T)
enrich<-as.data.frame(enrich)
enrich<-enrich[which(enrich$Category %in% map$id_map),]

enrich<-enrich[which(enrich$pheno==toupper(pheno)),]

dim(enrich)
# [1] 495  20

## only retain those singnificant
enrich<-enrich[which(enrich$Enrichment_p<=(0.05)),]
dim(enrich)
# [1] 491  20


# LOAD THE INDEPENDENT SIGNALS - see ~/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_38_2_define_independent_signals_using_as_ref_cojo_agnostic_results.R
all<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants.tsv.gz",sep=""))
all<-as.data.frame(all)
dim(all)
# [1] 631 182


# vec<-c("IIBDGC_GWAS_index_variant","MarkerName_eqtl_lead","gene","Beta_eqlt","SE_eqtl","Pvalue_eqtl","PP.H4.abf","cohort","study_label","condition_name",map$id_map)
# dat<-as.data.frame(matrix(nrow=1,ncol=length(vec)))
# colnames(dat)<-vec
# rm(vec)

dat_list <- vector("list", nrow(all))
for (i in seq_len(nrow(all))) {

    # IIBDGC_GWAS_index_variant in LD with local GWAS lead:
    tmp1<-coloc_results[coloc_results$IIBDGC_GWAS_index_variant.x==all$MarkerName[i] & coloc_results$R2_IIBDGC_GWAS_index_variant_MarkerName_eqtl_lead>=0.6, c("IIBDGC_GWAS_index_variant.x","qtl_lead","gene","gene_name","Beta_eqtl","SE_eqtl","Pvalue_eqtl","PP.H4.abf","cohort","study_label","condition_name","R2_IIBDGC_GWAS_index_variant_MarkerName_eqtl_lead")]

    # IIBDGC_GWAS_index_variant in LD with eQTL lead:
    tmp2<-coloc_results[coloc_results$IIBDGC_GWAS_index_variant.y==all$MarkerName[i] & coloc_results$R2_IIBDGC_GWAS_index_variant_MarkerName_gwas_local_lead>=0.6, c("IIBDGC_GWAS_index_variant.y","qtl_lead","PP.H4.abf","Beta_eqtl","SE_eqtl","Pvalue_eqtl","cohort","gene","gene_name","study_label","condition_name","R2_IIBDGC_GWAS_index_variant_MarkerName_gwas_local_lead")]

    colnames(tmp1)[c(1, 12)] <- c("IIBDGC_GWAS_index_variant", "R2")
    colnames(tmp2)[c(1, 12)] <- c("IIBDGC_GWAS_index_variant", "R2")

    tmp <- rbind(tmp1, tmp2)
    dat_list[[i]] <- tmp[!duplicated(tmp), ]

}
dat <- do.call(rbind, dat_list)
rm(dat_list)

# relabel pqtl id:
dat$gene_name[which(is.na(dat$gene_name))]<-dat$gene[which(is.na(dat$gene_name))]


dat1<-dat[,c("IIBDGC_GWAS_index_variant","gene_name")]
dat1<-dat1[!duplicated(dat1),]

final_list <- vector("list", nrow(dat1))
for (i in seq_len(nrow(dat1))) {

    tmp <- dat[dat$IIBDGC_GWAS_index_variant==dat1$IIBDGC_GWAS_index_variant[i] & dat$gene_name==dat1$gene_name[i], ]

    tmp_enrich <- as.data.frame(t(enrich[enrich$Category %in% tmp$condition_name, c("Category","Coefficient_zscore")]))
    colnames(tmp_enrich) <- tmp_enrich[1, ]
    tmp_enrich <- tmp_enrich[-1, , drop=FALSE]

    row <- cbind(dat1[i, ], tmp_enrich)
    final_list[[i]] <- as.data.frame(row)

}
dat_final <- rbind.fill(final_list)
rm(final_list)

dim(dat_final)
# [1] 652 410

dat_final$score<-do.call(pmax,c(dat_final[,3:ncol(dat_final)],na.rm=T))
dat_final$score<-as.numeric(dat_final$score)

summary(dat_final$score)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.3957  0.5777  2.0593  1.7415  2.9672  3.9476

fwrite(dat_final,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_matrices/matrix_",pheno,"_zscore_genes_tissue_noFM.tsv.gz"),col.names=T,row.names=F,quote=F,sep="\t")


q("no")