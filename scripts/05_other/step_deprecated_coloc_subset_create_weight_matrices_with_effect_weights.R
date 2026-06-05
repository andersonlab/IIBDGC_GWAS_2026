# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

# # singularity exec iibdgc_postprocess_10_singularity.sif
# MEM=4000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q cpu-interactive R

library(data.table)
library(rtracklayer)
library(stringr)
library(rtracklayer)
library(qvalue)
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
coloc_results<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_filtered_by_gwas_r2_0.6_allpheno_with_beta.tsv",sep=""),head=T)
coloc_results<-coloc_results[which(coloc_results$pheno_coloc==pheno),]

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
enrich<-enrich[which(enrich$Enrichment_p<=(0.05)),]
dim(enrich)


## only retain those singificant

# LOAD THE INDEPENDENT SIGNALS - see ~/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_38_2_define_independent_signals_using_as_ref_cojo_agnostic_results.R
all<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_with_fm_signals_with_multiancestry_signals.tsv.gz",sep=""))
all<-as.data.frame(all)
dim(all)
# [1] 579 218

# vec<-c("IIBDGC_GWAS_index_variant","MarkerName_eqtl_lead","gene","Beta_eqlt","SE_eqtl","Pvalue_eqtl","PP.H4.abf","cohort","study_label","condition_name",map$id_map)
# dat<-as.data.frame(matrix(nrow=1,ncol=length(vec)))
# colnames(dat)<-vec
# rm(vec)

for (i in c(1:nrow(all))) {

    # IIBDGC_GWAS_index_variant in LD with local GWAS lead:
    tmp1<-coloc_results[which(coloc_results$IIBDGC_GWAS_index_variant.x==all$MarkerName[i]),c("IIBDGC_GWAS_index_variant.x","qtl_lead","gene","gene_name","Beta_eqtl","SE_eqtl","Pvalue_eqtl","PP.H4.abf","cohort","study_label","condition_name")]

    # IIBDGC_GWAS_index_variant in LD with eQTL lead:
    tmp2<-coloc_results[which(coloc_results$IIBDGC_GWAS_index_variant.y==all$MarkerName[i]),c("IIBDGC_GWAS_index_variant.y","qtl_lead","PP.H4.abf","Beta_eqtl","SE_eqtl","Pvalue_eqtl","cohort","gene","gene_name","study_label","condition_name")]

    colnames(tmp1)[1]<-"IIBDGC_GWAS_index_variant"
    colnames(tmp2)[1]<-"IIBDGC_GWAS_index_variant"

    tmp<-rbind(tmp1,tmp2)
    tmp<-tmp[!duplicated(tmp),]

    if (i==1) {
        dat<-tmp
    } else {
        dat<-rbind(dat,tmp)
    }

}

# relabel pqtl id:
dat$gene_name[which(is.na(dat$gene_name))]<-dat$gene[which(is.na(dat$gene_name))]

dat1<-dat[,c("IIBDGC_GWAS_index_variant","gene_name")]
dat1<-dat1[!duplicated(dat1),]


vec<-c(colnames(dat1),map$id_map)
dat_final<-as.data.frame(matrix(ncol=length(vec),nrow=0))
colnames(dat_final)<-vec

for (i in 1:nrow(dat1)) {


    tmp<-dat[which(dat$IIBDGC_GWAS_index_variant==dat1$IIBDGC_GWAS_index_variant[i] & dat$gene_name==dat1$gene_name[i]),]

    tmp_enrich<-as.data.frame(t(enrich[which(enrich$Category %in% tmp$condition_name),c("Category","Coefficient_zscore")]))
    colnames(tmp_enrich)<-tmp_enrich[1,]
    tmp_enrich<-tmp_enrich[-1,,drop=F]

    tmp<-dat1[i,]
    tmp<-cbind(tmp,tmp_enrich)
    tmp<-as.data.frame(tmp)

    dat_final<-rbind.fill(dat_final,tmp)
}

dat_final$score<-do.call(pmax,c(dat_final[,3:ncol(dat_final)],na.rm=T))
dat_final$score<-as.numeric(dat_final$score)

summary(dat_final$score)


if (pheno=="ibd") {
    dat_final<-dat_final[which(dat_final$IIBDGC_GWAS_index_variant %in% all$MarkerName[which(all$pheno %in% c("IBD_unsaturated","IBD_saturated","CD","UC"))]),]
    all<-all[,c("MarkerName","BETA_ibd_eur_tier_2","SE_ibd_eur_tier_2")]
} else if (pheno=="cd") {
    dat_final<-dat_final[which(dat_final$IIBDGC_GWAS_index_variant %in% all$MarkerName[which(all$pheno %in% c("IBD_unsaturated","IBD_saturated","CD"))]),]
    all<-all[,c("MarkerName","BETA_cd_eur_tier_2","SE_cd_eur_tier_2")]
} else if (pheno=="uc") {
    dat_final<-dat_final[which(dat_final$IIBDGC_GWAS_index_variant %in% all$MarkerName[which(all$pheno %in% c("IBD_unsaturated","IBD_saturated","UC"))]),]
    all<-all[,c("MarkerName","BETA_uc_eur_tier_2","SE_uc_eur_tier_2")]
}

colnames(all)<-c("MarkerName","BETA","SE")
all$variant_score<-abs(all$BETA/all$SE)

### weight those by beta/se effect sizes:

dat_final<-merge(dat_final,all,by.x="IIBDGC_GWAS_index_variant",by.y="MarkerName",all.x=T)
dat_final$score_tissue_variant<-dat_final$score*dat_final$variant_score

fwrite(dat_final,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_matrices/matrix_",pheno,"_zscore_genes_tissue_variant.tsv.gz"),col.names=T,row.names=F,quote=F,sep="\t")


# Write output to test perfornamnce of the socre

dat_final<-dat_final[!grepl("ENSG",dat_final$gene_name),]
dat_final<-dat_final[which(!is.na(dat_final$score_tissue_variant)),]


dat_final<-dat_final %>% mutate(group = ntile(score_tissue_variant, 2))

for (i in 1:2) {

    write.table(dat_final[which(dat_final$group==i),"gene_name",drop=F],paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/top_",i,"_",pheno,"_zscore_genes_tissue_variant"),col.names=F,row.names=F,sep="\t")

}

q("no")