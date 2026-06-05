# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
 
# # singularity exec iibdgc_postprocess_10_singularity.sif
# MEM=8000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q cpu-interactive R

library(data.table)
library(rtracklayer)
library(stringr)
library(qvalue)
library(ggplot2)
library(ggarrange)

rm(list=ls())

# args<-commandArgs(trailingOnly=TRUE)
pheno<-c("ibd","cd","uc")

# for testing purposes

path_gwas<-"/path/to/ibdgwas/IIBDGC/"
path_eqtl<-"/path/to/project"


# LOAD THE MAP FOR COLOCALIZATION RESULTS - SEE ~/git/snakemake_colocalisation/scripts/create_maps_for_IIBDGC.R
map<-fread("/path/to/project",head=T)
map<-as.data.frame(map)
map<-map[which(map$cohort!="edQTL_GTEX"),]
map<-map[which(map$quant_method %in% c("aptamer","ge","microarray")),]

# LOAD THE INDEPENDENT SIGNALS - see ~/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_38_2_define_independent_signals_using_as_ref_cojo_agnostic_results.R
all<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_with_fm_signals_with_multiancestry_signals_with_pheno.tsv.gz",sep=""))
all<-as.data.frame(all)
dim(all)

# LOAD THE COLOCALIZATION RESULTS - SEE ~/git/IIBDGC_GWAS/scripts/other/evaluate_direction_effect_QTL_risk_variant.R

for (chr in c(1:22)) {
    for (ph in pheno) {
        tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/tmp_coloc_filtered_by_gwas_r2_0.6_allpheno_with_beta_with_wald_ratio_",chr,"_",ph,"_only_gwas_lead.tsv"),head=T)
        if (ph=="ibd" & chr==1) {
            dat<-tmp
        }else{
            dat<-rbind(dat,tmp,fill=T)
        }
    }
}
dat<-as.data.frame(dat)
table(dat$cohort[!grepl("ENSG",dat$gene)])
# pQTL_decode  pQTL_sparc 
#          27           6

dat$gene[!grepl("ENSG",dat$gene)]
dat$gene_name[which(dat$gene_name=="")]<-dat$gene[which(dat$gene_name=="")]

table(dat$pheno_coloc)
#   cd  ibd   uc 
# 6189 7028 5614

coloc_results<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_filtered_by_gwas_r2_0.6_allpheno_with_beta.tsv",sep=""))
coloc_results<-as.data.frame(coloc_results)

dat$X<-paste(dat$condition_name,dat$gene,dat$MarkerName_eqtl_lead,sep="_")
coloc_results$X<-paste(coloc_results$condition_name,coloc_results$gene,coloc_results$MarkerName_eqtl_lead,sep="_")

dim(dat[which(!dat$X %in% coloc_results$X),])
# 0

dim(coloc_results[which(!coloc_results$X %in% dat$X),])
# 0

length(names(table(dat$gene_name)))
# [1] 667


# for each gene, plot the wald_test results:
summary(dat$n_variants_in_ld_for_wr)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    1.00   15.00   42.00   70.85   99.00  516.00

list_genes<-names(table(dat$gene_name))
length(list_genes)
# [1] 667

dat<-merge(dat,map[,c("id_map","tissue_label","cell_label","condition_label","tissue_cell_condition_label","quant_method")],by.x="condition_name",by.y="id_map")
dim(dat)
# [1] 18831    25


dim(dat)
# [1] 18831    25

dim(dat[which(is.na(dat$n_variants_in_ld_for_wr)),])
# [1] 2835   25

dim(dat[which(!is.na(dat$n_variants_in_ld_for_wr)),])
# [1] 15996    27

2835/18831
# 0.1505496

# 1-0.1505496
# [1] 0.8494504

table(dat$study_label[which(is.na(dat$n_variants_in_ld_for_wr))])
# Allegbe_Harris_2025           BLUEPRINT        Fairfax_2014             GENCORD 
#                2490                   3                   1                   1 
#                GTEx             Hu_2021          Lepik_2017       Panousis_2025 
#                 164                  11                  36                 123 
#          Perez_2022      Schmiedel_2018             TwinsUK 
#                   1                   3                   2 

### compare with mean_beta from variants in LD:
table(dat$n_variants_in_ld_for_wr)
#     1 
# 15996 

for (chr in c(1:22)) {
    for (ph in pheno) {
        tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/tmp_coloc_filtered_by_gwas_r2_0.6_allpheno_with_beta_with_wald_ratio_",chr,"_",ph,".tsv"),head=T)
        if (ph=="ibd" & chr==1) {
            dat_r2<-tmp
        }else{
            dat_r2<-rbind(dat_r2,tmp,fill=T)
        }
    }
}
dat_r2<-as.data.frame(dat_r2)
dat_r2$gene_name[which(dat_r2$gene_name=="")]<-dat_r2$gene[which(dat_r2$gene_name=="")]

dat_r2<-merge(dat_r2,map[,c("id_map","tissue_label","cell_label","condition_label","tissue_cell_condition_label","quant_method")],by.x="condition_name",by.y="id_map")
dim(dat_r2)

dat_r2$X<-paste(dat_r2$study_label,dat_r2$condition_name,dat_r2$gene,dat_r2$MarkerName_eqtl_lead,dat_r2$pheno_coloc,dat_r2$PP.H4.abf,dat_r2$R2,dat_r2$IIBDGC_GWAS_index_variant,sep="_")
dat$X<-paste(dat$study_label,dat$condition_name,dat$gene,dat$MarkerName_eqtl_lead,dat$pheno_coloc,dat$PP.H4.abf,dat$R2,dat$IIBDGC_GWAS_index_variant,sep="_")

dat$X[which(!dat$X %in% dat_r2$X)]
dat_r2$X[which(!dat_r2$X %in% dat$X)]


tmp<-merge(dat[,c("X","wald_ratio","n_variants_in_ld_for_wr")],dat_r2[,c("X","mean_wald_ratio","n_variants_in_ld_for_wr","betas_wald_ratio")],by="X")
cor.test(tmp$wald_ratio,tmp$mean_wald_ratio)
#       cor 
# 0.9654505 , p-value < 2.2e-16

dim(tmp[which(!is.na(tmp$wald_ratio) & !is.na(tmp$mean_wald_ratio)),])
# [1] 15996     5

dim(tmp[which( (tmp$wald_ratio>0 & tmp$mean_wald_ratio<0) |  (tmp$wald_ratio<0 & tmp$mean_wald_ratio>0)),])
# [1] 21  5

# clearly signals with + and negative variants


q("no")
