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

rm(list=ls())

args<-commandArgs(trailingOnly=TRUE)

chr<-args[1]
ph<-args[2]

# for testing purposes
# chr<-1
# ph<-"ibd"


path_gwas<-"/path/to/ibdgwas/IIBDGC/"
path_eqtl<-"/path/to/project"


# LOAD THE INDEPENDENT SIGNALS - see ~/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_38_2_define_independent_signals_using_as_ref_cojo_agnostic_results.R
all<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants.tsv.gz",sep=""))
all<-as.data.frame(all)
dim(all)

all<-all[which(all$chr==chr),]
dim(all)

# LOAD THE COLOCALIZATION RESULTS - SEE ~/git/IIBDGC_GWAS/scripts/other/coloc_subset_results_by_ld_leads.R
coloc_results<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_filtered_by_gwas_r2_0.6_allpheno_noFM_with_beta.tsv",sep=""),head=T)
coloc_results<-as.data.frame(coloc_results)
coloc_results<-coloc_results[which(coloc_results$chr==chr),]
coloc_results<-coloc_results[which(coloc_results$pheno_coloc==ph),]
dim(coloc_results)


# COMBINE BOTH
for (i in c(1:nrow(all))) {

    # IIBDGC_GWAS_index_variant in LD with local GWAS lead:
    tmp1<-coloc_results[which(coloc_results$IIBDGC_GWAS_index_variant.x==all$MarkerName[i] & coloc_results$R2_IIBDGC_GWAS_index_variant_MarkerName_eqtl_lead>=0.6),c("IIBDGC_GWAS_index_variant.x","pheno_coloc","qtl_lead","gene","gene_name","Beta_eqtl","SE_eqtl","Pvalue_eqtl","PP.H4.abf","cohort","study_label","condition_name","R2_IIBDGC_GWAS_index_variant_MarkerName_eqtl_lead")]

    # IIBDGC_GWAS_index_variant in LD with eQTL lead:
    tmp2<-coloc_results[which(coloc_results$IIBDGC_GWAS_index_variant.y==all$MarkerName[i] & coloc_results$R2_IIBDGC_GWAS_index_variant_MarkerName_gwas_local_lead>=0.6),c("IIBDGC_GWAS_index_variant.y","pheno_coloc","qtl_lead","PP.H4.abf","Beta_eqtl","SE_eqtl","Pvalue_eqtl","cohort","gene","gene_name","study_label","condition_name","R2_IIBDGC_GWAS_index_variant_MarkerName_gwas_local_lead")]

    colnames(tmp1)[1]<-"IIBDGC_GWAS_index_variant"
    colnames(tmp2)[1]<-"IIBDGC_GWAS_index_variant"

    colnames(tmp1)[13]<-"R2"
    colnames(tmp2)[13]<-"R2"

    tmp<-rbind(tmp1,tmp2)
    tmp<-tmp[!duplicated(tmp),]

    if (i==1) {
        dat<-tmp
    } else {
        dat<-rbind(dat,tmp)
    }

}

dat<-as.data.frame(dat)
dat$Y<-paste(dat$condition_name,dat$pheno_coloc,sep="_")
dat$MarkerName_eqtl_lead<-gsub("_",":",dat$qtl_lead)

dat$chr<-gsub("chr","",dat$MarkerName_eqtl_lead)
dat$chr<-as.numeric(gsub(":.*","",dat$chr))

dat$position<-gsub(":[A-Z]*:[A-Z]*$","",dat$MarkerName_eqtl_lead)
dat$position<-as.numeric(gsub("chr[0-9]{1,2}:","",dat$position))
dim(dat)

# exclude chr12:110259525:G:A - not present in genotyping data
dat<-dat[which(dat$IIBDGC_GWAS_index_variant!="chr12:110259525:G:A"),]
dat$n_variants_in_ld_for_wr<-NA
dat$mean_wald_ratio<-NA
dat$sd_wald_ratio<-NA
dat$betas_wald_ratio<-NA

for (i in 1:nrow(dat))  {

    file_name<-paste0("/path/to/project",dat$cohort[i],"/nominal/",dat$condition_name[i],"/",dat$condition_name[i],".tsv.gz")

    if (dat$cohort[i]=="IBDverse") {
        file_name<-paste0("/path/to/project",dat$condition_name[i],"/",dat$condition_name[i],".tsv.gz")
    } else if (dat$cohort[i]=="hu_2021") {
        file_name<-paste0("/path/to/project",dat$cohort[i],"/nominal/",dat$condition_name[i],"/",dat$condition_name[i],".tsv.gz")
    }

    if (!file.exists(file_name)) {
        print(i)
    }

    param<-GRanges(c(dat$chr[i]), IRanges((dat$position[i]-2000000):(dat$position[i]+2000000)))
    tbx<-Rsamtools::TabixFile(file_name)

    res <- Rsamtools::scanTabix(tbx, param=param)

    tmp <- Map(function(elt) {
        read.csv(textConnection(elt), sep="\t", header=F)
        }, res)
    tmp<-as.data.frame(tmp[1])
    colnames(tmp)<-c("variant","r2","pvalue","molecular_trait_object_id","molecular_trait_id","maf",
                 "gene_id","median_tpm","beta","se","an","ac","chromosome","position","ref","alt","type","rsid")

    tmp$gene<-gsub("\\.[0-9]{1,2}$","",tmp$molecular_trait_object_id)
    tmp<-tmp[which(tmp$gene %in% dat$gene[i]),]
    
    tmp$eqtl_MarkerName<-gsub("_",":",tmp$variant)
    tmp$eqtl_MarkerName<-gsub("chr","",tmp$eqtl_MarkerName)
    tmp$eqtl_MarkerName<-paste0("chr",tmp$eqtl_MarkerName)

    # load ld:
    ld<-fread(paste0(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr",chr,"_",dat$IIBDGC_GWAS_index_variant[i],"_list_653_index_ld_r2_0.6.tags"),head=F)
    colnames(ld)<-"MarkerName"
    
    tmp<-tmp[which(tmp$eqtl_MarkerName %in% ld$MarkerName),]
    # print(dim(tmp))

    # load GWAS results:
    gwas<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",ph,"/",chr,"_",ph,"_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz",sep=""),head=T)
    gwas<-gwas[which(gwas$MarkerName %in% ld$MarkerName),c("MarkerName","BETA","SE")]
    colnames(gwas)[2:3]<-paste(colnames(gwas)[2:3],"gwas",sep="_")

    tmp<-merge(tmp,gwas,by.x="eqtl_MarkerName",by.y="MarkerName")

    tmp$wald_ratio<-tmp$beta/tmp$BETA_gwas

    dat$n_variants_in_ld_for_wr[i]<-nrow(tmp)
    dat$mean_wald_ratio[i]<-mean(tmp$wald_ratio)
    dat$sd_wald_ratio[i]<-sd(tmp$wald_ratio)
    dat$betas_wald_ratio[i]<-paste(tmp$wald_ratio,collapse="|")

}

fwrite(dat,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/tmp_coloc_filtered_by_gwas_r2_0.6_allpheno_with_beta_with_wald_ratio_",chr,"_",ph,".tsv"),
col.names=T,row.names=F,quote=F,sep="\t")

q("no")