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

path_gwas<-"/path/to/ibdgwas/IIBDGC/"
pheno<-c("ibd","cd","uc")

# LOAD THE INDEPENDENT SIGNALS - see ~/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_38_2_define_independent_signals_using_as_ref_cojo_agnostic_results.R
all<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_with_fm_signals_with_multiancestry_signals.tsv.gz",sep=""))
all<-as.data.frame(all)
dim(all)
# [1] 579 211

# six of them to exclude: no significant at all in EUR tier1 only signif when tier 2 is included:
ids_exclude<-c("chr2:18792705:C:CT","chr9:19972835:AG:A","chr4:42396369:G:T","chr4:42635348:ATATC:A","chr15:34731239:A:C","chr10:125821742:T:G")  
all<-all[which(!all$MarkerName %in% ids_exclude),]
dim(all)
# [1] 573 211


for (ph in pheno) {

    tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_matrices/matrix_",ph,"_zscore_genes_tissue.tsv.gz"),head=T)

    tmp$pheno<-ph

    if(ph=="ibd") {
        dat_final<-tmp
    } else {
        dat_final<-rbind(dat_final,tmp)
    }

}

dim(dat_final)
# [1] 1358  499

table(dat_final$pheno)
#  cd ibd  uc 
# 449 498 411 

dat_final<-as.data.frame(dat_final)

for (ph in pheno) {

    print(ph)

    if (ph=="ibd") {

        tmp1<-dat_final[which(dat_final$pheno %in% c("ibd") & dat_final$IIBDGC_GWAS_index_variant %in% all$MarkerName[which(all$phenotype %in% c("IBD_unsaturated","IBD_saturated","CD","UC"))]),c("IIBDGC_GWAS_index_variant","gene_name","score")]
        tmp1<-tmp1[order(tmp1$gene_name,decreasing=T),]
        tmp1<-tmp1[order(tmp1$score,decreasing=T),]
        tmp1<-tmp1[!duplicated(tmp1$gene_name),]

        # tmp2<-dat_final[which(dat_final$pheno %in% c("cd","uc") & dat_final$IIBDGC_GWAS_index_variant %in% all$MarkerName[which(all$phenotype %in% c("IBD_unsaturated","IBD_saturated","CD","UC"))]),c("IIBDGC_GWAS_index_variant","gene_name","score")]
        # tmp2<-tmp2[which(!tmp2$gene_name %in% tmp1$gene_name),]
        # tmp2<-tmp2[order(tmp2$score,decreasing=T),]
        # tmp2<-tmp2[!duplicated(tmp2$gene_name),]

        # tmp<-rbind(tmp1,tmp2)
        tmp<-tmp1

    } else if (ph=="cd") {

        tmp1<-dat_final[which(dat_final$pheno %in% c("cd") & dat_final$IIBDGC_GWAS_index_variant %in% all$MarkerName[which(all$phenotype %in% c("IBD_unsaturated","IBD_saturated","CD"))]),c("IIBDGC_GWAS_index_variant","gene_name","score")]
        tmp1<-tmp1[order(tmp1$gene_name,decreasing=T),]
        tmp1<-tmp1[order(tmp1$score,decreasing=T),]
        tmp1<-tmp1[!duplicated(tmp1$gene_name),]

        # tmp2<-dat_final[which(dat_final$pheno %in% c("ibd") & dat_final$IIBDGC_GWAS_index_variant %in% all$MarkerName[which(all$phenotype %in% c("IBD_unsaturated","IBD_saturated","CD"))]),c("IIBDGC_GWAS_index_variant","gene_name","score")]
        # tmp2<-tmp2[which(!tmp2$gene_name %in% tmp1$gene_name),]
        # tmp2<-tmp2[order(tmp2$score,decreasing=T),]
        # tmp2<-tmp2[!duplicated(tmp2$gene_name),]

        # tmp<-rbind(tmp1,tmp2)
        tmp<-tmp1

    } else if (ph=="uc") {

        tmp1<-dat_final[which(dat_final$pheno %in% c("uc") & dat_final$IIBDGC_GWAS_index_variant %in% all$MarkerName[which(all$phenotype %in% c("IBD_unsaturated","IBD_saturated","UC"))]),c("IIBDGC_GWAS_index_variant","gene_name","score")]
        tmp1<-tmp1[order(tmp1$gene_name,decreasing=T),]
        tmp1<-tmp1[order(tmp1$score,decreasing=T),]
        tmp1<-tmp1[!duplicated(tmp1$gene_name),]

        # tmp2<-dat_final[which(dat_final$pheno %in% c("ibd") & dat_final$IIBDGC_GWAS_index_variant %in% all$MarkerName[which(all$phenotype %in% c("IBD_unsaturated","IBD_saturated","UC"))]),c("IIBDGC_GWAS_index_variant","gene_name","score")]
        # tmp2<-tmp2[which(!tmp2$gene_name %in% tmp1$gene_name),]
        # tmp2<-tmp2[order(tmp2$score,decreasing=T),]
        # tmp2<-tmp2[!duplicated(tmp2$gene_name),]

        # tmp<-rbind(tmp1,tmp2)
        tmp<-tmp1
    }

    rm(tmp1,tmp2)

    tmp<-tmp[,c("gene_name","score"),drop=F]
    tmp<-tmp[order(tmp$score,decreasing=T),]
    tmp<-tmp[!duplicated(tmp$gene_name),,drop=F]

    print(dim(tmp))

    # tmp<-tmp[!grepl("ENSG",tmp$gene_name),]
    tmp<-tmp[which(!is.na(tmp$score)),]

    tmp<-tmp %>% mutate(group = ntile(score, 2))


    for (i in 1:2) {

        print(nrow(tmp[which(tmp$group==i),]))
        write.table(tmp[which(tmp$group==i),"gene_name",drop=F],paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/top_",i,"_",ph,"_zscore_genes_tissue"),col.names=F,row.names=F,sep="\t")

    }

    write.table(tmp[,"gene_name",drop=F],paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/allgenes_",ph,"_zscore_genes_tissue"),col.names=F,row.names=F,sep="\t")
    print(nrow(tmp))

}
    
# ONLY COLOC BY PHENO
# [1] "ibd"
# [1] 463   2
# [1] 232
# [1] 231
# [1] 463
# [1] "cd"
# [1] 415   2
# [1] 208
# [1] 207
# [1] 415
# [1] "uc"
# [1] 371   2
# [1] 186
# [1] 185
# [1] 371


# [1] "ibd"
# [1] 667   2
# [1] 334
# [1] 333
# [1] 667
# [1] "cd"
# [1] 537   2
# [1] 269
# [1] 268
# [1] 537
# [1] "uc"
# [1] 479   2
# [1] 240
# [1] 239
# [1] 479


q("no")