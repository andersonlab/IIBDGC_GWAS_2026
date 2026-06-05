# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

# # singularity exec iibdgc_postprocess_10_singularity.sif
# MEM=12000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q cpu-interactive R

library(data.table)
library(rtracklayer)
library(stringr)
library(rtracklayer)
library(qvalue)

rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"
path_eqtl<-"/path/to/project"


# LOAD THE COLOCALIZATION RESULTS - SEE ~/git/IIBDGC_GWAS/scripts/other/coloc_subset_results_by_ld_leads.R
coloc_results<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_filtered_by_gwas_r2_0.6_allpheno.tsv",sep=""),head=T)
coloc_results<-as.data.frame(coloc_results)

# unique set of gene and eqtl variants:
index_col<-coloc_results[,c("condition_name","MarkerName_eqtl_lead","cohort")]
index_col<-index_col[!duplicated(index_col),]
dim(index_col)
# [1] 12986    3

index_col$chr<-gsub("chr","",index_col$MarkerName_eqtl_lead)
index_col$chr<-as.numeric(gsub(":.*","",index_col$chr))

index_col$position<-gsub(":[A-Z]*:[A-Z]*$","",index_col$MarkerName_eqtl_lead)
index_col$position<-as.numeric(gsub("chr[0-9]{1,2}:","",index_col$position))
dim(index_col)
# [1] 12986    5



for (i in 1:nrow(index_col))  {

    file_name<-paste0("/path/to/project",index_col$cohort[i],"/nominal/",index_col$condition_name[i],"/",index_col$condition_name[i],".tsv.gz")

    if (index_col$cohort[i]=="IBDverse") {
        file_name<-paste0("/path/to/project",index_col$condition_name[i],"/",index_col$condition_name[i],".tsv.gz")
    } else if (index_col$cohort[i]=="hu_2021") {
        file_name<-paste0("/path/to/project",index_col$cohort[i],"/nominal/",index_col$condition_name[i],"/",index_col$condition_name[i],".tsv.gz")
    }

    if (!file.exists(file_name)) {
        print(i)
    }

    param<-GRanges(c(index_col$chr[i]), IRanges((index_col$position[i]-100000):(index_col$position[i]+100000)))
    tbx<-Rsamtools::TabixFile(file_name)

    res <- Rsamtools::scanTabix(tbx, param=param)

    tmp <- Map(function(elt) {
        read.csv(textConnection(elt), sep="\t", header=F)
        }, res)
    tmp<-as.data.frame(tmp[1])
    colnames(tmp)<-c("variant","r2","pvalue","molecular_trait_object_id","molecular_trait_id","maf",
                 "gene_id","median_tpm","beta","se","an","ac","chromosome","position","ref","alt","type","rsid")

    ld_var<-fread(paste0(path_gwas,"post_imputation/2022/analysis/ld_data/eur_tier2/",index_col$cohort[i],"/",index_col$cohort[i],"_",index_col$condition_name[i],"_allchr_list_tier_2_ld_r2_0.25.bed"),head=T)


    tmp$X<-paste(index_col$condition_name[i],tmp$gene_id,tmp$variant,sep="_")

    if (index_col$cohort[i]=="pQTL_sparc") {
            tmp$X<-paste0(index_col$condition_name[i],"_",tmp$gene_id,"_chr",tmp$variant)
    }

    tmp1.1<-tmp[which(tmp$X %in% coloc_results$X),]

    tmp$X<-paste(index_col$condition_name[i],tmp$molecular_trait_object_id,tmp$variant,sep="_")
    if (index_col$cohort[i]=="pQTL_sparc") {
            tmp$X<-paste0(index_col$condition_name[i],"_",tmp$molecular_trait_object_id,"_chr",tmp$variant)
    }

    tmp1.2<-tmp[which(tmp$X %in% coloc_results$X),]

    tmp<-rbind(tmp1.1,tmp1.2)
    tmp<-tmp[!duplicated(tmp),]

    if(i==1) {
        dat<-tmp
    } else {
        dat<-rbind(dat,tmp)
    }

    rm(tmp,tmp1.1,tmp1.2)

    if (i %in% c(100,200,300)) {
        print(i)
    }

}


dim(coloc_results[!coloc_results$X %in% dat$X,])
# [1] 0   38

dat$Beta_eqtl<-dat$beta
dat$SE_eqtl<-dat$se
dat$Pvalue_eqtl<-dat$pvalue


dim(coloc_results)
# [1] 20659    38
coloc_results<-merge(coloc_results,dat[,c("X","Beta_eqtl","SE_eqtl","Pvalue_eqtl")],by="X",all.x=T)
dim(coloc_results)
# [1] 20659    41

write.table(coloc_results,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_filtered_by_gwas_r2_0.6_allpheno_with_beta.tsv",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")

q("no")



##############################################################################################################################
##############################################################################################################################
