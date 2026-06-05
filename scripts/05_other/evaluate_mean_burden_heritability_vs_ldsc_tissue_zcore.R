# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# # singularity exec iibdgc_postprocess_10_singularity.sif

# MEM=8000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q cpu-interactive R

library(data.table)
library(stringr)
library(plyr)
library(dplyr)
library(ggalluvial)
library(ggpubr)
library(MASS)

rm(list=ls())

# provide PHENOTYPE
args = commandArgs(trailingOnly=TRUE)

pheno<-c("ibd","cd","uc")


path_gwas<-"/path/to/ibdgwas/IIBDGC/"
path_eqtl<-"/path/to/project"

# LOAD THE INDEPENDENT SIGNALS - see ~/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_38_2_define_independent_signals_using_as_ref_cojo_agnostic_results.R
all<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_with_fm_signals_with_multiancestry_signals_with_pheno.tsv.gz",sep=""))
all<-as.data.frame(all)
dim(all)

cols<-rev(c("#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d"))


for (ph in (pheno)) {

    print(ph)

    dat_final<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_matrices/matrix_",ph,"_zscore_genes_tissue.tsv.gz"),head=T)
    dat_final<-dat_final[,c("IIBDGC_GWAS_index_variant","gene_name","score")]
    dim(dat_final)
    # [1] 498   3

    dat_final<-merge(dat_final,all[,c("MarkerName","phenotype")],by.x="IIBDGC_GWAS_index_variant",by.y="MarkerName",all.x=T)

    if (ph=="cd") {
        dat_final<-dat_final[which(dat_final$phenotype %in% c("IBD_unsaturated","IBD_saturated","CD")),]
    } else if (ph=="uc") {
        dat_final<-dat_final[which(dat_final$phenotype %in% c("IBD_unsaturated","IBD_saturated","UC")),]
    }

    # for each signal linked to more than one gene, does the rank from ldsc matches the burden heritability rank?

    dat_final<-dat_final[duplicated(dat_final$IIBDGC_GWAS_index_variant),]
    gene_reference <- readRDS(paste0("/path/to/project",ph,"_nsyn_am.rds"))
    gene_reference<-gene_reference[which(gene_reference$geneid %in% dat_final$gene_name),]
    dim(gene_reference)
    # [1] 174  51


    dat_final<-dat_final[which(dat_final$gene_name %in% gene_reference$geneid),]
    dim(dat_final)
    # [1] 174  51

    dat_final<-dat_final[duplicated(dat_final$IIBDGC_GWAS_index_variant),]
    dim(dat_final)
    # [1] 90  3

    gene_reference<-gene_reference[which(gene_reference$geneid %in% dat_final$gene_name),]
    dim(gene_reference)
    # [1] 87 51


    dat_final<-merge(dat_final,gene_reference[,c("geneid","r2")],by.x="gene_name",by.y="geneid")

    dups<-dat_final$IIBDGC_GWAS_index_variant[duplicated(dat_final$IIBDGC_GWAS_index_variant)]
    length(dups)
    # [1] 45

    dat_final$rank_tissue_zcore<-NA
    dat_final$rank_r2<-NA

    for (i in 1:length(dups)) {

        tmp<-dat_final[which(dat_final$IIBDGC_GWAS_index_variant %in% dups[i]),]

        tmp<-tmp[order(tmp$score,decreasing=T),]
        tmp$rank_tissue_zcore<-seq(1:nrow(tmp))

        tmp<-tmp[order(tmp$r2,decreasing=T),]
        tmp$rank_r2<-seq(1:nrow(tmp))

        if(i==1) {
            datrank<-tmp
        } else {
            datrank<-rbind(datrank,tmp)
        }

    }

    datrank$ID<-paste(datrank$IIBDGC_GWAS_index_variant,datrank$gene_name,sep="_")
    datrank<-datrank[!duplicated(datrank$ID),]

    print(nrow(datrank))
    print(nrow(datrank[which(datrank$rank_r2==datrank$rank_tissue_zcore)]))
    print(nrow(datrank[which(datrank$rank_r2==1)]))
    print(nrow(datrank[which(datrank$rank_r2==datrank$rank_tissue_zcore & datrank$rank_r2==1)]))

    xtab<-table(datrank$rank_tissue_zcore,datrank$rank_r2)
    print(ph)
    print(xtab)
    pval<-fisher.test(xtab,workspace = 6e8,simulate.p.value=TRUE,B=500000)$p.value
    print(pval)
    # model1<-polr(formula = as.factor(rank_tissue_zcore) ~ r2, data = datrank, Hess = TRUE)
    # table<-coef(summary(model1))[1,]
    # p<-pnorm(abs(table[3]),lower.tail = F,log.p = TRUE)*2
    # round(p, 4) too low pval

    # not of use with > 4 categorical 
    # tidy(model1, exponentiate = TRUE,p.values = TRUE)

    # model 
    # library(ordinal)
    # mod<-clmm(as.factor(rank_tissue_zcore) ~ r2 + 1|IIBDGC_GWAS_index_variant, data=datrank, Hess=T)

    datrank<-datrank[order(datrank$IIBDGC_GWAS_index_variant),]
    datrank$ID<-paste(datrank$IIBDGC_GWAS_index_variant,datrank$gene_name,sep="_")
    datrank_1<-datrank[,c("ID","rank_tissue_zcore")]
    datrank_1$class<-"tissue_zcore"
    datrank_2<-datrank[,c("ID","rank_r2")]
    datrank_2$class<-"burden_r2"

    colnames(datrank_1)[2]<-"rank"
    colnames(datrank_2)[2]<-"rank"

    datrank_alluvial<-rbind(datrank_1,datrank_2)

    datrank_alluvial$rank<-paste0("rank",datrank_alluvial$rank)

    is_lodes_form(datrank_alluvial, key = "class", value = "rank", id = "ID")

    datrank_alluvial$class<-factor(datrank_alluvial$class,levels=c("tissue_zcore","burden_r2"))

    p<-ggplot(datrank_alluvial, aes(alluvium = ID, x = class, stratum = rank)) + 
    geom_alluvium(color = "black") +
    geom_stratum( color = "black",aes(fill=rank))  + 
    # Vanilla GGplot here onwards
    ggtitle(paste("Gene ranks",ph,"\nFisher exact test pval:",formatC(pval,format="e",digits=2)))+
    scale_y_discrete() +
    ylab("Gene variant pair")+
    theme_bw()+
    scale_fill_manual(values = cols) +
    theme(axis.text = element_text(size = 7))

    assign(ph,p)

}

p<-ggarrange(ibd,cd,uc,ncol=3)


ggsave(
  "~/git/IIBDGC_GWAS/plots/paper_figures/Figure_burden_heritability_tissue_zscore_allgenes_rank.pdf",
  p,
  width = 120,
  height = 90,
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)

q("no")