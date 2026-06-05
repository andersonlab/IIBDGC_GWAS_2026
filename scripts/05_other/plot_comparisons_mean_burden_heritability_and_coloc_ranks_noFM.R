# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# PLOT COMPARISONS MEAN BURDEN HERITABILITY AND COLOC:


# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=16000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -q cpu-interactive -G your_hpc_group R

library(data.table)
library(stringr)
library(plyr)
library(dplyr)
library(ggalluvial)
library(ggpubr)
library(MASS)
library(ggplot2)


#################################################################################################################################################
### CREATE PANEL 1: (see original test in ~/git/IIBDGC_GWAS/scripts/other/evaluate_mean_burden_heritability_vs_ldsc_tissue_zcore.R)

rm(list=ls())

pheno<-c("ibd","cd","uc")


path_gwas<-"/path/to/ibdgwas/IIBDGC/"
path_eqtl<-"/path/to/project"

# LOAD THE INDEPENDENT SIGNALS - see ~/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_38_2_define_independent_signals_using_as_ref_cojo_agnostic_results.R
all<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants.tsv.gz",sep=""))
all<-as.data.frame(all)
dim(all)
# [1] 631 182

cols<-rev(c("#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d"))


for (ph in (pheno)) {

    dat_final<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_matrices/matrix_",ph,"_zscore_genes_tissue_noFM.tsv.gz"),head=T)
    dat_final<-as.data.frame(dat_final)
    dat_final<-dat_final[,c("IIBDGC_GWAS_index_variant","gene_name","score")]
    dim(dat_final)
    # [1] 708   3

    dat_final<-dat_final[which(dat_final$IIBDGC_GWAS_index_variant %in% all$MarkerName),]

    dat_final<-merge(dat_final,all[,c("MarkerName","phenotype")],by.x="IIBDGC_GWAS_index_variant",by.y="MarkerName",all.x=T)

    # if (ph=="cd") {
    #     dat_final<-dat_final[which(dat_final$phenotype %in% c("IBD_unsaturated","IBD_saturated","CD")),]
    # } else if (ph=="uc") {
    #     dat_final<-dat_final[which(dat_final$phenotype %in% c("IBD_unsaturated","IBD_saturated","UC")),]
    # } 

    # for each signal linked to more than one gene, does the rank from ldsc matches the burden heritability rank?

    dat_final<-dat_final[duplicated(dat_final$IIBDGC_GWAS_index_variant),]
    
    gene_reference <- readRDS(paste0("/path/to/project",ph,"_nsyn_am.rds"))
    class(gene_reference)
    gene_reference<-gene_reference[which(gene_reference$geneid %in% dat_final$gene_name),]
    dim(gene_reference)
    # [1] 253  53

    dat_final<-dat_final[which(dat_final$gene_name %in% gene_reference$geneid),]
    dim(dat_final)
    # [1] 261   4

    dat_final<-dat_final[duplicated(dat_final$IIBDGC_GWAS_index_variant),]
    dim(dat_final)
    # [1] 127  3

    gene_reference<-gene_reference[which(gene_reference$geneid %in% dat_final$gene_name),]
    dim(gene_reference)
    # [1] 125  53


    dat_final<-merge(dat_final,gene_reference[,c("geneid","r2")],by.x="gene_name",by.y="geneid")

    dups<-dat_final$IIBDGC_GWAS_index_variant[duplicated(dat_final$IIBDGC_GWAS_index_variant)]

    dups<-dups[!duplicated(dups)]
    length(dups)
    # [1] 34

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

    xtab<-table(datrank$rank_tissue_zcore,datrank$rank_r2)
    print(ph)
    print(xtab)
    pval<-fisher.test(xtab,workspace = 6e8,simulate.p.value=TRUE,B=500000)$p.value
    print(pval)

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

    datrank_alluvial$rank<-as.factor(datrank_alluvial$rank)
    datrank_alluvial$rank<-factor(datrank_alluvial$rank,levels=paste0(rep("rank",10),seq(1:10)))

    p<-ggplot(datrank_alluvial, aes(alluvium = ID, x = class, stratum = rank)) + 
    geom_alluvium(color = "black") +
    geom_stratum( color = "black",aes(fill=rank))  + 
    # Vanilla GGplot here onwards
    ggtitle(paste("Gene ranks",toupper(ph),"\nFisher exact test pval:",formatC(pval,format="e",digits=2)))+
    scale_y_discrete() +
    ylab("Gene variant pair") +
    theme_bw() + 
    theme(axis.title=element_blank(),legend.position="bottom",legend.title = element_blank(),panel.border = element_blank(),plot.margin = unit(c(1,1,1,1), "cm")) +
    scale_fill_manual(values = cols) 

    assign(ph,p)

}

# [1] "ibd"
   
#      1  2  3  4  5  6  7  8  9
#   1 13  9  4  1  0  2  0  0  0
#   2 11 13  4  0  0  0  0  1  0
#   3  1  5  2  3  0  0  0  0  0
#   4  3  0  0  0  0  0  0  0  1
#   5  0  1  0  0  1  0  0  0  0
#   6  0  1  0  0  1  0  0  0  0
#   7  0  0  1  0  0  0  0  0  0
#   8  0  0  0  0  0  0  1  0  0
#   9  1  0  0  0  0  0  0  0  0
# [1] 0.002487995
# [1] "cd"
   
#      1  2  3  4  5  6  7  8  9
#   1 15 10  2  0  1  0  0  0  1
#   2  8 16  4  1  0  0  0  0  0
#   3  4  2  3  1  0  0  1  0  0
#   4  2  0  0  2  0  0  0  0  0
#   5  0  0  0  0  0  2  0  0  0
#   6  0  0  2  0  0  0  0  0  0
#   7  0  1  0  0  0  0  0  0  0
#   8  0  0  0  0  1  0  0  0  0
#   9  0  0  0  0  0  0  0  1  0
# [1] 1.799996e-05
# [1] "uc"
   
#      1  2  3  4  5  6  7  8  9
#   1 14 11  2  0  0  1  0  1  0
#   2 10 12  5  2  0  0  0  0  0
#   3  2  3  4  1  0  1  0  0  0
#   4  1  1  0  0  2  0  0  0  0
#   5  0  2  0  0  0  0  0  0  0
#   6  1  0  0  0  0  0  1  0  0
#   7  0  0  0  1  0  0  0  0  0
#   8  0  0  0  0  0  0  0  0  1
#   9  1  0  0  0  0  0  0  0  0
# [1] 0.003795992

p1<-ggarrange(cd,ibd,uc,ncol=3,common.legend=T,legend="right")




###########################################################################################################################################################################################################
### panel 2 - data generated in ~/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_53_combine_all_effector_gene_evidence_vs2.R

pheno <- c("ibd", "cd", "uc")
tests <- c(
    "genes_coloc", "genes_coloc_rank1", "genes_closest", "genes_closest_coloc_signals",
    "genes_closest_no_coloc_signals", "genes_otar_gentropy_0.8", "genes_otar_gentropy_0.5",
    "genes_coloc_score_0.8", "genes_coloc_score_0.7", "genes_coloc_score_0.6", "genes_coloc_score_0.5",
    "genes_coloc_mr_exo", "genes_coloc_rank1_mr_exo", "genes_closest_no_coding",
    "genes_coloc_rank1_no_coding", "genes_coloc_no_coding",
    "genes_closest_no_coloc_signals_no_coding", "genes_closest_coloc_signals_no_coding"
)

group_levels <- c(
    "genes_coloc", "genes_coloc_no_coding",
    "genes_coloc_rank1", "genes_coloc_rank1_no_coding",
    "genes_closest", "genes_closest_no_coding",
    "genes_closest_coloc_signals", "genes_closest_no_coloc_signals_no_coding",
    "genes_closest_no_coloc_signals", "genes_closest_coloc_signals_no_coding",
    "genes_otar_gentropy_0.8", "genes_otar_gentropy_0.5",
    "genes_coloc_score_0.8", "genes_coloc_score_0.7", "genes_coloc_score_0.6", "genes_coloc_score_0.5",
    "genes_coloc_mr_exo", "genes_coloc_rank1_mr_exo"
)

# without NOD2 — load all combinations and bind at once
dat <- rbindlist(lapply(pheno, function(ph) {
    rbindlist(lapply(tests, function(test) {
        tmp <- fread(paste0(
            path_gwas, "post_imputation/2022/analysis/metaanalysis/coloc_results/test/",
            test, "_", ph, "_nod2_exclussion_list_output"
        ), head = T)
        tmp$pheno <- ph
        tmp$group <- test
        tmp$nod2  <- "excluded"
        tmp
    }))
}))


dat <- as.data.frame(dat)
dat$group <- factor(dat$group, levels = group_levels)

mu <- dat[dat$permutate == "Target_genes", ]

table(mu$group, mu$nod2)

# Vectorised p-value computation via data.table join
dat_dt  <- as.data.table(dat)
mu_keys <- as.data.table(mu)[, .(pheno, nod2, group, target_r2 = mean_r2)]

pvals <- dat_dt[permutate != "Target_genes"][
    mu_keys, on = c("pheno", "nod2", "group"),
    .(pvalue = sum(mean_r2 >= target_r2) / 1e6),
    by = .EACHI
][, .(pheno, nod2, group, pvalue)]

mu <- merge(mu, pvals, by = c("pheno", "nod2", "group"), all.x = TRUE)

mu$logpvalue <- -log10(mu$pvalue)
mu$logpvalue[is.infinite(mu$logpvalue)] <- -log10(1 / 1e6)

plot_groups <- c("genes_coloc_rank1_mr_exo", "genes_closest", "genes_coloc", "genes_coloc_rank1")

xmin <- min(dat$mean_r2[dat$group %in% plot_groups])
xmax <- max(dat$mean_r2[dat$group %in% plot_groups])

panel_plots <- lapply(plot_groups, function(grp) {
    tmp    <- dat[dat$group == grp & dat$nod2 == "excluded" & dat$permutate != "Target_genes", ]
    mu_tmp <- mu[mu$group == grp, ]
    mu_tmp$pvalue_label <- ifelse(mu_tmp$pvalue == 0, "<1E-6", as.character(mu_tmp$pvalue))

    ggplot(tmp, aes(x = mean_r2)) +
        geom_histogram(color = "black", fill = "white") +
        geom_vline(data = mu_tmp, aes(xintercept = mean_r2),
                   linetype = "dotted", color = "blue", linewidth = 1.5) +
        geom_text(data = mu_tmp,
                  mapping = aes(x = 1.8E-6, y = 150000,
                                label = paste("Pvalue\n", pvalue_label))) +
        facet_grid(~ pheno) +
        xlim(xmin, xmax) +
        ggtitle(grp) +
        theme_bw()
})
names(panel_plots) <- plot_groups

p2 <- panel_plots[["genes_coloc_rank1_mr_exo"]]
p3 <- panel_plots[["genes_closest"]]
p4 <- panel_plots[["genes_coloc"]]
p5 <- panel_plots[["genes_coloc_rank1"]]



p<-ggarrange(p1,p2,p3,p4,p5,nrow=5,heights=c(6,2,2,2,2),align=c("v"),labels=c("a","b","c","d","e"))

ggsave(
  "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_burden_heritability_tissue_zscore_allgenes_rank.png",
  p,
  width = 180,
  height = 250,
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)

q("no")