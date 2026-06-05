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

rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

#################################################################################################################################################
### CREATE PANEL 1: (see original test in ~/git/IIBDGC_GWAS/scripts/other/evaluate_mean_burden_heritability_vs_ldsc_tissue_zcore.R)


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

    xtab<-table(datrank$rank_tissue_zcore,datrank$rank_r2)
    pval<-fisher.test(xtab,workspace = 6e8,simulate.p.value=TRUE,B=500000)$p.value

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
    datrank_alluvial$rank<-factor(datrank_alluvial$rank,levels=paste0(rep("rank",9),seq(1:9)))

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

p1<-ggarrange(cd,ibd,uc,ncol=3,common.legend=T,legend="right")




###########################################################################################################################################################################################################
### panel 2 - data generated in ~/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_53_combine_all_effector_gene_evidence_vs2.R

pheno<-c("ibd","cd","uc")
tests<-c("genes_coloc","genes_coloc_rank1","genes_closest","genes_closest_coloc_signals",
"genes_closest_no_coloc_signals","genes_otar_gentropy_0.8","genes_otar_gentropy_0.5",
"genes_coloc_score_0.8","genes_coloc_score_0.7","genes_coloc_score_0.6","genes_coloc_score_0.5",
"genes_coloc_mr_exo","genes_coloc_rank1_mr_exo","genes_closest_no_coding","genes_coloc_rank1_no_coding","genes_coloc_no_coding",
"genes_closest_no_coloc_signals_no_coding","genes_closest_coloc_signals_no_coding")

# without NOD2
for (ph in pheno) {

    for (test in tests) {

        tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/",test,"_",ph,"_nod2_exclussion_list_output"),head=T)
        tmp$pheno<-ph
        tmp$group<-test
        tmp$nod2<-"excluded"

        if (ph=="ibd" & test=="genes_coloc") {
            dat<-tmp
        } else {
            dat<-rbind(dat,tmp)
        }
    }

}


dat<-as.data.frame(dat)
dat$group<-as.factor(dat$group)

dat$group<-factor(dat$group,levels=c("genes_coloc","genes_coloc_no_coding",
"genes_coloc_rank1","genes_coloc_rank1_no_coding",
"genes_closest","genes_closest_no_coding",
"genes_closest_coloc_signals","genes_closest_no_coloc_signals_no_coding",
"genes_closest_no_coloc_signals","genes_closest_coloc_signals_no_coding",
"genes_otar_gentropy_0.8","genes_otar_gentropy_0.5","genes_coloc_score_0.8","genes_coloc_score_0.7","genes_coloc_score_0.6","genes_coloc_score_0.5",
"genes_coloc_mr_exo","genes_coloc_rank1_mr_exo"))



mu <- dat[which(dat$permutate=="Target_genes"),]
mu$pvalue<-NA


table(mu$group,mu$nod2)

for (ph in pheno) {
    for (i in 1:length(levels(mu$group))) {
        for (inc in c("excluded","included")) {
            mu$pvalue[which(mu$pheno==ph & mu$nod2==inc & mu$group==levels(mu$group)[i])]<-nrow(dat[which(dat$pheno==ph & dat$nod2==inc & dat$group==levels(mu$group)[i] & dat$permutate!="Target_genes" & dat$mean_r2>=mu$mean_r2[which(mu$pheno==ph & mu$nod2==inc & mu$group==levels(mu$group)[i])]),])/1000000
        }
    }
}


mu$logpvalue<-(-log10(mu$pvalue))
mu$logpvalue[which(mu$logpvalue==Inf)]<-(-log10(1/1000000))


list<-c("genes_coloc_rank1_mr_exo","genes_closest")

xmin<-min(dat$mean_r2[which(dat$group %in% list)])
xmax<-max(dat$mean_r2[which(dat$group %in% list)])

for (i in 1:length(list)) {

    tmp<-dat[which(dat$group==list[i] & dat$nod2=="excluded" & dat$mean_r2!="Target_genes"),]
    mu_tmp<-mu[which(mu$group==list[i]),]

    p<-ggplot(tmp, aes(x=mean_r2)) +
    geom_histogram(color="black", fill="white") + theme_bw() + geom_vline(data=mu_tmp,aes(xintercept = mean_r2), linetype="dotted", 
                color = "blue", size=1.5) + xlim(xmin,xmax) + facet_grid( ~ pheno)

    # save tmp plots for inicial visualization 
    # ggsave(
    #     paste0("~/git/IIBDGC_GWAS/plots/tmp_figures/",list[i],".png"),
    #     p,
    #     width = 120,
    #     height = 120,
    #     dpi = 300,
    #     units = c("mm"),
    #     limitsize = T,scale=2
    #     )

    assign(list[i],p)

    rm(tmp)

}

p2<-genes_coloc_rank1_mr_exo
p3<-genes_closest


# p2<-ggplot(mu[which(mu$group %in% c("genes_closest","genes_coloc_rank1_mr_exo")),], 
# aes(x=-log10(mean_r2),y=group,color=logpvalue,size=logpvalue)) +
#  geom_point() + facet_grid( ~ pheno) + theme(axis.title.y=element_blank(),) + labs(x = "Mean burden heritability (log(r2))") + 
# #  labs(color='-log10(pval)',size='-log10(pval)') + 
# #  lims(colour = c(0,8),size = c(0,8)) + 
# scale_color_continuous(limits=c(0,6), breaks=seq(0,6, by=1),name="-log10(pval)") +
# # guides(color= guide_legend(), size=guide_legend()) +
# # scale_size_continuous(limits=c(0,6), breaks=seq(0,6, by=1),name="-log10(pval)") + 
# theme(plot.margin = unit(c(1,1,1,0.5), "cm")) + 
# scale_y_discrete(labels=c("genes_closest" = "Closest Gene", "genes_coloc_rank1_mr_exo" = "Combined Coloc rank1, Exonic, MR")) 
# + xlim(0,6)




p<-ggarrange(p1,p2,p3,nrow=3,heights=c(6,2,2),align=c("v"),labels=c("a","b","c"))

ggsave(
  "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_burden_heritability_tissue_zscore_allgenes_rank.png",
  p,
  width = 150,
  height = 150,
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)


### plot a supplementary table:



# chekc out OTAR
# ggplot(mu[which(mu$group %in% c("genes_coloc_rank1_mr_exo","genes_coloc_rank1","genes_otar_gentropy_0.8","genes_otar_gentropy_0.5","genes_closest")),], aes(x=mean_r2,y=group,color=logpvalue,size=logpvalue)) +
#  geom_point() + facet_grid( ~ pheno) + theme(axis.title.y=element_blank(),) + labs(x = "Mean burden heritability (r2)") + 
# scale_color_continuous(limits=c(0,6), breaks=seq(0,6, by=1),name="-log10(pval)") +
# guides(color= guide_legend(), size=guide_legend()) +
# scale_size_continuous(limits=c(0,6), breaks=seq(0,6, by=1),name="-log10(pval)") + theme(plot.margin = unit(c(1,1,1,0.5), "cm"))





# p4<-ggplot(mu, aes(x=mean_r2,y=group,color=logpvalue,size=logpvalue)) +
#  geom_point() + facet_grid( ~ pheno) + theme(axis.title.y=element_blank(),) + labs(x = "Mean burden heritability (r2)") + 
# #  labs(color='-log10(pval)',size='-log10(pval)') + 
# #  lims(colour = c(0,8),size = c(0,8)) + 
# scale_color_continuous(limits=c(0,6), breaks=seq(0,6, by=1),name="-log10(pval)") +
# guides(color= guide_legend(), size=guide_legend()) +
# scale_size_continuous(limits=c(0,6), breaks=seq(0,6, by=1),name="-log10(pval)") + theme(plot.margin = unit(c(1,1,1,0.5), "cm"))

# ggsave(
#   "~/git/IIBDGC_GWAS/plots/metaanalysis/all_zscore_allgenes_rank.png",
#   p4,
#   width = 140,
#   height = 80,
#   dpi = 400,
#   units = c("mm"),
#   limitsize = T,scale=2
# )

# p5<-ggplot(mu[which(mu$group %in% c("genes_coloc","genes_coloc_rank1","genes_closest_coloc_signals","genes_closest_no_coloc_signals")),], aes(x=mean_r2,y=group,color=logpvalue,size=logpvalue)) +
#  geom_point() + facet_grid( ~ pheno) + theme(axis.title.y=element_blank(),) + labs(x = "Mean burden heritability (r2)") + 
# #  labs(color='-log10(pval)',size='-log10(pval)') + 
# #  lims(colour = c(0,8),size = c(0,8)) + 
# scale_color_continuous(limits=c(0,6), breaks=seq(0,6, by=1),name="-log10(pval)") +
# guides(color= guide_legend(), size=guide_legend()) +
# scale_size_continuous(limits=c(0,6), breaks=seq(0,6, by=1),name="-log10(pval)") + theme(plot.margin = unit(c(1,1,1,0.5), "cm")) +
# xlim(-13,0)

# p6<-ggplot(mu[which(mu$group %in% c("genes_coloc_no_coding","genes_coloc_rank1_no_coding","genes_closest_coloc_signals_no_coding","genes_closest_no_coloc_signals_no_coding")),], aes(x=mean_r2,y=group,color=logpvalue,size=logpvalue)) +
#  geom_point() + facet_grid( ~ pheno) + theme(axis.title.y=element_blank(),) + labs(x = "Mean burden heritability (r2)") + 
# #  labs(color='-log10(pval)',size='-log10(pval)') + 
# #  lims(colour = c(0,8),size = c(0,8)) + 
# scale_color_continuous(limits=c(0,6), breaks=seq(0,6, by=1),name="-log10(pval)") +
# guides(color= guide_legend(), size=guide_legend()) +
# scale_size_continuous(limits=c(0,6), breaks=seq(0,6, by=1),name="-log10(pval)") + theme(plot.margin = unit(c(1,1,1,0.5), "cm")) +
# xlim(-13,0)



# #  + 
# # scale_y_discrete(labels=c("genes_coloc" = "Colocalization Genes", "genes_coloc_rank1" = "Colocalization Genes, rank1", "genes_closest_no_coding" = "Closest Gene (no coding signals)"))

# p56<-ggarrange(p5,p6,nrow=2)


# ggsave(
#   "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_burden_heritability_tissue_zscore_allgenes_rank.png",
#   p3,
#   width = 140,
#   height = 40,
#   dpi = 400,
#   units = c("mm"),
#   limitsize = T,scale=2
# )


q("no")