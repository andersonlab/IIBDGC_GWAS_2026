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

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

pheno<-c("ibd","cd","uc")

for (ph in pheno) {
    
    for (i in 1:2) {
        tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/top_",i,"_",ph,"_zscore_genes_tissue_variant_output"),head=T)
        tmp$group<-i
        tmp$pheno<-ph
        if(i==1 & ph=="ibd") {
            dat<-tmp
        } else  {
            dat<-rbind(dat,tmp)
        }
    }
}

dat<-as.data.frame(dat)
dat$group_f<-as.character(dat$group)
dat$group_f[which(dat$group==1)]<-"low tissue*variant zscore"
dat$group_f[which(dat$group==2)]<-"high tissue*variant zscore"


mu <- dat[which(dat$permutate=="Target_genes"),]
head(mu)
mu$pvalue<-NA

for (ph in pheno) {
    for (i in 1:2) {
        mu$pvalue[which(mu$pheno==ph & mu$group==i)]<-nrow(dat[which(dat$group==i & dat$pheno==ph & dat$permutate!="Target_genes" & dat$mean_r2>=mu$mean_r2[which(mu$pheno==ph & mu$group==i)]),])/nrow(dat[which(dat$group==i & dat$pheno==ph),])
    }
}

print(mu)
#            permutate      mean_r2 group pheno                    group_f
# 1       Target_genes 2.728086e-06     1   ibd  low tissue*variant zscore
# 500002  Target_genes 3.035888e-06     2   ibd high tissue*variant zscore
# 1000003 Target_genes 2.480170e-06     1    cd  low tissue*variant zscore
# 1500004 Target_genes 3.334592e-06     2    cd high tissue*variant zscore
# 2000005 Target_genes 2.915787e-06     1    uc  low tissue*variant zscore
# 2500006 Target_genes 3.114773e-06     2    uc high tissue*variant zscore
#              pvalue
# 1       0.008893982
# 500002  0.040905918
# 1000003 0.117319765
# 1500004 0.025209950
# 2000005 0.134373731
# 2500006 0.001751996

p<-ggplot(dat, aes(x=mean_r2)) + 
  geom_histogram(color="black", fill="white") +
    geom_vline(data=mu,aes(xintercept=mean_r2),linetype="dashed") +
    geom_text(data=mu,aes(label=formatC(pvalue,digits=3,format="e"),hjust = -0.5),y=80000) + 
    facet_grid(group_f ~ pheno)

ggsave(
  "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_coloc_weights_by_tissue_variant_zscore.png",
  p,
  width = 120,
  height = 60,
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)
