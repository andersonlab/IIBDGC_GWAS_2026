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
library(ggplot2)

rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

# split by group:

pheno<-c("ibd","cd","uc")

for (ph in pheno) {
    
    for (i in 1:2) {
        tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/top_",i,"_",ph,"_zscore_genes_tissue_output"),head=T)
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
dat$group_f[which(dat$group==1)]<-"low tissue zscore"
dat$group_f[which(dat$group==2)]<-"high tissue zscore"


mu <- dat[which(dat$permutate=="Target_genes"),]
head(mu)
mu$pvalue<-NA

for (ph in pheno) {
    for (i in 1:2) {
        mu$pvalue[which(mu$pheno==ph & mu$group==i)]<-nrow(dat[which(dat$group==i & dat$pheno==ph & dat$permutate!="Target_genes" & dat$mean_r2>=mu$mean_r2[which(mu$pheno==ph & mu$group==i)]),])/nrow(dat[which(dat$group==i & dat$pheno==ph),])
    }
}

print(mu)
#            permutate      mean_r2 group pheno            group_f       pvalue
# 1       Target_genes 2.708317e-06     1   ibd  low tissue zscore 0.0209939580
# 500002  Target_genes 3.051525e-06     2   ibd high tissue zscore 0.0285759428
# 1000003 Target_genes 2.632746e-06     1    cd  low tissue zscore 0.0491539017
# 1500004 Target_genes 3.186143e-06     2    cd high tissue zscore 0.0428079144
# 2000005 Target_genes 2.842964e-06     1    uc  low tissue zscore 0.1594576811
# 2500006 Target_genes 3.184249e-06     2    uc high tissue zscore 0.0003119994

# p<-ggplot(dat, aes(x=mean_r2)) + 
#   geom_histogram(color="black", fill="white") +
#     geom_vline(data=mu,aes(xintercept=mean_r2),linetype="dashed") +
#     geom_text(data=mu,aes(label=formatC(pvalue,digits=3,format="e"),hjust = -0.5),y=80000) + 
#     facet_grid(group_f ~ pheno)

# ggsave(
#   "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_coloc_weights_by_tissue_zscore.png",
#   p,
#   width = 120,
#   height = 60,
#   dpi = 400,
#   units = c("mm"),
#   limitsize = T,scale=2
# )


##########################
# COMBINED

pheno<-c("ibd","cd","uc")

for (ph in pheno) {
        tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/allgenes_",ph,"_zscore_genes_tissue_output"),head=T)
        tmp$pheno<-ph
        tmp$group<-3
        tmp$group_f<-"combined"
        dat<-rbind(dat,tmp)
}


dat<-as.data.frame(dat)
dat$group_f<-factor(dat$group_f,levels=c("high tissue zscore","low tissue zscore","combined"))

mu <- dat[which(dat$permutate=="Target_genes"),]
head(mu)
mu$pvalue<-NA

for (ph in pheno) {
    for (i in 1:3) {
        mu$pvalue[which(mu$pheno==ph & mu$group==i)]<-nrow(dat[which(dat$group==i & dat$pheno==ph & dat$permutate!="Target_genes" & dat$mean_r2>=mu$mean_r2[which(mu$pheno==ph & mu$group==i)]),])/nrow(dat[which(dat$group==i & dat$pheno==ph),])
    }
}
print(mu)
#            permutate      mean_r2 group pheno            group_f       pvalue
# 1       Target_genes 2.708317e-06     1   ibd  low tissue zscore 0.0209939580
# 500002  Target_genes 3.051525e-06     2   ibd high tissue zscore 0.0285759428
# 1000003 Target_genes 2.632746e-06     1    cd  low tissue zscore 0.0491539017
# 1500004 Target_genes 3.186143e-06     2    cd high tissue zscore 0.0428079144
# 2000005 Target_genes 2.842964e-06     1    uc  low tissue zscore 0.1594576811
# 2500006 Target_genes 3.184249e-06     2    uc high tissue zscore 0.0003119994
# 3000007 Target_genes 2.884623e-06     3   ibd           combined 0.0036559927
# 3500008 Target_genes 2.910012e-06     3    cd           combined 0.0113299773
# 4000009 Target_genes 3.019588e-06     3    uc           combined 0.0089839820


p<-ggplot(dat, aes(x=mean_r2)) + 
  geom_histogram(color="black", fill="white") +
    geom_vline(data=mu,aes(xintercept=mean_r2),linetype="dashed") +
    geom_text(data=mu,aes(label=formatC(pvalue,digits=3,format="e"),hjust = -0.5),y=80000) + 
    facet_grid(group_f ~ pheno)

ggsave(
  "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_coloc_weights_by_tissue_zscore_allgenes.png",
  p,
  width = 120,
  height = 90,
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)



##########################
# CLOSEST GENE

pheno<-c("ibd","cd","uc")

for (ph in pheno) {
        tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/coloc_results/test/closest_gene_",ph,"_output"),head=T)
        tmp$pheno<-ph
        tmp$group<-4
        tmp$group_f<-"closest"
        dat<-rbind(dat,tmp)
}


dat<-as.data.frame(dat)
dat$group_f<-factor(dat$group_f,levels=c("high tissue zscore","low tissue zscore","combined","closest"))

mu <- dat[which(dat$permutate=="Target_genes"),]
head(mu)
mu$pvalue<-NA

for (ph in pheno) {
    for (i in 1:4) {
        mu$pvalue[which(mu$pheno==ph & mu$group==i)]<-nrow(dat[which(dat$group==i & dat$pheno==ph & dat$permutate!="Target_genes" & dat$mean_r2>=mu$mean_r2[which(mu$pheno==ph & mu$group==i)]),])/nrow(dat[which(dat$group==i & dat$pheno==ph),])
    }
}
print(mu)
#            permutate      mean_r2 group pheno            group_f       pvalue
# 1       Target_genes 2.708317e-06     1   ibd  low tissue zscore 0.0209939580
# 500002  Target_genes 3.051525e-06     2   ibd high tissue zscore 0.0285759428
# 1000003 Target_genes 2.632746e-06     1    cd  low tissue zscore 0.0491539017
# 1500004 Target_genes 3.186143e-06     2    cd high tissue zscore 0.0428079144
# 2000005 Target_genes 2.842964e-06     1    uc  low tissue zscore 0.1594576811
# 2500006 Target_genes 3.184249e-06     2    uc high tissue zscore 0.0003119994
# 3000007 Target_genes 2.884623e-06     3   ibd           combined 0.0036559927
# 3500008 Target_genes 2.910012e-06     3    cd           combined 0.0113299773
# 4000009 Target_genes 3.019588e-06     3    uc           combined 0.0089839820


p<-ggplot(dat, aes(x=mean_r2)) + 
  geom_histogram(color="black", fill="white") +
    geom_vline(data=mu,aes(xintercept=mean_r2),linetype="dashed") +
    geom_text(data=mu,aes(label=formatC(pvalue,digits=3,format="e"),hjust = -0.5),y=80000) + 
    facet_grid(group_f ~ pheno)

ggsave(
  "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_coloc_weights_by_tissue_zscore_allgenes.png",
  p,
  width = 120,
  height = 90,
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=2
)
