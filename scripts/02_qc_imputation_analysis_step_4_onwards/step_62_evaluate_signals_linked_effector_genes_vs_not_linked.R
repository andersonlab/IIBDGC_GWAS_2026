# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
############################################################################################################################
############################################################################################################################


# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=6000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q yesterday R

library(data.table)
library(ggplot2)
library(colorspace)
library(dplyr)
library(tibble)
library(ggrepel)
library(ggpubr)
library(rtracklayer)
rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

all<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants_join_conditinal_model.tsv.gz",sep=""))
dim(all)
# [1] 619 205

list_genes<-fread(paste0(path_gwas,"post_imputation/2022/analysis/final_tables/list_effector_genes_coloc_mr_exonic_with_coloc_rank1_noFM_with_otar_evidence_with_ligand_receptor_with_protein_complexes_with_monogenic_with_gps.tsv.gz"))
vec<-unlist(strsplit(list_genes$independent_index_variants,";"))


summary(all$A2FREQ_CONTROLS_ibd_eur_tier_2)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001065 0.108627 0.266899 0.311333 0.488174 0.984240 

all$maf<-pmin(all$A2FREQ_CONTROLS_ibd_eur_tier_2,1-all$A2FREQ_CONTROLS_ibd_eur_tier_2)

summary(all$maf)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001065 0.104442 0.248914 0.242480 0.380063 0.498348 

summary(all$maf[which(all$MarkerName %in% vec)])
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001065 0.146346 0.265632 0.265557 0.388178 0.498348

summary(all$maf[which(!all$MarkerName %in% vec)])
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001187 0.057000 0.202559 0.214475 0.359032 0.497298


all$effector_gene_nominated<-0
all$effector_gene_nominated[which(all$MarkerName %in% vec)]<-1

all$effector_gene_nominated<-as.numeric(all$effector_gene_nominated)

table(all$effector_gene_nominated)
#   n   y 
# 264 355

summary(glm(effector_gene_nominated ~ maf,data = all,family = "binomial"))
# Coefficients:
#             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -0.3339     0.1479  -2.258    0.024 *  
# maf           1.9263     0.5372   3.586 0.000336 ***

all$effector_gene_nominated<-as.character(all$effector_gene_nominated)
wilcox.test(maf ~ effector_gene_nominated, data=all) 
# data:  maf by effector_gene_nominated
# W = 42738, p-value = 0.0003679

# ggplot(all,aes(x=effector_gene_nominated,y=maf)) +
# geom_boxplot()

my_comparisons <- list( c("0", "1"))


p1<-ggplot(all,aes(x=effector_gene_nominated,y=maf)) +
geom_boxplot() + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",method.args = list(alternative = "two.sided")) + theme_bw() + ggtitle("MAF lead variants\nwith/without effector genes nominated\n")


### explore this for the coloc only:
vec_coloc<-unlist(strsplit(list_genes$independent_index_variants[which(list_genes$coloc==1)],";"))
vec_no_coloc<-unlist(strsplit(list_genes$independent_index_variants[which(list_genes$coloc!=1)],";"))


all$effector_gene_nominated_coloc<-0
all$effector_gene_nominated_coloc[which(all$MarkerName %in% vec_coloc)]<-1
all$effector_gene_nominated_coloc[which(all$MarkerName %in% vec_no_coloc)]<-NA

p2<-ggplot(all[which(!is.na(all$effector_gene_nominated_coloc))],aes(x=effector_gene_nominated,y=maf)) +
geom_boxplot() + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",method.args = list(alternative = "two.sided")) + theme_bw() + ggtitle("MAF Lead variants\nwith/without effector genes nominated\n(Colocalization only)")


p12<-ggarrange(p1,p2)



### add imputation info:

# EUR tier2:

chr<-as.numeric(unlist(labels(table(all$chr))))
chr<-chr[which(!is.na(chr))]

pheno<-c("ibd","cd","uc")

for (i in 1:length(pheno)) {
  
  print(pheno[i])
  
  for (j in 1:length(chr)) {
    
    tmp1<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr[j],"_",pheno[i],"_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz",sep=""),head=T)
    tmp1<-tmp1[which(tmp1$MarkerName %in% all$MarkerName),]

    tmp1<-tmp1[,c("MarkerName","INFO")]
    colnames(tmp1)<-c("MarkerName","INFO")
    
    if(j==1) {
      tmp<-tmp1
    } else {
      tmp<-rbind(tmp,tmp1)
    }
    
    rm(tmp1)
    
  }

  colnames(tmp)[2]<-paste(colnames(tmp)[2],pheno[i],"eur_tier_2",sep="_")
  
  if(i==1) {
    dat<-tmp
  } else {
    dat<-merge(dat,tmp,by="MarkerName",sort=F)
  }
  rm(tmp)
}

all<-merge(all,dat,by="MarkerName",all.x=T)



p3<-ggplot(all,aes(x=effector_gene_nominated,y=INFO_ibd_eur_tier_2)) +
geom_boxplot() + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",method.args = list(alternative = "two.sided")) + theme_bw() + ggtitle("INFO lead variants\nwith/without effector genes nominated\n") + ylim(0.5,1.1)

p4<-ggplot(all[which(!is.na(all$effector_gene_nominated_coloc))],aes(x=effector_gene_nominated,y=INFO_ibd_eur_tier_2)) +
geom_boxplot() + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",method.args = list(alternative = "two.sided")) + theme_bw() + ggtitle("INFO Lead variants\nwith/without effector genes nominated\n(Colocalization only)") + ylim(0.5,1.1)

p34<-ggarrange(p3,p4)


#  NEFF rate

p5<-ggplot(all,aes(x=effector_gene_nominated,y=rate_Neff_ibd_eur_tier_2)) +
geom_boxplot() + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",method.args = list(alternative = "two.sided")) + theme_bw() + ggtitle("Neff rate lead variants\nwith/without effector genes nominated\n") + ylim(0.3,1.1)

p6<-ggplot(all[which(!is.na(all$effector_gene_nominated_coloc))],aes(x=effector_gene_nominated,y=rate_Neff_ibd_eur_tier_2)) +
geom_boxplot() + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",method.args = list(alternative = "two.sided")) + theme_bw() + ggtitle("Neff rate Lead variants\nwith/without effector genes nominated\n(Colocalization only)") + ylim(0.3,1.1)

p56<-ggarrange(p5,p6)

#### ARE THOSE LOCATED FURTHER AWAY FROM TSS?

# ADD CLOSEST TSS:

# merge with gene ID:
gtf<-rtracklayer::import(paste(path_gwas,'post_imputation/2022/analysis/metaanalysis/annotation/gencode.v44.annotation.gtf.gz',sep=""))
gene<-as.data.frame(gtf)

rm(gtf)
gene1<-gene[which(gene$type=="transcript" & gene$transcript_type=="protein_coding"),]
dim(gene1)
# [1] 89067    27

gene1$tss<-NA
gene1$tss[which(gene1$strand=="+")]<-gene1$start
gene1$tss[which(gene1$strand=="-")]<-gene1$end

all$closest_gene_tss<-NA
all$closest_gene_tss_pos<-NA

for (i in 1:nrow(all)) {

  tmp1<-gene1[which( (gene1$seqnames==paste("chr",all$chr[i],sep=""))),]

  tmp1$distance_to_tss<-NA
  tmp1$distance_to_tss<-abs(tmp1$tss-all$pos[i])

  tmp1<-tmp1[which(tmp1$distance_to_tss==min(tmp1$distance_to_tss)),c("gene_name","tss","distance_to_tss")]
  tmp1<-tmp1[!duplicated(tmp1),]
  
  if(nrow(tmp1)>0) {
    all$closest_gene_tss[i]<-paste(tmp1$gene_name,collapse="|")
    all$closest_gene_tss_pos[i]<-tmp1$distance_to_tss[1]
  }

  rm(tmp1)
}



p7<-ggplot(all,aes(x=effector_gene_nominated,y=closest_gene_tss_pos)) +
geom_boxplot() + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",method.args = list(alternative = "two.sided")) + theme_bw() + ggtitle("Distance closest TSS lead variants\nwith/without effector genes nominated\n")

p8<-ggplot(all[which(!is.na(all$effector_gene_nominated_coloc))],aes(x=effector_gene_nominated,y=closest_gene_tss_pos)) +
geom_boxplot() + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",method.args = list(alternative = "two.sided")) + theme_bw() + ggtitle("Distance closest TSS lead variants\nwith/without effector genes nominated\n(Colocalization only)")

p78<-ggarrange(p7,p8)

### effect size

all$zscore<-all$BETA_ibd_eur_tier_2/all$SE_ibd_eur_tier_2
all$zscore[which(all$pheno=="CD")]<-all$BETA_cd_eur_tier_2[which(all$pheno=="CD")]/all$SE_cd_eur_tier_2[which(all$pheno=="CD")]
all$zscore[which(all$pheno=="UC")]<-all$BETA_uc_eur_tier_2[which(all$pheno=="UC")]/all$SE_uc_eur_tier_2[which(all$pheno=="UC")]



p9<-ggplot(all,aes(x=effector_gene_nominated,y=zscore)) +
geom_boxplot() + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",method.args = list(alternative = "two.sided")) + theme_bw() + ggtitle("Effect zscore lead variants\nwith/without effector genes nominated\n")

p10<-ggplot(all[which(!is.na(all$effector_gene_nominated_coloc))],aes(x=effector_gene_nominated,y=zscore)) +
geom_boxplot() + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",method.args = list(alternative = "two.sided")) + theme_bw() + ggtitle("Effect zscore lead variants\nwith/without effector genes nominated\n(Colocalization only)")

p910<-ggarrange(p9,p10)



p<-ggarrange(p12,p34,p56,p78,p910,ncol=2,nrow=3,labels=c("a","b","c","d","e"))


ggsave(
  paste0("~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_lead_variants_maf_info_neff_by_effector_genes_nomination.png"),
  p,
  width = 250,
  height = 200,
  dpi = 600,
  units = c("mm"),
  limitsize = T,scale=2
)

q("no")