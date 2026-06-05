# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# retrieve a curated list of receptors and ligands, and explore how many of these pairs are in our list of nearby genes/genes in locus - later explore likely target genes

# database used as ref: http://tcm.zju.edu.cn/celltalkdb/
# paper describing database: https://academic.oup.com/bib/article-abstract/22/4/bbaa269/5955941?redirectedFrom=fulltext&login=true#no-access-message

# version downloaded:
#
# File	Description	Statistics	Download times	Latest Download	License	Size
# human_lr_pair.txt Human ligand-receptor interaction pairs	3,399 lines, 10 columns.	2624 times	2022-10-04 18:24	GPL-3.0	355.15 Kb

# Setup (run in shell before starting R):
# # singularity exec iibdgc_postprocess_10_singularity.sif
# MEM=6000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -q cpu-interactive -G your_hpc_group R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(stringr)

path<-"/path/to/ibdgwas/IIBDGC/"


# LOAD THE INDEPENDENT SIGNALS - see ~/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_38_2_define_independent_signals_using_as_ref_cojo_agnostic_results.R
all<-fread(paste(path,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants.tsv.gz",sep=""))
all<-as.data.frame(all)
dim(all)
# [1] 619 211


reli<-read.table("/path/to/ibdgwas/IIBDGC/resources/celltalkdb/human_lr_pair.txt",head=T)

# eqtl/coloc genes
all_genes<-fread(paste0(path,"post_imputation/2022/analysis/final_tables/list_effector_genes_coloc_mr_exonic_with_coloc_rank1_noFM_with_otar_evidence_with_coloc_mr_direction.tsv.gz"))
all_genes<-all_genes$list_genes

reli$ligand_ibd_mr_coloc_exonic_coding<-0
reli$receptor_ibd_mr_coloc_exonic_coloc_coding<-0

reli$ligand_ibd_mr_coloc_exonic_coding[which(reli$ligand_gene_symbol %in% all_genes)]<-1
reli$receptor_ibd_mr_coloc_exonic_coloc_coding[which(reli$receptor_gene_symbol %in% all_genes)]<-1

table(reli$ligand_ibd_mr_coloc_exonic_coding,reli$receptor_ibd_mr_coloc_exonic_coloc_coding)
  #     0    1
  # 0 3022  228
  # 1  125   23

reli[which(reli$ligand_ibd_mr_coloc_exonic_coding==1 & reli$receptor_ibd_mr_coloc_exonic_coloc_coding==1),c("lr_pair")]
#  [1] "ICAM5_ITGAL"       "CCL2_CCR2"         "CCL2_CCR5"        
#  [4] "CCL13_CCR2"        "CCL13_CCR5"        "VTN_ITGAV"        
#  [7] "VTN_PLAUR"         "VTN_TNFRSF11B"     "VCAM1_ITGA4"      
# [10] "CXCL5_CXCR2"       "ADAM17_MUC1"       "ICAM4_ITGAL"      
# [13] "ICAM4_ITGA4"       "GNAS_ADCY7"        "TNFSF15_TNFRSF25" 
# [16] "VTN_ITGB8"         "PLAU_ITGAV"        "PDGFB_ITGAV"      
# [19] "ICAM4_ITGAV"       "ADAM15_ITGAV"      "CXCL5_CXCR1"      
# [22] "TNFSF11_TNFRSF11B" "PLAU_PLAUR"   


### probability of getting same N genes as per coloc being 

gtf<-rtracklayer::import(paste(path,'post_imputation/2022/analysis/metaanalysis/annotation/gencode.v44.annotation.gtf.gz',sep=""))
gene<-as.data.frame(gtf)
gene<-gene[which(!gene$seqnames %in% c("chrX","chrY","chrM")),]
gene<-gene[which(gene$type=="exon"),c("gene_name","gene_type","gene_id")]

gene<-gene[!duplicated(gene),]
gene<-gene[which(gene$gene_type=="protein_coding"),]

#######################################################
# get random N genes and check how many are dups:

nrow(gene)
# [1] 19111


iter<-as.data.frame(matrix(ncol=4,nrow=500000))
colnames(iter)<-c("N_iter","N_pairs_receptor_ligand","N_ligand","N_receptor")

# Pre-extract vectors to avoid repeated data.frame column access inside loop
ligand_vec   <- reli$ligand_gene_symbol
receptor_vec <- reli$receptor_gene_symbol
gene_names   <- gene[[1]]
n_sample     <- length(all_genes)

for (i in 1:500000) {

  gene_tmp     <- sample(gene_names, n_sample, replace=FALSE)
  ligand_hit   <- ligand_vec  %in% gene_tmp
  receptor_hit <- receptor_vec %in% gene_tmp

  iter$N_iter[i]                <- i
  iter$N_pairs_receptor_ligand[i] <- sum(ligand_hit & receptor_hit)
  iter$N_ligand[i]              <- sum(ligand_hit)
  iter$N_receptor[i]            <- sum(receptor_hit)

}

table(iter$N_pairs_receptor_ligand)

#### Ligand - Receptor pairs:
# pvalue
nrow(iter[which(iter$N_pairs_receptor_ligand>=nrow(reli[which(reli$ligand_ibd_mr_coloc_exonic_coding==1 & reli$receptor_ibd_mr_coloc_exonic_coloc_coding==1),])),])/nrow(iter)
# [1] 2.2e-05

# OR
nrow(reli[which(reli$ligand_ibd_mr_coloc_exonic_coding==1 & reli$receptor_ibd_mr_coloc_exonic_coloc_coding==1),])/mean(iter$N_pairs_receptor_ligand)
# [1] 5.913274


#### Ligands:
# pvalue
nrow(iter[which(iter$N_ligand>=nrow(reli[which(reli$ligand_ibd_mr_coloc_exonic_coding==1),])),])/nrow(iter)
# [1] 0.142376

# OR
nrow(reli[which(reli$ligand_ibd_mr_coloc_exonic_coding==1),])/mean(iter$N_ligand)
# [1] 1.288257


#### Receptors:
# pvalue
nrow(iter[which(iter$N_receptor>=nrow(reli[which(reli$receptor_ibd_mr_coloc_exonic_coloc_coding==1),])),])/nrow(iter)
# [1]0.005014

# OR
nrow(reli[which(reli$receptor_ibd_mr_coloc_exonic_coloc_coding==1),])/mean(iter$N_receptor)
# [1]2.17783


reli[which(reli$ligand_ibd_mr_coloc_exonic_coding==1 & reli$receptor_ibd_mr_coloc_exonic_coloc_coding==1),c("lr_pair")]
#  [1] "ICAM5_ITGAL"       "CCL2_CCR2"         "CCL2_CCR5"        
#  [4] "CCL13_CCR2"        "CCL13_CCR5"        "VTN_ITGAV"        
#  [7] "VTN_PLAUR"         "VTN_TNFRSF11B"     "VCAM1_ITGA4"      
# [10] "CXCL5_CXCR2"       "ADAM17_MUC1"       "ICAM4_ITGAL"      
# [13] "ICAM4_ITGA4"       "GNAS_ADCY7"        "TNFSF15_TNFRSF25" 
# [16] "VTN_ITGB8"         "PLAU_ITGAV"        "PDGFB_ITGAV"      
# [19] "ICAM4_ITGAV"       "ADAM15_ITGAV"      "CXCL5_CXCR1"      
# [22] "TNFSF11_TNFRSF11B" "PLAU_PLAUR"


reli[which(reli$ligand_ibd_mr_coloc_exonic_coding==1 & reli$receptor_ibd_mr_coloc_exonic_coloc_coding==0 & reli$receptor_ibd_region==1),]


reli[which(reli$ligand_ibd_mr_coloc_exonic_coding==0 & reli$ligand_ibd_region==1 & reli$receptor_ibd_mr_coloc_exonic_coloc_coding==1),]

# get receptors whose ligands activate downstreak JAK/STAT
# used as ref A. Salas, C. Hernandez-Rocha, M. Duijvestein, W. Faubion, D. McGovern, S. Vermeire, et al.
# Nat Rev Gastroenterol Hepatol 2020 Vol. 17 Issue 6 Pages 323-337
# Accession Number: 32203403 DOI: 10.1038/s41575-020-0273-0

# CSF3 encodes for G-CSF
# CSF2 encodes for GM-CSF
# GH1 encodes for GH
# LEP encodes for leptyn

# ligands<-c("IL2","IL4","IL7","IL9","IL15","IL21","IL12","IL23","IL10","IL19","IL20","IL22","IL26","EPO","TPO","CSF3","CSF2","GH1","LEP",
#            "IL3","IL5","IL6","IL11","IL13","OSM","LIF",
#            "IFNA1","IFNA10","IFNA13","IFNA14","IFNA16","IFNA17","IFNA2","IFNA21","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8",
#            "IFNB1","IFNG")

# reli$jak_stat<-0
# reli$jak_stat[which(reli$ligand_gene_symbol %in% ligands)]<-1

# reli[which(reli$jak_stat==1 & (reli$ligand_ibd_mr_coloc_exonic_coding==1 | reli$receptor_ibd_mr_coloc_exonic_coloc_coding==1)),
#      c("lr_pair","ligand_gene_symbol","receptor_gene_symbol",
#        "ligand_ibd_mr_coloc_exonic_coding","receptor_ibd_mr_coloc_exonic_coloc_coding","ligand_ibd_region","receptor_ibd_region")]


write.table(reli,"~/git/IIBDGC_GWAS/plots/paper_tables/Source_data_human_lr_pair_edited.txt",
            col.names=T,row.names=F,quote=F,sep="\t")


# eqtl/coloc genes
list_genes<-fread(paste0(path,"post_imputation/2022/analysis/final_tables/list_effector_genes_coloc_mr_exonic_with_coloc_rank1_noFM_with_otar_evidence_with_coloc_mr_direction.tsv.gz"))

list_genes$ligand_gene<-0
list_genes$ligand_gene[which(list_genes$list_genes %in% reli$ligand_gene_symbol)]<-1
list_genes$receptor_gene<-0
list_genes$receptor_gene[which(list_genes$list_genes %in% reli$receptor_gene_symbol)]<-1

table(list_genes$ligand_gene)
#   0   1 
# 691  26 
table(list_genes$receptor_gene)
#   0   1 
# 671  47

write.table(list_genes,paste0(path,"post_imputation/2022/analysis/final_tables/list_effector_genes_coloc_mr_exonic_with_coloc_rank1_noFM_with_otar_evidence_with_coloc_mr_direction_with_ligand_receptor.tsv.gz"),
            col.names=T,row.names=F,quote=F,sep="\t")

q("no")





##################

# plot these and /path/to/user/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_57_define_receptor_ligands_pairs.R outputs:
# ~/git/IIBDGC_GWAS/scripts/other/plot_ligand_receptor_pairs_and_protein_complexes.R

