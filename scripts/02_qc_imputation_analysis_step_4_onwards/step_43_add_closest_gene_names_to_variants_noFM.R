# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

# add genes in region to the region file:

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=10000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(stringr)
library(rtracklayer)

rm(list=ls())
path_gwas<-"/path/to/ibdgwas/IIBDGC/"

# LOAD THE INDEPENDENT SIGNALS - see /path/to/user/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_41_consolidate_independent_signals_noFM.R
all<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_no_subset_final_2.tsv.gz",sep=""))
all<-as.data.frame(all)
dim(all)
# [1] 619 173

all$min<-gsub("^[0-9]{1,2}_","",all$updated_region)
all$max<-as.numeric(gsub("^[0-9]*_","",all$min))
all$min<-as.numeric(gsub("_[0-9]*$","",all$min))

all$genes_coding<-NA
all$closest_gene<-NA

gtf<-rtracklayer::import(paste(path_gwas,'post_imputation/2022/analysis/metaanalysis/annotation/gencode.v44.annotation.gtf.gz',sep=""))
gene<-as.data.frame(gtf)
rm(gtf)
gene1<-gene[which(gene$type=="exon"),]


for (i in 1:nrow(all)) {
  
  tmp<-gene[which( (gene$seqnames==paste("chr",all$chr[i],sep="")) &
                     ((gene$start<=all$min[i] & gene$end>=all$min[i]) | 
                        (gene$start<=all$max[i] & gene$end>=all$max[i]) | 
                        (gene$start>=all$min[i] & gene$end<=all$max[i])) )
            ,c("gene_name","gene_type","start","end","type","strand")]
  tmp<-tmp[!duplicated(tmp$gene_name),]
  tmp<-tmp[which(tmp$gene_type=="protein_coding"),]
  
  all$genes_coding[i]<-paste(tmp$gene_name[!duplicated(tmp$gene_name)],collapse="|")
  
  if(nrow(tmp)>0) {
    tmp$tss<-NA
    tmp$tss[which(tmp$strand=="+")]<-tmp$start[which(tmp$strand=="+")]
    tmp$tss[which(tmp$strand=="-")]<-tmp$end[which(tmp$strand=="-")]
    
    tmp$dist<-abs(tmp$tss-all$pos[i])
    tmp<-tmp[which(tmp$dist==min(tmp$dist)),]
  }
  
  tmp1<-gene1[which( (gene1$seqnames==paste("chr",all$chr[i],sep="")) &
                       (gene1$start<=all$pos[i] & gene1$end>=all$pos[i]) )
              ,c("gene_name","gene_type","start","end","type","strand")]
  
  tmp1<-tmp1[which(tmp1$gene_type=="protein_coding"),]
  
  if(nrow(tmp1)>0) {
    all$closest_gene[i]<-paste(tmp1$gene_name[!duplicated(tmp1$gene_name)],collapse="|")
  } else if(nrow(tmp)>0) {
    all$closest_gene[i]<-paste(tmp$gene_name[!duplicated(tmp$gene_name)],collapse="|")
  }
  rm(tmp)
}

fwrite(all[,c("MarkerName","closest_gene")],paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_closest_gene.tsv.gz",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

q("no")

