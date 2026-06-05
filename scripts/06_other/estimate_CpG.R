# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# R script to estimate N CpG in +-50Kb from a SNP

library(plyr)
library(ggplot2)
library(data.table)

rm(list=ls())

args<-commandArgs()

n<-Sys.getenv("LSB_JOBINDEX")
print(n)

chr<-args[6]
print(chr)

path<-"/path/to/ibdgwas/IIBDGC/"

keep<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bed",sep=""),head=F)

tmp1<-fread(paste(path,"post_imputation/analysis/metaanalysis/annotation/cpgislands/chr",chr,"_list_union_variants_cd_uc_ibd_metaanalysis_plus_minus_50kb.bed",sep=""),head=F)
tmp1$nCpG<-0
tmp1<-tmp1[which(tmp1$V4 %in% keep$V4),]
print(nrow(tmp1))

tmp2<-fread(paste(path,"post_imputation/analysis/metaanalysis/annotation/cpgislands/chr",chr,"_list_union_variants_cd_uc_ibd_metaanalysis_plus_minus_50kb_cpgIslandExt_merged.bed",sep=""),head=F)
tmp2<-tmp2[which(tmp2$V4 %in% keep$V4),]
print(nrow(tmp2))

rm(keep)

# estimate approx N CpG = ( Number bd overlap (V15) * (Percentage of island that is CpG (V12) /100) ) / 2
tmp2$ncpg_ov<- (tmp2$V15 * (tmp2$V12/100) ) / 2

ids<-tmp2[!duplicated(tmp2$V4),"V4"]
ids<-as.data.frame(ids)
print(nrow(ids))

for(i in 1:nrow(ids)) {

  tmp<-tmp2[which(tmp2$V4==ids$V4[i]),]

  if(nrow(tmp)>0) {
    tmp1$nCpG[which(tmp1$V4==ids$V4[i])]<-sum(tmp$ncpg_ov)
  }else{
     next
  }
}

write.table(tmp1,paste(path,"post_imputation/analysis/metaanalysis/annotation/cpgislands/chr",chr,"_list_union_variants_cd_uc_ibd_metaanalysis_plus_minus_50kb_withCpG.bed",sep=""),
              col.names=T,row.names=F,quote = F,sep="\t")

q("no")
