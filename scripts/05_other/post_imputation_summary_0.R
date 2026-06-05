# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# GET LIST OF DIRECTLY GENTOYPED VARIANTS TO BE EXCLUDED

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(dplyr)
library(stringr)
library(viridis)
library(ggExtra)
library(scales)

rm(list=ls())

path<-"/path/to/ibdgwas/IIBDGC/"
args<-commandArgs()

study<-args[6]
print(study)

chr<-args[7]
print(chr)

# ancestry<-c("eur","noneur")
ancestry<-c("eur")

for (j in 1:length(ancestry)) {

  if(chr==23) {
    df<-fread(paste(path,"imputed/",study,"/2022/",ancestry[j],"/chrX.info.gz",sep=""),head=T)
  } else {
    df<-fread(paste(path,"imputed/",study,"/2022/",ancestry[j],"/chr",chr,".info.gz",sep=""),head=T)
  }
  
  df<-as.data.frame(df)
  
  df$chr<-chr
  df$Rsq<-as.numeric(df$Rsq)
  df$EmpR<-as.numeric(df$EmpR)
  df$EmpRsq<-as.numeric(df$EmpRsq)
 
  # variants neg EmpR to flip:
  df_flip<-df[which(df$EmpR<= -0.5),]
  
  df_flip$SNP<-gsub(":","_",df_flip$SNP)
  df_flip$SNP<-sub("_",":",df_flip$SNP)
  df_flip$SNP<-sub("chrX","chr23",df_flip$SNP)
  df_flip$SNP<-sub("chr","",df_flip$SNP)
  
  write.table(df_flip[,"SNP",drop=F],paste(path,"imputed/",study,"/qc/2022/chr",chr,"_list_variants_negative_EmpR_toflip",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  
  # variants low EmpRsq to remove:
  df_remove<-df[which(df$MAF > 0 & df$EmpRsq<=0.50),]
  
  df_remove$SNP<-gsub(":","_",df_remove$SNP)
  df_remove$SNP<-sub("_",":",df_remove$SNP)
  df_remove$SNP<-sub("chrX","chr23",df_remove$SNP)
  df_remove$SNP<-sub("chr","",df_remove$SNP)
  df_remove<-df_remove[which(!df_remove$SNP %in% df_flip$SNP),]
  
  write.table(df_remove[,"SNP",drop=F],paste(path,"imputed/",study,"/qc/2022/chr",chr,"_list_variants_low_EmpRsq_toexclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

}

q("no")


