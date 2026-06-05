# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# get imputation summary per chromosome

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

for (i in 1:23) {
  
  if(i==23) {
    tmp_file<-paste(path,"imputed/",study,"/2022/eur/chrX.info.gz",sep="")
  } else {
    tmp_file<-paste(path,"imputed/",study,"/2022/eur/chr",i,".info.gz",sep="")
  }
  
  if (file.exists(tmp_file)) {
    
    tmp<-fread(tmp_file,head=T)
    
    print(i)
    
    tmp<-as.data.frame(tmp)
    
    tmp$chr<-i
    tmp$EmpRsq<-as.numeric(tmp$EmpRsq)
    
    tmp<-tmp[which(tmp$MAF > 0 & tmp$EmpRsq<=0.50),]
  }
  
  if(i==1) {
    df<-tmp
  }else{
    df<-rbind(df,tmp)
  }
  rm(tmp)
}

write.table(df,paste(path,"imputed/",study,"/2022/eur/plots/",study,"_eur_low_EmpRsq_after_first_roudn_exclussions.tsv",sep=""),col.names=T,row.names=F,quote=F,sep="\t")


q("no")

