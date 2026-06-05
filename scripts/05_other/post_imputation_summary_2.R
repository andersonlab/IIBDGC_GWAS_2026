# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

# GET SUMMARY DATA FOR ALL CHROMOSOMES TOGETHER

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

# ancestry<-c("eur","noneur")
ancestry<-c("eur")

for (j in 1:length(ancestry)) {
  
  for (i in 1:23) {
    
    if(i==23) {
      tmp<-fread(paste(path,"imputed/",study,"/",ancestry[j],"/chrX.info.gz",sep=""),head=T)
    } else {
      tmp<-fread(paste(path,"imputed/",study,"/",ancestry[j],"/chr",i,".info.gz",sep=""),head=T)
    }
    
    print(i)
    
    tmp<-as.data.frame(tmp)
    
    tmp$chr<-i
    tmp$Rsq<-as.numeric(tmp$Rsq)
    
    ##################################################
    # create table summary MAF for all chrs:

    sum_tmp_maf<-tmp %>%
      group_by(MAF_range = cut(MAF, breaks = c(0,0.0001,0.001,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.501), 
                               right = F)) %>% 
      summarise(MAF_n = n())
    
    sum_tmp_maf<-as.data.frame(sum_tmp_maf)
    colnames(sum_tmp_maf)[2]<-paste("N_chr",i,sep="")
    sum_tmp_maf$MAF_range<-as.character(sum_tmp_maf$MAF_range)
    sum_tmp_maf$MAF_range[nrow(sum_tmp_maf)]<-'[0.45,0.50]'
    
    if(i==1){
      maf_dist<-sum_tmp_maf
    } else {
      maf_dist<-merge(maf_dist,sum_tmp_maf,by="MAF_range",sort=F)
    }

    ##################################################
    # create table summary MAF for variants Rsq>=0.4 across all chrs:
    
    tmp1<-tmp[which(tmp$Rsq>=0.4),]
    sum_tmp_maf2<-tmp1 %>%
      group_by(MAF_range = cut(MAF, breaks = c(0,0.0001,0.001,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.501), 
                               right = F)) %>% 
      summarise(MAF_n = n())
    
    sum_tmp_maf2<-as.data.frame(sum_tmp_maf2)
    colnames(sum_tmp_maf2)[2]<-paste("N_chr",i,sep="")
    sum_tmp_maf2$MAF_range<-as.character(sum_tmp_maf2$MAF_range)
    sum_tmp_maf2$MAF_range[nrow(sum_tmp_maf2)]<-'[0.45,0.50]'
    
    if(i==1){
      maf_dist_2<-sum_tmp_maf2
    } else {
      maf_dist_2<-merge(maf_dist_2,sum_tmp_maf2,by="MAF_range",sort=F)
    }
    
    ##################################################
    # create table Rsqr for all chrs:

    sum_tmp_rsq<-tmp %>%
      group_by(Rsq_range = cut(Rsq, breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.01), 
                               right = F)) %>% 
      summarise(Rsq_n = n())
    
    sum_tmp_rsq<-as.data.frame(sum_tmp_rsq)
    sum_tmp_rsq<-sum_tmp_rsq[which(!is.na(sum_tmp_rsq$Rsq_range)),]
    
    colnames(sum_tmp_rsq)[2]<-paste("N_chr",i,sep="")
    sum_tmp_rsq$Rsq_range<-as.character(sum_tmp_rsq$Rsq_range)
    sum_tmp_rsq$Rsq_range[nrow(sum_tmp_rsq)]<-'[0.9,1.0]'
    
    
    if(i==1){
      rsq_dist<-sum_tmp_rsq
    } else {
      rsq_dist<-merge(rsq_dist,sum_tmp_rsq,by="Rsq_range",sort=F)
    }
    
    ##################################################
    # create table Genotyped vs Imputed for all chrs:

    x<-table(tmp$Genotyped)
    x<-as.data.frame.matrix(t(as.matrix(x)))
    x$chr<-i

    if (ncol(x)==3) {
      x$Typed_Only<-0
    }

    if(i==1){
      gtImp<-x
    } else {
      gtImp<-rbind(gtImp,x)
    }

    ##################################################

    tmp<-tmp[which(tmp$Genotyped=="Genotyped"),]
    
    if(i==1) {
      df<-tmp
    }else{
      df<-rbind(tmp,df)
    }
    
    rm(tmp,tmp1)
    
  }
  
  df$EmpR<-as.numeric(df$EmpR)
  df$EmpRsq<-as.numeric(df$EmpRsq)

  
  maf_dist$all_chr<-NA
  maf_dist$all_chr<-rowSums(maf_dist[,2:(ncol(maf_dist)-1)])
  
  rsq_dist$all_chr<-NA
  rsq_dist$all_chr<-rowSums(rsq_dist[,2:(ncol(rsq_dist)-1)])
  
  write.table(maf_dist,paste(path,"imputed/",study,"/",ancestry[j],"/plots/",study,"_",ancestry[j],"_maf_distribution.tsv",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  write.table(maf_dist_2,paste(path,"imputed/",study,"/",ancestry[j],"/plots/",study,"_",ancestry[j],"_maf_distribution_variants_rsq_0.4.tsv",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  write.table(rsq_dist,paste(path,"imputed/",study,"/",ancestry[j],"/plots/",study,"_",ancestry[j],"_rsq_distribution.tsv",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  write.table(gtImp,paste(path,"imputed/",study,"/",ancestry[j],"/plots/",study,"_",ancestry[j],"_genotyped_imputed_variants_distribution.tsv",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

  
  ##################################################
  # plot EmpR and EmpRsq  vs MAF for all chrs (similar to Kyle plots)
  
  p1<-ggplot(df, aes(x = MAF, y = EmpR)) + ggtitle(study) + ylim(-1,1) + ylab("Empirical R") + xlab("MAF") +
    geom_point(alpha = 0.3,shape=20,size=1)
  
  p2<-ggplot(df, aes(x = MAF, y = EmpRsq)) + ggtitle(study) + ylim(0,1) + ylab("Empirical Rsq") + xlab("MAF") +
    geom_point(alpha = 0.3,shape=20,size=1)
  
  pdf(paste(path,"imputed/",study,"/",ancestry[j],"/plots/",study,"_",ancestry[j],"_maf_EmpRsq.pdf",sep=""),width = 14, height = 5)
  print(ggarrange(ggMarginal(p1, type="density",margins = 'y',colour="grey",fill = "grey"),
                  ggMarginal(p2, type="density",margins = 'y',colour="grey",fill = "grey"),ncol=2))
  dev.off()
  
  rm(p1,p2)
  
  ##################################################
  # EmpRsqr mean bins by MAF:

  scientific <- function(x){
    ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
  }
  
  sum_df<-df %>%
    group_by(MAF_range = cut(MAF, breaks = c(0,5E-6,2.5E-5,1.25E-4,2.5E-4,5E-4,7.5E-4,1E-3,2.5E-3,5E-3,7.5E-3,1E-2,5E-2,2.5E-2,1E-1,2E-1,3E-1,4E-1,1E0), 
                             right = F)) %>% 
    summarise(EmpRsq_mean = mean(EmpRsq),
              EmpRsq_n = n())
  
  sum_df<-as.data.frame(sum_df)
  sum_df$group<-"Genotyped variants"
  sum_df$tmp_maf<-gsub("\\[","",sum_df$MAF_range)
  sum_df$tmp_maf<-gsub(")","",sum_df$tmp_maf)
  
  xx<-as.data.frame(str_split(sum_df$tmp_maf,",",simplify=T))
  xx$V1<-as.numeric(as.character(xx$V1))
  xx$V2<-as.numeric(as.character(xx$V2))
  
  sum_df$MAF<-NA
  for(xxx in 1:nrow(xx)) {
    sum_df$MAF[xxx]<-sum(xx$V1[xxx],xx$V2[xxx])/2
  }
  
  p1<-ggplot(sum_df, aes(x = -log10(MAF), y = EmpRsq_mean, group = group, color = group)) + ylim(0,1) + ggtitle(study) + ylab("Mean Empirical RSQ") +
    geom_point() + geom_line() + expand_limits(x=c(-log10(1E-5),-log10(1E0)) ) + scale_x_continuous(name ="MAF",trans='reverse', 
                                                                                                    breaks= -log10(c(1E-5,1E-4,1E-3,1E-2,1E-1,1E0)) ,
                                                                                                    labels=scientific(c(1E-5,1E-4,1E-3,1E-2,1E-1,1E0)))
  
  pdf(paste(path,"imputed/",study,"/",ancestry[j],"/plots/",study,"_",ancestry[j],"_maf_EmpRsq_mean.pdf",sep=""),width = 7, height = 4)
  print(p1)
  dev.off()
  
  write.table(sum_df,paste(path,"imputed/",study,"/",ancestry[j],"/plots/",study,"_",ancestry[j],"_maf_EmpRsq_mean.tsv",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  
  # variants neg EmpR to flip:
  df_flip<-df[which(df$EmpR<= -0.5),]
  
  df_flip$SNP<-gsub(":","_",df_flip$SNP)
  df_flip$SNP<-sub("_",":",df_flip$SNP)
  df_flip$SNP<-sub("chrX","chr23",df_flip$SNP)
  df_flip$SNP<-sub("chr","",df_flip$SNP)
  
  write.table(df_flip[,"SNP",drop=F],paste(path,"imputed/",study,"/qc/list_variants_negative_EmpR_toflip",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  
  
  
  # variants low EmpRsq to remove:
  df_remove<-df[which(df$MAF > 0 & df$EmpRsq<=0.25),]
  
  df_remove$SNP<-gsub(":","_",df_remove$SNP)
  df_remove$SNP<-sub("_",":",df_remove$SNP)
  df_remove$SNP<-sub("chrX","chr23",df_remove$SNP)
  df_remove$SNP<-sub("chr","",df_remove$SNP)
  df_remove<-df_remove[which(!df_remove$SNP %in% df_flip$SNP),]
  
  write.table(df_remove[,"SNP",drop=F],paste(path,"imputed/",study,"/qc/list_variants_low_EmpRsq_toexclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

}


  
q("no")


