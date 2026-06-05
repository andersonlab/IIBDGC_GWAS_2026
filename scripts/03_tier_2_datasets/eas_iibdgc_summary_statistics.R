# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# Retrieve data from Liu et al. 2022
path_gwas="/path/to/ibdgwas/IIBDGC/"

  
  
cd ${path_gwas}post_imputation/2022/analysis/stage_2_summary_statistics/eas_iibdgc/
tar -xzvf liu-2022-east-asian-gwas.tar.gz 

cd liu-2022-east-asian-gwas/

less readme.sst.txt

# # version 12202022 
# 
# EAS fixed-effect meta-analysis summary statistics:

# 7,372 CD and 15,456 Control
# ibd_SiKJ_meta_CD.TBL.meta.txt.gz
# 6,862 UC and 15,456 Control
# ibd_SiKJ_meta_UC.TBL.meta.txt.gz
# 14,393 IBD and 15,456 Control
# ibd_SiKJ_meta_IBD.TBL.meta.txt.gz
# 
# ...
# 
# Format:
# MarkerName: variant (chr:pos:ref:alt) (GRCh38/hg38)
# Allele1: reference allele
# Allele2: effect allele
# Freq1: mean effect allele frequency across cohorts
# FreqSE: standard error of effect allele frequency across cohorts
# MinFreq: minimum of effect allele frequency across cohorts
# MaxFreq: maximum of effect allele frequency across cohorts
# Effect: effect size
# StdErr: standard error
# P-value: P value

# Direction: Direction of effect in each cohort
# (order for EAS: SHA1, ICH1, KOR1, JPN1) - infer sample size per SNP:



# MEM=15000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1

library(data.table)
path_gwas<-"/path/to/ibdgwas/IIBDGC/"
phenotypes<-c("ibd","cd","uc")


# See paper Fig. 1 to find out where N comes from:

total_N_all<-as.data.frame(matrix(ncol=4,nrow=12))
colnames(total_N_all)<-c("study","n_cases","n_controls","pheno")
total_N_all$study<-rep(c("SHA1","ICH1","KOR1","JPN1"),3)
total_N_all$pheno<-c(rep("ibd",4),rep("cd",4),rep("uc",4))

study<-c("SHA1","ICH1","KOR1","JPN1")
total_N_all$n_cases<-c(5088,2735,3188,3382,2552,1611,1619,1590,2400,1124,1569,1769)
total_N_all$n_controls<-rep(c(6279,3724,4419,1034),3)

for (pheno in phenotypes) {
  
  print(pheno)
  
  dat<-fread(paste(path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/eas_iibdgc/liu-2022-east-asian-gwas/summary-stats/ibd_EAS_SiKJ_meta_",toupper(pheno),".TBL.txt.gz",sep=""),head=T)
  tmp<-as.data.frame(matrix(unlist(strsplit(dat$Direction,"")),nrow=nrow(dat),ncol=4,byrow=T))
  
  dat<-cbind(dat,tmp)
  rm(tmp)
  
  dat$N<-NA
  dat$N_CASES<-NA
  dat$N_CONTROLS<-NA
  
  dat<-as.data.frame(dat)
  
  if(pheno=="cd") {
    
    dat$N<-sum(7372,15456)
    dat$N_CASES<-7372
    dat$N_CONTROLS<-15456
    

   } else if (pheno=="ibd") {
    
    dat$N<-sum(14393,15456)
    dat$N_CASES<-14393
    dat$N_CONTROLS<-15456
 
  } else if (pheno=="uc") {
    
    dat$N<-sum(6862,15456)
    dat$N_CASES<-6862
    dat$N_CONTROLS<-15456

  }

  # exclude samples from studies that do not contribute to meta-analysis
  for (i in c(1:4)) {
    
    total_N<-total_N_all[which(total_N_all$study==study[i] & total_N_all$pheno==pheno),]
    
    dat$N[which(dat[,17+i]=="?")]<-dat$N[which(dat[,17+i]=="?")]-sum(total_N[,2:3])
    dat$N_CASES[which(dat[,17+i]=="?")]<-dat$N_CASES[which(dat[,17+i]=="?")]-sum(total_N[,2])
    dat$N_CONTROLS[which(dat[,17+i]=="?")]<-dat$N_CONTROLS[which(dat[,17+i]=="?")]-sum(total_N[,3])
    
  }
  
  dat$CHISQ<-(dat$Effect/dat$StdErr)^2
  
  dat$A1FREQ_CASES<-NA
  dat$A1FREQ_CONTROLS<-NA
  dat$INFO<-1
  dat$TEST<-"add"
  dat$EXTRA<-NA
  
  dat$ALLELE1<-toupper(dat$Allele2)
  dat$ALLELE0<-toupper(dat$Allele1)
  
  dat$LOG10P<- -log10(dat$'P-value')

  dat<-dat[,c("CHR","BP","MarkerName","ALLELE0","ALLELE1","Freq1","A1FREQ_CASES","A1FREQ_CONTROLS","INFO","N","N_CASES","N_CONTROLS","TEST",
              "Effect","StdErr","CHISQ","LOG10P","EXTRA")]
  colnames(dat)<-c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","A1FREQ_CASES","A1FREQ_CONTROLS","INFO","N","N_CASES","N_CONTROLS","TEST",
                   "BETA","SE","CHISQ","LOG10P","EXTRA")
  
  fwrite(dat,paste(path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/eas_iibdgc/",pheno,"/ibd_EAS_SiKJ_meta_",pheno,".regenie",sep=""),
         col.names=T,row.names=F,sep=" ",quote=F,na="NA")
  
  rm(dat,total_N,tmp)
  
}

q("no")








