# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# get N excluded variant for empRsq >0.5


# unzip all files:

studies=(belgium_inf1_old_gwas belgium_inf2_old_gwas)
path_gwas=/path/to/ibdgwas/IIBDGC/
for i in ${studies[@]}; do
echo ${i}
ls -la ${path_gwas}imputed/${i}/eur/old/
done

for i in ${studies[@]}; do
cd ${path_gwas}imputed/${i}/eur/old/
for chr in {1..22}; do unzip -P PkSw6GlBg5VGuc chr_${chr}.zip;done
unzip -P PkSw6GlBg5VGuc chr_X.zip
done


#######################################################


studies=(swedish_uc_old_gwas niddk_old_gwas)
path_gwas=/path/to/ibdgwas/IIBDGC/
for i in ${studies[@]}; do
echo ${i}
cd ${path_gwas}imputed/${i}/eur/old/
for chr in {1..22}; do unzip -P PkSw6GlBg5VGuc chr_${chr}.zip;done
unzip -P PkSw6GlBg5VGuc chr_X.zip
done

#######################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
studies=(chop_old_gwas german_affy6_old_gwas)

for i in ${studies[@]}; do
echo ${i}
cd ${path_gwas}imputed/${i}/eur/old/
for chr in {1..22}; do unzip -P PkSw6GlBg5VGuc chr_${chr}.zip;done
unzip -P PkSw6GlBg5VGuc chr_X.zip
done


#######################################################

studies=(norway_affy6_old_gwas cedars_370k_old_gwas)
path_gwas=/path/to/ibdgwas/IIBDGC/
for i in ${studies[@]}; do
echo ${i}
ls -la ${path_gwas}imputed/${i}/eur/old/
done

for i in ${studies[@]}; do
cd ${path_gwas}imputed/${i}/eur/old/
for chr in {1..22}; do unzip -P PkSw6GlBg5VGuc chr_${chr}.zip;done
unzip -P PkSw6GlBg5VGuc chr_X.zip
done


#######################################################

studies=( cedars_610k_old_gwas cedars_omni_old_gwas)
path_gwas=/path/to/ibdgwas/IIBDGC/
for i in ${studies[@]}; do
echo ${i}
ls -la ${path_gwas}imputed/${i}/eur/old/
done

for i in ${studies[@]}; do
cd ${path_gwas}imputed/${i}/eur/old/
  for chr in {1..22}; do unzip -P PkSw6GlBg5VGuc chr_${chr}.zip;done
unzip -P PkSw6GlBg5VGuc chr_X.zip
done

#######################################################

studies=(australia_omniexome gwas1 gwas2)
path_gwas=/path/to/ibdgwas/IIBDGC/
  for i in ${studies[@]}; do
echo ${i}
ls -la ${path_gwas}imputed/${i}/eur/old/
  done

for i in ${studies[@]}; do
cd ${path_gwas}imputed/${i}/eur/old/
  for chr in {1..22}; do unzip -P PkSw6GlBg5VGuc chr_${chr}.zip;done
unzip -P PkSw6GlBg5VGuc chr_X.zip
done

#######################################################

studies=(all_hce)
path_gwas=/path/to/ibdgwas/IIBDGC/
for i in ${studies[@]}; do
echo ${i}
ls -la ${path_gwas}imputed/${i}/eur/
done

for i in ${studies[@]}; do
cd ${path_gwas}imputed/${i}/eur/
for chr in {1..22}; do unzip -P PkSw6GlBg5VGuc chr_${chr}.zip;done
unzip -P PkSw6GlBg5VGuc chr_X.zip
done

#######################################################

studies=(pittsburgh_gsa italy_gsa)
path_gwas=/path/to/ibdgwas/IIBDGC/
for i in ${studies[@]}; do
echo ${i}
ls -la ${path_gwas}imputed/${i}/eur/old/
done

for i in ${studies[@]}; do
cd ${path_gwas}imputed/${i}/eur/old/
for chr in {1..22}; do unzip -P PkSw6GlBg5VGuc chr_${chr}.zip;done
unzip -P PkSw6GlBg5VGuc chr_X.zip
done

#######################################################

studies=(spain_gsa)
path_gwas=/path/to/ibdgwas/IIBDGC/
for i in ${studies[@]}; do
echo ${i}
ls -la ${path_gwas}imputed/${i}/eur/old/
done

for i in ${studies[@]}; do
cd ${path_gwas}imputed/${i}/eur/old/
for chr in {1..22}; do unzip -P PkSw6GlBg5VGuc chr_${chr}.zip;done
unzip -P PkSw6GlBg5VGuc chr_X.zip
done

#######################################################

studies=(niddk_feinstein_gsa)
path_gwas=/path/to/ibdgwas/IIBDGC/
for i in ${studies[@]}; do
echo ${i}
ls -la ${path_gwas}imputed/${i}/eur/
done

for i in ${studies[@]}; do
cd ${path_gwas}imputed/${i}/eur/
for chr in {1..22}; do unzip -P PkSw6GlBg5VGuc chr_${chr}.zip;done
unzip -P PkSw6GlBg5VGuc chr_X.zip
done

#######################################################


studies=(basque_gsa lithuania_gsa)
path_gwas=/path/to/ibdgwas/IIBDGC/
  for i in ${studies[@]}; do
echo ${i}
ls -la ${path_gwas}imputed/${i}/eur/
  done

for i in ${studies[@]}; do
cd ${path_gwas}imputed/${i}/eur/
for chr in {1..22}; do unzip -P PkSw6GlBg5VGuc chr_${chr}.zip;done
unzip -P PkSw6GlBg5VGuc chr_X.zip
done

#######################################################

#######################################################


studies=(prism_nfe_gwas finland_illugwas)
path_gwas=/path/to/ibdgwas/IIBDGC/
  for i in ${studies[@]}; do
echo ${i}
ls -la ${path_gwas}imputed/${i}/eur/
  done

for i in ${studies[@]}; do
cd ${path_gwas}imputed/${i}/eur/
  for chr in {1..22}; do unzip -P PkSw6GlBg5VGuc chr_${chr}.zip;done
unzip -P PkSw6GlBg5VGuc chr_X.zip
done

#######################################################

studies=(netherlands_gsa)
path_gwas=/path/to/ibdgwas/IIBDGC/
  for i in ${studies[@]}; do
echo ${i}
ls -la ${path_gwas}imputed/${i}/eur/
  done

for i in ${studies[@]}; do
cd ${path_gwas}imputed/${i}/eur/
  for chr in {1..22}; do unzip -P PkSw6GlBg5VGuc chr_${chr}.zip;done
unzip -P PkSw6GlBg5VGuc chr_X.zip
done

#######################################################

studies=(slovenia_gsa)
path_gwas=/path/to/ibdgwas/IIBDGC/
  for i in ${studies[@]}; do
echo ${i}
ls -la ${path_gwas}imputed/${i}/eur/
  done

for i in ${studies[@]}; do
cd ${path_gwas}imputed/${i}/eur/
  for chr in {1..22}; do unzip -P PkSw6GlBg5VGuc chr_${chr}.zip;done
unzip -P PkSw6GlBg5VGuc chr_X.zip
done

#######################################################

studies=(sweden_gsa)
path_gwas=/path/to/ibdgwas/IIBDGC/
  for i in ${studies[@]}; do
echo ${i}
ls -la ${path_gwas}imputed/${i}/eur/
  done

for i in ${studies[@]}; do
cd ${path_gwas}imputed/${i}/eur/
  for chr in {1..22}; do unzip -P PkSw6GlBg5VGuc chr_${chr}.zip;done
unzip -P PkSw6GlBg5VGuc chr_X.zip
done



studies=(australia_omniexome gwas1 gwas2 all_hce pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas chop_old_gwas german_affy6_old_gwas norway_affy6_old_gwas niddk_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa)


##### 


MEM=6000
bsub -Is -m "modern_hardware" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1


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

# no german_illu
study<-c("australia_omniexome","gwas1","gwas2"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"chop_old_gwas","german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas"
           ,"swedish_uc_old_gwas","mccauley_gsa","ccfa_gsa","cedars_gsa","all_hce")

length(study)
# [1] 34


data<-matrix(nrow=length(study),ncol=17)
data<-as.data.frame(data)
colnames(data)<-c("study","N_genotyped"
                  ,"N_excluded_EmpRsq_0.25","N_excluded_EmpRsq_0.25_maf_below_0.01","N_excluded_EmpRsq_0.25_maf_below_0.05","Mean_maf_excluded_EmpRsq_0.25","Median_maf_excluded_EmpRsq_0.25","Min_maf_excluded_EmpRsq_0.25","Max_maf_excluded_EmpRsq_0.25"
                  ,"N_excluded_EmpRsq_0.50","N_excluded_EmpRsq_0.50_maf_below_0.01","N_excluded_EmpRsq_0.50_maf_below_0.05","Mean_maf_excluded_EmpRsq_0.50","Median_maf_excluded_EmpRsq_0.50","Min_maf_excluded_EmpRsq_0.50","Max_maf_excluded_EmpRsq_0.50"
                  ,"N_flip")

for (j in c(31:length(study))) {
  
  for (i in 1:23) {
    
    if(dir.exists(paste(path,"imputed/",study[j],"/eur/old/",sep=""))) {
      if(i==23) {
        tmp<-fread(paste(path,"imputed/",study[j],"/eur/old/chrX.info.gz",sep=""),head=T)
      } else {
        tmp<-fread(paste(path,"imputed/",study[j],"/eur/old/chr",i,".info.gz",sep=""),head=T)
      }
    } else {
      if(i==23) {
        tmp<-fread(paste(path,"imputed/",study[j],"/eur/chrX.info.gz",sep=""),head=T)
      } else {
        tmp<-fread(paste(path,"imputed/",study[j],"/eur/chr",i,".info.gz",sep=""),head=T)
      }
    }

    
    print(i)
    tmp<-tmp[which(tmp$Genotyped=="Genotyped"),]
    tmp<-as.data.frame(tmp)
    
    tmp$chr<-i

    if(i==1) {
      df<-tmp
    }else{
      df<-rbind(tmp,df)
    }
    
    rm(tmp,tmp1)
    
  }
  
  df$EmpRsq<-as.numeric(df$EmpRsq)
  df$EmpR<-as.numeric(df$EmpR)
  
  
  data$study[j]<-study[j]
  data$N_genotyped[j]<-nrow(df)
  
  data$N_excluded_EmpRsq_0.25[j]<-nrow(df[which(df$EmpRsq<=0.25),])
  data$N_excluded_EmpRsq_0.25_maf_below_0.01[j]<-nrow(df[which(df$EmpRsq<=0.25 & df$MAF<0.01),])
  data$N_excluded_EmpRsq_0.25_maf_below_0.05[j]<-nrow(df[which(df$EmpRsq<=0.25 & df$MAF<0.05),])
  data$Mean_maf_excluded_EmpRsq_0.25[j]<-mean(df$MAF[which(df$EmpRsq<=0.25)])
  data$Median_maf_excluded_EmpRsq_0.25[j]<-median(df$MAF[which(df$EmpRsq<=0.25)])
  data$Min_maf_excluded_EmpRsq_0.25[j]<-min(df$MAF[which(df$EmpRsq<=0.25)])
  data$Max_maf_excluded_EmpRsq_0.25[j]<-max(df$MAF[which(df$EmpRsq<=0.25)])
  
  data$N_excluded_EmpRsq_0.50[j]<-nrow(df[which(df$EmpRsq<=0.50),])
  data$N_excluded_EmpRsq_0.50_maf_below_0.01[j]<-nrow(df[which(df$EmpRsq<=0.50 & df$MAF<0.01),])
  data$N_excluded_EmpRsq_0.50_maf_below_0.05[j]<-nrow(df[which(df$EmpRsq<=0.50 & df$MAF<0.05),])
  data$Mean_maf_excluded_EmpRsq_0.50[j]<-mean(df$MAF[which(df$EmpRsq<=0.50)])
  data$Median_maf_excluded_EmpRsq_0.50[j]<-median(df$MAF[which(df$EmpRsq<=0.50)])
  data$Min_maf_excluded_EmpRsq_0.50[j]<-min(df$MAF[which(df$EmpRsq<=0.50)])
  data$Max_maf_excluded_EmpRsq_0.50[j]<-max(df$MAF[which(df$EmpRsq<=0.50)])
  
  data$N_flip[j]<-nrow(df[which(df$EmpR<= -0.5),])
  
  print(study[j])
  print(data[j,])
  
  # variants neg EmpR to flip:
  df_flip<-df[which(df$EmpR<= -0.5),]
  
  df_flip$SNP<-gsub(":","_",df_flip$SNP)
  df_flip$SNP<-sub("_",":",df_flip$SNP)
  df_flip$SNP<-sub("chrX","chr23",df_flip$SNP)
  df_flip$SNP<-sub("chr","",df_flip$SNP)
  
  write.table(df_flip[,"SNP",drop=F],paste(path,"imputed/",study[j],"/qc/list_variants_negative_EmpR_toflip",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  
  
  # variants low EmpRsq to remove:
  df_remove<-df[which(df$EmpRsq<=0.25),]
  
  df_remove$SNP<-gsub(":","_",df_remove$SNP)
  df_remove$SNP<-sub("_",":",df_remove$SNP)
  df_remove$SNP<-sub("chrX","chr23",df_remove$SNP)
  df_remove$SNP<-sub("chr","",df_remove$SNP)
  df_remove<-df_remove[which(!df_remove$SNP %in% df_flip$SNP),]
  
  write.table(df_remove[,"SNP",drop=F],paste(path,"imputed/",study[j],"/qc/list_variants_0.25_EmpRsq_toexclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  
  # variants low EmpRsq to remove:
  df_remove<-df[which(df$EmpRsq<=0.50),]
  
  df_remove$SNP<-gsub(":","_",df_remove$SNP)
  df_remove$SNP<-sub("_",":",df_remove$SNP)
  df_remove$SNP<-sub("chrX","chr23",df_remove$SNP)
  df_remove$SNP<-sub("chr","",df_remove$SNP)
  df_remove<-df_remove[which(!df_remove$SNP %in% df_flip$SNP),]
  
  write.table(df_remove[,"SNP",drop=F],paste(path,"imputed/",study[j],"/qc/list_variants_0.5_EmpRsq_toexclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  
  
}


data$percentage_excluded_0.25<-NA
data$percentage_excluded_0.25<-100*(data$N_excluded_EmpRsq_0.25/data$N_genotyped)

data$percentage_excluded_0.50<-NA
data$percentage_excluded_0.50<-100*(data$N_excluded_EmpRsq_0.50/data$N_genotyped)


write.table(data,paste(path,"imputed/plots/summary_variants_excluded_by_EmpRsq_percohort.tsv",sep=""),col.names=T,row.names=F,sep="\t",quote=F)


q("no")




cp /path/to/ibdgwas/IIBDGC/imputed/plots/summary_variants_excluded_by_EmpRsq_percohort.tsv ~/tmp_plots/






