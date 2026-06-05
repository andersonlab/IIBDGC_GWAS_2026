# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# RUN POST IMPUTATION QC:

############################################################################
# 1.- CREATE RAINFALL PLOTS PER CHROMOSOME FOR GENOTYPED AND IMPUTED DATA: #
############################################################################

# see post_imputation_summary_1.R


###################################################
# 2.- CREATE IMPUTATION QUALITY REPORTS PER STUDY #
###################################################

path_gwas="/path/to/ibdgwas/IIBDGC/"
studies=(australia_omniexome gwas1 gwas2 all_hce pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas chop_old_gwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas niddk_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa)

for i in ${studies[@]}
do
wc -l ${path_gwas}imputed/${i}/qc/list_variants_low_EmpRsq_toexclude
done


for i in ${studies[@]}
do
bash ${path_gwas}scripts/post_imputation_summary_2.sh ${i}
done
# Job <915028..915061> is submitted to queue <normal>.


##### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(dplyr)
library(stringr)
library(viridis)
library(ggExtra)
library(scales)
library(reshape2)

rm(list=ls())

path<-"/path/to/ibdgwas/IIBDGC/"

study<-c("australia_omniexome","gwas1","gwas2","all_hce"
         ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
         ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
         ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
         ,"prism_nfe_gwas","finland_illugwas"
         ,"chop_old_gwas","german_affy6_old_gwas","norway_affy6_old_gwas"
         ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
         ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas"
         ,"swedish_uc_old_gwas","mccauley_gsa","ccfa_gsa","cedars_gsa")

# ancestry<-c("eur","noneur")
ancestry<-c("eur")

for (j in 1:length(ancestry)) {
  
  for (i in 1:length(study)) {
    
    print(study[i])
    
    
    # create summary of number of variants imputed MAF ≥ 0.001 and Rsq ≥ 0.4
    
    maf.file<-paste(path,"imputed/",study[i],"/",ancestry[j],"/plots/",study[i],"_",ancestry[j],"_maf_distribution_variants_rsq_0.4.tsv",sep="")
    
    if(file.exists(maf.file)) {
      
      maf_dist_tmp<-read.table(maf.file,head=T)
      maf_dist_tmp[,'N_allchr']<-rowSums(maf_dist_tmp[,2:(ncol(maf_dist_tmp)-1)])
      
      maf_dist_tmp<-maf_dist_tmp[,c(1,ncol(maf_dist_tmp)),drop=F]
      maf_dist_tmp$study<-study[i]
      
      if(i==1){
        maf_dist<-maf_dist_tmp
      } else {
        maf_dist<-rbind(maf_dist,maf_dist_tmp)
      }
      
    }
    
    
    # Combine all EmpRsq graphs:
    
    empr.file<-paste(path,"imputed/",study[i],"/",ancestry[j],"/plots/",study[i],"_",ancestry[j],"_maf_EmpRsq_mean.tsv",sep="")
    
    if(file.exists(empr.file)) {
      
      emp_tmp<-read.table(empr.file,head=T,sep="\t")
      emp_tmp$study<-study[i]
      emp_tmp<-emp_tmp[,c("MAF_range","study","EmpRsq_mean","EmpRsq_n","MAF")]
      
      if(i==1){
        emp_dist<-emp_tmp
      } else {
        emp_dist<-rbind(emp_dist,emp_tmp)
      }

    }

  }

}


maf_dist$study[which(maf_dist$study=="all_hce")]<-"UK_gwas3"
maf_dist$study[which(maf_dist$study=="gwas1")]<-"UK_gwas1"
maf_dist$study[which(maf_dist$study=="gwas2")]<-"UK_gwas2"

maf_dist$MAF_range<-factor(maf_dist$MAF_range, levels = c("[0,0.0001)","[0.0001,0.001)","[0.001,0.01)","[0.01,0.05)","[0.05,0.1)",
                                                          "[0.1,0.15)","[0.15,0.2)","[0.2,0.25)","[0.25,0.3)","[0.3,0.35)","[0.35,0.4)",
                                                          "[0.4,0.45)","[0.45,0.50]"))
maf_dist<-maf_dist[which(maf_dist$MAF_range %in% c("[0.001,0.01)","[0.01,0.05)","[0.05,0.1)",
                                                   "[0.1,0.15)","[0.15,0.2)","[0.2,0.25)","[0.25,0.3)","[0.3,0.35)","[0.35,0.4)",
                                                   "[0.4,0.45)","[0.45,0.50]")),]


pdf(paste(path,"imputed/plots/maf_distribution_per_study_variants_rsq_0.4_maf_0.001.pdf",sep=""),width = 20, height = 8)
ggplot(data=maf_dist, aes(x=MAF_range,y=N_allchr,fill=factor(study))) +
  geom_bar(position="dodge",stat="identity") + xlab("MAF") + ylab("N variants") +
  ggtitle("Post-imputation MAF distribution\nVariants Rsq >= 0.4 MAF >=0.001")+ scale_y_continuous(labels = comma)
dev.off()
system(paste("cp ",path,"imputed/plots/maf_distribution_per_study_variants_rsq_0.4_maf_0.001.pdf ~/tmp_plots/",sep=""))



emp_dist$study[which(emp_dist$study=="all_hce")]<-"UK_gwas3"
emp_dist$study[which(emp_dist$study=="gwas1")]<-"UK_gwas1"
emp_dist$study[which(emp_dist$study=="gwas2")]<-"UK_gwas2"

emp_dist<-emp_dist[which(emp_dist$MAF_range %in% c("[0.001,0.0025)","[0.0025,0.005)","[0.005,0.0075)","[0.0075,0.01)"
                                                   ,"[0.01,0.025)","[0.025,0.05)","[0.05,0.1)","[0.1,0.2)","[0.2,0.3)"
                                                   ,"[0.3,0.4)","[0.4,1)")),]

emp_dist$MAF_range<-factor(emp_dist$MAF_range, levels = c("[0.001,0.0025)","[0.0025,0.005)","[0.005,0.0075)","[0.0075,0.01)"
                                                          ,"[0.01,0.025)","[0.025,0.05)","[0.05,0.1)","[0.1,0.2)","[0.2,0.3)"
                                                          ,"[0.3,0.4)","[0.4,1)"))

pdf(paste(path,"imputed/plots/maf_EmpRsq_mean_all_studies.pdf",sep=""),width = 11, height = 6)
ggplot(emp_dist, aes(x = -log10(MAF), y = EmpRsq_mean, group = study, color = study)) + ylim(0,1) + ggtitle("Genotyped Variants Mean EmpRsq\nVariants MAF >=0.001") + ylab("Mean Empirical RSQ") +
  geom_point() + geom_line() + expand_limits(x=c(-log10(1E-3),-log10(1E0)) ) + scale_x_continuous(name ="MAF",trans='reverse', 
                                                                                                  breaks= -log10(c(1E-3,1E-2,1E-1,1E0)) ,
                                                                                                  labels=scientific(c(1E-3,1E-2,1E-1,1E0)))
dev.off()
system(paste("cp ",path,"imputed/plots/maf_EmpRsq_mean_all_studies.pdf ~/tmp_plots/",sep=""))


for (j in 1:length(ancestry)) {
  for (i in 1:length(study)) {
    system(paste("cp ",path,"imputed/",study[i],"/",ancestry[j],"/plots/",study[i],"_",ancestry[j],"_maf_EmpRsq.pdf ~/tmp_plots/",sep=""))
  }
}




dat<-matrix(ncol=2,nrow=length(study))
dat<-as.data.frame(dat)
colnames(dat)<-c("study","N_variants")
dat$study<-levels(as.factor(maf_dist$study))

for (i in 1:length(study)) {
  dat$N_variants[i]<-sum(maf_dist$N_allchr[which(maf_dist$study==dat$study[i])])
}

pdf(paste(path,"imputed/plots/n_variants_per_study_variants_rsq_0.4_maf_0.001.pdf",sep=""),width = 20, height = 8)
ggplot(data=dat, aes(x=study,y=N_variants,fill=factor(study))) +
  geom_bar(position="dodge",stat="identity") + xlab("Study") + ylab("N variants") +
  ggtitle("Post-imputation MAF distribution\nVariants Rsq >= 0.4 MAF >=0.001") + scale_y_continuous(labels = comma)
dev.off()
system(paste("cp ",path,"imputed/plots/n_variants_per_study_variants_rsq_0.4_maf_0.001.pdf ~/tmp_plots/",sep=""))



# # syncronize plots
# rsync -a user-server:/path/to/user/tmp_plots/ ~/tmp_plots


###################################################
# 3.- CREATE IMPUTATION QUALITY REPORTS PER STUDY #
###################################################

# AFTER EXCLUDING LOW (<0.25) RSQ VARIANTS, POST IMPUTATION

path_gwas="/path/to/ibdgwas/IIBDGC/"
studies=(australia_omniexome gwas1 gwas2 all_hce pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas chop_old_gwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas niddk_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa)

for i in ${studies[@]}
do
ls -la ${path_gwas}imputed/${i}/qc/list_variants_low_EmpRsq_toexclude
done


for i in ${studies[@]}
do
bash ${path_gwas}scripts/post_imputation_summary_2_update.sh ${i}
done


