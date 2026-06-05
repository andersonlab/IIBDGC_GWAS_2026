# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# 
# # singularity exec iibdgc_postprocess_10_singularity.sif
# MEM="8000"
# bsub -J"ph" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 16 \
# -e ${path_gwas}post_imputation/analysis/metaanalysis/log/stderr_def_associated_ph_index \
# -o ${path_gwas}post_imputation/analysis/metaanalysis/log/stdout_def_associated_ph_index \
# "R CMD BATCH ~/git/IIBDGC_GWAS/scripts/other/define_associated_phenotype_index_variants_mc16.R"

# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G humgen-priority -q yesterday -n 4 R 

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(nnet)
library(multinomRob)
library(plyr)
library(doMC)

rm(list=ls())

registerDoMC(16)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

args = commandArgs(trailingOnly=TRUE)

chr<-args[1]

# keep only samples two largest arrays, including small ones and adjust by array does not allow to converge the model:
array<-c("humancoreexome","gsa")

# "chr10:105644358:TA:T_T" "chr10:125821742:T:G_G"

rm(dat)
for (chr in chr) {

  for (i in 1:length(array)) {

    # file_raw<-paste(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/",array[i],"_chr",chr,"_subset_new_and_old_index_gwas_variants_plus_credible_sets.raw",sep="")
    file_raw<-paste(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/",array[i],"_chr",chr,"_list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_with_fm_signals_with_multiancestry_signals.raw",sep="")
    
    if(file.exists(file_raw)) {

      tmp<-fread(file_raw,head=T,check.names=F)
      tmp$array<-array[i]
      
      if(i==1){
        tmp1<-tmp
      } else {

        # ids<-intersect(colnames(tmp1),colnames(tmp))
        # tmp<-tmp[,ids]
        # tmp1<-tmp1[,ids]

        tmp1<-rbind(tmp1,tmp,fill=T)

      }
      rm(tmp)
    }
  }

  if(file.exists(file_raw)) {
    if(!exists("dat")){
      dat<-tmp1
    }else{
      dat<-merge(dat,tmp1[,c(2,7:(ncol(tmp1)-1))],by="IID")
    }
    rm(tmp1)
  } else {
    next
  }
  
}

dim(dat)
# [1] 78561   906
dat<-as.data.frame(dat)

ids1<-colnames(dat)[7:ncol(dat)]

ids1<-ids1[which(!ids1 %in% c("chr1:66934185:C:T_T","array"))]
length(ids1)
# [1] 898

all_ids<-c("IID","FID","PAT","MAT","SEX","PHENOTYPE","array",as.character(ids1))
dat<-dat[,all_ids]
dim(dat)
# [1] 78561   905


## add phenotype data:

for (i in 1:length(array)) {
  
  print(array[i])
  
  tmp1<-read.table(paste(path_gwas,"post_imputation/2022/",array[i],"/phenotype_data/",array[i],"_all_studies_merged_eur_all_covariate_sex_PCs",sep=""),head=T)
  tmp2<-read.table(paste(path_gwas,"post_imputation/2022/",array[i],"/phenotype_data/",array[i],"_all_studies_merged_eur_all_phenotype_ibd",sep=""),head=T)
  tmp3<-read.table(paste(path_gwas,"post_imputation/2022/",array[i],"/phenotype_data/",array[i],"_all_studies_merged_eur_all_phenotype_cd",sep=""),head=T)
  tmp4<-read.table(paste(path_gwas,"post_imputation/2022/",array[i],"/phenotype_data/",array[i],"_all_studies_merged_eur_all_phenotype_uc",sep=""),head=T)
  
  tmp<-merge(tmp1[,c(1:6)],tmp2[,c(2:3)],by="IID",all.x=T)
  tmp<-merge(tmp,tmp3[,c(2:3)],by="IID",all.x=T)
  tmp<-merge(tmp,tmp4[,c(2:3)],by="IID",all.x=T)
  
  if(i==1){
    pheno<-tmp
  } else {
    pheno<-rbind(pheno,tmp)
  }
}

dat<-merge(dat[,c(1,7:ncol(dat))],pheno[,c(1,3:ncol(pheno))],by="IID",all.y=T)

# exclude relatives for this test:
rel<-read.table(paste(path_gwas,"pre_imputation/QC/pca_1000gp/list_eur_related_3rd_degree_samples",sep=""),head=T)
dat<-dat[which(!dat$IID %in% rel$IID),]
dim(dat)
# [1] 58913    9788

# exclude jewish ancestry samples:
eur_nonjw<-read.table(paste(path_gwas,"pre_imputation/QC/pca_1000gp/list_eur_nonjewish_ancestry_samples_noDuplicates",sep=""),head=F)
dat<-dat[which(dat$IID %in% eur_nonjw$V2),]
dim(dat)
# [1] 53878    9788

# get new pheno variable

dat$status_ibd<-NA
dat$status_ibd[which(dat$ibd==0)]<-0
dat$status_ibd[which(dat$cd==1)]<-1
dat$status_ibd[which(dat$uc==1)]<-2

dat<-dat[!is.na(dat$status_ibd),]
table(dat$status_ibd,useNA="ifany")

table(dat$status_ibd,dat$cd,useNA="ifany")
table(dat$status_ibd,dat$uc,useNA="ifany")

# get pheno in format for multinomrob
dat$ctr<-NA
dat$cd<-NA
dat$uc<-NA

dat$ctr <- ifelse(dat$status_ibd=="0",1, 0)
dat$cd <- ifelse(dat$status_ibd=="1", 1, 0)
dat$uc <- ifelse(dat$status_ibd=="2", 1, 0)

table(dat$status_ibd,dat$ctr,useNA="ifany")
table(dat$status_ibd,dat$cd,useNA="ifany")
table(dat$status_ibd,dat$uc,useNA="ifany")

# get list of variants
ids<-colnames(dat)[grep("chr",colnames(dat))]
ids<-c("IID","array","status_ibd","cd","uc","ctr","PC1","PC2","PC3",as.character(ids))
dat<-dat[,ids]

table(dat$status_ibd,dat$ctr)
table(dat$status_ibd,dat$cd)
table(dat$status_ibd,dat$uc)


mod<-foreach (n=10:ncol(dat),.combine=rbind) %dopar% {
  
  # print(n)
  tmp<-dat[,c(2,4:9,n)]
  print(colnames(tmp)[ncol(tmp)])

  colnames(tmp)[8]<-"var"

  # to handle SNPs with data only from one array
  tmp <- na.omit(tmp)
  
  tmp$array<-as.factor(tmp$array)
  if (length(levels(tmp$array))==2) {

    # model with no variant included
    model.0<-multinomRob(list(cd ~ 1 + array + PC1 + PC2 + PC3, uc ~ 1 + array + PC1 + PC2 + PC3, ctr ~ 0),data=tmp, MLEonly=TRUE)
    
    # Crohns disease specific model: UC = 0
    model.1<-multinomRob(list(cd ~ var + array + PC1 + PC2 + PC3, uc ~ 1 + array + PC1 + PC2 + PC3, ctr ~ 0),data=tmp, MLEonly=TRUE)
    
    # Ulcerative colitis disease specific model: CD = 0
    model.2<-multinomRob(list(cd ~ 1 + array + PC1 + PC2 + PC3, uc ~ var + array + PC1 + PC2 + PC3, ctr ~ 0),data=tmp, MLEonly=TRUE)
    
    # IBD unsaturated model, equal constraints CD = UC = ibd
    model.3<-multinomRob(list(cd ~ var + array + PC1 + PC2 + PC3, uc ~ var + array + PC1 + PC2 + PC3, ctr ~ 0),data=tmp,equality=list(list(cd~var + 0,uc~var + 0)), MLEonly=TRUE)
    
    # IBD saturated model
    model.4<-multinomRob(list(cd ~ var + array + PC1 + PC2 + PC3, uc ~ var + array + PC1 + PC2 + PC3, ctr ~ 0),data=tmp, MLEonly=TRUE)

  } else if (length(levels(tmp$array))==1) {

    # model with no variant included
    model.0<-multinomRob(list(cd ~ 1 + PC1 + PC2 + PC3, uc ~ 1 + PC1 + PC2 + PC3, ctr ~ 0),data=tmp, MLEonly=TRUE)
    
    # Crohns disease specific model: UC = 0
    model.1<-multinomRob(list(cd ~ var + PC1 + PC2 + PC3, uc ~ 1 + PC1 + PC2 + PC3, ctr ~ 0),data=tmp, MLEonly=TRUE)
    
    # Ulcerative colitis disease specific model: CD = 0
    model.2<-multinomRob(list(cd ~ 1 + PC1 + PC2 + PC3, uc ~ var + PC1 + PC2 + PC3, ctr ~ 0),data=tmp, MLEonly=TRUE)
    
    # IBD unsaturated model, equal constraints CD = UC = ibd
    model.3<-multinomRob(list(cd ~ var + PC1 + PC2 + PC3, uc ~ var + PC1 + PC2 + PC3, ctr ~ 0),data=tmp,equality=list(list(cd~var + 0,uc~var + 0)), MLEonly=TRUE)
    
    # IBD saturated model
    model.4<-multinomRob(list(cd ~ var + PC1 + PC2 + PC3, uc ~ var + PC1 + PC2 + PC3, ctr ~ 0),data=tmp, MLEonly=TRUE)

  }
    
  cbind(as.character(colnames(dat)[n]),model.0$deviance,model.1$deviance,model.2$deviance,model.3$deviance,model.4$deviance)
  
}

write.table(mod,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/chr",chr,"_tier1_EUR_array_data_fine_mapping_tier_2_unsupervised_conditional_old_index_subphenotype_models_adjusted_by_array_mc16_with_PCs.tsv.gz",sep=""),
            col.names=T,row.names=F,sep="\t",quote=F)
            
colnames(mod)<-c("var_id","deviance_model_0","deviance_model_1","deviance_model_2","deviance_model_3","deviance_model_4")

mod<-as.data.frame(mod)
mod$deviance_model_0<-as.numeric(as.character(mod$deviance_model_0))
mod$deviance_model_1<-as.numeric(as.character(mod$deviance_model_1))
mod$deviance_model_2<-as.numeric(as.character(mod$deviance_model_2))
mod$deviance_model_3<-as.numeric(as.character(mod$deviance_model_3))
mod$deviance_model_4<-as.numeric(as.character(mod$deviance_model_4))

mod$pvalue_m1_m4<-pchisq(mod$deviance_model_1-mod$deviance_model_4,df=1,lower.tail=F)
mod$pvalue_m2_m4<-pchisq(mod$deviance_model_2-mod$deviance_model_4,df=1,lower.tail=F)
mod$pvalue_m3_m4<-pchisq(mod$deviance_model_3-mod$deviance_model_4,df=1,lower.tail=F)

mod$pvalue_m1_m0<-pchisq(mod$deviance_model_0-mod$deviance_model_1,df=1,lower.tail=F)
mod$pvalue_m2_m0<-pchisq(mod$deviance_model_0-mod$deviance_model_2,df=1,lower.tail=F)
mod$pvalue_m3_m0<-pchisq(mod$deviance_model_0-mod$deviance_model_3,df=1,lower.tail=F)

mod$model_1<-NA

for (i in 1:nrow(mod)) {
  
  if(mod$pvalue_m1_m4[i]<0.05 & mod$pvalue_m2_m4[i]<0.05 & mod$pvalue_m3_m4[i]<0.05) {
    mod$model_1[i]<-"IBD_saturated"
  } else {
    mod$model_1[i]<-colnames(mod[i,which(mod[i,]==min(mod[i,3:5])),drop=F])
  }
}


table(mod$model_1)

mod$phenotype<-NA
mod$phenotype[which(mod$model_1=="IBD_saturated")]<-"IBD_saturated"
mod$phenotype[which(mod$model_1=="deviance_model_1")]<-"CD"
mod$phenotype[which(mod$model_1=="deviance_model_2")]<-"UC"
mod$phenotype[which(mod$model_1=="deviance_model_3")]<-"IBD_unsaturated"

write.table(mod,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/chr",chr,"_tier1_EUR_array_data_fine_mapping_tier_2_unsupervised_conditional_old_index_subphenotype_models_adjusted_by_array_by_pcs_mc16_with_PCs.tsv.gz",sep=""),
            col.names=T,row.names=F,sep="\t",quote=F)

q("no")