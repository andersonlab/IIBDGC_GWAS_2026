# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
############################################################################################################################
############################################################################################################################

# singularity exec iibdgc_postprocess_10_singularity.sif

# 1.- create a new definition of region (updated wiht Liu etc)

MEM=5000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)

rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

pheno<-c("cd","ibd","uc")

rm(dat,all_dat)
for (ph in pheno) {

       print(ph)

       for (chr in c(1:23)) {

              print(chr)

              file_tmp<-paste0(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/",chr,"_",ph,"_independent_index_variants_eur_tier_2_ld_allarrays_no_conditioning.jma.cojo")

              if (chr==23) {
                file_tmp<-paste0(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/X_",ph,"_independent_index_variants_eur_tier_2_ld_allarrays_no_conditioning.jma.cojo")
              }
              
              if(file.exists(file_tmp)) {
                tmp<-fread(file_tmp,head=T)

                tmp<-tmp[,c("SNP","bJ","bJ_se","pJ")]
                colnames(tmp)<-c("MarkerName","BETA_cond","SE_cond","Pvalue_cond")
                colnames(tmp)[2:ncol(tmp)]<-paste(colnames(tmp)[2:ncol(tmp)],ph,sep="_")

                if (!exists("dat")) {
                     dat<-tmp
                } else {
                    dat<-rbind(dat,tmp)
                }
                rm(tmp)
              }

              rm(file_tmp)
       }

       if(!exists("all_dat")) {
              all_dat<-dat
       } else {
              all_dat<-merge(all_dat,dat,by="MarkerName",all=T)
       }
       rm(dat)

}

dat<-all_dat
fwrite(dat,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")

dat$chr<-gsub(":.*","",dat$MarkerName)
dat$chr<-gsub("chr*","",dat$chr)

chr<-as.numeric(unlist(labels(table(dat$chr))))
chr<-chr[which(!is.na(chr))]

dat<-as.data.frame(dat)

dim(dat)
# [1] 1080   11


##############################################################
# results from EUR tier 1

for (i in 1:length(pheno)) {
  
  print(pheno[i])
  
  for (j in 1:length(chr)) {
    
    tmp1<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr[j],"_",pheno[i],"_meta_eur_tier_1_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff_with_rsid_dbsnp154.txt.gz",sep=""),head=T)
    tmp1<-tmp1[which(tmp1$MarkerName %in% dat$MarkerName),]

    if(chr[j]=="X") {
      tmp1$dbsnp154<-NA
    }
    
    if(j==1) {
      tmp<-tmp1
    } else {
      tmp<-rbind(tmp,tmp1)
    }
    
    rm(tmp1)
  
  }
  
  if(i==1) {
    tmp<-tmp[,c("MarkerName","dbsnp154","Position_b38","A1","A2","BETA","SE","P-value","Direction_ed","HetISq",
                "HetChiSq","HetDf","HetPVal","rate_Neff")]
    colnames(tmp)[2:4]<-paste(colnames(tmp)[2:4],"eur_tier_1",sep="_")
    colnames(tmp)[5:ncol(tmp)]<-paste(colnames(tmp)[5:ncol(tmp)],pheno[i],"eur_tier_1",sep="_")
    all<-tmp
  } else {
    tmp<-tmp[,c("MarkerName","BETA","SE","P-value","Direction_ed","HetISq",
                "HetChiSq","HetDf","HetPVal","rate_Neff")]
    colnames(tmp)[2:ncol(tmp)]<-paste(colnames(tmp)[2:ncol(tmp)],pheno[i],"eur_tier_1",sep="_")
    all<-merge(all,tmp,by="MarkerName",sort=F,all=T)
  }
  rm(tmp)
}

dat<-as.data.frame(dat)
all<-as.data.frame(all)

dat<-merge(dat,all,by="MarkerName",all.x=T)


fwrite(dat,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")

rm(all)


##############################################################
# results from EUR tier 2

for (i in 1:length(pheno)) {
  
  print(pheno[i])
  
  for (j in 1:length(chr)) {
    
    tmp1<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr[j],"_",pheno[i],"_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz",sep=""),head=T)
    tmp1<-tmp1[which(tmp1$MarkerName %in% dat$MarkerName),]
    
    if(j==1) {
      tmp<-tmp1
    } else {
      tmp<-rbind(tmp,tmp1)
    }
    
    rm(tmp1)
    
  }
  
  if(i==1) {
    tmp<-tmp[,c("MarkerName","Position_b38","A1","A2","BETA","SE","P-value","Direction_ed","HetISq",
                "HetChiSq","HetDf","HetPVal","rate_Neff")]
    colnames(tmp)[2:4]<-paste(colnames(tmp)[2:4],"eur_tier_2",sep="_")
    colnames(tmp)[5:ncol(tmp)]<-paste(colnames(tmp)[5:ncol(tmp)],pheno[i],"eur_tier_2",sep="_")
    all<-tmp
  } else {
    tmp<-tmp[,c("MarkerName","BETA","SE","P-value","Direction_ed","HetISq",
                "HetChiSq","HetDf","HetPVal","rate_Neff")]
    colnames(tmp)[2:ncol(tmp)]<-paste(colnames(tmp)[2:ncol(tmp)],pheno[i],"eur_tier_2",sep="_")
    all<-merge(all,tmp,by="MarkerName",sort=F,all=T)
  }
  rm(tmp)
}

dim(dat)
dim(all)


all<-as.data.frame(all)

dat<-merge(dat,all,by="MarkerName",all.x=T)

      
fwrite(dat,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")
rm(all)


##############################################################
### add results from meta SAS

for (i in 1:length(pheno)) {
  
  print(pheno[i])
  
  for (j in 1:length(chr)) {
    
    tmp1<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr[j],"_",pheno[i],"_meta_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz",sep=""),head=T)
    tmp1<-tmp1[which(tmp1$MarkerName %in% dat$MarkerName),]
    
    if(j==1) {
      tmp<-tmp1
    } else {
      tmp<-rbind(tmp,tmp1)
    }
    rm(tmp1)
  }
  
  if(i==1) {
    tmp<-tmp[,c("MarkerName","Position_b38","A1","A2","BETA","SE","P-value","Direction_ed","HetISq",
                "HetChiSq","HetDf","HetPVal","rate_Neff")]
    colnames(tmp)[2:4]<-paste(colnames(tmp)[2:4],"sas",sep="_")
    colnames(tmp)[5:ncol(tmp)]<-paste(colnames(tmp)[5:ncol(tmp)],pheno[i],"sas",sep="_")
    all<-tmp
  } else {
    tmp<-tmp[,c("MarkerName","BETA","SE","P-value","Direction_ed","HetISq",
                "HetChiSq","HetDf","HetPVal","rate_Neff")]
    colnames(tmp)[2:ncol(tmp)]<-paste(colnames(tmp)[2:ncol(tmp)],pheno[i],"sas",sep="_")
    all<-merge(all,tmp,by="MarkerName",sort=F,all=T)
  }
  rm(tmp)
}

all<-as.data.frame(all)
dat<-merge(dat,all,by="MarkerName",all.x=T)

fwrite(dat,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")
rm(all)


##############################################################
### add results from meta IIBDGC_EUR_tier2 with EAS (liu) AND SAS


for (i in 1:length(pheno)) {
  
  print(pheno[i])
  
  for (j in 1:length(chr)) {
    
    tmp1<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr[j],"_",pheno[i],"_meta_eur_eas_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz",sep=""),head=T)
    tmp1<-tmp1[which(tmp1$MarkerName %in% dat$MarkerName),]
    
    if(j==1) {
      tmp<-tmp1
    } else {
      tmp<-rbind(tmp,tmp1)
    }
    
    rm(tmp1)
    
  }
  
  if(i==1) {
    tmp<-tmp[,c("MarkerName","Position_b38","A1","A2","BETA","SE","P-value","Direction_ed","HetISq",
                "HetChiSq","HetDf","HetPVal","rate_Neff")]
    colnames(tmp)[2:4]<-paste(colnames(tmp)[2:4],"eur_tier_2_eas_sas",sep="_")
    colnames(tmp)[5:ncol(tmp)]<-paste(colnames(tmp)[5:ncol(tmp)],pheno[i],"eur_tier_2_eas_sas",sep="_")
    all<-tmp
  } else {
    tmp<-tmp[,c("MarkerName","BETA","SE","P-value","Direction_ed","HetISq",
                "HetChiSq","HetDf","HetPVal","rate_Neff")]
    colnames(tmp)[2:ncol(tmp)]<-paste(colnames(tmp)[2:ncol(tmp)],pheno[i],"eur_tier_2_eas_sas",sep="_")
    all<-merge(all,tmp,by="MarkerName",sort=F,all=T)
  }
  rm(tmp)
}

all<-as.data.frame(all)

dat<-merge(dat,all,by="MarkerName",all.x=T)
dim(dat)
# [1] 1080  132

fwrite(dat,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")

rm(all)

q("no")


# MEM=35000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"
pheno<-c("cd","ibd","uc")

dat<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency.tsv.gz",sep=""))
dat<-as.data.frame(dat)

# add results from MR-mega

for (i in 1:length(pheno)) {
  
  print(pheno[i])
  
  tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/mrmega/output_files/",pheno[i],"/mrmega_eur_eas_sas_tier_2_",pheno[i],"_allchr_pc1.result",sep=""))
  tmp<-tmp[which(tmp$MarkerName %in% dat$MarkerName),]
  
  if(i==1) {
    tmp<-tmp[,c("MarkerName","Position","NEA","EA","P-value_association","chisq_ancestry_het","P-value_ancestry_het","P-value_residual_het")]
    colnames(tmp)[2:4]<-paste(colnames(tmp)[2:4],"eur_tier_2_eas_sas_mrmega",sep="_")
    colnames(tmp)[5:ncol(tmp)]<-paste(colnames(tmp)[5:ncol(tmp)],pheno[i],"eur_tier_2_eas_sas_mrmega",sep="_")
    all<-tmp
  } else {
    tmp<-tmp[,c("MarkerName","Position","NEA","EA","P-value_association","chisq_ancestry_het","P-value_ancestry_het","P-value_residual_het")]
    colnames(tmp)[2:ncol(tmp)]<-paste(colnames(tmp)[2:ncol(tmp)],pheno[i],"eur_tier_2_eas_sas_mrmega",sep="_")
    all<-merge(all,tmp,by="MarkerName",sort=F,all=T)
  }
  rm(tmp)
}

all<-as.data.frame(all)

dat<-merge(dat,all,by="MarkerName",all.x=T)

fwrite(dat,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")

("no")

##############################################################################

# singularity exec iibdgc_postprocess_10_singularity.sif


MEM=25000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)

rm(list=ls())


path_gwas<-"/path/to/ibdgwas/IIBDGC/"
pheno<-c("cd","ibd","uc")

dat<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency.tsv.gz",sep=""))
dat<-as.data.frame(dat)

###############################
# add results from MR-mega

for (i in 1:length(pheno)) {
  
  print(pheno[i])
  
  tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/mrmega/output_files/",pheno[i],"/mrmega_eur_eas_sas_tier_2_",pheno[i],"_allchr_pc1.result",sep=""))
  tmp<-tmp[which(tmp$MarkerName %in% dat$MarkerName),]
  
  if(i==1) {
    tmp<-tmp[,c("MarkerName","Position","NEA","EA","P-value_association","chisq_ancestry_het","P-value_ancestry_het","P-value_residual_het")]
    colnames(tmp)[2:4]<-paste(colnames(tmp)[2:4],"eur_tier_2_eas_sas_mrmega",sep="_")
    colnames(tmp)[5:ncol(tmp)]<-paste(colnames(tmp)[5:ncol(tmp)],pheno[i],"eur_tier_2_eas_sas_mrmega",sep="_")
    all<-tmp
  } else {
    tmp<-tmp[,c("MarkerName","Position","NEA","EA","P-value_association","chisq_ancestry_het","P-value_ancestry_het","P-value_residual_het")]
    colnames(tmp)[2:ncol(tmp)]<-paste(colnames(tmp)[2:ncol(tmp)],pheno[i],"eur_tier_2_eas_sas_mrmega",sep="_")
    all<-merge(all,tmp,by="MarkerName",sort=F)
  }
  rm(tmp)
}

all<-as.data.frame(all)

dat<-merge(dat,all,by="MarkerName",all.x=T)
dim(dat)
# [1] [1] 1050  153


fwrite(dat,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")

q("no")


#####################################################################################################

# singularity exec iibdgc_postprocess_10_singularity.sif

# 1.- create a new definition of region (updated wiht Liu etc)

MEM=5000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)

rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

pheno<-c("cd","ibd","uc")

dat<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency.tsv.gz",sep=""))
dat<-as.data.frame(dat)
dim(dat)
# [1] 1080  153

## assign the variants to loci:

all<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_new_and_old_index_gwas_variants_conditional_on_known_index_gwas_fm_2023_with_iibdgc_tier1_tier2_eas_sas_summary_stats_plus_unsupervised_conditional_plus_new_fm_plus_pheno.tsv",sep=""))
all<-as.data.frame(all)
dim(all)

all<-all[,c("updated_region"),drop=F]
all<-all[which(all$updated_region!=""),,drop=F]

all<-all[!duplicated(all),,drop=F]
dim(all)
# [1] 422   1

all$chr<-gsub("_.*","",all$updated_region)
head(all)

all$start<-gsub("^[0-9]{1,2}_","",all$updated_region)
head(all)

all$end<-as.numeric(gsub("^[0-9]*_","",all$start))
all$start<-as.numeric(gsub("_[0-9]*$","",all$start))
head(all)

dat$pos<-gsub("chr[0-9]{1,2}:","",dat$MarkerName)
dat$pos<-as.numeric(gsub(":[A-Z]*:[A-Z]*$","",dat$pos))

dat$updated_region<-NA
for (i in 1:nrow(all)) {

  dat$updated_region[which(dat$chr==all$chr[i] & dat$pos>=all$start[i] & dat$pos<=all$end[i])]<-all$updated_region[i]

}

rm(all)

dim(dat[which(is.na(dat$updated_region)),])
# [1] [1]  16 176

table(dat$chr[which(is.na(dat$updated_region))])
#  2  4  9 12 14 16 18 
#  1  7  2  3  1  1  1 

summary(dat$rate_Neff_cd_eur_tier_2_eas_sas[which(is.na(dat$updated_region))])

# exclude those with Neff<0.5 - rerunning new cojo - excluding those with rate Neff <0.5
dat<-dat[which(dat$rate_Neff_cd_eur_tier_2_eas_sas>=0.5 | dat$rate_Neff_uc_eur_tier_2_eas_sas>=0.5 | dat$rate_Neff_ibd_eur_tier_2_eas_sas>=0.5),] 

dim(dat[which(is.na(dat$updated_region)),])
# [1]  15 155

dim(dat)
# [1] 1077  155

dim(dat[duplicated(dat$MarkerName),])
# [1]   0 155


# add the known signals to identify any LD with those:
# retain the old loci, and the new loci that do not come from conditional eur_tier_2
old<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_new_and_old_index_gwas_variants_conditional_on_known_index_gwas_fm_2023_with_iibdgc_tier1_tier2_eas_sas_summary_stats_plus_unsupervised_conditional.tsv",sep=""))
old<-as.data.frame(old)
dim(old)
# [1] 744 259
old$old_region_signal_name[which(old$MarkerName=="chr17:40619224:C:T")]<-"Cordero_region259_noFM|Liu_region328_noFM"
old$'FM study'[which(old$MarkerName=="chr17:40619224:C:T")]<-"Cordero|Liu"
old$old_region_signal_name[which(old$MarkerName=="chr8:10915904:A:G")]<-"Cordero_region245_noFM|Cordero_region256_noFM"
old$'FM study'[which(old$MarkerName=="chr8:10915904:A:G")]

old<-old[!duplicated(old$MarkerName),]
dim(old)
# [1] 742 259

table(old$class)
# (old)   new   old 
#    31   152   559

table(old$class[which(!old$updated_region %in% dat$updated_region)])
# (old)   new   old 
#     1    30    48

dim(old[which(old$class=="new" & (!old$updated_region %in% dat$updated_region) & (old$"P-value_ibd_eur_tier_2_eas_sas"<=5E-8 | old$"P-value_cd_eur_tier_2_eas_sas"<=5E-8 | old$"P-value_uc_eur_tier_2_eas_sas"<=5E-8)),])
# [1]  21 259


old<-old[which( (old$class_signal %in% c("old","(old)")) | (old$class=="new" & !(old$updated_region %in% dat$updated_region))),]
table(old$class,old$class_signal)
  #       (old) new_signal_new_region new_signal_new_region_to_inspect old
  # (old)    31                     0                                0   0
  # new       0                    24                                6   0
  # old       0                     0                                0 389


old$class_signal[which(old$class_signal=="new_signal_new_region")]<-"new_signal_new_region_cojo_informed"
old$class_signal[which(old$class_signal=="new_signal_new_region_to_inspect")]<-"new_signal_new_region_to_inspect_cojo_informed"
table(old$class,old$class_signal)
  #       (old) new_signal_new_region_cojo_informed
  # (old)    31                                   0
  # new       0                                  24
  # old       0                                   0
       
  #       new_signal_new_region_to_inspect_cojo_informed old
  # (old)                                              0   0
  # new                                                6   0
  # old                                                0 389

list_ids<-c(dat$MarkerName,old$MarkerName)
length(list_ids)
# [1] 1525

list_ids<-list_ids[!duplicated(list_ids)]
length(list_ids)
# [1] 1452

write.table(list_ids,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_variants_for_ld.tsv",sep=""),col.names=T,row.names=F,quote=F,sep="\t")


colnames(dat)[!colnames(dat) %in% colnames(old)]

old$BETA_cond_cd<-NA
old$SE_cond_cd<-NA
old$Pvalue_cond_cd<-NA
old$BETA_cond_ibd<-NA
old$SE_cond_ibd<-NA
old$Pvalue_cond_ibd<-NA
old$BETA_cond_uc<-NA
old$SE_cond_uc<-NA
old$Pvalue_cond_uc<-NA

dat$class<-"new"
dat$class[which(dat$updated_region %in% old$updated_region[which(old$class %in% c("old","(old)"))])]<-"old"

dat$class_signal<-NA

dat$class_signal[which(dat$MarkerName %in% old$MarkerName[which(old$class_signal=="old")])]<-"old|new_cojo_unsupervised"
dat$class_signal[which(dat$MarkerName %in% old$MarkerName[which(old$class_signal=="(old)")])]<-"old|new_cojo_unsupervised"

old<-old[which(!old$MarkerName %in% dat$MarkerName),]
old<-old[,colnames(dat)]

dat<-rbind(old,dat)
table(dat$class_signal,useNA="ifany")
#                                          (old) 
#                                             19 
#            new_signal_new_region_cojo_informed 
#                                             24 
# new_signal_new_region_to_inspect_cojo_informed 
#                                              6 
#                                            old 
#                                            326 
#                      old|new_cojo_unsupervised 
#                                             73 
#                                           <NA> 
#                                           1004 

fwrite(dat,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency_with_old_data.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")

q("no")

####################################

# estimate LD between these variants:


path_gwas=/path/to/ibdgwas/IIBDGC/
ph=ibd

# note this will creat hardcalls, and those will be used for ld estimation, all variants good imputation, it should not be biased...

MEM=2000
for chr in {1..22} 
do
bsub -J"cojo_i" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-o ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_ld1_stdout \
-e ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_ld1_stderr \
"/path/to/software/username/./plink2 \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_subset_included_in_${ph}_analysis \
--extract ${path_gwas}post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_variants_for_ld.tsv \
--make-bed --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_list_tier_2_unsupervised_conditional_variants_for_ld"
done

for chr in {1..22} 
do
echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_ld1_stdout | grep "Successfully"
done

for chr in {1..22} 
do
bsub -J"cojo_i" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-o ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_ld2_stdout \
-e ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_ld2_stderr \
"/path/to/software/username/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_list_tier_2_unsupervised_conditional_variants_for_ld \
--r2 --ld-window-kb '1000000000' --ld-window '1000000000' --ld-window-r2 '0' --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_list_tier_2_unsupervised_conditional_variants_for_ld"
done

for chr in {1..22} 
do
echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_ld2_stdout | grep "Successfully"
done



####################################

# MEM=5000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G humgen-priority -q normal R

rm(list=ls())

library(data.table)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

all<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency_with_old_data.tsv.gz",sep=""))
all<-as.data.frame(all)
dim(all)
# [1] 1452  157

# exclude Cordero_region243_noFM
# all<-all[which(all$MarkerName!="chr1:151281565:T:A"),]
# dim(all)
# # [1] 1451  157

all<-all[!duplicated(all),]
dim(all)
# [1] 1452  157

dim(all[which(all$updated_region==""),c("chr","pos","MarkerName","class_signal")])
# [1] 18  5 - 


all[which(all$updated_region==""),c("chr","pos","MarkerName","class")]
#      chr       pos                       MarkerName class_signal
# 2641  12   8922106                chr12:8922106:G:A         <NA> - new loci1, no overlap - new locus 12_8422106_10257536
# 2661  12   9736259                chr12:9736259:G:A         <NA> - new loci1, no overlap - new locus 12_8422106_10257536
# 267   12   9757536                chr12:9736259:G:A         <NA> - new loci1, no overlap - new locus 12_8422106_10257536
# 2961  14 106697560              chr14:106697560:G:T         <NA> - new loci2, no overlap - new locus 14_106197560_107197560
# 3931  16  69627624 chr16:69627624:TCAGGCTGCATAGTG:T         <NA> - new loci3, no overlap  - new locus 16_69127624_70127624
# 459   18  22486192               chr18:22486192:G:A         <NA> - new loci4, no overlap - new locus 18_21986192_22986192
# 5651   2 230982974               chr2:230982974:G:A         <NA> - 2_229716690_230796312 - new intervals: 230482974_231482974
# 734    4  42396369                chr4:42396369:G:T         <NA> - new loci5, no overlap - new locus 4_41896369_42896369
# 7351   4  42413003                chr4:42413003:C:T         <NA> - new loci5, no overlap - new locus 4_41896369_42896369
# 736    4  42459885                chr4:42459885:G:A         <NA> - new loci5, no overlap - new locus 4_41896369_42896369
# 7371   4  42516358               chr4:42516358:GA:G         <NA> - new loci5, no overlap - new locus 4_41896369_42896369
# 738    4  42635348            chr4:42635348:ATATC:A         <NA> - new loci5, no overlap - new locus 4_41896369_42896369
# 7391   4  42689519                chr4:42689519:A:G         <NA> - new loci5, no overlap - new locus 4_41896369_42896369
# 1076   9  94254276                chr9:94254276:T:C         <NA> - new loci6, no overlap - new locus 9_93754276_94754276
# 1077   9  94268783                chr9:94268783:T:C         <NA> - new loci6, no overlap - new locus 9_93754276_94754276



# manually re-do the new loci:
all$updated_region[which(all$MarkerName %in% c("chr12:8922106:G:A","chr12:9736259:G:A","chr12:9757536:C:A"))]<-"12_8422106_10257536"
all$updated_region[which(all$MarkerName %in% c("chr14:106697560:G:T"))]<-"14_106197560_107197560"
all$updated_region[which(all$MarkerName %in% c("chr16:69627624:TCAGGCTGCATAGTG:T"))]<-"16_69127624_70127624"
all$updated_region[which(all$MarkerName %in% c("chr18:22486192:G:A"))]<-"18_21986192_22986192"
all$updated_region[which(all$MarkerName %in% c("chr4:42396369:G:T","chr4:42413003:C:T","chr4:42459885:G:A","chr4:42516358:GA:G","chr4:42635348:ATATC:A","chr4:42689519:A:G"))]<-"4_41896369_42896369"
all$updated_region[which(all$MarkerName %in% c("chr9:94254276:T:C","chr9:94268783:T:C"))]<-"9_93754276_94754276"

# update intervals of old loci:
all$updated_region[which(all$MarkerName %in% c("chr2:230982974:G:A"))]<-"2_229716690_231482974"
all$updated_region[which(all$updated_region=="2_229716690_230796312")]<-"2_229716690_231482974"

all<-all[which(all$updated_region!=""),] # some old regions by cordero


# check for duplicated variants
dups<-all[duplicated(all$MarkerName),"MarkerName"]
length(dups)
# [1] 0

for (chr in 1:22) {
  ld.tmp<-read.table(paste0(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr",chr,"_list_tier_2_unsupervised_conditional_variants_for_ld.ld"),head=T)
  ld.tmp<-ld.tmp[which(ld.tmp$R2>=0.1),]

  if(chr==1) {
    ld<-ld.tmp
  } else {
    ld<-rbind(ld,ld.tmp)
  }
  rm(ld.tmp)
}


ld<-ld[which( (ld$SNP_A %in% all$MarkerName) & (ld$SNP_B %in% all$MarkerName)),]
all$variants_in_ld_r2_0.1<-NA

all<-all[order(all$chr,all$pos,decreasing=F),]

for (i in 1:nrow(ld)) {

  tmp<-all[which(all$MarkerName %in% c(as.character(ld$SNP_A[i]),as.character(ld$SNP_B[i]))),]

  if (nrow(tmp)>1) {

    # add LD:
    all$variants_in_ld_r2_0.1[which(all$MarkerName %in% c(as.character(ld$SNP_A[i]),as.character(ld$SNP_B[i])))]<-paste(c(as.character(ld$SNP_A[i]),as.character(ld$SNP_B[i]),ld$R2[i]),collapse="|")
    
    # flag as ld_old those in LD with known signal:
    
    # retain only the most signif if both are new:

    if (length(tmp$class_signal[which(tmp$class_signal=="")])==length(tmp$class_signal)) {
      min_pval<-min(c(tmp$Pvalue_cond_cd,tmp$Pvalue_cond_uc,tmp$Pvalue_cond_ibd),na.rm=T)
      new_marker<-tmp$MarkerName[which(tmp$Pvalue_cond_cd==min_pval | tmp$Pvalue_cond_uc==min_pval | tmp$Pvalue_cond_ibd==min_pval)]
      nonew_marker<-tmp$MarkerName[which(!tmp$MarkerName %in% new_marker)]

      all$class_signal[which(all$MarkerName %in% new_marker)]<-"new_cojo_unsupervised"
      all$class_signal[which(all$MarkerName %in% nonew_marker)]<-"new_cojo_unsupervised_in_ld_new_cojo_unsupervised"
      rm(min_pval,new_marker,nonew_marker)

    } else {
      all$class_signal[which(all$MarkerName %in% tmp$MarkerName[which(tmp$class_signal=="")])]<-"new_cojo_unsupervised_in_ld_with_old"
    }
  }

}


all$class_signal[which(all$class_signal=="new_signal_new_region_cojo_informed")]<-"new_signal_new_region_cojo_supervised"
all$class_signal[which(all$class_signal=="new_signal_new_region_to_inspect_cojo_informed")]<-"new_signal_new_region_to_inspect_cojo_supervised"

table(all$class_signal)
                                               
#                                               312 
#                                             (old) 
#                                                19 
#                             new_cojo_unsupervised 
#                                               148 
# new_cojo_unsupervised_in_ld_new_cojo_unsupervised 
#                                               148 
#              new_cojo_unsupervised_in_ld_with_old 
#                                               396 
#             new_signal_new_region_cojo_supervised 
#                                                24 
#  new_signal_new_region_to_inspect_cojo_supervised 
#                                                 6 
                        #                       old 
                        #                       323 
                        # old|new_cojo_unsupervised 
                        #                        73 


tmp<-all[which(all$class_signal==""),]
tmp[which(is.na(tmp$Pvalue_cond_cd) & is.na(tmp$Pvalue_cond_uc) & is.na(tmp$Pvalue_cond_ibd)),]
# 0


# all come from conditional unsupervised:
all$class_signal[which(all$class_signal=="")]<-"new_cojo_unsupervised"

# all come from conditional unsupervised in LD with known signals:
all$class_signal[which(all$class_signal=="(old)")]<-"new_cojo_supervised_in_ld_with_old"

# exclude those previously labelled as to inspect:
all<-all[which(all$class_signal!="new_signal_new_region_to_inspect_cojo_supervised"),]
dim(all)
# [1] 1443  158

fwrite(all,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")


reg<-all[,c("updated_region","class")]

reg$class[which(reg$class=="(old)")]<-"old"
reg<-reg[,c("updated_region","class")]
reg<-reg[!duplicated(reg),]
dim(reg)
# [1] 423   2

# how many loci:
old<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_new_and_old_index_gwas_variants_conditional_on_known_index_gwas_fm_2023_with_iibdgc_tier1_tier2_eas_sas_summary_stats_plus_unsupervised_conditional.tsv",sep=""))
old<-as.data.frame(old)
old<-old[which(old$updated_region!="")]

# # exclude those to inspect:
# old<-old[which(old$class_signal!="new_signal_new_region_to_inspect"),]
table(old$class)
# (old)   new   old 
#    31   152   561 

old$class[which(old$class=="(old)")]<-"old"

old<-old[,c("updated_region","class")]
old<-old[!duplicated(old),]
dim(old)
# [1] 424   2

# exclude the ones we agreed in the call:
old<-old[which(!old$updated_region %in% c("6_24782949_28504154","6_29357657_37435948","X_616740_2923444","")),]
dim(old)
# [1] 420   2

old[which(old$updated_region=="7_69936720_70961219"),]
#          updated_region class
# 683 7_69936720_70961219   new
# 684 7_69936720_70961219   old

old<-old[!duplicated(old$updated_region),]
dim(old)
# [1] 419   2

old[which(!old$updated_region %in% reg$updated_region),]
#            updated_region class
# 412 2_229716690_230796312   old - the one that changed intervals

dim(reg)
# [1] [1] 423   2

reg[which(!reg$updated_region %in% old$updated_region),]
#              updated_region class
# 714     12_8422106_10257536   new
# 746  14_106197560_107197560   new
# 843    16_69127624_70127624   new
# 909    18_21986192_22986192   new
# 1013  2_229716690_231482974   old
# 1015  2_229716690_231482974   new
# 1184    4_41896369_42896369   new
# 1305    6_24782949_28504154   new
# 380     6_29357657_37435948   old
# 1526    9_93754276_94754276   new

# exclude the three we agreed in the call:
reg<-reg[which(!reg$updated_region %in% c("6_24782949_28504154","6_29357657_37435948","X_616740_2923444","")),]
dim(reg)
# [1] 421   2

# exclude the duplicated:
reg<-reg[!duplicated(reg$updated_region),]
dim(reg)
# [1] 420   2

reg[which(!reg$updated_region %in% old$updated_region),]
#              updated_region class
# 714     12_8422106_10257536   new
# 746  14_106197560_107197560   new
# 843    16_69127624_70127624   new
# 909    18_21986192_22986192   new
# 1013  2_229716690_231482974   old - interval updated
# 1184    4_41896369_42896369   new
# 1526    9_93754276_94754276   new

dim(reg)
# [1] 420   2

fwrite(reg,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_regions_420_cojo_supervised_plus_unsupervised.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")


###########################################################################
# to count the number of new signals - not accounted for known signals:

dim(all)
# [1] 1443  158

# keep only the 420 regions:
all<-all[which(all$updated_region %in% reg$updated_region),]
dim(all)
# [1] 1374  158

length(names(table(all$updated_region)))
# [1] 420 - OK

table(all$class_signal)

dim(all[which( !(all$class_signal %in% c("new_cojo_unsupervised_in_ld_with_old","new_cojo_unsupervised_in_ld_new_cojo_unsupervised","old","(old)","old|new_cojo_unsupervised","new_cojo_supervised_in_ld_with_old")) ),])
# [1] 458 158

table(all$class_signal[which( !(all$class_signal %in% c("new_cojo_unsupervised_in_ld_with_old","new_cojo_unsupervised_in_ld_new_cojo_unsupervised","old","(old)","old|new_cojo_unsupervised","new_cojo_supervised_in_ld_with_old")) )])
#                            new_cojo_unsupervised 
#                                              434 
#            new_signal_new_region_cojo_supervised 
#                                               24 

# none of them come from new_cojo_unsupervised signals
all[which(all$class_signal %in% c("new_signal_new_region_cojo_supervised")),1:10]

all$updated_region[which(all$class_signal %in% c("new_signal_new_region_cojo_supervised"))]
# [1]  [1] "1_15325338_16449012"    "1_45712246_46717401"    "10_44514562_45514563"  
#  [4] "10_71205481_72205481"   "10_105121605_106144359" "11_10145868_11167716"  
#  [7] "11_124471689_125471689" "12_30154261_31154262"   "12_124048792_125048793"
# [10] "14_23646226_24646227"   "14_38616772_39690790"   "14_101322326_102322326"
# [13] "15_44302938_45302939"   "15_100444058_101448378" "2_99645086_100677415"  
# [16] "21_34850264_35850265"   "22_22409672_23409673"   "4_6539650_7542654"     
# [19] "4_53695959_54727865"    "5_3922049_4922050"      "5_88197056_89330037"   
# [22] "6_87193208_88227520"    "7_48450632_49522286"    "9_29200044_30202045" 

# any of those with FM results?

# regions with FM data:
fm_reg<-fread("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/fine_mapping/fine_mapping_region_reclassification.csv")
dim(fm_reg)
# [1] 419   4 - to be udpated to 420
fm_reg$region_label[which(fm_reg$region_label=="2_229716690_230796312")]<-"2_229716690_231482974"
table(fm_reg$reason_noCS)

table(fm_reg$reason_noCS[which(fm_reg$region_label %in% all$updated_region[which(all$class_signal %in% c("new_signal_new_region_cojo_supervised"))])])
                                #  . classC1_nosignificant_after_filter 
                                #  3                                 21 


length(names(table(all$updated_region)))
# [1] 420

table(all$class_signal[which( !(all$class_signal %in% c("new_cojo_unsupervised_in_ld_with_old","new_cojo_unsupervised_in_ld_new_cojo_unsupervised","old","(old)","old|new_cojo_unsupervised","new_cojo_supervised_in_ld_with_old")) )])
                # new_cojo_unsupervised new_signal_new_region_cojo_supervised 
                #                   434                                    24


dim(all[which(all$class_signal=="new_cojo_unsupervised"),])
# [1] 434 158

# across how many regions:
tmp<-all[which(all$class_signal=="new_cojo_unsupervised"),c("updated_region","class")]
tmp<-tmp[!duplicated(tmp),]
tmp$class[which(tmp$updated_region=="2_229716690_231482974")]<-"old"
tmp[duplicated(tmp$updated_region),]
table(tmp$class)
# new old 
# 120 145 


# all new_cojo_unsupervised have significant unsupervised p-values
dim(all[which(all$class_signal=="new_cojo_unsupervised" & is.na(all$Pvalue_cond_cd) & is.na(all$Pvalue_cond_uc) & is.na(all$Pvalue_cond_ibd)),1:10])
# [1]  0 10

# new_signal_new_region_cojo_supervised - are those singificant in EUR tier 1:
tmp<-all[which(all$class_signal=="new_signal_new_region_cojo_supervised"),]

dim(tmp[which(tmp$"P-value_cd_eur_tier_1"<=5E-8 | tmp$"P-value_uc_eur_tier_1"<=5E-8 | tmp$"P-value_ibd_eur_tier_1"<=5E-8),])
# [1]  2 158

dim(tmp[which(tmp$"P-value_cd_eur_tier_2"<=5E-8 | tmp$"P-value_uc_eur_tier_2"<=5E-8 | tmp$"P-value_ibd_eur_tier_2"<=5E-8),])
# [1]  22 158

dim(tmp[which(tmp$"P-value_cd_eur_tier_2_eas_sas"<=5E-8 | tmp$"P-value_uc_eur_tier_2_eas_sas"<=5E-8 | tmp$"P-value_ibd_eur_tier_2_eas_sas"<=5E-8),])
# [1]  24 158

tmp$updated_region[!duplicated(tmp$updated_region)]
#  [1] "1_15325338_16449012"    "1_45712246_46717401"    "10_44514562_45514563"  
#  [4] "10_71205481_72205481"   "10_105121605_106144359" "11_10145868_11167716"  
#  [7] "11_124471689_125471689" "12_30154261_31154262"   "12_124048792_125048793"
# [10] "14_23646226_24646227"   "14_38616772_39690790"   "14_101322326_102322326"
# [13] "15_44302938_45302939"   "15_100444058_101448378" "2_99645086_100677415"  
# [16] "21_34850264_35850265"   "22_22409672_23409673"   "4_6539650_7542654"     
# [19] "4_53695959_54727865"    "5_3922049_4922050"      "5_88197056_89330037"   
# [22] "6_87193208_88227520"    "7_48450632_49522286"    "9_29200044_30202045"

# add labels for the old signals:

old<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_new_and_old_index_gwas_variants_conditional_on_known_index_gwas_fm_2023_with_iibdgc_tier1_tier2_eas_sas_summary_stats_plus_unsupervised_conditional_frequency.tsv.gz",sep=""))
old<-as.data.frame(old)

old<-old[,c("old_region_signal_name","FM study","MarkerName")]
old<-old[which(old$old_region_signal_name!=""),]

# for duplicated known markers keep only one
vec<-old$MarkerName[which(duplicated(old$MarkerName))]
old[which(old$MarkerName %in% vec),]
#     old_region_signal_name FM study         MarkerName
# 307 Cordero_region259_noFM  Cordero chr17:40619224:C:T
# 308     Liu_region328_noFM      Liu chr17:40619224:C:T
# 691 Cordero_region256_noFM  Cordero  chr8:10915904:A:G
# 692 Cordero_region245_noFM  Cordero  chr8:10915904:A:G

old$old_region_signal_name[which(old$MarkerName=="chr17:40619224:C:T")]<-"Cordero_region259_noFM|Liu_region328_noFM"
old$'FM study'[which(old$MarkerName=="chr17:40619224:C:T")]<-"Cordero|Liu"
old$old_region_signal_name[which(old$MarkerName=="chr8:10915904:A:G")]<-"Cordero_region245_noFM|Cordero_region256_noFM"

old<-old[!duplicated(old),]

old[which(!old$MarkerName %in% all$MarkerName),]
#                            old_region_signal_name FM study         MarkerName
# 109                        Cordero_region258_noFM  Cordero chr10:48916679:G:A
# 629                        Cordero_region252_noFM  Cordero  chr6:29636347:G:A
# 630                         Sakaue_region260_noFM   Sakaue  chr6:32247045:T:G
# 631                     Goyette_region100_signal1  Goyette  chr6:32644620:A:G
# 632                        Cordero_region253_noFM  Cordero  chr6:33908305:C:T
# 634                        Cordero_region254_noFM  Cordero  chr6:34550494:G:A
# 635                        Cordero_region255_noFM  Cordero  chr6:35089554:G:A
# 691 Cordero_region245_noFM|Cordero_region256_noFM  Cordero  chr8:10915904:A:G
# 714                        Cordero_region257_noFM  Cordero  chr8:78765490:C:A

old<-old[which(old$MarkerName %in% all$MarkerName),]

dim(all)
# [1] 1374  158
all<-merge(old,all,by="MarkerName",all=T)
dim(all)
# [1] 1374  160

table(all$class)
# (old)   new   old 
#    19   218  1137 
all$class[which(all$class=="(old)")]<-"old"
table(all$class)
#  new  old 
#  218 1156 

fwrite(all,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")


new<-all[which(is.na(all$old_region_signal_name)),]
dim(new)
# [1] 996 160

table(new$class_signal)
#                new_cojo_supervised_in_ld_with_old 
#                                                19 
#                             new_cojo_unsupervised 
#                                               434 
# new_cojo_unsupervised_in_ld_new_cojo_unsupervised 
#                                               137 
#              new_cojo_unsupervised_in_ld_with_old 
#                                               370 
#             new_signal_new_region_cojo_supervised 
#                                                24 
#                         old|new_cojo_unsupervised 
#                                                12 
table(new$class)
# new old 
# 218 778

new$class[which(new$class)]
fwrite(new,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_new_signals.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")


# create two separate tables from this, one with known regions and signals, a second one with new signals at this stage:

# known signals
dim(all[which(all$old_region_signal_name!=""),])
# [1] 378 160

old<-all[which(all$old_region_signal_name!=""),]
table(old$class_signal)
                      # old old|new_cojo_unsupervised 
                      # 316                        61

# out of these we used:
dim(old[which(!is.na(old$Position_b38_eur_tier_1)),])
# [1] 370 160

table(old$class)
# old 
# 378 

fwrite(old,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_old_regions_and_signals.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")

q("no")


##################

# old regions that no longer are required:

# MEM=5000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

rm(list=ls())

library(data.table)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

old<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_new_and_old_index_gwas_variants_conditional_on_known_index_gwas_fm_2023_with_iibdgc_tier1_tier2_eas_sas_summary_stats_plus_unsupervised_conditional.tsv",sep=""))
old<-as.data.frame(old)
old<-old[which(old$updated_region!="")]

# # exclude those to inspect:
# old<-old[which(old$class_signal!="new_signal_new_region_to_inspect"),]
table(old$class)
# (old)   new   old 
#    31   152   561 

old$class[which(old$class=="(old)")]<-"old"

old<-old[,c("updated_region","class")]
old<-old[!duplicated(old),]
dim(old)
# [1] 424   2

# exclude the ones we agreed in the call:
old<-old[which(!old$updated_region %in% c("6_24782949_28504154","6_29357657_37435948","X_616740_2923444","")),]
dim(old)
# [1] 420   2

old[which(old$updated_region=="7_69936720_70961219"),]
#          updated_region class
# 683 7_69936720_70961219   new
# 684 7_69936720_70961219   old

old<-old[!duplicated(old$updated_region),]
dim(old)
# [1] 419   2

reg<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_regions_420_cojo_supervised_plus_unsupervised.tsv.gz",sep=""))

reg[which(!reg$updated_region %in% old$updated_region),]
          #  updated_region  class
# 1:    12_8422106_10257536    new
# 2: 14_106197560_107197560    new
# 3:   16_69127624_70127624    new
# 4:   18_21986192_22986192    new
# 5:  2_229716690_231482974    old
# 6:    4_41896369_42896369    new
# 7:    9_93754276_94754276    new


old[which(!old$updated_region %in% reg$updated_region),]
#            updated_region class
# 106  10_42619246_43753122   new
# 202  12_92237438_93240527   new
# 215  13_53996833_54996834   new
# 412 2_229716690_230796312   old
# 465  21_18171986_19171986   new
# 577     5_2840882_3846018   new

q("no")
