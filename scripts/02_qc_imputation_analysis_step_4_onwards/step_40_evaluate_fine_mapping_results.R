# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
############################################################################################################################
############################################################################################################################

# two different versions of FM resutls stored here:
# /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/fine_mapping/

# latest:
# tier1_EUR_array_data_fine_mapping_results_08_16_2025_*

# 1.- ADD THE NEW FINE MAPPING RESULTS FROM RUI:

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=22000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

all<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal.tsv.gz",sep=""))
all<-as.data.frame(all)
dim(all)
# [1] 1374 160

fm<-fread("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/fine_mapping/tier1_EUR_array_data_fine_mapping_results_08_16_2025_list_of_variants.csv")

##############################
# merge with FM

fm$X<-paste(fm$region_label,fm$CS_index,sep="_")

# 1.- if lead variant (or variant in credible set) in index table, link by that variant:
tmp<-fm[which( (fm$MarkerName %in% all$MarkerName)),]
dim(tmp)
# [1] 433  57
length(names(table(tmp$X)))
# [1] 237 - duplicated

dups<-tmp$MarkerName[duplicated(tmp$MarkerName)]
dups
# [1] "chr1:197547956:T:A" "chr1:197732149:G:C" "chr2:102322576:C:T"
# [4] "chr4:3447925:G:A"   "chr5:132336964:T:C" "chr5:132435113:C:T"
# [7] "chr11:76582714:G:A" "chr16:28534266:G:C" "chr16:28569012:T:C"

tmp<-tmp[!duplicated(tmp$MarkerName),]
dim(tmp)
# [1] 424  57

all$Y<-paste(all$updated_region,all$MarkerName,sep="_")
tmp$Y<-paste(tmp$region_label,tmp$MarkerName,sep="_")

all1<-merge(all,tmp,by="Y")
dim(all1)
# [1] 422 218


# 2.- replace the old lead variant by a new one - in those cases where they do not fully align:
tmp2<-tmp[which(tmp$lead=="FALSE"),]
dim(tmp2)
# [1] 299  58

dat1<-fm[which(fm$X %in% tmp2$X),]
dat1<-dat1[which(dat1$lead==TRUE),]
dim(dat1)
# [1] 169  53


# 3.- for any other signal, merge by lead from FM
dat2<-fm[which(!fm$X %in% tmp$X),]

dat2<-dat2[which(dat2$lead==TRUE),]
dim(dat2)
# [1] 33 57


# combine both to extract summary stats:
dat<-rbind(dat1,dat2)
dim(dat)
# [1] 202  53


chr<-as.numeric(gsub("chr","",unlist(labels(table(dat$chr)))))
chr<-chr[which(!is.na(chr))]

pheno<-c("cd","uc","ibd")

###################################################################################################################################################################################

# extract sumamry statistics for those 202 variants:

##############################################################
# results from EUR tier 1

for (i in 1:length(pheno)) {
  
  print(pheno[i])
  
  for (j in 1:length(chr)) {
    
    tmp1<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr[j],"_",pheno[i],"_meta_eur_tier_1_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff_with_rsid_dbsnp154.txt.gz",sep=""),head=T)
    tmp1<-tmp1[which(tmp1$MarkerName %in% dat$MarkerName),]
    
    if(j==1) {
      tmp<-tmp1
    } else {
      tmp<-rbind(tmp,tmp1)
    }
    
    rm(tmp1)
  
  }
  
  if(i==1) {
    tmp<-tmp[,c("MarkerName","dbsnp154","Position_b38","A1","A2","BETA","SE","P-value","Direction_ed","HetISq",
                "HetChiSq","HetDf","HetPVal","rate_total_sample","avgA2FREQ","avgA2FREQ_CASES","avgA2FREQ_CONTROLS","Neff","rate_Neff","N_CASES","N_CONTROLS")]
    colnames(tmp)[2:4]<-paste(colnames(tmp)[2:4],"eur_tier_1",sep="_")
    colnames(tmp)[5:ncol(tmp)]<-paste(colnames(tmp)[5:ncol(tmp)],pheno[i],"eur_tier_1",sep="_")
    all<-tmp
  } else {
    tmp<-tmp[,c("MarkerName","BETA","SE","P-value","Direction_ed","HetISq",
                "HetChiSq","HetDf","HetPVal","rate_total_sample","avgA2FREQ","avgA2FREQ_CASES","avgA2FREQ_CONTROLS","Neff","rate_Neff","N_CASES","N_CONTROLS")]
    colnames(tmp)[2:ncol(tmp)]<-paste(colnames(tmp)[2:ncol(tmp)],pheno[i],"eur_tier_1",sep="_")
    all<-merge(all,tmp,by="MarkerName",sort=F)
  }
  rm(tmp)
}



dat<-as.data.frame(dat)
dat<-merge(dat,all,by="MarkerName",all.x=T)
rm(all)
write.table(dat,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/tmp_list_all_variants_to_define_association_signals.tsv.gz",sep=""),col.names=T,row.names=F,quote=F,sep="\t")


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
                "HetChiSq","HetDf","HetPVal","rate_total_sample","avgA2FREQ","avgA2FREQ_CASES","avgA2FREQ_CONTROLS","Neff","rate_Neff","N_CASES","N_CONTROLS")]
    colnames(tmp)[2:4]<-paste(colnames(tmp)[2:4],"eur_tier_2",sep="_")
    colnames(tmp)[5:ncol(tmp)]<-paste(colnames(tmp)[5:ncol(tmp)],pheno[i],"eur_tier_2",sep="_")
    all<-tmp
  } else {
    tmp<-tmp[,c("MarkerName","BETA","SE","P-value","Direction_ed","HetISq",
                "HetChiSq","HetDf","HetPVal","rate_total_sample","avgA2FREQ","avgA2FREQ_CASES","avgA2FREQ_CONTROLS","Neff","rate_Neff","N_CASES","N_CONTROLS")]
    colnames(tmp)[2:ncol(tmp)]<-paste(colnames(tmp)[2:ncol(tmp)],pheno[i],"eur_tier_2",sep="_")
    all<-merge(all,tmp,by="MarkerName",sort=F)
  }
  rm(tmp)
}

dim(dat)
#[1] 248  84
dim(all)
#[1] 248  31

all<-as.data.frame(all)
dat<-merge(dat,all,by="MarkerName",all.x=T)
rm(all)

write.table(dat,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/tmp_list_all_variants_to_define_association_signals.tsv.gz",sep=""),col.names=T,row.names=F,quote=F,sep="\t")



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
                "HetChiSq","HetDf","HetPVal","rate_total_sample","avgA2FREQ","avgA2FREQ_CASES","avgA2FREQ_CONTROLS","Neff","rate_Neff","N_CASES","N_CONTROLS")]
    colnames(tmp)[2:4]<-paste(colnames(tmp)[2:4],"sas",sep="_")
    colnames(tmp)[5:ncol(tmp)]<-paste(colnames(tmp)[5:ncol(tmp)],pheno[i],"sas",sep="_")
    all<-tmp
  } else {
    tmp<-tmp[,c("MarkerName","BETA","SE","P-value","Direction_ed","HetISq",
                "HetChiSq","HetDf","HetPVal","rate_total_sample","avgA2FREQ","avgA2FREQ_CASES","avgA2FREQ_CONTROLS","Neff","rate_Neff","N_CASES","N_CONTROLS")]
    colnames(tmp)[2:ncol(tmp)]<-paste(colnames(tmp)[2:ncol(tmp)],pheno[i],"sas",sep="_")
    all<-merge(all,tmp,by="MarkerName",sort=F)
  }
  rm(tmp)
}

dim(all)
#[1] 200  49
dim(dat)
#[1]248 150

all<-as.data.frame(all)
dat<-merge(dat,all,by="MarkerName",all.x=T)
rm(all)

write.table(dat,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/tmp_list_all_variants_to_define_association_signals.tsv.gz",sep=""),col.names=T,row.names=F,quote=F,sep="\t")



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
                "HetChiSq","HetDf","HetPVal","rate_total_sample","avgA2FREQ","avgA2FREQ_CASES","avgA2FREQ_CONTROLS","Neff","rate_Neff","N_CASES","N_CONTROLS")]
    colnames(tmp)[2:4]<-paste(colnames(tmp)[2:4],"eur_tier_2_eas_sas",sep="_")
    colnames(tmp)[5:ncol(tmp)]<-paste(colnames(tmp)[5:ncol(tmp)],pheno[i],"eur_tier_2_eas_sas",sep="_")
    all<-tmp
  } else {
    tmp<-tmp[,c("MarkerName","BETA","SE","P-value","Direction_ed","HetISq",
                "HetChiSq","HetDf","HetPVal","rate_total_sample","avgA2FREQ","avgA2FREQ_CASES","avgA2FREQ_CONTROLS","Neff","rate_Neff","N_CASES","N_CONTROLS")]
    colnames(tmp)[2:ncol(tmp)]<-paste(colnames(tmp)[2:ncol(tmp)],pheno[i],"eur_tier_2_eas_sas",sep="_")
    all<-merge(all,tmp,by="MarkerName",sort=F)
  }
  rm(tmp)
}

all<-as.data.frame(all)

dim(all)
# [1] 248  49
dim(dat)
# [1] 248 198

dat<-merge(dat,all,by="MarkerName",all.x=T)
rm(all)
write.table(dat,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/tmp_list_all_variants_to_define_association_signals.tsv.gz",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

# dat<-fread("~/tmp.tsv.gz")
# dat<-as.data.frame(dat)

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

dim(all)
# [1] 248  22
dim(dat)
# [1] 248 246

all<-as.data.frame(all)
dat<-merge(dat,all,by="MarkerName",all.x=T)

rm(all)

write.table(dat,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/tmp_list_all_variants_to_define_association_signals.tsv.gz",sep=""),col.names=T,row.names=F,quote=F,sep="\t")


q("no")

###############

# merge all data together:

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=5000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

dat<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/tmp_list_all_variants_to_define_association_signals.tsv.gz",sep=""))
dat<-as.data.frame(dat)
dim(dat)
# [1] 202 279

all<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal.tsv.gz",sep=""))
all<-as.data.frame(all)
dim(all)
# [1] 1374  160

dim(all[duplicated(all$MarkerName),])
#0 

length(names(table(all$updated_region)))
# [1] 420

fm<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/fine_mapping/tier1_EUR_array_data_fine_mapping_results_08_16_2025_list_of_variants.csv",sep=""))
fm<-as.data.frame(fm)
dim(fm)
# 7221

dim(fm[duplicated(fm$MarkerName),])
# [1]  59 56

dim(fm[duplicated(fm),])
# [1]  0 56

# # rename temporarily the updated region:
# fm$region_label[which(fm$region_label=="2_229716690_230796312")]<-"2_229716690_231482974"


# fm file includes regions later excluded (post FM) - update accordingly:
fm_reg<-fread("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/fine_mapping/fine_mapping_region_reclassification.csv")
dim(fm_reg)
# [1] 419   4

fm_reg$region_label[which(fm_reg$region_label=="2_229716690_230796312")]<-"2_229716690_231482974"
table(fm_reg$reason_noCS)
  #                                .              classA_no_significant 
  #                              170                                  1 
  #                    classB_purity classC1_nosignificant_after_filter 
  #                                1                                231 
  # classC2_significant_after_filter 
  #                               16 

dim(fm)
# [1] 7221   51
fm<-fm[which(fm$region_label %in% fm_reg$region_label[which(fm_reg$reason_noCS %in% c("."))]),]
dim(fm)
# [1] 7221   51


lead<-fm[which(fm$lead==TRUE),"MarkerName"]
lead<-lead[!duplicated(lead)]
length(lead)
# [1] 268

fm$X<-paste(fm$region_label,fm$CS_index,sep="_")


##############################
# merge with FM

# 1.- if lead variant (or variant in credible set) in index table, link by that variant:
tmp<-fm[which( (fm$MarkerName %in% all$MarkerName)),]
dim(tmp)
# [1] 433  53

length(names(table(tmp$X)))
# [1] 237 - 9 duplicated

dups<-tmp$MarkerName[dno_conditional_significant_eur_tier2uplicated(tmp$MarkerName)]
dups
# [1] [1] "chr1:197547956:T:A" "chr1:197732149:G:C" "chr2:102322576:C:T"
# [4] "chr4:3447925:G:A"   "chr5:132336964:T:C" "chr5:132435113:C:T"
# [7] "chr11:76582714:G:A" "chr16:28534266:G:C" "chr16:28569012:T:C"

# keep only fm lead (when >1)


for (i in 1:length(dups)) {
  if (nrow(tmp[which(tmp$MarkerName==dups[i] & tmp$lead==TRUE),])>0) {
    print(dups[i])
    tmp$MarkerName[which(tmp$MarkerName==dups[i] & tmp$lead!=TRUE)]<-paste0(tmp$MarkerName[which(tmp$MarkerName==dups[i] & tmp$lead!=TRUE)],"_rm")
  } else {
    tmp$MarkerName[which(tmp$MarkerName==dups[i])][1]<-paste0(tmp$MarkerName[which(tmp$MarkerName==dups[i])][1],"_rm")
  }
}

dim(tmp)
# [1] 433  53
tmp<-tmp[!grepl("_rm",tmp$MarkerName),]
dim(tmp)
# [1] 424  53

tmp$MarkerName[duplicated(tmp$MarkerName)]

all$Y<-paste(all$updated_region,all$MarkerName,sep="_")
tmp$Y<-paste(tmp$region_label,tmp$MarkerName,sep="_")

dim(tmp[which(tmp$Y %in% all$Y),])
# [1] 424  54

tmp<-tmp[,colnames(tmp)[!colnames(tmp) %in% c("chr","MarkerName","N")]]

dim(all)
# [1] 1374  161

all<-merge(all,tmp,by="Y",all.x=T)
dim(all)
# [1] 1374 210

dim(all[duplicated(all$MarkerName),])
# [1]   0 210

all$chr<-gsub(":[0-9]*:[A-Z]*:[A-Z]*$","",all$MarkerName)
all$chr<-gsub("chr","",all$chr)

# rename class_signal for those wiht same index from new FM:
all$class_signal[which(!is.na(all$trait_CS))]<-paste0(all$class_signal[which(!is.na(all$trait_CS))],"|new_fm_credible_set")
all$class_signal[which(all$class_signal=="|new_fm_credible_set")]<-"new_fm_credible_set"
table(all$class_signal)
#                                    new_cojo_supervised_in_ld_with_old 
#                                                                    16 
#                new_cojo_supervised_in_ld_with_old|new_fm_credible_set 
#                                                                     3 
#                                                 new_cojo_unsupervised 
#                                                                   360 
#                     new_cojo_unsupervised_in_ld_new_cojo_unsupervised 
#                                                                   101 
# new_cojo_unsupervised_in_ld_new_cojo_unsupervised|new_fm_credible_set 
#                                                                    36 
#                                  new_cojo_unsupervised_in_ld_with_old 
#                                                                   195 
#              new_cojo_unsupervised_in_ld_with_old|new_fm_credible_set 
#                                                                   175 
#                             new_cojo_unsupervised|new_fm_credible_set 
#                                                                    74 
#                                 new_signal_new_region_cojo_supervised 
#                                                                    21 
#             new_signal_new_region_cojo_supervised|new_fm_credible_set 
#                                                                     3 
#                                                                   old 
#                                                                   233 
#                                             old|new_cojo_unsupervised 
#                                                                    24 
#                         old|new_cojo_unsupervised|new_fm_credible_set 
#                                                                    49 
#                                               old|new_fm_credible_set 
#                                                                    84 



# rbind the new signals that do not match any previous index variant, either because they are independent signals, or they are the new index
# colnames(all)
# colnames(dat)

dim(dat)
# 202 278

dim(dat[duplicated(dat$MarkerName),])
# [1]   0 278

dat$class_signal<-"new_fm_credible_set"
dat$updated_region<-dat$region_label

dat$pos<-dat$position

# relabel the class:
dat<-dat[,colnames(dat)[!colnames(dat) %in% "class"]]

tmp<-all[,c("updated_region","class")]
tmp$class[which(tmp$class=="(old)")]<-"old"
tmp[duplicated(tmp$updated_region),]
#            updated_region class
# 766 2_229716690_231482974   new

tmp[which(tmp$updated_region=="2_229716690_231482974"),]
#            updated_region class
# 763 2_229716690_231482974   old
# 766 2_229716690_231482974   new

tmp<-tmp[!duplicated(tmp$updated_region),]


tmp<-tmp[!duplicated(tmp),]


dat<-merge(dat,tmp[,c("updated_region","class")],by="updated_region",all.x=T)

cols_to_add<-colnames(all)[!colnames(all) %in% colnames(dat)]
dat[,cols_to_add] <- NA

colnames(dat)[!colnames(dat) %in% colnames(all)]

dat<-dat[,colnames(all)]
dim(dat)
# [1] 202 210

# rename the re-labelled updated region:
dat$updated_region[which(dat$updated_region=="2_229716690_230796312")]<-"2_229716690_231482974"

# all should be in the 420 regions we are keeping:
reg<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_regions_420_cojo_supervised_plus_unsupervised.tsv.gz",sep=""))

dat<-dat[which(dat$updated_region %in% reg$updated_region),]
dim(dat)
# [1] 202 211

dim(dat[duplicated(dat$MarkerName),])
# [1]   0 215

dim(dat[which(dat$MarkerName %in% all$MarkerName),])
# [1]  59 215

dat<-dat[which(!dat$MarkerName %in% all$MarkerName),]


all<-rbind(all,dat)

all<-all[,colnames(all)[which(!colnames(all) %in% c("Y","X"))]]
dim(all)
# [1] 1517  208

# remove X, V5 ... column:
all<-all[,colnames(all)[!colnames(all) %in% c("X","V5","min","max","regions_combined","region_class")],]

length(names(table(all$updated_region)))
# [1] 420

all$class[which(all$updated_region=="2_229716690_231482974")]<-"old"

all$region_fm_2025<-NA
all$region_fm_2025[which(all$updated_region %in% fm_reg$region_label[which(fm_reg$reason_noCS %in% c("classC1_nosignificant_after_filter","classA_no_significant"))])]<-"N_filtered_out_preFM"
all$region_fm_2025[which(all$updated_region %in% fm_reg$region_label[which(fm_reg$reason_noCS %in% c("classB_purity","classC2_significant_after_filter"))])]<-"N_filtered_out_postFM"
all$region_fm_2025[which(all$updated_region %in% fm_reg$region_label[which(fm_reg$reason_noCS %in% c("."))])]<-"Y"

table(all$region_fm_2025,useNA="ifany")

table(all$updated_region[which(is.na(all$region_fm_2025))],all$class_signal[which(is.na(all$region_fm_2025))])
  #                        in_ld_new_cojo_unsupervised new_cojo_unsupervised
  # 12_8422106_10257536                              1                     2
  # 14_106197560_107197560                           0                     1
  # 16_69127624_70127624                             0                     1
  # 18_21986192_22986192                             0                     1
  # 4_41896369_42896369                              1                     5
  # 9_93754276_94754276                              1                     1

all$region_fm_2025[which(is.na(all$region_fm_2025))]<-"N_filtered_out_preFM"

table(all$region_fm_2025,useNA="ifany")
# N_filtered_out_postFM  N_filtered_out_preFM                     Y 
#                    65                   438                  1014 


# identify the fine mapping leads:
length(lead)
# [1] 268

length(lead[which(lead %in% all$MarkerName)])
# [1] 268

table(all$class_signal[which(all$MarkerName %in% lead)])
#                new_cojo_supervised_in_ld_with_old|new_fm_credible_set 
#                                                                     1 
# new_cojo_unsupervised_in_ld_new_cojo_unsupervised|new_fm_credible_set 
#                                                                     4 
#              new_cojo_unsupervised_in_ld_with_old|new_fm_credible_set 
#                                                                    43 
#                             new_cojo_unsupervised|new_fm_credible_set 
#                                                                    35 
#                                                   new_fm_credible_set 
#                                                                   143 
#             new_signal_new_region_cojo_supervised|new_fm_credible_set 
#                                                                     2 
#                         old|new_cojo_unsupervised|new_fm_credible_set 
#                                                                    30 
#                                               old|new_fm_credible_set 
#                                                                    10

for (i in 1:length(lead)) {

  all$class_signal[which(all$MarkerName==lead[i])]<-gsub("new_fm_credible_set","new_fm_credible_set_lead",all$class_signal[which(all$MarkerName==lead[i])])

}

table(all$class_signal[which(all$MarkerName %in% lead)])
#                new_cojo_supervised_in_ld_with_old|new_fm_credible_set_lead 
#                                                                          1 
# new_cojo_unsupervised_in_ld_new_cojo_unsupervised|new_fm_credible_set_lead 
#                                                                          4 
#              new_cojo_unsupervised_in_ld_with_old|new_fm_credible_set_lead 
#                                                                         43 
#                             new_cojo_unsupervised|new_fm_credible_set_lead 
#                                                                         35 
#                                                   new_fm_credible_set_lead 
#                                                                        143 
#             new_signal_new_region_cojo_supervised|new_fm_credible_set_lead 
#                                                                          2 
#                         old|new_cojo_unsupervised|new_fm_credible_set_lead 
#                                                                         30 
#                                               old|new_fm_credible_set_lead 
#                                                                         10 

dim(all)
# [1] 1517  213

fwrite(all,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_with_fm_signals.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")

list_ids<-all$MarkerName
list_ids<-list_ids[!duplicated(list_ids)]
length(list_ids)
# [1] 1517

write.table(list_ids,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_variants_for_ld_with_fm_with_old_data.tsv",sep=""),col.names=T,row.names=F,quote=F,sep="\t")


q("no")

#############################################################################

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
--extract ${path_gwas}post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_variants_for_ld_with_fm_with_old_data.tsv \
--make-bed --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_list_tier_2_unsupervised_conditional_variants_for_ld_with_fm_with_old_data"
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
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_list_tier_2_unsupervised_conditional_variants_for_ld_with_fm_with_old_data \
--r2 --ld-window-kb '1000000000' --ld-window '1000000000' --ld-window-r2 '0' --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_list_tier_2_unsupervised_conditional_variants_for_ld_with_fm_with_old_data"
done

for chr in {1..22} 
do
echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_ld2_stdout | grep "Successfully"
done


#################################################################################################################################################

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=5000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

all<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_with_fm_signals.tsv.gz",sep=""))
all<-as.data.frame(all)
dim(all)
# [1] 1517  160


# for regions with new FM data - retain just the FM signals
table(all$class_signal,all$region_fm_2025,useNA="ifany")

# for the 170 regions successfully FM
length(names(table(all$updated_region[which(all$region_fm_2025=="Y")])))
# [1] 170


#########################################
# REGIONS THAT HAVE BEEN FINE MAPPED:

fm<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/fine_mapping/tier1_EUR_array_data_fine_mapping_results_06_01_2025_list_of_variants.csv",sep=""))
fm<-as.data.frame(fm)
dim(fm)
# 7221

lead<-fm[which(fm$lead==TRUE),"MarkerName"]
lead<-lead[!duplicated(lead)]
length(lead)
# [1] 268
rm(fm)

all$class_signal_postfm<-NA

all1<-all[which((all$region_fm_2025=="Y" & all$class_signal %in% names(table(all$class_signal[which(all$MarkerName %in% lead)])))),]
length(names(table(all1$updated_region)))
# [1] 170

# keep only the ones that coincide with the marker lead
dim(all1[which(all1$lead==TRUE),])
# [1] [1] 268 213

all1<-all1[which(all1$lead==TRUE),]

table(all1$class_signal)
#                new_cojo_supervised_in_ld_with_old|new_fm_credible_set_lead 
#                                                                          1 
# new_cojo_unsupervised_in_ld_new_cojo_unsupervised|new_fm_credible_set_lead 
#                                                                          4 
#              new_cojo_unsupervised_in_ld_with_old|new_fm_credible_set_lead 
#                                                                         43 
#                             new_cojo_unsupervised|new_fm_credible_set_lead 
#                                                                         35 
#                                                   new_fm_credible_set_lead 
#                                                                        143 
#             new_signal_new_region_cojo_supervised|new_fm_credible_set_lead 
#                                                                          2 
#                         old|new_cojo_unsupervised|new_fm_credible_set_lead 
#                                                                         30 
#                                               old|new_fm_credible_set_lead 
#                                                                         10 

all1$class_signal_postfm<-"new_fm_credible_set_lead"
table(all1$class_signal_postfm)
# new_fm_credible_set_lead 
#                      268

# old known signals
all1_tmp<-all[which(all$old_region_signal_name!=""),]
dim(all1_tmp)
# [1] 378 213

dim(all1_tmp[which(all1_tmp$MarkerName %in% all1$MarkerName),])
# [1]  38 213

all1_tmp<-all1_tmp[which(!all1_tmp$MarkerName %in% all1$MarkerName),]
ids<-c(all1_tmp$MarkerName,all1$MarkerName)
ids<-ids[!duplicated(ids)]


for (chr in 1:22) {
  ld.tmp<-read.table(paste0(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr",chr,"_list_tier_2_unsupervised_conditional_variants_for_ld_with_fm_with_old_data.ld"),head=T)
  ld.tmp<-ld.tmp[which(ld.tmp$R2>=0.1),]

  if(chr==1) {
    ld<-ld.tmp
  } else {
    ld<-rbind(ld,ld.tmp)
  }
  rm(ld.tmp)
}


ld<-ld[which( (ld$SNP_A %in% ids) & (ld$SNP_B %in% ids)),]
all1$variants_in_ld_r2_0.1<-NA

all1$class[which(all1$class=="(old)")]<-"old"
table(all1$class)
# new old 
#  21 247 


for (i in 1:nrow(all1)) {

  tmp<-ld[which(ld$SNP_A %in% all1$MarkerName[i] | (ld$SNP_B %in% all1$MarkerName[i])),]
  tmp<-tmp[which(tmp$SNP_A %in% all1_tmp$MarkerName | (tmp$SNP_B %in% all1_tmp$MarkerName)),]

  tmp$ld<-paste(as.character(tmp$SNP_A),as.character(tmp$SNP_B),tmp$R2,sep="|")

  if (nrow(tmp)>0) {

    # print(i)
    # print(all1_tmp[which(all1_tmp$MarkerName %in% c(tmp$SNP_A,tmp$SNP_B)),"old_region_signal_name"])
    # add LD:
    all1$variants_in_ld_r2_0.1[i]<-paste(tmp$ld,collapse=";")

  }

}

# rename those signals:
all1$class_signal_postfm[!is.na(all1$variants_in_ld_r2_0.1)]<-"new_fm_credible_set_lead_known_signal"
all1$class_signal_postfm[which(all1$old_region_signal_name!="")]<-"new_fm_credible_set_lead_known_signal"
table(all1$class_signal_postfm)
#  new_fm_credible_set_lead new_fm_credible_set_lead_known_signal 
#                        76                                   192 

all1$class_signal_postfm[which(all1$class_signal_postfm=="new_fm_credible_set_lead")]<-"new_fm_credible_set_lead_new_signal"
table(all1$class_signal_postfm)
# new_fm_credible_set_lead_known_signal   new_fm_credible_set_lead_new_signal 
#                                   192                                    76

table(all1$class_signal_postfm,all1$class)
  #                                       new old
  # new_fm_credible_set_lead_known_signal   0 192
  # new_fm_credible_set_lead_new_signal    21  55

table(all1$class_signal_postfm,all1$class_signal)

length(names(table(all1$updated_region)))
# [1] 170


# how many signals we have on those regions from unsupervised conditional analyses?
### for comparison purposes - how many signals would have in the FM regions pre-fm:

tmp<-all[which(all$updated_region %in% all1$updated_region),]
tmp<-tmp[grep("new_cojo_unsupervised",tmp$class_signal),]
tmp<-tmp[!grepl("new_cojo_unsupervised_in_ld_new_cojo_unsupervised",tmp$class_signal),]
table(tmp$class_signal[grep("new_cojo_unsupervised",tmp$class_signal) & !grepl("new_cojo_unsupervised_in_ld_with_old",tmp$class_signal)])

nrow(tmp[grep("new_cojo_unsupervised",tmp$class_signal) & !grepl("new_cojo_unsupervised_in_ld_with_old",tmp$class_signal),])
# [1] 291


#########################################
# REGIONS THAT HAVE not BEEN FINE MAPPED:

# for the remaining regions - retain any signal from new_cojo_unsupervised
all2<-all[which(!all$updated_region %in% all1$updated_region),]
length(names(table(all2$updated_region)))
# [1] 250

table(all2$class_signal[which(all2$class!="old")])
#                             new_cojo_unsupervised  # report as new
#                                               130 
# new_cojo_unsupervised_in_ld_new_cojo_unsupervised  # do not report
#                                                33 
#              new_cojo_unsupervised_in_ld_with_old  # report as known
#                                                 1 
#             new_signal_new_region_cojo_supervised  # do not report
#                                                21 

table(all2$class_signal[which(all2$class=="old")])

#                new_cojo_supervised_in_ld_with_old   # do not report
#                                                 8 
#                             new_cojo_unsupervised   # report as new
#                                                79 
# new_cojo_unsupervised_in_ld_new_cojo_unsupervised   # do not report
#                                                24 
#              new_cojo_unsupervised_in_ld_with_old   # report as known
#                                                69 
#                                               old   # do not report
#                                               131 
#                         old|new_cojo_unsupervised   # report as known
#                                                 7

all2_tmp<-all2[which(!all2$class_signal %in% c("new_cojo_unsupervised_in_ld_new_cojo_unsupervised","new_cojo_supervised_in_ld_with_old","new_signal_new_region_cojo_supervised","old")),]
dim(all2_tmp)
# [1] 286 214

length(names(table(all2[which(!all2$updated_region %in% all2_tmp$updated_region),"updated_region"])))
# [1] 63 - 63 regions with no signals

# 22 have genome wide signals but do not come from unsupervised FM

table(all2$old_region_signal_name[which(!all2$updated_region %in% all2_tmp$updated_region)])
#                           delange_region142_signal1 delange_region154_signal1 
#                        22                         1                         1 
#  delange_region28_signal1   Huang_region226_signal1        Liu_region268_noFM 
#                         1                         1                         1 
#        Liu_region274_noFM        Liu_region276_noFM        Liu_region277_noFM 
#                         1                         1                         1 
# ....


table(all2$class_signal)
#           in_ld_new_cojo_unsupervised                        in_ld_with_old 
#                                    57                                    72 
#    new_cojo_supervised_in_ld_with_old                 new_cojo_unsupervised 
#                                     8                                   209 
# new_signal_new_region_cojo_supervised                                   old 
#                                    21                                   131 
#             old|new_cojo_unsupervised 
#                                     5 


all2<-all2[which(all2$class_signal %in% c("new_cojo_unsupervised","old|new_cojo_unsupervised","new_cojo_unsupervised","new_cojo_unsupervised_in_ld_with_old")),]

table(all2$class_signal_postfm,useNA="ifany")
# <NA> 
#  286 

all2$class_signal_postfm[which(all2$class_signal %in% c("old|new_cojo_unsupervised","new_cojo_unsupervised_in_ld_with_old"))]<-"new_cojo_unsupervised_known_signal"
all2$class_signal_postfm[which(all2$class_signal %in% c("new_cojo_unsupervised"))]<-"new_cojo_unsupervised_new_signal"
table(all2$class_signal_postfm)
# new_cojo_unsupervised_known_signal   new_cojo_unsupervised_new_signal 
#                                 77                                209

table(all2$class_signal,all2$class)
  #                                    new old
  # new_cojo_unsupervised                130  79
  # new_cojo_unsupervised_in_ld_with_old   1  69
  # old|new_cojo_unsupervised              0   7


table(all2$class_signal_postfm,all2$class)                           
  #                                    new old
  # new_cojo_unsupervised_known_signal   1  76
  # new_cojo_unsupervised_new_signal   130  79

dim(all2[which(all2$updated_region %in% all1$updated_region),])
# [1]   0 214



#####################################################################################################################################################
# how many regions with signals from each category:


reg<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_regions_420_cojo_supervised_plus_unsupervised.tsv.gz",sep=""))

dim(reg[which(reg$updated_region %in% all2$updated_region),])
# [1] 187   2

dim(reg[which(reg$updated_region %in% all1$updated_region),])
# [1] 170   2

dim(reg[which(!(reg$updated_region %in% all1$updated_region) & !(reg$updated_region %in% all2$updated_region)),])
# [1] 63   2


all3<-all[which(all$updated_region %in% reg$updated_region[which(!(reg$updated_region %in% all1$updated_region) & !(reg$updated_region %in% all2$updated_region))]),]


###############################################################################################
# table to report final signals from our analyses:


all<-rbind(all1,all2)

table(all$class_signal_postfm,useNA="ifany")
#    new_cojo_unsupervised_known_signal      new_cojo_unsupervised_new_signal 
#                                    77                                   209 
# new_fm_credible_set_lead_known_signal   new_fm_credible_set_lead_new_signal 
#                                   192                                    76 

table(all$class,all$region_fm_2025,useNA="ifany")
  #     N_filtered_out_postFM N_filtered_out_preFM   Y
  # new                     6                  124  21
  # old                    22                   62 247


length(names(table(all$updated_region)))
# [1] 357

length(names(table(all$updated_region[which(all$region_fm_2025=="Y")])))
# [1] 170
length(names(table(all$updated_region[which(all$region_fm_2025 %in% c("N_filtered_out_postFM","N_filtered_out_preFM"))])))
# [1] 187

length(names(table(all$updated_region[which(all$region_fm_2025 %in% c("N_filtered_out_postFM"))])))
# [1] 17
length(names(table(all$updated_region[which(all$region_fm_2025 %in% c("N_filtered_out_preFM"))])))
# [1] 170

dim(reg[which(!reg$updated_region %in% all$updated_region)])
# [1] 63   2

# fm file includes regions later excluded (post FM) - update accordingly:
fm_reg<-fread("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/fine_mapping/fine_mapping_region_reclassification.csv")
dim(fm_reg)
fm_reg$region_label[which(fm_reg$region_label=="2_229716690_230796312")]<-"2_229716690_231482974"
table(fm_reg$reason_noCS)
  #                                .              classA_no_significant 
  #                              170                                  1 
  #                    classB_purity classC1_nosignificant_after_filter 
  #                                1                                231 
  # classC2_significant_after_filter 
  #                               16 


table(fm_reg$reason_noCS[!fm_reg$region_label %in% all$updated_region])
  #            classA_no_significant classC1_nosignificant_after_filter 
  #                                1                                106 
  # classC2_significant_after_filter 
  #                                1 

dim(all[which(all$region_fm_2025=="Y"),])
# [1] 268 213

dim(all[which(all$region_fm_2025!="Y"),])
# [1] 286 213

table(all$class_signal_postfm,useNA="ifany")
#    new_cojo_unsupervised_known_signal      new_cojo_unsupervised_new_signal 
#                                    77                                   209 
# new_fm_credible_set_lead_known_signal   new_fm_credible_set_lead_new_signal 
#                                   192                                    76 

all$chr<-gsub("chr","",all$chr)

table(all$phenotype,useNA="ifany")
            #  CD   IBD_saturated IBD_unsaturated              UC            <NA> 
            # 100             143             220              76              15 

all[which(is.na(all$phenotype)),c("MarkerName","class_signal")]
#                 MarkerName
# 159      chr1:66934185:C:T - model does not converge
# 451    chr14:105753554:C:T - not in both gsa and hce
# 487     chr15:44802938:T:C - not in both gsa and hce
# 559     chr16:50537867:C:T - very low freq in gsa or hce
# 1204      chr5:4422049:C:T - not in both gsa and hce
# 216    chr10:125821742:T:G - not in both gsa and hce
# 306    chr11:117998788:C:T - very low freq in gsa or hce
# 483     chr15:34731239:A:C - not in both gsa and hce
# 816     chr2:18792705:C:CT - very low freq in gsa or hce
# 1097     chr4:42396369:G:T
# 1101 chr4:42635348:ATATC:A
# 1364      chr7:4827597:G:C
# 1365      chr7:4881854:T:C
# 1367      chr7:4981470:C:G
# 1377      chr7:5139807:G:A
#                                                        class_signal
# 159                  new_cojo_unsupervised|new_fm_credible_set_lead
# 451                                        new_fm_credible_set_lead
# 487  new_signal_new_region_cojo_supervised|new_fm_credible_set_lead
# 559                                        new_fm_credible_set_lead
# 1204 new_signal_new_region_cojo_supervised|new_fm_credible_set_lead
# 216                                           new_cojo_unsupervised
# 306                                           new_cojo_unsupervised
# 483                                           new_cojo_unsupervised
# 816                                           new_cojo_unsupervised
# 1097                                          new_cojo_unsupervised
# 1101                                          new_cojo_unsupervised
# 1364                                          new_cojo_unsupervised
# 1365                                          new_cojo_unsupervised
# 1367                                          new_cojo_unsupervised
# 1377                                          new_cojo_unsupervised

fwrite(all,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_with_fm_signals_no_overlap.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")

all_keep<-all

############################################
# annotate the remaining regions

all<-all3

table(all$class_signal)
all$class_signal[which(all$class_signal=="new_signal_new_region_cojo_supervised")]<-"new_cojo_supervised"
table(all$class_signal)


all$"P-value_cd_eur_tier_1"<-as.numeric(all$"P-value_cd_eur_tier_1")
all$"P-value_uc_eur_tier_1"<-as.numeric(all$"P-value_uc_eur_tier_1")
all$"P-value_ibd_eur_tier_1"<-as.numeric(all$"P-value_ibd_eur_tier_1")


all$"P-value_cd_eur_tier_2"<-as.numeric(all$"P-value_cd_eur_tier_2")
all$"P-value_uc_eur_tier_2"<-as.numeric(all$"P-value_uc_eur_tier_2")
all$"P-value_ibd_eur_tier_2"<-as.numeric(all$"P-value_ibd_eur_tier_2")


all$"P-value_cd_eur_tier_2_eas_sas"<-as.numeric(all$"P-value_cd_eur_tier_2_eas_sas")
all$"P-value_uc_eur_tier_2_eas_sas"<-as.numeric(all$"P-value_uc_eur_tier_2_eas_sas")
all$"P-value_ibd_eur_tier_2_eas_sas"<-as.numeric(all$"P-value_ibd_eur_tier_2_eas_sas")

all$"P-value_association_cd_eur_tier_2_eas_sas_mrmega"<-as.numeric(all$"P-value_association_cd_eur_tier_2_eas_sas_mrmega")
all$"P-value_association_uc_eur_tier_2_eas_sas_mrmega"<-as.numeric(all$"P-value_association_uc_eur_tier_2_eas_sas_mrmega")
all$"P-value_association_ibd_eur_tier_2_eas_sas_mrmega"<-as.numeric(all$"P-value_association_ibd_eur_tier_2_eas_sas_mrmega")


all$genomewide_significance_eur_tier_1<-"N"
all$genomewide_significance_eur_tier_1[which(all$"P-value_cd_eur_tier_1"<5E-8 | all$"P-value_uc_eur_tier_1"<5E-8 | all$"P-value_ibd_eur_tier_1"<5E-8)]<-"Y"
table(all$genomewide_significance_eur_tier_1)
#   N 
#  64

all$genomewide_significance_eur_tier_2<-"N"
all$genomewide_significance_eur_tier_2[which(all$"P-value_cd_eur_tier_2"<5E-8 | all$"P-value_uc_eur_tier_2"<5E-8 | all$"P-value_ibd_eur_tier_2"<5E-8)]<-"Y"
table(all$genomewide_significance_eur_tier_2)
#  N  Y 
# 61  3 

all$genomewide_significance_eur_eas_sas_tier_2<-"N"
all$genomewide_significance_eur_eas_sas_tier_2[which(all$"P-value_cd_eur_tier_2_eas_sas"<5E-8 | all$"P-value_uc_eur_tier_2_eas_sas"<5E-8 | all$"P-value_ibd_eur_tier_2_eas_sas"<5E-8 |
all$"P-value_association_cd_eur_tier_2_eas_sas_mrmega"<5E-8 | all$"P-value_association_uc_eur_tier_2_eas_sas_mrmega"<5E-8 | all$"P-value_association_ibd_eur_tier_2_eas_sas_mrmega"<5E-8)]<-"Y"
table(all$genomewide_significance_eur_eas_sas_tier_2)
#  N  Y 
# 38 26  

table(all$genomewide_significance_eur_tier_2,all$genomewide_significance_eur_tier_1,useNA="ifany")
  #    N
  # N 61
  # Y  3


table(all$genomewide_significance_eur_tier_2,all$genomewide_significance_eur_eas_sas_tier_2,useNA="ifany")
  #    N  Y
  # N 38 23
  # Y  0  3

### create one flag per region:

all$region_genomewide_significance_eur_tier_1<-"N"
all$region_genomewide_significance_eur_tier_2<-"N"
all$region_genomewide_significance_eur_eas_sas_tier_2<-"N"

regions<-all$updated_region
regions<-regions[!duplicated(regions)]
length(regions)
# [1] 63

for (i in 1:length(regions)) {

  tmp<-all[which(all$updated_region==regions[i]),]

  if (nrow(tmp[which(tmp$genomewide_significance_eur_tier_1=="Y"),])>0) {
    all$region_genomewide_significance_eur_tier_1[which(all$updated_region==regions[i])]<-"Y"
  } 

  if (nrow(tmp[which(tmp$genomewide_significance_eur_tier_2=="Y"),])>0) {
    all$region_genomewide_significance_eur_tier_2[which(all$updated_region==regions[i])]<-"Y"
  }

  if (nrow(tmp[which(tmp$genomewide_significance_eur_eas_sas_tier_2=="Y"),])>0) {
    all$region_genomewide_significance_eur_eas_sas_tier_2[which(all$updated_region==regions[i])]<-"Y"
  }

}

# regions with significant variants in eur tier 1
length(names(table(all$updated_region[which(all$genomewide_significance_eur_tier_1=="Y")])))
# [1] 0

# regions with significant variants in eur tier 2
length(names(table(all$updated_region[which(all$genomewide_significance_eur_tier_2=="Y")])))
# [1] 3

# regions with significant variants in eur_eas_sas tier 2
length(names(table(all$updated_region[which(all$region_genomewide_significance_eur_eas_sas_tier_2=="Y")])))
# [1] 25

all$class_signal

fwrite(all,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_new_and_old_index_gwas_variants_regions_no_conditional_significant_eur_tier2.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")



######################################################
# recapture those that reach genome wide singificance in multiancestry data

all_recap<-all[which(all$genomewide_significance_eur_eas_sas_tier_2=="Y"),colnames(all_keep)]
table(all_recap$class_signal)
              #  new_cojo_supervised new_cojo_supervised_in_ld_with_old 
              #                   19                                  1 
              #                  old 
              #                    6 
# remove new_cojo_supervised_in_ld_with_old
all_recap<-all_recap[which(all_recap$class_signal!="new_cojo_supervised_in_ld_with_old"),]


all_recap$class_signal_postfm[which(all_recap$class_signal=="old")]<-"new_cojo_supervised_gw_significant_multiancestry_known_signal"
all_recap$class_signal_postfm[which(all_recap$class_signal=="new_cojo_supervised")]<-"new_cojo_supervised_gw_significant_multiancestry_new_signal"
table(all_recap$class_signal_postfm)
# new_cojo_supervised_gw_significant_multiancestry_known_signal 
#                                                             6 
#   new_cojo_supervised_gw_significant_multiancestry_new_signal 
#                                                            19



for (chr in 1:22) {
  ld.tmp<-read.table(paste0(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr",chr,"_list_tier_2_unsupervised_conditional_variants_for_ld_with_fm_with_old_data.ld"),head=T)
  ld.tmp<-ld.tmp[which(ld.tmp$R2>=0.1),]

  if(chr==1) {
    ld<-ld.tmp
  } else {
    ld<-rbind(ld,ld.tmp)
  }
  rm(ld.tmp)
}


ld<-ld[which( (ld$SNP_A %in% all_recap$MarkerName) | (ld$SNP_B %in% all_recap$MarkerNamed)),]
# 0





all_recap<-rbind(all_recap,all_keep)

table(all_recap$class_signal)

#                                                        new_cojo_supervised 
#                                                                         19 
#                new_cojo_supervised_in_ld_with_old|new_fm_credible_set_lead 
#                                                                          1 
#                                                      new_cojo_unsupervised 
#                                                                        209 
# new_cojo_unsupervised_in_ld_new_cojo_unsupervised|new_fm_credible_set_lead 
#                                                                          4 
#                                       new_cojo_unsupervised_in_ld_with_old 
#                                                                         70 
#              new_cojo_unsupervised_in_ld_with_old|new_fm_credible_set_lead 
#                                                                         43 
#                             new_cojo_unsupervised|new_fm_credible_set_lead 
#                                                                         35 
#                                                   new_fm_credible_set_lead 
#                                                                        143 
#             new_signal_new_region_cojo_supervised|new_fm_credible_set_lead 
#                                                                          2 
#                                                                        old 
#                                                                          6 
#                                                  old|new_cojo_unsupervised 
#                                                                          7 
#                         old|new_cojo_unsupervised|new_fm_credible_set_lead 
#                                                                         30 
#                                               old|new_fm_credible_set_lead 
#                                                                         10 


table(all_recap$class_signal_postfm,all_recap$class,useNA="ifany") 
  #                                                               new old
  # new_cojo_supervised_gw_significant_multiancestry_known_signal   0   6
  # new_cojo_supervised_gw_significant_multiancestry_new_signal    19   0
  # new_cojo_unsupervised_known_signal                              1  76
  # new_cojo_unsupervised_new_signal                              130  79
  # new_fm_credible_set_lead_known_signal                           0 192
  # new_fm_credible_set_lead_new_signal                            21  55


# double check that no singal is in LD with each other:

fwrite(all_recap,paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_with_fm_signals_with_multiancestry_signals.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")

q("no")


###################################################################################################################################
# out of all the old index variants - how many are represented (tag) at different r2 thresholds by the new index variants:


# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=8000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

old<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_old_regions_and_signals.tsv.gz"))
old<-as.data.frame(old)
dim(old)
# [1] 378 160

new<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_with_fm_signals_no_overlap.tsv.gz"))
new<-as.data.frame(new)
dim(new)
# [1] 554 211

ids<-c(old$MarkerName,new$MarkerName)
ids<-ids[!duplicated(ids)]

for (chr in 1:22) {
  ld.tmp<-read.table(paste0(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr",chr,"_list_tier_2_unsupervised_conditional_variants_for_ld_with_fm_with_old_data.ld"),head=T)
  ld.tmp<-ld.tmp[which(ld.tmp$R2>=0.1),]

  if(chr==1) {
    ld<-ld.tmp
  } else {
    ld<-rbind(ld,ld.tmp)
  }
  rm(ld.tmp)
}


ld<-ld[which( (ld$SNP_A %in% ids) & (ld$SNP_B %in% ids)),]

old$final_index_variants_in_ld_r2_0.1<-NA
old$final_index_variants_in_ld_r2_0.2<-NA
old$final_index_variants_in_ld_r2_0.4<-NA
old$final_index_variants_in_ld_r2_0.5<-NA
old$final_index_variants_in_ld_r2_0.6<-NA
old$final_index_variants_in_ld_r2_0.7<-NA
old$final_index_variants_in_ld_r2_0.8<-NA
old$final_index_variants_in_ld_r2_0.9<-NA

r2_thresholds<-c(0.1,0.2,0.4,0.5,0.6,0.7,0.8,0.9)

for (i in 1:nrow(old)) {

  if (old$MarkerName[i] %in% new$MarkerName) {
    for (r2 in r2_thresholds) {
      old[i,paste0("final_index_variants_in_ld_r2_",r2)]<-paste(as.character(old$MarkerName[i]),as.character(old$MarkerName[i]),1,sep="|")
    }
  } else {
      tmp<-ld[which(ld$SNP_A %in% old$MarkerName[i] | (ld$SNP_B %in% old$MarkerName[i])),]
    tmp<-tmp[which(tmp$SNP_A %in% new$MarkerName | (tmp$SNP_B %in% new$MarkerName)),]

    tmp$ld<-paste(as.character(tmp$SNP_A),as.character(tmp$SNP_B),tmp$R2,sep="|")

    if (nrow(tmp)>0) {

      for (r2 in r2_thresholds) {
        old[i,paste0("final_index_variants_in_ld_r2_",r2)]<-paste(tmp$ld[which(tmp$R2>=r2)],collapse=";")
      }

    }


  }

}

old[old==""]<-NA

dim(old)
# [1] 378 167

dim(old[which(is.na(old$final_index_variants_in_ld_r2_0.1)),])
# [1] 132 167
dim(old[which(is.na(old$final_index_variants_in_ld_r2_0.5)),])
# [1] 177 167
dim(old[which(is.na(old$final_index_variants_in_ld_r2_0.6)),])
# [1] 195 167
dim(old[which(is.na(old$final_index_variants_in_ld_r2_0.7)),])
# [1] 209 167

table(old$'FM study')
    # Cordero Cordero|Liu     delange       Huang         Liu        noFM 
    #       9           1          66         139          81          80 
    #  Sakaue 
    #       2 

table(old$'FM study'[which(is.na(old$final_index_variants_in_ld_r2_0.7))])
    # Cordero Cordero|Liu     delange       Huang         Liu        noFM 
    #       8           1          27          54          61          56 
    #  Sakaue 
    #       2 

old<-old[,c("MarkerName","old_region_signal_name","FM study","updated_region",
"final_index_variants_in_ld_r2_0.1","final_index_variants_in_ld_r2_0.2","final_index_variants_in_ld_r2_0.4","final_index_variants_in_ld_r2_0.5","final_index_variants_in_ld_r2_0.6","final_index_variants_in_ld_r2_0.7","final_index_variants_in_ld_r2_0.8","final_index_variants_in_ld_r2_0.9",
"P-value_ibd_eur_tier_2","P-value_cd_eur_tier_2","P-value_uc_eur_tier_2")]

old<-merge(old,new[,c("MarkerName","class_signal")],by="MarkerName",all.x=T)


## add fine-mapping flag:

# fm file includes regions later excluded (post FM) - update accordingly:
fm_reg<-fread("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/fine_mapping/fine_mapping_region_reclassification.csv")
dim(fm_reg)
fm_reg$region_label[which(fm_reg$region_label=="2_229716690_230796312")]<-"2_229716690_231482974"
table(fm_reg$reason_noCS)

old$region_fm_2025<-NA
old$region_fm_2025[which(old$updated_region %in% fm_reg$region_label[which(fm_reg$reason_noCS %in% c("classC1_nosignificant_after_filter","classA_no_significant"))])]<-"N_filtered_out_preFM"
old$region_fm_2025[which(old$updated_region %in% fm_reg$region_label[which(fm_reg$reason_noCS %in% c("classB_purity","classC2_significant_after_filter"))])]<-"N_filtered_out_postFM"
old$region_fm_2025[which(old$updated_region %in% fm_reg$region_label[which(fm_reg$reason_noCS %in% c("."))])]<-"Y"

table(old$class_signal,old$region_fm_2025,useNA="ifany")
  #                                                    N_filtered_out_postFM
  # old|new_cojo_unsupervised                                              1
  # old|new_cojo_unsupervised|new_fm_credible_set_lead                     0
  # old|new_fm_credible_set_lead                                           0
  # <NA>                                                                  13
                                                    
  #                                                    N_filtered_out_preFM   Y
  # old|new_cojo_unsupervised                                             6   0
  # old|new_cojo_unsupervised|new_fm_credible_set_lead                    0  28
  # old|new_fm_credible_set_lead                                          0  10
  # <NA>                                                                118 202


old$tagged_at_r2_0.7<-"N"
old$tagged_at_r2_0.7[which(!is.na(old$final_index_variants_in_ld_r2_0.7))]<-"Y"
table(old$tagged_at_r2_0.7)
#   N   Y 
# 209 169 

old$tagged_at_r2_0.1<-"N"
old$tagged_at_r2_0.1[which(!is.na(old$final_index_variants_in_ld_r2_0.1))]<-"Y"
table(old$tagged_at_r2_0.1)
#   N   Y 
# 132 246 

old$tagged_at_r2_0.5<-"N"
old$tagged_at_r2_0.5[which(!is.na(old$final_index_variants_in_ld_r2_0.5))]<-"Y"
table(old$tagged_at_r2_0.5)
#   N   Y 
# 177 201

table(old$region_fm_2025)
# N_filtered_out_postFM  N_filtered_out_preFM                     Y 
                  #  14                   124                   240

table(old$tagged_at_r2_0.7,old$region_fm_2025)
  #   N_filtered_out_postFM N_filtered_out_preFM   Y
  # N                     6                   98 105
  # Y                     8                   26 135

table(old$tagged_at_r2_0.1,old$region_fm_2025)
  #   N_filtered_out_postFM N_filtered_out_preFM   Y
  # N                     1                   71  60
  # Y                    13                   53 180


table(old$region_fm_2025,old$tagged_at_r2_0.1)



# how many regions, how many fine mapped:
tmp<-old[,c("updated_region","region_fm_2025")]
tmp<-tmp[!duplicated(tmp),]

dim(tmp)
# [1] 276   2
length(names(table(old$updated_region)))
# [1] 276

table(old$'FM study',old$tagged_at_r2_0.1)
  #               N   Y
  # Cordero       5   4
  # Cordero|Liu   1   0
  # delange      14  52
  # Huang        24 115
  # Liu          49  32
  # noFM         37  43
  # Sakaue        2   0

table(old$'FM study',old$tagged_at_r2_0.7)
  #              N  Y
  # Cordero      8  1
  # Cordero|Liu  1  0
  # delange     27 39
  # Huang       54 85
  # Liu         61 20
  # noFM        56 24
  # Sakaue       2  0


table(old$'FM study',old$tagged_at_r2_0.7)[,1]/table(old$'FM study')
table(old$'FM study',old$tagged_at_r2_0.1)[,1]/table(old$'FM study')

# all decrease proportionally ~0.5

# are the main signals replicated?
tmp1<-old[grep("_signal1",old$old_region_signal_name),]
table(tmp1$'FM study',tmp1$tagged_at_r2_0.7)



# none of the regions where we do not nominate independent signals shoudl have a variant in ld:
ns<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_new_and_old_index_gwas_variants_regions_no_conditional_significant_eur_tier2.tsv.gz",sep=""))
old[which(old$updated_region %in% ns$updated_region),]
# NA - OK


# add the summary stats from de Lange:
pheno<-c("cd","uc","ibd")

for (ph in pheno) {

  # delange:

  tmp<-fread(paste0(path_gwas,"summary_files/summary_stats_",ph,"_delange_sorted_noheader_chr_pos_b37_pval_lifted_hg38.gz",sep=""),head=F)
  tmp<-tmp[,1:5]
  colnames(tmp)<-c("chr","position_b38","end","ID_b37","pvalue")
  tmp<-tmp[,c("chr","position_b38","ID_b37","pvalue")]

  tmp<-tmp[,c("ID_b37","chr","position_b38","pvalue")]

  colnames(tmp)[ncol(tmp)]<-paste(colnames(tmp)[ncol(tmp)],"_delange_",ph,sep="")

  if (!exists("old_delange")) {
    old_delange<-tmp
  } else {
    old_delange<-merge(old_delange,tmp[,c(1,4)],by="ID_b37",all=T)
  }
}

old_delange$ref<-gsub("[0-9]{1,2}:[0-9]*_","",old_delange$ID_b37)
old_delange$alt<-gsub(".*_","",old_delange$ref)
old_delange$ref<-gsub("_.*","",old_delange$ref)

old_delange$MarkerName_1<-paste(old_delange$chr,old_delange$position_b38,old_delange$ref,old_delange$alt,sep=":")
old_delange$MarkerName_2<-paste(old_delange$chr,old_delange$position_b38,old_delange$alt,old_delange$ref,sep=":")

old_delange1<-old_delange[which(old_delange$MarkerName_1 %in% old$MarkerName),]
old_delange2<-old_delange[which(old_delange$MarkerName_2 %in% old$MarkerName),]

dim(old_delange1)
# [1] 346  10
dim(old_delange2)
# [1]  1 10

old_delange1<-old_delange1[,c("ID_b37","chr","position_b38","pvalue_delange_cd","pvalue_delange_uc","pvalue_delange_ibd","ref","alt","MarkerName_1")]
old_delange2<-old_delange2[,c("ID_b37","chr","position_b38","pvalue_delange_cd","pvalue_delange_uc","pvalue_delange_ibd","ref","alt","MarkerName_2")]

colnames(old_delange1)<-c("ID_b37","chr","position_b38","pvalue_delange_cd","pvalue_delange_uc","pvalue_delange_ibd","ref","alt","MarkerName")
colnames(old_delange2)<-c("ID_b37","chr","position_b38","pvalue_delange_cd","pvalue_delange_uc","pvalue_delange_ibd","ref","alt","MarkerName")

old_delange<-rbind(old_delange1,old_delange2)

old<-merge(old,old_delange[,c("MarkerName","pvalue_delange_cd","pvalue_delange_uc","pvalue_delange_ibd")],by="MarkerName",all.x=T)

old$gwas_signif_delange<-"N"
old$gwas_signif_delange[which(old$pvalue_delange_cd<=5E-8 | old$pvalue_delange_ibd<=5E-8 | old$pvalue_delange_uc<=5E-8)]<-"Y"

table(old$tagged_at_r2_0.7,old$gwas_signif_delange)
  #     N   Y
  # N 154  55
  # Y  61 108

table(old$tagged_at_r2_0.7[which(old$'FM study'=="delange")],old$gwas_signif_delange[which(old$'FM study'=="delange")])
  #    N  Y
  # N  0 27
  # Y  1 38



# add the summary stats from de Liu (EUR):
pheno<-c("cd","uc","ibd")

rm(old_liu)
for (ph in pheno) {

  # liu liu gao:

  tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/eas_iibdgc/liu-2022-east-asian-gwas/summary-stats/ibd_EAS_EUR_SiKJEF_meta_",toupper(ph),".TBL.txt.gz",sep=""),head=T)
  tmp<-tmp[,c("MarkerName","Allele1","Allele2","CHR","BP","P-value")]
  tmp<-tmp[which(tmp$MarkerName %in% old$MarkerName),]
  
  tmp<-tmp[,c("MarkerName","CHR","BP","P-value")]
  colnames(tmp)<-c("MarkerName","chr","position_b38","pvalue")
  
  colnames(tmp)[ncol(tmp)]<-paste(colnames(tmp)[ncol(tmp)],"_liu_",ph,sep="")

  if (!exists("old_liu")) {
    old_liu<-tmp
  } else {
    old_liu<-merge(old_liu,tmp[,c(1,4)],by="MarkerName",all=T)
  }
}

dim(old_liu)
# [1] 327   6

old<-merge(old,old_liu[,c("MarkerName","pvalue_liu_cd","pvalue_liu_uc","pvalue_liu_ibd")],by="MarkerName",all.x=T)


old$gwas_signif_liu<-"N"
old$gwas_signif_liu[which(old$pvalue_liu_cd<=5E-8 | old$pvalue_liu_uc<=5E-8 | old$pvalue_liu_ibd<=5E-8)]<-"Y"

table(old$gwas_signif_liu)
#   N   Y 
# 127 251 

fwrite(old,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_old_regions_and_signals_with_ld_with_final_new.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")

q("no")

#######################################################################

# investigate not tagged variants:


# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=10000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
library(ggplot2)
library(ggpubr)

rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"


old<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_old_regions_and_signals_with_ld_with_final_new.tsv.gz",sep=""))
old$min_pval_delange<-pmin(old$pvalue_delange_cd,old$pvalue_delange_uc,old$pvalue_delange_ibd)
old$min_pval_liu<-pmin(old$pvalue_liu_cd,old$pvalue_liu_uc,old$pvalue_liu_ibd)


old$"P-value_ibd_eur_tier_2"<-as.numeric(old$"P-value_ibd_eur_tier_2")
old$"P-value_cd_eur_tier_2"<-as.numeric(old$"P-value_cd_eur_tier_2")
old$"P-value_uc_eur_tier_2"<-as.numeric(old$"P-value_uc_eur_tier_2")

old$min_pval_iibdgc<-as.numeric(pmin(old$"P-value_ibd_eur_tier_2",old$"P-value_cd_eur_tier_2",old$"P-value_uc_eur_tier_2"))

old$gwas_signif_iibdgc<-"N"
old$gwas_signif_iibdgc[which(old$"P-value_ibd_eur_tier_2"<=5E-8 | old$"P-value_cd_eur_tier_2"<=5E-8 | old$"P-value_uc_eur_tier_2"<=5E-8)]<-"Y"

table(old$gwas_signif_iibdgc,old$gwas_signif_delange)
  #     N   Y
  # N 112  24
  # Y 103 139

old$min_pval_iibdgc[which(old$min_pval_iibdgc<=1E-300)]<-1E-300

old$FM_study<-old$'FM study'

ggplot(old, aes(x=-log10(min_pval_delange),y=-log10(min_pval_iibdgc),color=tagged_at_r2_0.5)) + 
  geom_point() + xlim(0,max(-log10(old$min_pval_delange),-log10(old$min_pval_iibdgc),na.rm=T)) +
  ylim(0,max(-log10(old$min_pval_delange),-log10(old$min_pval_iibdgc),na.rm=T)) + facet_grid(gwas_signif_delange ~ gwas_signif_iibdgc)


old$signal<-gsub(".*_signal","",old$old_region_signal_name)
old$signal<-gsub("Liu_region.*","1",old$signal)
old$signal<-gsub("Cordero_region.*","1",old$signal)
old$signal<-gsub("Sakaue_region.*","1",old$signal)
old$signal<-gsub("Sakaue_region.*","1",old$signal)

table(old$tagged_at_r2_0.5,old$signal)

ggplot(old, aes(x=-log10(min_pval_delange),color=tagged_at_r2_0.5)) + 
  geom_density()


table(old$gwas_signif_delange,old$gwas_signif_iibdgc)
      N   Y
  N 112 103
  Y  24 139

table(old$gwas_signif_delange[which(old$tagged_at_r2_0.5=="N")],old$gwas_signif_iibdgc[which(old$tagged_at_r2_0.5=="N")])
  #     N   Y
  # N 101  35
  # Y  23  18


table(old$gwas_signif_delange[which(old$tagged_at_r2_0.5=="Y")],old$gwas_signif_iibdgc[which(old$tagged_at_r2_0.5=="Y")])
  #     N   Y
  # N  11  68
  # Y   1 121


table(old$gwas_signif_delange,old$FM_study)
  #   Cordero Cordero|Liu delange Huang Liu noFM Sakaue
  # N       9           1       1    52  71   79      2
  # Y       0           0      65    87  10    1      0

table(old$FM_study,old$tagged_at_r2_0.5)


table(old$FM_study,old$tagged_at_r2_0.5)[,2]/table(old$FM_study)
  #   Cordero Cordero|Liu     delange       Huang         Liu        noFM 
  # 0.4444444   0.0000000   0.6515152   0.7122302   0.2962963   0.3875000 
  #    Sakaue 
  # 0.0000000 

# exclude those regions where we are not assignig independent GWAS signals:


###############################################################
# De lange:
dl<-old[which(old$FM_study=="delange"),]

table(dl$tagged_at_r2_0.5)
#  N  Y 
# 23 43 

dl[which(dl$tagged_at_r2_0.5=="N"),]

table(dl$tagged_at_r2_0.5,dl$gwas_signif_iibdgc)
  #    N  Y
  # N 15  8
  # Y  1 42


# INDEX VARIANTS no GWAS signif in new analyses - not tagged by new signals
dl_plot<-ggplot(dl, aes(x=-log10(min_pval_delange),y=-log10(min_pval_iibdgc),color=tagged_at_r2_0.5)) + 
  geom_point() + xlim(0,max(-log10(dl$min_pval_delange),-log10(dl$min_pval_iibdgc),na.rm=T)) +
  ylim(0,max(-log10(dl$min_pval_delange),-log10(dl$min_pval_iibdgc),na.rm=T)) + facet_grid(gwas_signif_delange ~ gwas_signif_iibdgc) + 
  ylab("Genome-wide significant IIBDGC 2025 (Y/N)") + xlab("Genome-wide significant de Lange (Y/N)")



###############################################################
# Liu:
li<-old[which(old$FM_study=="Liu"),]
dim(li)
table(li$tagged_at_r2_0.5)
#  N  Y 
# 57 24 


table(li$tagged_at_r2_0.5,li$gwas_signif_iibdgc)
  #    N  Y
  # N 42 15
  # Y  2 22


# INDEX VARIANTS no GWAS signif in new analyses - not tagged by new signals
li_plot<-ggplot(li, aes(x=-log10(min_pval_liu),y=-log10(min_pval_iibdgc),color=tagged_at_r2_0.5)) + 
  geom_point() + xlim(0,max(-log10(li$min_pval_liu),-log10(li$min_pval_iibdgc),na.rm=T)) +
  ylim(0,max(-log10(li$min_pval_liu),-log10(li$min_pval_iibdgc),na.rm=T)) + facet_grid(gwas_signif_liu ~ gwas_signif_iibdgc) + 
  ylab("Genome-wide significant IIBDGC 2025 (Y/N)") + xlab("Genome-wide significant Liu Liu Gao (Y/N)")


p<-ggarrange(dl_plot,li_plot,ncol=2,common.legend=T)

ggsave(paste0("~/git/IIBDGC_GWAS/plots/replicated_non_replicated_ibd_known_signals.pdf",sep=""),
    p,
    width = 14,
    height = 7,
    dpi = 300,
    units = c("in"),
    limitsize = T
)

