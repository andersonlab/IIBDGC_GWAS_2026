# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
##########################################################################################################
#  1.4 EXCLUDE RELATED SAMPLES GENOTYPED WITH DIFFERENT GENOTYPING ARRAYS (WORK JUST WITH EUR SAMPLES):


# /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

array<-c("illumina370","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")

illumina370<-c("belgium_inf1_old_gwas","belgium_inf2_old_gwas","swedish_uc_old_gwas","niddk_old_gwas")
affymetrix6<-c("german_affy6_old_gwas","norway_affy6_old_gwas","gwas2")
humanomniexpress<-c("australia_omniexome")
affymetrix500<-c("gwas1")
humancoreexome<-c("all_hce")
humanomni1<-c("pittsburgh_gsa")
quad610<-c("spain_gsa")
gsa<-c("italy_gsa","kiel_austria_sibdcs_gsa"
       ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
       ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
       ,"mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa","franke_gsa","hyams_protect_gsa",
       "lewis_sparc_gsa","mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
       "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
       "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
       "xavier_share_gsa")
illuminaexome<-c("prism_nfe_gwas","helmsley_prism_gsa","helmsley_xavier_prism_gsa")

studies<-c("all_hce","niddk_old_gwas","australia_omniexome","gwas1","gwas2","pittsburgh_gsa","spain_gsa","italy_gsa",
           "kiel_austria_sibdcs_gsa","netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa",
           "niddk_feinstein_gsa","basque_gsa","lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa",
           "belgium_vermeire_gsa","prism_nfe_gsa","prism_nfe_gwas","finland_illugwas","german_affy6_old_gwas",
           "norway_affy6_old_gwas","belgium_inf1_old_gwas","belgium_inf2_old_gwas",
           "cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas",
           "swedish_uc_old_gwas",
           "mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")


ii<-1
ancestry<-"eur"

for (i in 1:length(array)) {
  
  tmp<-read.table(paste(path,"post_imputation/2022/",array[i],"/genotyped_data/",array[i],"_all_studies_merged_",ancestry[ii],".fam",sep=""),head=F)
  
  tmp$array<-as.character(array[i])
  tmp$V1<-as.character(tmp$V1)
  tmp$V2<-as.character(tmp$V2)
  
  if(i==1){
    dat<-tmp
  }else{
    dat<-rbind(dat,tmp)
  }
}

table(dat$array)
# affymetrix500      affymetrix6              gsa   humancoreexome 
# 4652            11087            81334            22588 
# humanomni1 humanomniexpress      illumina370    illuminaexome 
# 2701             1245             5633             2557 
# quad610 
# 3396 



# keep only nonDuplicated samples:
for (j in 1:length(studies)) {
  
  fam_tmp<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy_interstudy.fam",sep=""),head=F)
  
  fam_tmp$study<-studies[j]
  if(j==1) {
    fam<-fam_tmp
  }else{
    fam<-rbind(fam,fam_tmp)
  }
  
}

dim(fam)
dim(fam[which(fam$V1 %in% dat$V1),])
# [1] 107079      7

dim(fam[which(!fam$V1 %in% dat$V1),])
# [1] 7898    7
table(fam$study[which(!fam$V1 %in% dat$V1)])
# non- eur subset


dat<-merge(dat,fam[,c("V2","study")],by="V2")
dim(dat)
# [1] 107079      8

dat<-dat[which(dat$V2 %in% fam$V2),]
dim(dat)
# [1] 107079     7

table(dat$array)
# affymetrix500      affymetrix6              gsa   humancoreexome 
#          4648             8246            61440            21243 
# humanomni1 humanomniexpress      illumina370    illuminaexome 
#       2053             1238             3452             1363 
# quad610 
#    3396 


# keep only the samples in the nonDuplicated list:


# double check all samples included in iibdgc relatedness exercise, see 'evaluate_relatedness_IIBDGC.R'
famall<-read.table(paste(path,"pre_imputation/QC/relatedness/iibdgc_merged_step25.fam",sep=""),head=F)

dim(famall)
# [1] 116872     6

dim(dat)
# [1] 107079     8

dim(dat[which(!dat$V2 %in% famall$V2),])
# [1] 1931     8

table(dat$study[which(!dat$V2 %in% famall$V2)])
# ccfa_gsa 
# 1931  # OK ccfa_gsa completely replaced by lewis_spark_gsa


colnames(dat)[6]<-"pheno"
colnames(dat)[5]<-"sex"

kin<-read.table(paste(path,"pre_imputation/QC/relatedness/iibdgc_merged_step25_plink.kin0",sep=""),head=F)
colnames(kin)<-c("FID1","IID1","FID2","IID2","NSNP","HETHET","IBS0","KINSHIP")

kin<-kin[which( (kin$IID1 %in% dat$V2) & (kin$IID2 %in% dat$V2)),]

table(cut(kin$KINSHIP,breaks=c(1,0.354,0.177)))
# (0.177,0.354]     (0.354,1]
#          3610             1 

# set the first relatedness threshold based on distribution (see 'evaluate_relatedness_IIBDGC.R')
system(paste("cp ",path,"pre_imputation/QC/relatedness/histogram_iibdgc_merged_plink.pdf ~/tmp_plots/",sep=""))

rel<-merge(kin,dat[,c("V2","array","sex","pheno")],by.x="FID2",by.y="V2",all.x=T,sort=F)
colnames(rel)[(ncol(rel)-3):ncol(rel)]<-paste(colnames(rel)[(ncol(rel)-3):ncol(rel)],"FID2",sep="_")

rel<-merge(rel,dat[,c("V2","array","sex","pheno")],by.x="FID1",by.y="V2",all.x=T,sort=F)
colnames(rel)[(ncol(rel)-3):ncol(rel)]<-paste(colnames(rel)[(ncol(rel)-3):ncol(rel)],"FID1",sep="_")

# exclude related samples in same array:

dim(rel)
# [1] 3611   14
rel<-rel[which(rel$array_FID1!=rel$array_FID2),]
dim(rel)
# [1] 274  14

table(rel$array_FID1,rel$array_FID2)

#                affymetrix500 affymetrix6 gsa humancoreexome humanomni1
# affymetrix500              0          18   1             26          0
# affymetrix6                9           0   6             19          0
# gsa                        0          33   0              0          9
# humancoreexome             0           1   0              0          0
# humanomni1                 0           0   1              0          0
# illumina370                0           0  74              0          0
# illuminaexome              0           0  18              0          0
# 
#                illumina370 illuminaexome quad610
# affymetrix500            0             0       0
# affymetrix6              0             0       0
# gsa                     45            12       1
# humancoreexome           0             0       0
# humanomni1               1             0       0
# illumina370              0             0       0
# illuminaexome            0             0       0




rel_ids<-c(as.character(rel$FID1),as.character(rel$FID2))
length(rel_ids)
#[1] 548

table(length(rel_ids)==(nrow(rel)*2))
# TRUE
# 1

rel_ids<-rel_ids[!duplicated(rel_ids)]
length(rel_ids)
#[1] 523

dim(rel)
# [1] 274    14

rm(data_remove)
for (iii in 1:length(rel_ids)) {
  
  tmp1<-rel[which(rel$FID1==rel_ids[iii]),]
  tmp2<-rel[which(rel$FID2==rel_ids[iii]),]
  colnames(tmp2)[1:2]<-colnames(tmp2)[2:1]
  
  tmp<-rbind(tmp1,tmp2)
  
  # if(nrow(tmp)>1) {
  #   print(iii)
  # }
  
  ids_tmp<-c(as.character(tmp$FID1),as.character(tmp$FID2))
  ids_tmp<-ids_tmp[!duplicated(ids_tmp)]
  
  #keep same order as in rel_ids and in all
  ids_tmp<-ids_tmp[match(rel_ids[which(rel_ids %in% ids_tmp)],ids_tmp)]
  
  
  #keep same order as in rel_ids and in all
  ids_tmp<-ids_tmp[match(rel_ids[which(rel_ids %in% ids_tmp)],ids_tmp)]
  
  data<-as.data.frame(matrix(ncol=1,nrow=length(ids_tmp)))
  data$V1<-ids_tmp
  data<-merge(data,dat[,c("V1","sex","pheno","array")],by="V1",all.x=T,sort=F)
  # data$group<-iii
  
  
  if (any(data$array %in% c("gsa"))){
    
    data<-data[which(data$array!="gsa"),]
    
    if (!exists("data_remove")) {
      data_remove<-data
    } else {
      data_remove<-rbind(data_remove,data)
    }
    
  } else if (!any(data$array %in% c("gsa")) & any(data$array %in% c("affymetrix6"))) {
    
    data<-data[which(data$array!="affymetrix6"),]
    
    if (!exists("data_remove")) {
      data_remove<-data
    } else {
      data_remove<-rbind(data_remove,data)
    }
    
  } else if (!any(data$array %in% c("gsa")) & !any(data$array %in% c("affymetrix6")) & any(data$array %in% c("humancoreexome"))) {
    
    data<-data[which(data$array!="humancoreexome"),]
    
    if (!exists("data_remove")) {
      data_remove<-data
    } else {
      data_remove<-rbind(data_remove,data)
    }
    
  } else if (!any(data$array %in% c("gsa")) & !any(data$array %in% c("affymetrix6")) & !any(data$array %in% c("humancoreexome")) & any(data$array %in% c("affymetrix500"))) {
    
    data<-data[which(data$array!="affymetrix500"),]
    
    if (!exists("data_remove")) {
      data_remove<-data
    } else {
      data_remove<-rbind(data_remove,data)
    }
    
  } else if (!any(data$array %in% c("gsa")) & !any(data$array %in% c("affymetrix6")) & !any(data$array %in% c("humancoreexome"))
             & !any(data$array %in% c("affymetrix500")) & any(data$array %in% c("humanomni1")) ) {
    
    data<-data[which(data$array!="humanomni1"),]
    
    if (!exists("data_remove")) {
      data_remove<-data
    } else {
      data_remove<-rbind(data_remove,data)
    }
    
  } else if (!any(data$array %in% c("gsa")) & !any(data$array %in% c("affymetrix6")) & !any(data$array %in% c("humancoreexome"))
             & !any(data$array %in% c("affymetrix500")) & !any(data$array %in% c("humanomni1")) & any(data$array %in% c("quad610"))) {
    
    data<-data[which(data$array!="quad610"),]
    
    if (!exists("data_remove")) {
      data_remove<-data
    } else {
      data_remove<-rbind(data_remove,data)
    }
    
  } else if (!any(data$array %in% c("gsa")) & !any(data$array %in% c("affymetrix6")) & !any(data$array %in% c("humancoreexome"))
             & !any(data$array %in% c("affymetrix500")) & !any(data$array %in% c("humanomni1")) & !any(data$array %in% c("quad610"))
             & any(data$array %in% c("illumina550"))) {
    
    data<-data[which(data$array!="illumina550"),]
    
    if (!exists("data_remove")) {
      data_remove<-data
    } else {
      data_remove<-rbind(data_remove,data)
    }
    
  } else if (!any(data$array %in% c("gsa")) & !any(data$array %in% c("affymetrix6")) & !any(data$array %in% c("humancoreexome"))
             & !any(data$array %in% c("affymetrix500")) & !any(data$array %in% c("humanomni1")) & !any(data$array %in% c("quad610"))
             & !any(data$array %in% c("illumina550")) & any(data$array %in% c("illumina370")) ) {
    
    data<-data[which(data$array!="illumina370"),]
    
    if (!exists("data_remove")) {
      data_remove<-data
    } else {
      data_remove<-rbind(data_remove,data)
    }
    
  }  else if (!any(data$array %in% c("gsa")) & !any(data$array %in% c("affymetrix6")) & !any(data$array %in% c("humancoreexome"))
              & !any(data$array %in% c("affymetrix500")) & !any(data$array %in% c("humanomni1")) & !any(data$array %in% c("quad610"))
              & !any(data$array %in% c("illumina550")) & !any(data$array %in% c("illumina370")) & any(data$array %in% c("illuminaexome"))) {
    
    data<-data[which(data$array!="illuminaexome"),]
    
    if (!exists("data_remove")) {
      data_remove<-data
    } else {
      data_remove<-rbind(data_remove,data)
    }
    
  }
  
}

dim(data_remove)
#[1] 528   5
# create phenotype files

data_remove<-data_remove[!duplicated(data_remove$V1),]
dim(data_remove)
# [1] 254   4

table(data_remove$array)
# affymetrix500    affymetrix6 humancoreexome     humanomni1    illumina370 
# 54             38             20             10            101 
# illuminaexome        quad610 
# 30              1


data_remove$study<-NA

for (i in 1:length(array)){
  
  print(array[i])
  study<-get(array[i])
  
  print(study)
  
  for (j in 1:length(study)) {
    tmp<-read.table(paste(path,"pre_imputation/QC/",study[j],"/",study[j],"_postqc_preimp_",ancestry[ii],".fam",sep=""),head=F)
    data_remove$study[which(data_remove$V1 %in% tmp$V1)]<-study[j]
  }
}

write.table(data_remove,paste(path,"post_imputation/2022/analysis/list_related_samples_across_arrays_toexclude_",ancestry[ii],sep=""),col.names=T,row.names=F,quote=F,sep="\t")

