# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
#### /software/R-4.3.1/bin/R

###################################
# list of variants across UKIBDGC #
###################################

library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("gwas1","gwas2","all_hce")


for (i in 2:length(cohorts)){
  if(i==2) {
    a<-read.table(paste(path,"pre_imputation/QC/",cohorts[1],"/",cohorts[1],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    common<-intersect(a$V2,b$V2)
    print(cohorts[1])
    print(cohorts[i])
    print(length(common))
    rm(a)
  } else {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    common<-intersect(common,b$V2)
    print(cohorts[i])
    print(length(common))
  }
}

# [1] "gwas1"
# [1] "gwas2"
# [1] 403645
# [1] "all_hce"
# [1] 36884

# keep only those in autosomal chromosomes
table(b[which(b$V2 %in% common),"V1"])
#

common_2<-b[which(b$V2 %in% common),]
common_2<-common_2[which(common_2$V1 %in% seq(1:22)),]
dim(common_2)
# [1] 36239     6

write.table(common_2[,2,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_ukibdgc",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

##################################
# list of variants across German #
##################################

library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("kiel_austria_sibdcs_gsa","german_illu550_old_gwas","german_affy6_old_gwas")

for (i in 2:length(cohorts)){
  if(i==2) {
    a<-read.table(paste(path,"pre_imputation/QC/",cohorts[1],"/",cohorts[1],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    common<-intersect(a$V2,b$V2)
    print(cohorts[1])
    print(cohorts[i])
    print(length(common))
    rm(a)
  } else {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    common<-intersect(common,b$V2)
    print(cohorts[i])
    print(length(common))
  }
}

# [1] "kiel_austria_sibdcs_gsa"
# [1] "german_illu550_old_gwas"
# [1] 100154
# [1] "german_affy6_old_gwas"
# [1] 29032

# keep only those in autosomal chromosomes
common_2<-b[which(b$V2 %in% common),]
common_2<-common_2[which(common_2$V1 %in% seq(1:22)),]
dim(common_2)
# [1] 28009    6

write.table(common_2[,2,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_german",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


#################################
# list of variants across NIDDK #
#################################

library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("niddk_broad_gsa","niddk_feinstein_gsa","niddk_old_gwas")


for (i in 2:length(cohorts)){
  if(i==2) {
    a<-read.table(paste(path,"pre_imputation/QC/",cohorts[1],"/",cohorts[1],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    common<-intersect(a$V2,b$V2)
    print(cohorts[1])
    print(cohorts[i])
    print(length(common))
    rm(a)
  } else {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    common<-intersect(common,b$V2)
    print(cohorts[i])
    print(length(common))
  }
}

# [1] "niddk_broad_gsa"
# [1] "niddk_feinstein_gsa"
# [1] 565523
# [1] "niddk_old_gwas"
# [1] 61889


# keep only those in autosomal chromosomes
common_2<-b[which(b$V2 %in% common),]
common_2<-common_2[which(common_2$V1 %in% seq(1:22)),]
dim(common_2)
# [1] [1] 59626    6

write.table(common_2[,2,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_niddk",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#################################
# list of variants across Spain #
#################################

library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("spain_gsa","basque_gsa")


for (i in 2:length(cohorts)){
  if(i==2) {
    a<-read.table(paste(path,"pre_imputation/QC/",cohorts[1],"/",cohorts[1],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    common<-intersect(a$V2,b$V2)
    print(cohorts[1])
    print(cohorts[i])
    print(length(common))
    rm(a)
  } else {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    common<-intersect(common,b$V2)
    print(cohorts[i])
    print(length(common))
  }
}

# [1] "spain_gsa"
# [1] "basque_gsa"
# [1] 117245

# keep only those in autosomal chromosomes
common_2<-b[which(b$V2 %in% common),]
common_2<-common_2[which(common_2$V1 %in% seq(1:22)),]
dim(common_2)
# [1] 106588    6

write.table(common_2[,2,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_spain",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

###################################
# list of variants across Belgium #
###################################

library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("belgium_inf1_old_gwas","belgium_inf2_old_gwas","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa")


for (i in 2:length(cohorts)){
  if(i==2) {
    a<-read.table(paste(path,"pre_imputation/QC/",cohorts[1],"/",cohorts[1],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    common<-intersect(a$V2,b$V2)
    print(cohorts[1])
    print(cohorts[i])
    print(length(common))
    rm(a)
  } else {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    common<-intersect(common,b$V2)
    print(cohorts[i])
    print(length(common))
  }
}
# [1] "belgium_inf1_old_gwas"
# [1] "belgium_inf2_old_gwas"
# [1] 283815
# [1] "belgium_louis_gsa"
# [1] 61008
# [1] "belgium_franchimont_gsa"
# [1] 60689
# [1] "belgium_vermeire_gsa"
# [1] 60652

# keep only those in autosomal chromosomes
common_2<-b[which(b$V2 %in% common),]
common_2<-common_2[which(common_2$V1 %in% seq(1:22)),]
dim(common_2)
# [1] 58696    6

write.table(common_2[,2,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_belgium",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


#################################
# list of variants across PRISM #
#################################

library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("prism_nfe_gsa","prism_nfe_gwas")


for (i in 2:length(cohorts)){
  if(i==2) {
    a<-read.table(paste(path,"pre_imputation/QC/",cohorts[1],"/",cohorts[1],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    common<-intersect(a$V2,b$V2)
    print(cohorts[1])
    print(cohorts[i])
    print(length(common))
    rm(a)
  } else {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    common<-intersect(common,b$V2)
    print(cohorts[i])
    print(length(common))
  }
}

# [1] "prism_nfe_gsa"
# [1] "prism_nfe_gwas"
# [1] 124092

# keep only those in autosomal chromosomes
common_2<-b[which(b$V2 %in% common),]
common_2<-common_2[which(common_2$V1 %in% seq(1:22)),]
dim(common_2)
# [1] 119802    6

write.table(common_2[,2,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_prism",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


##################################
# list of variants across SWEDEN #
##################################


library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("swedish_uc_old_gwas","sweden_gsa")


for (i in 2:length(cohorts)){
  if(i==2) {
    a<-read.table(paste(path,"pre_imputation/QC/",cohorts[1],"/",cohorts[1],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    common<-intersect(a$V2,b$V2)
    print(cohorts[1])
    print(cohorts[i])
    print(length(common))
    rm(a)
  } else {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    common<-intersect(common,b$V2)
    print(cohorts[i])
    print(length(common))
  }
}

# [1] "prism_nfe_gsa"
# [1] "prism_nfe_gwas"
# [1] 62075

# keep only those in autosomal chromosomes
common_2<-b[which(b$V2 %in% common),]
common_2<-common_2[which(common_2$V1 %in% seq(1:22)),]
dim(common_2)
# [1] 62075    6

write.table(common_2[,2,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_sweden",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



##################################
# list of variants across CEDARS #
##################################

library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","cedars_gsa")


for (i in 2:length(cohorts)){
  if(i==2) {
    a<-read.table(paste(path,"pre_imputation/QC/",cohorts[1],"/",cohorts[1],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    common<-intersect(a$V2,b$V2)
    print(cohorts[1])
    print(cohorts[i])
    print(length(common))
    rm(a)
  } else {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.bim",sep=""),head=F)
    common<-intersect(common,b$V2)
    print(cohorts[i])
    print(length(common))
  }
}

# [1] "cedars_370k_old_gwas"
# [1] "cedars_610k_old_gwas"
# [1] 337664
# [1] "cedars_omni_old_gwas"
# [1] 161119
# [1] "cedars_gsa"
# [1] 46849


# keep only those in autosomal chromosomes
common_2<-b[which(b$V2 %in% common),]
common_2<-common_2[which(common_2$V1 %in% seq(1:22)),]
dim(common_2)
# [1] 44870      6

write.table(common_2[,2,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_cedars",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



###########################################################################################################################################################
###########################################################################################################################################################


# after running intra-centre dup checking run following:

##################################
# list of variants across IIBDGC #
##################################

# create final file for all studies:

path_gwas=/path/to/ibdgwas/IIBDGC/
  
# studies=(australia_omniexome pittsburgh_gsa spain_gsa italy_gsa netherlands_gsa slovenia_gsa sweden_gsa basque_gsa lithuania_gsa finland_illugwas chop_old_gwas norway_affy6_old_gwas)
# 
# for i in ${studies[@]}
# do
# ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy.bim
# done
# 
# for i in ${studies[@]}
# do
# /path/to/software/plink_linux_x86_64_20181202/./plink \
# --bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample \
# --make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy
# done

#### /software/R-4.3.1/bin/R

library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"

# cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
#            ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
#            ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
#            ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
#            ,"prism_nfe_gwas","finland_illugwas"
#            ,"chop_old_gwas","german_illu550_old_gwas","german_affy6_old_gwas","norway_affy6_old_gwas"
#            ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_cd_old_gwas") # illu 317

# no german_illu, no chop
cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","swedish_uc_old_gwas"
           ,"mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

length(cohorts)
# [1] 59


for (i in 2:length(cohorts)){
  if(i==2) {
    a<-read.table(paste(path,"pre_imputation/QC/",cohorts[1],"/",cohorts[1],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy.bim",sep=""),head=F)
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy.bim",sep=""),head=F)
    common<-intersect(a$V2,b$V2)
    print(cohorts[1])
    print(cohorts[i])
    print(length(common))
    rm(a)
  } else {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy.bim",sep=""),head=F)
    common<-intersect(common,b$V2)
    print(cohorts[i])
    print(length(common))
  }
}

###########
# [1] "australia_omniexome"
# [1] "gwas1"
# [1] 105118
# [1] "gwas2"
# [1] 99829
# [1] "all_hce"
# [1] 35672
# [1] "pittsburgh_gsa"
# [1] 32721
# [1] "spain_gsa"
# [1] 16404
# [1] "italy_gsa"
# [1] 9274
# [1] "kiel_austria_sibdcs_gsa"
# [1] 9000
# [1] "netherlands_gsa"
# [1] 9000
# [1] "slovenia_gsa"
# [1] 8981
# [1] "sweden_gsa"
# [1] 8971
# [1] "niddk_broad_gsa"
# [1] 8957
# [1] "niddk_feinstein_gsa"
# [1] 8924
# [1] "basque_gsa"
# [1] 8924
# [1] "prism_nfe_gsa"
# [1] 8920
# [1] "lithuania_gsa"
# [1] 8914
# [1] "belgium_louis_gsa"
# [1] 8914
# [1] "belgium_franchimont_gsa"
# [1] 8911
# [1] "belgium_vermeire_gsa"
# [1] 8911
# [1] "prism_nfe_gwas"
# [1] 7828
# [1] "finland_illugwas"
# [1] 7826
# [1] "german_affy6_old_gwas"
# [1] 7259
# [1] "norway_affy6_old_gwas"
# [1] 6964
# [1] "belgium_inf1_old_gwas"
# [1] 4398
# [1] "belgium_inf2_old_gwas"
# [1] 4281
# [1] "niddk_old_gwas"
# [1] 4223
# [1] "cedars_370k_old_gwas"
# [1] 4220
# [1] "cedars_610k_old_gwas"
# [1] 4220
# [1] "cedars_omni_old_gwas"
# [1] 3280
# [1] "swedish_uc_old_gwas"
# [1] 3190
# [1] "mccauley_gsa"
# [1] 3190
# [1] "ccfa_gsa"
# [1] 3190
# [1] "cedars_gsa"
# [1] 3190
# [1] "mccauley_gsa"
# [1] 3190
# [1] "ccfa_gsa"
# [1] 3190
# [1] "cedars_gsa"
# [1] 3190
# [1] "bernstein_gsa"
# [1] 3190
# [1] "farkkila_gsa"
# [1] 3188
# [1] "franchimont_gsa"
# [1] 3188
# [1] "franke_gsa"
# [1] 3188
# [1] "helmsley_prism_gsa"
# [1] 3188
# [1] "helmsley_xavier_prism_gsa"
# [1] 3188
# [1] "hyams_protect_gsa"
# [1] 3188
# [1] "lewis_sparc_gsa"
# [1] 3188
# [1] "mccauley_new_gsa"
# [1] 3188
# [1] "mcgovern_gsa"
# [1] 3188
# [1] "moayyedi_imagine_gsa"
# [1] 3188
# [1] "newberry_share_gsa"
# [1] 3188
# [1] "niddk_cho_gsa"
# [1] 3188
# [1] "niddk_duerr_gsa"
# [1] 3188
# [1] "niddk_rioux_gsa"
# [1] 3188
# [1] "niddk_silverberg_gsa"
# [1] 3188
# [1] "palotie_hus_gsa"
# [1] 3188
# [1] "pekow_share_gsa"
# [1] 3188
# [1] "rioux_igenomed_gsa"
# [1] 3188
# [1] "sands_msccr_gsa"
# [1] 3188
# [1] "stampfer_gsa"
# [1] 3186
# [1] "vermeire_gsa"
# [1] 3186
# [1] "weersma_gsa"
# [1] 3186
# [1] "xavier_prism_gsa"
# [1] 3186
# [1] "xavier_share_gsa"
# [1] 3186
###########

# keep only those in autosomal chromosomes
common_2<-b[which(b$V2 %in% common),]
common_2<-common_2[which(common_2$V1 %in% seq(1:22)),]
dim(common_2)
# [1] 3186 6


write.table(common_2[,2,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_iibdgc",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


###########################################################################################################################################################

rm(b,common)

# to run PCA with 1000GP:

table(common_2$V1)
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
# 263 293 186 190 212 216 160 211 156 183 171 106 131  94  74 105  60 123  52  78 
# 21  22 
# 56  66 

for (i in c(1:22)) {

  tmp<-fread(paste(path,"pre_imputation/QC/1000gp/1000GP_chr",i,"_b37_edited.bim",sep=""),head=F)
  tmp<-tmp[which(tmp$V2 %in% common_2$V2),]

  print(paste("chr",i))
  print(nrow(tmp))

  if (i==1) {
    bb<-tmp
    rm(tmp)
  } else {
    bb<-rbind(bb,tmp)
    rm(tmp)
  }
}

table(bb$V1)
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
# 263 293 186 190 212 216 160 211 156 183 171 106 131  94  74 105  60 123  52  78 
# 21  22 
# 56  66 

dim(bb)
# [1] 3186    6

# all variants also in 1000GP

###########################################################################################################################################################

# NOT NEEDED ANYMORE

# # to run PCA with UKBB 
# 
# # /software/R-4.3.1/bin/R
# 
# library(data.table)
# 
# path<-"/path/to/ibdgwas/IIBDGC/"
# 
# var<-read.table(paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_iibdgc",sep=""),head=F)
# ukbb<-read.table("/path/to/project",head=F)
# dim(ukbb)
# # [1] 556174      6 OK
# 
# 
# dim(ukbb[which(ukbb$V2 %in% var$V1),])
# # [1] 825   6
# 
# # final check:
# 
# tmp<-var
# tmp$XX<-gsub("_.*","",tmp$V1)
# 
# ukbb$XX<-gsub("_.*","",ukbb$V2)
# 
# dim(ukbb[which(ukbb$XX %in% tmp$XX),])
# # [1] 825   7 # OK
# 
# write.table(var[which(var$V1 %in% ukbb$V2),,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_iibdgc_ukbb",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#################################################################################################################################################

# create a list between old broad/gsa batches and new (includes gsa and Illumina Infinium Exome)

#### /software/R-4.3.1/bin/R

library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"

# use list of post_QC variants first:

cohorts<-c("basque_gsa","belgium_franchimont_gsa","belgium_louis_gsa","belgium_vermeire_gsa","ccfa_gsa","cedars_gsa",
           "italy_gsa","kiel_austria_sibdcs_gsa","lithuania_gsa","mccauley_gsa","netherlands_gsa",
           "niddk_broad_gsa","niddk_feinstein_gsa","prism_nfe_gsa","slovenia_gsa","sweden_gsa",
           "prism_nfe_gwas","finland_illugwas",
           "helmsley_prism_gsa","helmsley_xavier_prism_gsa",
           "bernstein_gsa","farkkila_gsa","franchimont_gsa","franke_gsa",
           "hyams_protect_gsa","lewis_sparc_gsa","mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa",
           "newberry_share_gsa","niddk_cho_gsa","niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa",
           "palotie_hus_gsa","pekow_share_gsa","rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa",
           "vermeire_gsa","weersma_gsa","xavier_prism_gsa","xavier_share_gsa")

length(cohorts)
# [1] 43

list_good_gsa_illuminainfexome<-read.table(paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_prism",sep=""))

for (i in 1:length(cohorts)){
  if(i==1) {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip.bim",sep=""),head=F)
    common<-intersect(list_good_gsa_illuminainfexome$V1,b$V2)
    print(cohorts[i])
    print(length(common))
    rm(b)
  } else {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip.bim",sep=""),head=F)
    common<-intersect(common,b$V2)
    print(cohorts[i])
    print(length(common))
  }
}

# keep only those in autosomal chromosomes
common_2<-b[which(b$V2 %in% common),]
common_2<-common_2[which(common_2$V1 %in% seq(1:22)),]
dim(common_2)
# [1] 44870      6

write.table(common_2[,2,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_old_new_broad_batch",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#################################################################################################################################################

# create a list between old broad/gsa batches and new 
# (includes gsa and Illumina Infinium Exome) - NOW WITH VARIANTS THAT PASS QC

#### /software/R-4.3.1/bin/R

library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"

# use list of post_QC variants first:

cohorts<-c("basque_gsa","belgium_franchimont_gsa","belgium_louis_gsa","belgium_vermeire_gsa","ccfa_gsa","cedars_gsa",
           "italy_gsa","kiel_austria_sibdcs_gsa","lithuania_gsa","mccauley_gsa","netherlands_gsa",
           "niddk_broad_gsa","niddk_feinstein_gsa","prism_nfe_gsa","slovenia_gsa","sweden_gsa",
           "prism_nfe_gwas","finland_illugwas",
           "helmsley_prism_gsa","helmsley_xavier_prism_gsa",
           "bernstein_gsa","farkkila_gsa","franchimont_gsa","franke_gsa",
           "hyams_protect_gsa","lewis_sparc_gsa","mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa",
           "newberry_share_gsa","niddk_cho_gsa","niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa",
           "palotie_hus_gsa","pekow_share_gsa","rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa",
           "vermeire_gsa","weersma_gsa","xavier_prism_gsa","xavier_share_gsa")

length(cohorts)
# [1] 43

list_good_gsa_illuminainfexome<-read.table(paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_prism",sep=""))

for (i in 1:length(cohorts)){
  if(i==1) {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim",sep=""),head=F)
    common<-intersect(list_good_gsa_illuminainfexome$V1,b$V2)
    print(cohorts[i])
    print(length(common))
    rm(b)
  } else {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim",sep=""),head=F)
    common<-intersect(common,b$V2)
    print(cohorts[i])
    print(length(common))
  }
}

# keep only those in autosomal chromosomes
common_2<-b[which(b$V2 %in% common),]
common_2<-common_2[which(common_2$V1 %in% seq(1:22)),]
dim(common_2)
# [1] 108295      6

write.table(common_2[,2,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_old_new_broad_batch_2",sep=""),col.names=F,row.names=F,quote=F,sep="\t")




####################################################################################################################################################
####################################################################################################################################################

# STEP 25 UPDATES

###################################
# list of variants across UKIBDGC #
###################################

library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("gwas1","gwas2","all_hce")


for (i in 2:length(cohorts)){
  if(i==2) {
    a<-read.table(paste(path,"pre_imputation/QC/",cohorts[1],"/",cohorts[1],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(a$V2,b$V2)
    print(cohorts[1])
    print(cohorts[i])
    print(length(common))
    rm(a)
  } else {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(common,b$V2)
    print(cohorts[i])
    print(length(common))
  }
}

# [1] "gwas1"
# [1] "gwas2"
# [1] 399424
# [1] "all_hce"
# [1] 36537

# keep only those in autosomal chromosomes
table(b[which(b$V2 %in% common),"V1"])
#

common_2<-b[which(b$V2 %in% common),]
common_2<-common_2[which(common_2$V1 %in% seq(1:22)),]
dim(common_2)
# [1] 35897     6

write.table(common_2[,2,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_step25_ukibdgc",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

##################################
# list of variants across German #
##################################

library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("kiel_austria_sibdcs_gsa","german_affy6_old_gwas","franke_gsa")

for (i in 2:length(cohorts)){
  if(i==2) {
    a<-read.table(paste(path,"pre_imputation/QC/",cohorts[1],"/",cohorts[1],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(a$V2,b$V2)
    print(cohorts[1])
    print(cohorts[i])
    print(length(common))
    rm(a)
  } else {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(common,b$V2)
    print(cohorts[i])
    print(length(common))
  }
}

# [1] "kiel_austria_sibdcs_gsa"
# [1] "german_affy6_old_gwas"
# [1] 70450
# [1] "franke_gsa"
# [1] 70171

# keep only those in autosomal chromosomes
common_2<-b[which(b$V2 %in% common),]
common_2<-common_2[which(common_2$V1 %in% seq(1:22)),]
dim(common_2)
# [1] 67623    6

write.table(common_2[,2,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_step25_german",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


#################################
# list of variants across NIDDK #
#################################

library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("niddk_broad_gsa","niddk_feinstein_gsa","niddk_old_gwas","niddk_cho_gsa","niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa")


for (i in 2:length(cohorts)){
  if(i==2) {
    a<-read.table(paste(path,"pre_imputation/QC/",cohorts[1],"/",cohorts[1],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(a$V2,b$V2)
    print(cohorts[1])
    print(cohorts[i])
    print(length(common))
    rm(a)
  } else {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(common,b$V2)
    print(cohorts[i])
    print(length(common))
  }
}

# [1] "niddk_broad_gsa"
# [1] "niddk_feinstein_gsa"
# [1] 564283
# [1] "niddk_old_gwas"
# [1] 61855
# [1] "niddk_cho_gsa"
# [1] 61776
# [1] "niddk_duerr_gsa"
# [1] 61770
# [1] "niddk_rioux_gsa"
# [1] 61769
# [1] "niddk_silverberg_gsa"
# [1] 61761


# keep only those in autosomal chromosomes
common_2<-b[which(b$V2 %in% common),]
common_2<-common_2[which(common_2$V1 %in% seq(1:22)),]
dim(common_2)
# [1] [1] 59499    6

write.table(common_2[,2,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_step25_niddk",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


#################################
# list of variants across Spain #
#################################

library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("spain_gsa","basque_gsa")


for (i in 2:length(cohorts)){
  if(i==2) {
    a<-read.table(paste(path,"pre_imputation/QC/",cohorts[1],"/",cohorts[1],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(a$V2,b$V2)
    print(cohorts[1])
    print(cohorts[i])
    print(length(common))
    rm(a)
  } else {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(common,b$V2)
    print(cohorts[i])
    print(length(common))
  }
}

# [1] "spain_gsa"
# [1] "basque_gsa"
# [1] 117176

# keep only those in autosomal chromosomes
common_2<-b[which(b$V2 %in% common),]
common_2<-common_2[which(common_2$V1 %in% seq(1:22)),]
dim(common_2)
# [1] 106519    6

write.table(common_2[,2,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_step25_spain",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

###################################
# list of variants across Belgium #
###################################

library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("belgium_inf1_old_gwas","belgium_inf2_old_gwas","belgium_louis_gsa","belgium_franchimont_gsa"
           ,"belgium_vermeire_gsa","franchimont_gsa","vermeire_gsa")


for (i in 2:length(cohorts)){
  if(i==2) {
    a<-read.table(paste(path,"pre_imputation/QC/",cohorts[1],"/",cohorts[1],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(a$V2,b$V2)
    print(cohorts[1])
    print(cohorts[i])
    print(length(common))
    rm(a)
  } else {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(common,b$V2)
    print(cohorts[i])
    print(length(common))
  }
}
# [1] "belgium_inf1_old_gwas"
# [1] "belgium_inf2_old_gwas"
# [1] 283815
# [1] "belgium_louis_gsa"
# [1] 60998
# [1] "belgium_franchimont_gsa"
# [1] 60652
# [1] "belgium_vermeire_gsa"
# [1] 60599
# [1] "franchimont_gsa"
# [1] 60590
# [1] "vermeire_gsa"
# [1] 60579

# keep only those in autosomal chromosomes
common_2<-b[which(b$V2 %in% common),]
common_2<-common_2[which(common_2$V1 %in% seq(1:22)),]
dim(common_2)
# [1] 58624    6

write.table(common_2[,2,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_step25_belgium",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


#################################
# list of variants across PRISM #
#################################

library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("prism_nfe_gsa","prism_nfe_gwas","helmsley_prism_gsa","helmsley_xavier_prism_gsa","xavier_prism_gsa","xavier_share_gsa")


for (i in 2:length(cohorts)){
  if(i==2) {
    a<-read.table(paste(path,"pre_imputation/QC/",cohorts[1],"/",cohorts[1],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(a$V2,b$V2)
    print(cohorts[1])
    print(cohorts[i])
    print(length(common))
    rm(a)
  } else {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(common,b$V2)
    print(cohorts[i])
    print(length(common))
  }
}

# [1] "prism_nfe_gsa"
# [1] "prism_nfe_gwas"
# [1] 124061
# [1] "helmsley_prism_gsa"
# [1] 124061
# [1] "helmsley_xavier_prism_gsa"
# [1] 124030
# [1] "xavier_prism_gsa"
# [1] 123983
# [1] "xavier_share_gsa"
# [1] 123709

# keep only those in autosomal chromosomes
common_2<-b[which(b$V2 %in% common),]
common_2<-common_2[which(common_2$V1 %in% seq(1:22)),]
dim(common_2)
# [1] 119440    6

write.table(common_2[,2,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_step25_prism",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


##################################
# list of variants across SWEDEN #
##################################


library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("swedish_uc_old_gwas","sweden_gsa")


for (i in 2:length(cohorts)){
  if(i==2) {
    a<-read.table(paste(path,"pre_imputation/QC/",cohorts[1],"/",cohorts[1],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(a$V2,b$V2)
    print(cohorts[1])
    print(cohorts[i])
    print(length(common))
    rm(a)
  } else {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(common,b$V2)
    print(cohorts[i])
    print(length(common))
  }
}

# [1] "swedish_uc_old_gwas"
# [1] "sweden_gsa"
# [1] 62063

# keep only those in autosomal chromosomes
common_2<-b[which(b$V2 %in% common),]
common_2<-common_2[which(common_2$V1 %in% seq(1:22)),]
dim(common_2)
# [1] 62063    6

write.table(common_2[,2,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_step25_sweden",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



##################################
# list of variants across CEDARS #
##################################

library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","cedars_gsa","mcgovern_gsa")


for (i in 2:length(cohorts)){
  if(i==2) {
    a<-read.table(paste(path,"pre_imputation/QC/",cohorts[1],"/",cohorts[1],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(a$V2,b$V2)
    print(cohorts[1])
    print(cohorts[i])
    print(length(common))
    rm(a)
  } else {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(common,b$V2)
    print(cohorts[i])
    print(length(common))
  }
}

# [1] "cedars_370k_old_gwas"
# [1] "cedars_610k_old_gwas"
# [1] 337609
# [1] "cedars_omni_old_gwas"
# [1] 161080
# [1] "cedars_gsa"
# [1] 46819
# [1] "mcgovern_gsa"
# [1] 46774


# keep only those in autosomal chromosomes
common_2<-b[which(b$V2 %in% common),]
common_2<-common_2[which(common_2$V1 %in% seq(1:22)),]
dim(common_2)
# [1] 44796      6

write.table(common_2[,2,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_step25_cedars",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#######################################
# list of variants across NETHERLANDS #
#######################################

library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("netherlands_gsa","weersma_gsa")


for (i in 2:length(cohorts)){
  if(i==2) {
    a<-read.table(paste(path,"pre_imputation/QC/",cohorts[1],"/",cohorts[1],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(a$V2,b$V2)
    print(cohorts[1])
    print(cohorts[i])
    print(length(common))
    rm(a)
  } else {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(common,b$V2)
    print(cohorts[i])
    print(length(common))
  }
}

# [1] "netherlands_gsa"
# [1] "weersma_gsa"
# [1] 559129


# keep only those in autosomal chromosomes
common_2<-b[which(b$V2 %in% common),]
common_2<-common_2[which(common_2$V1 %in% seq(1:22)),]
dim(common_2)
# [1] 545037      6

write.table(common_2[,2,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_step25_netherlands",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

####################################
# list of variants across MCCAULEY #
####################################

library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("mccauley_new_gsa","mccauley_gsa")


for (i in 2:length(cohorts)){
  if(i==2) {
    a<-read.table(paste(path,"pre_imputation/QC/",cohorts[1],"/",cohorts[1],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(a$V2,b$V2)
    print(cohorts[1])
    print(cohorts[i])
    print(length(common))
    rm(a)
  } else {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(common,b$V2)
    print(cohorts[i])
    print(length(common))
  }
}

# [1] "mccauley_new_gsa"
# [1] "mccauley_gsa"
# [1] 556134

# keep only those in autosomal chromosomes
common_2<-b[which(b$V2 %in% common),]
common_2<-common_2[which(common_2$V1 %in% seq(1:22)),]
dim(common_2)
# [1] 542092      6

write.table(common_2[,2,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_step25_mccauley",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

####################################
# list of variants across CCFA #
####################################

library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("ccfa_gsa","lewis_sparc_gsa")


for (i in 2:length(cohorts)){
  if(i==2) {
    a<-read.table(paste(path,"pre_imputation/QC/",cohorts[1],"/",cohorts[1],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(a$V2,b$V2)
    print(cohorts[1])
    print(cohorts[i])
    print(length(common))
    rm(a)
  } else {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(common,b$V2)
    print(cohorts[i])
    print(length(common))
  }
}

# [1] "ccfa_gsa"
# [1] "lewis_sparc_gsa"
# [1] 598973

# keep only those in autosomal chromosomes
common_2<-b[which(b$V2 %in% common),]
common_2<-common_2[which(common_2$V1 %in% seq(1:22)),]
dim(common_2)
# [1] 583991      6

write.table(common_2[,2,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_step25_ccfa",sep=""),col.names=F,row.names=F,quote=F,sep="\t")



###########################################################################################################################################################
###########################################################################################################################################################


# after running intra-centre dup checking run following:

##################################
# list of variants across IIBDGC #
##################################

# create final file for all studies:

path_gwas=/path/to/ibdgwas/IIBDGC/
  
#### /software/R-4.3.1/bin/R

library(data.table)

path<-"/path/to/ibdgwas/IIBDGC/"

# no german_illu, no chop
cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","swedish_uc_old_gwas"
           ,"mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

length(cohorts)
# [1] 58


for (i in 2:length(cohorts)){
  if(i==2) {
    a<-read.table(paste(path,"pre_imputation/QC/",cohorts[1],"/",cohorts[1],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(a$V2,b$V2)
    print(cohorts[1])
    print(cohorts[i])
    print(length(common))
    rm(a)
  } else {
    b<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.bim",sep=""),head=F)
    common<-intersect(common,b$V2)
    print(cohorts[i])
    print(length(common))
  }
}

###########
# [1] "australia_omniexome"
# [1] "gwas1"
# [1] 104721
# [1] "gwas2"
# [1] 98793
# [1] "all_hce"
# [1] 35334
# [1] "pittsburgh_gsa"
# [1] 32433
# [1] "spain_gsa"
# [1] 16274
# [1] "italy_gsa"
# [1] 9205
# [1] "kiel_austria_sibdcs_gsa"
# [1] 8932
# [1] "netherlands_gsa"
# [1] 8932
# [1] "slovenia_gsa"
# [1] 8913
# [1] "sweden_gsa"
# [1] 8903
# [1] "niddk_broad_gsa"
# [1] 8889
# [1] "niddk_feinstein_gsa"
# [1] 8854
# [1] "basque_gsa"
# [1] 8853
# [1] "prism_nfe_gsa"
# [1] 8849
# [1] "lithuania_gsa"
# [1] 8838
# [1] "belgium_louis_gsa"
# [1] 8838
# [1] "belgium_franchimont_gsa"
# [1] 8834
# [1] "belgium_vermeire_gsa"
# [1] 8834
# [1] "prism_nfe_gwas"
# [1] 7764
# [1] "finland_illugwas"
# [1] 7762
# [1] "german_affy6_old_gwas"
# [1] 7184
# [1] "norway_affy6_old_gwas"
# [1] 6860
# [1] "belgium_inf1_old_gwas"
# [1] 4344
# [1] "belgium_inf2_old_gwas"
# [1] 4228
# [1] "niddk_old_gwas"
# [1] 4170
# [1] "cedars_370k_old_gwas"
# [1] 4167
# [1] "cedars_610k_old_gwas"
# [1] 4167
# [1] "cedars_omni_old_gwas"
# [1] 3242
# [1] "swedish_uc_old_gwas"
# [1] 3152
# [1] "mccauley_gsa"
# [1] 3152
# [1] "ccfa_gsa"
# [1] 3152
# [1] "cedars_gsa"
# [1] 3152
# [1] "bernstein_gsa"
# [1] 3152
# [1] "farkkila_gsa"
# [1] 3150
# [1] "franchimont_gsa"
# [1] 3150
# [1] "franke_gsa"
# [1] 3150
# [1] "helmsley_prism_gsa"
# [1] 3150
# [1] "helmsley_xavier_prism_gsa"
# [1] 3150
# [1] "hyams_protect_gsa"
# [1] 3150
# [1] "lewis_sparc_gsa"
# [1] 3150
# [1] "mccauley_new_gsa"
# [1] 3150
# [1] "mcgovern_gsa"
# [1] 3150
# [1] "moayyedi_imagine_gsa"
# [1] 3150
# [1] "newberry_share_gsa"
# [1] 3150
# [1] "niddk_cho_gsa"
# [1] 3150
# [1] "niddk_duerr_gsa"
# [1] 3150
# [1] "niddk_rioux_gsa"
# [1] 3150
# [1] "niddk_silverberg_gsa"
# [1] 3150
# [1] "palotie_hus_gsa"
# [1] 3150
# [1] "pekow_share_gsa"
# [1] 3150
# [1] "rioux_igenomed_gsa"
# [1] 3150
# [1] "sands_msccr_gsa"
# [1] 3150
# [1] "stampfer_gsa"
# [1] 3148
# [1] "vermeire_gsa"
# [1] 3148
# [1] "weersma_gsa"
# [1] 3148
# [1] "xavier_prism_gsa"
# [1] 3148
# [1] "xavier_share_gsa"
# [1] 3148
###########

# keep only those in autosomal chromosomes
common_2<-b[which(b$V2 %in% common),]
common_2<-common_2[which(common_2$V1 %in% seq(1:22)),]
dim(common_2)
# [1] 3186 6


write.table(common_2[,2,drop=F],paste(path,"pre_imputation/QC/relatedness/list_common_variants_among_cohorts_iibdgc_step25",sep=""),col.names=F,row.names=F,quote=F,sep="\t")


###########################################################################################################################################################

rm(b,common)

# to run PCA with 1000GP:

table(common_2$V1)
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
# 263 293 186 190 212 216 160 211 156 183 171 106 131  94  74 105  60 123  52  78 
# 21  22 
# 56  66 

for (i in c(1:22)) {
  
  tmp<-fread(paste(path,"pre_imputation/QC/1000gp/1000GP_chr",i,"_b37_edited.bim",sep=""),head=F)
  tmp<-tmp[which(tmp$V2 %in% common_2$V2),]
  
  print(paste("chr",i))
  print(nrow(tmp))
  
  if (i==1) {
    bb<-tmp
    rm(tmp)
  } else {
    bb<-rbind(bb,tmp)
    rm(tmp)
  }
}

table(bb$V1)
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
# 258 290 186 186 211 213 157 208 155 179 169 101 130  94  74 105  60 122  52  77 
# 21  22 
# 56  65 


dim(bb)
# [1] 3148    6





