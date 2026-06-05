# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# EVALUATE RELATEDNESS BETWEEN cedars SAMPLES

###########################
# 1.- SELECT THE VARIANTS #
###########################

# SEE SCRIPT:
# create_list_common_variants_among_cohorts

########################
# 2.- EXTRACT VARIANTS #
########################

# hpc-server

path=/path/to/ibdgwas/IIBDGC/
  
studies=(cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas cedars_gsa)
  
path_variants=${path}pre_imputation/QC/relatedness/list_common_variants_among_cohorts_cedars
  
for i in ${studies[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample \
--allow-no-sex \
--extract ${path_variants} \
--make-bed \
--out ${path}pre_imputation/QC/relatedness/${i}_subset_cedars
done
  
# estimate missingness for later step
for i in ${studies[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample \
--allow-no-sex \
--missing \
--out ${path}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample
done
  
  
wc -l ${path}pre_imputation/QC/relatedness/*_subset_cedars.fam
# 605 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_370k_old_gwas_subset_cedars.fam
# 889 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_610k_old_gwas_subset_cedars.fam
# 3076 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_gsa_subset_cedars.fam
# 1219 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_omni_old_gwas_subset_cedars.fam
# 5788 total


wc -l ${path}pre_imputation/QC/relatedness/*_subset_cedars.bim
# 44870  OK
  
##############################################################################################################################################
  
#######################
# 3.- MERGE ALL FILES #
#######################

#### /software/R-4.3.1/bin/R
  
library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
  
path<-"/path/to/ibdgwas/IIBDGC/"
  
cohorts<-c("cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","cedars_gsa")
  
cohorts<-cohorts[-1]
  
dat<-matrix(ncol=3,nrow=length(cohorts))
dat<-as.data.frame(dat)
  
for (i in 1:nrow(dat)){
  dat[i,1]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_cedars.bed",sep="")
  dat[i,2]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_cedars.bim",sep="")
  dat[i,3]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_cedars.fam",sep="")
}
  
write.table(dat,paste(path,"pre_imputation/QC/relatedness/list_cohorts_tomerge_cedars.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  
#############################
  
# hpc-server
  
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/relatedness/cedars_370k_old_gwas_subset_cedars \
--allow-no-sex \
--merge-list ${path}pre_imputation/QC/relatedness/list_cohorts_tomerge_cedars.txt \
--make-bed --out ${path}pre_imputation/QC/relatedness/cedars_merged
  
# 44870 variants and 4841 people pass filters and QC.
# Among remaining phenotypes, 4837 are cases and 4 are controls.
  
  
############################
# 3.- ESTIMATE RELATEDNESS #
############################

# The current release is KING version 2.2.4 (released on October 11, 2019). 

# --related provides integrative, fast, and accurate inference for close relationships. This option is highly recommended, 
# especially when dealing with biobank-level datasets. Integration of the IBD segment inference furthur improves the inference accuracy. 
# In our tests, --related has been successfully applied to datasets consisting of ~10 million samples. When "--rplot" is specified, 
# several relationship plots are generated automatically. --related --degree 2 specifies that only related pairs
# (up to the 2nd-degree in this case) between families are included in the output. Specifically all pairs across families with 
# a kinship coefficient less than 0.0884 will be excluded from the output. More details of 
# --related analysis are available in the INTEGRATED RELATIONSHIP INFERENCE section later in this tutorial.

# --duplicate implements the fastest (and accurate) algorithm to identify duplicates or MZ twins. The running time is in seconds, 
# unless the number of samples is > 1,000,000 in which case a few minutes may be needed. In our tests, --duplicate has been successfully 
# applied to datasets consisting of ~10 million samples. One potential application of the duplicate analysis is to identify duplicates 
# accross different studies, in which case multiple datasets can be read in conveniently as shown in GENERAL INPUT FILES section.

# The amount of computer memory required by KING analysis is modest, at ~N ✕ M / 4 (where N is the number of samples and M is the number of SNPs)
# plus a small percentage of overhead cost. E.g., for a dataset consisting of 100,000 samples each genotyped at 1,000,000 SNPs, 
# the required memory size is ~25GB.


# (156413*2713/4)/1E+6
# [1] 106.0871

path=/path/to/ibdgwas/IIBDGC/
MEM=1000

bsub -J"kg" -M"$MEM" -n8 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path}pre_imputation/QC/relatedness/log/stderr_king_cedars \
-o ${path}pre_imputation/QC/relatedness/log/stdout_king_cedars \
"/path/to/software/./king \
-b ${path}pre_imputation/QC/relatedness/cedars_merged.bed \
--related --cpus 8 --prefix ${path}pre_imputation/QC/relatedness/cedars_merged"
# Job <983092> is submitted to queue <normal>.

less ${path}pre_imputation/QC/relatedness/log/stdout_king_cedars
# Relationship summary (total relatives: 0 by pedigree, 470 by inference)
# MZ      PO      FS      2nd
# =====================================================
# Inference     357     332     102     0




#### R

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","cedars_gsa")

for (i in 1:length(cohorts)) {
  tmp<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],
                        "_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.fam",sep=""),head=F)
  sample_miss<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],
                                "_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.imiss",sep=""),head=T)
  tmp<-merge(tmp,sample_miss[,c("IID","F_MISS")],by.x="V2",by.y="IID",all.x=T)
  tmp$cohort<-as.character(cohorts[i])
  tmp$V1<-as.character(tmp$V1)
  tmp$V2<-as.character(tmp$V2)
  
  if(cohorts[i]=="german_illu550_old_gwas") {
    tmp$V1<-paste(tmp$V1,"_illu550",sep="")
    tmp$V2<-paste(tmp$V2,"_illu550",sep="")
  }
  
  if(i==1){
    dat<-tmp
  }else{
    dat<-rbind(dat,tmp)
  }
}

table(dat$cohort)
# cedars_370k_old_gwas cedars_610k_old_gwas           cedars_gsa 
#                  605                  889                 3076 
# cedars_omni_old_gwas 
#                 1219 

# double check all samples included in exercise
famall<-read.table(paste(path,"pre_imputation/QC/relatedness/cedars_merged.fam",sep=""),head=F)

dim(famall)
# [1] 5789     6
dim(dat)
# [1] 5789     8
dim(dat[which(dat$V1 %in% famall$V1),])
# [1] 5789     8

colnames(dat)[6]<-"pheno"
aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary)
#                Group.1  x.1  x.2
# 1 cedars_370k_old_gwas    0  605
# 2 cedars_610k_old_gwas    0  889
# 3           cedars_gsa  947 2128
# 4 cedars_omni_old_gwas    4 1215

colnames(dat)[5]<-"sex"
aggregate(as.factor(dat$sex), by=list(dat$cohort), FUN=summary)
#                Group.1  x.0  x.1  x.2
# 1 cedars_370k_old_gwas    1  324  280
# 2 cedars_610k_old_gwas    3  465  421
# 3           cedars_gsa    9 1481 1586
# 4 cedars_omni_old_gwas    0  594  625


table(dat$V1==dat$V2)
# TRUE 
# 5788

kin<-read.table(paste(path,"pre_imputation/QC/relatedness/cedars_merged.kin0",sep=""),head=T)
table(kin$InfType)
# Dup/MZ     FS     PO 
#    357    102    332 

table(cut(kin$Kinship,breaks=c(1,0.354,0.177)))
# (0.177,0.354]     (0.354,1] 
# 434           357 

dup<-kin[which(kin$Kinship>0.354 & kin$InfType=="Dup/MZ"),c("FID1","FID2","Kinship")]
summary(dup$Kinship)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4938  0.4993  0.4996  0.4994  0.4997  0.5000 

dup<-merge(dup,dat[,c("V1","cohort","sex","pheno","F_MISS")],by.x="FID2",by.y="V1",all.x=T,sort=F)
colnames(dup)[(ncol(dup)-3):ncol(dup)]<-paste(colnames(dup)[(ncol(dup)-3):ncol(dup)],"FID2",sep="_")

dup<-merge(dup,dat[,c("V1","cohort","sex","pheno","F_MISS")],by.x="FID1",by.y="V1",all.x=T,sort=F)
colnames(dup)[(ncol(dup)-3):ncol(dup)]<-paste(colnames(dup)[(ncol(dup)-3):ncol(dup)],"FID1",sep="_")

x<-table(dup$cohort_FID1,dup$cohort_FID2)
x
#                      cedars_370k_old_gwas cedars_610k_old_gwas
# cedars_370k_old_gwas                    0                    0
# cedars_gsa                             52                  122
# 
#                      cedars_omni_old_gwas
# cedars_370k_old_gwas                    1
# cedars_gsa                            182

xx<-dup[,c("cohort_FID1","cohort_FID2")]
yy<-dup[,c("cohort_FID2","cohort_FID1")]
colnames(yy)<-colnames(xx)

dim(xx)
# [1] 357    2

xx<-rbind(xx,yy)
x<-table(xx$cohort_FID1,xx$cohort_FID2)
x
#                       cedars_370k_old_gwas cedars_610k_old_gwas cedars_gsa
# cedars_370k_old_gwas                    0                    0         52
# cedars_610k_old_gwas                    0                    0        122
# cedars_gsa                             52                  122          0
# cedars_omni_old_gwas                    1                    0        182
# 
#                      cedars_omni_old_gwas
# cedars_370k_old_gwas                    1
# cedars_610k_old_gwas                    0
# cedars_gsa                            182
# cedars_omni_old_gwas                    0


write.table(x,paste(path,"pre_imputation/QC/relatedness/cedars_dup_summary_table",sep=""),col.names=T,row.names=T,sep="\t",quote=F)


table(dup$sex_FID1,dup$sex_FID2)
#     1   2
# 1 166   0
# 2   0 191

dup<-dup[,c("FID1","FID2","Kinship","cohort_FID1","pheno_FID1","F_MISS_FID1","cohort_FID2","pheno_FID2","F_MISS_FID2")]
table(dup$pheno_FID1,dup$pheno_FID2)
#     1   2
# 1   2   0
# 2   0 355


colnames(dat)[6]<-"pheno"
aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary)
# Group.1  x.1  x.2
# 1 cedars_370k_old_gwas    0  605
# 2 cedars_610k_old_gwas    0  889
# 3           cedars_gsa  947 2128
# 4 cedars_omni_old_gwas    4 1215


###############################

dup_ids<-c(as.character(dup$FID1),as.character(dup$FID2))
length(dup_ids)
table(length(dup_ids)==(nrow(dup)*2))
dup_ids<-dup_ids[!duplicated(dup_ids)]
length(dup_ids)
dim(dup)
# [1] 357    9


# CHECK ALSO THAT BOTH SAMPLES HAVE SAME PHENOTYPE AND GENDER, 
# - keep all new gwas samples

rm(data_remove,data_inconsist)
for (i in 1:length(dup_ids)) {
  
  tmp1<-dup[which(dup$FID1==dup_ids[i]),]
  tmp2<-dup[which(dup$FID2==dup_ids[i]),]
  colnames(tmp2)[1:2]<-colnames(tmp2)[2:1]
  
  tmp<-rbind(tmp1,tmp2)
  
  ids_tmp<-c(as.character(tmp$FID1),as.character(tmp$FID2))
  ids_tmp<-ids_tmp[!duplicated(ids_tmp)]
  
  #keep same order as in dup_ids and in all
  ids_tmp<-ids_tmp[match(dup_ids[which(dup_ids %in% ids_tmp)],ids_tmp)]
  
  # get number of possible combinations
  n_possible_combinations<-nrow(permutations(length(ids_tmp), 2, letters[1:length(ids_tmp)]))/2
  
  
  # DOUBLE CHEK THAT WE GET ALL EXPECTED PAIRS
  
  if (nrow(dup[which((dup$FID1 %in% ids_tmp) | (dup$FID2 %in% ids_tmp)),])==n_possible_combinations ) {
    
    data<-as.data.frame(matrix(ncol=1,nrow=length(ids_tmp)))
    data$V1<-ids_tmp
    data<-merge(data,dat[,c("V1","sex","pheno","cohort","F_MISS")],by="V1",all.x=T,sort=F)
    
    
    # OPTION A: samples have different pheno, remove all:
    
    if ( (dim(table(data$sex))!=1 | dim(table(data$pheno))!=1) ) {
      
      if(!exists("data_remove")) {
        data_remove<-data
      } else {
        data_remove<-rbind(data_remove,data)
      }
      
      # keep this also as independent table:
      if(!exists("data_inconsist")){
        jj<-1
        data_inconsist<-data
        data_inconsist$group<-jj
      } else {
        jj<-jj+1
        data$group<-jj
        data_inconsist<-rbind(data_inconsist,data)
      }
    } 
    
    # OPTION B: samples have same pheno
    
    else {
      
      # B.1: keep cedars_gsa:
      
      if ( (any(data$cohort %in% c("cedars_gsa")) )) {
        
        keep_sample<-data$V1[which(data$cohort %in% c("cedars_gsa"))]
        
        # if more than one new pair in new gwas, keep sample with largest call rate:
        if (length(keep_sample)>1) {
          keep_sample<-data[which(data$cohort %in% c("cedars_gsa")),]
          keep_sample<-keep_sample$V1[which(keep_sample$F_MISS==min(keep_sample$F_MISS))][1]# to deal with >1 with same F_miss
        }
        
        if(!exists("data_remove")) {
          data_remove<-data[which(!data$V1 %in% keep_sample),]
        } else {
          data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
        }
        
      } 
      
      # B.2: any other option, keep sample with largest call rate:
      
      else {
        
        keep_sample<-data$V1[which(data$F_MISS==min(data$F_MISS))][1]# to deal with >1 with same F_miss
        
        if(!exists("data_remove")) {
          data_remove<-data[which(!data$V1 %in% keep_sample),]
        } else {
          data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
        }
      }
      
    }
  }
  
  else {
    
    # not all combinations of duplicated pairs found, show list of IDs, to manually inspect issues
    print(i)
    print(paste("Number of expected combinations: ",n_possible_combinations,sep=""))
    print(paste("Number of observed combinations: ",nrow(dup[which((dup$FID1 %in% ids_tmp) | (dup$FID2 %in% ids_tmp)),]),sep=""))
    print(ids_tmp)
    
  }
  
}


dim(data_remove)
# [1] 714   5

data_remove<-data_remove[which(!duplicated(data_remove$V1)),]
dim(data_remove)
# [1] 357  5

table(data_remove$cohort)
# cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas 
#                   52                  122                  183 

data_inconsist<-data_inconsist[which(!duplicated(data_inconsist$V1)),]
dim(data_inconsist)
# [1] 0  6
# data_inconsist


data_remove_tmp<-data_remove[,c(1,1)]
colnames(data_remove_tmp)<-c("FID","IID")

write.table(data_remove,paste(path,"pre_imputation/QC/relatedness/cedars_data_remove",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
write.table(data_remove_tmp,paste(path,"pre_imputation/QC/relatedness/list_cedars_data_remove.tsv",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
system(paste("cp ",path,"pre_imputation/QC/relatedness/list_cedars_data_remove.tsv ~/files_iibdgc/",sep=""))

# write.table(data_inconsist,paste(path,"pre_imputation/QC/relatedness/cedars_data_inconsist",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
write.table(aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary),paste(path,"pre_imputation/QC/relatedness/cedars_case_control_pre_per_study_relatedness",sep=""),col.names=T,row.names=F,sep="\t",quote=F)

write.table(dup,paste(path,"pre_imputation/QC/relatedness/all_cedars_duplicated_pairs.tsv",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
system(paste("cp ",path,"pre_imputation/QC/relatedness/all_cedars_duplicated_pairs.tsv ~/files_iibdgc/",sep=""))

####################################################################################################################################################

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","cedars_gsa")

for (i in 1:length(cohorts)) {
  tmp<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],
                        "_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy.fam",sep=""),head=F)
  tmp$cohort<-as.character(cohorts[i])
  tmp$V1<-as.character(tmp$V1)
  tmp$V2<-as.character(tmp$V2)
  
  if(i==1){
    dat<-tmp
  }else{
    dat<-rbind(dat,tmp)
  }
}

table(dat$cohort)
# cedars_370k_old_gwas cedars_610k_old_gwas           cedars_gsa 
#                  553                  767                 3075 
# cedars_omni_old_gwas 
#                 1036

colnames(dat)[6]<-"pheno"
aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary)
#                Group.1  x.1  x.2
# 1 cedars_370k_old_gwas    0  552
# 2 cedars_610k_old_gwas    0  767
# 3           cedars_gsa    0 2127
# 4 cedars_omni_old_gwas    4 1035

colnames(dat)[5]<-"sex"
aggregate(as.factor(dat$sex), by=list(dat$cohort), FUN=summary)
               Group.1  x.0  x.1  x.2
# 1 cedars_370k_old_gwas    1  297  254
# 2 cedars_610k_old_gwas    3  403  361
# 3           cedars_gsa    6 1030 1091
# 4 cedars_omni_old_gwas    0  517  522

write.table(aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary),paste(path,"pre_imputation/QC/relatedness/cedars_case_control_post_per_study_relatedness",sep=""),col.names=T,row.names=F,sep="\t",quote=F)

##############################
# remove intermediate files:

ls -la ${path}pre_imputation/QC/relatedness/*_subset_cedars*
rm ${path}pre_imputation/QC/relatedness/*_subset_cedars*
  
  
  
  