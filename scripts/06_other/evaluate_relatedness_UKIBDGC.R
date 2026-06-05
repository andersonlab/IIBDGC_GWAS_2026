# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# EVALUATE RELATEDNESS BETWEEN UKIBDGC SAMPLES

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
  
studies=(gwas1 gwas2 all_hce)

path_variants=${path}pre_imputation/QC/relatedness/list_common_variants_among_cohorts_ukibdgc

for i in ${studies[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample \
--allow-no-sex \
--extract ${path_variants} \
--make-bed \
--out ${path}pre_imputation/QC/relatedness/${i}_subset_uk
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
  
  
wc -l ${path}pre_imputation/QC/relatedness/*_subset_uk.fam
# 22614 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/all_hce_subset_uk.fam
# 4679 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/gwas1_subset_uk.fam
# 7777 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/gwas2_subset_uk.fam
# 35070 total
  
wc -l ${path}pre_imputation/QC/relatedness/*_subset_uk.bim
# 36239  OK
  
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
  
cohorts<-c("gwas1","gwas2","all_hce")

cohorts<-cohorts[-1]
  
dat<-matrix(ncol=3,nrow=length(cohorts))
dat<-as.data.frame(dat)
  
for (i in 1:nrow(dat)){
  dat[i,1]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_uk.bed",sep="")
  dat[i,2]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_uk.bim",sep="")
  dat[i,3]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_uk.fam",sep="")
}

write.table(dat,paste(path,"pre_imputation/QC/relatedness/list_cohorts_tomerge_uk.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#############################

# hpc-server
  
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/relatedness/gwas1_subset_uk \
--allow-no-sex \
--merge-list ${path}pre_imputation/QC/relatedness/list_cohorts_tomerge_uk.txt \
--make-bed --out ${path}pre_imputation/QC/relatedness/ukibdgc_merged
  
# 336239 variants and 35070 people pass filters and QC.
# Among remaining phenotypes, 16275 are cases and 18795 are controls.


##############################################################################################################################################
  
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
  
  
# (36877*34904/4)/1E+6
# [1] 321.7887
  
path=/path/to/ibdgwas/IIBDGC/
MEM=1000
  
bsub -J"kg" -M"$MEM" -n8 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path}pre_imputation/QC/relatedness/log/stderr_king_ukibdgc \
-o ${path}pre_imputation/QC/relatedness/log/stdout_king_ukibdgc \
"/path/to/software/./king \
-b ${path}pre_imputation/QC/relatedness/ukibdgc_merged.bed \
--related --cpus 8 --prefix ${path}pre_imputation/QC/relatedness/ukibdgc_merged"
# Job <20399> is submitted to queue <normal>.


#### R

path<-"/path/to/ibdgwas/IIBDGC/"
  
cohorts<-c("gwas1","gwas2","all_hce")
  
  
for (i in 1:length(cohorts)) {
  tmp<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],
                          "_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.fam",sep=""),head=F)
  sample_miss<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],
                                  "_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.imiss",sep=""),head=T)
  tmp<-merge(tmp,sample_miss[,c("IID","F_MISS")],by.x="V2",by.y="IID",all.x=T)
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
# all_hce   gwas1   gwas2 
# 22614    4679    7777 
  
# double check all samples included in exercise
famall<-read.table(paste(path,"pre_imputation/QC/relatedness/ukibdgc_merged.fam",sep=""),head=F)
  
dim(famall)
# [1] 35070     6
dim(dat)
# [1] 35070     8
dim(dat[which(dat$V1 %in% famall$V1),])
# [1] 35070     8
  
colnames(dat)[6]<-"pheno"
aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary)
#   Group.1   x.1   x.2
# 1 all_hce 10447 12167
# 2   gwas1  2932  1747
# 3   gwas2  5416  2361

colnames(dat)[5]<-"sex"
aggregate(as.factor(dat$sex), by=list(dat$cohort), FUN=summary)
#   Group.1   x.0   x.1   x.2
# 1 all_hce    83 10427 12104
# 2   gwas1     2  2122  2555
# 3   gwas2    78  3831  3868

table(dat$V1==dat$V2)
# TRUE 
# 35070

kin<-read.table(paste(path,"pre_imputation/QC/relatedness/ukibdgc_merged.kin0",sep=""),head=T)
table(kin$InfType)
# 2nd    3rd    4th Dup/MZ     FS     PO     UN 
# 78     55      9   2835    206    990     11  

table(cut(kin$Kinship,breaks=c(1,0.354,0.177)))
# (0.177,0.354]     (0.354,1] 
#          1300          2877  

summary(kin$Kinship[which(kin$Kinship>0.354 & kin$InfType!="Dup/MZ")])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.3546  0.3648  0.3756  0.3806  0.3977  0.4140 

dup<-kin[which(kin$Kinship>0.354 & kin$InfType=="Dup/MZ"),c("FID1","FID2","Kinship")]
dim(dup)
# [1] 2835    3


dup<-merge(dup,dat[,c("V1","cohort","sex","pheno","F_MISS")],by.x="FID2",by.y="V1",all.x=T,sort=F)
colnames(dup)[(ncol(dup)-3):ncol(dup)]<-paste(colnames(dup)[(ncol(dup)-3):ncol(dup)],"FID2",sep="_")
  
dup<-merge(dup,dat[,c("V1","cohort","sex","pheno","F_MISS")],by.x="FID1",by.y="V1",all.x=T,sort=F)
colnames(dup)[(ncol(dup)-3):ncol(dup)]<-paste(colnames(dup)[(ncol(dup)-3):ncol(dup)],"FID1",sep="_")
  
x<-table(dup$cohort_FID1,dup$cohort_FID2)
x
#          all_hce gwas1 gwas2
#  all_hce       1   123    56
#  gwas1         0     0  2535
#  gwas2        10   110     0

xx<-dup[,c("cohort_FID1","cohort_FID2")]
yy<-dup[,c("cohort_FID2","cohort_FID1")]
colnames(yy)<-colnames(xx)
  
dim(xx)
# [1] 2835    2
  
xx<-rbind(xx,yy)
x<-table(xx$cohort_FID1,xx$cohort_FID2)
x
#         all_hce gwas1 gwas2
# all_hce       2   123    66
# gwas1       123     0  2645
# gwas2        66  2645     0

write.table(x,paste(path,"pre_imputation/QC/relatedness/ukibdgc_dup_summary_table",sep=""),col.names=T,row.names=T,sep="\t",quote=F)

dup[which(dup$cohort_FID1=="all_hce" & dup$cohort_FID2=="all_hce"),]
#                        FID1                     FID2 Kinship cohort_FID2
# 28 331743_B05_IBD_EX5371959 367989_C03_IBD_Ex5963356   0.474     all_hce
#    sex_FID2 pheno_FID2 F_MISS_FID2 cohort_FID1 sex_FID1 pheno_FID1 F_MISS_FID1
# 28        1          2     0.02913     all_hce        1          2   0.0003665

vec<-c("331743_B05_IBD_EX5371959","367989_C03_IBD_Ex5963356")
kin[which(kin$FID1 %in% vec & kin$FID1 %in% vec),]
#                          FID1                      ID1                     FID2
# 1057 331743_B05_IBD_EX5371959 331743_B05_IBD_EX5371959 367989_C03_IBD_Ex5963356
#                           ID2 N_SNP HetHet IBS0 HetConc HomIBS0 Kinship IBD1Seg
# 1057 367989_C03_IBD_Ex5963356 36149 0.3595    0   0.905       0   0.474  0.7492
#      IBD2Seg PropIBD InfType
# 1057  0.1829  0.5575  Dup/MZ

# GWAS3 only
# kin[which(kin$FID1 %in% vec & kin$FID1 %in% vec),]
# #                          FID1                      ID1                     FID2
# # 1257 331743_B05_IBD_EX5371959 331743_B05_IBD_EX5371959 367989_C03_IBD_Ex5963356
# #                           ID2 N_SNP HetHet  IBS0 HetConc HomIBS0 Kinship
# # 1257 367989_C03_IBD_Ex5963356 420997 0.1874 8e-04  0.7838  0.0156  0.4279
# #      IBD1Seg IBD2Seg PropIBD InfType
# # 1257  0.8105  0.1444  0.5497      FS


# duplicated pairs in GWAS1 and GWAS2:
table(dup[which(dup$cohort_FID1 %in% c("gwas1","gwas2") & dup$cohort_FID2 %in% c("gwas1","gwas2")),"pheno_FID1"],
      dup[which(dup$cohort_FID1 %in% c("gwas1","gwas2") & dup$cohort_FID2 %in% c("gwas1","gwas2")),"pheno_FID2"])
#      1    2
# 1 2640    0
# 2    0    5

# duplicated pairs in GWAS1 and GWAS3:
table(dup[which(dup$cohort_FID1 %in% c("gwas1","all_hce") & dup$cohort_FID2 %in% c("gwas1","all_hce")),"pheno_FID1"],
      dup[which(dup$cohort_FID1 %in% c("gwas1","all_hce") & dup$cohort_FID2 %in% c("gwas1","all_hce")),"pheno_FID2"])
#     1   2
# 1   2   0
# 2   2 120

# duplicated pairs in GWAS2 and GWAS3:
table(dup[which(dup$cohort_FID1 %in% c("gwas2","all_hce") & dup$cohort_FID2 %in% c("gwas2","all_hce")),"pheno_FID1"],
      dup[which(dup$cohort_FID1 %in% c("gwas2","all_hce") & dup$cohort_FID2 %in% c("gwas2","all_hce")),"pheno_FID2"])
#    1  2
# 1  2  0
# 2  2 63



table(dup$sex_FID1,dup$sex_FID2)
#      0    1    2
# 1   27 1360    0
# 2    0    0 1448
(dup[which(dup$sex_FID2=="0" & dup$sex_FID1==1),])
# all gwas2 sex =0, gwas1 sex = 1


dup<-dup[,c("FID1","FID2","Kinship","cohort_FID1","pheno_FID1","F_MISS_FID1","cohort_FID2","pheno_FID2","F_MISS_FID2")]
table(dup$pheno_FID1,dup$pheno_FID2)
#      1    2
# 1 2644    0
# 2    4  187
  
dup_ids<-c(as.character(dup$FID1),as.character(dup$FID2))
length(dup_ids)
table(length(dup_ids)==(nrow(dup)*2))
dup_ids<-dup_ids[!duplicated(dup_ids)]
length(dup_ids)
dim(dup)
# [1] 2835    9


# CHECK ALSO THAT BOTH SAMPLES HAVE SAME PHENOTYPE AND GENDER, 
# - keep all new gwas samples
  
# - for ukbb-ukbb duplicated pairs:
#    - keep sample with pheno data when both samples have the same phenotype, otherwise exclude both
  

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
      
      # make an exeption to deal with sex = 0
      if ( any(data$sex==0) & any(data$sex %in% c(1,2)) & dim(table(data$pheno))==1 ) {

        keep_sample<-data$V1[which(data$sex!=0)]
        
        if(!exists("data_remove")) {
          data_remove<-data[which(!data$V1 %in% keep_sample),]
        } else {
          data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
        }
        
      } else {
        if (!exists("data_remove")) {
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
    } 
    
    # OPTION B: samples have same pheno
    
    else {
      
      # B.1: samples in gwas1 and gwas2, keep samples in gwas1, to keep ca:ctr ratio:
      
      if ( any(data$cohort %in% c("gwas1")) & any(data$cohort %in% c("gwas2")) & !any(data$cohort %in% c("all_hce")) ) {
        
        keep_sample<-data$V1[which(!data$cohort %in% c("gwas2"))]
        
        # if more than one new pair in new gwas, keep sample with largest call rate:
        if (length(keep_sample)>1) {
          keep_sample<-data[which(!data$cohort %in% c("gwas2")),]
          keep_sample<-keep_sample$V1[which(keep_sample$F_MISS==min(keep_sample$F_MISS))][1]# to deal with >1 with same F_miss
        }
        
        if(!exists("data_remove")) {
          data_remove<-data[which(!data$V1 %in% keep_sample),]
        } else {
          data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
        }
      }
      
      
      # B.2: gwas1 or gwas2 AND gwas3, keep gwas1 or gwas2 to keep ca/ctr ratio
      
      else if ( (any(data$cohort %in% c("gwas1")) | any(data$cohort %in% c("gwas2")) ) & any(data$cohort %in% c("all_hce")) ) {
          
        keep_sample<-data$V1[which(!data$cohort %in% c("all_hce"))]
          
        # if more than one new pair in good array, keep sample with largest call rate:
        if (length(keep_sample)>1) {
          keep_sample<-data[which(!data$cohort %in% c("all_hce")),]
          keep_sample<-keep_sample$V1[which(keep_sample$F_MISS==min(keep_sample$F_MISS))][1]# to deal with >1 with same F_miss
        }
          
        if (!exists("data_remove")) {
          data_remove<-data[which(!data$V1 %in% keep_sample),]
        } else {
          data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
        }
        
      }
      
      # B.3: gwas1 or gwas2 or gwas3 only, keep better sample
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
# [1] 5677    5

data_remove<-data_remove[which(!duplicated(data_remove$V1)),]
dim(data_remove)
# [1] 2836    5

table(data_remove$cohort)
# all_hce   gwas1   gwas2 
#     188       3    2645 



data_inconsist<-data_inconsist[which(!duplicated(data_inconsist$V1)),]
dim(data_inconsist)
# [1] 7  6
data_inconsist
#                                    V1 sex pheno  cohort    F_MISS group
# 1             368340_E11_IBD_N5965129   2     2 all_hce 0.0002602     1
# 2                         sample_id   2     1   gwas1 0.0045190     1
# 3                        WTCCCT511185   2     1   gwas2 0.0024860     1
# 7  urn:wtsi:435974_H12_ibdgwas6295779   2     2 all_hce 0.0005837     3
# 8                         sample_id   2     1   gwas1 0.0013280     3
# 9  urn:wtsi:435979_B11_ibdgwas6296281   2     2 all_hce 0.0028530     4
# 10                       WTCCCT542888   2     1   gwas2 0.0003903     4

data_remove_tmp<-data_remove[,c(1,1)]
colnames(data_remove_tmp)<-c("FID","IID")

write.table(data_remove,paste(path,"pre_imputation/QC/relatedness/ukibdgc_data_remove",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
write.table(data_remove_tmp,paste(path,"pre_imputation/QC/relatedness/list_ukiibdgc_data_remove",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
write.table(data_inconsist,paste(path,"pre_imputation/QC/relatedness/ukibdgc_data_inconsist",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
write.table(aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary),paste(path,"pre_imputation/QC/relatedness/ukibdgc_case_control_pre_per_study_relatedness",sep=""),col.names=T,row.names=F,sep="\t",quote=F)


####################################################################################################################################################

# path<-"/path/to/ibdgwas/IIBDGC/"
# 
# cohorts<-c("gwas1","gwas2","all_hce")
# 
# 
# for (i in 1:length(cohorts)) {
#   tmp<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],
#                         "_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy.fam",sep=""),head=F)
#   tmp$cohort<-as.character(cohorts[i])
#   tmp$V1<-as.character(tmp$V1)
#   tmp$V2<-as.character(tmp$V2)
#   
#   if(i==1){
#     dat<-tmp
#   }else{
#     dat<-rbind(dat,tmp)
#   }
# }
# 
# table(dat$cohort)
# # all_hce   gwas1   gwas2 
# # 22426    4676    5132 
# 
# colnames(dat)[6]<-"pheno"
# aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary)
# #   Group.1   x.1   x.2
# # 1 all_hce 10444 11982
# # 2   gwas1  2929  1747
# # 3   gwas2  2776  2356
# 
# colnames(dat)[5]<-"sex"
# # aggregate(as.factor(dat$sex), by=list(dat$cohort), FUN=summary)
# #   Group.1   x.0   x.1   x.2
# # 1 all_hce    83 10350 11993
# # 2   gwas1     2  2121  2553
# # 3   gwas2    51  2550  2531
# 
# 
# write.table(aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary),paste(path,"pre_imputation/QC/relatedness/ukibdgc_case_control_post_per_study_relatedness",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
