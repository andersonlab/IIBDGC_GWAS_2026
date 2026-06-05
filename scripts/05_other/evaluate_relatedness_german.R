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
  
studies=(kiel_austria_sibdcs_gsa german_illu550_old_gwas german_affy6_old_gwas)
  
path_variants=${path}pre_imputation/QC/relatedness/list_common_variants_among_cohorts_german
  
for i in ${studies[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample \
--extract ${path_variants} \
--allow-no-sex \
--make-bed \
--out ${path}pre_imputation/QC/relatedness/${i}_subset_german
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

wc -l ${path}pre_imputation/QC/relatedness/*_subset_german.fam
# 2787 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/german_affy6_old_gwas_subset_german.fam
# 1621 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/german_illu550_old_gwas_subset_german.fam
# 13906 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/kiel_austria_sibdcs_gsa_subset_german.fam
# 18314 total

wc -l ${path}pre_imputation/QC/relatedness/*_subset_german.bim
# 28009  OK

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

cohorts<-c("kiel_austria_sibdcs_gsa","german_illu550_old_gwas","german_affy6_old_gwas")

cohorts<-cohorts[-1]

dat<-matrix(ncol=3,nrow=length(cohorts))
dat<-as.data.frame(dat)

for (i in 1:nrow(dat)){
  dat[i,1]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_german.bed",sep="")
  dat[i,2]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_german.bim",sep="")
  dat[i,3]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_german.fam",sep="")
}

write.table(dat,paste(path,"pre_imputation/QC/relatedness/list_cohorts_tomerge_german.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#############################

# edit german_illu550_old_gwas fam files so IDs are not duplicated (all of them in affy6) and confirm all samples are in affy6, and merge again:

fam<-read.table(paste(path,"pre_imputation/QC/relatedness/german_illu550_old_gwas_subset_german.fam",sep=""),head=F)
fam$V1<-paste(fam$V1,"_illu550",sep="")
fam$V2<-paste(fam$V2,"_illu550",sep="")
write.table(fam,paste(path,"pre_imputation/QC/relatedness/german_illu550_old_gwas_subset_german.fam",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

#############################

# hpc-server

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/relatedness/kiel_austria_sibdcs_gsa_subset_german \
--allow-no-sex \
--merge-list ${path}pre_imputation/QC/relatedness/list_cohorts_tomerge_german.txt \
--make-bed --out ${path}pre_imputation/QC/relatedness/german_merged

# 28013 variants and 18314 people pass filters and QC.
# Among remaining phenotypes, 10968 are cases and 7346 are controls.



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


# (28014*18108/4)/1E+6
# [1] 126.8194

path=/path/to/ibdgwas/IIBDGC/
MEM=1000

bsub -J"kg" -M"$MEM" -n8 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path}pre_imputation/QC/relatedness/log/stderr_king_german \
-o ${path}pre_imputation/QC/relatedness/log/stdout_king_german \
"/path/to/software/./king \
-b ${path}pre_imputation/QC/relatedness/german_merged.bed \
--related --cpus 8 --prefix ${path}pre_imputation/QC/relatedness/german_merged"
# Job <257841> is submitted to queue <normal>.


# #### R
# 
# path<-"/path/to/ibdgwas/IIBDGC/"
# 
# cohorts<-c("kiel_austria_sibdcs_gsa","german_illu550_old_gwas","german_affy6_old_gwas")
# 
# for (i in 1:length(cohorts)) {
#   tmp<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],
#                         "_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.fam",sep=""),head=F)
#   sample_miss<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],
#                                 "_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample.imiss",sep=""),head=T)
#   tmp<-merge(tmp,sample_miss[,c("IID","F_MISS")],by.x="V2",by.y="IID",all.x=T)
#   tmp$cohort<-as.character(cohorts[i])
#   tmp$V1<-as.character(tmp$V1)
#   tmp$V2<-as.character(tmp$V2)
#   
#   if(cohorts[i]=="german_illu550_old_gwas") {
#     tmp$V1<-paste(tmp$V1,"_illu550",sep="")
#     tmp$V2<-paste(tmp$V2,"_illu550",sep="")
#   }
#   
#   if(i==1){
#     dat<-tmp
#   }else{
#     dat<-rbind(dat,tmp)
#   }
# }
# 
# table(dat$cohort)
# # german_affy6_old_gwas german_illu550_old_gwas kiel_austria_sibdcs_gsa 
# # 787                    1621                   13906
# 
# # double check all samples included in exercise
# famall<-read.table(paste(path,"pre_imputation/QC/relatedness/german_merged.fam",sep=""),head=F)
# 
# dim(famall)
# # [1] 18314     6
# dim(dat)
# # [1] 18314     8
# dim(dat[which(dat$V1 %in% famall$V1),])
# # [1] 18314     8
# 
# colnames(dat)[6]<-"pheno"
# aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary)
# #                   Group.1  x.1  x.2
# # 1   german_affy6_old_gwas 1750 1037
# # 2 german_illu550_old_gwas 1141  480
# # 3 kiel_austria_sibdcs_gsa 4455 9451
# 
# colnames(dat)[5]<-"sex"
# aggregate(as.factor(dat$sex), by=list(dat$cohort), FUN=summary)
# # Group.1  x.0  x.1  x.2
# # 1   german_affy6_old_gwas    5 1358 1424
# # 2 german_illu550_old_gwas   10  636  975
# # 3 kiel_austria_sibdcs_gsa   12 6896 6998
# 
# table(dat$V1==dat$V2)
# # TRUE 
# # 18314
# 
# kin<-read.table(paste(path,"pre_imputation/QC/relatedness/german_merged.kin0",sep=""),head=T)
# table(kin$InfType)
# # 2nd Dup/MZ     FS     PO 
# # 6   1732    193    419 
# 
# table(cut(kin$Kinship,breaks=c(1,0.354,0.177)))
# # (0.177,0.354]     (0.354,1] 
# #           597          1732 
# 
# dup<-kin[which(kin$Kinship>0.354),c("FID1","FID2","Kinship")]
# summary(dup$Kinship)
# # #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 0.4377  0.4989  0.4994  0.4989  0.4999  0.5000 
# 
# dup<-merge(dup,dat[,c("V1","cohort","sex","pheno","F_MISS")],by.x="FID2",by.y="V1",all.x=T,sort=F)
# colnames(dup)[(ncol(dup)-3):ncol(dup)]<-paste(colnames(dup)[(ncol(dup)-3):ncol(dup)],"FID2",sep="_")
# 
# dup<-merge(dup,dat[,c("V1","cohort","sex","pheno","F_MISS")],by.x="FID1",by.y="V1",all.x=T,sort=F)
# colnames(dup)[(ncol(dup)-3):ncol(dup)]<-paste(colnames(dup)[(ncol(dup)-3):ncol(dup)],"FID1",sep="_")
# 
# x<-table(dup$cohort_FID1,dup$cohort_FID2)
# x
# # german_affy6_old_gwas german_illu550_old_gwas
# # german_affy6_old_gwas                       0                    1111
# # german_illu550_old_gwas                     5                       0
# # kiel_austria_sibdcs_gsa                    73                     391
# # 
# # kiel_austria_sibdcs_gsa
# # german_affy6_old_gwas                        82
# # german_illu550_old_gwas                      70
# # kiel_austria_sibdcs_gsa                       0
# 
# 
# xx<-dup[,c("cohort_FID1","cohort_FID2")]
# yy<-dup[,c("cohort_FID2","cohort_FID1")]
# colnames(yy)<-colnames(xx)
# 
# dim(xx)
# # [1] 1732    2
# 
# xx<-rbind(xx,yy)
# x<-table(xx$cohort_FID1,xx$cohort_FID2)
# x
# # german_affy6_old_gwas german_illu550_old_gwas
# # german_affy6_old_gwas                       0                    1116
# # german_illu550_old_gwas                  1116                       0
# # kiel_austria_sibdcs_gsa                   155                     461
# # 
# # kiel_austria_sibdcs_gsa
# # german_affy6_old_gwas                       155
# # german_illu550_old_gwas                     461
# # kiel_austria_sibdcs_gsa                       0
# 
# 
# write.table(x,paste(path,"pre_imputation/QC/relatedness/german_dup_summary_table",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
# 
# 
# 
# table(dup$sex_FID1,dup$sex_FID2)
# #     0   1   2
# # 0   1   0   0
# # 1   6 737   0
# # 2   4   0 984
# 
# dup<-dup[,c("FID1","FID2","Kinship","cohort_FID1","pheno_FID1","F_MISS_FID1","cohort_FID2","pheno_FID2","F_MISS_FID2")]
# table(dup$pheno_FID1,dup$pheno_FID2)
# #      1    2
# # 1 1266    2
# # 2    0  464
# 
# 
# # duplicated pairs in Illu and Affy:
# table(dup[which(dup$cohort_FID1 %in% c("german_illu550_old_gwas","german_affy6_old_gwas") & dup$cohort_FID2 %in% c("german_illu550_old_gwas","german_affy6_old_gwas")),"pheno_FID1"],
#       dup[which(dup$cohort_FID1 %in% c("german_illu550_old_gwas","german_affy6_old_gwas") & dup$cohort_FID2 %in% c("german_illu550_old_gwas","german_affy6_old_gwas")),"pheno_FID2"])
# #      1    2
# # 1 1115    1
# 
# # duplicated pairs in illu and kiel:
# table(dup[which(dup$cohort_FID1 %in% c("german_illu550_old_gwas","kiel_austria_sibdcs_gsa") & dup$cohort_FID2 %in% c("german_illu550_old_gwas","kiel_austria_sibdcs_gsa")),"pheno_FID1"],
#       dup[which(dup$cohort_FID1 %in% c("german_illu550_old_gwas","kiel_austria_sibdcs_gsa") & dup$cohort_FID2 %in% c("german_illu550_old_gwas","kiel_austria_sibdcs_gsa")),"pheno_FID2"])
# #     1   2
# # 1  70   0
# # 2   0 391
# 
# # duplicated pairs in affy6 and kiel:
# table(dup[which(dup$cohort_FID1 %in% c("german_affy6_old_gwas","kiel_austria_sibdcs_gsa") & dup$cohort_FID2 %in% c("german_affy6_old_gwas","kiel_austria_sibdcs_gsa")),"pheno_FID1"],
#       dup[which(dup$cohort_FID1 %in% c("german_affy6_old_gwas","kiel_austria_sibdcs_gsa") & dup$cohort_FID2 %in% c("german_affy6_old_gwas","kiel_austria_sibdcs_gsa")),"pheno_FID2"])
# #    1  2
# # 1 81  1
# # 2  0 73
# 
# 
# # how many Illu cases in Kiel?
# test<-dup[which(dup$cohort_FID1 %in% c("german_illu550_old_gwas","kiel_austria_sibdcs_gsa") & dup$cohort_FID2 %in% c("german_illu550_old_gwas","kiel_austria_sibdcs_gsa")),]
# test1<-test[which( (test$cohort_FID1 %in% c("german_illu550_old_gwas")) & (test$pheno_FID1==2) ),c("FID1")]
# test2<-test[which( (test$cohort_FID2 %in% c("german_illu550_old_gwas")) & (test$pheno_FID2==2) ),c("FID2")]
# testall<-c(test1,test2)
# testall<-testall[!duplicated(testall)]
# length(testall)
# # [1] 391
# 
# # how many Illu controls in Affy6?
# test<-dup[which(dup$cohort_FID1 %in% c("german_illu550_old_gwas","german_affy6_old_gwas") & dup$cohort_FID2 %in% c("german_illu550_old_gwas","german_affy6_old_gwas")),]
# test1<-test[which( (test$cohort_FID1 %in% c("german_illu550_old_gwas")) & (test$pheno_FID1==1) ),c("FID1")]
# test2<-test[which( (test$cohort_FID2 %in% c("german_illu550_old_gwas")) & (test$pheno_FID2==1) ),c("FID2")]
# testall<-c(test1,test2)
# testall<-testall[!duplicated(testall)]
# length(testall)
# # [1] 1104
# 
# ##############################
# # ARE ALL ILLU550 IN AFFY6.0
# 
# aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary)
# #                   Group.1  x.1  x.2
# # 1   german_affy6_old_gwas 1749 1033
# # 2 german_illu550_old_gwas 1131  480
# # 3 kiel_austria_sibdcs_gsa 4454 9261
# 
# 
# ilu<-dat[which(dat$cohort=="german_illu550_old_gwas"),]
# affy<-dat[which(dat$cohort=="german_affy6_old_gwas"),]
# kiel<-dat[which(dat$cohort=="kiel_austria_sibdcs_gsa"),]
# 
# ilu$ID<-gsub("_illu550","",ilu$V1)
# 
# dim(ilu)
# # [1] 1621    9
# dim(ilu[which(ilu$ID %in% affy$V1),])
# # [1] 1108    9
# dim(ilu[which(ilu$ID %in% kiel$V1),])
# # [1] 0 9
# 
# 
# table(ilu[which(ilu$ID %in% affy$V1),"pheno"])
# # 1 
# # 1108
# 
# table(affy$pheno)
# # 1    2 
# # 1750 1037 
# 
# table(ilu$pheno)
# # 1    2 
# # 1141  480 
# 
# ###############################
# 
# 
# dup_ids<-c(as.character(dup$FID1),as.character(dup$FID2))
# length(dup_ids)
# table(length(dup_ids)==(nrow(dup)*2))
# dup_ids<-dup_ids[!duplicated(dup_ids)]
# length(dup_ids)
# dim(dup)
# # [1] 1732    9
# 
# 
# # CHECK ALSO THAT BOTH SAMPLES HAVE SAME PHENOTYPE AND GENDER, 
# # - keep all new gwas samples
# 
# # - for ukbb-ukbb duplicated pairs:
# #    - keep sample with pheno data when both samples have the same phenotype, otherwise exclude both
# 
# 
# for (i in 1:length(dup_ids)) {
#   
#   tmp1<-dup[which(dup$FID1==dup_ids[i]),]
#   tmp2<-dup[which(dup$FID2==dup_ids[i]),]
#   colnames(tmp2)[1:2]<-colnames(tmp2)[2:1]
#   
#   tmp<-rbind(tmp1,tmp2)
#   
#   ids_tmp<-c(as.character(tmp$FID1),as.character(tmp$FID2))
#   ids_tmp<-ids_tmp[!duplicated(ids_tmp)]
#   
#   #keep same order as in dup_ids and in all
#   ids_tmp<-ids_tmp[match(dup_ids[which(dup_ids %in% ids_tmp)],ids_tmp)]
#   
#   # get number of possible combinations
#   n_possible_combinations<-nrow(permutations(length(ids_tmp), 2, letters[1:length(ids_tmp)]))/2
#   
#   
#   # DOUBLE CHEK THAT WE GET ALL EXPECTED PAIRS
#   
#   if (nrow(dup[which((dup$FID1 %in% ids_tmp) | (dup$FID2 %in% ids_tmp)),])==n_possible_combinations ) {
#     
#     data<-as.data.frame(matrix(ncol=1,nrow=length(ids_tmp)))
#     data$V1<-ids_tmp
#     data<-merge(data,dat[,c("V1","sex","pheno","cohort","F_MISS")],by="V1",all.x=T,sort=F)
#     
#     
#     # OPTION A: samples have different pheno, remove all:
#     
#     if ( (dim(table(data$sex))!=1 | dim(table(data$pheno))!=1) ) {
#       
#       # make an exeption to deal with sex = 0
#       if ( any(data$sex==0) & any(data$sex %in% c(1,2)) & dim(table(data$pheno))==1 ) {
#         
#         keep_sample<-data$V1[which(data$sex!=0)]
#         
#         if(!exists("data_remove")) {
#           data_remove<-data[which(!data$V1 %in% keep_sample),]
#         } else {
#           data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
#         }
#         
#       } else {
#         if (!exists("data_remove")) {
#           data_remove<-data
#         } else {
#           data_remove<-rbind(data_remove,data)
#         }
#         
#         # keep this also as independent table:
#         if(!exists("data_inconsist")){
#           jj<-1
#           data_inconsist<-data
#           data_inconsist$group<-jj
#         } else {
#           jj<-jj+1
#           data$group<-jj
#           data_inconsist<-rbind(data_inconsist,data)
#         }
#       }
#     } 
#     
#     
#     # OPTION B: samples have same pheno
#     
#     else {
#       
#       # B.1: samples in german_affy6_old_gwas and german_illu550_old_gwas, keep samples in german_affy6_old_gwas, to keep ca:ctr ratio:
#       
#       if ( any(data$cohort %in% c("german_illu550_old_gwas")) & any(data$cohort %in% c("german_affy6_old_gwas")) & !any(data$cohort %in% c("kiel_austria_sibdcs_gsa")) ) {
#         
#         keep_sample<-data$V1[which(!data$cohort %in% c("german_illu550_old_gwas"))]
#         
#         # if more than one new pair in new gwas, keep sample with largest call rate:
#         if (length(keep_sample)>1) {
#           keep_sample<-data[which(!data$cohort %in% c("german_illu550_old_gwas")),]
#           keep_sample<-keep_sample$V1[which(keep_sample$F_MISS==min(keep_sample$F_MISS))][1]# to deal with >1 with same F_miss
#         }
#         
#         if(!exists("data_remove")) {
#           data_remove<-data[which(!data$V1 %in% keep_sample),]
#         } else {
#           data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
#         }
#       }
#       
#       
#       # B.2: german_affy6_old_gwas or german_illu550_old_gwas AND kiel_austria_sibdcs_gsa, keep german_affy6_old_gwas or german_illu550_old_gwas to keep ca/ctr ratio
#       
#       else if ( (any(data$cohort %in% c("german_illu550_old_gwas")) | any(data$cohort %in% c("german_affy6_old_gwas")) ) & any(data$cohort %in% c("kiel_austria_sibdcs_gsa")) ) {
#         
#         keep_sample<-data$V1[which(!data$cohort %in% c("kiel_austria_sibdcs_gsa"))]
#         
#         # if more than one new pair in good array, keep sample with largest call rate:
#         if (length(keep_sample)>1) {
#           keep_sample<-data[which(!data$cohort %in% c("kiel_austria_sibdcs_gsa")),]
#           keep_sample<-keep_sample$V1[which(keep_sample$F_MISS==min(keep_sample$F_MISS))][1]# to deal with >1 with same F_miss
#         }
#         
#         if(!exists("data_remove")) {
#           data_remove<-data[which(!data$V1 %in% keep_sample),]
#         } else {
#           data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
#         }
#         
#       }
#       
#       # B.3: german_affy6_old_gwas or german_illu550_old_gwas or kiel_austria_sibdcs_gsa only, keep better sample
#       else {
#         
#         keep_sample<-data$V1[which(data$F_MISS==min(data$F_MISS))][1]# to deal with >1 with same F_miss
#         
#         if(!exists("data_remove")) {
#           data_remove<-data[which(!data$V1 %in% keep_sample),]
#         } else {
#           data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
#         }
#       }
#     }
#   }
#   
#   else {
#     
#     # not all combinations of duplicated pairs found, show list of IDs, to manually inspect issues
#     print(i)
#     print(paste("Number of expected combinations: ",n_possible_combinations,sep=""))
#     print(paste("Number of observed combinations: ",nrow(dup[which((dup$FID1 %in% ids_tmp) | (dup$FID2 %in% ids_tmp)),]),sep=""))
#     print(ids_tmp)
#     
#   }
#   
# }
# 
# dim(data_remove)
# # [1] 3467    5
# 
# data_remove<-data_remove[which(!duplicated(data_remove$V1)),]
# dim(data_remove)
# # [1] 1666    5
# 
# table(data_remove$cohort)
# # german_affy6_old_gwas german_illu550_old_gwas kiel_austria_sibdcs_gsa 
# #                    67                    1050                     549
# 
# 
# 
# data_inconsist<-data_inconsist[which(!duplicated(data_inconsist$V1)),]
# dim(data_inconsist)
# # [1] 3  6
# data_inconsist
# #                      V1 sex pheno                  cohort    F_MISS group
# # 1     sample_id   1     1   german_affy6_old_gwas 0.0036120     1
# # 2             sample_id   1     2 kiel_austria_sibdcs_gsa 0.0116400     1
# # 3 sample_id   1     2 german_illu550_old_gwas 0.0002305     1
# 
# data_remove$V1<-gsub("_illu550","",data_remove$V1)
# 
# 
# ## edit the list of samples to remove so we can keep good Ca:ctr ratio:
# 
# data_remove_tmp<-data_remove[which(data_remove$cohort=="german_illu550_old_gwas"),]
# dim(data_remove_tmp)
# # [1] 1050    5
# 
# data_remove_tmp<-data_remove_tmp[which(data_remove_tmp$V1 %in% affy$V1),]
# dim(data_remove_tmp)
# # [1] 1033    5
# 
# # select a number of random elements to keep
# list_ids_to_keep_illu<-sample(data_remove_tmp$V1, 450)
# 
# table(data_remove$cohort)
# # german_affy6_old_gwas german_illu550_old_gwas kiel_austria_sibdcs_gsa 
# # 67                    1050                     549 
# data_remove$cohort[which(data_remove$V1 %in% list_ids_to_keep_illu)]<-"german_affy6_old_gwas"
# table(data_remove$cohort)
# # german_affy6_old_gwas german_illu550_old_gwas kiel_austria_sibdcs_gsa 
# #                   517                     600                     549 
# 
# data_remove_tmp1<-data_remove[which(data_remove$cohort %in% c("german_affy6_old_gwas","kiel_austria_sibdcs_gsa")),c(1,1)]
# colnames(data_remove_tmp1)<-c("FID","IID")
# 
# data_remove_tmp2<-data_remove[which(data_remove$cohort %in% c("german_illu550_old_gwas")),c(1,1)]
# colnames(data_remove_tmp2)<-c("FID","IID")
# 
# 
# write.table(data_remove,paste(path,"pre_imputation/QC/relatedness/german_data_remove",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
# write.table(data_remove_tmp1,paste(path,"pre_imputation/QC/relatedness/list_german_affy_kiel_data_remove",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
# write.table(data_remove_tmp2,paste(path,"pre_imputation/QC/relatedness/list_german_illu_data_remove",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
# 
# write.table(data_inconsist,paste(path,"pre_imputation/QC/relatedness/german_data_inconsist",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
# write.table(aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary),paste(path,"pre_imputation/QC/relatedness/german_case_control_pre_per_study_relatedness",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
# 
# 
# ####################################################################################################################################################
# 
# path<-"/path/to/ibdgwas/IIBDGC/"
# 
# cohorts<-c("kiel_austria_sibdcs_gsa","german_illu550_old_gwas","german_affy6_old_gwas")
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
# # german_affy6_old_gwas german_illu550_old_gwas kiel_austria_sibdcs_gsa 
# #                 2270                    1021                   13178 
# 
# colnames(dat)[6]<-"pheno"
# aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary)
# #                   Group.1  x.1  x.2
# # 1   german_affy6_old_gwas 1233 1037
# # 2 german_illu550_old_gwas  542  479
# # 3 kiel_austria_sibdcs_gsa 4371 8807
# 
# colnames(dat)[5]<-"sex"
# aggregate(as.factor(dat$sex), by=list(dat$cohort), FUN=summary)
# # 
# 
# write.table(aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary),paste(path,"pre_imputation/QC/relatedness/german_case_control_post_per_study_relatedness",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
# 
# 
# ##############################################################################################################################################

### WHAT HAPPENS IF WE EXCLUDE ILLU550 - final option

#### R

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("kiel_austria_sibdcs_gsa","german_affy6_old_gwas")

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
# german_affy6_old_gwas kiel_austria_sibdcs_gsa 
# 2787                   13907

# double check all samples included in exercise
famall<-read.table(paste(path,"pre_imputation/QC/relatedness/german_merged.fam",sep=""),head=F)

dim(famall)
# [1] 18315     6
dim(dat)
# [1] 16694     8
dim(dat[which(dat$V1 %in% famall$V1),])
# [1] 16694     8

colnames(dat)[6]<-"pheno"
aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary)
#                   Group.1  x.1  x.2
# 1   german_affy6_old_gwas 1750 1037
# 2 kiel_austria_sibdcs_gsa 4455 9452

colnames(dat)[5]<-"sex"
aggregate(as.factor(dat$sex), by=list(dat$cohort), FUN=summary)
#
table(dat$V1==dat$V2)
# TRUE 
# 16694

kin<-read.table(paste(path,"pre_imputation/QC/relatedness/german_merged.kin0",sep=""),head=T)
table(kin$InfType)

# 2nd Dup/MZ     FS     PO 
#  5   1733    194    419

table(cut(kin$Kinship,breaks=c(1,0.354,0.177)))
# (0.177,0.354]     (0.354,1] 
# 618          1733

dup<-kin[which(kin$Kinship>0.354 & kin$InfType=="Dup/MZ"),c("FID1","FID2","Kinship")]

dup<-merge(dup,dat[,c("V1","cohort","sex","pheno","F_MISS")],by.x="FID2",by.y="V1",all.x=T,sort=F)
colnames(dup)[(ncol(dup)-3):ncol(dup)]<-paste(colnames(dup)[(ncol(dup)-3):ncol(dup)],"FID2",sep="_")

dup<-merge(dup,dat[,c("V1","cohort","sex","pheno","F_MISS")],by.x="FID1",by.y="V1",all.x=T,sort=F)
colnames(dup)[(ncol(dup)-3):ncol(dup)]<-paste(colnames(dup)[(ncol(dup)-3):ncol(dup)],"FID1",sep="_")


dup<-dup[which(dup$cohort_FID1!="german_illu550_old_gwas" & dup$cohort_FID2!="german_illu550_old_gwas"),]


x<-table(dup$cohort_FID1,dup$cohort_FID2)
x
#                         german_affy6_old_gwas kiel_austria_sibdcs_gsa
# german_affy6_old_gwas                       0                      82
# kiel_austria_sibdcs_gsa                    73                       1

xx<-dup[,c("cohort_FID1","cohort_FID2")]
yy<-dup[,c("cohort_FID2","cohort_FID1")]
colnames(yy)<-colnames(xx)

dim(xx)
# [1] 156    2

xx<-rbind(xx,yy)
x<-table(xx$cohort_FID1,xx$cohort_FID2)
x
#                         german_affy6_old_gwas kiel_austria_sibdcs_gsa
# german_affy6_old_gwas                       0                     155
# kiel_austria_sibdcs_gsa                   155                       2


write.table(x,paste(path,"pre_imputation/QC/relatedness/german_dup_summary_table_2",sep=""),col.names=T,row.names=T,sep="\t",quote=F)



table(dup$sex_FID1,dup$sex_FID2)
#    0  1  2
# 0  1  0  0
# 1  0 89  0
# 2  0  0 66

dup<-dup[,c("FID1","FID2","Kinship","cohort_FID1","pheno_FID1","F_MISS_FID1","cohort_FID2","pheno_FID2","F_MISS_FID2")]
table(dup$pheno_FID1,dup$pheno_FID2)
#    1  2
# 1 81  1
# 2  0 74


dup_ids<-c(as.character(dup$FID1),as.character(dup$FID2))
length(dup_ids)
table(length(dup_ids)==(nrow(dup)*2))
dup_ids<-dup_ids[!duplicated(dup_ids)]
length(dup_ids)
dim(dup)
# [1] 156    9


# CHECK ALSO THAT BOTH SAMPLES HAVE SAME PHENOTYPE AND GENDER, 
# - keep all new gwas samples

# - for ukbb-ukbb duplicated pairs:
#    - keep sample with pheno data when both samples have the same phenotype, otherwise exclude both


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
      
      # print("OLD GWAS - NEW GWAS pairs with different sex/pheno:")
      # print(data)
      
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
      
      # B.1: keep kiel_austria_sibdcs_gsa:
      
      keep_sample<-data$V1[which(data$cohort %in% c("kiel_austria_sibdcs_gsa"))]
        
        # if more than one new pair in new gwas, keep sample with largest call rate:
        if (length(keep_sample)>1) {
          keep_sample<-data[which(data$cohort %in% c("kiel_austria_sibdcs_gsa")),]
          keep_sample<-keep_sample$V1[which(keep_sample$F_MISS==min(keep_sample$F_MISS))][1]# to deal with >1 with same F_miss
        }
        
        if(!exists("data_remove")) {
          data_remove<-data[which(!data$V1 %in% keep_sample),]
        } else {
          data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
        }
      }
    }
  }
}

dim(data_remove)
# [1] 314    5

data_remove<-data_remove[which(!duplicated(data_remove$V1)),]
dim(data_remove)
# [1] 157   5

table(data_remove$cohort)
# german_affy6_old_gwas kiel_austria_sibdcs_gsa 
#                    155                      2 



data_inconsist<-data_inconsist[which(!duplicated(data_inconsist$V1)),]
dim(data_inconsist)
# [1] 2  6
data_inconsist
#                      V1 sex pheno                  cohort    F_MISS group
# 1     sample_id   1     1   german_affy6_old_gwas 0.0036120     1
# 2             sample_id   1     2 kiel_austria_sibdcs_gsa 0.0116400     1


data_remove_tmp<-data_remove[,c(1,1)]
colnames(data_remove_tmp)<-c("FID","IID")


write.table(data_remove,paste(path,"pre_imputation/QC/relatedness/german_data_remove_2",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
write.table(data_remove_tmp,paste(path,"pre_imputation/QC/relatedness/list_german_affy_kiel_data_remove_2",sep=""),col.names=T,row.names=F,sep="\t",quote=F)

write.table(data_inconsist,paste(path,"pre_imputation/QC/relatedness/german_data_inconsist_2",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
write.table(aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary),paste(path,"pre_imputation/QC/relatedness/german_case_control_pre_per_study_relatedness_2",sep=""),col.names=T,row.names=F,sep="\t",quote=F)

####################################################################################################################################################


####################################################################################################################################################

# path<-"/path/to/ibdgwas/IIBDGC/"
# 
# cohorts<-c("kiel_austria_sibdcs_gsa","german_affy6_old_gwas")
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
# # german_affy6_old_gwas kiel_austria_sibdcs_gsa 
# # 2632                   13905 
# 
# colnames(dat)[6]<-"pheno"
# aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary)
# #                   Group.1  x.1  x.2
# # 1   german_affy6_old_gwas 1668  964
# # 2 kiel_austria_sibdcs_gsa 4455 9450
# 
# colnames(dat)[5]<-"sex"
# aggregate(as.factor(dat$sex), by=list(dat$cohort), FUN=summary)
# #                  Group.1  x.0  x.1  x.2
# # 1   german_affy6_old_gwas    4 1270 1358
# # 2 kiel_austria_sibdcs_gsa   12 6895 6998
# 
# 
# write.table(aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary),paste(path,"pre_imputation/QC/relatedness/german_case_control_post_per_study_relatedness_2",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
# 

## goign for second option - easier later analysis





