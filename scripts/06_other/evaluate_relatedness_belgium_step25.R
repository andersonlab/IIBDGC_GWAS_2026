# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# EVALUATE RELATEDNESS BETWEEN BELGIUM SAMPLES

###########################
# 1.- SELECT THE VARIANTS #
###########################

# SEE SCRIPT:
# create_list_common_variants_among_cohorts

########################
# 2.- EXTRACT VARIANTS #
########################

# hpc-server

path_gwas=/path/to/ibdgwas/IIBDGC/
  
studies=(belgium_inf1_old_gwas belgium_inf2_old_gwas belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa franchimont_gsa vermeire_gsa)
path_variants=${path_gwas}pre_imputation/QC/relatedness/list_common_variants_among_cohorts_step25_belgium
  
for i in ${studies[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample \
--allow-no-sex \
--extract ${path_variants} \
--make-bed \
--out ${path_gwas}pre_imputation/QC/relatedness/${i}_subset_belgium_step25
done
  
# estimate missingness for later step
for i in ${studies[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample \
--allow-no-sex \
--missing \
--out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample
done


wc -l ${path_gwas}pre_imputation/QC/relatedness/*_subset_belgium_step25.fam
# 1484 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_franchimont_gsa_subset_belgium_step25.fam
# 1406 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_inf1_old_gwas_subset_belgium_step25.fam
# 271 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_inf2_old_gwas_subset_belgium_step25.fam
# 1505 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_louis_gsa_subset_belgium_step25.fam
# 3965 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_vermeire_gsa_subset_belgium_step25.fam
# 2598 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/franchimont_gsa_subset_belgium_step25.fam
# 4645 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/vermeire_gsa_subset_belgium_step25.fam
# 15874 total

wc -l ${path_gwas}pre_imputation/QC/relatedness/*_subset_belgium_step25.bim
# 58624  OK


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

cohorts<-c("belgium_inf1_old_gwas","belgium_inf2_old_gwas","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa","franchimont_gsa","vermeire_gsa")

cohorts<-cohorts[-1]

dat<-matrix(ncol=3,nrow=length(cohorts))
dat<-as.data.frame(dat)

for (i in 1:nrow(dat)){
  dat[i,1]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_belgium_step25.bed",sep="")
  dat[i,2]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_belgium_step25.bim",sep="")
  dat[i,3]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_belgium_step25.fam",sep="")
}

write.table(dat,paste(path,"pre_imputation/QC/relatedness/list_cohorts_tomerge_belgium_step25.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#############################

# hpc-server

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/relatedness/belgium_inf1_old_gwas_subset_belgium_step25 \
--allow-no-sex \
--merge-list ${path_gwas}pre_imputation/QC/relatedness/list_cohorts_tomerge_belgium_step25.txt \
--make-bed --out ${path_gwas}pre_imputation/QC/relatedness/belgium_merged_step25

# 58624 variants and 15874 people pass filters and QC.
# Among remaining phenotypes, 10698 are cases and 5176 are controls..



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


# (58635*11341/4)/1E+6
# [1] 166.2449

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=1000

bsub -J"kg" -M"$MEM" -n8 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/relatedness/logs/stderr_king_belgium_step25 \
-o ${path_gwas}pre_imputation/QC/relatedness/logs/stdout_king_belgium_step25 \
"/path/to/software/./king \
-b ${path_gwas}pre_imputation/QC/relatedness/belgium_merged_step25.bed \
--related --cpus 8 --prefix ${path_gwas}pre_imputation/QC/relatedness/belgium_merged_step25"
# Job <67672> is submitted to queue <normal>.

####################################################################################################################################################

####################################################################################################################################################

#### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)


path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("belgium_inf1_old_gwas","belgium_inf2_old_gwas","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa","franchimont_gsa","vermeire_gsa")

for (i in 1:length(cohorts)) {
  tmp<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],
                        "_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.fam",sep=""),head=F)
  sample_miss<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],
                                "_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.imiss",sep=""),head=T)
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
# belgium_franchimont_gsa   belgium_inf1_old_gwas   belgium_inf2_old_gwas 
# 1484                    1406                     271 
# belgium_louis_gsa    belgium_vermeire_gsa         franchimont_gsa 
# 1505                    3965                    2598 
# vermeire_gsa 
# 4645 

# double check all samples included in exercise
famall<-read.table(paste(path,"pre_imputation/QC/relatedness/belgium_merged_step25.fam",sep=""),head=F)
updated_fam<-fread(paste(path,"pheno/gwas-mega-2-core_updated_sex_gender_phenotype_Dec22.txt.gz",sep=""),head=T)
famall<-merge(famall[,c(1:4)],updated_fam[,c(1,9,8)],by="V2",all.x=T,sort=F)
rm(updated_fam)

dim(famall)
# [1] 15874     6
dim(dat)
# [1] 15874     8
dim(dat[which(dat$V1 %in% famall$V1),])
# [1] 15874     8

# add updated phenotype and sex:


colnames(dat)[6]<-"pheno"
aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary)
#                   Group.1  x.1  x.2
# 1 belgium_franchimont_gsa  591  893
# 2   belgium_inf1_old_gwas  893  513
# 3   belgium_inf2_old_gwas  110  161
# 4       belgium_louis_gsa  596  909
# 5    belgium_vermeire_gsa  811 3154
# 6         franchimont_gsa 1368 1230
# 7            vermeire_gsa  807 3838

colnames(dat)[5]<-"sex"
aggregate(as.factor(dat$sex), by=list(dat$cohort), FUN=summary)

table(dat$V1==dat$V2)
# TRUE 
# 15874

kin<-read.table(paste(path,"pre_imputation/QC/relatedness/belgium_merged_step25.kin0",sep=""),head=T)
table(kin$InfType)
# Dup/MZ     FS     PO 
#   5925    222    473  

table(cut(kin$Kinship,breaks=c(1,0.354,0.177)))
# (0.177,0.354]     (0.354,1] 
#           695          5925

dup<-kin[which(kin$Kinship>0.354 & kin$InfType=="Dup/MZ"),c("FID1","FID2","Kinship")]
dim(dup)
# [1] 5925   3

dup<-merge(dup,dat[,c("V1","cohort","sex","pheno","F_MISS")],by.x="FID2",by.y="V1",all.x=T,sort=F)
colnames(dup)[(ncol(dup)-3):ncol(dup)]<-paste(colnames(dup)[(ncol(dup)-3):ncol(dup)],"FID2",sep="_")

dup<-merge(dup,dat[,c("V1","cohort","sex","pheno","F_MISS")],by.x="FID1",by.y="V1",all.x=T,sort=F)
colnames(dup)[(ncol(dup)-3):ncol(dup)]<-paste(colnames(dup)[(ncol(dup)-3):ncol(dup)],"FID1",sep="_")

x<-table(dup$cohort_FID1,dup$cohort_FID2)
x
# belgium_franchimont_gsa belgium_inf1_old_gwas
# belgium_franchimont_gsa                       0                     3
# belgium_inf1_old_gwas                        10                     0
# belgium_inf2_old_gwas                        77                     0
# belgium_louis_gsa                             6                   336
# belgium_vermeire_gsa                          0                     3
# franchimont_gsa                               0                     4
# vermeire_gsa                                  0                     3
# 
# belgium_inf2_old_gwas belgium_vermeire_gsa
# belgium_franchimont_gsa                     0                    7
# belgium_inf1_old_gwas                       0                    0
# belgium_inf2_old_gwas                       0                    1
# belgium_louis_gsa                          34                    2
# belgium_vermeire_gsa                        0                    0
# franchimont_gsa                             0                    0
# vermeire_gsa                                0                    0
# 
# franchimont_gsa vermeire_gsa
# belgium_franchimont_gsa            1480            8
# belgium_inf1_old_gwas                10            0
# belgium_inf2_old_gwas                79            1
# belgium_louis_gsa                     7            2
# belgium_vermeire_gsa                  7         3836
# franchimont_gsa                       0            9
# vermeire_gsa                          0            0


xx<-dup[,c("cohort_FID1","cohort_FID2")]
yy<-dup[,c("cohort_FID2","cohort_FID1")]
colnames(yy)<-colnames(xx)

dim(xx)
# [1] 5925   2

xx<-rbind(xx,yy)
x<-table(xx$cohort_FID1,xx$cohort_FID2)
x
#                         belgium_franchimont_gsa belgium_inf1_old_gwas
# belgium_franchimont_gsa                       0                    13
# belgium_inf1_old_gwas                        13                     0
# belgium_inf2_old_gwas                        77                     0
# belgium_louis_gsa                             6                   336
# belgium_vermeire_gsa                          7                     3
# franchimont_gsa                            1480                    14
# vermeire_gsa                                  8                     3
# 
#                         belgium_inf2_old_gwas belgium_louis_gsa
# belgium_franchimont_gsa                    77                 6
# belgium_inf1_old_gwas                       0               336
# belgium_inf2_old_gwas                       0                34
# belgium_louis_gsa                          34                 0
# belgium_vermeire_gsa                        1                 2
# franchimont_gsa                            79                 7
# vermeire_gsa                                1                 2
# 
#                         belgium_vermeire_gsa franchimont_gsa vermeire_gsa
# belgium_franchimont_gsa                    7            1480            8
# belgium_inf1_old_gwas                      3              14            3
# belgium_inf2_old_gwas                      1              79            1
# belgium_louis_gsa                          2               7            2
# belgium_vermeire_gsa                       0               7         3836
# franchimont_gsa                            7               0            9
# vermeire_gsa                            3836               9            0

write.table(x,paste(path,"pre_imputation/QC/relatedness/belgium_dup_summary_table_step25.tsv",sep=""),col.names=T,row.names=T,sep="\t",quote=F)

table(dup$sex_FID1,dup$sex_FID2)
#      0    1    2
# 0   11    0    1
# 1    1 2728    0
# 2    0    0 3184

dup<-dup[,c("FID1","FID2","Kinship","cohort_FID1","pheno_FID1","F_MISS_FID1","cohort_FID2","pheno_FID2","F_MISS_FID2")]
table(dup$pheno_FID1,dup$pheno_FID2)
#      1    2
# 1 1481    0
# 2    0 4444


colnames(dat)[6]<-"pheno"
aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary)
#                   Group.1  x.1  x.2
# 1 belgium_franchimont_gsa  591  893
# 2   belgium_inf1_old_gwas  893  513
# 3   belgium_inf2_old_gwas  110  161
# 4       belgium_louis_gsa  596  909
# 5    belgium_vermeire_gsa  811 3154
# 6         franchimont_gsa 1368 1230
# 7            vermeire_gsa  807 3838

###############################

dup_ids<-c(as.character(dup$FID1),as.character(dup$FID2))
length(dup_ids)
table(length(dup_ids)==(nrow(dup)*2))
dup_ids<-dup_ids[!duplicated(dup_ids)]
length(dup_ids)
dim(dup)
# [1] 5925   9


# CHECK ALSO THAT BOTH SAMPLES HAVE SAME PHENOTYPE AND GENDER, 
# - keep all new gwas samples

# - for ukbb-ukbb duplicated pairs:
#    - keep sample with pheno data when both samples have the same phenotype, otherwise exclude both

old_gwas<-c("belgium_inf1_old_gwas","belgium_inf2_old_gwas")
new_gwas<-c("belgium_franchimont_gsa","belgium_louis_gsa","belgium_vermeire_gsa","franchimont_gsa","vermeire_gsa")

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
    
    # make an exception to deal with sex = 0; allowing that to happen
    if ( (dim(table(data$sex))!=1 | dim(table(data$pheno))!=1) & all(data$sex %in% c(1,2))) {
      
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
      
      # OPTION B: samples have same pheno
      
    } else {
      
      # B.1: if one sample from franchimont_gsa, keep franchimont_gsa (update of belgium_franchimont) or vermeire
      
      if (any(data$cohort %in% c("franchimont_gsa","vermeire_gsa"))) {
        
        keep_sample<-data$V1[which(data$cohort %in% c("franchimont_gsa","vermeire_gsa"))]
        
        # if more than one new pair in new gwas, keep sample with largest call rate:
        if (length(keep_sample)>1) {
          keep_sample<-data[which(data$cohort %in% c("franchimont_gsa","vermeire_gsa")),]
          keep_sample<-keep_sample$V1[which(keep_sample$F_MISS==min(keep_sample$F_MISS))][1]# to deal with >1 with same F_miss
        }
        
        if(!exists("data_remove")) {
          data_remove<-data[which(!data$V1 %in% keep_sample),]
        } else {
          data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
        }
        
        
      # B.2: samples in new and old arrays, keep only samples in new arrays:
        
      } else if (any(data$cohort %in% old_gwas) & any(data$cohort %in% new_gwas)) {
        
        keep_sample<-data$V1[which(!data$cohort %in% old_gwas)]
        
        # if more than one new pair in new gwas, keep sample with largest call rate :
        if (length(keep_sample)>1) {
          keep_sample<-data[which(!data$cohort %in% old_gwas),]
          keep_sample<-keep_sample$V1[which(keep_sample$F_MISS==min(keep_sample$F_MISS))][1]# to deal with >1 with same F_miss
        }
        
        if(!exists("data_remove")) {
          data_remove<-data[which(!data$V1 %in% keep_sample),]
        } else {
          data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
        }
        
      # B.3: samples in new OR old arrays, keep only samples with best call rate:
        
      } else {
        
        keep_sample<-data$V1[which(data$F_MISS==min(data$F_MISS))][1]# to deal with >1 with same F_miss
        
        if(!exists("data_remove")) {
          data_remove<-data[which(!data$V1 %in% keep_sample),]
        } else {
          data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
        }
        
      }
    }
  } else {
    
    # not all combinations of duplicated pairs found, show list of IDs, to manually inspect issues
    print(i)
    print(paste("Number of expected combinations: ",n_possible_combinations,sep=""))
    print(paste("Number of observed combinations: ",nrow(dup[which((dup$FID1 %in% ids_tmp) | (dup$FID2 %in% ids_tmp)),]),sep=""))
    print(ids_tmp)
    
  }
  
}

dim(data_remove)
# [1] 11850    5

data_remove<-data_remove[which(!duplicated(data_remove$V1)),]
dim(data_remove)
# [1] 5799   5

table(data_remove$cohort)
# belgium_franchimont_gsa   belgium_inf1_old_gwas   belgium_inf2_old_gwas 
# 1480                     350                     113 
# belgium_louis_gsa    belgium_vermeire_gsa         franchimont_gsa 
# 10                    3837                       5 
# vermeire_gsa 
# 4 

aggregate(as.factor(dat$pheno[which(!dat$V1 %in% data_remove$V1)]), by=list(dat$cohort[which(!dat$V1 %in% data_remove$V1)]), FUN=summary)
#                   Group.1  x.1  x.2
# 1 belgium_franchimont_gsa    1    3
# 2   belgium_inf1_old_gwas  833  223
# 3   belgium_inf2_old_gwas   76   82
# 4       belgium_louis_gsa  596  899
# 5    belgium_vermeire_gsa   14  114
# 6         franchimont_gsa 1368 1225
# 7            vermeire_gsa  807 3834

head(dup[which(!dup$FID1 %in% data_remove$V1 & !dup$FID2 %in% data_remove$V1),])
# 0 OK

data_inconsist<-data_inconsist[which(!duplicated(data_inconsist$V1)),]
dim(data_inconsist)
# [1] 0  6
# data_inconsist


data_remove_tmp<-data_remove[,c(1,1)]
colnames(data_remove_tmp)<-c("FID","IID")

write.table(data_remove,paste(path,"pre_imputation/QC/relatedness/belgium_data_remove_step25",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
write.table(data_remove_tmp,paste(path,"pre_imputation/QC/relatedness/list_belgium_data_remove_step25",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
write.table(data_inconsist,paste(path,"pre_imputation/QC/relatedness/belgium_data_inconsist_step25",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
write.table(aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary),paste(path,"pre_imputation/QC/relatedness/belgium_case_control_pre_per_study_relatedness_step25",sep=""),col.names=T,row.names=F,sep="\t",quote=F)

####################################################################################################################################################

