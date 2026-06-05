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

path_variants=${path_gwas}pre_imputation/QC/relatedness/list_common_variants_among_cohorts_step25_ukibdgc

for i in ${studies[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample \
--allow-no-sex \
--extract ${path_variants} \
--make-bed \
--out ${path_gwas}pre_imputation/QC/relatedness/${i}_subset_uk_step25
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


wc -l ${path_gwas}pre_imputation/QC/relatedness/*_subset_uk_step25.fam
# 22498 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/all_hce_subset_uk_step25.fam
# 4658 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/gwas1_subset_uk_step25.fam
# 7772 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/gwas2_subset_uk_step25.fam
# 34928 total

wc -l ${path_gwas}pre_imputation/QC/relatedness/*_subset_uk_step25.bim
# 35897 # OK


##############################################################################################################################################
  
#######################
# 3.- MERGE ALL FILES #
#######################
  
# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=8000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q cpu-interactive R

  
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
  dat[i,1]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_uk_step25.bed",sep="")
  dat[i,2]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_uk_step25.bim",sep="")
  dat[i,3]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_uk_step25.fam",sep="")
}

write.table(dat,paste(path,"pre_imputation/QC/relatedness/list_cohorts_tomerge_uk_step25.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#############################

# hpc-server
  
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/relatedness/gwas1_subset_uk_step25 \
--allow-no-sex \
--merge-list ${path_gwas}pre_imputation/QC/relatedness/list_cohorts_tomerge_uk_step25.txt \
--make-bed --out ${path_gwas}pre_imputation/QC/relatedness/ukibdgc_merged_step25
  
# 35897 variants and 34928 people pass filters and QC.
# Among remaining phenotypes, 16168 are cases and 18760 are controls.


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
  
  
# (35897*34928/4)/1E+6
# [1] 313.4526
  
path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=1000
  
bsub -J"kg" -M"$MEM" -n8 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/relatedness/logs/stderr_king_ukibdgc_step25 \
-o ${path_gwas}pre_imputation/QC/relatedness/logs/stdout_king_ukibdgc_step25 \
"/path/to/software/./king \
-b ${path_gwas}pre_imputation/QC/relatedness/ukibdgc_merged_step25.bed \
--related --cpus 8 --prefix ${path_gwas}pre_imputation/QC/relatedness/ukibdgc_merged_step25"


# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=8000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q cpu-interactive R


library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"
  
cohorts<-c("gwas1","gwas2","all_hce")
  
  
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
# all_hce   gwas1   gwas2 
#  22498    4658    7772 
  
# double check all samples included in exercise
famall<-read.table(paste(path,"pre_imputation/QC/relatedness/ukibdgc_merged_step25.fam",sep=""),head=F)
updated_fam<-fread(paste(path,"pheno/gwas-mega-2-core_updated_sex_gender_phenotype_Dec22.txt.gz",sep=""),head=T)
famall<-merge(famall[,c(1:4)],updated_fam[,c(1,9,8)],by="V2",all.x=T,sort=F)
rm(updated_fam)

dim(famall)
# [1] 34928     6
dim(dat)
# [1] 34928     8
dim(dat[which(dat$V1 %in% famall$V1),])
# [1] 34928     8
  
colnames(dat)[6]<-"pheno"
aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary)
#   Group.1   x.1   x.2
# 1 all_hce 10425 12073
# 2   gwas1  2922  1736
# 3   gwas2  5413  2359

colnames(dat)[5]<-"sex"
aggregate(as.factor(dat$sex), by=list(dat$cohort), FUN=summary)
#   Group.1   x.0   x.1   x.2
# 1 all_hce    27 10402 12069
# 2   gwas1     2  2111  2545
# 3   gwas2    78  3831  3863

table(dat$V1==dat$V2)
# TRUE 
# 35070

kin<-read.table(paste(path,"pre_imputation/QC/relatedness/ukibdgc_merged_step25.kin0",sep=""),head=T)
table(kin$InfType)
# Dup/MZ     FS     PO     UN 
# 2830    206    935      1 

table(cut(kin$Kinship,breaks=c(1,0.354,0.177)))
# (0.177,0.354]     (0.354,1] 
#          1141          2831

summary(kin$Kinship[which(kin$Kinship>0.354 & kin$InfType!="Dup/MZ")])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.414   0.414   0.414   0.414   0.414   0.414 


rel<-kin[which(kin$Kinship<=0.354 & kin$Kinship>0.177),c("FID1","FID2","Kinship")]

rel_ids<-c(as.character(rel$FID1),as.character(rel$FID2))
length(rel_ids)
# [1] 2282
table(length(rel_ids)==(nrow(rel)*2))
# TRUE 
# 1 


rel<-merge(rel,dat[,c("V1","cohort","sex","pheno","F_MISS")],by.x="FID2",by.y="V1",all.x=T,sort=F)
colnames(rel)[(ncol(rel)-3):ncol(rel)]<-paste(colnames(rel)[(ncol(rel)-3):ncol(rel)],"FID2",sep="_")
  
rel<-merge(rel,dat[,c("V1","cohort","sex","pheno","F_MISS")],by.x="FID1",by.y="V1",all.x=T,sort=F)
colnames(rel)[(ncol(rel)-3):ncol(rel)]<-paste(colnames(rel)[(ncol(rel)-3):ncol(rel)],"FID1",sep="_")
  
x<-table(rel$cohort_FID1,rel$cohort_FID2)
x
#           all_hce gwas1 gwas2
#   all_hce    1057    26    27
#   gwas1         0     0     9
#   gwas2         2    19     1

xx<-rel[,c("cohort_FID1","cohort_FID2")]
yy<-rel[,c("cohort_FID2","cohort_FID1")]
colnames(yy)<-colnames(xx)
  
dim(xx)
# [1] 1141    2


xx<-rbind(xx,yy)
x<-table(xx$cohort_FID1,xx$cohort_FID2)
x
#           all_hce gwas1 gwas2
#   all_hce    2114    26    29
#   gwas1        26     0    28
#   gwas2        29    28     2



rm(data_remove)

for (i in 1:length(rel_ids)) {
  
  tmp1<-rel[which(rel$FID1==rel_ids[i]),]
  tmp2<-rel[which(rel$FID2==rel_ids[i]),]
  colnames(tmp2)[1:2]<-colnames(tmp2)[2:1]
  
  tmp<-rbind(tmp1,tmp2)
  
  ids_tmp<-c(as.character(tmp$FID1),as.character(tmp$FID2))
  ids_tmp<-ids_tmp[!duplicated(ids_tmp)]
  
  #keep same order as in rel_ids and in all
  ids_tmp<-ids_tmp[match(rel_ids[which(rel_ids %in% ids_tmp)],ids_tmp)]
  
  # remove the ones with smaller call rate
    
    data<-as.data.frame(matrix(ncol=1,nrow=length(ids_tmp)))
    data$V1<-ids_tmp
    data<-merge(data,dat,by="V1",all.x=T,sort=F)

    if (dim(table(data$pheno))==2) {
      keep_sample<-data$V1[which(data$pheno==2)][1]
    } else {
      keep_sample<-data$V1[which(data$F_MISS==min(data$F_MISS))][1]
    }
    
    

    if(!exists("data_remove")) {
      data_remove<-data[which(!data$V1 %in% keep_sample),]
    } else {
      data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
    }
  
}

data_remove<-data_remove[!duplicated(data_remove),]
dim(data_remove)
# [1] 1030    8

nrow(rel[which(rel$FID1 %in% data_remove$V1 | rel$FID2 %in% data_remove$V1),])
# [1] 1141

nrow(rel)
# [1] 1141

fwrite(data_remove,paste(path,"pre_imputation/QC/relatedness/ukiibdgc_list_first_degree_samples_to_remove.tsv.gz",sep=""),col.names=T,row.names=F,sep="\t",quote=F)



dup<-kin[which(kin$Kinship>0.354 & kin$InfType=="Dup/MZ"),c("FID1","FID2","Kinship")]
dim(dup)
# [1] 2830    3


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
# all_hce       0   123    65
# gwas1       123     0  2642
# gwas2        65  2642     0

write.table(x,paste(path,"pre_imputation/QC/relatedness/ukibdgc_dup_summary_table_step25.tsv",sep=""),col.names=T,row.names=T,sep="\t",quote=F)

# duplicated pairs in GWAS1 and GWAS2:
table(dup[which(dup$cohort_FID1 %in% c("gwas1","gwas2") & dup$cohort_FID2 %in% c("gwas1","gwas2")),"pheno_FID1"],
      dup[which(dup$cohort_FID1 %in% c("gwas1","gwas2") & dup$cohort_FID2 %in% c("gwas1","gwas2")),"pheno_FID2"])
#      1    2
# 1 2637    0
# 2    0    5

# duplicated pairs in GWAS1 and GWAS3:
table(dup[which(dup$cohort_FID1 %in% c("gwas1","all_hce") & dup$cohort_FID2 %in% c("gwas1","all_hce")),"pheno_FID1"],
      dup[which(dup$cohort_FID1 %in% c("gwas1","all_hce") & dup$cohort_FID2 %in% c("gwas1","all_hce")),"pheno_FID2"])
#     1   2
# 1   2   0
# 2   2 119

# duplicated pairs in GWAS2 and GWAS3:
table(dup[which(dup$cohort_FID1 %in% c("gwas2","all_hce") & dup$cohort_FID2 %in% c("gwas2","all_hce")),"pheno_FID1"],
      dup[which(dup$cohort_FID1 %in% c("gwas2","all_hce") & dup$cohort_FID2 %in% c("gwas2","all_hce")),"pheno_FID2"])
#    1  2
# 1  2  0
# 2  2 61



table(dup$sex_FID1,dup$sex_FID2)
#      0    1    2
# 1   27 1358    0
# 2    0    0 1445

(dup[which(dup$sex_FID2=="0" & dup$sex_FID1==1),])
# all gwas2 sex =0, gwas1 sex = 1


dup<-dup[,c("FID1","FID2","Kinship","cohort_FID1","pheno_FID1","F_MISS_FID1","cohort_FID2","pheno_FID2","F_MISS_FID2")]
table(dup$pheno_FID1,dup$pheno_FID2)
#      1    2
# 1 2641    0
# 2    4  185
  
dup_ids<-c(as.character(dup$FID1),as.character(dup$FID2))
length(dup_ids)
table(length(dup_ids)==(nrow(dup)*2))
dup_ids<-dup_ids[!duplicated(dup_ids)]
length(dup_ids)
dim(dup)
# [1] 2830    9


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
    # make an exception to deal with sex = 0; allowing that to happen
    if ( (dim(table(data$sex))!=1 | dim(table(data$pheno))!=1) & all(data$sex %in% c(1,2))) {
      
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
# [1] 5667    5

data_remove<-data_remove[which(!duplicated(data_remove$V1)),]
dim(data_remove)
# [1] 2831    5

table(data_remove$cohort)
# all_hce   gwas1   gwas2 
#     186       3    2642 

head(dup[which(!dup$FID1 %in% data_remove$V1 & !dup$FID2 %in% data_remove$V1),])
# 0 OK

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

write.table(data_remove,paste(path,"pre_imputation/QC/relatedness/ukibdgc_data_remove_step25",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
write.table(data_remove_tmp,paste(path,"pre_imputation/QC/relatedness/list_ukiibdgc_data_remove_step25",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
write.table(data_inconsist,paste(path,"pre_imputation/QC/relatedness/ukibdgc_data_inconsist_step25",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
write.table(aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary),paste(path,"pre_imputation/QC/relatedness/ukibdgc_case_control_pre_per_study_relatedness_step25",sep=""),col.names=T,row.names=F,sep="\t",quote=F)

q("no")


cp /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/ukibdgc_data_remove_step25 /path/to/project
cp /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/ukiibdgc_list_first_degree_samples_to_remove.tsv.gz /path/to/project