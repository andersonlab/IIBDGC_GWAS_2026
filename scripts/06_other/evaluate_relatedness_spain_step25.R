# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# EVALUATE RELATEDNESS BETWEEN SPAIN SAMPLES

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
  
studies=(spain_gsa basque_gsa)
  
path_variants=${path_gwas}pre_imputation/QC/relatedness/list_common_variants_among_cohorts_step25_spain
  
for i in ${studies[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample \
--allow-no-sex \
--extract ${path_variants} \
--make-bed \
--out ${path_gwas}pre_imputation/QC/relatedness/${i}_subset_spain_step25
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


wc -l ${path_gwas}pre_imputation/QC/relatedness/*_subset_spain_step25.fam
# 1481 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/basque_gsa_subset_spain_step25.fam
# 3408 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/spain_gsa_subset_spain_step25.fam
# 4889 total

wc -l ${path_gwas}pre_imputation/QC/relatedness/*_subset_spain_step25.bim
# 106519  OK

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

cohorts<-c("spain_gsa","basque_gsa")

cohorts<-cohorts[-1]

dat<-matrix(ncol=3,nrow=length(cohorts))
dat<-as.data.frame(dat)

for (i in 1:nrow(dat)){
  dat[i,1]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_spain_step25.bed",sep="")
  dat[i,2]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_spain_step25.bim",sep="")
  dat[i,3]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_spain_step25.fam",sep="")
}

write.table(dat,paste(path,"pre_imputation/QC/relatedness/list_cohorts_tomerge_spain_step25.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#############################

# hpc-server

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/relatedness/spain_gsa_subset_spain_step25 \
--allow-no-sex \
--merge-list ${path_gwas}pre_imputation/QC/relatedness/list_cohorts_tomerge_spain_step25.txt \
--make-bed --out ${path_gwas}pre_imputation/QC/relatedness/spain_merged_step25

# 106519 variants and 4889 people pass filters and QC.
# Among remaining phenotypes, 2463 are cases and 2426 are controls.



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


# (106519*4889/4)/1E+6
# [1] 130.1928

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=1000

bsub -J"kg" -M"$MEM" -n8 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/relatedness/logs/stderr_king_spain_step25 \
-o ${path_gwas}pre_imputation/QC/relatedness/logs/stdout_king_spain_step25 \
"/path/to/software/./king \
-b ${path_gwas}pre_imputation/QC/relatedness/spain_merged_step25.bed \
--related --cpus 8 --prefix ${path_gwas}pre_imputation/QC/relatedness/spain_merged_step25"
# Job <66743> is submitted to queue <normal>.

less ${path_gwas}pre_imputation/QC/relatedness/log/stdout_king_spain_step25
# Relationship summary (total relatives: 0 by pedigree, 25 by inference)
#                MZ      PO      FS      2nd
# =====================================================
#   Inference     0       9       15      0



#### R

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("spain_gsa","basque_gsa")

for (i in 1:length(cohorts)) {
  tmp<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],
                        "_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.fam",sep=""),head=F)
  sample_miss<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],
                                "_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample.imiss",sep=""),head=T)
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
# basque_gsa  spain_gsa 
# 1481       3408 

# double check all samples included in exercise
famall<-read.table(paste(path,"pre_imputation/QC/relatedness/spain_merged_step25.fam",sep=""),head=F)
updated_fam<-fread(paste(path,"pheno/gwas-mega-2-core_updated_sex_gender_phenotype_Dec22.txt.gz",sep=""),head=T)
famall<-merge(famall[,c(1:4)],updated_fam[,c(1,9,8)],by="V2",all.x=T,sort=F)
rm(updated_fam)


dim(famall)
# [1] 4889     6
dim(dat)
# [1] 4889     8
dim(dat[which(dat$V1 %in% famall$V1),])
# [1] 4889     8

colnames(dat)[6]<-"pheno"
aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary)
#      Group.1  x.1  x.2
# 1 basque_gsa  960  521
# 2  spain_gsa 1466 1942

colnames(dat)[5]<-"sex"
aggregate(as.factor(dat$sex), by=list(dat$cohort), FUN=summary)
#      Group.1  x.0  x.1  x.2
# 1 basque_gsa    6  905  570
# 2  spain_gsa   56 1861 1491

table(dat$V1==dat$V2)
# TRUE 
# 4889

kin<-read.table(paste(path,"pre_imputation/QC/relatedness/spain_merged_step25.kin0",sep=""),head=T)
table(kin$InfType)
# FS PO 
# 15  9

table(cut(kin$Kinship,breaks=c(1,0.354,0.177)))
# (0.177,0.354]     (0.354,1] 
#  24             0 

write.table(aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary),paste(path,"pre_imputation/QC/relatedness/spain_case_control_pre_per_study_relatedness_step25",sep=""),col.names=T,row.names=F,sep="\t",quote=F)


##############################
# remove intermediate files:

ls -la ${path_gwas}pre_imputation/QC/relatedness/*_subset_spain_step25*
rm ${path_gwas}pre_imputation/QC/relatedness/*_subset_spain_step25*
