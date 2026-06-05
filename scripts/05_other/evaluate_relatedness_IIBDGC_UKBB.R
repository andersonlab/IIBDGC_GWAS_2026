# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# EVALUATE RELATEDNESS BETWEEN IIBDGC SAMPLES

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
  
studies=(australia_omniexome gwas1 gwas2 all_hce pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa chop_old_gwas german_illu550_old_gwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas niddk_cd_old_gwas)
  
for i in ${studies[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample \
--extract ${path}pre_imputation/QC/relatedness/list_common_variants_among_cohorts_iibdgc \
--make-bed \
--out ${path}pre_imputation/QC/relatedness/${i}_subset
done

# estimate missingness for later step
for i in ${studies[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample \
--missing \
--out ${path}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample
done


wc -l ${path}pre_imputation/QC/relatedness/*_subset.fam
# 22528 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/all_hce_subset.fam
# 1287 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/australia_omniexome_subset.fam
# 1406 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_inf1_old_gwas_subset.fam
# 271 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_inf2_old_gwas_subset.fam
# 8505 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/chop_old_gwas_subset.fam
# 2782 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/german_affy6_old_gwas_subset.fam
# 1611 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/german_illu550_old_gwas_subset.fam
# 4677 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/gwas1_subset.fam
# 7699 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/gwas2_subset.fam
# 979 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/italy_gsa_subset.fam
# 13715 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/kiel_austria_sibdcs_gsa_subset.fam
# 4548 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/netherlands_gsa_subset.fam
# 5383 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_broad_gsa_subset.fam
# 1722 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_cd_old_gwas_subset.fam
# 7961 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_feinstein_gsa_subset.fam
# 549 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/norway_affy6_old_gwas_subset.fam
# 2739 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/pittsburgh_gsa_subset.fam
# 257 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/slovenia_gsa_subset.fam
# 3387 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/spain_gsa_subset.fam
# 1334 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/sweden_gsa_subset.fam
# 93340 total


wc -l ${path}pre_imputation/QC/relatedness/*_subset.bim
# 4608  OK

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

cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa"
           ,"chop_old_gwas","german_illu550_old_gwas","german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_cd_old_gwas") # illu 317

# needs to be QCed:
# "swedish_uc_old_gwas"
# "niddk_uc_old_gwas"

# edit german_illu550_old_gwas fam files so IDs are not duplicated (all of them in affy6) and confirm all samples are in affy6, and merge again:

fam<-read.table(paste(path,"pre_imputation/QC/relatedness/german_illu550_old_gwas_subset.fam",sep=""),head=F)
fam$V1<-paste(fam$V1,"_illu550",sep="")
fam$V2<-paste(fam$V2,"_illu550",sep="")
write.table(fam,paste(path,"pre_imputation/QC/relatedness/german_illu550_old_gwas_subset.fam",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

#####

cohorts<-cohorts[-1]

dat<-matrix(ncol=3,nrow=length(cohorts))
dat<-as.data.frame(dat)

for (i in 1:nrow(dat)){
  dat[i,1]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset.bed",sep="")
  dat[i,2]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset.bim",sep="")
  dat[i,3]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset.fam",sep="")
}

write.table(dat,paste(path,"pre_imputation/QC/relatedness/list_cohorts_tomerge.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#############################

# hpc-server

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/relatedness/australia_omniexome_subset \
--merge-list ${path}pre_imputation/QC/relatedness/list_cohorts_tomerge.txt \
--make-bed --out ${path}pre_imputation/QC/relatedness/iibdgc_merged

# 4608 variants and 93340 people pass filters and QC.
# Among remaining phenotypes, 50016 are cases and 43324 are controls.


# remove intermediate files:
ls -la ${path}pre_imputation/QC/relatedness/*_subset*
rm ${path}pre_imputation/QC/relatedness/*_subset*
  
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


# (93340*4608/4)/1E+6
# [1] 107.5277

path=/path/to/ibdgwas/IIBDGC/
MEM=200

bsub -J"kg" -M"$MEM" -n8 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path}pre_imputation/QC/relatedness/log/stderr_king_iibdgc \
-o ${path}pre_imputation/QC/relatedness/log/stdout_king_iibdgc \
"/path/to/software/./king \
-b ${path}pre_imputation/QC/relatedness/iibdgc_merged.bed \
--related --cpus 8 --prefix ${path}pre_imputation/QC/relatedness/iibdgc_merged"
# Job <428640> is submitted to queue <normal>.



################################

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa"
           ,"chop_old_gwas","german_illu550_old_gwas","german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_cd_old_gwas")


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
# all_hce     australia_omniexome   belgium_inf1_old_gwas 
# 22528                    1287                    1406 
# belgium_inf2_old_gwas           chop_old_gwas   german_affy6_old_gwas 
# 271                    8505                    2782 
# german_illu550_old_gwas                   gwas1                   gwas2 
# 1611                    4677                    7699 
# italy_gsa kiel_austria_sibdcs_gsa         netherlands_gsa 
# 979                   13715                    4548 
# niddk_broad_gsa       niddk_cd_old_gwas     niddk_feinstein_gsa 
# 5383                    1722                    7961 
# norway_affy6_old_gwas          pittsburgh_gsa            slovenia_gsa 
# 549                    2739                     257 
# spain_gsa              sweden_gsa 
# 3387                    1334 

# double check all samples included in exercise
famall<-read.table(paste(path,"pre_imputation/QC/relatedness/iibdgc_merged.fam",sep=""),head=F)

dim(famall)
# [1] 93340     6

dim(dat)
# [1] 93340     8

dim(dat[which(dat$V1 %in% famall$V1),])
# [1] 93340     8

colnames(dat)[6]<-"pheno"

aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary)
#                    Group.1   x.1   x.2
# 1                  all_hce 10430 12098
# 2      australia_omniexome   612   675
# 3    belgium_inf1_old_gwas   889   517
# 4    belgium_inf2_old_gwas   111   160
# 5            chop_old_gwas  6097  2408
# 6    german_affy6_old_gwas  1749  1033
# 7  german_illu550_old_gwas  1131   480
# 8                    gwas1  2930  1747
# 9                    gwas2  5360  2339
# 10               italy_gsa   385   594
# 11 kiel_austria_sibdcs_gsa  4454  9261
# 12         netherlands_gsa   547  4001
# 13         niddk_broad_gsa  1011  4372
# 14       niddk_cd_old_gwas   931   791
# 15     niddk_feinstein_gsa  2351  5610
# 16   norway_affy6_old_gwas   281   268
# 17          pittsburgh_gsa  1485  1254
# 18            slovenia_gsa   171    86
# 19               spain_gsa  1468  1919
# 20              sweden_gsa   931   403

colnames(dat)[5]<-"sex"

aggregate(as.factor(dat$sex), by=list(dat$cohort), FUN=summary)
#                    Group.1   x.1   x.2
# 1                  all_hce 10426 12102
# 2      australia_omniexome   594   693
# 3    belgium_inf1_old_gwas   991   415
# 4    belgium_inf2_old_gwas   137   134
# 5            chop_old_gwas  4479  4026
# 6    german_affy6_old_gwas  1358  1424
# 7  german_illu550_old_gwas   636   975
# 8                    gwas1  2122  2555
# 9                    gwas2  3831  3868
# 10               italy_gsa   528   451
# 11 kiel_austria_sibdcs_gsa  6806  6909
# 12         netherlands_gsa  1895  2653
# 13         niddk_broad_gsa  2805  2578
# 14       niddk_cd_old_gwas   838   884
# 15     niddk_feinstein_gsa  4013  3948
# 16   norway_affy6_old_gwas   304   245
# 17          pittsburgh_gsa  1195  1544
# 18            slovenia_gsa    55   202
# 19               spain_gsa  1883  1504
# 20              sweden_gsa   773   561

table(dat$V1==dat$V2)
# TRUE 
# 93340


kin<-read.table(paste(path,"pre_imputation/QC/relatedness/iibdgc_merged.kin0",sep=""),head=T)
# No informative IBD segments.
# Inference will be based on kinship estimation only.

table(cut(kin$Kinship,breaks=c(1,0.354,0.177)))
# (0.177,0.354]     (0.354,1] 
#          2798          6972 

dup<-kin[which(kin$Kinship>0.354),c("FID1","FID2","Kinship")]

dup<-merge(dup,dat[,c("V1","cohort","sex","pheno","F_MISS")],by.x="FID2",by.y="V1",all.x=T,sort=F)
colnames(dup)[(ncol(dup)-3):ncol(dup)]<-paste(colnames(dup)[(ncol(dup)-3):ncol(dup)],"FID2",sep="_")

dup<-merge(dup,dat[,c("V1","cohort","sex","pheno","F_MISS")],by.x="FID1",by.y="V1",all.x=T,sort=F)
colnames(dup)[(ncol(dup)-3):ncol(dup)]<-paste(colnames(dup)[(ncol(dup)-3):ncol(dup)],"FID1",sep="_")

x<-table(dup$cohort_FID1,dup$cohort_FID2)
x
#                         all_hce chop_old_gwas german_affy6_old_gwas
# all_hce                       2           104                     0
# chop_old_gwas                 2             0                     0
# german_affy6_old_gwas         0             0                     0
# german_illu550_old_gwas       0             0                     5
# gwas1                         0             0                     0
# gwas2                        10             0                     0
# kiel_austria_sibdcs_gsa       0             0                    72
# netherlands_gsa               0             0                     0
# niddk_broad_gsa               0             0                     0
# niddk_cd_old_gwas             1            59                     0
# niddk_feinstein_gsa           0             0                     0
# pittsburgh_gsa                0             0                     0
# 
#                         german_illu550_old_gwas gwas1 gwas2 italy_gsa
# all_hce                                       0   123    56         0
# chop_old_gwas                                 0    27    15        11
# german_affy6_old_gwas                      1101     0     0         0
# german_illu550_old_gwas                       0     0     0         0
# gwas1                                         0     0  2508         0
# gwas2                                         0   110     0         0
# kiel_austria_sibdcs_gsa                     391     0     0         0
# netherlands_gsa                               0     0     0         0
# niddk_broad_gsa                               0     0     0         0
# niddk_cd_old_gwas                             0     0     0         0
# niddk_feinstein_gsa                           0     0     0         0
# pittsburgh_gsa                                0     0     0         0
# 
#                         kiel_austria_sibdcs_gsa niddk_broad_gsa
# all_hce                                      16               0
# chop_old_gwas                                 1             217
# german_affy6_old_gwas                        82               0
# german_illu550_old_gwas                      70               0
# gwas1                                         0               0
# gwas2                                         0               0
# kiel_austria_sibdcs_gsa                       0               0
# netherlands_gsa                               1               0
# niddk_broad_gsa                               1               0
# niddk_cd_old_gwas                             0             411
# niddk_feinstein_gsa                           0               0
# pittsburgh_gsa                                0               0
# 
#                         niddk_cd_old_gwas niddk_feinstein_gsa pittsburgh_gsa
# all_hce                                 0                   0              0
# chop_old_gwas                           1                 112              0
# german_affy6_old_gwas                   0                   0              0
# german_illu550_old_gwas                 0                   0              0
# gwas1                                   0                   0              0
# gwas2                                   0                   0              0
# kiel_austria_sibdcs_gsa                 0                   0              1
# netherlands_gsa                         0                   0              0
# niddk_broad_gsa                        41                 283            507
# niddk_cd_old_gwas                       0                 117            254
# niddk_feinstein_gsa                     4                   0            234
# pittsburgh_gsa                         22                   0              0


xx<-dup[,c("cohort_FID1","cohort_FID2")]
yy<-dup[,c("cohort_FID2","cohort_FID1")]
colnames(yy)<-colnames(xx)
  
dim(xx)
# [1] 6972    2

xx<-rbind(xx,yy)

x<-table(xx$cohort_FID1,xx$cohort_FID2)

write.table(x,paste(path,"pre_imputation/QC/relatedness/iibdgc_dup_summary_table",sep=""),col.names=T,row.names=T,sep="\t",quote=F)





table(dup$sex_FID1,dup$sex_FID2)
#      1    2
# 1 3336    0
# 2    0 3336


dup<-dup[,c("FID1","FID2","Kinship","cohort_FID1","pheno_FID1","F_MISS_FID1","cohort_FID2","pheno_FID2","F_MISS_FID2")]

table(dup$pheno_FID1,dup$pheno_FID2)
#      1    2
# 1 4272   14
# 2    8 2678


dup_ids<-c(as.character(dup$FID1),as.character(dup$FID2))
length(dup_ids)
#[1] 13944
table(length(dup_ids)==(nrow(dup)*2))
# TRUE 
# 1 

dup_ids<-dup_ids[!duplicated(dup_ids)]
length(dup_ids)
#[1] 12791

dim(dup)
# [1] 6972    9

# CHECK ALSO THAT BOTH SAMPLES HAVE SAME PHENOTYPE AND GENDER, 
# - keep all new gwas samples

# - for ukbb-ukbb duplicated pairs:
#    - keep sample with pheno data when both samples have the same phenotype, otherwise exclude both

old_gwas<-c("australia_omniexome","gwas1","gwas2","all_hce"
            ,"chop_old_gwas","german_illu550_old_gwas","german_affy6_old_gwas","norway_affy6_old_gwas"
            ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_cd_old_gwas")

new_gwas<-c("pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa")

good_array<-c("chop_old_gwas","german_illu550_old_gwas","german_affy6_old_gwas","norway_affy6_old_gwas","australia_omniexome","gwas1","gwas2","all_hce")

bad_array<-c("belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_cd_old_gwas")

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
    
    # OPTION B: samples have same pheno:

    else {

      # A.1: samples in new and old arrays, keep only samples in new arrays:

      if (any(data$cohort %in% old_gwas) & any(data$cohort %in% new_gwas)) {

        keep_sample<-data$V1[which(!data$cohort %in% old_gwas)]

        # if more than one new pair in new gwas, keep sample with largest call rate:
        if (length(keep_sample)>1) {
          keep_sample<-data[which(!data$cohort %in% old_gwas),]
          keep_sample<-keep_sample$V1[which(keep_sample$F_MISS==min(keep_sample$F_MISS))][1]# to deal with >1 with same F_miss
        }

        if(!exists("data_remove")) {
          data_remove<-data[which(!data$V1 %in% keep_sample),]
        } else {
          data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
        }
      }

      
      # OPTION B: old_gwas duplicated pairs only
      
      else if (all(data$cohort %in% old_gwas)) {
      
        # B.1: if in good and in bad array, keep sample in the best array:
        
        if (any(data$cohort %in% good_array) & any(data$cohort %in% bad_array)) {

          keep_sample<-data$V1[which(data$cohort %in% good_array)]
          
          # if more than one new pair in good array, keep sample with largest call rate:
          if (length(keep_sample)>1) {
            keep_sample<-data[which(data$cohort %in% good_array),]
            keep_sample<-keep_sample$V1[which(keep_sample$F_MISS==min(keep_sample$F_MISS))][1]# to deal with >1 with same F_miss
          }
          
          if(!exists("data_remove")) {
            data_remove<-data[which(!data$V1 %in% keep_sample),]
          } else {
            data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
          }
        }
        
        # B.2: if all in good OR all in bad array, keep sample in the best array:
        
        else if (all(data$cohort %in% good_array) | all(data$cohort %in% bad_array)) {
          
          keep_sample<-data$V1[which(data$F_MISS==min(data$F_MISS))][1]# to deal with same F_miss
          
          if(!exists("data_remove")) {
            data_remove<-data[which(!data$V1 %in% keep_sample),]
          } else {
            data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
          }
          
        }
      }
      
      # OPTION B: new_gwas duplicated pairs only, keep sample with largest call rate: 
      
      else if (all(data$cohort %in% new_gwas)) {
        
        keep_sample<-data$V1[which(data$F_MISS==min(data$F_MISS))]
        
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

# # check these:
# [1] 738
# [1] "Number of expected combinations: 1"
# [1] "Number of observed combinations: 2"
# [1] "305235_C04_usgwas5502691" "367426_C12_IBD_N5786918" 
# [1] 739
# [1] "Number of expected combinations: 1"
# [1] "Number of observed combinations: 2"
# [1] "305236_A05_usgwas5502808" "367426_C12_IBD_N5786918" 
# [1] 7132
# [1] "Number of expected combinations: 3"
# [1] "Number of observed combinations: 2"
# [1] "305235_C04_usgwas5502691" "305236_A05_usgwas5502808"
# [3] "367426_C12_IBD_N5786918" 


########################################################################

# [1] 741
# [1] "Number of expected combinations: 1"
# [1] "Number of observed combinations: 2"
# [1] "305235_C04_usgwas5502691" "367426_C12_IBD_N5786918" 
n_possible_combinations
# [1] 1
length(ids_tmp)
# [1] 2
dup[which((dup$FID1 %in% ids_tmp) | (dup$FID2 %in% ids_tmp)),]
#                         FID1                    FID2 Kinship cohort_FID1
# 941 305235_C04_usgwas5502691 367426_C12_IBD_N5786918  0.3588       gwas3
# 942 305236_A05_usgwas5502808 367426_C12_IBD_N5786918  0.3554       gwas3
#     pheno_FID1 F_MISS_FID1 cohort_FID2 pheno_FID2 F_MISS_FID2
# 941          1     0.01793       gwas3          2     0.03168
# 942          1     0.01677       gwas3          2     0.03168

# lack this combination:
kin[which((kin$FID1=="305235_C04_usgwas5502691" & kin$FID2=="305236_A05_usgwas5502808") | (kin$FID2=="305235_C04_usgwas5502691" & kin$FID1=="305236_A05_usgwas5502808")),]
# FID1                      ID1                     FID2
# 1771 305235_C04_usgwas5502691 305235_C04_usgwas5502691 305236_A05_usgwas5502808
# ID2 N_SNP HetHet   IBS0 HetConc HomIBS0 Kinship
# 1771 305236_A05_usgwas5502808  4215 0.4695 0.0043  0.5559  0.1374  0.3507

# OK third combination close to limit

########################################################################

# [1] 742
# [1] "Number of expected combinations: 1"
# [1] "Number of observed combinations: 2"
# [1] "305236_A05_usgwas5502808" "367426_C12_IBD_N5786918" 

n_possible_combinations
# [1] 1
dup[which((dup$FID1 %in% ids_tmp) | (dup$FID2 %in% ids_tmp)),]
#                         FID1                    FID2 Kinship cohort_FID1
# 941 305235_C04_usgwas5502691 367426_C12_IBD_N5786918  0.3588       gwas3
# 942 305236_A05_usgwas5502808 367426_C12_IBD_N5786918  0.3554       gwas3
#     pheno_FID1 F_MISS_FID1 cohort_FID2 pheno_FID2 F_MISS_FID2
# 941          1     0.01793       gwas3          2     0.03168
# 942          1     0.01677       gwas3          2     0.03168

# involves same group

########################################################################

# [1] 7395
# [1] "Number of expected combinations: 3"
# [1] "Number of observed combinations: 2"
# [1] "305235_C04_usgwas5502691" "305236_A05_usgwas5502808"
# [3] "367426_C12_IBD_N5786918" 

vec<-c("305235_C04_usgwas5502691","305236_A05_usgwas5502808","367426_C12_IBD_N5786918")

kin[which((kin$FID1 %in% vec) & (kin$FID2 %in% vec)),c(1,3,5:9)]
#                          FID1                     FID2 N_SNP HetHet   IBS0
# 1771 305235_C04_usgwas5502691 305236_A05_usgwas5502808  4215 0.4695 0.0043
# 1787 305235_C04_usgwas5502691  367426_C12_IBD_N5786918  4031 0.4738 0.0022
# 1820 305236_A05_usgwas5502808  367426_C12_IBD_N5786918  4036 0.4698 0.0025
#      HetConc HomIBS0
# 1771  0.5559  0.1374
# 1787  0.5712  0.1084
# 1820  0.5636  0.1111



dim(data_remove)
# [1] 13976     5

data_remove<-data_remove[which(!duplicated(data_remove$V1)),]

dim(data_remove)
# [1] 6599    5

table(data_remove$cohort)
# all_hce           chop_old_gwas   german_affy6_old_gwas 
# 46                     440                     888 
# german_illu550_old_gwas                   gwas1                   gwas2 
# 766                    2296                     522 
# italy_gsa kiel_austria_sibdcs_gsa         netherlands_gsa 
# 2                       3                       1 
# niddk_broad_gsa       niddk_cd_old_gwas     niddk_feinstein_gsa 
# 646                     662                      60 
# pittsburgh_gsa 
# 267 


table(data_remove$pheno)
# 1    2 
# 4169 2430 

y<-table(data_remove$cohort,data_remove$pheno)
#                            1    2
# all_hce                    1   16
# chop_old_gwas              6  473
# german_affy6_old_gwas    816   72
# german_illu550_old_gwas  374  391
# gwas1                   2153  120
# gwas2                    466   53
# italy_gsa                  0    2
# kiel_austria_sibdcs_gsa    0   19
# netherlands_gsa            0    1
# niddk_broad_gsa           81  565
# niddk_cd_old_gwas        259  403
# niddk_feinstein_gsa        1   59
# pittsburgh_gsa            11  256


data_remove_tmp<-data_remove[,c(1,1)]
colnames(data_remove_tmp)<-c("FID","IID")

write.table(y,paste(path,"pre_imputation/QC/relatedness/iibdgc_table_data_remove",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
write.table(data_remove,paste(path,"pre_imputation/QC/relatedness/iibdgc_data_remove",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
write.table(data_remove_tmp,paste(path,"pre_imputation/QC/relatedness/list_iibdgc_data_remove",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
write.table(data_inconsist,paste(path,"pre_imputation/QC/relatedness/iibdgc_data_inconsist",sep=""),col.names=T,row.names=T,sep="\t",quote=F)


#### manually inspect duplicated samples that, given where they come from, is suspicious that they are duplicates:

# 16 GWAS3-Kiel

dup[which((dup$cohort_FID1 %in% c("kiel_austria_sibdcs_gsa","all_hce")) & (dup$cohort_FID2 %in% c("kiel_austria_sibdcs_gsa","all_hce"))),]
# FID1                    FID2 Kinship cohort_FID1
# 938     305235_C04_usgwas5502691 367426_C12_IBD_N5786918  0.3572     all_hce
# 939     305236_A05_usgwas5502808 367426_C12_IBD_N5786918  0.3541     all_hce
# 996   333817_E03_1794STDY5084767              sample_id  0.5000     all_hce
# 1028  333830_C12_IBD_20135522173              sample_id  0.5000     all_hce
# 1068 351574_A05_SC_COLORS5537568              sample_id  0.5000     all_hce
# 1069 351574_B05_SC_COLORS5537569              sample_id  0.5000     all_hce
# 1070 351574_B06_SC_COLORS5537577              sample_id  0.5000     all_hce
# 1073 351574_F04_SC_COLORS5537565              sample_id  0.5000     all_hce
# 1074 351574_C06_SC_COLORS5537578              sample_id  0.5000     all_hce
# 1075 351574_F05_SC_COLORS5537573              sample_id  0.5000     all_hce
# 1076 351574_E05_SC_COLORS5537572              sample_id  0.5000     all_hce
# 1077 351574_C05_SC_COLORS5537570              sample_id  0.5000     all_hce
# 1078 351574_D05_SC_COLORS5537571              sample_id  0.5000     all_hce
# 1079 351574_D06_SC_COLORS5537579              sample_id  0.5000     all_hce
# 1081 351574_G05_SC_COLORS5537574              sample_id  0.5000     all_hce
# 1082 351574_H04_SC_COLORS5537567              sample_id  0.4999     all_hce
# 1083 351574_H05_SC_COLORS5537575              sample_id  0.5000     all_hce
# 1084 351574_G04_SC_COLORS5537566              sample_id  0.5000     all_hce
# pheno_FID1 F_MISS_FID1             cohort_FID2 pheno_FID2 F_MISS_FID2
# 938           1   1.615e-02                 all_hce          2   0.0309900
# 939           1   1.505e-02                 all_hce          2   0.0309900
# 996           2   3.846e-04 kiel_austria_sibdcs_gsa          2   0.0011900
# 1028          2   1.380e-04 kiel_austria_sibdcs_gsa          2   0.0139000
# 1068          2   6.335e-05 kiel_austria_sibdcs_gsa          2   0.0006199
# 1069          2   5.430e-05 kiel_austria_sibdcs_gsa          2   0.0005192
# 1070          2   3.394e-05 kiel_austria_sibdcs_gsa          2   0.0009908
# 1073          2   5.204e-05 kiel_austria_sibdcs_gsa          2   0.0007594
# 1074          2   1.041e-04 kiel_austria_sibdcs_gsa          2   0.0007894
# 1075          2   5.430e-05 kiel_austria_sibdcs_gsa          2   0.0007417
# 1076          2   9.729e-05 kiel_austria_sibdcs_gsa          2   0.0005369
# 1077          2   5.656e-05 kiel_austria_sibdcs_gsa          2   0.0008477
# 1078          2   8.824e-05 kiel_austria_sibdcs_gsa          2   0.0006799
# 1079          2   5.204e-05 kiel_austria_sibdcs_gsa          2   0.0009060
# 1081          2   6.335e-05 kiel_austria_sibdcs_gsa          2   0.0006075
# 1082          2   1.516e-04 kiel_austria_sibdcs_gsa          2   0.0009590
# 1083          2   8.824e-05 kiel_austria_sibdcs_gsa          2   0.0006111
# 1084          2   7.240e-05 kiel_austria_sibdcs_gsa          2   0.0013920




dup[which((dup$cohort_FID1 %in% c("niddk_cd_old_gwas","all_hce")) & (dup$cohort_FID2 %in% c("niddk_cd_old_gwas","all_hce"))),]
# FID1                       FID2 Kinship
# 243                   130903 333830_B07_IBD_20135522131  0.5000
# 938 305235_C04_usgwas5502691    367426_C12_IBD_N5786918  0.3572
# 939 305236_A05_usgwas5502808    367426_C12_IBD_N5786918  0.3541
# cohort_FID1 pheno_FID1 F_MISS_FID1 cohort_FID2 pheno_FID2 F_MISS_FID2
# 243 niddk_cd_old_gwas          2    0.004148     all_hce          2   6.788e-05
# 938           all_hce          1    0.016150     all_hce          2   3.099e-02
# 939           all_hce          1    0.015050     all_hce          2   3.099e-02


dup[which((dup$cohort_FID1 %in% c("kiel_austria_sibdcs_gsa","netherlands_gsa")) & (dup$cohort_FID2 %in% c("kiel_austria_sibdcs_gsa","netherlands_gsa"))),]
# FID1       FID2 Kinship     cohort_FID1 pheno_FID1 F_MISS_FID1
# 3867 sample_id sample_id     0.5 netherlands_gsa          2    0.001611
# cohort_FID2 pheno_FID2 F_MISS_FID2
# 3867 kiel_austria_sibdcs_gsa          2   0.0007541



##########################
# two GWAS3-GWAS3 pairs:

dup[which((dup$cohort_FID1 %in% c("all_hce")) & (dup$cohort_FID2 %in% c("all_hce"))),]
# FID1                    FID2 Kinship cohort_FID1
# 938 305235_C04_usgwas5502691 367426_C12_IBD_N5786918  0.3572     all_hce
# 939 305236_A05_usgwas5502808 367426_C12_IBD_N5786918  0.3541     all_hce
# pheno_FID1 F_MISS_FID1 cohort_FID2 pheno_FID2 F_MISS_FID2
# 938          1     0.01615     all_hce          2     0.03099
# 939          1     0.01505     all_hce          2     0.03099


# again part of the trio above: 367426_C12_IBD_N5786918 - 305235_C04_usgwas5502691 - 305236_A05_usgwas5502808
kin[which((kin$FID1 %in% vec) & (kin$FID2 %in% vec)),c(1,3,5:9)]
#                          FID1                     FID2 N_SNP HetHet   IBS0
# 1771 305235_C04_usgwas5502691 305236_A05_usgwas5502808  4215 0.4695 0.0043
# 1787 305235_C04_usgwas5502691  367426_C12_IBD_N5786918  4031 0.4738 0.0022
# 1820 305236_A05_usgwas5502808  367426_C12_IBD_N5786918  4036 0.4698 0.0025
#      HetConc HomIBS0
# 1771  0.5559  0.1374
# 1787  0.5712  0.1084
# 1820  0.5636  0.1111


##############################################################################################################################################

#################################################
# 4.- DOUBLE CHECK NO DUPLICATED SAMPLES REMAIN #
#################################################


studies=(australia_omniexome gwas1 gwas2 all_hce pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa chop_old_gwas german_illu550_old_gwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas niddk_cd_old_gwas)

for i in ${studies[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/relatedness/${i}_subset \
--remove ${path}pre_imputation/QC/relatedness/list_iibdgc_data_remove \
--make-bed \
--out ${path}pre_imputation/QC/relatedness/${i}_subset_nodup
done

wc -l ${path}pre_imputation/QC/relatedness/*_subset_nodup.fam
# 22511 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/all_hce_subset_nodup.fam
# 1287 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/australia_omniexome_subset_nodup.fam
# 1406 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_inf1_old_gwas_subset_nodup.fam
# 271 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_inf2_old_gwas_subset_nodup.fam
# 8026 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/chop_old_gwas_subset_nodup.fam
# 1894 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/german_affy6_old_gwas_subset_nodup.fam
# 846 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/german_illu550_old_gwas_subset_nodup.fam
# 2404 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/gwas1_subset_nodup.fam
# 7180 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/gwas2_subset_nodup.fam
# 977 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/italy_gsa_subset_nodup.fam
# 13696 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/kiel_austria_sibdcs_gsa_subset_nodup.fam
# 4547 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/netherlands_gsa_subset_nodup.fam
# 4737 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_broad_gsa_subset_nodup.fam
# 1060 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_cd_old_gwas_subset_nodup.fam
# 7901 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_feinstein_gsa_subset_nodup.fam
# 549 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/norway_affy6_old_gwas_subset_nodup.fam
# 2472 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/pittsburgh_gsa_subset_nodup.fam
# 257 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/slovenia_gsa_subset_nodup.fam
# 3387 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/spain_gsa_subset_nodup.fam
# 1334 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/sweden_gsa_subset_nodup.fam
# 86741 total

86742-93340
# [1] -6599 OK


#####

cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa"
           ,"chop_old_gwas","german_illu550_old_gwas","german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_cd_old_gwas") # illu 317

cohorts<-cohorts[-1]

dat<-matrix(ncol=3,nrow=length(cohorts))
dat<-as.data.frame(dat)

for (i in 1:nrow(dat)){
  dat[i,1]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_nodup.bed",sep="")
  dat[i,2]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_nodup.bim",sep="")
  dat[i,3]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_nodup.fam",sep="")
}

write.table(dat,paste(path,"pre_imputation/QC/relatedness/list_cohorts_tomerge_2.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#############################

# hpc-server

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/relatedness/australia_omniexome_subset_nodup \
--merge-list ${path}pre_imputation/QC/relatedness/list_cohorts_tomerge_2.txt \
--make-bed --out ${path}pre_imputation/QC/relatedness/iibdgc_merged_nodup

# 4608 variants and 86741 people pass filters and QC.
# Among remaining phenotypes, 47586 are cases and 39155 are controls.

  
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


# (93340*4608/4)/1E+6
# [1] 107.5277

path=/path/to/ibdgwas/IIBDGC/
MEM=200

bsub -J"kg" -M"$MEM" -n8 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path}pre_imputation/QC/relatedness/log/stderr_king_iibdgc_2 \
-o ${path}pre_imputation/QC/relatedness/log/stdout_king_iibdgc_2 \
"/path/to/software/./king \
-b ${path}pre_imputation/QC/relatedness/iibdgc_merged_nodup.bed \
--related --cpus 8 --prefix ${path}pre_imputation/QC/relatedness/iibdgc_merged_nodup"
# Job <517806> is submitted to queue <normal>.


################

kin<-read.table(paste(path,"pre_imputation/QC/relatedness/iibdgc_merged_nodup.kin0",sep=""),head=T)
table(cut(kin$Kinship,breaks=c(1,0.354,0.177)))
# (0.177,0.354]     (0.354,1] 
# 2520             2 

dup<-kin[which(kin$Kinship>0.354),c("FID1","FID2","Kinship")]

cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa"
           ,"chop_old_gwas","german_illu550_old_gwas","german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_cd_old_gwas")


for (i in 1:length(cohorts)) {
  
  file.tmp<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_nodup.fam",sep="")
  if(file.exists(file.tmp)){
    tmp<-read.table(file.tmp,head=F)
  }else{
    tmp<-read.table(paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset.fam",sep=""),head=F)
  }
  

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
  
  rm(tmp,file.tmp)
}


colnames(dat)[6]<-"pheno"
colnames(dat)[5]<-"sex"
dup<-merge(dup,dat[,c("V1","cohort","sex","pheno","F_MISS")],by.x="FID2",by.y="V1",all.x=T,sort=F)
colnames(dup)[(ncol(dup)-3):ncol(dup)]<-paste(colnames(dup)[(ncol(dup)-3):ncol(dup)],"FID2",sep="_")

dup<-merge(dup,dat[,c("V1","cohort","sex","pheno","F_MISS")],by.x="FID1",by.y="V1",all.x=T,sort=F)
colnames(dup)[(ncol(dup)-3):ncol(dup)]<-paste(colnames(dup)[(ncol(dup)-3):ncol(dup)],"FID1",sep="_")


aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary)
#                    Group.1   x.1   x.2
# 1                  all_hce 10429 12053
# 2      australia_omniexome   612   675
# 3    belgium_inf1_old_gwas   889   517
# 4    belgium_inf2_old_gwas   111   160
# 5            chop_old_gwas  6091  1974
# 6    german_affy6_old_gwas   933   961
# 7  german_illu550_old_gwas   756    89
# 8                    gwas1   777  1604
# 9                    gwas2  4894  2283
# 10               italy_gsa   385   592
# 11 kiel_austria_sibdcs_gsa  4454  9258
# 12         netherlands_gsa   547  4000
# 13         niddk_broad_gsa   930  3807
# 14       niddk_cd_old_gwas   672   388
# 15     niddk_feinstein_gsa  2350  5551
# 16   norway_affy6_old_gwas   281   268
# 17          pittsburgh_gsa  1474   998
# 18            slovenia_gsa   171    86
# 19               spain_gsa  1468  1919
# 20              sweden_gsa   931   403





######## DOUBLE CHECK REMAINING DUPLICATES:
dup
# FID1                    FID2 Kinship cohort_FID2 sex_FID2
# 1 305235_C04_usgwas5502691 367426_C12_IBD_N5786918  0.3572     all_hce        2
# 2 305236_A05_usgwas5502808 367426_C12_IBD_N5786918  0.3541     all_hce        2
# pheno_FID2 F_MISS_FID2 cohort_FID1 sex_FID1 pheno_FID1 F_MISS_FID1
# 1          2     0.03099     all_hce        2          1     0.01615
# 2          2     0.03099     all_hce        2          1     0.01505

# OK

