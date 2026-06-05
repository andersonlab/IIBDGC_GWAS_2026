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
  
studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)
  
for i in ${studies[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy \
--allow-no-sex \
--extract ${path}pre_imputation/QC/relatedness/list_common_variants_among_cohorts_iibdgc \
--make-bed \
--out ${path}pre_imputation/QC/relatedness/${i}_subset
done

# estimate missingness for later step
for i in ${studies[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy \
--allow-no-sex \
--missing \
--out ${path}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy
done


ls -la ${path}pre_imputation/QC/relatedness/*_subset.fam | wc -l 
# 34 OK

for i in ${studies[@]}
do wc -l ${path}pre_imputation/QC/relatedness/${i}_subset.fam
done
# 22426 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/all_hce_subset.fam
# 1572 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_old_gwas_subset.fam
# 1298 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/australia_omniexome_subset.fam
# 4676 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/gwas1_subset.fam
# 5132 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/gwas2_subset.fam
# 2725 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/pittsburgh_gsa_subset.fam
# 3443 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/spain_gsa_subset.fam
# 949 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/italy_gsa_subset.fam
# 13905 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/kiel_austria_sibdcs_gsa_subset.fam
# 4559 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/netherlands_gsa_subset.fam
# 262 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/slovenia_gsa_subset.fam
# 1379 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/sweden_gsa_subset.fam
# 5124 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_broad_gsa_subset.fam
# 7962 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_feinstein_gsa_subset.fam
# 1491 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/basque_gsa_subset.fam
# 2210 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/lithuania_gsa_subset.fam
# 1506 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_louis_gsa_subset.fam
# 1488 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_franchimont_gsa_subset.fam
# 3980 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_vermeire_gsa_subset.fam
# 466 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/prism_nfe_gsa_subset.fam
# 756 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/prism_nfe_gwas_subset.fam
# 444 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/finland_illugwas_subset.fam
# 2632 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/german_affy6_old_gwas_subset.fam
# 550 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/norway_affy6_old_gwas_subset.fam
# 1067 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_inf1_old_gwas_subset.fam
# 161 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_inf2_old_gwas_subset.fam
# 553 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_370k_old_gwas_subset.fam
# 767 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_610k_old_gwas_subset.fam
# 1036 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_omni_old_gwas_subset.fam
# 1248 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/swedish_uc_old_gwas_subset.fam
# 781 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/mccauley_gsa_subset.fam
# 2176 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/ccfa_gsa_subset.fam
# 3076 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_gsa_subset.fam
# 481 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/bernstein_gsa_subset.fam
# 68 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/farkkila_gsa_subset.fam
# 2727 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/franchimont_gsa_subset.fam
# 860 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/franke_gsa_subset.fam
# 760 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/helmsley_prism_gsa_subset.fam
# 1293 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/helmsley_xavier_prism_gsa_subset.fam
# 418 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/hyams_protect_gsa_subset.fam
# 2857 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/lewis_sparc_gsa_subset.fam
# 1612 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/mccauley_new_gsa_subset.fam
# 5960 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/mcgovern_gsa_subset.fam
# 1128 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/moayyedi_imagine_gsa_subset.fam
# 861 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/newberry_share_gsa_subset.fam
# 1746 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_cho_gsa_subset.fam
# 1817 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_duerr_gsa_subset.fam
# 913 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_rioux_gsa_subset.fam
# 2330 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_silverberg_gsa_subset.fam
# 878 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/palotie_hus_gsa_subset.fam
# 633 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/pekow_share_gsa_subset.fam
# 176 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/rioux_igenomed_gsa_subset.fam
# 1422 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/sands_msccr_gsa_subset.fam
# 1468 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/stampfer_gsa_subset.fam
# 4683 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/vermeire_gsa_subset.fam
# 709 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/weersma_gsa_subset.fam
# 687 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/xavier_prism_gsa_subset.fam
# 693 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/xavier_share_gsa_subset.fam

for i in ${studies[@]}
do
echo ${i} && head -10 ${path}pre_imputation/QC/relatedness/${i}_subset.fam
done
 


for i in ${studies[@]}
do wc -l ${path}pre_imputation/QC/relatedness/${i}_subset.bim
done
wc -l ${path}pre_imputation/QC/relatedness/*_subset.bim
# 3190  OK

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
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas"
           ,"swedish_uc_old_gwas","mccauley_gsa","ccfa_gsa","cedars_gsa")

length(cohorts)
# [1] 33 OK - no chop

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
--allow-no-sex \
--merge-list ${path}pre_imputation/QC/relatedness/list_cohorts_tomerge.txt \
--make-bed --out ${path}pre_imputation/QC/relatedness/iibdgc_merged

# 3190 variants and 101800 people pass filters and QC.
# Among remaining phenotypes, 62927 are cases and 38873 are controls.


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

# (3110*107598/4)/1E+6
# [1] 83.65744

path=/path/to/ibdgwas/IIBDGC/
MEM=3500

bsub -J"kg" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal -n 2 \
-e ${path}pre_imputation/QC/relatedness/log/stderr_king_iibdgc_plink2_noparallel \
-o ${path}pre_imputation/QC/relatedness/log/stdout_king_iibdgc_plink2_noparallel \
"/path/to/software/./plink2 \
--bfile ${path}pre_imputation/QC/relatedness/iibdgc_merged \
--make-king-table --king-table-filter 0.177 \
--memory $MEM \
--out ${path}pre_imputation/QC/relatedness/iibdgc_merged_plink"
# Job <984725> is submitted to queue <normal>.


################################

##### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas"
           ,"swedish_uc_old_gwas","mccauley_gsa","ccfa_gsa","cedars_gsa")

length(cohorts)
# [1] 33

for (i in 1:length(cohorts)) {
  
  tmp<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],
                        "_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy.fam",sep=""),head=F)

  
  sample_miss<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],
                                "_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy.imiss",sep=""),head=T)
  
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
# all_hce     australia_omniexome              basque_gsa 
# 22426                    1298                    1491 
# belgium_franchimont_gsa   belgium_inf1_old_gwas   belgium_inf2_old_gwas 
# 1488                    1067                     161 
# belgium_louis_gsa    belgium_vermeire_gsa                ccfa_gsa 
# 1506                    3980                    2176 
# cedars_370k_old_gwas    cedars_610k_old_gwas              cedars_gsa 
# 553                     767                    3076 
# cedars_omni_old_gwas        finland_illugwas   german_affy6_old_gwas 
# 1036                     444                    2632 
# gwas1                   gwas2               italy_gsa 
# 4676                    5132                     949 
# kiel_austria_sibdcs_gsa           lithuania_gsa            mccauley_gsa 
# 13905                    2210                     781 
# netherlands_gsa         niddk_broad_gsa     niddk_feinstein_gsa 
# 4559                    5124                    7962 
# niddk_old_gwas   norway_affy6_old_gwas          pittsburgh_gsa 
# 1572                     550                    2725 
# prism_nfe_gsa          prism_nfe_gwas            slovenia_gsa 
# 466                     756                     262 
# spain_gsa              sweden_gsa     swedish_uc_old_gwas 
# 3443                    1379                    1248 

# double check all samples included in exercise
famall<-read.table(paste(path,"pre_imputation/QC/relatedness/iibdgc_merged.fam",sep=""),head=F)

dim(famall)
# [1] 101800     6

dim(dat)
# [1] 101800     8

dim(dat[which(dat$V1 %in% famall$V1),])
# [1] 101800     8


colnames(dat)[6]<-"pheno"

aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary)
#                    Group.1   x.1   x.2
# 1                  all_hce 10444 11982
# 2      australia_omniexome   612   686
# 3               basque_gsa   966   525
# 4  belgium_franchimont_gsa   592   896
# 5    belgium_inf1_old_gwas   839   228
# 6    belgium_inf2_old_gwas    77    84
# 7        belgium_louis_gsa   599   907
# 8     belgium_vermeire_gsa   813  3167
# 9                 ccfa_gsa     0  2176
# 10    cedars_370k_old_gwas     0   553
# 11    cedars_610k_old_gwas     0   767
# 12              cedars_gsa   947  2129
# 13    cedars_omni_old_gwas     2  1034
# 14        finland_illugwas     0   444
# 15   german_affy6_old_gwas  1668   964
# 16                   gwas1  2929  1747
# 17                   gwas2  2776  2356
# 18               italy_gsa   374   575
# 19 kiel_austria_sibdcs_gsa  4455  9450
# 20           lithuania_gsa  1139  1071
# 21            mccauley_gsa     0   781
# 22         netherlands_gsa   548  4011
# 23         niddk_broad_gsa  1012  4112
# 24     niddk_feinstein_gsa  2350  5612
# 25          niddk_old_gwas   728   844
# 26   norway_affy6_old_gwas   282   268
# 27          pittsburgh_gsa  1468  1257
# 28           prism_nfe_gsa    33   433
# 29          prism_nfe_gwas   257   499
# 30            slovenia_gsa   175    87
# 31               spain_gsa  1482  1961
# 32              sweden_gsa   979   400
# 33     swedish_uc_old_gwas   327   921

colnames(dat)[5]<-"sex"

aggregate(as.factor(dat$sex), by=list(dat$cohort), FUN=summary)
#                    Group.1   x.0   x.1   x.2
# 1                  all_hce    83 10350 11993
# 2      australia_omniexome     1   602   695
# 3               basque_gsa     7   912   572
# 4  belgium_franchimont_gsa     2   642   844
# 5    belgium_inf1_old_gwas    11   848   208
# 6    belgium_inf2_old_gwas     1    92    68
# 7        belgium_louis_gsa    11   650   845
# 8     belgium_vermeire_gsa    26  1904  2050
# 9                 ccfa_gsa     2   997  1177
# 10    cedars_370k_old_gwas     1   298   254
# 11    cedars_610k_old_gwas     3   403   361
# 12              cedars_gsa     9  1481  1586
# 13    cedars_omni_old_gwas     0   516   520
# 14        finland_illugwas     0   295   149
# 15   german_affy6_old_gwas     4  1270  1358
# 16                   gwas1     2  2121  2553
# 17                   gwas2    51  2550  2531
# 18               italy_gsa     3   515   431
# 19 kiel_austria_sibdcs_gsa   193  6804  6908
# 20           lithuania_gsa     0  1174  1036
# 21            mccauley_gsa     2   400   379
# 22         netherlands_gsa    11  1895  2653
# 23         niddk_broad_gsa    16  2651  2457
# 24     niddk_feinstein_gsa    13  4001  3948
# 25          niddk_old_gwas     0   780   792
# 26   norway_affy6_old_gwas     1   304   245
# 27          pittsburgh_gsa     8  1183  1534
# 28           prism_nfe_gsa     4   216   246
# 29          prism_nfe_gwas     0   367   389
# 30            slovenia_gsa     5    55   202
# 31               spain_gsa    56  1883  1504
# 32              sweden_gsa     0   797   582
# 33     swedish_uc_old_gwas     0   688   560

table(dat$V1==dat$V2)
# TRUE 
# 101800


kin<-read.table(paste(path,"pre_imputation/QC/relatedness/iibdgc_merged_plink.kin0",sep=""),head=F)
colnames(kin)<-c("FID1","IID1","FID2","IID2","NSNP","HETHET","IBS0","KINSHIP")

table(cut(kin$KINSHIP,breaks=c(1,0.354,0.177)))
# (0.177,0.354]     (0.354,1] 
#        16924          3553 

pdf(paste(path,"pre_imputation/QC/relatedness/histogram_iibdgc_merged_plink.pdf",sep=""),height = 7,width = 14)
h <- hist(kin$KINSHIP, breaks = "FD", plot = FALSE)
ggplot(kin, aes(x=KINSHIP)) + geom_histogram(breaks = h$breaks) + geom_vline(xintercept=0.45, linetype="dashed", color = "red")
dev.off()

system(paste("cp ",path,"pre_imputation/QC/relatedness/histogram_iibdgc_merged_plink.pdf ~/tmp_plots/",sep=""))

dup<-kin[which(kin$KINSHIP>0.45),c("FID1","FID2","KINSHIP")]
summary(dup$KINSHIP)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.4771  0.5000  0.5000  0.4999  0.5000  0.5000 


# https://en.wikipedia.org/wiki/Freedman–Diaconis_rule

dup<-merge(dup,dat[,c("V1","cohort","sex","pheno","F_MISS")],by.x="FID2",by.y="V1",all.x=T,sort=F)
colnames(dup)[(ncol(dup)-3):ncol(dup)]<-paste(colnames(dup)[(ncol(dup)-3):ncol(dup)],"FID2",sep="_")

dup<-merge(dup,dat[,c("V1","cohort","sex","pheno","F_MISS")],by.x="FID1",by.y="V1",all.x=T,sort=F)
colnames(dup)[(ncol(dup)-3):ncol(dup)]<-paste(colnames(dup)[(ncol(dup)-3):ncol(dup)],"FID1",sep="_")

x<-table(dup$cohort_FID1,dup$cohort_FID2)
x



xx<-dup[,c("cohort_FID1","cohort_FID2")]
yy<-dup[,c("cohort_FID2","cohort_FID1")]
colnames(yy)<-colnames(xx)
  
dim(xx)
# [1] 3199    2

xx<-rbind(xx,yy)

x<-table(xx$cohort_FID1,xx$cohort_FID2)

write.table(x,paste(path,"pre_imputation/QC/relatedness/iibdgc_dup_summary_table",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
system(paste("cp ",path,"pre_imputation/QC/relatedness/iibdgc_dup_summary_table ~/tmp_plots/",sep=""))

table(dup$sex_FID1,dup$sex_FID2)
#      0    1    2
# 0    2    1    3
# 1    0 1587    0
# 2    7    0 1599


dup<-dup[,c("FID1","FID2","KINSHIP","cohort_FID1","pheno_FID1","F_MISS_FID1","cohort_FID2","pheno_FID2","F_MISS_FID2")]

table(dup$pheno_FID1,dup$pheno_FID2)
#      1    2
# 1  332    1
# 2    1 2865

dup_ids<-c(as.character(dup$FID1),as.character(dup$FID2))
length(dup_ids)
#[1] 6398
table(length(dup_ids)==(nrow(dup)*2))
# TRUE 
# 1 

dup_ids<-dup_ids[!duplicated(dup_ids)]
length(dup_ids)
#[1] 6230

dim(dup)
# [1] 3199    9

# CHECK ALSO THAT BOTH SAMPLES HAVE SAME PHENOTYPE AND GENDER, 
# - keep all new gwas samples

# - for ukbb-ukbb duplicated pairs:
#    - keep sample with pheno data when both samples have the same phenotype, otherwise exclude both

new_gwas<-c("spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa","all_hce"
            ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa"
            ,"basque_gsa","prism_nfe_gsa","lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
            ,"mccauley_gsa","ccfa_gsa","cedars_gsa")

old_gwas<-c("pittsburgh_gsa","australia_omniexome","gwas1","gwas2"
            ,"german_affy6_old_gwas","norway_affy6_old_gwas"
            ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas",
            "prism_nfe_gwas","finland_illugwas"
            ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","swedish_uc_old_gwas")

# classify old gwas into good and bad:

good_array<-c("pittsburgh_gsa","australia_omniexome","german_illu550_old_gwas","german_affy6_old_gwas","norway_affy6_old_gwas","australia_omniexome"
              ,"gwas1","gwas2","cedars_610k_old_gwas","cedars_omni_old_gwas")
bad_array<-c("belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas","prism_nfe_gwas","finland_illugwas","chop_old_gwas","cedars_370k_old_gwas"
             ,"swedish_uc_old_gwas")

cohorts[which(!cohorts %in% c(old_gwas,new_gwas))] 
# 0

old_gwas[which(!old_gwas %in% c(good_array,bad_array))] 
# 0

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
    
    
    # OPTION A: samples have different pheno, remove all (already checked before all samples have same sex:
    
    if ( (dim(table(data$pheno))!=1) ) {
      
      # # make an exeption to deal with sex = 0
      # if ( any(data$sex==0) & any(data$sex %in% c(1,2)) & dim(table(data$pheno))==1 ) {
      #   
      #   keep_sample<-data$V1[which(data$sex!=0)]
      #   
      #   if(!exists("data_remove")) {
      #     data_remove<-data[which(!data$V1 %in% keep_sample),]
      #   } else {
      #     data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
      #   }
      #   
      # } else 
      if (any(data$cohort %in% c("chop_old_gwas"))){
        
        data<-data[which(data$cohort=="chop_old_gwas"),]
        
        if (!exists("data_remove")) {
          data_remove<-data
        } else {
          data_remove<-rbind(data_remove,data)
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
    
    else {
    
    # OPTION B: samples have same pheno:
    
      # B.1: rest of samples, if samples in new and old arrays, keep only samples in new arrays:
      
      if (any(data$cohort %in% old_gwas) & any(data$cohort %in% new_gwas)) {
        
        keep_sample<-data$V1[which(!data$cohort %in% old_gwas)]
        
        # if more than one new pair in new gwas, keep sample with largest call rate, make exception for CCFA and McCauley:
        
        if (length(keep_sample)>1) {
          
          if ( all(data$cohort %in% c("ccfa_gsa","cedars_gsa")) | all(data$cohort %in% c("ccfa_gsa","mccauley_gsa"))) {
            
            keep_sample<-data$V1[which(data$F_MISS==min(data$F_MISS))][1]# to deal with >1 with same F_miss
            
          } 
          
          else if ( any(data$cohort %in% c("ccfa_gsa","mccauley_gsa","cedars_gsa")) ) {
            
            keep_sample<-data[which(!(data$cohort %in% old_gwas) & !(data$cohort %in% c("ccfa_gsa","mccauley_gsa","cedars_gsa"))),]
            keep_sample<-keep_sample$V1[which(keep_sample$F_MISS==min(keep_sample$F_MISS))][1]# to deal with >1 with same F_miss
            
          } 
          else {
            
            keep_sample<-data[which(!data$cohort %in% old_gwas),]
            keep_sample<-keep_sample$V1[which(keep_sample$F_MISS==min(keep_sample$F_MISS))][1]# to deal with >1 with same F_miss
            
          }
          
        }
        
        if (!exists("data_remove")) {
          data_remove<-data[which(!data$V1 %in% keep_sample),]
        } else {
          data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
        }
  
      } 
      
      # B.2: old_gwas duplicated pairs only
      
      else if (all(data$cohort %in% old_gwas)) {
        
        # B.2.1: if in good and in bad array, keep sample in the best array:
        
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
  
        # B.2.2: if all in good OR all in bad array, keep sample in the best array:
        
        else if (all(data$cohort %in% good_array) | all(data$cohort %in% bad_array)) {
          
          keep_sample<-data$V1[which(data$F_MISS==min(data$F_MISS))][1]# to deal with same F_miss
          
          if(!exists("data_remove")) {
            data_remove<-data[which(!data$V1 %in% keep_sample),]
          } else {
            data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
          }
          
        }
      }
      
      # B.3: new_gwas duplicated pairs only, keep sample with largest call rate: 
      
      else if (all(data$cohort %in% new_gwas)) {
        
        # B.3.1: Samples in ccfa (no ctr) and other new cohort, keep other:
        
        if ( all(data$cohort %in% c("ccfa_gsa","cedars_gsa")) | all(data$cohort %in% c("ccfa_gsa","mccauley_gsa"))) {
          
          keep_sample<-data$V1[which(data$F_MISS==min(data$F_MISS))][1]# to deal with >1 with same F_miss
          
          if(!exists("data_remove")) {
            data_remove<-data[which(!data$V1 %in% keep_sample),]
          } else {
            data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
          }
          
        } else if ( any(data$cohort %in% c("ccfa_gsa","mccauley_gsa","cedars_gsa")) ) {
          
          keep_sample<-data$V1[which(!data$cohort %in% c("ccfa_gsa","mccauley_gsa","cedars_gsa"))]
          
          # if more than one new pair in new gwas, keep sample with largest call rate:
          if (length(keep_sample)>1) {
            
            keep_sample<-data[which(!data$cohort %in% c("ccfa_gsa","mccauley_gsa","cedars_gsa")),]
            keep_sample<-keep_sample$V1[which(keep_sample$F_MISS==min(keep_sample$F_MISS))][1]# to deal with >1 with same F_miss
            
          }
          
          if(!exists("data_remove")) {
            data_remove<-data[which(!data$V1 %in% keep_sample),]
          } else {
            data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
          }
        } 
        
        # B.3.2: Samples in kiel and gwas3 OR belgium:
        
        else if ( all( any(data$cohort=="all_hce" | data$cohort=="belgium_vermeire_gsa") & any(data$cohort %in% c("kiel_austria_sibdcs_gsa")) )  ) {
          
          keep_sample<-data$V1[which(!data$cohort %in% c("kiel_austria_sibdcs_gsa"))]
          
          # if more than one new pair in new gwas, keep sample with largest call rate:
          if (length(keep_sample)>1) {
            keep_sample<-data[which(!data$cohort %in% c("kiel_austria_sibdcs_gsa")),]
            keep_sample<-keep_sample$V1[which(keep_sample$F_MISS==min(keep_sample$F_MISS))][1]# to deal with >1 with same F_miss
          }
          
          if(!exists("data_remove")) {
            data_remove<-data[which(!data$V1 %in% keep_sample),]
          } else {
            data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
          }
        } 
        
        # B.3.3: Other keep sample with largest call rate:
        
        else {
          
          keep_sample<-data$V1[which(data$F_MISS==min(data$F_MISS))]
          
          if(!exists("data_remove")) {
            data_remove<-data[which(!data$V1 %in% keep_sample),]
          } else {
            data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
          }
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
#[1] 6408  4

data_remove<-data_remove[!duplicated(data_remove$V1),]
dim(data_remove)
# [1] 3148    5
# double check which sample per pair removed in Kiel_Uk pairs:

table(data_remove$cohort)
# belgium_vermeire_gsa                ccfa_gsa    cedars_370k_old_gwas 
#                   1                     112                     100 
# cedars_610k_old_gwas              cedars_gsa    cedars_omni_old_gwas 
#                  300                    1273                     363 
# german_affy6_old_gwas kiel_austria_sibdcs_gsa            mccauley_gsa 
#                     1                      35                      12 
# netherlands_gsa         niddk_broad_gsa     niddk_feinstein_gsa 
#               4                       2                       1 
# niddk_old_gwas          pittsburgh_gsa           prism_nfe_gsa 
#            226                     695                       2 
# prism_nfe_gwas 
#             21 

# 
# table(dup[which(dup$cohort_FID1 %in% c("ccfa_gsa","mccauley_gsa") | dup$cohort_FID2 %in% c("ccfa_gsa","mccauley_gsa")),"cohort_FID1"],
#       dup[which(dup$cohort_FID1 %in% c("ccfa_gsa","mccauley_gsa") | dup$cohort_FID2 %in% c("ccfa_gsa","mccauley_gsa")),"cohort_FID2"])
# #                      ccfa_gsa chop_old_gwas mccauley_gsa niddk_broad_gsa
# # ccfa_gsa                    0            19            0              56
# # cedars_370k_old_gwas        1             0            1               0
# # cedars_610k_old_gwas        1             0            0               0
# # cedars_gsa                  4             0            5               0
# # cedars_omni_old_gwas        0             0            1               0
# # mccauley_gsa                0             2            0               1
# # 
# #                      niddk_feinstein_gsa niddk_old_gwas pittsburgh_gsa
# # ccfa_gsa                              47              3             45
# # cedars_370k_old_gwas                   0              0              0
# # cedars_610k_old_gwas                   0              0              0
# # cedars_gsa                             0              0              0
# # cedars_omni_old_gwas                   0              0              0
# # mccauley_gsa                           8              0              0
# # 
# #                       prism_nfe_gsa prism_nfe_gwas
# # ccfa_gsa                         9             13
# # cedars_370k_old_gwas             0              0
# # cedars_610k_old_gwas             0              0
# # cedars_gsa                       0              0
# # cedars_omni_old_gwas             0              0
# # mccauley_gsa                     0              2
# 
# tmp<-data_remove[which(data_remove$cohort=="niddk_broad_gsa"),]
# dup[which(dup$FID1 %in% tmp$V1 | dup$FID2 %in% tmp$V1),]


#### create a file to report inconsitend data:

dim(data_inconsist)
# [1] 8 5
data_inconsist<-data_inconsist[!duplicated(data_inconsist$V1),]
dim(data_inconsist)
# [1] 4  6

data_inconsist
#           V1 sex pheno          cohort    F_MISS group
# 1 sample_id   1     1  pittsburgh_gsa 0.0006155     1
# 2 sample_id   1     2 niddk_broad_gsa 0.0011710     1
# 3 sample_id   1     2  pittsburgh_gsa 0.0006611     2
# 4     160284   1     1  niddk_old_gwas 0.0065020     2


data_inconsist$cohort[which(data_inconsist$cohort=="pittsburgh_gsa")]<-"pittsburgh"


# which niddk_old_gwas is sample "160284" from?
tmp1<-read.table(paste(path,"pre_imputation/QC/niddk_cd_old_gwas/niddk_cd_old_gwas_hg17.fam",sep=""),head=F)
tmp2<-read.table(paste(path,"pre_imputation/QC/niddk_uc_old_gwas/niddk_uc_old_gwas_hg18.fam",sep=""),head=F)

tmp1[which(tmp1$V1=="160284"),]
# V1     V2 V3 V4 V5 V6
# 638 160284 160284  0  0  1  1

tmp2[which(tmp2$V1=="160284"),]
# [1] V1 V2 V3 V4 V5 V6
# <0 rows> (or 0-length row.names)

data_inconsist$cohort[which(data_inconsist$cohort=="niddk_old_gwas")]<-"niddk_cd_old_gwas"


### save files:

data_remove<-data_remove[,c(1,1)]
colnames(data_remove)<-c("FID","IID")


# these two are not anymore in data inconsist because sample_id has been excluded from CEDARS GSA
# 5 sample_id   2     1          cedars_gsa 0.0005495     3
# 6 sample_id   2     2 niddk_feinstein_gsa 0.0002400     3


# but from Talin's message, we should also exclude sample_id:
data_remove[nrow(data_remove)+1,]<-c("sample_id","sample_id")

write.table(data_remove,paste(path,"pre_imputation/QC/relatedness/list_iibdgc_data_remove",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

colnames(data_inconsist)[1]<-c("IIBDGC_id")
write.table(data_inconsist[,c("IIBDGC_id","sex","pheno","cohort","group")],paste(path,"pre_imputation/QC/relatedness/all_iibdgc_list_duplicated_samples_inconsistent_phenotype.tsv",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")

# create for Talin and kyle list with duplicates between NIDDK and cedars:

dup1<-dup[which((dup$cohort_FID1 %in% c("cedars_gsa")) | (dup$cohort_FID2 %in% c("cedars_gsa"))),]
dup1<-dup1[which((dup1$cohort_FID1 %in% c("niddk_feinstein_gsa","niddk_broad_gsa")) | (dup1$cohort_FID2 %in% c("niddk_feinstein_gsa","niddk_broad_gsa"))),]

write.table(dup1,paste(path,"pre_imputation/QC/relatedness/all_iibdgc_list_duplicated_pairs_inconsistent_phenotype_cedars_gsa_niddk.tsv",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")
system(paste("cp ",path,"pre_imputation/QC/relatedness/all_iibdgc_list_duplicated_pairs_inconsistent_phenotype_cedars_gsa_niddk.tsv ~/files_iibdgc/",sep=""))




# #### manually inspect duplicated samples that, given where they come from, is suspicious that they are duplicates:
# 
# # 16 GWAS3-Kiel
# 
# dup[which((dup$cohort_FID1 %in% c("kiel_austria_sibdcs_gsa","all_hce")) & (dup$cohort_FID2 %in% c("kiel_austria_sibdcs_gsa","all_hce"))),]
# # FID1                    FID2 Kinship cohort_FID1
# # 938     305235_C04_usgwas5502691 367426_C12_IBD_N5786918  0.3572     all_hce
# # 939     305236_A05_usgwas5502808 367426_C12_IBD_N5786918  0.3541     all_hce
# # 996   333817_E03_1794STDY5084767              sample_id  0.5000     all_hce
# # 1028  333830_C12_IBD_20135522173              sample_id  0.5000     all_hce
# # 1068 351574_A05_SC_COLORS5537568              sample_id  0.5000     all_hce
# # 1069 351574_B05_SC_COLORS5537569              sample_id  0.5000     all_hce
# # 1070 351574_B06_SC_COLORS5537577              sample_id  0.5000     all_hce
# # 1073 351574_F04_SC_COLORS5537565              sample_id  0.5000     all_hce
# # 1074 351574_C06_SC_COLORS5537578              sample_id  0.5000     all_hce
# # 1075 351574_F05_SC_COLORS5537573              sample_id  0.5000     all_hce
# # 1076 351574_E05_SC_COLORS5537572              sample_id  0.5000     all_hce
# # 1077 351574_C05_SC_COLORS5537570              sample_id  0.5000     all_hce
# # 1078 351574_D05_SC_COLORS5537571              sample_id  0.5000     all_hce
# # 1079 351574_D06_SC_COLORS5537579              sample_id  0.5000     all_hce
# # 1081 351574_G05_SC_COLORS5537574              sample_id  0.5000     all_hce
# # 1082 351574_H04_SC_COLORS5537567              sample_id  0.4999     all_hce
# # 1083 351574_H05_SC_COLORS5537575              sample_id  0.5000     all_hce
# # 1084 351574_G04_SC_COLORS5537566              sample_id  0.5000     all_hce
# # pheno_FID1 F_MISS_FID1             cohort_FID2 pheno_FID2 F_MISS_FID2
# # 938           1   1.615e-02                 all_hce          2   0.0309900
# # 939           1   1.505e-02                 all_hce          2   0.0309900
# # 996           2   3.846e-04 kiel_austria_sibdcs_gsa          2   0.0011900
# # 1028          2   1.380e-04 kiel_austria_sibdcs_gsa          2   0.0139000
# # 1068          2   6.335e-05 kiel_austria_sibdcs_gsa          2   0.0006199
# # 1069          2   5.430e-05 kiel_austria_sibdcs_gsa          2   0.0005192
# # 1070          2   3.394e-05 kiel_austria_sibdcs_gsa          2   0.0009908
# # 1073          2   5.204e-05 kiel_austria_sibdcs_gsa          2   0.0007594
# # 1074          2   1.041e-04 kiel_austria_sibdcs_gsa          2   0.0007894
# # 1075          2   5.430e-05 kiel_austria_sibdcs_gsa          2   0.0007417
# # 1076          2   9.729e-05 kiel_austria_sibdcs_gsa          2   0.0005369
# # 1077          2   5.656e-05 kiel_austria_sibdcs_gsa          2   0.0008477
# # 1078          2   8.824e-05 kiel_austria_sibdcs_gsa          2   0.0006799
# # 1079          2   5.204e-05 kiel_austria_sibdcs_gsa          2   0.0009060
# # 1081          2   6.335e-05 kiel_austria_sibdcs_gsa          2   0.0006075
# # 1082          2   1.516e-04 kiel_austria_sibdcs_gsa          2   0.0009590
# # 1083          2   8.824e-05 kiel_austria_sibdcs_gsa          2   0.0006111
# # 1084          2   7.240e-05 kiel_austria_sibdcs_gsa          2   0.0013920


##############################################################################################################################################

#################################################
# 4.- DOUBLE CHECK NO DUPLICATED SAMPLES REMAIN #
#################################################

path=/path/to/ibdgwas/IIBDGC/
  
studies=(australia_omniexome gwas1 gwas2 all_hce pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas chop_old_gwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas niddk_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa)

for i in ${studies[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/relatedness/${i}_subset \
--remove ${path}pre_imputation/QC/relatedness/list_iibdgc_data_remove \
--make-bed \
--out ${path}pre_imputation/QC/relatedness/${i}_subset_nodup
done

wc -l ${path}pre_imputation/QC/relatedness/*_subset_nodup.fam
# 22426 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/all_hce_subset_nodup.fam
# 1298 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/australia_omniexome_subset_nodup.fam
# 1491 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/basque_gsa_subset_nodup.fam
# 1489 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_franchimont_gsa_subset_nodup.fam
# 1070 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_inf1_old_gwas_subset_nodup.fam
# 161 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_inf2_old_gwas_subset_nodup.fam
# 1490 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_louis_gsa_subset_nodup.fam
# 3979 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_vermeire_gsa_subset_nodup.fam
# 2064 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/ccfa_gsa_subset_nodup.fam
# 452 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_370k_old_gwas_subset_nodup.fam
# 467 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_610k_old_gwas_subset_nodup.fam
# 1802 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_gsa_subset_nodup.fam
# 673 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_omni_old_gwas_subset_nodup.fam
# 7895 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/chop_old_gwas_subset_nodup.fam
# 444 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/finland_illugwas_subset_nodup.fam
# 2631 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/german_affy6_old_gwas_subset_nodup.fam
# 4676 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/gwas1_subset_nodup.fam
# 5132 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/gwas2_subset_nodup.fam
# 949 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/italy_gsa_subset_nodup.fam
# 13870 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/kiel_austria_sibdcs_gsa_subset_nodup.fam
# 2210 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/lithuania_gsa_subset_nodup.fam
# 769 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/mccauley_gsa_subset_nodup.fam
# 4555 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/netherlands_gsa_subset_nodup.fam
# 5122 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_broad_gsa_subset_nodup.fam
# 7960 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_feinstein_gsa_subset_nodup.fam
# 1284 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_old_gwas_subset_nodup.fam
# 550 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/norway_affy6_old_gwas_subset_nodup.fam
# 2030 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/pittsburgh_gsa_subset_nodup.fam
# 464 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/prism_nfe_gsa_subset_nodup.fam
# 735 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/prism_nfe_gwas_subset_nodup.fam
# 262 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/slovenia_gsa_subset_nodup.fam
# 3443 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/spain_gsa_subset_nodup.fam
# 1334 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/sweden_gsa_subset_nodup.fam
# 1248 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/swedish_uc_old_gwas_subset_nodup.fam
# 106425 total

110318-(3892+1)
# [1] 106425 OK


########################################################################################################################
# create a new merged set:

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
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"chop_old_gwas","german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas"
           ,"swedish_uc_old_gwas","mccauley_gsa","ccfa_gsa","cedars_gsa")

length(cohorts)
# [1] 34

#####

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
--allow-no-sex \
--merge-list ${path}pre_imputation/QC/relatedness/list_cohorts_tomerge_2.txt \
--make-bed --out ${path}pre_imputation/QC/relatedness/iibdgc_merged_nodup
# 3111 variants and 106425 people pass filters and QC.
# Among remaining phenotypes, 61780 are cases and 44645 are controls.

#############################
# 
# path<-"/path/to/ibdgwas/IIBDGC/"
# 
# cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
#            ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
#            ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
#            ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
#            ,"prism_nfe_gwas","finland_illugwas"
#            ,"chop_old_gwas","german_affy6_old_gwas","norway_affy6_old_gwas"
#            ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
#            ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas"
#            ,"swedish_uc_old_gwas","mccauley_gsa","ccfa_gsa","cedars_gsa")
# 
# length(cohorts)
# # [1] 34
# 
# for (i in 1:length(cohorts)) {
#   
#   tmp<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],
#                         "_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy.fam",sep=""),head=F)
#   
#   
#   sample_miss<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],
#                                 "_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy.imiss",sep=""),head=T)
#   
#   tmp<-merge(tmp,sample_miss[,c("IID","F_MISS")],by.x="V2",by.y="IID",all.x=T)
#   
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
# # all_hce     australia_omniexome              basque_gsa 
# # 22426                    1298                    1491 
# # belgium_franchimont_gsa   belgium_inf1_old_gwas   belgium_inf2_old_gwas 
# # 1489                    1070                     161 
# # belgium_louis_gsa    belgium_vermeire_gsa                ccfa_gsa 
# # 1490                    3980                    2176 
# # cedars_370k_old_gwas    cedars_610k_old_gwas              cedars_gsa 
# # 553                     767                    3075 
# # cedars_omni_old_gwas           chop_old_gwas        finland_illugwas 
# # 1036                    8576                     444 
# # german_affy6_old_gwas                   gwas1                   gwas2 
# # 2632                    4676                    5132 
# # italy_gsa kiel_austria_sibdcs_gsa           lithuania_gsa 
# # 949                   13905                    2210 
# # mccauley_gsa         netherlands_gsa         niddk_broad_gsa 
# # 781                    4559                    5124 
# # niddk_feinstein_gsa          niddk_old_gwas   norway_affy6_old_gwas 
# # 7962                    1572                     550 
# # pittsburgh_gsa           prism_nfe_gsa          prism_nfe_gwas 
# # 2725                     466                     756 
# # slovenia_gsa               spain_gsa              sweden_gsa 
# # 262                    3443                    1334 
# # swedish_uc_old_gwas 
# # 1248 
# 
# 
# # double check all samples included in exercise
# famall<-read.table(paste(path,"pre_imputation/QC/relatedness/iibdgc_merged_nodup.fam",sep=""),head=F)
# 
# dim(famall)
# # [1] 105670     6
# 
# dim(dat)
# # [1] 109371     8
# 
# dat<-dat[which(dat$V1 %in% famall$V1),]
# dim(dat)
# # [1] 105670     8
# 
# colnames(dat)[6]<-"pheno"
# 
# aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary)
# # Group.1   x.1   x.2
# # 1                  all_hce 10444 11982
# # 2      australia_omniexome   612   686
# # 3               basque_gsa   966   525
# # 4  belgium_franchimont_gsa   592   897
# # 5    belgium_inf1_old_gwas   839   231
# # 6    belgium_inf2_old_gwas    77    84
# # 7        belgium_louis_gsa   595   895
# # 8     belgium_vermeire_gsa   813  3166
# # 9                 ccfa_gsa     0  2064
# # 10    cedars_370k_old_gwas     0   452
# # 11    cedars_610k_old_gwas     0   467
# # 12              cedars_gsa   762  1046
# # 13    cedars_omni_old_gwas     2   671
# # 14           chop_old_gwas  6158  1737
# # 15        finland_illugwas     0   444
# # 16   german_affy6_old_gwas  1668   963
# # 17                   gwas1  2929  1747
# # 18                   gwas2  2776  2356
# # 19               italy_gsa   374   575
# # 20 kiel_austria_sibdcs_gsa  4455  9415
# # 21           lithuania_gsa  1139  1071
# # 22            mccauley_gsa     0   769
# # 23         netherlands_gsa   548  4007
# # 24         niddk_broad_gsa  1012  4111
# # 25     niddk_feinstein_gsa  2352  5606
# # 26          niddk_old_gwas   674   610
# # 27   norway_affy6_old_gwas   282   268
# # 28          pittsburgh_gsa  1378   652
# # 29           prism_nfe_gsa    33   431
# # 30          prism_nfe_gwas   257   478
# # 31            slovenia_gsa   175    87
# # 32               spain_gsa  1482  1961
# # 33              sweden_gsa   931   403
# # 34     swedish_uc_old_gwas   327   921
# 
# 
# 
# 
# kin<-read.table(paste(path,"pre_imputation/QC/relatedness/iibdgc_merged_plink.kin0",sep=""),head=F)
# colnames(kin)<-c("FID1","IID1","FID2","IID2","NSNP","HETHET","IBS0","KINSHIP")
# 
# table(cut(kin$KINSHIP,breaks=c(1,0.45,0.177)))
# # (0.177,0.354]     (0.354,1] 
# #        15923         2795  
# 
# kin_nodup<-kin[which(kin$FID1 %in% dat$V1 & kin$FID2 %in% dat$V1),]
# 
# table(cut(kin_nodup$KINSHIP,breaks=c(1,0.45,0.177)))
# # (0.177,0.45]     (0.45,1] 
# # 16135            0 


