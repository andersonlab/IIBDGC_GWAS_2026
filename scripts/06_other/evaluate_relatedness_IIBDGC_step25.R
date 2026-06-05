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
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy \
--allow-no-sex \
--extract ${path}pre_imputation/QC/relatedness/list_common_variants_among_cohorts_iibdgc_step25 \
--make-bed \
--out ${path}pre_imputation/QC/relatedness/${i}_subset_step25
done

# estimate missingness for later step
for i in ${studies[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy \
--allow-no-sex \
--missing \
--out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy 
done


ls -la ${path}pre_imputation/QC/relatedness/*_subset_step25.fam | wc -l 
# 58 OK

for i in ${studies[@]}
do wc -l ${path}pre_imputation/QC/relatedness/${i}_subset_step25.fam
done

############
# 22312 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/all_hce_subset_step25.fam
# 1157 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_old_gwas_subset_step25.fam
# 1291 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/australia_omniexome_subset_step25.fam
# 4655 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/gwas1_subset_step25.fam
# 5130 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/gwas2_subset_step25.fam
# 2709 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/pittsburgh_gsa_subset_step25.fam
# 3408 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/spain_gsa_subset_step25.fam
# 940 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/italy_gsa_subset_step25.fam
# 13830 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/kiel_austria_sibdcs_gsa_subset_step25.fam
# 4514 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/netherlands_gsa_subset_step25.fam
# 261 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/slovenia_gsa_subset_step25.fam
# 1374 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/sweden_gsa_subset_step25.fam
# 826 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_broad_gsa_subset_step25.fam
# 7376 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_feinstein_gsa_subset_step25.fam
# 1481 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/basque_gsa_subset_step25.fam
# 2204 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/lithuania_gsa_subset_step25.fam
# 1495 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_louis_gsa_subset_step25.fam
# 4 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_franchimont_gsa_subset_step25.fam
# 128 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_vermeire_gsa_subset_step25.fam
# 426 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/prism_nfe_gsa_subset_step25.fam
# 491 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/prism_nfe_gwas_subset_step25.fam
# 442 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/finland_illugwas_subset_step25.fam
# 2596 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/german_affy6_old_gwas_subset_step25.fam
# 544 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/norway_affy6_old_gwas_subset_step25.fam
# 1056 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_inf1_old_gwas_subset_step25.fam
# 158 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_inf2_old_gwas_subset_step25.fam
# 443 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_370k_old_gwas_subset_step25.fam
# 519 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_610k_old_gwas_subset_step25.fam
# 734 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_omni_old_gwas_subset_step25.fam
# 1242 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/swedish_uc_old_gwas_subset_step25.fam
# 5 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/mccauley_gsa_subset_step25.fam
############# 2166 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/ccfa_gsa_subset_step25.fam ##### exclude!!!!!
# 46 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_gsa_subset_step25.fam
# 479 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/bernstein_gsa_subset_step25.fam
# 68 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/farkkila_gsa_subset_step25.fam
# 2593 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/franchimont_gsa_subset_step25.fam
# 767 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/franke_gsa_subset_step25.fam
# 91 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/helmsley_prism_gsa_subset_step25.fam
# 961 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/helmsley_xavier_prism_gsa_subset_step25.fam
# 413 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/hyams_protect_gsa_subset_step25.fam
# 2831 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/lewis_sparc_gsa_subset_step25.fam
# 1603 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/mccauley_new_gsa_subset_step25.fam
# 5910 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/mcgovern_gsa_subset_step25.fam
# 1119 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/moayyedi_imagine_gsa_subset_step25.fam
# 856 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/newberry_share_gsa_subset_step25.fam
# 1438 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_cho_gsa_subset_step25.fam
# 1263 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_duerr_gsa_subset_step25.fam
# 634 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_rioux_gsa_subset_step25.fam
# 1913 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_silverberg_gsa_subset_step25.fam
# 868 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/palotie_hus_gsa_subset_step25.fam
# 631 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/pekow_share_gsa_subset_step25.fam
# 170 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/rioux_igenomed_gsa_subset_step25.fam
# 1405 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/sands_msccr_gsa_subset_step25.fam
# 1464 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/stampfer_gsa_subset_step25.fam
# 4641 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/vermeire_gsa_subset_step25.fam
# 13 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/weersma_gsa_subset_step25.fam
# 283 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/xavier_prism_gsa_subset_step25.fam
# 661 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/xavier_share_gsa_subset_step25.fam
##############

### ccfa study completely replaced by lewis:
rm /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/ccfa_gsa_subset_step25.*

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do wc -l ${path}pre_imputation/QC/relatedness/${i}_subset_step25.bim
done
wc -l ${path}pre_imputation/QC/relatedness/*_subset.bim
# 3148  OK

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

cohorts<-c("all_hce","niddk_old_gwas","australia_omniexome","gwas1","gwas2","pittsburgh_gsa","spain_gsa","italy_gsa",
           "kiel_austria_sibdcs_gsa","netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa",
           "niddk_feinstein_gsa","basque_gsa","lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa",
           "belgium_vermeire_gsa","prism_nfe_gsa","prism_nfe_gwas","finland_illugwas","german_affy6_old_gwas",
           "norway_affy6_old_gwas","belgium_inf1_old_gwas","belgium_inf2_old_gwas",
           "cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas",
           "swedish_uc_old_gwas",
           "mccauley_gsa",
           # "ccfa_gsa",
           "cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

length(cohorts)
# [1] 57 OK - no ccfa

#####

cohorts<-cohorts[-1]

dat<-matrix(ncol=3,nrow=length(cohorts))
dat<-as.data.frame(dat)

for (i in 1:nrow(dat)){
  dat[i,1]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_step25.bed",sep="")
  dat[i,2]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_step25.bim",sep="")
  dat[i,3]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_step25.fam",sep="")
}

write.table(dat,paste(path,"pre_imputation/QC/relatedness/list_cohorts_tomerge_step25.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#############################

# hpc-server

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/relatedness/all_hce_subset_step25 \
--allow-no-sex \
--merge-list ${path}pre_imputation/QC/relatedness/list_cohorts_tomerge_step25.txt \
--make-bed --out ${path}pre_imputation/QC/relatedness/iibdgc_merged_step25

# 3148 variants and 116872 people pass filters and QC.
# Among remaining phenotypes, 74403 are cases and 42469 are controls.

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

# (3148*116819/4)/1E+6
# [1] 91

path=/path/to/ibdgwas/IIBDGC/
MEM=3500

bsub -J"kg" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal -n 2 \
-e ${path}pre_imputation/QC/relatedness/logs/stderr_king_iibdgc_plink2_noparallel_step25 \
-o ${path}pre_imputation/QC/relatedness/logs/stdout_king_iibdgc_plink2_noparallel_step25 \
"/path/to/software/./plink2 \
--bfile ${path}pre_imputation/QC/relatedness/iibdgc_merged_step25 \
--make-king-table --king-table-filter 0.177 \
--memory $MEM \
--out ${path}pre_imputation/QC/relatedness/iibdgc_merged_step25_plink"
# Job <984725> is submitted to queue <normal>.


################################

##### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("all_hce","niddk_old_gwas","australia_omniexome","gwas1","gwas2","pittsburgh_gsa","spain_gsa","italy_gsa",
           "kiel_austria_sibdcs_gsa","netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa",
           "niddk_feinstein_gsa","basque_gsa","lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa",
           "belgium_vermeire_gsa","prism_nfe_gsa","prism_nfe_gwas","finland_illugwas","german_affy6_old_gwas",
           "norway_affy6_old_gwas","belgium_inf1_old_gwas","belgium_inf2_old_gwas",
           "cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas",
           "swedish_uc_old_gwas",
           "mccauley_gsa",
           # "ccfa_gsa",
           "cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

length(cohorts)

for (i in 1:length(cohorts)) {
  
  tmp<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],
                        "_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy.fam",sep=""),head=F)

  
  sample_miss<-read.table(paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],
                                "_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe_het_sexcheck_step24_step25_nodupsample_perstudy.imiss",sep=""),head=T)
  
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
# all_hce       australia_omniexome                basque_gsa 
# 22312                      1291                      1481 
# belgium_franchimont_gsa     belgium_inf1_old_gwas     belgium_inf2_old_gwas 
# 4                      1056                       158 
# belgium_louis_gsa      belgium_vermeire_gsa             bernstein_gsa 
# 1495                       128                       479 
# cedars_370k_old_gwas      cedars_610k_old_gwas                cedars_gsa 
# 443                       519                        46 
# cedars_omni_old_gwas              farkkila_gsa          finland_illugwas 
# 734                        68                       442 
# franchimont_gsa                franke_gsa     german_affy6_old_gwas 
# 2593                       767                      2596 
# gwas1                     gwas2        helmsley_prism_gsa 
# 4655                      5130                        91 
# helmsley_xavier_prism_gsa         hyams_protect_gsa                 italy_gsa 
# 961                       413                       940 
# kiel_austria_sibdcs_gsa           lewis_sparc_gsa             lithuania_gsa 
# 13830                      2831                      2204 
# mccauley_gsa          mccauley_new_gsa              mcgovern_gsa 
# 5                      1603                      5910 
# moayyedi_imagine_gsa           netherlands_gsa        newberry_share_gsa 
# 1119                      4514                       856 
# niddk_broad_gsa             niddk_cho_gsa           niddk_duerr_gsa 
# 826                      1438                      1263 
# niddk_feinstein_gsa            niddk_old_gwas           niddk_rioux_gsa 
# 7376                      1157                       634 
# niddk_silverberg_gsa     norway_affy6_old_gwas           palotie_hus_gsa 
# 1913                       544                       868 
# pekow_share_gsa            pittsburgh_gsa             prism_nfe_gsa 
# 631                      2709                       426 
# prism_nfe_gwas        rioux_igenomed_gsa           sands_msccr_gsa 
# 491                       170                      1405 
# slovenia_gsa                 spain_gsa              stampfer_gsa 
# 261                      3408                      1464 
# sweden_gsa       swedish_uc_old_gwas              vermeire_gsa 
# 1374                      1242                      4641 
# weersma_gsa          xavier_prism_gsa          xavier_share_gsa 
# 13                       283                       661

# double check all samples included in exercise
famall<-read.table(paste(path,"pre_imputation/QC/relatedness/iibdgc_merged_step25.fam",sep=""),head=F)
updated_fam<-fread(paste(path,"pheno/gwas-mega-2-core_updated_sex_gender_phenotype_Dec22.txt.gz",sep=""),head=T)
famall<-merge(famall[,c(1:4)],updated_fam[,c(1,9,8)],by="V2",all.x=T,sort=F)
rm(updated_fam)


dim(famall)
# [1] 116872     6

dim(dat)
# [1] 116872     8

dim(dat[which(dat$V1 %in% famall$V1),])
# [1] 116872     8



colnames(dat)[6]<-"pheno"

aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary)
#                      Group.1   x.1   x.2
# 1                    all_hce 10422 11890
# 2        australia_omniexome   608   683
# 3                 basque_gsa   960   521
# 4    belgium_franchimont_gsa     1     3
# 5      belgium_inf1_old_gwas   833   223
# 6      belgium_inf2_old_gwas    76    82
# 7          belgium_louis_gsa   596   899
# 8       belgium_vermeire_gsa    14   114
# 9              bernstein_gsa    10   469
# 10      cedars_370k_old_gwas     0   443
# 11      cedars_610k_old_gwas     0   519
# 12                cedars_gsa     9    37
# 13      cedars_omni_old_gwas     2   732
# 14              farkkila_gsa     0    68
# 15          finland_illugwas     0   442
# 16           franchimont_gsa  1368  1225
# 17                franke_gsa   413   354
# 18     german_affy6_old_gwas  1642   954
# 19                     gwas1  2919  1736
# 20                     gwas2  2776  2354
# 21        helmsley_prism_gsa    36    55
# 22 helmsley_xavier_prism_gsa   158   803
# 23         hyams_protect_gsa     0   413
# 24                 italy_gsa   368   572
# 25   kiel_austria_sibdcs_gsa  4439  9391
# 26           lewis_sparc_gsa     0  2831
# 27             lithuania_gsa  1138  1066
# 28              mccauley_gsa     0     5
# 29          mccauley_new_gsa   235  1368
# 30              mcgovern_gsa  1168  4742
# 31      moayyedi_imagine_gsa   251   868
# 32           netherlands_gsa   548  3966
# 33        newberry_share_gsa     0   856
# 34           niddk_broad_gsa    30   796
# 35             niddk_cho_gsa   443   995
# 36           niddk_duerr_gsa   437   826
# 37       niddk_feinstein_gsa  2333  5043
# 38            niddk_old_gwas   659   498
# 39           niddk_rioux_gsa   313   321
# 40      niddk_silverberg_gsa   158  1755
# 41     norway_affy6_old_gwas   279   265
# 42           palotie_hus_gsa     0   868
# 43           pekow_share_gsa     0   631
# 44            pittsburgh_gsa  1457  1252
# 45             prism_nfe_gsa    31   395
# 46            prism_nfe_gwas   217   274
# 47        rioux_igenomed_gsa     0   170
# 48           sands_msccr_gsa   302  1103
# 49              slovenia_gsa   175    86
# 50                 spain_gsa  1466  1942
# 51              stampfer_gsa  1028   436
# 52                sweden_gsa   979   395
# 53       swedish_uc_old_gwas   327   915
# 54              vermeire_gsa   807  3834
# 55               weersma_gsa     8     5
# 56          xavier_prism_gsa    30   253
# 57          xavier_share_gsa     0   661

colnames(dat)[5]<-"sex"

aggregate(as.factor(dat$sex), by=list(dat$cohort), FUN=summary)
#                      Group.1   x.0   x.1   x.2
# 1                    all_hce    27 10326 11959
# 2        australia_omniexome     1   598   692
# 3                 basque_gsa     6   905   570
# 4    belgium_franchimont_gsa     0     1     3
# 5      belgium_inf1_old_gwas    11   839   206
# 6      belgium_inf2_old_gwas     1    91    66
# 7          belgium_louis_gsa     7   648   840
# 8       belgium_vermeire_gsa     3    60    65
# 9              bernstein_gsa     1   188   290
# 10      cedars_370k_old_gwas     1   234   208
# 11      cedars_610k_old_gwas     2   268   249
# 12                cedars_gsa     1    28    17
# 13      cedars_omni_old_gwas     0   379   355
# 14              farkkila_gsa    17    25    26
# 15          finland_illugwas     0   294   148
# 16           franchimont_gsa     0  1202  1391
# 17                franke_gsa    16   448   303
# 18     german_affy6_old_gwas     2  1247  1347
# 19                     gwas1     2  2110  2543
# 20                     gwas2    51  2551  2528
# 21        helmsley_prism_gsa     1    51    39
# 22 helmsley_xavier_prism_gsa     4   426   531
# 23         hyams_protect_gsa     0   210   203
# 24                 italy_gsa     3   509   428
# 25   kiel_austria_sibdcs_gsa   203  6755  6872
# 26           lewis_sparc_gsa     6  1274  1551
# 27             lithuania_gsa     0  1170  1034
# 28              mccauley_gsa     0     1     4
# 29          mccauley_new_gsa    11   823   769
# 30              mcgovern_gsa     7  2899  3004
# 31      moayyedi_imagine_gsa     5   486   628
# 32           netherlands_gsa    10  1876  2628
# 33        newberry_share_gsa     5   364   487
# 34           niddk_broad_gsa     0   465   361
# 35             niddk_cho_gsa     5   737   696
# 36           niddk_duerr_gsa     9   624   630
# 37       niddk_feinstein_gsa     7  3699  3670
# 38            niddk_old_gwas     0   582   575
# 39           niddk_rioux_gsa     2   290   342
# 40      niddk_silverberg_gsa    13   985   915
# 41     norway_affy6_old_gwas     1   301   242
# 42           palotie_hus_gsa     8   410   450
# 43           pekow_share_gsa     2   314   315
# 44            pittsburgh_gsa     7  1179  1523
# 45             prism_nfe_gsa     3   195   228
# 46            prism_nfe_gwas     0   230   261
# 47        rioux_igenomed_gsa     0    81    89
# 48           sands_msccr_gsa   702   403   300
# 49              slovenia_gsa     5    55   201
# 50                 spain_gsa    56  1861  1491
# 51              stampfer_gsa     2   235  1227
# 52                sweden_gsa     0   795   579
# 53       swedish_uc_old_gwas     0   684   558
# 54              vermeire_gsa    20  2242  2379
# 55               weersma_gsa     0     2    11
# 56          xavier_prism_gsa     1   152   130
# 57          xavier_share_gsa     6   327   328


table(dat$V1==dat$V2)
# TRUE 
# 116872


kin<-read.table(paste(path,"pre_imputation/QC/relatedness/iibdgc_merged_step25_plink.kin0",sep=""),head=F)
colnames(kin)<-c("FID1","IID1","FID2","IID2","NSNP","HETHET","IBS0","KINSHIP")

table(cut(kin$KINSHIP,breaks=c(1,0.354,0.177)))
# (0.177,0.354]     (0.354,1] 
#          4535          4104


pdf(paste(path,"pre_imputation/QC/relatedness/histogram_iibdgc_merged_step25_plink.pdf",sep=""),height = 7,width = 14)
h <- hist(kin$KINSHIP, breaks = "FD", plot = FALSE)
ggplot(kin, aes(x=KINSHIP)) + geom_histogram(breaks = h$breaks) + geom_vline(xintercept=0.45, linetype="dashed", color = "red")
dev.off()

dup<-kin[which(kin$KINSHIP>0.45),c("FID1","FID2","KINSHIP")]
summary(dup$KINSHIP)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4836  0.5000  0.5000  0.4999  0.5000  0.5000


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
# [1] 4103    2

xx<-rbind(xx,yy)

x<-table(xx$cohort_FID1,xx$cohort_FID2)

write.table(x,paste(path,"pre_imputation/QC/relatedness/iibdgc_dup_summary_table_step25.tsv",sep=""),col.names=T,row.names=T,sep="\t",quote=F)

table(dup$sex_FID1,dup$sex_FID2)
#      0    1    2
# 0    3  113  109
# 1    2 1940    0
# 2    8    0 1928


dup<-dup[,c("FID1","FID2","KINSHIP","cohort_FID1","pheno_FID1","F_MISS_FID1","cohort_FID2","pheno_FID2","F_MISS_FID2")]

table(dup$pheno_FID1,dup$pheno_FID2)
#      1    2
# 1  440   37
# 2    4 3622

dup_ids<-c(as.character(dup$FID1),as.character(dup$FID2))
length(dup_ids)
#[1] 8206
table(length(dup_ids)==(nrow(dup)*2))
# TRUE 
# 1 

dup_ids<-dup_ids[!duplicated(dup_ids)]
length(dup_ids)
#[1] 7966

dim(dup)
# [1] 4103    9

# CHECK ALSO THAT BOTH SAMPLES HAVE SAME PHENOTYPE AND GENDER, 
# - keep all new gwas samples

# - for ukbb-ukbb duplicated pairs:
#    - keep sample with pheno data when both samples have the same phenotype, otherwise exclude both

new_gwas<-c("spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa","all_hce"
            ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa"
            ,"basque_gsa","prism_nfe_gsa","lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
            ,"mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
            "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
            "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
            "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
            "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
            "xavier_share_gsa")

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
      
      # B.1: rest of samples, if samples in new and old arrays, keep only samples in new arrays:
      
      # B.1: if one sample from niddk_old - exclude this one first, but check if any other from broad or feinstein
      
      if (any(data$cohort %in% c("mcgovern_gsa"))) {
        
        keep_sample<-data$V1[which(data$cohort=="mcgovern_gsa")]
        
        # if more than one new pair in new gwas, keep sample with largest call rate:
        if (length(data$V1[which(data$cohort=="mcgovern_gsa")])>1) {
          keep_sample<-data$V1[which(data$cohort=="mcgovern_gsa"),]
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
        
      # if more than one new pair in new gwas, keep sample with largest call rate, make exception for CCFA and McCauley:
        
        if (length(keep_sample)>1) {
          keep_sample<-data[which(!data$cohort %in% old_gwas),]
          keep_sample<-keep_sample$V1[which(keep_sample$F_MISS==min(keep_sample$F_MISS))][1]# to deal with >1 with same F_miss
         }
        
        if (!exists("data_remove")) {
         data_remove<-data[which(!data$V1 %in% keep_sample),]
         } else {
          data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
        }
  
        # B.2: old_gwas duplicated pairs only
      } else if (all(data$cohort %in% old_gwas)) {
        
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
        # B.2.2: if all in good OR all in bad array, keep sample with lowest missingness
        } else {
          
          keep_sample<-data$V1[which(data$F_MISS==min(data$F_MISS))][1]# to deal with same F_miss
          
          if(!exists("data_remove")) {
            data_remove<-data[which(!data$V1 %in% keep_sample),]
          } else {
               data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
          }
        }
  
        # B.3: new_gwas duplicated pairs only, keep sample with largest call rate: 
      
        } else if ( all( any(data$cohort=="all_hce" | data$cohort=="belgium_vermeire_gsa") & any(data$cohort %in% c("kiel_austria_sibdcs_gsa")) )  ) {
          
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
        
        # B.3.3: Other keep sample with largest call rate:
        
        } else {
          
          keep_sample<-data$V1[which(data$F_MISS==min(data$F_MISS))]
          
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
#[1] 8248  4

data_remove<-data_remove[!duplicated(data_remove$V1),]
dim(data_remove)
# [1] 4046    5
# double check which sample per pair removed in Kiel_Uk pairs:

sort(table(data_remove$cohort),decreasing=T)
# niddk_feinstein_gsa            pittsburgh_gsa           niddk_broad_gsa 
# 1301                       652                       641 
# cedars_omni_old_gwas           sands_msccr_gsa      cedars_610k_old_gwas 
# 266                       243                       196 
# niddk_cho_gsa            niddk_old_gwas           lewis_sparc_gsa 
# 149                       140                       106 
# cedars_370k_old_gwas           niddk_duerr_gsa   kiel_austria_sibdcs_gsa 
# 61                        47                        35 
# mcgovern_gsa      moayyedi_imagine_gsa          mccauley_new_gsa 
# 28                        26                        22 
# finland_illugwas                cedars_gsa           pekow_share_gsa 
# 19                        17                        16 
# prism_nfe_gwas          xavier_share_gsa             bernstein_gsa 
# 13                        13                        12 
# newberry_share_gsa        rioux_igenomed_gsa      niddk_silverberg_gsa 
# 10                         9                         8 
# helmsley_xavier_prism_gsa             prism_nfe_gsa              vermeire_gsa 
# 7                         7                         5 
# netherlands_gsa         hyams_protect_gsa           niddk_rioux_gsa 
# 4                         3                         3 
# german_affy6_old_gwas           palotie_hus_gsa 
# 1                         1 


head(dup[which(!dup$FID1 %in% data_remove$V1 & !dup$FID2 %in% data_remove$V1),])
# 0 OK



#### create a file to report inconsitend data:
dim(data_inconsist)
# [1] 152 5
data_inconsist<-data_inconsist[!duplicated(data_inconsist$V1),]
dim(data_inconsist)
# [1] 76  6
data_inconsist
#

data_inconsist$cohort[which(data_inconsist$cohort=="pittsburgh_gsa")]<-"pittsburgh"


### save files:

data_remove_tmp<-data_remove[,c(1,1)]
colnames(data_remove_tmp)<-c("FID","IID")

write.table(data_remove,paste(path,"pre_imputation/QC/relatedness/all_iibdgc_data_remove_step25",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
write.table(data_remove_tmp,paste(path,"pre_imputation/QC/relatedness/list_all_iibdgc_data_remove_step25",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
write.table(data_inconsist,paste(path,"pre_imputation/QC/relatedness/all_iibdgc_data_inconsist_step25",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
write.table(aggregate(as.factor(dat$pheno), by=list(dat$cohort), FUN=summary),paste(path,"pre_imputation/QC/relatedness/all_iibdgc_case_control_pre_per_study_relatedness_step25",sep=""),col.names=T,row.names=F,sep="\t",quote=F)


##############################################################################################################################################
##############################################################################################################################################

## get summary by cohort and per ancestry:

##### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("all_hce","niddk_old_gwas","australia_omniexome","gwas1","gwas2","pittsburgh_gsa","spain_gsa","italy_gsa",
           "kiel_austria_sibdcs_gsa","netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa",
           "niddk_feinstein_gsa","basque_gsa","lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa",
           "belgium_vermeire_gsa","prism_nfe_gsa","prism_nfe_gwas","finland_illugwas","german_affy6_old_gwas",
           "norway_affy6_old_gwas","belgium_inf1_old_gwas","belgium_inf2_old_gwas",
           "cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas",
           "swedish_uc_old_gwas",
           "mccauley_gsa",
           # "ccfa_gsa",
           "cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")


data_remove,paste(path,"pre_imputation/QC/relatedness/all_iibdgc_data_remove_step25",sep="")



##############################################################################################################################################
##############################################################################################################################################

# combine all phenotype inconsistencies into one file:

#### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

data_inconsist<-fread(paste(path,"pre_imputation/QC/relatedness/all_iibdgc_data_inconsist_step25",sep=""),head=T,sep="\t")
data_inconsist<-as.data.frame(data_inconsist)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("all_hce","niddk_old_gwas","australia_omniexome","gwas1","gwas2","pittsburgh_gsa","spain_gsa","italy_gsa",
           "kiel_austria_sibdcs_gsa","netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa",
           "niddk_feinstein_gsa","basque_gsa","lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa",
           "belgium_vermeire_gsa","prism_nfe_gsa","prism_nfe_gwas","finland_illugwas","german_affy6_old_gwas",
           "norway_affy6_old_gwas","belgium_inf1_old_gwas","belgium_inf2_old_gwas",
           "cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas",
           "swedish_uc_old_gwas",
           "mccauley_gsa","ccfa_gsa",
           "cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")


# intra-study:
for (j in 1:length(cohorts)) {
  file_inconsist<-paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_list_duplicated_samples_inconsistent_phenotype_step25",sep="")
  
  if(file.exists(file_inconsist)) {
    tmp<-read.table(file_inconsist,head=T)
    tmp$cohort<-cohorts[j]
    tmp<-tmp[,colnames(data_inconsist)]

    tmp$group<-max(data_inconsist$group)+(tmp$group)
    data_inconsist<-rbind(data_inconsist,tmp)
    rm(tmp)
  }
  rm(file_inconsist)
}



centers<-c("ukibdgc","german","niddk","belgium","prism","sweden","cedars","netherlands","mccauley","ccfa")


for (j in 1:length(centers)){
  
  file_inconsist<-paste(path,"pre_imputation/QC/relatedness/",centers[j],"_data_inconsist_step25",sep="")

  if(file.exists(file_inconsist)) {
    tmp<-read.table(file_inconsist,head=T)
    tmp<-tmp[,colnames(data_inconsist)]

    tmp$group<-max(data_inconsist$group)+(tmp$group)
    data_inconsist<-rbind(data_inconsist,tmp)
    rm(tmp)
  }
  rm(file_inconsist)
}

length(table(data_inconsist$group))
# [1] 1484

groups<-names(table(data_inconsist$group))

data_inconsist$groups_2<-NA
for (i in 1:length(groups)) {
  data_inconsist$groups_2[which(data_inconsist$group==groups[i])]<-i
  
}

data_inconsist<-data_inconsist[,c(1:5,7)]
colnames(data_inconsist)[6]<-"group"

write.table(data_inconsist,paste(path,"pre_imputation/QC/relatedness/list_iibdgc_data_inconsistent_phenotype_Dec22.tsv.gz",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")
  
  
  
