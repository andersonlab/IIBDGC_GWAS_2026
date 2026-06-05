# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# UK-IIBGBD:

# See how variants were defined in:

# create_list_common_variants_among_cohorts.R
# after talking with Carl, generate PC with 1000GP samples and then project the IIBDGC

##################################################
# 0.- CREATE FILES INCLUDING DUPLICATED SAMPLES: #
##################################################

path=/path/to/ibdgwas/IIBDGC/
  
studies=(australia_omniexome gwas1 gwas2 all_hce pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas chop_old_gwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas niddk_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa)

for i in ${studies[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy \
--allow-no-sex \
--extract ${path}pre_imputation/QC/relatedness/list_common_variants_among_cohorts_iibdgc \
--make-bed \
--out ${path}pre_imputation/QC/relatedness/${i}_subset_withDuplicates
done


#######################
# 0.1 MERGE ALL FILES

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
# [1] 34 OK
  

  
#####
  
cohorts<-cohorts[-1]
  
dat<-matrix(ncol=3,nrow=length(cohorts))
dat<-as.data.frame(dat)
  
for (i in 1:nrow(dat)){
  dat[i,1]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_withDuplicates.bed",sep="")
  dat[i,2]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_withDuplicates.bim",sep="")
  dat[i,3]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_withDuplicates.fam",sep="")
}
  
write.table(dat,paste(path,"pre_imputation/QC/relatedness/list_cohorts_tomerge_withDuplicates.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  
#############################
  
# hpc-server
  
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/relatedness/australia_omniexome_subset_withDuplicates \
--allow-no-sex \
--merge-list ${path}pre_imputation/QC/relatedness/list_cohorts_tomerge_withDuplicates.txt \
--make-bed --out ${path}pre_imputation/QC/relatedness/iibdgc_merged_withDuplicates
# 3111 variants and 117923 people pass filters and QC.
# Among remaining phenotypes, 69577 are cases and 48346 are controls.

wc -l ${path}pre_imputation/QC/relatedness/*_subset_withDuplicates.fam
# 23890 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/all_hce_subset_withDuplicates.fam
# 1306 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/australia_omniexome_subset_withDuplicates.fam
# 1508 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/basque_gsa_subset_withDuplicates.fam
# 1501 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_franchimont_gsa_subset_withDuplicates.fam
# 1417 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_inf1_old_gwas_subset_withDuplicates.fam
# 272 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_inf2_old_gwas_subset_withDuplicates.fam
# 1502 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_louis_gsa_subset_withDuplicates.fam
# 3993 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/belgium_vermeire_gsa_subset_withDuplicates.fam
# 2177 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/ccfa_gsa_subset_withDuplicates.fam
# 605 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_370k_old_gwas_subset_withDuplicates.fam
# 889 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_610k_old_gwas_subset_withDuplicates.fam
# 3078 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_gsa_subset_withDuplicates.fam
# 1221 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/cedars_omni_old_gwas_subset_withDuplicates.fam
# 8610 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/chop_old_gwas_subset_withDuplicates.fam
# 446 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/finland_illugwas_subset_withDuplicates.fam
# 2823 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/german_affy6_old_gwas_subset_withDuplicates.fam
# 4680 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/gwas1_subset_withDuplicates.fam
# 7778 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/gwas2_subset_withDuplicates.fam
# 953 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/italy_gsa_subset_withDuplicates.fam
# 14403 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/kiel_austria_sibdcs_gsa_subset_withDuplicates.fam
# 2231 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/lithuania_gsa_subset_withDuplicates.fam
# 782 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/mccauley_gsa_subset_withDuplicates.fam
# 4581 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/netherlands_gsa_subset_withDuplicates.fam
# 5455 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_broad_gsa_subset_withDuplicates.fam
# 8177 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_feinstein_gsa_subset_withDuplicates.fam
# 2760 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/niddk_old_gwas_subset_withDuplicates.fam
# 550 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/norway_affy6_old_gwas_subset_withDuplicates.fam
# 2725 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/pittsburgh_gsa_subset_withDuplicates.fam
# 466 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/prism_nfe_gsa_subset_withDuplicates.fam
# 818 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/prism_nfe_gwas_subset_withDuplicates.fam
# 264 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/slovenia_gsa_subset_withDuplicates.fam
# 3443 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/spain_gsa_subset_withDuplicates.fam
# 1355 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/sweden_gsa_subset_withDuplicates.fam
# 1264 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/relatedness/swedish_uc_old_gwas_subset_withDuplicates.fam
# 117923 total

###########################################################################################################################################################################
  
#########################################
# 1.- REDEFINE LIST OF COMMON VARIANTS: #
#########################################

# SEE evaluate_relatedness_IIBDGC.R to see how ${path}pre_imputation/QC/relatedness/iibdgc_merged was created

path=/path/to/ibdgwas/IIBDGC/

  
# 1.1 EXCLUDE VARIANTS IN HIGH LD REGIONS:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/relatedness/iibdgc_merged_nodup \
--allow-no-sex \
--exclude /nfs/team152/carl/protocols/QC-paper/high-LD-regions.txt --range \
--make-bed --out ${path}pre_imputation/QC/pca_1000gp/iibdgc_merged_noHighLD
# 3077 variants and 106425 people pass filters and QC.
# Among remaining phenotypes, 61780 are cases and 44645 are controls.



# 1.2 PRUNE VARIANTS:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/pca_1000gp/iibdgc_merged_noHighLD \
--allow-no-sex \
--indep-pairwise 50 5 0.2 \
--make-bed --out ${path}pre_imputation/QC/pca_1000gp/iibdgc_merged_noHighLD

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/pca_1000gp/iibdgc_merged_noHighLD \
--extract ${path}pre_imputation/QC/pca_1000gp/iibdgc_merged_noHighLD.prune.in \
--make-bed --out ${path}pre_imputation/QC/pca_1000gp/iibdgc_merged_noHighLD_pruned
# 2989 variants and 106425 people pass filters and QC.
# Among remaining phenotypes, 61780 are cases and 44645 are controls.

# 1.3 EXCLUDE ANY ASSOCIATED VARIANT:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/pca_1000gp/iibdgc_merged_noHighLD_pruned \
--assoc --out ${path}pre_imputation/QC/pca_1000gp/iibdgc_merged_noHighLD_pruned

#### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(viridis)

path<-"/path/to/ibdgwas/IIBDGC/"

as<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_merged_noHighLD_pruned.assoc",sep=""),head=T)

summary(as$P)
# 
nrow(as[which(as$P<=1E-4),])
# [1] 666

write.table(as[which(as$P<=1E-4),"SNP",drop=F],paste(path,"pre_imputation/QC/pca_1000gp/list_associated_variants_toexclude",sep=""),col.names=F
            ,row.names=F,quote=F,sep="\t")

nrow(as[which(as$P>1E-4),])
# [1] 2323
write.table(as[which(as$P>1E-4),"SNP",drop=F],paste(path,"pre_imputation/QC/pca_1000gp/list_noassociated_variants_tokeep",sep=""),col.names=F
            ,row.names=F,quote=F,sep="\t")

####

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/pca_1000gp/iibdgc_merged_noHighLD_pruned \
--exclude ${path}pre_imputation/QC/pca_1000gp/list_associated_variants_toexclude \
--make-bed --out ${path}pre_imputation/QC/pca_1000gp/iibdgc_merged_noHighLD_pruned_noassoclogP4
# 2323 variants and 106425 people pass filters and QC.
# Among remaining phenotypes, 61780 are cases and 44645 are controls.



###########################################################################################################################################################################

######################################
# 2.- SUBSET LIST OF COMMON VARIANTS #
######################################

#################################################
# 2.1- FROM IIBDGC FILES WITH DUP SAMPLES:

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/relatedness/iibdgc_merged_withDuplicates \
--extract ${path}pre_imputation/QC/pca_1000gp/list_noassociated_variants_tokeep \
--make-bed --out ${path}pre_imputation/QC/pca_1000gp/iibdgc_merged_withDuplicates_subset_iibdgc
# 2323 variants and 117923 people pass filters and QC.
# Among remaining phenotypes, 69577 are cases and 48346 are controls.

#################################################
# 2.1- FROM 1000GP FILES AND MERGE 

for chr in {1..22}; do /path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path}pre_imputation/QC/1000gp/1000GP_chr${chr}_b37.bed \
--bim ${path}pre_imputation/QC/1000gp/1000GP_chr${chr}_b37_edited.bim \
--fam ${path}pre_imputation/QC/1000gp/1000GP_chr${chr}_b37.fam \
--extract ${path}pre_imputation/QC/pca_1000gp/list_noassociated_variants_tokeep \
--make-bed --out ${path}pre_imputation/QC/1000gp/1000GP_chr${chr}_b37_subset_iibdgc;done

#### R

# no chrx or chry
dat<-matrix(ncol=3,nrow=21)
dat<-as.data.frame(dat)
for (i in 1:nrow(dat)){
  dat[i,1]<-paste(path,"pre_imputation/QC/1000gp/1000GP_chr",i+1,"_b37_subset_iibdgc.bed",sep="")
  dat[i,2]<-paste(path,"pre_imputation/QC/1000gp/1000GP_chr",i+1,"_b37_subset_iibdgc.bim",sep="")
  dat[i,3]<-paste(path,"pre_imputation/QC/1000gp/1000GP_chr",i+1,"_b37_subset_iibdgc.fam",sep="")
}

write.table(dat,paste(path,"pre_imputation/QC/1000gp/list_1000GP_files_subset_iibdgc_tomerge.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/1000gp/1000GP_chr1_b37_subset_iibdgc \
--merge-list ${path}pre_imputation/QC/1000gp/list_1000GP_files_subset_iibdgc_tomerge.txt \
--allow-no-sex \
--make-bed --out ${path}pre_imputation/QC/pca_1000gp/1000GP_all_b37_subset_iibdgc
# 2323 variants and 2504 people pass filters and QC.
# Note: No phenotypes present.

# remove intermediate files:
rm ${path}pre_imputation/QC/1000gp/1000GP_chr*_b37_subset*
  
###################################################################################################################################################
# # EXPLORE PROJECTIONS - pre plink2 instructions
# 
# # 
# # Yes, it's likely that something like PLINK 1.9's --pca-clusters/--pca-cluster-names will eventually make it into 2.0.  
# # With that said, PCA projection is actually already supported, the workflow is just a bit more convoluted for now.
# # 
# # Step 1: Export allele frequencies and PCA variant weights from your reference dataset.  E.g.
# # plink2 --bfile hapmap --freq --pca var-wts --out pca_hapmap
# # 
# # Step 2: Use --score to compute the necessary dot products with the variant weights.  E.g.
# # plink2 --bfile mydata --read-freq pca_hapmap.afreq --score pca_hapmap.eigenvec.var 2 3 header-read no-mean-imputation variance-normalize 
# # --score-col-nums 5-14 --out pca_proj_mydata
# # 
# # These PCs will be scaled a bit differently from pca_hapmap.eigenvec (you need to multiply or divide the PCs by sqrt(eigenvalue) to put 
# them on the same scale).  One way to make the PCs directly comparable is to also run step 2 on the original dataset.
# # 

# # EXPLORE PROJECTIONS - plink2 instructions
# 
# PCA projection with --score
# Since --score's new 'variance-standardize' modifier applies the same transformation to G as --pca does, --score can now 
# execute the vector-matrix multiply corresponding to PCA projection.
# 
# The following command exports PCs to project onto, along with the allele frequencies needed to calibrate the 'variance-standardize' operation:
# 
# plink2 --pfile ref_data \
#        --freq counts \
#        --pca allele-wts \
#        --out ref_pcs
# 
# You can then project onto those PCs with
# 
# plink2 --pfile new_data \
#        --read-freq ref_pcs.acount \
#        --score ref_pcs.eigenvec.allele 2 5 header-read no-mean-imputation \
#                variance-standardize \
#        --score-col-nums 6-15 \
#        --out new_projection
# 
# Note that these PCs will be scaled a bit differently from ref_data.eigenvec; you need to multiply or divide the PCs by a multiple of sqrt(eigenvalue)
# to put them on the same scale.


################################
# 3.- RUN PCA AND PROJECT DATA #
################################

#########################################################
## 3.1.- ESTIMATE PC, EXPORT PCS, AND ALLELE FREQUENCIES:

# The following command exports PCs to project onto, along with the allele frequencies needed to calibrate the 'variance-standardize' operation:

path=/path/to/ibdgwas/IIBDGC/
MEM=1000 # next time only this

bsub -J"bcftools" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal -n 10 \
-e ${path}pre_imputation/QC/pca_1000gp/logs/stderr_pca_1000gp \
-o ${path}pre_imputation/QC/pca_1000gp/logs/stdout_pca_1000gp \
"/path/to/software/./plink2  \
--bfile ${path}pre_imputation/QC/pca_1000gp/1000GP_all_b37_subset_iibdgc \
--freq counts \
--pca allele-wts \
--out ${path}pre_imputation/QC/pca_1000gp/1000GP_all_b37_subset_iibdgc_ref_pcs \
--threads 8 --allow-no-sex --memory $MEM"
# Job <993391> is submitted to queue <normal>.

# less ${path}pre_imputation/QC/pca_1000gp/logs/stdout_pca_1000gp
# 2504 samples (0 females, 0 males, 2504 ambiguous; 2504 founders) loaded from
# /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/1000GP_all_b37_su
# bset_iibdgc.fam.
# 2323 variants loaded from
# /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/1000GP_all_b37
# _subset_iibdgc.bim.
# Note: No phenotype data present.
# Calculating allele frequencies... done.
# --freq counts: 99%^M--freq counts: Allele counts (founders only) written to
# /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/1000GP_all_b37
# _subset_iibdgc_ref_pcs.acount
# .
# Constructing GRM: done.
# Extracting eigenvalues and eigenvectors... done.
# --pca: Allele weights written to
# /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/1000GP_all_b37
# _subset_iibdgc_ref_pcs.eigenvec.allele
# .
# --pca: Eigenvectors written to
# /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/1000GP_all_b37
# _subset_iibdgc_ref_pcs.eigenvec
# , and eigenvalues written to
# /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/1000GP_all_b37
# _subset_iibdgc_ref_pcs.eigenval


#########################################################
## 3.2 MERGE 1000GP AND IIBDGC (WITH DUPLICATED SAMPLES)

# plink1.9, 2.0 has not bmerge available yet

path=/path/to/ibdgwas/IIBDGC/
MEM=800 

bsub -J"bcftools" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path}pre_imputation/QC/pca_1000gp/logs/stderr_pca_1000gp_merge_iibdgc \
-o ${path}pre_imputation/QC/pca_1000gp/logs/stdout_pca_1000gp_merge_iibdgc \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/pca_1000gp/1000GP_all_b37_subset_iibdgc \
--bmerge ${path}pre_imputation/QC/pca_1000gp/iibdgc_merged_withDuplicates_subset_iibdgc.bed \
${path}pre_imputation/QC/pca_1000gp/iibdgc_merged_withDuplicates_subset_iibdgc.bim \
${path}pre_imputation/QC/pca_1000gp/iibdgc_merged_withDuplicates_subset_iibdgc.fam \
--allow-no-sex \
--make-bed --out ${path}pre_imputation/QC/pca_1000gp/1000GP_all_b37_plus_iibdgc_subset_iibdgc"
# Job <993399> is submitted to queue <normal>.

less ${path}pre_imputation/QC/pca_1000gp/logs/stdout_pca_1000gp_merge_iibdgc
# 2323 variants and 120427 people pass filters and QC.
# Among remaining phenotypes, 69577 are cases and 48346 are controls.  (2504
#                                                                       phenotypes are missing.)


#########################################################
## 3.3.- PROJECT ONTO THOSE PCS 1000GP AND IIBDGC

path=/path/to/ibdgwas/IIBDGC/
MEM=1000 # next time only this

bsub -J"bcftools" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal -n 10 \
-e ${path}pre_imputation/QC/pca_1000gp/logs/stderr_pca_1000gp_iibdgc_pca \
-o ${path}pre_imputation/QC/pca_1000gp/logs/stdout_pca_1000gp_iibdgc_pca \
"/path/to/software/./plink2 \
--bfile ${path}pre_imputation/QC/pca_1000gp/1000GP_all_b37_plus_iibdgc_subset_iibdgc \
--read-freq ${path}pre_imputation/QC/pca_1000gp/1000GP_all_b37_subset_iibdgc_ref_pcs.acount \
--score ${path}pre_imputation/QC/pca_1000gp/1000GP_all_b37_subset_iibdgc_ref_pcs.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
--score-col-nums 6-15 \
--out ${path}pre_imputation/QC/pca_1000gp/1000GP_all_b37_plus_iibdgc_subset_iibdgc_new_projection"
# Job <993401> is submitted to queue <normal>.

less ${path}pre_imputation/QC/pca_1000gp/logs/stdout_pca_1000gp_iibdgc_pca
# 120427 samples (59607 females, 57869 males, 2951 ambiguous; 120427 founders)
# --score: Results written to
# /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/1000GP_all_b37_plus_iibdgc_subset_iibdgc_new_projection.sscore



########## /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(viridis)

path<-"/path/to/ibdgwas/IIBDGC/"

# no german_illu
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


for (i in 1:length(cohorts)){
  a<-read.table(paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_withDuplicates.fam",sep=""),head=F)
  a$cohort<-cohorts[i]
  if(i==1){
    fam<-a
  }else{
    fam<-rbind(a,fam)
  }
}
dim(fam)
# [1] 117923      7

table(fam$cohort)

samples1000<-read.table("/path/to/project",head=T)
table(samples1000$pop,samples1000$super_pop)
#      AFR AMR EAS EUR SAS
# ACB  96   0   0   0   0
# ASW  61   0   0   0   0
# BEB   0   0   0   0  86
# CDX   0   0  93   0   0
# CEU   0   0   0  99   0
# CHB   0   0 103   0   0
# CHS   0   0 105   0   0
# CLM   0  94   0   0   0
# ESN  99   0   0   0   0
# FIN   0   0   0  99   0
# GBR   0   0   0  91   0
# GIH   0   0   0   0 103
# GWD 113   0   0   0   0
# IBS   0   0   0 107   0
# ITU   0   0   0   0 102
# JPT   0   0 104   0   0
# KHV   0   0  99   0   0
# LWK  99   0   0   0   0
# MSL  85   0   0   0   0
# MXL   0  64   0   0   0
# PEL   0  85   0   0   0
# PJL   0   0   0   0  96
# PUR   0 104   0   0   0
# STU   0   0   0   0 102
# TSI   0   0   0 107   0
# YRI 108   0   0   0   0

fam<-fam[,c(1,7)]
colnames(fam)[1]<-c("sample")
fam$super_pop<-"IIBDGC"
fam$pop<-NA

samples1000$cohort<-"1000GP"
samples1000<-samples1000[,c("sample","cohort","super_pop","pop")]
fam<-rbind(fam,samples1000)
dim(fam)
# [1] 120427      4


pca<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/1000GP_all_b37_plus_iibdgc_subset_iibdgc_new_projection.sscore",sep=""),head=F)
colnames(pca)<-c("FID","IID","PHENO1","ALLELE_CT","NAMED_ALLELE_DOSAGE_SUM","PC1_AVG","PC2_AVG","PC3_AVG","PC4_AVG","PC5_AVG","PC6_AVG"
                 ,"PC7_AVG","PC8_AVG","PC9_AVG","PC10_AVG")
dim(pca)
# [1] 120427     12


pca<-merge(pca,fam[,c("sample","cohort","super_pop","pop")],by.x="IID",by.y="sample",all.x=T)
dim(pca)
# [1] 120427     18

table(pca$super_pop)
# AFR    AMR    EAS    EUR IIBDGC    SAS 
# 661    347    504    503 117923    489 

# TOTAL VARIANCE EXPLAINED BY EACH PC:

eigenval<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/1000GP_all_b37_subset_iibdgc_ref_pcs.eigenval",sep=""),head=F)

eigenval$var_exp<-NA
for (i in 1:nrow(eigenval)){
  eigenval$var_exp[i]<-(eigenval$V1[i] / sum(eigenval$V1))*100
}

eigenval
# V1   var_exp
# 1  220.44900 56.149531
# 2   90.95560 23.166874
# 3   29.99140  7.638969
# 4   23.36370  5.950858
# 5    5.88810  1.499730
# 6    5.12991  1.306615
# 7    4.62689  1.178493
# 8    4.26146  1.085416
# 9    3.98774  1.015699
# 10   3.95678  1.007813

# PC1 EUR - non EUR
# PC2 AFR  - non AFR
# PC3 north EUR -South EUR
# PC4 EAS - SAS
# PC5 Finish - non-Finish
# PC6: AMR - non-AMR
# PC7 - IIBDGC outliers (from more than 1 cohort)
# PC8 transition EUR - EAS - SAS/AFR - AMR

cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"chop_old_gwas","german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas"
           ,"swedish_uc_old_gwas","mccauley_gsa","ccfa_gsa","cedars_gsa")

pca$study<-NA

pca$study[which(pca$cohort=="australia_omniexome")]<-"Australia"
pca$study[which(pca$cohort %in% c("gwas1","gwas2","all_hce"))]<-"UK"
pca$study[which(pca$cohort %in% c("spain_gsa","basque_gsa"))]<-"Spain"
pca$study[which(pca$cohort %in% c("kiel_austria_sibdcs_gsa","german_affy6_old_gwas"))]<-"Germany"                 
pca$study[which(pca$cohort %in% c("italy_gsa"))]<-"Italy"   
pca$study[which(pca$cohort %in% c("netherlands_gsa"))]<-"Netherlands"                     
pca$study[which(pca$cohort %in% c("slovenia_gsa"))]<-"Slovenia"                     
pca$study[which(pca$cohort %in% c("sweden_gsa","swedish_uc_old_gwas"))]<-"Sweden"                    
pca$study[which(pca$cohort %in% c("finland_illugwas"))]<-"Finland"
pca$study[which(pca$cohort %in% c("chop_old_gwas"))]<-"CHOP"
pca$study[which(pca$cohort %in% c("norway_affy6_old_gwas"))]<-"Norway"
pca$study[which(pca$cohort %in% c("belgium_inf1_old_gwas","belgium_inf2_old_gwas","belgium_louis_gsa","belgium_franchimont_gsa"
                                  ,"belgium_vermeire_gsa"))]<-"Belgium"
pca$study[which(pca$cohort %in% c("lithuania_gsa"))]<-"Lithuania"
pca$study[which(pca$cohort %in% c("prism_nfe_gwas","prism_nfe_gsa","cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","cedars_gsa",
                                  "niddk_broad_gsa","niddk_feinstein_gsa","niddk_old_gwas","pittsburgh_gsa","mccauley_gsa","ccfa_gsa"))]<-"USA"


table(pca$study,useNA="ifany")
# Australia     Belgium        CHOP     Finland     Germany       Italy 
# 1306        8685        8610         446       17226         953 
# Lithuania Netherlands      Norway    Slovenia       Spain      Sweden 
# 2231        4581         550         264        4951        2619 
# UK         USA        <NA> 
#   36348       29153        2504 


####################
# 4.- PLOT RESULTS #
####################

######################################################
# 4.1 PLOT 1000GP ALONE, COLOUR BY SUPER POPULATION:

# creat plots PC1 to PC6, 1000GP only
pca$pop<-as.factor(pca$pop)
pca$pop<-factor(pca$pop, levels=c("CEU","FIN","GBR","IBS","TSI",
                                  "CDX","CHB","CHS","JPT","KHV",
                                  "BEB","GIH","ITU","PJL","STU",
                                  "ACB","ASW","ESN","GWD","LWK","MSL","YRI",
                                  "CLM","MXL","PEL","PUR"))

# EUR = "CEU","FIN","GBR","IBS","TSI"
# EAS = "CDX","CHB","CHS","JPT","KHV"
# SAS =  "BEB","GIH","ITU","PJL","STU"
# AFR = "ACB","ASW","ESN","GWD","LWK","MSL","YRI"
# AMF = "CLM","MXL","PEL","PUR"

pca$super_pop<-as.factor(pca$super_pop)
pca$super_pop<-factor(pca$super_pop, levels=c("EUR","AFR","EAS","SAS","AMR","IIBDGC"))

colnames(pca)<-gsub("_AVG","",colnames(pca))

pna<-qplot()+theme(
  panel.background = element_rect(fill = "transparent") # bg of the panel
  , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
  , panel.grid.major = element_blank() # get rid of major grid
  , panel.grid.minor = element_blank() # get rid of minor grid
  , legend.background = element_rect(fill = "transparent") # get rid of legend bg
  , axis.title=element_blank()
  , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
)

p11<-qplot(pca$PC1[which(pca$cohort=="1000GP")],pca$PC1[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p12<-qplot(pca$PC2[which(pca$cohort=="1000GP")],pca$PC1[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p13<-qplot(pca$PC3[which(pca$cohort=="1000GP")],pca$PC1[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p14<-qplot(pca$PC4[which(pca$cohort=="1000GP")],pca$PC1[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p15<-qplot(pca$PC5[which(pca$cohort=="1000GP")],pca$PC1[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p16<-qplot(pca$PC6[which(pca$cohort=="1000GP")],pca$PC1[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")

p21<-qplot(pca$PC1[which(pca$cohort=="1000GP")],pca$PC2[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p22<-qplot(pca$PC2[which(pca$cohort=="1000GP")],pca$PC2[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p23<-qplot(pca$PC3[which(pca$cohort=="1000GP")],pca$PC2[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p24<-qplot(pca$PC4[which(pca$cohort=="1000GP")],pca$PC2[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p25<-qplot(pca$PC5[which(pca$cohort=="1000GP")],pca$PC2[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p26<-qplot(pca$PC6[which(pca$cohort=="1000GP")],pca$PC2[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")

p31<-qplot(pca$PC1[which(pca$cohort=="1000GP")],pca$PC3[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p32<-qplot(pca$PC2[which(pca$cohort=="1000GP")],pca$PC3[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p33<-qplot(pca$PC3[which(pca$cohort=="1000GP")],pca$PC3[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p34<-qplot(pca$PC4[which(pca$cohort=="1000GP")],pca$PC3[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p35<-qplot(pca$PC5[which(pca$cohort=="1000GP")],pca$PC3[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p36<-qplot(pca$PC6[which(pca$cohort=="1000GP")],pca$PC3[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")

p41<-qplot(pca$PC1[which(pca$cohort=="1000GP")],pca$PC4[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p42<-qplot(pca$PC2[which(pca$cohort=="1000GP")],pca$PC4[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p43<-qplot(pca$PC3[which(pca$cohort=="1000GP")],pca$PC4[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p44<-qplot(pca$PC4[which(pca$cohort=="1000GP")],pca$PC4[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p45<-qplot(pca$PC5[which(pca$cohort=="1000GP")],pca$PC4[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p46<-qplot(pca$PC6[which(pca$cohort=="1000GP")],pca$PC4[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")

p51<-qplot(pca$PC1[which(pca$cohort=="1000GP")],pca$PC5[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p52<-qplot(pca$PC2[which(pca$cohort=="1000GP")],pca$PC5[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p53<-qplot(pca$PC3[which(pca$cohort=="1000GP")],pca$PC5[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p54<-qplot(pca$PC4[which(pca$cohort=="1000GP")],pca$PC5[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p55<-qplot(pca$PC5[which(pca$cohort=="1000GP")],pca$PC5[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p56<-qplot(pca$PC6[which(pca$cohort=="1000GP")],pca$PC5[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")

p61<-qplot(pca$PC1[which(pca$cohort=="1000GP")],pca$PC6[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p62<-qplot(pca$PC2[which(pca$cohort=="1000GP")],pca$PC6[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p63<-qplot(pca$PC3[which(pca$cohort=="1000GP")],pca$PC6[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p64<-qplot(pca$PC4[which(pca$cohort=="1000GP")],pca$PC6[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p65<-qplot(pca$PC5[which(pca$cohort=="1000GP")],pca$PC6[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p66<-qplot(pca$PC6[which(pca$cohort=="1000GP")],pca$PC6[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")


r1<-ggarrange(pna,p12,p13,p14,p15,p16,ncol=6,legend=c("none"))
r2<-ggarrange(p21,pna,p23,p24,p25,p26,ncol=6,legend=c("none"))
r3<-ggarrange(p31,p32,pna,p34,p35,p36,ncol=6,legend=c("none"))
r4<-ggarrange(p41,p42,p43,pna,p45,p36,ncol=6,legend=c("none"))
r5<-ggarrange(p51,p52,p53,p54,pna,p56,ncol=6,legend=c("none"))
r6<-ggarrange(p61,p62,p63,p64,p65,pna,ncol=6,common.legend = TRUE,legend=c("bottom"))
dev.off()


pdf(paste(path,"pre_imputation/QC/pca_1000gp/pca_PC1_PC2_PC3_PC5_PC6_3kvariants.pdf",sep=""),width=30,height=30)
ggarrange(r1,r2,r3,r4,r5,r6,nrow=6,common.legend = TRUE,legend=c("bottom"),heights = c(1,1,1,1,1,1.2))
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/pca_PC1_PC2_PC3_PC5_PC6_3kvariants.pdf ~/tmp_plots/",sep=""))



#######################################################################
# 4.1 PLOT 1000GP ALONE, AND SEPARATELY , COLOUR BY SUPER POPULATION:


# PC1 EUR - non EUR
# PC2 AFR  - non AFR
# PC3 north EUR -South EUR
# PC4 EAS - SAS
# PC5 Finish - non-Finish
# PC6: AMR - non-AMR
# PC7 - IIBDGC outliers (from more than 1 cohort)
# PC8 transition EUR - EAS - SAS/AFR - AMR

pca$study<-as.factor(pca$study)
pca$study<-factor(pca$study, levels=c("Spain","Italy","Slovenia","Belgium","Germany","UK","Netherlands","Lithuania","Norway","Sweden","Finland",
                                      "CHOP","USA","Australia"))

# PC1 - PC2

ylims<-c(min(pca$PC1),max(pca$PC1))
xlims<-c(min(pca$PC2),max(pca$PC2))

p0<-ggplot(pca[which(pca$cohort=="1000GP"),],aes(PC2,PC1)) +
  geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims)

p1<-ggplot(pca[which(pca$cohort!="1000GP"),],aes(PC2,PC1)) +
  geom_point(aes(color = study)) + 
  xlim(xlims) + ylim(ylims)


pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2.pdf",sep=""),width=14,height=6)
ggarrange(p0,p1,legend=c("right"),widths = c(1,1.1))
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2.pdf ~/tmp_plots/",sep=""))



# PC3 - PC4

ylims<-c(min(pca$PC3),max(pca$PC3))
xlims<-c(min(pca$PC4),max(pca$PC4))

p0<-ggplot(pca[which(pca$cohort=="1000GP"),],aes(PC4,PC3)) +
  geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims)

p1<-ggplot(pca[which(pca$cohort!="1000GP"),],aes(PC4,PC3)) +
  geom_point(aes(color = study)) + 
  xlim(xlims) + ylim(ylims)


pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC3_PC4.pdf",sep=""),width=14,height=6)
ggarrange(p0,p1,legend=c("right"),widths = c(1,1))
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC3_PC4.pdf ~/tmp_plots/",sep=""))


############################################################################################################################################
# INFER POPULATIONS:

pop<-c("AMR","AFR","EAS","SAS","EUR")

pca$inferred_population<-NA

for (i in 1:length(pop)) {
  
  pca$inferred_population[which(pca$PC1>=min(pca$PC1[which(pca$super_pop==pop[i])]) & pca$PC1<=max(pca$PC1[which(pca$super_pop==pop[i])]) &
                                  pca$PC2>=min(pca$PC2[which(pca$super_pop==pop[i])]) & pca$PC2<=max(pca$PC2[which(pca$super_pop==pop[i])]) &
                                  pca$PC3>=min(pca$PC3[which(pca$super_pop==pop[i])]) & pca$PC3<=max(pca$PC3[which(pca$super_pop==pop[i])]) &
                                  pca$PC4>=min(pca$PC4[which(pca$super_pop==pop[i])]) & pca$PC4<=max(pca$PC4[which(pca$super_pop==pop[i])]) &
                                  pca$PC5>=min(pca$PC5[which(pca$super_pop==pop[i])]) & pca$PC5<=max(pca$PC5[which(pca$super_pop==pop[i])]) &
                                  pca$PC6>=min(pca$PC6[which(pca$super_pop==pop[i])]) & pca$PC6<=max(pca$PC6[which(pca$super_pop==pop[i])]) &
                                  pca$PC7>=min(pca$PC7[which(pca$super_pop==pop[i])]) & pca$PC7<=max(pca$PC7[which(pca$super_pop==pop[i])]) &
                                  pca$PC8>=min(pca$PC8[which(pca$super_pop==pop[i])]) & pca$PC8<=max(pca$PC8[which(pca$super_pop==pop[i])]) &
                                  pca$PC9>=min(pca$PC9[which(pca$super_pop==pop[i])]) & pca$PC9<=max(pca$PC9[which(pca$super_pop==pop[i])]) &
                                  pca$PC10>=min(pca$PC10[which(pca$super_pop==pop[i])]) & pca$PC10<=max(pca$PC10[which(pca$super_pop==pop[i])]) )]<-pop[i]
  
}


table(pca$super_pop,pca$inferred_population,useNA="ifany")
#           AFR    AMR    EAS    EUR    SAS   <NA>
# EUR         0      0      0    503      0      0
# AFR       661      0      0      0      0      0
# EAS         0      0    504      0      0      0
# SAS         0      0      0      0    489      0
# AMR        72    252      0     23      0      0
# IIBDGC   1315   1701    639 112040    820   1408


# PC1 - PC2

ylims<-c(min(pca$PC1),max(pca$PC1))
xlims<-c(min(pca$PC2),max(pca$PC2))

p0<-ggplot(pca[which(pca$cohort=="1000GP"),],aes(PC2,PC1)) +
  geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims)

p1<-ggplot(pca[which(is.na(pca$inferred_population)),],aes(PC2,PC1)) +
  geom_point(aes(color = study)) + 
  xlim(xlims) + ylim(ylims)

# PC1 - PC3

ylims<-c(min(pca$PC1),max(pca$PC1))
xlims<-c(min(pca$PC3),max(pca$PC3))

p2<-ggplot(pca[which(pca$cohort=="1000GP"),],aes(PC3,PC1)) +
  geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims)

p3<-ggplot(pca[which(is.na(pca$inferred_population)),],aes(PC3,PC1)) +
  geom_point(aes(color = study)) + 
  xlim(xlims) + ylim(ylims)


c1<-ggarrange(p0,p2,nrow=2,common.legend = TRUE,legend=c("right"))
c2<-ggarrange(p1,p3,nrow=2,common.legend = TRUE,legend=c("right"))
dev.off()

pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2_PC3_noInferredEth.pdf",sep=""),width=14,height=11)
ggarrange(c1,c2,widths = c(1,1),ncol=2)
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2_PC3_noInferredEth.pdf ~/tmp_plots/",sep=""))


####################

# EXPAND LIMITS FOR EUR SUBSET OF SAMPLES

pc1max<-max(pca$PC1[which(pca$super_pop=="EUR")])+((max(pca$PC1[which(pca$super_pop=="EUR")])-min(pca$PC1[which(pca$super_pop=="EUR")]))*0.3)
pc1min<-min(pca$PC1[which(pca$super_pop=="EUR")])-((max(pca$PC1[which(pca$super_pop=="EUR")])-min(pca$PC1[which(pca$super_pop=="EUR")]))*0.3)

pc2max<-max(pca$PC2[which(pca$super_pop=="EUR")])+((max(pca$PC2[which(pca$super_pop=="EUR")])-min(pca$PC2[which(pca$super_pop=="EUR")]))*0.3)
pc2min<-min(pca$PC2[which(pca$super_pop=="EUR")])-((max(pca$PC2[which(pca$super_pop=="EUR")])-min(pca$PC2[which(pca$super_pop=="EUR")]))*0.3)

pc3max<-max(pca$PC3[which(pca$super_pop=="EUR")])+((max(pca$PC3[which(pca$super_pop=="EUR")])-min(pca$PC3[which(pca$super_pop=="EUR")]))*0.3)
pc3min<-min(pca$PC3[which(pca$super_pop=="EUR")])-((max(pca$PC3[which(pca$super_pop=="EUR")])-min(pca$PC3[which(pca$super_pop=="EUR")]))*0.3)


print(paste("PC1 expanded limits",pc1max,pc1min))
# [1] "PC1 expanded limits 0.12924582 0.01132678"
print(paste("PC2 expanded limits",pc2max,pc2min))
# "PC2 expanded limits -0.05502758 -0.17508102"
print(paste("PC3 expanded limits",pc3max,pc3min))
# [1] "PC3 expanded limits 0.03398294 -0.06531114"


# PC1 - PC2 - with new limits


ylims<-c(min(pca$PC1),max(pca$PC1))
xlims<-c(min(pca$PC2),max(pca$PC2))

p0<-ggplot(pca[which(pca$cohort=="1000GP"),],aes(PC2,PC1)) +
  geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims) + 
  geom_hline(yintercept=pc1max, linetype="dashed", color = "grey") + 
  geom_hline(yintercept=pc1min, linetype="dashed", color = "grey") + 
  geom_vline(xintercept=pc2max, linetype="dashed", color = "grey") + 
  geom_vline(xintercept=pc2min, linetype="dashed", color = "grey")


p1<-ggplot(pca[which(is.na(pca$inferred_population)),],aes(PC2,PC1)) +
  geom_point(aes(color = study)) + 
  xlim(xlims) + ylim(ylims) + 
  geom_hline(yintercept=pc1max, linetype="dashed", color = "grey") + 
  geom_hline(yintercept=pc1min, linetype="dashed", color = "grey") + 
  geom_vline(xintercept=pc2max, linetype="dashed", color = "grey") + 
  geom_vline(xintercept=pc2min, linetype="dashed", color = "grey")


# PC1 - PC3

ylims<-c(min(pca$PC1),max(pca$PC1))
xlims<-c(min(pca$PC3),max(pca$PC3))

p2<-ggplot(pca[which(pca$cohort=="1000GP"),],aes(PC3,PC1)) +
  geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims) + 
  geom_hline(yintercept=pc1max, linetype="dashed", color = "grey") + 
  geom_hline(yintercept=pc1min, linetype="dashed", color = "grey") + 
  geom_vline(xintercept=pc3max, linetype="dashed", color = "grey") + 
  geom_vline(xintercept=pc3min, linetype="dashed", color = "grey")

p3<-ggplot(pca[which(is.na(pca$inferred_population)),],aes(PC3,PC1)) +
  geom_point(aes(color = study)) + 
  xlim(xlims) + ylim(ylims) + 
  geom_hline(yintercept=pc1max, linetype="dashed", color = "grey") + 
  geom_hline(yintercept=pc1min, linetype="dashed", color = "grey") + 
  geom_vline(xintercept=pc3max, linetype="dashed", color = "grey") + 
  geom_vline(xintercept=pc3min, linetype="dashed", color = "grey")


c1<-ggarrange(p0,p2,nrow=2,common.legend = TRUE,legend=c("right"))
c2<-ggarrange(p1,p3,nrow=2,common.legend = TRUE,legend=c("right"))
dev.off()

pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2_PC3_noInferredEth_newLim.pdf",sep=""),width=14,height=11)
ggarrange(c1,c2,widths = c(1,1),ncol=2)
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2_PC3_noInferredEth_newLim.pdf ~/tmp_plots/",sep=""))


#### reclasify

pca$inferred_population[which(is.na(pca$inferred_population) & pca$PC1<=pc1max & pca$PC1>=pc1min & pca$PC2<=pc2max & pca$PC2>=pc2min
                              & pca$PC3<=pc3max & pca$PC3>=pc3min)]<-"EUR"



table(pca$super_pop,pca$inferred_population,useNA="ifany")
#           AFR    AMR    EAS    EUR    SAS   <NA>
# EUR         0      0      0    503      0      0
# AFR       661      0      0      0      0      0
# EAS         0      0    504      0      0      0
# SAS         0      0      0      0    489      0
# AMR        72    252      0     23      0      0
# IIBDGC   1315   1701    639 113112    820    336


# 
# # PC1 - PC2 - with new limits
# 
# 
# ylims<-c(min(pca$PC1),max(pca$PC1))
# xlims<-c(min(pca$PC2),max(pca$PC2))
# 
# p0<-ggplot(pca[which(pca$cohort=="1000GP"),],aes(PC2,PC1)) +
#   geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") +
#   xlim(xlims) + ylim(ylims) +
#   geom_hline(yintercept=pc1max, linetype="dashed", color = "grey") +
#   geom_hline(yintercept=pc1min, linetype="dashed", color = "grey") +
#   geom_vline(xintercept=pc2max, linetype="dashed", color = "grey") +
#   geom_vline(xintercept=pc2min, linetype="dashed", color = "grey")
# 
# 
# p1<-ggplot(pca[which(is.na(pca$inferred_population)),],aes(PC2,PC1)) +
#   geom_point(aes(color = study)) +
#   xlim(xlims) + ylim(ylims) +
#   geom_hline(yintercept=pc1max, linetype="dashed", color = "grey") +
#   geom_hline(yintercept=pc1min, linetype="dashed", color = "grey") +
#   geom_vline(xintercept=pc2max, linetype="dashed", color = "grey") +
#   geom_vline(xintercept=pc2min, linetype="dashed", color = "grey")
# 
# 
# # PC1 - PC3
# 
# ylims<-c(min(pca$PC1),max(pca$PC1))
# xlims<-c(min(pca$PC3),max(pca$PC3))
# 
# p2<-ggplot(pca[which(pca$cohort=="1000GP"),],aes(PC3,PC1)) +
#   geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") +
#   xlim(xlims) + ylim(ylims) +
#   geom_hline(yintercept=pc1max, linetype="dashed", color = "grey") +
#   geom_hline(yintercept=pc1min, linetype="dashed", color = "grey") +
#   geom_vline(xintercept=pc3max, linetype="dashed", color = "grey") +
#   geom_vline(xintercept=pc3min, linetype="dashed", color = "grey")
# 
# p3<-ggplot(pca[which(is.na(pca$inferred_population)),],aes(PC3,PC1)) +
#   geom_point(aes(color = study)) +
#   xlim(xlims) + ylim(ylims) +
#   geom_hline(yintercept=pc1max, linetype="dashed", color = "grey") +
#   geom_hline(yintercept=pc1min, linetype="dashed", color = "grey") +
#   geom_vline(xintercept=pc3max, linetype="dashed", color = "grey") +
#   geom_vline(xintercept=pc3min, linetype="dashed", color = "grey")
# 
# 
# c1<-ggarrange(p0,p2,nrow=2,common.legend = TRUE,legend=c("right"))
# c2<-ggarrange(p1,p3,nrow=2,common.legend = TRUE,legend=c("right"))
# dev.off()
# 
# pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2_PC3_noInferredEth_newLim_2.pdf",sep=""),width=14,height=11)
# ggarrange(c1,c2,widths = c(1,1),ncol=2)
# dev.off()
# system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2_PC3_noInferredEth_newLim_2.pdf ~/tmp_plots/",sep=""))
# 


####################

# EXPAND LIMITS FOR SAS SUBSET OF SAMPLES

pc1max_sas<-max(pca$PC1[which(pca$super_pop=="SAS")])+((max(pca$PC1[which(pca$super_pop=="SAS")])-min(pca$PC1[which(pca$super_pop=="SAS")]))*0.3)
pc1min_sas<-min(pca$PC1[which(pca$super_pop=="SAS")])-((max(pca$PC1[which(pca$super_pop=="SAS")])-min(pca$PC1[which(pca$super_pop=="SAS")]))*0.3)

pc2max_sas<-max(pca$PC2[which(pca$super_pop=="SAS")])+((max(pca$PC2[which(pca$super_pop=="SAS")])-min(pca$PC2[which(pca$super_pop=="SAS")]))*0.3)
pc2min_sas<-min(pca$PC2[which(pca$super_pop=="SAS")])-((max(pca$PC2[which(pca$super_pop=="SAS")])-min(pca$PC2[which(pca$super_pop=="SAS")]))*0.3)

pc3max_sas<-max(pca$PC3[which(pca$super_pop=="SAS")])+((max(pca$PC3[which(pca$super_pop=="SAS")])-min(pca$PC3[which(pca$super_pop=="SAS")]))*0.3)
pc3min_sas<-min(pca$PC3[which(pca$super_pop=="SAS")])-((max(pca$PC3[which(pca$super_pop=="SAS")])-min(pca$PC3[which(pca$super_pop=="SAS")]))*0.3)


print(paste("PC1 expanded limits",pc1max_sas,pc1min_sas))
# [1] "PC1 expanded limits 0.13417616 0.02434064"
print(paste("PC2 expanded limits",pc2max_sas,pc2min_sas))
# [1] "PC2 expanded limits 0.05007192 -0.12805192"
print(paste("PC3 expanded limits",pc3max_sas,pc3min_sas))
# [1] "PC3 expanded limits 0.14906881 0.01080849"

pca$inferred_population[which(is.na(pca$inferred_population) & pca$PC1<=pc1max_sas & pca$PC1>=pc1min_sas & pca$PC2<=pc2max_sas & pca$PC2>=pc2min_sas
                              & pca$PC3<=pc3max_sas & pca$PC3>=pc3min_sas)]<-"SAS"

table(pca$super_pop,pca$inferred_population,useNA="ifany")
#           AFR    AMR    EAS    EUR    SAS   <NA>
# EUR         0      0      0    503      0      0
# AFR       661      0      0      0      0      0
# EAS         0      0    504      0      0      0
# SAS         0      0      0      0    489      0
# AMR        72    252      0     23      0      0
# IIBDGC   1315   1701    639 113112    909    247

####################

# EXPAND LIMITS FOR EAS SUBSET OF SAMPLES

pc1max_eas<-max(pca$PC1[which(pca$super_pop=="EAS")])+((max(pca$PC1[which(pca$super_pop=="EAS")])-min(pca$PC1[which(pca$super_pop=="EAS")]))*0.3)
pc1min_eas<-min(pca$PC1[which(pca$super_pop=="EAS")])-((max(pca$PC1[which(pca$super_pop=="EAS")])-min(pca$PC1[which(pca$super_pop=="EAS")]))*0.3)

pc2max_eas<-max(pca$PC2[which(pca$super_pop=="EAS")])+((max(pca$PC2[which(pca$super_pop=="EAS")])-min(pca$PC2[which(pca$super_pop=="EAS")]))*0.3)
pc2min_eas<-min(pca$PC2[which(pca$super_pop=="EAS")])-((max(pca$PC2[which(pca$super_pop=="EAS")])-min(pca$PC2[which(pca$super_pop=="EAS")]))*0.3)

pc3max_eas<-max(pca$PC3[which(pca$super_pop=="EAS")])+((max(pca$PC3[which(pca$super_pop=="EAS")])-min(pca$PC3[which(pca$super_pop=="EAS")]))*0.3)
pc3min_eas<-min(pca$PC3[which(pca$super_pop=="EAS")])-((max(pca$PC3[which(pca$super_pop=="EAS")])-min(pca$PC3[which(pca$super_pop=="EAS")]))*0.3)


print(paste("PC1 expanded limits",pc1max_eas,pc1min_eas))
# [1] "PC1 expanded limits 0.16424364 0.08424556"
print(paste("PC2 expanded limits",pc2max_eas,pc2min_eas))
#[1] "PC2 expanded limits 0.206812 0.104236"
print(paste("PC3 expanded limits",pc3max_eas,pc3min_eas))
# [1] "PC3 expanded limits 0.06467009 -0.05431119"

pca$inferred_population[which(is.na(pca$inferred_population) & pca$PC1<=pc1max_eas & pca$PC1>=pc1min_eas & pca$PC2<=pc2max_eas & pca$PC2>=pc2min_eas
                              & pca$PC3<=pc3max_eas & pca$PC3>=pc3min_eas)]<-"EAS"

table(pca$super_pop,pca$inferred_population,useNA="ifany")
#           AFR    AMR    EAS    EUR    SAS   <NA>
# EUR         0      0      0    503      0      0
# AFR       661      0      0      0      0      0
# EAS         0      0    504      0      0      0
# SAS         0      0      0      0    489      0
# AMR        72    252      0     23      0      0
# IIBDGC   1315   1701    667 113112    909    219

####################

####################

# EXPAND LIMITS FOR AFR SUBSET OF SAMPLES

pc1max_afr<-max(pca$PC1[which(pca$super_pop=="AFR")])+((max(pca$PC1[which(pca$super_pop=="AFR")])-min(pca$PC1[which(pca$super_pop=="AFR")]))*0.3)
pc1min_afr<-min(pca$PC1[which(pca$super_pop=="AFR")])-((max(pca$PC1[which(pca$super_pop=="AFR")])-min(pca$PC1[which(pca$super_pop=="AFR")]))*0.3)

pc2max_afr<-max(pca$PC2[which(pca$super_pop=="AFR")])+((max(pca$PC2[which(pca$super_pop=="AFR")])-min(pca$PC2[which(pca$super_pop=="AFR")]))*0.3)
pc2min_afr<-min(pca$PC2[which(pca$super_pop=="AFR")])-((max(pca$PC2[which(pca$super_pop=="AFR")])-min(pca$PC2[which(pca$super_pop=="AFR")]))*0.3)

pc3max_afr<-max(pca$PC3[which(pca$super_pop=="AFR")])+((max(pca$PC3[which(pca$super_pop=="AFR")])-min(pca$PC3[which(pca$super_pop=="AFR")]))*0.3)
pc3min_afr<-min(pca$PC3[which(pca$super_pop=="AFR")])-((max(pca$PC3[which(pca$super_pop=="AFR")])-min(pca$PC3[which(pca$super_pop=="AFR")]))*0.3)


print(paste("PC1 expanded limits",pc1max_afr,pc1min_afr))
# [1] "PC1 expanded limits 0.2187623 -0.4029193"
print(paste("PC2 expanded limits",pc2max_afr,pc2min_afr))
# [1] "PC2 expanded limits 0.08841407 -0.07793137"
print(paste("PC3 expanded limits",pc3max_afr,pc3min_afr))
# [1] "PC3 expanded limits 0.09617305 -0.18353655"

pca$inferred_population[which(is.na(pca$inferred_population) & pca$PC1<=pc1max_afr & pca$PC1>=pc1min_afr & pca$PC2<=pc2max_afr & pca$PC2>=pc2min_afr
                              & pca$PC3<=pc3max_afr & pca$PC3>=pc3min_afr)]<-"AFR"

table(pca$super_pop,pca$inferred_population,useNA="ifany")
#           AFR    AMR    EAS    EUR    SAS   <NA>
# EUR         0      0      0    503      0      0
# AFR       661      0      0      0      0      0
# EAS         0      0    504      0      0      0
# SAS         0      0      0      0    489      0
# AMR        72    252      0     23      0      0
# IIBDGC   1485   1701    667 113112    909     49

####################

# EXPAND LIMITS FOR AMR SUBSET OF SAMPLES

pc1max_amr<-max(pca$PC1[which(pca$super_pop=="AMR")])+((max(pca$PC1[which(pca$super_pop=="AMR")])-min(pca$PC1[which(pca$super_pop=="AMR")]))*0.3)
pc1min_amr<-min(pca$PC1[which(pca$super_pop=="AMR")])-((max(pca$PC1[which(pca$super_pop=="AMR")])-min(pca$PC1[which(pca$super_pop=="AMR")]))*0.3)

pc2max_amr<-max(pca$PC2[which(pca$super_pop=="AMR")])+((max(pca$PC2[which(pca$super_pop=="AMR")])-min(pca$PC2[which(pca$super_pop=="AMR")]))*0.3)
pc2min_amr<-min(pca$PC2[which(pca$super_pop=="AMR")])-((max(pca$PC2[which(pca$super_pop=="AMR")])-min(pca$PC2[which(pca$super_pop=="AMR")]))*0.3)

pc3max_amr<-max(pca$PC3[which(pca$super_pop=="AMR")])+((max(pca$PC3[which(pca$super_pop=="AMR")])-min(pca$PC3[which(pca$super_pop=="AMR")]))*0.3)
pc3min_amr<-min(pca$PC3[which(pca$super_pop=="AMR")])-((max(pca$PC3[which(pca$super_pop=="AMR")])-min(pca$PC3[which(pca$super_pop=="AMR")]))*0.3)


print(paste("PC1 expanded limits",pc1max_amr,pc1min_amr))
# [1] "PC1 expanded limits 0.2300459 -0.2356709"
print(paste("PC2 expanded limits",pc2max_amr,pc2min_amr))
# [1] "PC2 expanded limits 0.14405716 -0.19010796"
print(paste("PC3 expanded limits",pc3max_amr,pc3min_amr))
# [1] "PC3 expanded limits 0.064479033 -0.295113623"

pca$inferred_population[which(is.na(pca$inferred_population) & pca$PC1<=pc1max_amr & pca$PC1>=pc1min_amr & pca$PC2<=pc2max_amr & pca$PC2>=pc2min_amr
                              & pca$PC3<=pc3max_amr & pca$PC3>=pc3min_amr)]<-"AMR"

table(pca$super_pop,pca$inferred_population,useNA="ifany")
#           AFR    AMR    EAS    EUR    SAS
# EUR         0      0      0    503      0
# AFR       661      0      0      0      0
# EAS         0      0    504      0      0
# SAS         0      0      0      0    489
# AMR        72    252      0     23      0
# IIBDGC   1485   1750    667 113112    909


pca$inferred_population<-as.factor(pca$inferred_population)
pca$inferred_population<-factor(pca$inferred_population, levels=c("EUR","AFR","EAS","SAS","AMR"))

###########################################
# some inferred AFR in AMR, reclassify:

table(pca[which(pca$PC3< -0.05),"super_pop"])
# EUR    AFR    EAS    SAS    AMR IIBDGC 
#   0      3      0      0    253    709   
table(pca[which(pca$PC3< -0.05),"inferred_population"])
# EUR AFR EAS SAS AMR 
# 29 371   0   0 565 

pca$inferred_population[which((pca$inferred_population=="AFR") & (pca$PC3< -0.05) )]<-"AMR"
table(pca$super_pop,pca$inferred_population,useNA="ifany")
#           EUR    AFR    EAS    SAS    AMR
# EUR       503      0      0      0      0
# AFR         0    658      0      0      3
# EAS         0      0    504      0      0
# SAS         0      0      0    489      0
# AMR        23      8      0      0    316
# IIBDGC 113112   1181    667    909   2054


table(pca[which(pca$PC1< -0.05),"super_pop"])
# EUR    AFR    EAS    SAS    AMR IIBDGC 
#   0    658      0      0      5   1037 
table(pca[which(pca$PC1> -0.05),"super_pop"])
# EUR    AFR    EAS    SAS    AMR IIBDGC 
# 503      3    504    489    342 116886

table(pca[which(pca$PC1< -0.05),"inferred_population"])
# EUR  AFR  EAS  SAS  AMR 
# 0 1591    0    0  109

pca$inferred_population[which((pca$inferred_population=="AMR") & (pca$PC1< -0.05) )]<-"AFR"
table(pca$super_pop,pca$inferred_population,useNA="ifany")
#           EUR    AFR    EAS    SAS    AMR
# EUR       503      0      0      0      0
# AFR         0    658      0      0      3
# EAS         0      0    504      0      0
# SAS         0      0      0    489      0
# AMR        23      9      0      0    315
# IIBDGC 113112   1289    667    909   1946


pca$inferred_population[which((pca$inferred_population=="AFR") & (pca$PC1> -0.05) )]<-"AMR"
table(pca$super_pop,pca$inferred_population,useNA="ifany")
#           EUR    AFR    EAS    SAS    AMR
# EUR       503      0      0      0      0
# AFR         0    658      0      0      3
# EAS         0      0    504      0      0
# SAS         0      0      0    489      0
# AMR        23      5      0      0    319
# IIBDGC 113112   1037    667    909   2198

# PLOT FINAL RESULTS

# PC1 - PC2 - with new limits

ylims<-c(min(pca$PC1),max(pca$PC1))
xlims<-c(min(pca$PC2),max(pca$PC2))

p0<-ggplot(pca[which(pca$cohort=="1000GP"),],aes(PC2,PC1)) +
  geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims) 


p1<-ggplot(pca[which(pca$cohort!="1000GP"),],aes(PC2,PC1)) +
  geom_point(aes(color = inferred_population)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims) 


# PC1 - PC3

ylims<-c(min(pca$PC1),max(pca$PC1))
xlims<-c(min(pca$PC3),max(pca$PC3))

p2<-ggplot(pca[which(pca$cohort=="1000GP"),],aes(PC3,PC1)) +
  geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims) 

p3<-ggplot(pca[which(pca$cohort!="1000GP"),],aes(PC3,PC1)) +
  geom_point(aes(color = inferred_population)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims) 


c1<-ggarrange(p0,p2,nrow=2,common.legend = TRUE,legend=c("right"))
c1<-annotate_figure(c1,top = text_grob("1000GP"))
c2<-ggarrange(p1,p3,nrow=2,common.legend = TRUE,legend=c("right"))
c2<-annotate_figure(c2,top = text_grob("IIBDGC"))
dev.off()

pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2_PC3_inferred_final.pdf",sep=""),width=14,height=11)
ggarrange(c1,c2,widths = c(1,1.1),ncol=2)
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2_PC3_inferred_final.pdf ~/tmp_plots/",sep=""))


#########

table(pca$cohort[which(pca$PHENO1=="2")],pca$inferred_population[which(pca$PHENO1=="2")],useNA="ifany")
#                           EUR   AFR   EAS   SAS   AMR
# all_hce                 12342   110   325   494   164
# australia_omniexome       659     1     4    14    13
# basque_gsa                518     0     0     0     7
# belgium_franchimont_gsa   880     2     0     2    22
# belgium_inf1_old_gwas     513     0     0     0     4
# belgium_inf2_old_gwas     159     0     0     0     2
# belgium_louis_gsa         897     2     0     0     7
# belgium_vermeire_gsa     3098    10     3     9    59
# ccfa_gsa                 1931   144    12    25    65
# cedars_370k_old_gwas      509    13     9    12    62
# cedars_610k_old_gwas      830    15     6     2    36
# cedars_gsa               1884    40    43    36   128
# cedars_omni_old_gwas     1104    36     9     9    59
# chop_old_gwas            2371     0     0     1    41
# finland_illugwas          443     1     0     0     2
# german_affy6_old_gwas    1045     0     0     0    10
# gwas1                    1740     0     0     0     7
# gwas2                    2353     0     0     0     8
# italy_gsa                 576     1     0     0     1
# kiel_austria_sibdcs_gsa  9597    29    12    52   129
# lithuania_gsa            1072     0     0     0    16
# mccauley_gsa              760     1     0     0    21
# netherlands_gsa          3818    38    13    42   120
# niddk_broad_gsa          4157    66    30    42   143
# niddk_feinstein_gsa      5186   140    45    77   305
# niddk_old_gwas           1808     1     0     0    16
# norway_affy6_old_gwas     266     0     0     0     2
# pittsburgh_gsa           1251     0     0     0     6
# prism_nfe_gsa             401    11     2     6    13
# prism_nfe_gwas            510    18     2    10    16
# slovenia_gsa               87     0     0     0     0
# spain_gsa                1957     0     0     0     4
# sweden_gsa                390     5     2     5    15
# swedish_uc_old_gwas       921     0     0     0     2

table(pca$cohort[which(pca$PHENO1=="1")],pca$inferred_population[which(pca$PHENO1=="1")],useNA="ifany")
#                           EUR   AFR   EAS   SAS   AMR
# all_hce                 10356    12     3    24    60
# australia_omniexome       592     1     9     2    11
# basque_gsa                981     0     0     0     2
# belgium_franchimont_gsa   586     0     0     0     9
# belgium_inf1_old_gwas     891     0     0     0     9
# belgium_inf2_old_gwas     111     0     0     0     0
# belgium_louis_gsa         592     2     0     0     2
# belgium_vermeire_gsa      807     2     0     0     5
# cedars_gsa                801     2    67    10    67
# cedars_omni_old_gwas        3     0     0     0     1
# chop_old_gwas            6104     0     0     2    91
# german_affy6_old_gwas    1759     3     1     0     5
# gwas1                    2923     0     0     1     9
# gwas2                    5400     0     0     0    17
# italy_gsa                 371     1     0     0     3
# kiel_austria_sibdcs_gsa  4536     5     4     2    37
# lithuania_gsa            1126     0     0     0    17
# netherlands_gsa           534     2     3     1    10
# niddk_broad_gsa           991     7     3     0    16
# niddk_feinstein_gsa      1781   303    45    13   282
# niddk_old_gwas            922     2     1     0    10
# norway_affy6_old_gwas     281     0     0     0     1
# pittsburgh_gsa           1464     1     0     0     3
# prism_nfe_gsa              28     0     1     1     3
# prism_nfe_gwas            210    10    13    15    14
# slovenia_gsa              177     0     0     0     0
# spain_gsa                1479     0     0     0     3
# sweden_gsa                934     0     0     0     4
# swedish_uc_old_gwas       339     0     0     0     2

write.table(pca,paste(path,"pre_imputation/QC/pca_1000gp/1000GP_all_b37_plus_iibdgc_subset_iibdgc_new_projection_edited_withDuplicates.sscore",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")

### WHERE SELF-REPORTED ASKENAZI ARE:

# # only new cohorts and cedars will have data for this:
# pheno<-read.csv("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes.csv",head=T)
# colnames(pheno)[13]<-"self_jewish"
# pheno_cd<-fread(paste(path,"pheno/list_cedars_old_gwas_ids_TH.txt",sep=""),head=T,sep="\t")
# 
# pheno<-rbind(pheno[,c("sample_id","self_jewish")],pheno_cd[,c("sample_id","self_jewish")])
# table(pheno$self_jewish)
# #               No    Unknown        Yes     Jewish Non-Jewish 
# # 21749      31383      15431       3431       1217       1513 
# 
# pheno$self_jewish<-as.character(pheno$self_jewish)
# pheno$self_jewish[which(!pheno$self_jewish %in% c("Yes","Jewish"))]<-"0"
# pheno$self_jewish[which(pheno$self_jewish %in% c("Yes","Jewish"))]<-"1"
# table(pheno$self_jewish)
# # 0     1 
# # 70076  4648
# 
# 3431+1217
# # [1] 4648
# 
# 
# pca<-merge(pca,pheno,by.x="FID",by.y="sample_id",all.x=T)
# table(pca$self_jewish,useNA="ifany")
#   #     0     1  <NA> 
#   # 51171  3860 49602 
# 
# table(pca$inferred_population,pca$self_jewish,useNA="ifany")
# #          0     1  <NA>
# # AFR    950     9   851
# # AMR    969    48   690
# # EAS    186     0   824
# # EUR  48627  3794 46242
# # SAS    246     0   916
# # <NA>   193     9    79
# 
# 
# 
# pca$self_jewish[which(is.na(pca$self_jewish))]<-"0"
# table(pca$self_jewish)
# # 0     1 
# # 100773  3860
# 
# table(pca$inferred_population,pca$self_jewish,useNA="ifany")
# #          0     1
# # AFR   1801     9
# # AMR   1659    48
# # EAS   1010     0
# # EUR  94869  3794
# # SAS   1162     0
# # <NA>   272     9
# 
# 
# p11<-qplot(pca$PC1[which(pca$cohort!="1000GP")],pca$PC1[which(pca$cohort!="1000GP")], data = pca[which(pca$cohort!="1000GP"),], colour = self_jewish)
# p12<-qplot(pca$PC2[which(pca$cohort!="1000GP")],pca$PC1[which(pca$cohort!="1000GP")], data = pca[which(pca$cohort!="1000GP"),], colour = self_jewish)
# p13<-qplot(pca$PC3[which(pca$cohort!="1000GP")],pca$PC1[which(pca$cohort!="1000GP")], data = pca[which(pca$cohort!="1000GP"),], colour = self_jewish)
# 
# p21<-qplot(pca$PC1[which(pca$cohort!="1000GP")],pca$PC2[which(pca$cohort!="1000GP")], data = pca[which(pca$cohort!="1000GP"),], colour = self_jewish)
# p22<-qplot(pca$PC2[which(pca$cohort!="1000GP")],pca$PC2[which(pca$cohort!="1000GP")], data = pca[which(pca$cohort!="1000GP"),], colour = self_jewish)
# p23<-qplot(pca$PC3[which(pca$cohort!="1000GP")],pca$PC2[which(pca$cohort!="1000GP")], data = pca[which(pca$cohort!="1000GP"),], colour = self_jewish)
# 
# p31<-qplot(pca$PC1[which(pca$cohort!="1000GP")],pca$PC3[which(pca$cohort!="1000GP")], data = pca[which(pca$cohort!="1000GP"),], colour = self_jewish)
# p32<-qplot(pca$PC2[which(pca$cohort!="1000GP")],pca$PC3[which(pca$cohort!="1000GP")], data = pca[which(pca$cohort!="1000GP"),], colour = self_jewish)
# p33<-qplot(pca$PC3[which(pca$cohort!="1000GP")],pca$PC3[which(pca$cohort!="1000GP")], data = pca[which(pca$cohort!="1000GP"),], colour = self_jewish)
# 
# r1<-ggarrange(pna,p12,p13,ncol=3,legend=c("none"))
# r2<-ggarrange(p21,pna,p23,ncol=3,legend=c("none"))
# r3<-ggarrange(p31,p32,pna,ncol=3,common.legend = TRUE,legend=c("bottom"))
# dev.off()
# 
# 
# pdf(paste(path,"pre_imputation/QC/pca_1000gp/pca_PC1_PC2_PC3_3kvariants_sr_jewish.pdf",sep=""),width=15,height=15)
# ggarrange(r1,r2,r3,nrow=3,common.legend = TRUE,legend=c("bottom"),heights = c(1,1,1,1,1,1.2))
# dev.off()
# system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/pca_PC1_PC2_PC3_3kvariants_sr_jewish.pdf ~/tmp_plots/",sep=""))


#########################################################
# 5.- EXTRACT EUR SAMPLES AND IDENTIFY ASHKENAZI SUBSET #
#########################################################

pca<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/1000GP_all_b37_plus_iibdgc_subset_iibdgc_new_projection_edited_withDuplicates.sscore",sep=""),head=T)

table(pca$super_pop,pca$inferred_population,useNA="ifany")
#           AFR    AMR    EAS    EUR    SAS
# AFR       658      3      0      0      0
# AMR         5    319      0     23      0
# EAS         0      0    504      0      0
# EUR         0      0      0    503      0
# IIBDGC   1037   2198    667 113112    909
# SAS         0      0      0      0    489



write.table(pca[which(pca$inferred_population=="EUR"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_eur_ancestry_samples_withDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

write.table(pca[which(pca$inferred_population=="AFR"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_afr_ancestry_samples_withDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

write.table(pca[which(pca$inferred_population=="EAS"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_eas_ancestry_samples_withDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

write.table(pca[which(pca$inferred_population=="SAS"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_sas_ancestry_samples_withDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

write.table(pca[which(pca$inferred_population=="AMR"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_amr_ancestry_samples_withDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

#### create lists with no duplicates:

fam<-read.table(paste(path,"pre_imputation/QC/relatedness/iibdgc_merged_nodup.fam",sep=""),head=F)

pca_nodup<-pca[which(pca$IID %in% fam$V2),]
table(pca_nodup$super_pop,pca_nodup$inferred_population,useNA="ifany")
#           AFR    AMR    EAS    EUR    SAS
# AFR         0      0      0      0      0
# AMR         0      0      0      0      0
# EAS         0      0      0      0      0
# EUR         0      0      0      0      0
# IIBDGC    957   1933    605 102134    796
# SAS         0      0      0      0      0

write.table(pca_nodup[which(pca_nodup$inferred_population=="EUR"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_eur_ancestry_samples_noDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

write.table(pca_nodup[which(pca_nodup$inferred_population=="AFR"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_afr_ancestry_samples_noDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

write.table(pca_nodup[which(pca_nodup$inferred_population=="EAS"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_eas_ancestry_samples_noDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

write.table(pca_nodup[which(pca_nodup$inferred_population=="SAS"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_sas_ancestry_samples_noDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")

write.table(pca_nodup[which(pca_nodup$inferred_population=="AMR"),c("FID","IID")],paste(path,"pre_imputation/QC/pca_1000gp/list_amr_ancestry_samples_noDuplicates",sep=""),
            col.names=F,row.names=F,quote=F,sep="\t")



table(pca$cohort[which(pca$PHENO1=="1")],pca$inferred_population[which(pca$PHENO1=="1")],useNA="ifany")
table(pca$cohort[which(pca$PHENO1=="1")],pca$inferred_population[which(pca$PHENO1=="2")],useNA="ifany")


table(pca_nodup$cohort[which(pca_nodup$PHENO1=="1")],pca_nodup$inferred_population[which(pca_nodup$PHENO1=="1")],useNA="ifany")
#                           AFR   AMR   EAS   EUR   SAS
# 1000GP                      0     0     0     0     0
# all_hce                    12    60     3 10345    24
# australia_omniexome         1    11     9   589     2
# basque_gsa                  0     2     0   964     0
# belgium_franchimont_gsa     0     9     0   583     0
# belgium_inf1_old_gwas       0     8     0   831     0
# belgium_inf2_old_gwas       0     0     0    77     0
# belgium_louis_gsa           2     2     0   591     0
# belgium_vermeire_gsa        2     5     0   806     0
# ccfa_gsa                    0     0     0     0     0
# cedars_370k_old_gwas        0     0     0     0     0
# cedars_610k_old_gwas        0     0     0     0     0
# cedars_gsa                  2    51    55   641     8
# cedars_omni_old_gwas        0     1     0     1     0
# chop_old_gwas               0    90     0  6066     2
# finland_illugwas            0     0     0     0     0
# german_affy6_old_gwas       3     5     1  1659     0
# gwas1                       0     9     0  2919     1
# gwas2                       0     9     0  2767     0
# italy_gsa                   1     3     0   370     0
# kiel_austria_sibdcs_gsa     5    35     4  4409     2
# lithuania_gsa               0    17     0  1122     0
# mccauley_gsa                0     0     0     0     0
# netherlands_gsa             2    10     3   532     1
# niddk_broad_gsa             7    16     3   986     0
# niddk_feinstein_gsa       291   268    44  1735    12
# niddk_old_gwas              2     5     1   666     0
# norway_affy6_old_gwas       0     1     0   281     0
# pittsburgh_gsa              1     2     0  1375     0
# prism_nfe_gsa               0     3     1    28     1
# prism_nfe_gwas             10    14    13   205    15
# slovenia_gsa                0     0     0   175     0
# spain_gsa                   0     3     0  1479     0
# sweden_gsa                  0     4     0   927     0
# swedish_uc_old_gwas         0     2     0   325     0

table(pca_nodup$cohort[which(pca_nodup$PHENO1=="2")],pca_nodup$inferred_population[which(pca_nodup$PHENO1=="2")],useNA="ifany")
#                           AFR   AMR   EAS   EUR   SAS
# 1000GP                      0     0     0     0     0
# all_hce                    95   139   321 10996   431
# australia_omniexome         1    13     4   654    14
# basque_gsa                  0     7     0   518     0
# belgium_franchimont_gsa     1    22     0   872     2
# belgium_inf1_old_gwas       0     3     0   228     0
# belgium_inf2_old_gwas       0     2     0    82     0
# belgium_louis_gsa           2     7     0   886     0
# belgium_vermeire_gsa       10    59     3  3085     9
# ccfa_gsa                  143    63    12  1821    25
# cedars_370k_old_gwas        5    42     3   395     7
# cedars_610k_old_gwas       12    14     1   440     0
# cedars_gsa                 22    64    17   927    15
# cedars_omni_old_gwas       33    30     2   600     6
# chop_old_gwas               0    29     0  1707     1
# finland_illugwas            1     2     0   441     0
# german_affy6_old_gwas       0    10     0   953     0
# gwas1                       0     7     0  1740     0
# gwas2                       0     8     0  2348     0
# italy_gsa                   1     1     0   573     0
# kiel_austria_sibdcs_gsa    27   120    12  9210    46
# lithuania_gsa               0    16     0  1055     0
# mccauley_gsa                1    21     0   747     0
# netherlands_gsa            38   120    13  3794    42
# niddk_broad_gsa            52   133    29  3859    37
# niddk_feinstein_gsa       138   300    45  5054    73
# niddk_old_gwas              1     3     0   606     0
# norway_affy6_old_gwas       0     2     0   266     0
# pittsburgh_gsa              0     2     0   650     0
# prism_nfe_gsa              11    13     2   399     6
# prism_nfe_gwas             17    15     2   435     9
# slovenia_gsa                0     0     0    87     0
# spain_gsa                   0     4     0  1957     0
# sweden_gsa                  5    15     2   376     5
# swedish_uc_old_gwas         0     2     0   919     0


