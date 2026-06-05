# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
######################################################
# 20 - EXCLUDE GSA VARIANTS FLAGGED BY GENOME STUDIO #
######################################################

# GSA array specific step

## /software/R-4.3.1/bin/R

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
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","swedish_uc_old_gwas"
           ,"mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

length(cohorts)
# [1] 58

cohorts_gsa<-c("italy_gsa","kiel_austria_sibdcs_gsa"
               ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
               ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
               ,"mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
               "franke_gsa","hyams_protect_gsa","lewis_sparc_gsa",
               "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
               "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
               "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
               "xavier_share_gsa")

for (j in 1:length(cohorts)) {
  
  print(cohorts[j])
  
  bim<-read.table(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck.bim",sep=""),head=F)
  
  if (cohorts[j] %in% cohorts_gsa) {
    fl<-read.table(paste(path,"from_kyle/list_variants_flagged_genome_studio_b38.tsv",sep=""),head=F,sep="\t")
    fl$ID1<-paste(fl$V2,":",fl$V3,"_",fl$V4,"_",fl$V5,sep="")
    fl$ID2<-paste(fl$V2,":",fl$V3,"_",fl$V5,"_",fl$V4,sep="")
    fl$X<-paste(fl$V2,":",fl$V3,sep="")
    
    fl$a1<-""
    fl$a1[which(fl$V4=="A")]<-"T"
    fl$a1[which(fl$V4=="T")]<-"A"
    fl$a1[which(fl$V4=="C")]<-"G"
    fl$a1[which(fl$V4=="G")]<-"C"
    
    fl$a2<-""
    fl$a2[which(fl$V5=="A")]<-"T"
    fl$a2[which(fl$V5=="T")]<-"A"
    fl$a2[which(fl$V5=="C")]<-"G"
    fl$a2[which(fl$V5=="G")]<-"C"
    
    fl$ID3<-paste(fl$V2,":",fl$V3,"_",fl$a1,"_",fl$a2,sep="")
    fl$ID4<-paste(fl$V2,":",fl$V3,"_",fl$a2,"_",fl$a1,sep="")
    
    
    excl<-bim[which( (bim$V2 %in% fl$ID1) | (bim$V2 %in% fl$ID2) | (bim$V2 %in% fl$ID3) | (bim$V2 %in% fl$ID4)),]
    print(paste("N variants to exclude:",nrow(excl),sep=""))
  }else{
    excl<-as.data.frame(matrix(ncol=1,nrow=0))
    colnames(excl)<-"V2"
  }
  
  write.table(excl[,"V2",drop=F],paste(path,"pre_imputation/QC/",cohorts[j],"/list_gsa_variants_flagged_to_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  rm(list=ls()[!ls() %in% c("path","j","cohorts","cohorts_gsa")]) 
}

#########
# [1] "australia_omniexome"
# [1] "gwas1"
# [1] "gwas2"
# [1] "all_hce"
# [1] "pittsburgh_gsa"
# [1] "spain_gsa"
# [1] "italy_gsa"
# [1] "N variants to exclude:40877"
# [1] "kiel_austria_sibdcs_gsa"
# [1] "N variants to exclude:21059"
# [1] "netherlands_gsa"
# [1] "N variants to exclude:45050"
# [1] "slovenia_gsa"
# [1] "N variants to exclude:34977"
# [1] "sweden_gsa"
# [1] "N variants to exclude:41923"
# [1] "niddk_broad_gsa"
# [1] "N variants to exclude:33799"
# [1] "niddk_feinstein_gsa"
# [1] "N variants to exclude:32590"
# [1] "basque_gsa"
# [1] "N variants to exclude:40110"
# [1] "prism_nfe_gsa"
# [1] "N variants to exclude:41241"
# [1] "lithuania_gsa"
# [1] "N variants to exclude:37422"
# [1] "belgium_louis_gsa"
# [1] "N variants to exclude:44437"
# [1] "belgium_franchimont_gsa"
# [1] "N variants to exclude:38891"
# [1] "belgium_vermeire_gsa"
# [1] "N variants to exclude:43728"
# [1] "prism_nfe_gwas"
# [1] "finland_illugwas"
# [1] "german_affy6_old_gwas"
# [1] "norway_affy6_old_gwas"
# [1] "belgium_inf1_old_gwas"
# [1] "belgium_inf2_old_gwas"
# [1] "niddk_old_gwas"
# [1] "cedars_370k_old_gwas"
# [1] "cedars_610k_old_gwas"
# [1] "cedars_omni_old_gwas"
# [1] "swedish_uc_old_gwas"
# [1] "mccauley_gsa"
# [1] "N variants to exclude:42531"
# [1] "ccfa_gsa"
# [1] "N variants to exclude:48464"
# [1] "cedars_gsa"
# [1] "N variants to exclude:42766"
# [1] "bernstein_gsa"
# [1] "N variants to exclude:37513"
# [1] "farkkila_gsa"
# [1] "N variants to exclude:35346"
# [1] "franchimont_gsa"
# [1] "N variants to exclude:41070"
# [1] "franke_gsa"
# [1] "N variants to exclude:43272"
# [1] "helmsley_prism_gsa"
# [1] "helmsley_xavier_prism_gsa"
# [1] "hyams_protect_gsa"
# [1] "N variants to exclude:45311"
# [1] "lewis_sparc_gsa"
# [1] "N variants to exclude:48509"
# [1] "mccauley_new_gsa"
# [1] "N variants to exclude:45631"
# [1] "mcgovern_gsa"
# [1] "N variants to exclude:44678"
# [1] "moayyedi_imagine_gsa"
# [1] "N variants to exclude:50194"
# [1] "newberry_share_gsa"
# [1] "N variants to exclude:48021"
# [1] "niddk_cho_gsa"
# [1] "N variants to exclude:45639"
# [1] "niddk_duerr_gsa"
# [1] "N variants to exclude:44781"
# [1] "niddk_rioux_gsa"
# [1] "N variants to exclude:43959"
# [1] "niddk_silverberg_gsa"
# [1] "N variants to exclude:46073"
# [1] "palotie_hus_gsa"
# [1] "N variants to exclude:44139"
# [1] "pekow_share_gsa"
# [1] "N variants to exclude:48644"
# [1] "rioux_igenomed_gsa"
# [1] "N variants to exclude:44546"
# [1] "sands_msccr_gsa"
# [1] "N variants to exclude:49597"
# [1] "stampfer_gsa"
# [1] "N variants to exclude:39214"
# [1] "vermeire_gsa"
# [1] "N variants to exclude:41534"
# [1] "weersma_gsa"
# [1] "N variants to exclude:47336"
# [1] "xavier_prism_gsa"
# [1] "N variants to exclude:43275"
# [1] "xavier_share_gsa"
# [1] "N variants to exclude:42306"
#########

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500

for i in ${studies[@]}
do
bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_exclude_gsa_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_exclude_gsa_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck \
--keep-allele-order --allow-no-sex \
--exclude ${path_gwas}pre_imputation/QC/${i}/list_gsa_variants_flagged_to_exclude \
--make-bed --out ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_compare_freq_3_${i} | grep -E "completed"
done

for i in ${studies[@]}
do echo ${i} && wc -l ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
done

#######
# all_hce
# 402621 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/all_hce/all_hce_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# niddk_old_gwas
# 296707 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_old_gwas/niddk_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# australia_omniexome
# 726852 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/australia_omniexome/australia_omniexome_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# gwas1
# 437006 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/gwas1/gwas1_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# gwas2
# 757166 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/gwas2/gwas2_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# pittsburgh_gsa
# 873299 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/pittsburgh_gsa/pittsburgh_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# spain_gsa
# 567253 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/spain_gsa/spain_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# italy_gsa
# 502511 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/italy_gsa/italy_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# kiel_austria_sibdcs_gsa
# 507134 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# netherlands_gsa
# 520390 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/netherlands_gsa/netherlands_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# slovenia_gsa
# 472856 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/slovenia_gsa/slovenia_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# sweden_gsa
# 495600 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/sweden_gsa/sweden_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# niddk_broad_gsa
# 535756 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# niddk_feinstein_gsa
# 548429 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# basque_gsa
# 503285 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/basque_gsa/basque_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# lithuania_gsa
# 496986 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/lithuania_gsa/lithuania_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# belgium_louis_gsa
# 498223 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_louis_gsa/belgium_louis_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# belgium_franchimont_gsa
# 495548 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_franchimont_gsa/belgium_franchimont_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# belgium_vermeire_gsa
# 558983 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_vermeire_gsa/belgium_vermeire_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# prism_nfe_gsa
# 488279 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# prism_nfe_gwas
# 239737 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/prism_nfe_gwas/prism_nfe_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# finland_illugwas
# 231862 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/finland_illugwas/finland_illugwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# german_affy6_old_gwas
# 766440 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/german_affy6_old_gwas/german_affy6_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# norway_affy6_old_gwas
# 708069 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# belgium_inf1_old_gwas
# 298229 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# belgium_inf2_old_gwas
# 286780 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# cedars_370k_old_gwas
# 333600 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# cedars_610k_old_gwas
# 574805 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# cedars_omni_old_gwas
# 603028 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# swedish_uc_old_gwas
# 293205 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# mccauley_gsa
# 494663 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mccauley_gsa/mccauley_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# ccfa_gsa
# 529246 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/ccfa_gsa/ccfa_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# cedars_gsa
# 543255 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/cedars_gsa/cedars_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# bernstein_gsa
# 423775 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/bernstein_gsa/bernstein_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# farkkila_gsa
# 435690 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/farkkila_gsa/farkkila_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# franchimont_gsa
# 510770 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/franchimont_gsa/franchimont_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# franke_gsa
# 495087 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/franke_gsa/franke_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# helmsley_prism_gsa
# 239932 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/helmsley_prism_gsa/helmsley_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# helmsley_xavier_prism_gsa
# 239968 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/helmsley_xavier_prism_gsa/helmsley_xavier_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# hyams_protect_gsa
# 504126 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/hyams_protect_gsa/hyams_protect_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# lewis_sparc_gsa
# 535128 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/lewis_sparc_gsa/lewis_sparc_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# mccauley_new_gsa
# 511576 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mccauley_new_gsa/mccauley_new_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# mcgovern_gsa
# 555799 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# moayyedi_imagine_gsa
# 518613 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/moayyedi_imagine_gsa/moayyedi_imagine_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# newberry_share_gsa
# 508084 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/newberry_share_gsa/newberry_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# niddk_cho_gsa
# 525238 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_cho_gsa/niddk_cho_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# niddk_duerr_gsa
# 514622 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_duerr_gsa/niddk_duerr_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# niddk_rioux_gsa
# 505968 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_rioux_gsa/niddk_rioux_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# niddk_silverberg_gsa
# 531796 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/niddk_silverberg_gsa/niddk_silverberg_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# palotie_hus_gsa
# 487940 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/palotie_hus_gsa/palotie_hus_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# pekow_share_gsa
# 504572 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/pekow_share_gsa/pekow_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# rioux_igenomed_gsa
# 480503 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# sands_msccr_gsa
# 524394 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/sands_msccr_gsa/sands_msccr_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# stampfer_gsa
# 509287 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/stampfer_gsa/stampfer_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# vermeire_gsa
# 520566 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/vermeire_gsa/vermeire_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# weersma_gsa
# 494166 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/weersma_gsa/weersma_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# xavier_prism_gsa
# 504358 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/xavier_prism_gsa/xavier_prism_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
# xavier_share_gsa
# 505923 /path/to/ibdgwas/IIBDGC//pre_imputation/QC/xavier_share_gsa/xavier_share_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim
#######


################################################################################################################################################################
