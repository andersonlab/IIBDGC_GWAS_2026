# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

#############################################################################
# 11.- REMOVE VARIANTS WITH SIGNIF DISCREPANCY IN CALL RATE BETWEEN STUDIES #
#############################################################################

###### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"
studies<-c("all_hce","australia_omniexome","kiel_austria_sibdcs_gsa")

for (j in 1:length(studies)) {
  
  print(studies[j])
  
  fam<-read.table(paste(path,
                        "pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam",sep=""),head=F)
  
  if(studies[j]=="all_hce") {
    cohorts<-c("gwas3","new_wave")
  } else if (studies[j]=="australia_omniexome") {
    cohorts<-c("australia_2012_omniexome","australia_2014_omniexome")
  } else if (studies[j]=="kiel_austria_sibdcs_gsa") {
    cohorts<-c("kiel_foc_gsa","kiel_eze_gsa","kiel_bc_gsa","kiel_ibd_gsa"
               ,"kiel_hlitm_gsa","austria_gsa","sibdcs_gsa")
  }
  
  for (i in 1:length(cohorts)){
    
    file_fam<-paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19.fam",sep="")
    
    if (studies[j]=="all_hce" & cohorts[i]=="gwas3") {
      file_fam<-paste(path,"pre_imputation/QC/",cohorts[i],"/",cohorts[i],"_hg19_tmp.fam",sep="")
    }
    
    if (studies[j]=="kiel_austria_sibdcs_gsa") {
      file_fam<-paste(path,"pre_imputation/QC/",studies[j],"/",cohorts[i],"_hg19.fam",sep="")
    }
    
    if (i==1) {
      a<-read.table(file_fam,head=F)
      a<-a[which(a$V1 %in% fam$V1),]
      a$study<-cohorts[i]
      
    }else{
      b<-read.table(file_fam,head=F)
      b<-b[which(b$V1 %in% fam$V1),]
      b$study<-cohorts[i]
      
      a<-rbind(a,b)
    }
  }
  table(as.character(a$V1)==as.character(fam$V1))
  
  a<-a[match(fam$V1,a$V1),]
  print(table(as.character(a$V1)==as.character(fam$V1)))
  
  for (i in 1:length(cohorts)) {
    tmp<-a[which(a$study==cohorts[i]),c(1,1),drop=F]
    colnames(tmp)<-c("FID","IID")
    write.table(tmp,paste(path,"pre_imputation/QC/",studies[j],"/list_ids_",cohorts[i],sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  }
}

# [1] "all_hce"
# TRUE 
# 23618 

# [1] "australia_omniexome"
# TRUE 
# 1306 

# [1] "kiel_austria_sibdcs_gsa"
# TRUE 
# 14403 


########################################

studies=(all_hce)
cohorts=(gwas3 new_wave)

#####

studies=(australia_omniexome)
cohorts=(australia_2012_omniexome australia_2014_omniexome)

#####

studies=(kiel_austria_sibdcs_gsa)
cohorts=(austria_gsa kiel_bc_gsa kiel_eze_gsa kiel_foc_gsa kiel_hlitm_gsa kiel_ibd_gsa sibdcs_gsa)

for i in ${studies[@]}
do wc -l ${path_gwas}pre_imputation/QC/${i}/list_ids_*
done

for i in ${studies[@]}
do for j in ${cohorts[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep \
--allow-no-sex \
--keep ${path_gwas}pre_imputation/QC/${i}/list_ids_${j} \
--missing \
--out ${path_gwas}pre_imputation/QC/${i}/tmp_${j}
done
done


### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"
studies<-c("all_hce","australia_omniexome","kiel_austria_sibdcs_gsa")

# distribution of missingness by array

for (j in 1:length(studies)) {
  
  print(studies[j])
  
  fam<-read.table(paste(path,
                        "pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam",sep=""),head=F)
  
  if(studies[j]=="all_hce") {
    cohorts<-c("gwas3","new_wave")
  } else if (studies[j]=="australia_omniexome") {
    cohorts<-c("australia_2012_omniexome","australia_2014_omniexome")
  } else if (studies[j]=="kiel_austria_sibdcs_gsa") {
    cohorts<-c("kiel_foc_gsa","kiel_eze_gsa","kiel_bc_gsa","kiel_ibd_gsa"
               ,"kiel_hlitm_gsa","austria_gsa","sibdcs_gsa")
  }
  
  for (i in 1:length(cohorts)){
    
    if(i==1){
      
      tmp<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/tmp_",cohorts[i],".lmiss",sep=""),head=T)
      a1<-as.data.frame(table(cut(tmp$F_MISS,breaks=c(-1,0.05,0.1,0.2,0.3,0.4),labels=c("0-0.05","0.05-0.1","0.1-0.2","0.2-0.3","0.3-0.4"))) )
      colnames(a1)<-c("numbers",cohorts[i])
      
    }else{
      
      tmp<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/tmp_",cohorts[i],".lmiss",sep=""),head=T)
      a2<-as.data.frame(table(cut(tmp$F_MISS,breaks=c(-1,0.05,0.1,0.2,0.3,0.4),labels=c("0-0.05","0.05-0.1","0.1-0.2","0.2-0.3","0.3-0.4"))) )
      colnames(a2)<-c("numbers",cohorts[i])
      
      a1<-merge(a1,a2,by="numbers",sort=F)
      
    }
    
  }
  
  print(a1)
  
  write.table(a1,paste(path,"pre_imputation/QC/",studies[j],"/table_breakdown_number_missing_variants_perstudy.csv",sep=""),col.names=T,row.names=F,quote=F,sep=",")
  
  
  
  for (i in 1:length(cohorts)){
    
    if (i==1) {
      
      tmp<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/tmp_",cohorts[i],".lmiss",sep=""),head=T)
      tmp$study<-cohorts[i]
      tmp<-tmp[,c("SNP","F_MISS","study")]
      
      all<-tmp
      
    } else {
      
      tmp<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/tmp_",cohorts[i],".lmiss",sep=""),head=T)
      tmp$study<-cohorts[i]
      tmp<-tmp[,c("SNP","F_MISS","study")]
      
      all<-rbind(all,tmp)
      
    }
    
  }
  
  print(table(all$study))
  
  variants<-all[which(all$F_MISS>0.1),"SNP",drop=F]
  variants<-variants[!duplicated(variants$SNP),,drop=F]
  print(nrow(variants))
  # [1] 14  1
  
  bim<-read.table(paste(path,
                        "pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim",sep=""),head=F)
  
  # table(bim$V1[which(bim$V2 %in% variants$SNP)])
  # dim(bim[which(bim$V2 %in% variants$SNP),])
  
  write.table(variants,paste(path,"pre_imputation/QC/",studies[j],"/list_variants_per_study_missing_0.1",sep=""),col.names=F,row.names=F,quote=F)

}


###############

# [1] "all_hce"
# numbers  gwas3 new_wave
# 1   0-0.05 436289   436244
# 2 0.05-0.1      3       36
# 3  0.1-0.2      0       10
# 4  0.2-0.3      0        2
# 5  0.3-0.4      0        0
# 
# gwas3 new_wave 
# 436292   436292 
# [1] 12


# [1] "australia_omniexome"
# numbers australia_2012_omniexome australia_2014_omniexome
# 1   0-0.05                   739114                   739812
# 2 0.05-0.1                      912                      212
# 3  0.1-0.2                        0                        2
# 4  0.2-0.3                        0                        0
# 5  0.3-0.4                        0                        0
# 
# australia_2012_omniexome australia_2014_omniexome 
# 740026                   740026 
# [1] 2


# [1] "kiel_austria_sibdcs_gsa"
# numbers kiel_foc_gsa kiel_eze_gsa kiel_bc_gsa kiel_ibd_gsa kiel_hlitm_gsa
# 1   0-0.05       564564       558745      558758       554635         563587
# 2 0.05-0.1         1530         7253        7339        11464           2512
# 3  0.1-0.2            5          101           2            0              0
# 4  0.2-0.3            0            0           0            0              0
# 5  0.3-0.4            0            0           0            0              0
# austria_gsa sibdcs_gsa
# 1      566081     566004
# 2          18         95
# 3           0          0
# 4           0          0
# 5           0          0
# 
# austria_gsa    kiel_bc_gsa   kiel_eze_gsa   kiel_foc_gsa kiel_hlitm_gsa 
# 566099         566099         566099         566099         566099 
# kiel_ibd_gsa     sibdcs_gsa 
# 566099         566099 
# [1] 106

###############

# keep FID and IID with same code:

###### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"
studies<-c("all_hce","niddk_old_gwas","australia_omniexome","gwas1","gwas2","pittsburgh_gsa","spain_gsa","italy_gsa",
           "kiel_austria_sibdcs_gsa","netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa",
           "niddk_feinstein_gsa","basque_gsa","lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa",
           "belgium_vermeire_gsa","prism_nfe_gsa","prism_nfe_gwas","finland_illugwas","german_affy6_old_gwas",
           "norway_affy6_old_gwas","belgium_inf1_old_gwas","belgium_inf2_old_gwas",
           "cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas",
           "swedish_uc_old_gwas",
           "mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

for (j in 1:length(studies)) {
  
  print(studies[j])
  
  fam<-read.table(paste(path,
                        "pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.fam",sep=""),head=F)
  print(table(fam$V1==fam$V2))
  
  if(studies[j]=="belgium_inf2_old_gwas") {
    fam$V1<-fam$V2
  }else {
    fam$V2<-fam$V1
  }
  
  write.table(fam,paste(path,
                        "pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_edited.fam",sep=""),
              col.names=F,row.names=F,quote=F,sep="\t")
}

###############

path_gwas=/path/to/ibdgwas/IIBDGC/
studies=(niddk_old_gwas gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_32_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_32_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bed \
--bim ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim \
--fam ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_edited.fam \
--allow-no-sex \
--make-bed \
--missing \
--out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy"
done


studies=(all_hce australia_omniexome kiel_austria_sibdcs_gsa)

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_32_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_32_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bed \
--bim ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep.bim \
--fam ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_edited.fam \
--allow-no-sex \
--exclude ${path_gwas}pre_imputation/QC/${i}/list_variants_per_study_missing_0.1 \
--make-bed \
--missing \
--out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy"
done

for i in ${studies[@]}
do echo ${i} && wc -l ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy.bim
done
# all_hce
# 436280 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/all_hce/all_hce_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy.bim
# australia_omniexome
# 740024 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/australia_omniexome/australia_omniexome_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy.bim
# kiel_austria_sibdcs_gsa
# 565993 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy.bim


studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_32_${i} | grep -E "completed"
done

###########
# all_hce
# Successfully completed.
# niddk_old_gwas
# Successfully completed.
# australia_omniexome
# Successfully completed.
# gwas1
# Successfully completed.
# gwas2
# Successfully completed.
# pittsburgh_gsa
# Successfully completed.
# spain_gsa
# Successfully completed.
# italy_gsa
# Successfully completed.
# kiel_austria_sibdcs_gsa
# Successfully completed.
# netherlands_gsa
# Successfully completed.
# slovenia_gsa
# Successfully completed.
# sweden_gsa
# Successfully completed.
# niddk_broad_gsa
# Successfully completed.
# niddk_feinstein_gsa
# Successfully completed.
# basque_gsa
# Successfully completed.
# lithuania_gsa
# Successfully completed.
# belgium_louis_gsa
# Successfully completed.
# belgium_franchimont_gsa
# Successfully completed.
# belgium_vermeire_gsa
# Successfully completed.
# prism_nfe_gsa
# Successfully completed.
# prism_nfe_gwas
# Successfully completed.
# finland_illugwas
# Successfully completed.
# german_affy6_old_gwas
# Successfully completed.
# norway_affy6_old_gwas
# Successfully completed.
# belgium_inf1_old_gwas
# Successfully completed.
# belgium_inf2_old_gwas
# Successfully completed.
# cedars_370k_old_gwas
# Successfully completed.
# cedars_610k_old_gwas
# Successfully completed.
# cedars_omni_old_gwas
# Successfully completed.
# swedish_uc_old_gwas
# Successfully completed.
# mccauley_gsa
# Successfully completed.
# ccfa_gsa
# Successfully completed.
# cedars_gsa
# Successfully completed.
# bernstein_gsa
# Successfully completed.
# farkkila_gsa
# Successfully completed.
# franchimont_gsa
# Successfully completed.
# franke_gsa
# Successfully completed.
# helmsley_prism_gsa
# Successfully completed.
# helmsley_xavier_prism_gsa
# Successfully completed.
# hyams_protect_gsa
# Successfully completed.
# lewis_sparc_gsa
# Successfully completed.
# mccauley_new_gsa
# Successfully completed.
# mcgovern_gsa
# Successfully completed.
# moayyedi_imagine_gsa
# Successfully completed.
# newberry_share_gsa
# Successfully completed.
# niddk_cho_gsa
# Successfully completed.
# niddk_duerr_gsa
# Successfully completed.
# niddk_rioux_gsa
# Successfully completed.
# niddk_silverberg_gsa
# Successfully completed.
# palotie_hus_gsa
# Successfully completed.
# pekow_share_gsa
# Successfully completed.
# rioux_igenomed_gsa
# Successfully completed.
# sands_msccr_gsa
# Successfully completed.
# stampfer_gsa
# Successfully completed.
# vermeire_gsa
# Successfully completed.
# weersma_gsa
# Successfully completed.
# xavier_prism_gsa
# Successfully completed.
# xavier_share_gsa
# Successfully completed.
###########


############################################################################################################################################
