# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

##############################################
# 6.- COMPARE ALLELE FREQUENCIES WITH 1000GP #
##############################################

#########################################################################
# 6.1 COMPARE ALLELE FREQUENCIES WITH 1000GP, GENERATE LIST OF VARIANTS:

path_gwas=/path/to/ibdgwas/IIBDGC/
  
studies=(australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)
  
for i in ${studies[@]}
do
cat ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup.bim | cut -f 2 > ${path_gwas}/pre_imputation/QC/${i}/list_variants_${i}_hg19_noind_posstr_nodup
done

for i in ${studies[@]}
do
echo ${i} && ls -la ${path_gwas}/pre_imputation/QC/${i}/list_variants_${i}_hg19_noind_posstr_nodup
done

##############################################
# 6.2 EXTRACT VARIANTS FROM 1000GP

# see script 1000gp_data_to_compare.R about how 1000GP data was generated

path_gwas=/path/to/ibdgwas/IIBDGC/
studies=(australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

MEM=5000

for i in ${studies[@]}
do
for chr in {1..24}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_6_${i}_${chr} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_6_${i}_${chr} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/1000gp/1000GP_EUR_chr${chr}_b37 \
--extract ${path_gwas}/pre_imputation/QC/${i}/list_variants_${i}_hg19_noind_posstr_nodup \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/1000GP_EUR_chr${chr}_b37_${i}_variants"
done
done
# Job <810587..811190> is submitted to queue <normal>.

for i in ${studies[@]}
do
echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/1000GP_EUR_chr*_b37_${i}_variants.bim | wc -l
done

##############
# NO CHRY DATA
# gwas1
# 23
# prism_nfe_gwas
# 23
# finland_illugwas
# 23
# belgium_inf1_old_gwas
# 23
# belgium_inf2_old_gwas
# 23
# helmsley_prism_gsa
# 23
# helmsley_xavier_prism_gsa
# 23
# NO CHRY OR CHRX DATA
# swedish_uc_old_gwas
# 22

# rest chrx and chrY
##############

##############################################
# 6.3 COMBINE FILES

##############################
### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

# "niddk_old_gwas","all_hce" - not included
path<-"/path/to/ibdgwas/IIBDGC/"

studies<-c("australia_omniexome","gwas1","gwas2","pittsburgh_gsa","spain_gsa","italy_gsa",
          "kiel_austria_sibdcs_gsa","netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa",
          "niddk_feinstein_gsa","basque_gsa","lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa",
          "belgium_vermeire_gsa","prism_nfe_gsa","prism_nfe_gwas","finland_illugwas","german_affy6_old_gwas",
          "norway_affy6_old_gwas","belgium_inf1_old_gwas","belgium_inf2_old_gwas",
          "cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","swedish_uc_old_gwas",
          "mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
          "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
          "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
          "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
          "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
          "xavier_share_gsa")

non_chrxy<-c("swedish_uc_old_gwas")
non_chry<-c("gwas1","prism_nfe_gwas","finland_illugwas","belgium_inf1_old_gwas","belgium_inf2_old_gwas",
"helmsley_prism_gsa","helmsley_xavier_prism_gsa")

for (j in 1:length(studies)) {
  
  if(studies[j] %in% non_chry){
    N_rows<-22
  }else if(studies[j] %in% non_chrxy) {
    N_rows<-21
  }else {
    N_rows<-23
  }
  
  dat<-matrix(ncol=3,nrow=N_rows)
  dat<-as.data.frame(dat)
  
  for (i in 1:nrow(dat)){
    dat[i,1]<-paste(path,"pre_imputation/QC/",studies[j],"/1000GP_EUR_chr",i+1,"_b37_",studies[j],"_variants.bed",sep="")
    dat[i,2]<-paste(path,"pre_imputation/QC/",studies[j],"/1000GP_EUR_chr",i+1,"_b37_",studies[j],"_variants.bim",sep="")
    dat[i,3]<-paste(path,"pre_imputation/QC/",studies[j],"/1000GP_EUR_chr",i+1,"_b37_",studies[j],"_variants.fam",sep="")
  }
  
  write.table(dat,paste(path,"pre_imputation/QC/",studies[j],"/list_1000GP_files_",studies[j],"_variants_tomerge.txt",sep=""),
              col.names=F,row.names=F,quote=F,sep="\t")
  
}

q("no")

##############################

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=2500

for i in ${studies[@]}
do
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_7_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_7_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/1000GP_EUR_chr1_b37_${i}_variants \
--merge-list ${path_gwas}pre_imputation/QC/${i}/list_1000GP_files_${i}_variants_tomerge.txt \
--allow-no-sex \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/1000GP_EUR_b37_${i}_variants"
done

for i in ${studies[@]}
do echo ${i} && tail -100 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_7_${i} | grep -E "completed|people pass"
done

##############################
# australia_omniexome
# 767234 variants and 503 people pass filters and QC.
# Successfully completed.
# gwas1
# 457782 variants and 503 people pass filters and QC.
# Successfully completed.
# gwas2
# 792837 variants and 503 people pass filters and QC.
# Successfully completed.
# pittsburgh_gsa
# 936134 variants and 503 people pass filters and QC.
# Successfully completed.
# spain_gsa
# 583256 variants and 503 people pass filters and QC.
# Successfully completed.
# italy_gsa
# 578601 variants and 503 people pass filters and QC.
# Successfully completed.
# kiel_austria_sibdcs_gsa
# 643055 variants and 503 people pass filters and QC.
# Successfully completed.
# netherlands_gsa
# 657510 variants and 503 people pass filters and QC.
# Successfully completed.
# slovenia_gsa
# 554910 variants and 503 people pass filters and QC.
# Successfully completed.
# sweden_gsa
# 572206 variants and 503 people pass filters and QC.
# Successfully completed.
# niddk_broad_gsa
# 615081 variants and 503 people pass filters and QC.
# Successfully completed.
# niddk_feinstein_gsa
# 631736 variants and 503 people pass filters and QC.
# Successfully completed.
# basque_gsa
# 573074 variants and 503 people pass filters and QC.
# Successfully completed.
# lithuania_gsa
# 576182 variants and 503 people pass filters and QC.
# Successfully completed.
# belgium_louis_gsa
# 595810 variants and 503 people pass filters and QC.
# Successfully completed.
# belgium_franchimont_gsa
# 580165 variants and 503 people pass filters and QC.
# Successfully completed.
# belgium_vermeire_gsa
# 648977 variants and 503 people pass filters and QC.
# Successfully completed.
# prism_nfe_gsa
# 570948 variants and 503 people pass filters and QC.
# Successfully completed.
# prism_nfe_gwas
# 243339 variants and 503 people pass filters and QC.
# Successfully completed.
# finland_illugwas
# 238611 variants and 503 people pass filters and QC.
# Successfully completed.
# german_affy6_old_gwas
# 886572 variants and 503 people pass filters and QC.
# Successfully completed.
# norway_affy6_old_gwas
# 842682 variants and 503 people pass filters and QC.
# Successfully completed.
# belgium_inf1_old_gwas
# 301706 variants and 503 people pass filters and QC.
# Successfully completed.
# belgium_inf2_old_gwas
# 289804 variants and 503 people pass filters and QC.
# Successfully completed.
# cedars_370k_old_gwas
# 338482 variants and 503 people pass filters and QC.
# Successfully completed.
# cedars_610k_old_gwas
# 583912 variants and 503 people pass filters and QC.
# Successfully completed.
# cedars_omni_old_gwas
# 612985 variants and 503 people pass filters and QC.
# Successfully completed.
# swedish_uc_old_gwas
# 296077 variants and 503 people pass filters and QC.
# Successfully completed.
# mccauley_gsa
# 567549 variants and 503 people pass filters and QC.
# Successfully completed.
# ccfa_gsa
# 605583 variants and 503 people pass filters and QC.
# Successfully completed.
# cedars_gsa
# 647961 variants and 503 people pass filters and QC.
# Successfully completed.
# bernstein_gsa
# 563551 variants and 503 people pass filters and QC.
# Successfully completed.
# farkkila_gsa
# 527169 variants and 503 people pass filters and QC.
# Successfully completed.
# franchimont_gsa
# 646368 variants and 503 people pass filters and QC.
# Successfully completed.
# franke_gsa
# 568298 variants and 503 people pass filters and QC.
# Successfully completed.
# helmsley_prism_gsa
# 243339 variants and 503 people pass filters and QC.
# Successfully completed.
# helmsley_xavier_prism_gsa
# 243309 variants and 503 people pass filters and QC.
# Successfully completed.
# hyams_protect_gsa
# 577956 variants and 503 people pass filters and QC.
# Successfully completed.
# lewis_sparc_gsa
# 611251 variants and 503 people pass filters and QC.
# Successfully completed.
# mccauley_new_gsa
# 586270 variants and 503 people pass filters and QC.
# Successfully completed.
# mcgovern_gsa
# 639753 variants and 503 people pass filters and QC.
# Successfully completed.
# moayyedi_imagine_gsa
# 609263 variants and 503 people pass filters and QC.
# Successfully completed.
# newberry_share_gsa
# 588149 variants and 503 people pass filters and QC.
# Successfully completed.
# niddk_cho_gsa
# 602279 variants and 503 people pass filters and QC.
# Successfully completed.
# niddk_duerr_gsa
# 589769 variants and 503 people pass filters and QC.
# Successfully completed.
# niddk_rioux_gsa
# 584504 variants and 503 people pass filters and QC.
# Successfully completed.
# niddk_silverberg_gsa
# 618032 variants and 503 people pass filters and QC.
# Successfully completed.
# palotie_hus_gsa
# 562617 variants and 503 people pass filters and QC.
# Successfully completed.
# pekow_share_gsa
# 578515 variants and 503 people pass filters and QC.
# Successfully completed.
# rioux_igenomed_gsa
# 550181 variants and 503 people pass filters and QC.
# Successfully completed.
# sands_msccr_gsa
# 600363 variants and 503 people pass filters and QC.
# Successfully completed.
# stampfer_gsa
# 647980 variants and 503 people pass filters and QC.
# Successfully completed.
# vermeire_gsa
# 649533 variants and 503 people pass filters and QC.
# Successfully completed.
# weersma_gsa
# 571212 variants and 503 people pass filters and QC.
# Successfully completed.
# xavier_prism_gsa
# 580947 variants and 503 people pass filters and QC.
# Successfully completed.
# xavier_share_gsa
# 580888 variants and 503 people pass filters and QC.
# Successfully completed.
##############################

# remove intermediate files:
for i in ${studies[@]}
do 
rm ${path_gwas}pre_imputation/QC/${i}/1000GP_EUR_chr*_b37_${i}_variants*
done
  
##############################################
# 6.4 ESTIMATE FREQUENCIES:

# 1000GP

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_8_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_8_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/1000GP_EUR_b37_${i}_variants \
--freq \
--allow-no-sex \
--out ${path_gwas}pre_imputation/QC/${i}/1000GP_EUR_b37_${i}_variants"
done

for i in ${studies[@]}
do echo ${i} && tail -100 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_8_${i} | grep -E "completed"
done

# Studies

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_9_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_9_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup \
--freq \
--allow-no-sex \
--out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup"
done

for i in ${studies[@]}
do echo ${i} && tail -100 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_9_${i} | grep -E "completed"
done

for i in ${studies[@]}
do echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup.frq \
&& ls -la ${path_gwas}pre_imputation/QC/${i}/1000GP_EUR_b37_${i}_variants.frq
done

### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)


# "niddk_old_gwas","all_hce" - not included
path<-"/path/to/ibdgwas/IIBDGC/"

studies<-c("australia_omniexome","gwas1","gwas2","pittsburgh_gsa","spain_gsa","italy_gsa",
           "kiel_austria_sibdcs_gsa","netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa",
           "niddk_feinstein_gsa","basque_gsa","lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa",
           "belgium_vermeire_gsa","prism_nfe_gsa","prism_nfe_gwas","finland_illugwas","german_affy6_old_gwas",
           "norway_affy6_old_gwas","belgium_inf1_old_gwas","belgium_inf2_old_gwas",
           "cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","swedish_uc_old_gwas",
           "mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

for (j in 1:length(studies)) {
  
  print(studies[j])
  gp<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/1000GP_EUR_b37_",studies[j],"_variants.frq",sep=""),head=T)
  g1<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup.frq",sep=""),head=T)
  
  print(paste("N variants in study:",nrow(g1)))
  print(paste("N variants from study in 1000GP:",nrow(gp)))
  #
  
  colnames(gp)[3:6]<-paste(colnames(gp)[3:6],"_gp",sep="")
  colnames(g1)[3:6]<-paste(colnames(g1)[3:6],"_g1",sep="")
  
  all<-merge(g1,gp[,2:6],by="SNP",all=T)
  
  check<-all[which(all$A1_g1!=all$A1_gp),]
  print(paste("N_variants different Ref allele:",nrow(check)))
  
  # keep only A/T C/G
  check<-check[which( (check$A1_g1=="G" & check$A2_g1=="C") | (check$A1_g1=="C" & check$A2_g1=="G") | (check$A1_g1=="A" & check$A2_g1=="T") | (check$A1_g1=="T" & check$A2_g1=="A")),]
  print(paste("N_variants to check:",nrow(check)))
  
  if(nrow(check)>0) {
    
    # LIST OF VARIANTS TO REMOVE, WE CANNOT REALLY BE SURE WHETHER THERE IS STRAND ISSUE OR NOT
    remove<-check[which(check$MAF_g1>=0.45),]
    print(paste("N_variants to remove (maf >=0.45):",nrow(remove)))
    
    # LIST OF VARIANTS TO FLIP:
    flip<-check[which(check$MAF_g1<0.45),]
    
    flip<-flip[order(flip$MAF_g1,decreasing=T),]
    
    pdf(paste(path,"pre_imputation/QC/",studies[j],"/plot_maf_AT_CG_",studies[j],"_1000gp.pdf",sep=""))
    plot(flip$MAF_g1,flip$MAF_gp)
    dev.off()
    
    remove_2<-flip[which(flip$MAF_g1>0.2 & flip$MAF_gp<0.1),]
    remove_3<-flip[which(flip$MAF_g1<0.1 & flip$MAF_gp>0.2),]
      
    remove<-rbind(remove,remove_2,remove_3)
    print(paste("N_variants to remove:",nrow(remove)))
    
    write.table(remove[,"SNP"],paste(path,"pre_imputation/QC/",studies[j],"/list_variants_to_remove_AT_CG",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
    
    flip<-flip[which(!flip$SNP %in% remove$SNP),]
    print(paste("N_variants to flip:",nrow(flip)))
    
    write.table(flip[,"SNP"],paste(path,"pre_imputation/QC/",studies[j],"/list_variants_to_flip_AT_CG",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
    
  }
  
  rm(list=ls()[!ls() %in% c("j","studies","path")])
}
q("no")

##############################

# [1] "australia_omniexome"
# [1] "N variants in study: 788230"
# [1] "N variants from study in 1000GP: 767234"
# [1] "N_variants different Ref allele: 12664"
# [1] "N_variants to check: 75"
# [1] "N_variants to remove (maf >=0.45): 67"
# [1] "N_variants to remove: 68"
# [1] "N_variants to flip: 7"

# [1] "gwas1"
# [1] "N variants in study: 459424"
# [1] "N variants from study in 1000GP: 457782"
# [1] "N_variants different Ref allele: 42771"
# [1] "N_variants to check: 37244"
# [1] "N_variants to remove (maf >=0.45): 2972"
# [1] "N_variants to remove: 2976"
# [1] "N_variants to flip: 34268"

# [1] "gwas2"
# [1] "N variants in study: 796326"
# [1] "N variants from study in 1000GP: 792837"
# [1] "N_variants different Ref allele: 66993"
# [1] "N_variants to check: 55861"
# [1] "N_variants to remove (maf >=0.45): 4742"
# [1] "N_variants to remove: 4772"
# [1] "N_variants to flip: 51089"

# [1] "pittsburgh_gsa"
# [1] "N variants in study: 945262"
# [1] "N variants from study in 1000GP: 936134"
# [1] "N_variants different Ref allele: 12437"
# [1] "N_variants to check: 297"
# [1] "N_variants to remove (maf >=0.45): 279"
# [1] "N_variants to remove: 285"
# [1] "N_variants to flip: 12"
# [1] "spain_gsa"
# [1] "N variants in study: 585728"
# [1] "N variants from study in 1000GP: 583256"
# [1] "N_variants different Ref allele: 12958"
# [1] "N_variants to check: 72"
# [1] "N_variants to remove (maf >=0.45): 66"
# [1] "N_variants to remove: 67"
# [1] "N_variants to flip: 5"

# [1] "italy_gsa"
# [1] "N variants in study: 585029"
# [1] "N variants from study in 1000GP: 578601"
# [1] "N_variants different Ref allele: 10315"
# [1] "N_variants to check: 144"
# [1] "N_variants to remove (maf >=0.45): 81"
# [1] "N_variants to remove: 81"
# [1] "N_variants to flip: 63"

# [1] "kiel_austria_sibdcs_gsa"
# [1] "N variants in study: 662049"
# [1] "N variants from study in 1000GP: 643055"
# [1] "N_variants different Ref allele: 6791"
# [1] "N_variants to check: 94"
# [1] "N_variants to remove (maf >=0.45): 54"
# [1] "N_variants to remove: 65"
# [1] "N_variants to flip: 29"

# [1] "netherlands_gsa"
# [1] "N variants in study: 684335"
# [1] "N variants from study in 1000GP: 657510"
# [1] "N_variants different Ref allele: 9470"
# [1] "N_variants to check: 154"
# [1] "N_variants to remove (maf >=0.45): 59"
# [1] "N_variants to remove: 103"
# [1] "N_variants to flip: 51"

# [1] "slovenia_gsa"
# [1] "N variants in study: 557385"
# [1] "N variants from study in 1000GP: 554910"
# [1] "N_variants different Ref allele: 8987"
# [1] "N_variants to check: 94"
# [1] "N_variants to remove (maf >=0.45): 63"
# [1] "N_variants to remove: 66"
# [1] "N_variants to flip: 28"

# [1] "sweden_gsa" - updated N cases ctr
# [1] "N variants in study: 574884"
# [1] "N variants from study in 1000GP: 572206"
# [1] "N_variants different Ref allele: 7642"
# [1] "N_variants to check: 91"
# [1] "N_variants to remove (maf >=0.45): 54"
# [1] "N_variants to remove: 60"
# [1] "N_variants to flip: 31"

# [1] "niddk_broad_gsa"
# [1] "N variants in study: 621811"
# [1] "N variants from study in 1000GP: 615081"
# [1] "N_variants different Ref allele: 6978"
# [1] "N_variants to check: 86"
# [1] "N_variants to remove (maf >=0.45): 52"
# [1] "N_variants to remove: 62"
# [1] "N_variants to flip: 24"

# [1] "niddk_feinstein_gsa"
# [1] "N variants in study: 644846"
# [1] "N variants from study in 1000GP: 631736"
# [1] "N_variants different Ref allele: 8749"
# [1] "N_variants to check: 108"
# [1] "N_variants to remove (maf >=0.45): 71"
# [1] "N_variants to remove: 78"
# [1] "N_variants to flip: 30"

# [1] "basque_gsa"
# [1] "N variants in study: 578033"
# [1] "N variants from study in 1000GP: 573074"
# [1] "N_variants different Ref allele: 8532"
# [1] "N_variants to check: 83"
# [1] "N_variants to remove (maf >=0.45): 69"
# [1] "N_variants to remove: 69"
# [1] "N_variants to flip: 14"

# [1] "lithuania_gsa"
# [1] "N variants in study: 579791"
# [1] "N variants from study in 1000GP: 576182"
# [1] "N_variants different Ref allele: 14689"
# [1] "N_variants to check: 1976"
# [1] "N_variants to remove (maf >=0.45): 117"
# [1] "N_variants to remove: 146"
# [1] "N_variants to flip: 1830"

# [1] "belgium_louis_gsa" - change in N
# [1] "N variants in study: 605456"
# [1] "N variants from study in 1000GP: 595810"
# [1] "N_variants different Ref allele: 11438"
# [1] "N_variants to check: 2328"
# [1] "N_variants to remove (maf >=0.45): 112"
# [1] "N_variants to remove: 147"
# [1] "N_variants to flip: 2181"

# [1] "belgium_franchimont_gsa"
# [1] "N variants in study: 584944"
# [1] "N variants from study in 1000GP: 580165"
# [1] "N_variants different Ref allele: 8907"
# [1] "N_variants to check: 2252"
# [1] "N_variants to remove (maf >=0.45): 122"
# [1] "N_variants to remove: 146"
# [1] "N_variants to flip: 2106"

# [1] "belgium_vermeire_gsa"
# [1] "N variants in study: 674335"
# [1] "N variants from study in 1000GP: 648977"
# [1] "N_variants different Ref allele: 10919"
# [1] "N_variants to check: 3385"
# [1] "N_variants to remove (maf >=0.45): 110"
# [1] "N_variants to remove: 175"
# [1] "N_variants to flip: 3210"

# [1] "prism_nfe_gsa"
# [1] "N variants in study: 573465"
# [1] "N variants from study in 1000GP: 570948"
# [1] "N_variants different Ref allele: 9532"
# [1] "N_variants to check: 1999"
# [1] "N_variants to remove (maf >=0.45): 106"
# [1] "N_variants to remove: 125"
# [1] "N_variants to flip: 1874"

# [1] "prism_nfe_gwas"
# [1] "N variants in study: 243983"
# [1] "N variants from study in 1000GP: 243339"
# [1] "N_variants different Ref allele: 6124"
# [1] "N_variants to check: 146"
# [1] "N_variants to remove (maf >=0.45): 15"
# [1] "N_variants to remove: 15"
# [1] "N_variants to flip: 131"

# [1] "finland_illugwas"
# [1] "N variants in study: 239244"
# [1] "N variants from study in 1000GP: 238611"
# [1] "N_variants different Ref allele: 11632"
# [1] "N_variants to check: 146"
# [1] "N_variants to remove (maf >=0.45): 15"
# [1] "N_variants to remove: 15"
# [1] "N_variants to flip: 131"

# [1] "german_affy6_old_gwas"
# [1] "N variants in study: 890363"
# [1] "N variants from study in 1000GP: 886572"
# [1] "N_variants different Ref allele: 79513"
# [1] "N_variants to check: 68290"
# [1] "N_variants to remove (maf >=0.45): 5114"
# [1] "N_variants to remove: 5165"
# [1] "N_variants to flip: 63125"
# [1] "norway_affy6_old_gwas"
# [1] "N variants in study: 846267"
# [1] "N variants from study in 1000GP: 842682"
# [1] "N_variants different Ref allele: 81497"
# [1] "N_variants to check: 64964"
# [1] "N_variants to remove (maf >=0.45): 5073"
# [1] "N_variants to remove: 5261"
# [1] "N_variants to flip: 59703"
# [1] "belgium_inf1_old_gwas"
# [1] "N variants in study: 302700"
# [1] "N variants from study in 1000GP: 301706"
# [1] "N_variants different Ref allele: 5253"
# [1] "N_variants to check: 0"
# [1] "belgium_inf2_old_gwas"
# [1] "N variants in study: 290730"
# [1] "N variants from study in 1000GP: 289804"
# [1] "N_variants different Ref allele: 6590"
# [1] "N_variants to check: 0"
# [1] "cedars_370k_old_gwas"
# [1] "N variants in study: 339974"
# [1] "N variants from study in 1000GP: 338482"
# [1] "N_variants different Ref allele: 9638"
# [1] "N_variants to check: 36"
# [1] "N_variants to remove (maf >=0.45): 28"
# [1] "N_variants to remove: 29"
# [1] "N_variants to flip: 7"
# [1] "cedars_610k_old_gwas"
# [1] "N variants in study: 586374"
# [1] "N variants from study in 1000GP: 583912"
# [1] "N_variants different Ref allele: 13553"
# [1] "N_variants to check: 87"
# [1] "N_variants to remove (maf >=0.45): 80"
# [1] "N_variants to remove: 81"
# [1] "N_variants to flip: 6"
# [1] "cedars_omni_old_gwas"
# [1] "N variants in study: 721651"
# [1] "N variants from study in 1000GP: 612985"
# [1] "N_variants different Ref allele: 14655"
# [1] "N_variants to check: 12"
# [1] "N_variants to remove (maf >=0.45): 1"
# [1] "N_variants to remove: 5"
# [1] "N_variants to flip: 7"
# [1] "swedish_uc_old_gwas"
# [1] "N variants in study: 297004"
# [1] "N variants from study in 1000GP: 296077"
# [1] "N_variants different Ref allele: 7614"
# [1] "N_variants to check: 0"
# [1] "mccauley_gsa"
# [1] "N variants in study: 570878"
# [1] "N variants from study in 1000GP: 567549"
# [1] "N_variants different Ref allele: 10857"
# [1] "N_variants to check: 2063"
# [1] "N_variants to remove (maf >=0.45): 107"
# [1] "N_variants to remove: 123"
# [1] "N_variants to flip: 1940"
# [1] "ccfa_gsa"
# [1] "N variants in study: 611031"
# [1] "N variants from study in 1000GP: 605583"
# [1] "N_variants different Ref allele: 9721"
# [1] "N_variants to check: 2523"
# [1] "N_variants to remove (maf >=0.45): 118"
# [1] "N_variants to remove: 141"
# [1] "N_variants to flip: 2382"
# [1] "cedars_gsa"
# [1] "N variants in study: 669011"
# [1] "N variants from study in 1000GP: 647961"
# [1] "N_variants different Ref allele: 13967"
# [1] "N_variants to check: 3151"
# [1] "N_variants to remove (maf >=0.45): 111"
# [1] "N_variants to remove: 154"
# [1] "N_variants to flip: 2997"
# [1] "bernstein_gsa"
# [1] "N variants in study: 565467"
# [1] "N variants from study in 1000GP: 563551"
# [1] "N_variants different Ref allele: 9031"
# [1] "N_variants to check: 1962"
# [1] "N_variants to remove (maf >=0.45): 112"
# [1] "N_variants to remove: 127"
# [1] "N_variants to flip: 1835"
# [1] "farkkila_gsa"
# [1] "N variants in study: 528039"
# [1] "N variants from study in 1000GP: 527169"
# [1] "N_variants different Ref allele: 16668"
# [1] "N_variants to check: 1577"
# [1] "N_variants to remove (maf >=0.45): 102"
# [1] "N_variants to remove: 117"
# [1] "N_variants to flip: 1460"
# [1] "franchimont_gsa"
# [1] "N variants in study: 670680"
# [1] "N variants from study in 1000GP: 646368"
# [1] "N_variants different Ref allele: 13269"
# [1] "N_variants to check: 3386"
# [1] "N_variants to remove (maf >=0.45): 110"
# [1] "N_variants to remove: 178"
# [1] "N_variants to flip: 3208"
# [1] "franke_gsa"
# [1] "N variants in study: 570827"
# [1] "N variants from study in 1000GP: 568298"
# [1] "N_variants different Ref allele: 9094"
# [1] "N_variants to check: 2002"
# [1] "N_variants to remove (maf >=0.45): 114"
# [1] "N_variants to remove: 129"
# [1] "N_variants to flip: 1873"
# [1] "helmsley_prism_gsa"
# [1] "N variants in study: 243983"
# [1] "N variants from study in 1000GP: 243339"
# [1] "N_variants different Ref allele: 5970"
# [1] "N_variants to check: 146"
# [1] "N_variants to remove (maf >=0.45): 14"
# [1] "N_variants to remove: 14"
# [1] "N_variants to flip: 132"
# [1] "helmsley_xavier_prism_gsa"
# [1] "N variants in study: 243952"
# [1] "N variants from study in 1000GP: 243309"
# [1] "N_variants different Ref allele: 5582"
# [1] "N_variants to check: 141"
# [1] "N_variants to remove (maf >=0.45): 9"
# [1] "N_variants to remove: 9"
# [1] "N_variants to flip: 132"
# [1] "hyams_protect_gsa"
# [1] "N variants in study: 580630"
# [1] "N variants from study in 1000GP: 577956"
# [1] "N_variants different Ref allele: 11515"
# [1] "N_variants to check: 2174"
# [1] "N_variants to remove (maf >=0.45): 123"
# [1] "N_variants to remove: 137"
# [1] "N_variants to flip: 2037"
# [1] "lewis_sparc_gsa"
# [1] "N variants in study: 619563"
# [1] "N variants from study in 1000GP: 611251"
# [1] "N_variants different Ref allele: 9790"
# [1] "N_variants to check: 2629"
# [1] "N_variants to remove (maf >=0.45): 120"
# [1] "N_variants to remove: 144"
# [1] "N_variants to flip: 2485"
# [1] "mccauley_new_gsa"
# [1] "N variants in study: 591040"
# [1] "N variants from study in 1000GP: 586270"
# [1] "N_variants different Ref allele: 11274"
# [1] "N_variants to check: 2394"
# [1] "N_variants to remove (maf >=0.45): 115"
# [1] "N_variants to remove: 135"
# [1] "N_variants to flip: 2259"
# [1] "mcgovern_gsa"
# [1] "N variants in study: 651675"
# [1] "N variants from study in 1000GP: 639753"
# [1] "N_variants different Ref allele: 12543"
# [1] "N_variants to check: 2852"
# [1] "N_variants to remove (maf >=0.45): 117"
# [1] "N_variants to remove: 157"
# [1] "N_variants to flip: 2695"
# [1] "moayyedi_imagine_gsa"
# [1] "N variants in study: 615197"
# [1] "N variants from study in 1000GP: 609263"
# [1] "N_variants different Ref allele: 8940"
# [1] "N_variants to check: 2423"
# [1] "N_variants to remove (maf >=0.45): 112"
# [1] "N_variants to remove: 131"
# [1] "N_variants to flip: 2292"
# [1] "newberry_share_gsa"
# [1] "N variants in study: 592395"
# [1] "N variants from study in 1000GP: 588149"
# [1] "N_variants different Ref allele: 10396"
# [1] "N_variants to check: 2307"
# [1] "N_variants to remove (maf >=0.45): 116"
# [1] "N_variants to remove: 130"
# [1] "N_variants to flip: 2177"
# [1] "niddk_cho_gsa"
# [1] "N variants in study: 608443"
# [1] "N variants from study in 1000GP: 602279"
# [1] "N_variants different Ref allele: 10938"
# [1] "N_variants to check: 2421"
# [1] "N_variants to remove (maf >=0.45): 112"
# [1] "N_variants to remove: 133"
# [1] "N_variants to flip: 2288"
# [1] "niddk_duerr_gsa"
# [1] "N variants in study: 595967"
# [1] "N variants from study in 1000GP: 589769"
# [1] "N_variants different Ref allele: 8804"
# [1] "N_variants to check: 2393"
# [1] "N_variants to remove (maf >=0.45): 118"
# [1] "N_variants to remove: 140"
# [1] "N_variants to flip: 2253"
# [1] "niddk_rioux_gsa"
# [1] "N variants in study: 589383"
# [1] "N variants from study in 1000GP: 584504"
# [1] "N_variants different Ref allele: 9362"
# [1] "N_variants to check: 2206"
# [1] "N_variants to remove (maf >=0.45): 109"
# [1] "N_variants to remove: 129"
# [1] "N_variants to flip: 2077"
# [1] "niddk_silverberg_gsa"
# [1] "N variants in study: 626523"
# [1] "N variants from study in 1000GP: 618032"
# [1] "N_variants different Ref allele: 10308"
# [1] "N_variants to check: 2550"
# [1] "N_variants to remove (maf >=0.45): 114"
# [1] "N_variants to remove: 136"
# [1] "N_variants to flip: 2414"
# [1] "palotie_hus_gsa"
# [1] "N variants in study: 565121"
# [1] "N variants from study in 1000GP: 562617"
# [1] "N_variants different Ref allele: 14357"
# [1] "N_variants to check: 1931"
# [1] "N_variants to remove (maf >=0.45): 102"
# [1] "N_variants to remove: 119"
# [1] "N_variants to flip: 1812"
# [1] "pekow_share_gsa"
# [1] "N variants in study: 580938"
# [1] "N variants from study in 1000GP: 578515"
# [1] "N_variants different Ref allele: 10436"
# [1] "N_variants to check: 2230"
# [1] "N_variants to remove (maf >=0.45): 125"
# [1] "N_variants to remove: 138"
# [1] "N_variants to flip: 2092"
# [1] "rioux_igenomed_gsa"
# [1] "N variants in study: 551520"
# [1] "N variants from study in 1000GP: 550181"
# [1] "N_variants different Ref allele: 10054"
# [1] "N_variants to check: 1805"
# [1] "N_variants to remove (maf >=0.45): 118"
# [1] "N_variants to remove: 127"
# [1] "N_variants to flip: 1678"
# [1] "sands_msccr_gsa"
# [1] "N variants in study: 603886"
# [1] "N variants from study in 1000GP: 600363"
# [1] "N_variants different Ref allele: 13784"
# [1] "N_variants to check: 2394"
# [1] "N_variants to remove (maf >=0.45): 117"
# [1] "N_variants to remove: 137"
# [1] "N_variants to flip: 2257"
# [1] "stampfer_gsa"
# [1] "N variants in study: 672168"
# [1] "N variants from study in 1000GP: 647980"
# [1] "N_variants different Ref allele: 13055"
# [1] "N_variants to check: 3385"
# [1] "N_variants to remove (maf >=0.45): 119"
# [1] "N_variants to remove: 183"
# [1] "N_variants to flip: 3202"
# [1] "vermeire_gsa"
# [1] "N variants in study: 674915"
# [1] "N variants from study in 1000GP: 649533"
# [1] "N_variants different Ref allele: 10895"
# [1] "N_variants to check: 3389"
# [1] "N_variants to remove (maf >=0.45): 117"
# [1] "N_variants to remove: 181"
# [1] "N_variants to flip: 3208"
# [1] "weersma_gsa"
# [1] "N variants in study: 573599"
# [1] "N variants from study in 1000GP: 571212"
# [1] "N_variants different Ref allele: 9768"
# [1] "N_variants to check: 2043"
# [1] "N_variants to remove (maf >=0.45): 110"
# [1] "N_variants to remove: 124"
# [1] "N_variants to flip: 1919"
# [1] "xavier_prism_gsa"
# [1] "N variants in study: 584058"
# [1] "N variants from study in 1000GP: 580947"
# [1] "N_variants different Ref allele: 9561"
# [1] "N_variants to check: 2197"
# [1] "N_variants to remove (maf >=0.45): 123"
# [1] "N_variants to remove: 141"
# [1] "N_variants to flip: 2056"

# [1] "xavier_share_gsa"
# [1] "N variants in study: 584645"
# [1] "N variants from study in 1000GP: 580888"
# [1] "N_variants different Ref allele: 9398"
# [1] "N_variants to check: 2197"
# [1] "N_variants to remove (maf >=0.45): 119"
# [1] "N_variants to remove: 141"
# [1] "N_variants to flip: 2056"


##############################

##########################################################################################
# 6.5 REMOVE VARIANTS WE CANNOT BE SURE ARE IN THE RIGHT STRAND, AND FLIP THE OTHERS

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_10_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_10_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup \
--exclude ${path_gwas}pre_imputation/QC/${i}/list_variants_to_remove_AT_CG \
--flip ${path_gwas}pre_imputation/QC/${i}/list_variants_to_flip_AT_CG \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip"
done


for i in ${studies[@]}
do echo ${i} && tail -100 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_10_${i} | grep -E "completed|flipped"
done

#####################
# australia_omniexome
# --flip: 7 SNPs flipped.
# Successfully completed.
# gwas1
# --flip: 34268 SNPs flipped.
# Successfully completed.
# gwas2
# --flip: 51089 SNPs flipped.
# Successfully completed.
# pittsburgh_gsa
# --flip: 12 SNPs flipped.
# Successfully completed.
# spain_gsa
# --flip: 5 SNPs flipped.
# Successfully completed.
# italy_gsa
# Successfully completed.
# --flip: 63 SNPs flipped.
# Successfully completed.
# kiel_austria_sibdcs_gsa
# Successfully completed.
# --flip: 29 SNPs flipped.
# Successfully completed.
# netherlands_gsa
# Successfully completed.
# --flip: 51 SNPs flipped.
# Successfully completed.
# slovenia_gsa
# Successfully completed.
# --flip: 28 SNPs flipped.
# Successfully completed.
# sweden_gsa
# Successfully completed.
# --flip: 31 SNPs flipped.
# Successfully completed.
# niddk_broad_gsa
# Successfully completed.
# --flip: 24 SNPs flipped.
# Successfully completed.
# niddk_feinstein_gsa
# Successfully completed.
# --flip: 30 SNPs flipped.
# Successfully completed.
# basque_gsa
# --flip: 14 SNPs flipped.
# Successfully completed.
# lithuania_gsa
# Successfully completed.
# --flip: 1830 SNPs flipped.
# Successfully completed.
# belgium_louis_gsa
# Successfully completed.
# --flip: 2181 SNPs flipped.
# Successfully completed.
# belgium_franchimont_gsa
# Successfully completed.
# --flip: 2106 SNPs flipped.
# Successfully completed.
# belgium_vermeire_gsa
# Successfully completed.
# --flip: 3210 SNPs flipped.
# Successfully completed.
# prism_nfe_gsa
# Successfully completed.
# --flip: 1874 SNPs flipped.
# Successfully completed.
# prism_nfe_gwas
# Successfully completed.
# --flip: 131 SNPs flipped.
# Successfully completed.
# finland_illugwas
# Successfully completed.
# --flip: 131 SNPs flipped.
# Successfully completed.
# german_affy6_old_gwas
# --flip: 63125 SNPs flipped.
# Successfully completed.
# norway_affy6_old_gwas
# --flip: 59703 SNPs flipped.
# Successfully completed.
# belgium_inf1_old_gwas
# --flip: 0 SNPs flipped.
# Successfully completed.
# belgium_inf2_old_gwas
# --flip: 0 SNPs flipped.
# Successfully completed.
# cedars_370k_old_gwas
# --flip: 7 SNPs flipped.
# Successfully completed.
# cedars_610k_old_gwas
# --flip: 6 SNPs flipped.
# Successfully completed.
# cedars_omni_old_gwas
# --flip: 7 SNPs flipped.
# Successfully completed.
# swedish_uc_old_gwas
# --flip: 0 SNPs flipped.
# Successfully completed.
# mccauley_gsa
# Successfully completed.
# --flip: 1940 SNPs flipped.
# Successfully completed.
# ccfa_gsa
# Successfully completed.
# --flip: 2382 SNPs flipped.
# Successfully completed.
# cedars_gsa
# Successfully completed.
# --flip: 2997 SNPs flipped.
# Successfully completed.
# bernstein_gsa
# Successfully completed.
# --flip: 1835 SNPs flipped.
# Successfully completed.
# farkkila_gsa
# Successfully completed.
# --flip: 1460 SNPs flipped.
# Successfully completed.
# franchimont_gsa
# Successfully completed.
# --flip: 3208 SNPs flipped.
# Successfully completed.
# franke_gsa
# Successfully completed.
# --flip: 1873 SNPs flipped.
# Successfully completed.
# helmsley_prism_gsa
# Successfully completed.
# --flip: 132 SNPs flipped.
# Successfully completed.
# helmsley_xavier_prism_gsa
# Successfully completed.
# --flip: 132 SNPs flipped.
# Successfully completed.
# hyams_protect_gsa
# Successfully completed.
# --flip: 2037 SNPs flipped.
# Successfully completed.
# lewis_sparc_gsa
# Successfully completed.
# --flip: 2485 SNPs flipped.
# Successfully completed.
# mccauley_new_gsa
# Successfully completed.
# --flip: 2259 SNPs flipped.
# Successfully completed.
# mcgovern_gsa
# Successfully completed.
# --flip: 2695 SNPs flipped.
# Successfully completed.
# moayyedi_imagine_gsa
# Successfully completed.
# --flip: 2292 SNPs flipped.
# Successfully completed.
# newberry_share_gsa
# Successfully completed.
# --flip: 2177 SNPs flipped.
# Successfully completed.
# niddk_cho_gsa
# Successfully completed.
# --flip: 2288 SNPs flipped.
# Successfully completed.
# niddk_duerr_gsa
# Successfully completed.
# --flip: 2253 SNPs flipped.
# Successfully completed.
# niddk_rioux_gsa
# Successfully completed.
# --flip: 2077 SNPs flipped.
# Successfully completed.
# niddk_silverberg_gsa
# Successfully completed.
# --flip: 2414 SNPs flipped.
# Successfully completed.
# palotie_hus_gsa
# Successfully completed.
# --flip: 1812 SNPs flipped.
# Successfully completed.
# pekow_share_gsa
# Successfully completed.
# --flip: 2092 SNPs flipped.
# Successfully completed.
# rioux_igenomed_gsa
# Successfully completed.
# --flip: 1678 SNPs flipped.
# Successfully completed.
# sands_msccr_gsa
# Successfully completed.
# --flip: 2257 SNPs flipped.
# Successfully completed.
# stampfer_gsa
# Successfully completed.
# --flip: 3202 SNPs flipped.
# Successfully completed.
# vermeire_gsa
# Successfully completed.
# --flip: 3208 SNPs flipped.
# Successfully completed.
# weersma_gsa
# Successfully completed.
# --flip: 1919 SNPs flipped.
# Successfully completed.
# xavier_prism_gsa
# Successfully completed.
# --flip: 2056 SNPs flipped.
# Successfully completed.
# xavier_share_gsa
# Successfully completed.
# --flip: 2056 SNPs flipped.
# Successfully completed.

########################################################################################################################################
