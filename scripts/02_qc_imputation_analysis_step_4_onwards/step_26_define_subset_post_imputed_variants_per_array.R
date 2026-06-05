# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
########################################################################################################
# 26.1 Define per array set of variants monomorphic in all studies - to decrease size of files

# no monomorphic variants in out put results. keep them all


########################################################################################################
# 26.2 Create a list of variants that do not pass HWE in EUR - to exclude



path_gwas=/path/to/ibdgwas/IIBDGC/

##### STUDIES WITH CONTROLS ONLY

path_gwas=/path/to/ibdgwas/IIBDGC/

studies_withctr=(kiel_austria_sibdcs_gsa)
MEM=40000

studies_withctr=(all_hce) 
MEM=60000

studies_withctr=(niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa netherlands_gsa slovenia_gsa 
                 sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa 
                 belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas german_affy6_old_gwas norway_affy6_old_gwas 
                 belgium_inf1_old_gwas belgium_inf2_old_gwas  
                 swedish_uc_old_gwas cedars_gsa bernstein_gsa franchimont_gsa franke_gsa 
                 helmsley_prism_gsa helmsley_xavier_prism_gsa mccauley_new_gsa mcgovern_gsa 
                 moayyedi_imagine_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa 
                 sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa 
                 xavier_prism_gsa)
MEM=30000

for i in ${studies_withctr[@]}
do 
for chr in {1..22}
do
bsub -J"kg1" -M"$MEM" -n 4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_26.1_${i}_${chr} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_26.1_${i}_${chr} \
"/path/to/software/./plink2 \
--vcf ${path_gwas}imputed/${i}/2022/eur/chr${chr}.dose.vcf.gz \
--double-id \
--keep ${path_gwas}pre_imputation/QC/${i}/list_controls2 \
--hardy \
--threads 4 \
--memory ${MEM} \
--out ${path_gwas}imputed/${i}/2022/eur/chr${chr}_hardy_ctr"
done
done



# CHRX - use plink1.9, plink2 assumes all in chr23 - haploid (no chr25-PAR labels in vcf?) - only in control females for non-PAR

###### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

tmp<-fread(paste(path,"pheno/gwas-mega-2-core_updated_sex_gender_Dec22.txt.gz",sep=""),head=T)

studies<-c("all_hce","niddk_old_gwas","australia_omniexome","gwas1","gwas2","pittsburgh_gsa","spain_gsa","italy_gsa",
           "kiel_austria_sibdcs_gsa","netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa",
           "niddk_feinstein_gsa","basque_gsa","lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa",
           "belgium_vermeire_gsa","prism_nfe_gsa","prism_nfe_gwas","finland_illugwas","german_affy6_old_gwas",
           "norway_affy6_old_gwas","belgium_inf1_old_gwas","belgium_inf2_old_gwas",
           "cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas",
           # "swedish_uc_old_gwas",
           "mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

tmp<-tmp[which(tmp$SNPSEX2==2),]


for (j in 1:length(studies)) {
  
  print(studies[j])
  tmp1<-fread(paste(path,"pre_imputation/QC/",studies[j],"/list_controls2",sep=""),head=F)

  if(nrow(tmp1)>0) {
    
    tmp_ctr_fe<-tmp1[which(tmp1$V1 %in% tmp$IID),]
    write.table(tmp_ctr_fe,paste(path,"pre_imputation/QC/",studies[j],"/list_female_controls2",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
    rm(tmp_ctr_fe)
    
    tmp_ca_fe<-tmp[which(! (tmp$IID %in% tmp1$V1) & (tmp$study==studies[j])),c("IID","FID")]
    write.table(tmp_ca_fe,paste(path,"pre_imputation/QC/",studies[j],"/list_female_cases2",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
    rm(tmp_ca_fe)
    
  } else {
    
    tmp_ca_fe<-tmp[which(! (tmp$IID %in% tmp1$V1) & (tmp$study==studies[j])),c("IID","FID")]
    write.table(tmp_ca_fe,paste(path,"pre_imputation/QC/",studies[j],"/list_female_cases2",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
    rm(tmp_ca_fe)
    
  }
  
  rm(tmp1)

}


############

for i in ${studies_withctr[@]}
do 
bsub -J"kg1" -M"$MEM" -n 4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_26.1_${i}_chrX \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_26.1_${i}_chrX \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--vcf ${path_gwas}imputed/${i}/2022/eur/chrX.dose.vcf.gz \
--double-id \
--keep ${path_gwas}pre_imputation/QC/${i}/list_female_controls2 \
--hardy \
--threads 4 \
--memory ${MEM} \
--out ${path_gwas}imputed/${i}/2022/eur/chrX_hardy_ctr"
done




##### STUDIES WITH CASES

studies=(niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa netherlands_gsa slovenia_gsa 
         sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa 
         belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas 
         belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas 
         swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa 
         helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa 
         moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa 
         palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa 
         xavier_prism_gsa xavier_share_gsa)
MEM=30000
for i in ${studies[@]}
do 
for chr in {1..22}
do
bsub -J"kg1" -M"$MEM" -n 4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_26.2_${i}_${chr} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_26.2_${i}_${chr} \
"/path/to/software/./plink2 \
--vcf ${path_gwas}imputed/${i}/2022/eur/chr${chr}.dose.vcf.gz \
--double-id \
--remove ${path_gwas}pre_imputation/QC/${i}/list_controls2 \
--hardy \
--threads 4 \
--memory ${MEM} \
--out ${path_gwas}imputed/${i}/2022/eur/chr${chr}_hardy_cases"
done
done



for i in ${studies[@]}
do 
bsub -J"kg1" -M"$MEM" -n 4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_26.2_${i}_chrX \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_26.2_${i}_chrX \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--vcf ${path_gwas}imputed/${i}/2022/eur/chrX.dose.vcf.gz \
--double-id \
--keep ${path_gwas}pre_imputation/QC/${i}/list_female_cases2 \
--hardy \
--threads 4 \
--memory ${MEM} \
--out ${path_gwas}imputed/${i}/2022/eur/chrX_hardy_cases"
done
done


studies_withctr=(kiel_austria_sibdcs_gsa all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa netherlands_gsa slovenia_gsa 
                 sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa 
                 belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas german_affy6_old_gwas norway_affy6_old_gwas 
                 belgium_inf1_old_gwas belgium_inf2_old_gwas  
                 swedish_uc_old_gwas cedars_gsa bernstein_gsa franchimont_gsa franke_gsa 
                 helmsley_prism_gsa helmsley_xavier_prism_gsa mccauley_new_gsa mcgovern_gsa 
                 moayyedi_imagine_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa 
                 sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa 
                 xavier_prism_gsa)

for i in ${studies_withctr[@]}
do echo ${i} && for chr in {1..22}
do echo ${chr} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_26.1_${i}_${chr} | grep "Successfully"
done
done

for i in ${studies_withctr[@]}
do echo ${i} && ls -la ${path_gwas}imputed/${i}/2022/eur/*_hardy_ctr.hardy | wc -l
done
done



studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do echo ${i} && for chr in {1..22}
do echo ${chr} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_26.2_${i}_${chr} | grep "Successfully"
done
done



#### chrX

for i in ${studies_withctr[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_26.1_${i}_chrX | grep "Successfully"
done

for i in ${studies_withctr[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_26.2_${i}_chrX | grep "Successfully"
done




####################
#### merge files

studies_withctr=(all_hce kiel_austria_sibdcs_gsa niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa netherlands_gsa slovenia_gsa 
                 sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa 
                 belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas german_affy6_old_gwas norway_affy6_old_gwas 
                 belgium_inf1_old_gwas belgium_inf2_old_gwas  
                 swedish_uc_old_gwas cedars_gsa bernstein_gsa franchimont_gsa franke_gsa 
                 helmsley_prism_gsa helmsley_xavier_prism_gsa mccauley_new_gsa mcgovern_gsa 
                 moayyedi_imagine_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa 
                 sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa 
                 xavier_prism_gsa)

MEM=250
for i in ${studies_withctr[@]}
do
bsub -J"kg1" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_26.3_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_26.3_${i}} \
"cat ${path_gwas}imputed/${i}/2022/eur/chr*_hardy_ctr.hardy | gzip > ${path_gwas}imputed/${i}/2022/eur/allchr_hardy_ctr.hardy.gz"
done

for i in ${studies_withctr[@]}
do 
ls -la ${path_gwas}imputed/${i}/2022/eur/allchr_hardy_ctr.hardy.gz
done



studies=(all_hce kiel_austria_sibdcs_gsa niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)
MEM=250
for i in ${studies[@]}
do 
bsub -J"kg1" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_26.4_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_26.4_${i}} \
"cat ${path_gwas}imputed/${i}/2022/eur/chr*_hardy_cases.hardy | gzip > ${path_gwas}imputed/${i}/2022/eur/allchr_hardy_cases.hardy.gz"
done

for i in ${studies[@]}
do 
ls -la ${path_gwas}imputed/${i}/2022/eur/allchr_hardy_cases.hardy.gz
done



# # create list with variants to exclude from analysis, only in controls if study has ctr, otherwise in cases with less stringent threshold

MEM=250
for i in ${studies_withctr[@]}
do
bsub -J"kg1" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_26.4_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_26.4_${i} \
"zcat ${path_gwas}imputed/${i}/2022/eur/allchr_hardy_ctr.hardy.gz | awk '\$10 < 1E-5' \
| gzip > ${path_gwas}imputed/${i}/2022/eur/list_variants_fail_hwe_ctr.gz"
done

MEM=250
for i in ${studies[@]}
do
bsub -J"kg1" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_26.5_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_26.5_${i} \
"zcat ${path_gwas}imputed/${i}/2022/eur/allchr_hardy_cases.hardy.gz | awk '\$10 < 1E-12' \
| gzip > ${path_gwas}imputed/${i}/2022/eur/list_variants_fail_hwe_cases.gz"
done



## CHRX create list with variants to exclude from analysis, only in controls if study has ctr, otherwise in cases with less stringent threshold

MEM=250
for i in ${studies_withctr[@]}
do
bsub -J"kg1" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_26.4_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_26.4_${i} \
"cat ${path_gwas}imputed/${i}/2022/eur/chrX_hardy_ctr.hwe | awk '\$9 < 1E-5' \
| gzip > ${path_gwas}imputed/${i}/2022/eur/list_variants_fail_hwe_chrX_ctr.hwe.gz"
done

MEM=250
for i in ${studies[@]}
do
bsub -J"kg1" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_26.5_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_26.5_${i} \
"cat ${path_gwas}imputed/${i}/2022/eur/chrX_hardy_cases.hwe | awk '\$9 < 1E-12' \
| gzip > ${path_gwas}imputed/${i}/2022/eur/list_variants_fail_hwe_chrX_cases.hwe.gz"
done




# double check that none of the excluded variants are in known IBD regions:

#### /software/R-4.3.1/bin/R

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


ibd<-read.table(paste(path,"post_imputation/analysis/metaanalysis/summary_results/list_index_variants_updated_with_NOD2_ADCY7_lifted_hg38_no_overlaps.bed",sep=""))
dim(ibd)
ibd$chr<-as.numeric(gsub("chr","",ibd$V1))

dat<-as.data.frame(matrix(ncol=4,nrow=length(studies)))
colnames(dat)<-c("studies","N_variants_excluded_hwe","N_variants_excluded_hwe_in_associated_regions","group")


for (i in 1:length(studies)) {
  
  print(studies[i])
  
  dat$studies[i]<-studies[i]
  
  file_cases<-paste(path,"imputed/",studies[i],"/2022/eur/list_variants_fail_hwe_cases.gz",sep="")
  if(file.exists(file_cases)) {
    all<-try(fread(paste(path,"imputed/",studies[i],"/2022/eur/list_variants_fail_hwe_cases.gz",sep="")),silent = TRUE)
    dat$group[i]<-"cases"
  }
  
  # overwrite for studies with ctr data, if that exists
  file_ctr<-paste(path,"imputed/",studies[i],"/2022/eur/list_variants_fail_hwe_ctr.gz",sep="")
  if(file.exists(file_ctr)) {
    all<-try(fread(file_ctr),silent = TRUE)
    dat$group[i]<-"ctr"
  }
  
  # get number of excluded variants
  if (!is.null(nrow(all))) {
    
    all$position<-gsub("chr[0-9]{1,2}:","",all$V2)
    all$position<-as.numeric(gsub(":.*","",all$position))
      
    for(ii in 1:nrow(ibd)) {
      if(ii==1) {
        all_assoc<-all[which(all$V1==ibd$chr[ii] & all$position>=ibd$V2[ii] & all$position<=ibd$V2[ii]),]
      }else{
        tmp<-all[which(all$V1==ibd$chr[ii] & all$position>=ibd$V2[ii] & all$position<=ibd$V3[ii]),]
        all_assoc<-rbind(all_assoc,tmp)
        rm(tmp)
      }
    }
    
    dat$N_variants_excluded_hwe[i]<-nrow(all)
    dat$N_variants_excluded_hwe_in_associated_regions[i]<-nrow(all_assoc)
    

    fwrite(all[,2],paste(path,"imputed/",studies[i],"/2022/eur/list_variants_post_imputation_hwe_toremove.tsv",sep=""),
             col.names=F,row.names=F,quote=F,sep="\t")
    rm(all,all_assoc)
    
  } else {
    dat$N_variants_excluded_hwe[i]<-0
    dat$N_variants_excluded_hwe_in_associated_regions[i]<-0
    
    all<-as.data.frame(matrix(ncol=1,nrow=0))
    fwrite(all,paste(path,"imputed/",studies[i],"/2022/eur/list_variants_post_imputation_hwe_toremove.tsv",sep=""),
           col.names=F,row.names=F,quote=F,sep="\t")
    rm(all)
  }
  
}

table(dat$N_variants_excluded_hwe_in_associated_regions)
# 0 
# 58

table(dat$group,dat$N_variants_excluded_hwe_in_associated_regions)
#        0  1  2  3  4  5  6  8 11 16 19 23 25 29 34 35 36 47 53 66 71 72 100
# cases 14  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0   0
# ctr    6  5  3  2  1  2  1  2  1  2  2  2  1  1  2  1  1  1  1  1  1  1   1
# 
#       116 176 184
# cases   0   0   0
# ctr     1   1   1


q("no")


####################
#                      studies N_variants_excluded_hwe
# 1                    all_hce                     643
# 2             niddk_old_gwas                     508
# 3        australia_omniexome                      20
# 4                      gwas1                     169
# 5                      gwas2                     362
# 6             pittsburgh_gsa                     394
# 7                  spain_gsa                     257
# 8                  italy_gsa                     112
# 9    kiel_austria_sibdcs_gsa                    1332
# 10           netherlands_gsa                      74
# 11              slovenia_gsa                     107
# 12                sweden_gsa                      60
# 13           niddk_broad_gsa                     231
# 14       niddk_feinstein_gsa                    1232
# 15                basque_gsa                     498
# 16             lithuania_gsa                     180
# 17         belgium_louis_gsa                      48
# 18   belgium_franchimont_gsa                      66
# 19      belgium_vermeire_gsa                     107
# 20             prism_nfe_gsa                      12
# 21            prism_nfe_gwas                      42
# 22          finland_illugwas                       0
# 23     german_affy6_old_gwas                     838
# 24     norway_affy6_old_gwas                      20
# 25     belgium_inf1_old_gwas                      30
# 26     belgium_inf2_old_gwas                       5
# 27      cedars_370k_old_gwas                       0
# 28      cedars_610k_old_gwas                       0
# 29      cedars_omni_old_gwas                      74
# 30       swedish_uc_old_gwas                      42
# 31              mccauley_gsa                      11
# 32                  ccfa_gsa                       0
# 33                cedars_gsa                     492
# 34             bernstein_gsa                       0
# 35              farkkila_gsa                       0
# 36           franchimont_gsa                    1212
# 37                franke_gsa                      37
# 38        helmsley_prism_gsa                      40
# 39 helmsley_xavier_prism_gsa                     119
# 40         hyams_protect_gsa                       0
# 41           lewis_sparc_gsa                       0
# 42          mccauley_new_gsa                      19
# 43              mcgovern_gsa                     534
# 44      moayyedi_imagine_gsa                      49
# 45        newberry_share_gsa                       0
# 46             niddk_cho_gsa                     128
# 47           niddk_duerr_gsa                     319
# 48           niddk_rioux_gsa                      56
# 49      niddk_silverberg_gsa                     498
# 50           palotie_hus_gsa                       0
# 51           pekow_share_gsa                       0
# 52        rioux_igenomed_gsa                       0
# 53           sands_msccr_gsa                      43
# 54              stampfer_gsa                     113
# 55              vermeire_gsa                     100
# 56               weersma_gsa                      58
# 57          xavier_prism_gsa                       6
# 58          xavier_share_gsa                       0
#    N_variants_excluded_hwe_in_associated_regions group
# 1                                             34   ctr
# 2                                             11   ctr
# 3                                              0   ctr
# 4                                             72   ctr
# 5                                            100   ctr
# 6                                             35   ctr
# 7                                             34   ctr
# 8                                             66   ctr
# 9                                            116   ctr
# 10                                             1   ctr
# 11                                             0   ctr
# 12                                            19   ctr
# 13                                            47   ctr
# 14                                           184   ctr
# 15                                            71   ctr
# 16                                            36   ctr
# 17                                             5   ctr
# 18                                            16   ctr
# 19                                             8   ctr
# 20                                             0   ctr
# 21                                             1   ctr
# 22                                             0 cases
# 23                                            53   ctr
# 24                                             1   ctr
# 25                                             1   ctr
# 26                                             0   ctr
# 27                                             0 cases
# 28                                             0 cases
# 29                                             0 cases
# 30                                             5   ctr
# 31                                             0 cases
# 32                                             0 cases
# 33                                            16   ctr
# 34                                             0   ctr
# 35                                             0 cases
# 36                                           176   ctr
# 37                                             2   ctr
# 38                                             3   ctr
# 39                                             3   ctr
# 40                                             0 cases
# 41                                             0 cases
# 42                                             0   ctr
# 43                                            19   ctr
# 44                                            23   ctr
# 45                                             0 cases
# 46                                             6   ctr
# 47                                            29   ctr
# 48                                             2   ctr
# 49                                            23   ctr
# 50                                             0 cases
# 51                                             0 cases
# 52                                             0 cases
# 53                                             2   ctr
# 54                                            25   ctr
# 55                                             8   ctr
# 56                                             4   ctr
# 57                                             1   ctr
# 58                                             0 cases
####################

studies=(all_hce kiel_austria_sibdcs_gsa niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do wc -l ${path_gwas}imputed/${i}/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
done


####################
# 643 /path/to/ibdgwas/IIBDGC/imputed/all_hce/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 1332 /path/to/ibdgwas/IIBDGC/imputed/kiel_austria_sibdcs_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 508 /path/to/ibdgwas/IIBDGC/imputed/niddk_old_gwas/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 20 /path/to/ibdgwas/IIBDGC/imputed/australia_omniexome/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 169 /path/to/ibdgwas/IIBDGC/imputed/gwas1/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 362 /path/to/ibdgwas/IIBDGC/imputed/gwas2/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 394 /path/to/ibdgwas/IIBDGC/imputed/pittsburgh_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 257 /path/to/ibdgwas/IIBDGC/imputed/spain_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 112 /path/to/ibdgwas/IIBDGC/imputed/italy_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 74 /path/to/ibdgwas/IIBDGC/imputed/netherlands_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 107 /path/to/ibdgwas/IIBDGC/imputed/slovenia_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 60 /path/to/ibdgwas/IIBDGC/imputed/sweden_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 231 /path/to/ibdgwas/IIBDGC/imputed/niddk_broad_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 1232 /path/to/ibdgwas/IIBDGC/imputed/niddk_feinstein_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 498 /path/to/ibdgwas/IIBDGC/imputed/basque_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 180 /path/to/ibdgwas/IIBDGC/imputed/lithuania_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 48 /path/to/ibdgwas/IIBDGC/imputed/belgium_louis_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 66 /path/to/ibdgwas/IIBDGC/imputed/belgium_franchimont_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 107 /path/to/ibdgwas/IIBDGC/imputed/belgium_vermeire_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 12 /path/to/ibdgwas/IIBDGC/imputed/prism_nfe_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 42 /path/to/ibdgwas/IIBDGC/imputed/prism_nfe_gwas/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 0 /path/to/ibdgwas/IIBDGC/imputed/finland_illugwas/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 838 /path/to/ibdgwas/IIBDGC/imputed/german_affy6_old_gwas/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 20 /path/to/ibdgwas/IIBDGC/imputed/norway_affy6_old_gwas/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 30 /path/to/ibdgwas/IIBDGC/imputed/belgium_inf1_old_gwas/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 5 /path/to/ibdgwas/IIBDGC/imputed/belgium_inf2_old_gwas/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 0 /path/to/ibdgwas/IIBDGC/imputed/cedars_370k_old_gwas/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 0 /path/to/ibdgwas/IIBDGC/imputed/cedars_610k_old_gwas/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 74 /path/to/ibdgwas/IIBDGC/imputed/cedars_omni_old_gwas/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 42 /path/to/ibdgwas/IIBDGC/imputed/swedish_uc_old_gwas/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 11 /path/to/ibdgwas/IIBDGC/imputed/mccauley_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 0 /path/to/ibdgwas/IIBDGC/imputed/ccfa_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 492 /path/to/ibdgwas/IIBDGC/imputed/cedars_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 0 /path/to/ibdgwas/IIBDGC/imputed/bernstein_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 0 /path/to/ibdgwas/IIBDGC/imputed/farkkila_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 1212 /path/to/ibdgwas/IIBDGC/imputed/franchimont_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 37 /path/to/ibdgwas/IIBDGC/imputed/franke_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 40 /path/to/ibdgwas/IIBDGC/imputed/helmsley_prism_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 119 /path/to/ibdgwas/IIBDGC/imputed/helmsley_xavier_prism_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 0 /path/to/ibdgwas/IIBDGC/imputed/hyams_protect_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 0 /path/to/ibdgwas/IIBDGC/imputed/lewis_sparc_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 19 /path/to/ibdgwas/IIBDGC/imputed/mccauley_new_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 534 /path/to/ibdgwas/IIBDGC/imputed/mcgovern_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 49 /path/to/ibdgwas/IIBDGC/imputed/moayyedi_imagine_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 0 /path/to/ibdgwas/IIBDGC/imputed/newberry_share_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 128 /path/to/ibdgwas/IIBDGC/imputed/niddk_cho_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 319 /path/to/ibdgwas/IIBDGC/imputed/niddk_duerr_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 56 /path/to/ibdgwas/IIBDGC/imputed/niddk_rioux_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 498 /path/to/ibdgwas/IIBDGC/imputed/niddk_silverberg_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 0 /path/to/ibdgwas/IIBDGC/imputed/palotie_hus_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 0 /path/to/ibdgwas/IIBDGC/imputed/pekow_share_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 0 /path/to/ibdgwas/IIBDGC/imputed/rioux_igenomed_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 43 /path/to/ibdgwas/IIBDGC/imputed/sands_msccr_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 113 /path/to/ibdgwas/IIBDGC/imputed/stampfer_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 100 /path/to/ibdgwas/IIBDGC/imputed/vermeire_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 58 /path/to/ibdgwas/IIBDGC/imputed/weersma_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 6 /path/to/ibdgwas/IIBDGC/imputed/xavier_prism_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
# 0 /path/to/ibdgwas/IIBDGC/imputed/xavier_share_gsa/2022/eur/list_variants_post_imputation_hwe_toremove.tsv
####################

#### SAME FOR CHRX

MEM=4000

bsub -Is -M"$MEM" -m "modern_hardware" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1


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


dat<-as.data.frame(matrix(ncol=4,nrow=length(studies)))
colnames(dat)<-c("studies","N_variants_excluded_hwe","group","N_variants_excluded_hwe_in_associated_regions")


for (i in 1:length(studies)) {
  
  print(studies[i])
  
  dat$studies[i]<-studies[i]
  
  file_cases<-paste(path,"imputed/",studies[i],"/2022/eur/list_variants_fail_hwe_chrX_cases.hwe.gz",sep="")

  if(file.exists(file_cases)) {
    all<-try(fread(paste(path,"imputed/",studies[i],"/2022/eur/list_variants_fail_hwe_chrX_cases.hwe.gz",sep="")),silent = TRUE)
    dat$group[i]<-"cases"
  }
  
  # overwrite for studies with ctr data, if that exists
  file_ctr<-paste(path,"imputed/",studies[i],"/2022/eur/list_variants_fail_hwe_chrX_ctr.hwe.gz",sep="")
  if(file.exists(file_ctr)) {
    all<-try(fread(file_ctr),silent = TRUE)
    dat$group[i]<-"ctr"
  }
  
  # get number of excluded variants
  if (!is.null(nrow(all))) {
    
    all$position<-gsub("chrX:","",all$V2)
    all$position<-as.numeric(gsub(":.*","",all$position))
    
    for(ii in 1:nrow(ibd)) {
      if(ii==1) {
        all_assoc<-all[which(all$V1==ibd$chr[ii] & all$position>=ibd$V2[ii] & all$position<=ibd$V2[ii]),]
      }else{
        tmp<-all[which(all$V1==ibd$chr[ii] & all$position>=ibd$V2[ii] & all$position<=ibd$V3[ii]),]
        all_assoc<-rbind(all_assoc,tmp)
        rm(tmp)
      }
    }
    
    dat$N_variants_excluded_hwe[i]<-nrow(all)
    dat$N_variants_excluded_hwe_in_associated_regions[i]<-nrow(all_assoc)
    
    
    fwrite(all[,2],paste(path,"imputed/",studies[i],"/2022/eur/list_variants_post_imputation_hwe_toremove_chrX.tsv",sep=""),
           col.names=F,row.names=F,quote=F,sep="\t")
    rm(all,all_assoc)
    
  } else {
    dat$N_variants_excluded_hwe[i]<-0
    dat$N_variants_excluded_hwe_in_associated_regions[i]<-0
    
    all<-as.data.frame(matrix(ncol=1,nrow=0))
    fwrite(all,paste(path,"imputed/",studies[i],"/2022/eur/list_variants_post_imputation_hwe_toremove_chrX.tsv",sep=""),
           col.names=F,row.names=F,quote=F,sep="\t")
    rm(all)
  }
  
}


dat
# studies N_variants_excluded_hwe group
# 1                    all_hce                       5   ctr
# 2             niddk_old_gwas                       0   ctr
# 3        australia_omniexome                       0   ctr
# 4                      gwas1                       0   ctr
# 5                      gwas2                       4   ctr
# 6             pittsburgh_gsa                       5   ctr
# 7                  spain_gsa                       1   ctr
# 8                  italy_gsa                       1   ctr
# 9    kiel_austria_sibdcs_gsa                      12   ctr
# 10           netherlands_gsa                       0   ctr
# 11              slovenia_gsa                       2   ctr
# 12                sweden_gsa                       2   ctr
# 13           niddk_broad_gsa                      36   ctr
# 14       niddk_feinstein_gsa                     244   ctr # only 1 in PAR1, 0 in PAR2
# 15                basque_gsa                       2   ctr
# 16             lithuania_gsa                       0   ctr
# 17         belgium_louis_gsa                       2   ctr
# 18   belgium_franchimont_gsa                       6   ctr
# 19      belgium_vermeire_gsa                       2   ctr
# 20             prism_nfe_gsa                       0   ctr
# 21            prism_nfe_gwas                       0   ctr
# 22          finland_illugwas                       0  <NA>
#   23     german_affy6_old_gwas                      22   ctr
# 24     norway_affy6_old_gwas                       0   ctr
# 25     belgium_inf1_old_gwas                       0   ctr
# 26     belgium_inf2_old_gwas                       1   ctr
# 27      cedars_370k_old_gwas                       0  <NA>
#   28      cedars_610k_old_gwas                       0  <NA>
#   29      cedars_omni_old_gwas                       0  <NA>
#   30       swedish_uc_old_gwas                       0   ctr
# 31              mccauley_gsa                       0  <NA>
#   32                  ccfa_gsa                       0  <NA>
#   33                cedars_gsa                       6   ctr
# 34             bernstein_gsa                       0   ctr
# 35              farkkila_gsa                       0  <NA>
#   36           franchimont_gsa                      10   ctr
# 37                franke_gsa                       0   ctr
# 38        helmsley_prism_gsa                       0   ctr
# 39 helmsley_xavier_prism_gsa                       0   ctr
# 40         hyams_protect_gsa                       0  <NA>
#   41           lewis_sparc_gsa                       0  <NA>
#   42          mccauley_new_gsa                       1   ctr
# 43              mcgovern_gsa                       2   ctr
# 44      moayyedi_imagine_gsa                       0   ctr
# 45        newberry_share_gsa                       0  <NA>
#   46             niddk_cho_gsa                       1   ctr
# 47           niddk_duerr_gsa                       8   ctr
# 48           niddk_rioux_gsa                       3   ctr
# 49      niddk_silverberg_gsa                       6   ctr
# 50           palotie_hus_gsa                       0  <NA>
#   51           pekow_share_gsa                       0  <NA>
#   52        rioux_igenomed_gsa                       0  <NA>
#   53           sands_msccr_gsa                       0   ctr
# 54              stampfer_gsa                       0   ctr
# 55              vermeire_gsa                       1   ctr
# 56               weersma_gsa                       0   ctr
# 57          xavier_prism_gsa                       0   ctr
# 58          xavier_share_gsa                       0  <NA>
