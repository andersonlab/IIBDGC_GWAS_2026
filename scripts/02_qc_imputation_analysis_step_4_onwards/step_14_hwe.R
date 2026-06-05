# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
############################################################################################################################################

################################
# 14 - HWE IN EUR ONLY SAMPLES #
################################

# samples genotyped in one batch
# evaluate HWE per batch and within EUR non-Jewish ancestry samples

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
           "mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

for (j in 1:length(cohorts)) {
  
  print(cohorts[j])
  
  tmp<-read.table(paste(path,
                        "pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy.fam",sep=""),head=F)
  colnames(tmp)<-c("FID","IID")
  write.table(tmp,paste(path,"pre_imputation/QC/",cohorts[j],"/list_ids_",cohorts[j],sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  
}

######

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_41_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_41_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample \
--allow-no-sex \
--keep ${path_gwas}pre_imputation/QC/pca_1000gp/list_eur_nonjewish_ancestry_samples_noDuplicates_perstudy \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/hwe_eur_nonjewish_autosomal_${i}"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_41_${i} | grep -E "completed"
done

for i in ${studies[@]}
do echo ${i} \
&& wc -l ${path_gwas}pre_imputation/QC/${i}/hwe_eur_nonjewish_autosomal_${i}.fam
done


for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_42_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_42_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/hwe_eur_nonjewish_autosomal_${i} \
--allow-no-sex \
--hardy \
--out ${path_gwas}pre_imputation/QC/${i}/hwe_eur_nonjewish_autosomal_${i}"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_42_${i} | grep -E "completed"
done


###

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_43_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_43_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample \
--allow-no-sex \
--filter-females --chr 23 25 \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/hwe_chr23_females_${i}"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_43_${i} | grep -E "completed"
done


for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_44_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_44_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/hwe_chr23_females_${i} \
--allow-no-sex \
--keep ${path_gwas}pre_imputation/QC/pca_1000gp/list_eur_nonjewish_ancestry_samples_noDuplicates --chr 23 25 \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/hwe_eur_nonjewish_chr23_females_${i}"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_44_${i} | grep -E "completed"
done


for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_45_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_45_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/hwe_eur_nonjewish_chr23_females_${i} \
--allow-no-sex \
--hardy \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/hwe_eur_nonjewish_chr23_females_${i}"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_45_${i} | grep -E "completed"
done


########## /software/R-4.3.1/bin/R

path<-"/path/to/ibdgwas/IIBDGC/"
cohorts<-c("all_hce","niddk_old_gwas","australia_omniexome","gwas1","gwas2","pittsburgh_gsa","spain_gsa","italy_gsa",
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

for (j in 1:length(cohorts)) {
  
  print(cohorts[j])

  file_hwe_aut<-paste(path,"pre_imputation/QC/",cohorts[j],"/hwe_eur_nonjewish_autosomal_",cohorts[j],".hwe",
                      sep="")
  file_hwe_23<-paste(path,"pre_imputation/QC/",cohorts[j],"/hwe_eur_nonjewish_chr23_females_",cohorts[j],".hwe",
                      sep="")
  
  if(file.exists(file_hwe_aut)) {
    
    aut_hwe<-read.table(file_hwe_aut,head=T)
    aut_hwe<-aut_hwe[which(!aut_hwe$CHR %in% c(23,25)),]
    dim(aut_hwe[which(aut_hwe$TEST=="AFF" & aut_hwe$P<=1E-12),])
    remove_1<-aut_hwe[which(aut_hwe$TEST=="UNAFF" & aut_hwe$P<=1E-5),]
    remove_2<-aut_hwe[which(aut_hwe$TEST=="AFF" & aut_hwe$P<=1E-12),]
    remove_2<-remove_2[which(!remove_2$SNP %in% remove_1$SNP),]
    
  }
  
  if(file.exists(file_hwe_23)){
    
    chr23_hwe<-read.table(file_hwe_23,head=T)
    remove_3<-chr23_hwe[which(chr23_hwe$TEST=="UNAFF" & chr23_hwe$P<=1E-5),]
    remove_4<-chr23_hwe[which(chr23_hwe$TEST=="AFF" & chr23_hwe$P<=1E-12),]
    remove_4<-remove_4[which(!remove_4$SNP %in% remove_3$SNP),]
    
    print(paste("N variants fail hwe in cases: ",sum(nrow(remove_2),nrow(remove_4))),sep="")
    
    to_remove<-rbind(remove_1,remove_2,remove_3,remove_4)
    to_remove<-to_remove[,c(1:2)]
    to_remove<-to_remove[!duplicated(to_remove),]
    print(dim(to_remove))
    # print(table(to_remove$CHR))
    
  } else {
    print(paste("N variants fail hwe in cases: ",sum(nrow(remove_2))))
    to_remove<-rbind(remove_1,remove_2)
    to_remove<-to_remove[,c(1:2)]
    to_remove<-to_remove[!duplicated(to_remove),]
    print(dim(to_remove))
    # print(table(to_remove$CHR))
  }
  # print(table(dim(to_remove)==dim(to_remove[!duplicated(to_remove$SNP),])))
  write.table(to_remove[,"SNP",drop=F],paste(path,"pre_imputation/QC/",cohorts[j],"/list_var_exclude_hwe",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  rm(list=ls()[!ls() %in% c("path","j","cohorts")])
}

##################
# [1] "all_hce"
# [1] "N variants fail hwe in cases:  402"
# [1] 1682    2
# [1] "niddk_old_gwas"
# [1] "N variants fail hwe in cases:  42"
# [1] 200   2
# [1] "australia_omniexome"
# [1] "N variants fail hwe in cases:  6"
# [1] 712   2
# [1] "gwas1"
# [1] "N variants fail hwe in cases:  12"
# [1] 1364    2
# [1] "gwas2"
# [1] "N variants fail hwe in cases:  7"
# [1] 10466     2
# [1] "pittsburgh_gsa"
# [1] "N variants fail hwe in cases:  0"
# [1] 687   2
# [1] "spain_gsa"
# [1] "N variants fail hwe in cases:  21"
# [1] 1448    2
# [1] "italy_gsa"
# [1] "N variants fail hwe in cases:  14"
# [1] 222   2
# [1] "kiel_austria_sibdcs_gsa"
# [1] "N variants fail hwe in cases:  272"
# [1] 1482    2
# [1] "netherlands_gsa"
# [1] "N variants fail hwe in cases:  673"
# [1] 1205    2
# [1] "slovenia_gsa"
# [1] "N variants fail hwe in cases:  0"
# [1] 123   2
# [1] "sweden_gsa"
# [1] "N variants fail hwe in cases:  1"
# [1] 764   2
# [1] "niddk_broad_gsa"
# [1] "N variants fail hwe in cases:  329"
# [1] 1010    2
# [1] "niddk_feinstein_gsa"
# [1] "N variants fail hwe in cases:  118"
# [1] 1249    2
# [1] "basque_gsa"
# [1] "N variants fail hwe in cases:  2"
# [1] 583   2
# [1] "lithuania_gsa"
# [1] "N variants fail hwe in cases:  14"
# [1] 2357    2
# [1] "belgium_louis_gsa"
# [1] "N variants fail hwe in cases:  30"
# [1] 948   2
# [1] "belgium_franchimont_gsa"
# [1] "N variants fail hwe in cases:  44"
# [1] 1243    2
# [1] "belgium_vermeire_gsa"
# [1] "N variants fail hwe in cases:  693"
# [1] 2210    2
# [1] "prism_nfe_gsa"
# [1] "N variants fail hwe in cases:  174"
# [1] 178   2
# [1] "prism_nfe_gwas"
# [1] "N variants fail hwe in cases:  15"
# [1] 132   2
# [1] "finland_illugwas"
# [1] "N variants fail hwe in cases:  40"
# [1] 40  2
# [1] "german_affy6_old_gwas"
# [1] "N variants fail hwe in cases:  624"
# [1] 6903    2
# [1] "norway_affy6_old_gwas"
# [1] "N variants fail hwe in cases:  569"
# [1] 5924    2
# [1] "belgium_inf1_old_gwas"
# [1] "N variants fail hwe in cases:  0"
# [1] 0 2
# [1] "belgium_inf2_old_gwas"
# [1] "N variants fail hwe in cases:  0"
# [1] 0 2
# [1] "cedars_370k_old_gwas"
# [1] "N variants fail hwe in cases:  19"
# [1] 19  2
# [1] "cedars_610k_old_gwas"
# [1] "N variants fail hwe in cases:  98"
# [1] 98  2
# [1] "cedars_omni_old_gwas"
# [1] "N variants fail hwe in cases:  363"
# [1] 363   2
# [1] "swedish_uc_old_gwas"
# [1] "N variants fail hwe in cases:  1"
# [1] 7 2
# [1] "mccauley_gsa"
# [1] "N variants fail hwe in cases:  318"
# [1] 318   2
# [1] "ccfa_gsa"
# [1] "N variants fail hwe in cases:  1270"
# [1] 1270    2
# [1] "cedars_gsa"
# [1] "N variants fail hwe in cases:  217"
# [1] 1509    2
# [1] "bernstein_gsa"
# [1] "N variants fail hwe in cases:  254"
# [1] 254   2
# [1] "farkkila_gsa"
# [1] "N variants fail hwe in cases:  8"
# [1] 8 2
# [1] "franchimont_gsa"
# [1] "N variants fail hwe in cases:  17"
# [1] 2193    2
# [1] "franke_gsa"
# [1] "N variants fail hwe in cases:  2"
# [1] 762   2
# [1] "helmsley_prism_gsa"
# [1] "N variants fail hwe in cases:  14"
# [1] 130   2
# [1] "helmsley_xavier_prism_gsa"
# [1] "N variants fail hwe in cases:  110"
# [1] 143   2
# [1] "hyams_protect_gsa"
# [1] "N variants fail hwe in cases:  247"
# [1] 247   2
# [1] "lewis_sparc_gsa"
# [1] "N variants fail hwe in cases:  1530"
# [1] 1530    2
# [1] "mccauley_new_gsa"
# [1] "N variants fail hwe in cases:  394"
# [1] 580   2
# [1] "mcgovern_gsa"
# [1] "N variants fail hwe in cases:  807"
# [1] 2365    2
# [1] "moayyedi_imagine_gsa"
# [1] "N variants fail hwe in cases:  190"
# [1] 667   2
# [1] "newberry_share_gsa"
# [1] "N variants fail hwe in cases:  554"
# [1] 554   2
# [1] "niddk_cho_gsa"
# [1] "N variants fail hwe in cases:  18"
# [1] 970   2
# [1] "niddk_duerr_gsa"
# [1] "N variants fail hwe in cases:  29"
# [1] 1262    2
# [1] "niddk_rioux_gsa"
# [1] "N variants fail hwe in cases:  0"
# [1] 902   2
# [1] "niddk_silverberg_gsa"
# [1] "N variants fail hwe in cases:  244"
# [1] 998   2
# [1] "palotie_hus_gsa"
# [1] "N variants fail hwe in cases:  714"
# [1] 714   2
# [1] "pekow_share_gsa"
# [1] "N variants fail hwe in cases:  393"
# [1] 393   2
# [1] "rioux_igenomed_gsa"
# [1] "N variants fail hwe in cases:  74"
# [1] 74  2
# [1] "sands_msccr_gsa"
# [1] "N variants fail hwe in cases:  179"
# [1] 428   2
# [1] "stampfer_gsa"
# [1] "N variants fail hwe in cases:  1"
# [1] 1648    2
# [1] "vermeire_gsa"
# [1] "N variants fail hwe in cases:  1073"
# [1] 2593    2
# [1] "weersma_gsa"
# [1] "N variants fail hwe in cases:  0"
# [1] 669   2
# [1] "xavier_prism_gsa"
# [1] "N variants fail hwe in cases:  246"
# [1] 282   2
# [1] "xavier_share_gsa"
# [1] "N variants fail hwe in cases:  345"
# [1] 345   2
##################


########


## noDuplicates

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_46_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_46_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy \
--allow-no-sex \
--exclude ${path_gwas}/pre_imputation/QC/${i}/list_var_exclude_hwe \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy_interstudy_hwe"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_46_${i} | grep -E "completed"
done


## withDuplicates

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_47_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_47_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy \
--allow-no-sex \
--exclude ${path_gwas}/pre_imputation/QC/${i}/list_var_exclude_hwe \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_hwe"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_47_${i} | grep -E "completed"
done
