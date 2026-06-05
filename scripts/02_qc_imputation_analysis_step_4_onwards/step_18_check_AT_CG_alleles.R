# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
####################################
# 18.- DOUBLE CHEK A/T C/G ALLELES #
####################################

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do 
bsub -J"TOPMed" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_b38_alleles_1_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_b38_alleles_1_${i} \
"awk -v OFS='\t' '{print \$2,\$6}' ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim > \
${path_gwas}pre_imputation/QC/${i}/list_variants_${i}_hg38_posstrandaligned_with_A1"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_b38_alleles_1_${i} | grep -E "completed"
done

for i in ${studies[@]}
do 
bsub -J"TOPMed" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_b38_alleles_1.1_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_b38_alleles_1.1_${i} \
"sort ${path_gwas}pre_imputation/QC/${i}/list_variants_${i}_hg38_posstrandaligned_with_A1 | uniq \
> ${path_gwas}pre_imputation/QC/${i}/list_variants_${i}_hg38_posstrandaligned_with_A1_ed"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_b38_alleles_1.1_${i} | grep -E "completed"
done
for i in ${studies[@]}
do echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/list_variants_${i}_hg38_posstrandaligned_with_A1_ed
done


# EDIT CONTROL LIST, SO FID AND IID ARE THE SAME - note that some studies do not have control set, but empty file was created
for i in ${studies[@]}
do 
bsub -J"TOPMed" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_b38_alleles_2_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_b38_alleles_2_${i} \
"awk -v OFS='\t' '{print \$2,\$2}' ${path_gwas}pre_imputation/QC/${i}/list_controls > ${path_gwas}pre_imputation/QC/${i}/list_controls2"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_b38_alleles_2_${i} | grep -E "completed"
done


### EUR NON-JEWISH CTR:

for i in ${studies[@]}
do 
bsub -J"TOPMed" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_b38_alleles_3_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_b38_alleles_3_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed \
--keep-allele-order --allow-no-sex \
--keep ${path_gwas}pre_imputation/QC/pca_1000gp/list_eur_nonjewish_ancestry_samples_withDuplicates \
--make-bed --out ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_eur"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_b38_alleles_3_${i} | grep -E "completed"
done


for i in ${studies[@]}
do 
bsub -J"TOPMed" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_b38_alleles_4_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_b38_alleles_4_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_eur \
--keep-allele-order --allow-no-sex \
--keep ${path_gwas}pre_imputation/QC/${i}/list_controls2 \
--freq counts \
--out ${path_gwas}/pre_imputation/QC/${i}/all_hce_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_eur_ctr"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_b38_alleles_4_${i} | grep -E "completed"
done

for i in ${studies[@]}
do echo ${i} && wc -l ${path_gwas}pre_imputation/QC/${i}/list_controls2
done


# resubmit for those with cases only - keeping all samples (all cases)
studies_noctr=(finland_illugwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas mccauley_gsa ccfa_gsa farkkila_gsa helmsley_prism_gsa hyams_protect_gsa lewis_sparc_gsa newberry_share_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa)

for i in ${studies_noctr[@]}
do echo ${i} && wc -l ${path_gwas}pre_imputation/QC/${i}/list_controls2
done

for i in ${studies_noctr[@]}
do 
bsub -J"TOPMed" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_b38_alleles_4_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_b38_alleles_4_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_eur \
--keep-allele-order --allow-no-sex \
--freq counts \
--out ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_eur_ctr"
done


for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_b38_alleles_4_${i} | grep -E "completed"
done

for i in ${studies[@]}
do echo ${i} && ls -la ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_eur_ctr.frq.counts
done


# GNOMAD:

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=200

for i in ${studies[@]}
do
bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_b38_alleles_5_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_b38_alleles_5_${i} \
"awk 'NR==FNR{vals[\$2];next} (\$1) in vals' ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim \
<(zcat ${path_gwas}resources/gnomad/gnomad_freq_edited.gz) > ${path_gwas}pre_imputation/QC/${i}/${i}_gnomad_variants"
done

# CONTINUE HERE!!

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_b38_alleles_5_${i} | grep -E "completed"
done

all_hce
niddk_old_gwas
australia_omniexome

## /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"


cohorts<-c("australia_omniexome")

cohorts<-c("gwas1","gwas2","all_hce","niddk_old_gwas","australia_omniexome"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","swedish_uc_old_gwas"
           ,"mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

length(cohorts)
# [1] 58

for (j in 35:length(cohorts)) {
  
  print(cohorts[j])
  
  gnomad<-read.table(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_gnomad_variants",sep=""),head=F)
  
  colnames(gnomad)<-c("SNP","CHROM","POS","REF","ALT","AF","AF_nfe","AF_afr","AF_amr","AF_eas","AF_sas","AF_asj")
  gnomad$AF_nfe<-as.numeric(as.character(gnomad$AF_nfe))
  
  vfreq<-read.table(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_eur_ctr.frq.counts",sep=""),head=T)
  vfreq$freq_alt<-vfreq$C1/(vfreq$C1+vfreq$C2)
  vfreq$freq_alt<-as.numeric(vfreq$freq_alt)
  
  print(paste("N variants not present in gnomad:",nrow(vfreq[which(!vfreq$SNP %in% gnomad$SNP),]),sep=""))

  all<-merge(vfreq,gnomad,by="SNP",sort=F)
  
  print(table(all$A2,all$REF))
  print(table(all$A1,all$ALT))

  remove1<-all[which( (all$freq_alt>=0.45 & all$freq_alt<=0.55) & 
                        ( (all$A1=="G" & all$A2=="C") | (all$A1=="C" & all$A2=="G") | 
                            (all$A1=="A" & all$A2=="T") | (all$A1=="T" & all$A2=="A") ) ),]
  
  # remove A/T G/C that cannot be evaluated:
  remove2<-vfreq[which(!(vfreq$SNP %in% gnomad$SNP) & ( (vfreq$A1=="G" & vfreq$A2=="C") | (vfreq$A1=="C" & vfreq$A2=="G") | 
                                                          (vfreq$A1=="A" & vfreq$A2=="T") | (vfreq$A1=="T" & vfreq$A2=="A") ) ),]

  remove<-rbind(remove1[,"SNP",drop=F],remove2[,"SNP",drop=F])
  
  bim<-read.table(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim",sep=""),head=F)
  bim$V2<-as.character(bim$V2)
  dups<-bim[duplicated(bim$V2),"V2"]
  
  if(length(dups)>0){
    for (i in 1:length(dups)){
      bim$V2[which(bim$V2==dups[i])][1]<-paste(dups[i],"_tmp",sep="")
    }
    
    write.table(bim,paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_edited.bim",sep=""),
                col.names=F,row.names=F,quote=F,sep="\t")
    remove3<-bim[grep("_tmp",bim$V2),]
    colnames(remove3)[2]<-"SNP"
    remove<-rbind(remove1[,"SNP",drop=F],remove2[,"SNP",drop=F],remove3[,"SNP",drop=F])
    print(paste("study with N dups:",nrow(remove3),sep=""))
  } else {
    bim<-read.table(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim",sep=""),head=F)
    write.table(bim,paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_edited.bim",sep=""),
                col.names=F,row.names=F,quote=F,sep="\t")
    remove<-rbind(remove1[,"SNP",drop=F],remove2[,"SNP",drop=F])
  }

  print(paste("N variants to remove: ",nrow(remove),sep=""))
  
  flip<-all[which( (all$freq_alt<0.45 | all$freq_alt>0.55) & 
                     ( (all$A1=="G" & all$A2=="C") | (all$A1=="C" & all$A2=="G") | 
                         (all$A1=="A" & all$A2=="T") | (all$A1=="T" & all$A2=="A") ) ),]
  flip<-flip[which( (flip$freq_alt>0.5 & flip$AF_nfe<0.5) | (flip$freq_alt<0.5 & flip$AF_nfe>0.5) ),]
  print(paste("N variants to flip: ",nrow(flip),sep=""))
  
  write.table(remove[,"SNP",drop=F],paste(path,"pre_imputation/QC/",cohorts[j],"/list_variants_to_remove_AT_CG_",cohorts[j],sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  write.table(flip[,"SNP",drop=F],paste(path,"pre_imputation/QC/",cohorts[j],"/list_variants_to_flip_AT_CG_",cohorts[j],sep=""),col.names=F,row.names=F,quote=F,sep="\t")
 
  rm(list=ls()[!ls() %in% c("path","j","cohorts")]) 

}


############

# [1] "netherlands_gsa"
# [1] "N variants not present in gnomad:117"
# 
# A      C      G      T
# A 124044      0      0      0
# C      0 159496      0      0
# G      0      0 158819      0
# T      0      0      0 124223
# 
# A      C      G      T
# A 158718      0      0      0
# C      0 124872      0      0
# G      0      0 124371      0
# T      0      0      0 158621
# [1] "study with N dups:1"
# [1] "N variants to remove: 162"
# [1] "N variants to flip: 1973"
# [1] "slovenia_gsa"
# [1] "N variants not present in gnomad:39"
# 
# A      C      G      T
# A 111526      0      0      0
# C      0 143122      0      0
# G      0      0 142715      0
# T      0      0      0 111575
# 
# A      C      G      T
# A 142710      0      0      0
# C      0 111878      0      0
# G      0      0 111713      0
# T      0      0      0 142637
# [1] "study with N dups:1"
# [1] "N variants to remove: 159"
# [1] "N variants to flip: 1388"
# [1] "sweden_gsa"
# [1] "N variants not present in gnomad:97"
# 
# A      C      G      T
# A 117735      0      0      0
# C      0 151740      0      0
# G      0      0 151098      0
# T      0      0      0 118030
# 
# A      C      G      T
# A 151156      0      0      0
# C      0 118410      0      0
# G      0      0 118036      0
# T      0      0      0 151001
# [1] "study with N dups:1"
# [1] "N variants to remove: 166"
# [1] "N variants to flip: 1679"
# [1] "niddk_broad_gsa"
# [1] "N variants not present in gnomad:221"
# 
# A      C      G      T
# A 121725      0      0      0
# C      0 163809      0      0
# G      0      0 163005      0
# T      0      0      0 121963
# 
# A      C      G      T
# A 162807      0      0      0
# C      0 122637      0      0
# G      0      0 122207      0
# T      0      0      0 162851
# [1] "study with N dups:1"
# [1] "N variants to remove: 172"
# [1] "N variants to flip: 2169"
# [1] "niddk_feinstein_gsa"
# [1] "N variants not present in gnomad:255"
# 
# A      C      G      T
# A 124254      0      0      0
# C      0 167052      0      0
# G      0      0 166218      0
# T      0      0      0 124394
# 
# A      C      G      T
# A 165671      0      0      0
# C      0 125441      0      0
# G      0      0 124931      0
# T      0      0      0 165875
# [1] "N variants to remove: 164"
# [1] "N variants to flip: 2931"
# [1] "basque_gsa"
# [1] "N variants not present in gnomad:83"
# 
# A      C      G      T
# A 119220      0      0      0
# C      0 153299      0      0
# G      0      0 152686      0
# T      0      0      0 119163
# 
# A      C      G      T
# A 152622      0      0      0
# C      0 119730      0      0
# G      0      0 119689      0
# T      0      0      0 152327
# [1] "N variants to remove: 190"
# [1] "N variants to flip: 2153"
# [1] "prism_nfe_gsa"
# [1] "N variants not present in gnomad:87"
# 
# A      C      G      T
# A 118126      0      0      0
# C      0 152483      0      0
# G      0      0 152052      0
# T      0      0      0 118339
# 
# A      C      G      T
# A 151932      0      0      0
# C      0 118820      0      0
# G      0      0 118446      0
# T      0      0      0 151802
# [1] "study with N dups:1"
# [1] "N variants to remove: 166"
# [1] "N variants to flip: 1745"
# [1] "lithuania_gsa"
# [1] "N variants not present in gnomad:143"
# 
# A      C      G      T
# A 116522      0      0      0
# C      0 151743      0      0
# G      0      0 151082      0
# T      0      0      0 116767
# 
# A      C      G      T
# A 150962      0      0      0
# C      0 117124      0      0
# G      0      0 116865      0
# T      0      0      0 151163
# [1] "study with N dups:1"
# [1] "N variants to remove: 132"
# [1] "N variants to flip: 1655"
# [1] "belgium_louis_gsa"
# [1] "N variants not present in gnomad:194"
# 
# A      C      G      T
# A 117665      0      0      0
# C      0 154488      0      0
# G      0      0 153880      0
# T      0      0      0 117899
# 
# A      C      G      T
# A 153750      0      0      0
# C      0 118401      0      0
# G      0      0 118053      0
# T      0      0      0 153728
# [1] "study with N dups:1"
# [1] "N variants to remove: 109"
# [1] "N variants to flip: 1890"
# [1] "belgium_franchimont_gsa"
# [1] "N variants not present in gnomad:124"
# 
# A      C      G      T
# A 116324      0      0      0
# C      0 151777      0      0
# G      0      0 151304      0
# T      0      0      0 116526
# 
# A      C      G      T
# A 151129      0      0      0
# C      0 116976      0      0
# G      0      0 116713      0
# T      0      0      0 151113
# [1] "study with N dups:1"
# [1] "N variants to remove: 118"
# [1] "N variants to flip: 1761"
# [1] "belgium_vermeire_gsa"
# [1] "N variants not present in gnomad:1411"
# 
# A      C      G      T
# A 124972      0      0      0
# C      0 176891      0      0
# G      0      0 176455      0
# T      0      0      0 125250
# 
# A      C      G      T
# A 176016      0      0      0
# C      0 126162      0      0
# G      0      0 125731      0
# T      0      0      0 175659
# [1] "study with N dups:1"
# [1] "N variants to remove: 219"
# [1] "N variants to flip: 2622"
# [1] "prism_nfe_gwas"
# [1] "N variants not present in gnomad:4"
# 
# A     C     G     T
# A 58479     0     0     0
# C     0 61699     0     0
# G     0     0 61562     0
# T     0     0     0 58120
# 
# A     C     G     T
# A 61462     0     0     0
# C     0 58285     0     0
# G     0     0 58389     0
# T     0     0     0 61724
# [1] "N variants to remove: 14"
# [1] "N variants to flip: 170"
# [1] "finland_illugwas"
# [1] "N variants not present in gnomad:4"
# 
# A     C     G     T
# A 56769     0     0     0
# C     0 59986     0     0
# G     0     0 59879     0
# T     0     0     0 56479
# 
# A     C     G     T
# A 59786     0     0     0
# C     0 56620     0     0
# G     0     0 56708     0
# T     0     0     0 59999
# [1] "N variants to remove: 16"
# [1] "N variants to flip: 165"
# [1] "german_affy6_old_gwas"
# [1] "N variants not present in gnomad:17"
# 
# A      C      G      T
# A 178429      0      0      0
# C      0 207513      0      0
# G      0      0 207351      0
# T      0      0      0 178524
# 
# A      C      G      T
# A 197384      0      0      0
# C      0 188401      0      0
# G      0      0 188497      0
# T      0      0      0 197535
# [1] "N variants to remove: 4777"
# [1] "N variants to flip: 54759"
# [1] "norway_affy6_old_gwas"
# [1] "N variants not present in gnomad:13"
# 
# A      C      G      T
# A 166190      0      0      0
# C      0 190860      0      0
# G      0      0 190298      0
# T      0      0      0 166174
# 
# A      C      G      T
# A 181118      0      0      0
# C      0 175287      0      0
# G      0      0 175652      0
# T      0      0      0 181465
# [1] "N variants to remove: 4918"
# [1] "N variants to flip: 50115"
# [1] "belgium_inf1_old_gwas"
# [1] "N variants not present in gnomad:4"
# 
# A     C     G     T
# A 72347     0     0     0
# C     0 76888     0     0
# G     0     0 76904     0
# T     0     0     0 72187
# 
# A     C     G     T
# A 77057     0     0     0
# C     0 72204     0     0
# G     0     0 72330     0
# T     0     0     0 76735
# [1] "N variants to remove: 0"
# [1] "N variants to flip: 0"
# [1] "belgium_inf2_old_gwas"
# [1] "N variants not present in gnomad:4"
# 
# A     C     G     T
# A 69508     0     0     0
# C     0 73958     0     0
# G     0     0 74129     0
# T     0     0     0 69326
# 
# A     C     G     T
# A 74293     0     0     0
# C     0 69292     0     0
# G     0     0 69542     0
# T     0     0     0 73794
# [1] "N variants to remove: 0"
# [1] "N variants to flip: 0"
# [1] "cedars_370k_old_gwas"
# [1] "N variants not present in gnomad:6"
# 
# A     C     G     T
# A 80983     0     0     0
# C     0 86010     0     0
# G     0     0 85956     0
# T     0     0     0 80870
# 
# A     C     G     T
# A 85967     0     0     0
# C     0 81137     0     0
# G     0     0 80938     0
# T     0     0     0 85777
# [1] "N variants to remove: 68"
# [1] "N variants to flip: 642"
# [1] "cedars_610k_old_gwas"
# [1] "N variants not present in gnomad:12"
# 
# A      C      G      T
# A 137944      0      0      0
# C      0 149478      0      0
# G      0      0 149158      0
# T      0      0      0 138631
# 
# A      C      G      T
# A 149109      0      0      0
# C      0 139003      0      0
# G      0      0 138195      0
# T      0      0      0 148904
# [1] "N variants to remove: 159"
# [1] "N variants to flip: 1322"
# [1] "cedars_omni_old_gwas"
# [1] "N variants not present in gnomad:2948"
# 
# A      C      G      T
# A 145209      0      0      0
# C      0 166176      0      0
# G      0      0 158468      0
# T      0      0      0 145293
# 
# A      C      G      T
# A 158415      0      0      0
# C      0 145529      0      0
# G      0      0 145009      0
# T      0      0      0 166193
# [1] "N variants to remove: 12"
# [1] "N variants to flip: 31"
# [1] "swedish_uc_old_gwas"
# [1] "N variants not present in gnomad:4"
# 
# A     C     G     T
# A 71064     0     0     0
# C     0 75471     0     0
# G     0     0 75611     0
# T     0     0     0 71160
# 
# A     C     G     T
# A 75745     0     0     0
# C     0 71199     0     0
# G     0     0 71025     0
# T     0     0     0 75337
# [1] "N variants to remove: 0"
# [1] "N variants to flip: 0"
# [1] "mccauley_gsa"
# [1] "N variants not present in gnomad:117"
# 
# A      C      G      T
# A 117571      0      0      0
# C      0 151875      0      0
# G      0      0 151398      0
# T      0      0      0 117597
# 
# A      C      G      T
# A 151302      0      0      0
# C      0 118089      0      0
# G      0      0 117874      0
# T      0      0      0 151176
# [1] "study with N dups:1"
# [1] "N variants to remove: 135"
# [1] "N variants to flip: 1696"
# [1] "ccfa_gsa"
# [1] "N variants not present in gnomad:224"
# 
# A      C      G      T
# A 124886      0      0      0
# C      0 164842      0      0
# G      0      0 164237      0
# T      0      0      0 125204
# 
# A      C      G      T
# A 164115      0      0      0
# C      0 125810      0      0
# G      0      0 125364      0
# T      0      0      0 163880
# [1] "study with N dups:1"
# [1] "N variants to remove: 134"
# [1] "N variants to flip: 2122"
# [1] "cedars_gsa"
# [1] "N variants not present in gnomad:301"
# 
# A      C      G      T
# A 126446      0      0      0
# C      0 168018      0      0
# G      0      0 167336      0
# T      0      0      0 126701
# 
# A      C      G      T
# A 167095      0      0      0
# C      0 127433      0      0
# G      0      0 126861      0
# T      0      0      0 167112
# [1] "study with N dups:1"
# [1] "N variants to remove: 142"
# [1] "N variants to flip: 2118"
# [1] "bernstein_gsa"
# [1] "N variants not present in gnomad:53"
# 
# A      C      G      T
# A 117523      0      0      0
# C      0 150903      0      0
# G      0      0 150188      0
# T      0      0      0 117636
# 
# A      C      G      T
# A 150148      0      0      0
# C      0 118051      0      0
# G      0      0 117776      0
# T      0      0      0 150275
# [1] "study with N dups:1"
# [1] "N variants to remove: 91"
# [1] "N variants to flip: 1636"
# [1] "farkkila_gsa"
# [1] "N variants not present in gnomad:31"
# 
# A      C      G      T
# A 106310      0      0      0
# C      0 133708      0      0
# G      0      0 133171      0
# T      0      0      0 106588
# 
# A      C      G      T
# A 133111      0      0      0
# C      0 106902      0      0
# G      0      0 106471      0
# T      0      0      0 133293
# [1] "study with N dups:1"
# [1] "N variants to remove: 110"
# [1] "N variants to flip: 1277"
# [1] "franchimont_gsa"
# [1] "N variants not present in gnomad:223"
# 
# A      C      G      T
# A 119112      0      0      0
# C      0 157830      0      0
# G      0      0 157131      0
# T      0      0      0 119273
# 
# A      C      G      T
# A 156969      0      0      0
# C      0 119873      0      0
# G      0      0 119549      0
# T      0      0      0 156955
# [1] "study with N dups:1"
# [1] "N variants to remove: 126"
# [1] "N variants to flip: 1950"
# [1] "franke_gsa"
# [1] "N variants not present in gnomad:99"
# 
# A      C      G      T
# A 117876      0      0      0
# C      0 152119      0      0
# G      0      0 151505      0
# T      0      0      0 118105
# 
# A      C      G      T
# A 151498      0      0      0
# C      0 118477      0      0
# G      0      0 118229      0
# T      0      0      0 151401
# [1] "study with N dups:1"
# [1] "N variants to remove: 123"
# [1] "N variants to flip: 1646"
# [1] "helmsley_prism_gsa"
# [1] "N variants not present in gnomad:4"
# 
# A     C     G     T
# A 58538     0     0     0
# C     0 61719     0     0
# G     0     0 61605     0
# T     0     0     0 58179
# 
# A     C     G     T
# A 61498     0     0     0
# C     0 58354     0     0
# G     0     0 58437     0
# T     0     0     0 61752
# [1] "N variants to remove: 13"
# [1] "N variants to flip: 171"
# [1] "helmsley_xavier_prism_gsa"
# [1] "N variants not present in gnomad:4"
# 
# A     C     G     T
# A 58549     0     0     0
# C     0 61737     0     0
# G     0     0 61617     0
# T     0     0     0 58179
# 
# A     C     G     T
# A 61517     0     0     0
# C     0 58355     0     0
# G     0     0 58454     0
# T     0     0     0 61756
# [1] "N variants to remove: 20"
# [1] "N variants to flip: 169"
# [1] "hyams_protect_gsa"
# [1] "N variants not present in gnomad:117"
# 
# A      C      G      T
# A 120145      0      0      0
# C      0 155497      0      0
# G      0      0 154904      0
# T      0      0      0 120219
# 
# A      C      G      T
# A 154845      0      0      0
# C      0 120723      0      0
# G      0      0 120390      0
# T      0      0      0 154807
# [1] "study with N dups:1"
# [1] "N variants to remove: 144"
# [1] "N variants to flip: 1782"
# [1] "lewis_sparc_gsa"
# [1] "N variants not present in gnomad:348"
# 
# A      C      G      T
# A 125713      0      0      0
# C      0 166903      0      0
# G      0      0 166370      0
# T      0      0      0 125973
# 
# A      C      G      T
# A 166206      0      0      0
# C      0 126658      0      0
# G      0      0 126204      0
# T      0      0      0 165891
# [1] "study with N dups:1"
# [1] "N variants to remove: 134"
# [1] "N variants to flip: 2208"
# [1] "mccauley_new_gsa"
# [1] "N variants not present in gnomad:185"
# 
# A      C      G      T
# A 121074      0      0      0
# C      0 158491      0      0
# G      0      0 158222      0
# T      0      0      0 121088
# 
# A      C      G      T
# A 157952      0      0      0
# C      0 121736      0      0
# G      0      0 121478      0
# T      0      0      0 157709
# [1] "study with N dups:1"
# [1] "N variants to remove: 133"
# [1] "N variants to flip: 1999"
# [1] "mcgovern_gsa"
# [1] "N variants not present in gnomad:513"
# 
# A      C      G      T
# A 128619      0      0      0
# C      0 173083      0      0
# G      0      0 172418      0
# T      0      0      0 128754
# 
# A      C      G      T
# A 172157      0      0      0
# C      0 129618      0      0
# G      0      0 129094      0
# T      0      0      0 172005
# [1] "study with N dups:1"
# [1] "N variants to remove: 153"
# [1] "N variants to flip: 2320"
# [1] "moayyedi_imagine_gsa"
# [1] "N variants not present in gnomad:115"
# 
# A      C      G      T
# A 124288      0      0      0
# C      0 160882      0      0
# G      0      0 160380      0
# T      0      0      0 124490
# 
# A      C      G      T
# A 160323      0      0      0
# C      0 124968      0      0
# G      0      0 124676      0
# T      0      0      0 160073
# [1] "study with N dups:1"
# [1] "N variants to remove: 126"
# [1] "N variants to flip: 1859"
# [1] "newberry_share_gsa"
# [1] "N variants not present in gnomad:91"
# 
# A      C      G      T
# A 121370      0      0      0
# C      0 157502      0      0
# G      0      0 156988      0
# T      0      0      0 121498
# 
# A      C      G      T
# A 156820      0      0      0
# C      0 122066      0      0
# G      0      0 121729      0
# T      0      0      0 156743
# [1] "study with N dups:1"
# [1] "N variants to remove: 126"
# [1] "N variants to flip: 1926"
# [1] "niddk_cho_gsa"
# [1] "N variants not present in gnomad:248"
# 
# A      C      G      T
# A 123409      0      0      0
# C      0 162918      0      0
# G      0      0 162395      0
# T      0      0      0 123544
# 
# A      C      G      T
# A 162124      0      0      0
# C      0 124161      0      0
# G      0      0 123797      0
# T      0      0      0 162184
# [1] "study with N dups:1"
# [1] "N variants to remove: 144"
# [1] "N variants to flip: 2021"
# [1] "niddk_duerr_gsa"
# [1] "N variants not present in gnomad:271"
# 
# A      C      G      T
# A 121081      0      0      0
# C      0 159433      0      0
# G      0      0 158850      0
# T      0      0      0 121308
# 
# A      C      G      T
# A 158710      0      0      0
# C      0 121811      0      0
# G      0      0 121544      0
# T      0      0      0 158607
# [1] "study with N dups:1"
# [1] "N variants to remove: 120"
# [1] "N variants to flip: 1993"
# [1] "niddk_rioux_gsa"
# [1] "N variants not present in gnomad:151"
# 
# A      C      G      T
# A 119752      0      0      0
# C      0 155971      0      0
# G      0      0 155602      0
# T      0      0      0 119908
# 
# A      C      G      T
# A 155428      0      0      0
# C      0 120432      0      0
# G      0      0 120034      0
# T      0      0      0 155339
# [1] "study with N dups:1"
# [1] "N variants to remove: 112"
# [1] "N variants to flip: 1796"
# [1] "niddk_silverberg_gsa"
# [1] "N variants not present in gnomad:291"
# 
# A      C      G      T
# A 124716      0      0      0
# C      0 165128      0      0
# G      0      0 164612      0
# T      0      0      0 124823
# 
# A      C      G      T
# A 164386      0      0      0
# C      0 125419      0      0
# G      0      0 125153      0
# T      0      0      0 164321
# [1] "study with N dups:1"
# [1] "N variants to remove: 142"
# [1] "N variants to flip: 2053"
# [1] "palotie_hus_gsa"
# [1] "N variants not present in gnomad:99"
# 
# A      C      G      T
# A 117240      0      0      0
# C      0 150699      0      0
# G      0      0 150025      0
# T      0      0      0 117397
# 
# A      C      G      T
# A 149975      0      0      0
# C      0 117776      0      0
# G      0      0 117535      0
# T      0      0      0 150075
# [1] "study with N dups:1"
# [1] "N variants to remove: 108"
# [1] "N variants to flip: 1594"
# [1] "pekow_share_gsa"
# [1] "N variants not present in gnomad:80"
# 
# A      C      G      T
# A 120866      0      0      0
# C      0 156561      0      0
# G      0      0 155946      0
# T      0      0      0 121133
# 
# A      C      G      T
# A 155861      0      0      0
# C      0 121670      0      0
# G      0      0 121186      0
# T      0      0      0 155789
# [1] "study with N dups:1"
# [1] "N variants to remove: 133"
# [1] "N variants to flip: 1861"
# [1] "rioux_igenomed_gsa"
# [1] "N variants not present in gnomad:52"
# 
# A      C      G      T
# A 115924      0      0      0
# C      0 147481      0      0
# G      0      0 146882      0
# T      0      0      0 116004
# 
# A      C      G      T
# A 146893      0      0      0
# C      0 116367      0      0
# G      0      0 116112      0
# T      0      0      0 146919
# [1] "study with N dups:1"
# [1] "N variants to remove: 109"
# [1] "N variants to flip: 1506"
# [1] "sands_msccr_gsa"
# [1] "N variants not present in gnomad:130"
# 
# A      C      G      T
# A 124805      0      0      0
# C      0 163102      0      0
# G      0      0 162493      0
# T      0      0      0 124933
# 
# A      C      G      T
# A 162310      0      0      0
# C      0 125618      0      0
# G      0      0 125145      0
# T      0      0      0 162260
# [1] "study with N dups:1"
# [1] "N variants to remove: 148"
# [1] "N variants to flip: 2023"
# [1] "stampfer_gsa"
# [1] "N variants not present in gnomad:291"
# 
# A      C      G      T
# A 118937      0      0      0
# C      0 156692      0      0
# G      0      0 156248      0
# T      0      0      0 118964
# 
# A      C      G      T
# A 156020      0      0      0
# C      0 119565      0      0
# G      0      0 119296      0
# T      0      0      0 155960
# [1] "study with N dups:1"
# [1] "N variants to remove: 128"
# [1] "N variants to flip: 1894"
# [1] "vermeire_gsa"
# [1] "N variants not present in gnomad:333"
# 
# A      C      G      T
# A 120878      0      0      0
# C      0 161014      0      0
# G      0      0 160610      0
# T      0      0      0 121038
# 
# A      C      G      T
# A 160336      0      0      0
# C      0 121700      0      0
# G      0      0 121305      0
# T      0      0      0 160199
# [1] "study with N dups:1"
# [1] "N variants to remove: 119"
# [1] "N variants to flip: 1986"
# [1] "weersma_gsa"
# [1] "N variants not present in gnomad:86"
# 
# A      C      G      T
# A 118929      0      0      0
# C      0 152694      0      0
# G      0      0 151958      0
# T      0      0      0 119136
# 
# A      C      G      T
# A 152021      0      0      0
# C      0 119593      0      0
# G      0      0 119198      0
# T      0      0      0 151905
# [1] "study with N dups:1"
# [1] "N variants to remove: 123"
# [1] "N variants to flip: 1673"
# [1] "xavier_prism_gsa"
# [1] "N variants not present in gnomad:99"
# 
# A      C      G      T
# A 119861      0      0      0
# C      0 155293      0      0
# G      0      0 154856      0
# T      0      0      0 120006
# 
# A      C      G      T
# A 154691      0      0      0
# C      0 120522      0      0
# G      0      0 120181      0
# T      0      0      0 154622
# [1] "study with N dups:1"
# [1] "N variants to remove: 129"
# [1] "N variants to flip: 1775"
# [1] "xavier_share_gsa"
# [1] "N variants not present in gnomad:167"
# 
# A      C      G      T
# A 119534      0      0      0
# C      0 155323      0      0
# G      0      0 154820      0
# T      0      0      0 119780
# 
# A      C      G      T
# A 154609      0      0      0
# C      0 120306      0      0
# G      0      0 119866      0
# T      0      0      0 154676
# [1] "study with N dups:1"
# [1] "N variants to remove: 113"
# [1] "N variants to flip: 1813"

############


##########################################################################################
# 18.2 REMOVE VARIANTS WE CANNOT BE SURE ARE IN THE RIGHT STRAND, AND FLIP THE OTHERS

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500

for i in ${studies[@]}
do
bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_b38_alleles_6_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_b38_alleles_6_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bed \
--bim ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_edited.bim \
--fam ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.fam \
--allow-no-sex \
--exclude ${path_gwas}pre_imputation/QC/${i}/list_variants_to_remove_AT_CG_${i} \
--flip ${path_gwas}pre_imputation/QC/${i}/list_variants_to_flip_AT_CG_${i} \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_b38_alleles_6_${i} | grep -E "completed"
done


for i in ${studies[@]}
do
bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_b38_alleles_7_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_b38_alleles_7_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip \
--allow-no-sex \
--a2-allele ${path_gwas}pre_imputation/QC/${i}/list_variants_${i}_hg38_posstrandaligned_with_A1_ed \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2"
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_b38_alleles_7_${i} | grep -E "completed"
done

