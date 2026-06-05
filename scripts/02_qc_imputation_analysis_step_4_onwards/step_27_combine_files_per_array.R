# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
######################################################################
# 0.- CREATE BED/BIM/FAM FILES WITH THE VARIANTS USED FOR IMPUTATION #
######################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
  
studies=(bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa 
         lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa 
         niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa 
         vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)
ancestry=(eur noneur)

MEM=250

for j in ${ancestry[@]}
do
for i in ${studies[@]}
do
bsub -J"pl" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stdout_post_server_qc_pre_regenie_1_${j}_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_post_server_qc_pre_regenie_1_${j}_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_${j}_nomonom \
--allow-no-sex \
--keep-allele-order \
--flip ${path_gwas}imputed/${i}/qc/2022/list_variants_negative_EmpR_toflip \
--exclude ${path_gwas}imputed/${i}/qc/2022/list_variants_to_exclude_2 \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_preimp_${j}"
done
done

# do the same with old studies pointing to old flip file
studies=(niddk_old_gwas australia_omniexome gwas2 pittsburgh_gsa
         kiel_austria_sibdcs_gsa netherlands_gsa 
         niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa
         belgium_vermeire_gsa
         prism_nfe_gsa prism_nfe_gwas finland_illugwas
         belgium_inf1_old_gwas belgium_inf2_old_gwas
         cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas
         swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa)
for j in ${ancestry[@]}
do
for i in ${studies[@]}
do
bsub -J"pl" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stdout_post_server_qc_pre_regenie_1_${j}_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_post_server_qc_pre_regenie_1_${j}_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_${j}_nomonom \
--allow-no-sex \
--keep-allele-order \
--flip ${path_gwas}imputed/${i}/qc/list_variants_negative_EmpR_toflip \
--exclude ${path_gwas}imputed/${i}/qc/2022/list_variants_to_exclude_2 \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_preimp_${j}"
done
done


# subset of studies with an final list of exclusions 
studies=(german_affy6_old_gwas italy_gsa all_hce slovenia_gsa gwas1 norway_affy6_old_gwas belgium_louis_gsa belgium_franchimont_gsa)
for j in ${ancestry[@]}
do
for i in ${studies[@]}
do
bsub -J"pl" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_post_server_qc_pre_regenie_1_${j}_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_post_server_qc_pre_regenie_1_${j}_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_${j}_nomonom \
--allow-no-sex \
--keep-allele-order \
--flip ${path_gwas}imputed/${i}/qc/list_variants_negative_EmpR_toflip \
--exclude ${path_gwas}imputed/${i}/qc/2022/list_variants_to_exclude_3 \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_preimp_${j}"
done
done


studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa 
         kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa 
         basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa 
         prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas 
         belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas 
         mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa 
         helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa 
         moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa 
         niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa 
         stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)
for i in ${studies[@]}
do echo ${i} && for j in ${ancestry[@]}
do echo ${j} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_post_server_qc_pre_regenie_1_${j}_${i} | grep "Successfully"
done
done



##########################################
# 1.- GROUP STUDIES BY GENOTYPING ARRRAY #
##########################################


#  1.1 CREATE ONE FOLDER PER GENOTYPING ARRAY

cd /path/to/ibdgwas/IIBDGC/
mkdir /path/to/ibdgwas/IIBDGC/post_imputation/2022/
  
cd /path/to/ibdgwas/IIBDGC/post_imputation/2022/
array=(illumina370 illumina550 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 gsa illuminaexome)
for i in ${array[@]}
do mkdir -p ${i}/genotyped_data/ && mkdir -p ${i}/imputed_data/ && mkdir -p ${i}/phenotype_data/
done


#  1.2 PER ARRAY, IDENTIFY SHARED VARIANTS AND MERGE FILES

# /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

array<-c("illumina370","illumina550","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")

illumina370<-c("belgium_inf1_old_gwas","belgium_inf2_old_gwas","swedish_uc_old_gwas","niddk_old_gwas")
affymetrix6<-c("german_affy6_old_gwas","norway_affy6_old_gwas","gwas2")
humanomniexpress<-c("australia_omniexome")
affymetrix500<-c("gwas1")
humancoreexome<-c("all_hce")
humanomni1<-c("pittsburgh_gsa")
quad610<-c("spain_gsa")
gsa<-c("italy_gsa","kiel_austria_sibdcs_gsa"
       ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
       ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
       ,"mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa","franke_gsa","hyams_protect_gsa",
       "lewis_sparc_gsa","mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
       "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
       "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
       "xavier_share_gsa")

illuminaexome<-c("prism_nfe_gwas","helmsley_prism_gsa","helmsley_xavier_prism_gsa","finland_illugwas")

array<-c("illuminaexome")


ancestry<-c("eur")

for (i in 1:length(array)){
  
  print(array[i])
  
  study<-get(array[i])
  print(study)
  
  for (ii in 1:length(ancestry)) {
    
    if (length(study)>1) {
      
      # slovenia does not have non-eur samples
      if(array[i]=="gsa" & ancestry[ii]=="noneur") {
        study<-study[which(study!="slovenia_gsa")]
      }
      
      print(ancestry[ii])
      for (j in 1:length(study)){
        
        # print(study[j])
        if(j==1) {
          bim<-read.table(paste(path,"pre_imputation/QC/",study[j],"/",study[j],"_postqc_preimp_",ancestry[ii],".bim",sep=""),head=F)
        } else {
          tmp<-read.table(paste(path,"pre_imputation/QC/",study[j],"/",study[j],"_postqc_preimp_",ancestry[ii],".bim",sep=""),head=F)
          bim<-bim[which(bim$V2 %in% tmp$V2),]
        }
      }
      
      write.table(bim[,2,drop=F],paste(path,"post_imputation/2022/",array[i],"/genotyped_data/list_variants_shared_across_studies_",ancestry[ii],sep=""),col.names=F,row.names=F,quote=F,sep="\t")
      print(paste("N shared variants ",array[i]," ",ancestry[ii],": ",nrow(bim),sep=""))
      rm(bim,tmp)
      
      # extract the shared variants
      for (j in 1:length(study)){
        system(paste("/path/to/software/plink_linux_x86_64_20181202/./plink --bfile ",path,"pre_imputation/QC/",study[j],"/",study[j],"_postqc_preimp_",ancestry[ii]," --allow-no-sex --extract ",path,"post_imputation/2022/",array[i],"/genotyped_data/list_variants_shared_across_studies_",ancestry[ii]," --make-bed --out ",path,"post_imputation/2022/",array[i],"/genotyped_data/",study[j],"_postqc_preimp_",ancestry[ii],"_shared_variants",sep=""),ignore.stdout = T,ignore.stderr = T)
        # system(paste("/path/to/software/plink_linux_x86_64_20181202/./plink --bfile ",path,"pre_imputation/QC/",study[j],"/",study[j],"_postqc_preimp_",ancestry[ii]," --allow-no-sex --extract ",path,"post_imputation/2022/",array[i],"/genotyped_data/list_variants_shared_across_studies_",ancestry[ii]," --make-bed --out ",path,"post_imputation/2022/",array[i],"/genotyped_data/",study[j],"_postqc_preimp_",ancestry[ii],"_shared_variants",sep=""))
      }
      
      
      # merge the studies
      study_tmp<-study[-1]
      dat<-matrix(ncol=3,nrow=length(study_tmp))
      dat<-as.data.frame(dat)
      
      for (jj in 1:length(study_tmp)){
        dat[jj,1]<-paste(path,"post_imputation/2022/",array[i],"/genotyped_data/",study_tmp[jj],"_postqc_preimp_",ancestry[ii],"_shared_variants.bed",sep="")
        dat[jj,2]<-paste(path,"post_imputation/2022/",array[i],"/genotyped_data/",study_tmp[jj],"_postqc_preimp_",ancestry[ii],"_shared_variants.bim",sep="")
        dat[jj,3]<-paste(path,"post_imputation/2022/",array[i],"/genotyped_data/",study_tmp[jj],"_postqc_preimp_",ancestry[ii],"_shared_variants.fam",sep="")
      }
      
      write.table(dat,paste(path,"post_imputation/2022/",array[i],"/genotyped_data/list_studies_to_merge_",ancestry[ii],sep=""),col.names=F,row.names=F,quote=F,sep="\t")
      
      system(paste("/path/to/software/plink_linux_x86_64_20181202/./plink --bfile ",path,"post_imputation/2022/",array[i],"/genotyped_data/",study[1],"_postqc_preimp_",ancestry[ii],"_shared_variants --allow-no-sex --merge-list ",path,"post_imputation/2022/",array[i],"/genotyped_data/list_studies_to_merge_",ancestry[ii]," --make-bed --out ",path,"post_imputation/2022/",array[i],"/genotyped_data/",array[i],"_all_studies_merged_",ancestry[ii],sep=""))
      
      # remove the intermediate files:
      for (j in 1:length(study)){
        system(paste("rm ",path,"post_imputation/2022/",array[i],"/genotyped_data/",study[j],"_postqc_preimp_",ancestry[ii],"_shared_variants.*",sep=""))
      }
      
    } else {
      system(paste("/path/to/software/plink_linux_x86_64_20181202/./plink --bfile ",path,"pre_imputation/QC/",study[1],"/",study[1],"_postqc_preimp_",ancestry[ii]," --allow-no-sex --make-bed --out ",path,"post_imputation/2022/",array[i],"/genotyped_data/",array[i],"_all_studies_merged_",ancestry[ii],sep=""))
    }
    
  }
  
}



##########################################################################################################
# 1.3 PRUNE THE DATA BY LD and MAF AHEAD OF RUNNING FIRST STAGE OF SAIGE:

# Note: regenie will throw an error if a low-variance SNP is included in the step 1 run. Hence, the user should run adequate 
# QC filtering prior to running regenie to identify and remove such SNPs.


# VIF is 1/(1-R^2)
# 1/(1-0.9)
# [1] 10

array=(illumina370 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 gsa illuminaexome)

for i in ${array[@]}
do 
# plink2 does not have --indep fx inplemented:
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}post_imputation/2022/${i}/genotyped_data/${i}_all_studies_merged_eur \
--keep-allele-order \
--maf 0.01 \
--make-bed --out ${path_gwas}post_imputation/2022/${i}/genotyped_data/${i}_all_studies_merged_eur_maf
done

for i in ${array[@]}
do 
# plink2 does not have --indep fx inplemented:
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}post_imputation/2022/${i}/genotyped_data/${i}_all_studies_merged_eur_maf \
--keep-allele-order \
--indep 50 5 10 \
--out ${path_gwas}post_imputation/2022/${i}/genotyped_data/${i}_all_studies_merged_eur
done

for i in ${array[@]}
do 
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}post_imputation/2022/${i}/genotyped_data/${i}_all_studies_merged_eur_maf \
--keep-allele-order \
--extract ${path_gwas}post_imputation/2022/${i}/genotyped_data/${i}_all_studies_merged_eur.prune.in \
--make-bed --out ${path_gwas}post_imputation/2022/${i}/genotyped_data/${i}_all_studies_merged_eur_pruned
done

# remove intermediate file:
for i in ${array[@]}
do 
rm ${path_gwas}post_imputation/2022/${i}/genotyped_data/${i}_all_studies_merged_eur_maf.bed 
done

for i in ${array[@]}
do 
rm ${path_gwas}post_imputation/2022/${i}/genotyped_data/${i}_all_studies_merged_eur_maf.bim 
done

for i in ${array[@]}
do 
rm ${path_gwas}post_imputation/2022/${i}/genotyped_data/${i}_all_studies_merged_eur_maf.fam 
done

for i in ${array[@]}
do 
rm ${path_gwas}post_imputation/2022/${i}/genotyped_data/${i}_all_studies_merged_eur_maf.log 
done

wc -l /path/to/ibdgwas/IIBDGC/post_imputation/2022/*/genotyped_data/*_all_studies_merged_eur_pruned.bim
# 192361 /path/to/ibdgwas/IIBDGC/post_imputation/2022/affymetrix500/genotyped_data/affymetrix500_all_studies_merged_eur_pruned.bim
# 264863 /path/to/ibdgwas/IIBDGC/post_imputation/2022/affymetrix6/genotyped_data/affymetrix6_all_studies_merged_eur_pruned.bim
# 263372 /path/to/ibdgwas/IIBDGC/post_imputation/2022/gsa/genotyped_data/gsa_all_studies_merged_eur_pruned.bim
# 161991 /path/to/ibdgwas/IIBDGC/post_imputation/2022/humancoreexome/genotyped_data/humancoreexome_all_studies_merged_eur_pruned.bim
# 352258 /path/to/ibdgwas/IIBDGC/post_imputation/2022/humanomni1/genotyped_data/humanomni1_all_studies_merged_eur_pruned.bim
# 312061 /path/to/ibdgwas/IIBDGC/post_imputation/2022/humanomniexpress/genotyped_data/humanomniexpress_all_studies_merged_eur_pruned.bim
# 191680 /path/to/ibdgwas/IIBDGC/post_imputation/2022/illumina370/genotyped_data/illumina370_all_studies_merged_eur_pruned.bim
# 166340 /path/to/ibdgwas/IIBDGC/post_imputation/2022/illuminaexome/genotyped_data/illuminaexome_all_studies_merged_eur_pruned.bim
# 298845 /path/to/ibdgwas/IIBDGC/post_imputation/2022/quad610/genotyped_data/quad610_all_studies_merged_eur_pruned.bim

wc -l /path/to/ibdgwas/IIBDGC/post_imputation/2022/*/genotyped_data/*_all_studies_merged_eur_pruned.fam
# 4652 /path/to/ibdgwas/IIBDGC/post_imputation/2022/affymetrix500/genotyped_data/affymetrix500_all_studies_merged_eur_pruned.fam
# 11087 /path/to/ibdgwas/IIBDGC/post_imputation/2022/affymetrix6/genotyped_data/affymetrix6_all_studies_merged_eur_pruned.fam
# 81334 /path/to/ibdgwas/IIBDGC/post_imputation/2022/gsa/genotyped_data/gsa_all_studies_merged_eur_pruned.fam
# 22588 /path/to/ibdgwas/IIBDGC/post_imputation/2022/humancoreexome/genotyped_data/humancoreexome_all_studies_merged_eur_pruned.fam
# 2701 /path/to/ibdgwas/IIBDGC/post_imputation/2022/humanomni1/genotyped_data/humanomni1_all_studies_merged_eur_pruned.fam
# 1245 /path/to/ibdgwas/IIBDGC/post_imputation/2022/humanomniexpress/genotyped_data/humanomniexpress_all_studies_merged_eur_pruned.fam
# 5633 /path/to/ibdgwas/IIBDGC/post_imputation/2022/illumina370/genotyped_data/illumina370_all_studies_merged_eur_pruned.fam
# 2997 /path/to/ibdgwas/IIBDGC/post_imputation/2022/illuminaexome/genotyped_data/illuminaexome_all_studies_merged_eur_pruned.fam
# 3396 /path/to/ibdgwas/IIBDGC/post_imputation/2022/quad610/genotyped_data/quad610_all_studies_merged_eur_pruned.fam
# 135633 total


# create exclusion criteria for SNPs in the ashkenazi samples datasets in gsa and illumina370 (only two arrays with enough N - see step 30)

path_gwas=/path/to/ibdgwas/IIBDGC/
array=(illumina370 gsa)
j=eur_jewish

for i in ${array[@]}
do 
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}post_imputation/2022/${i}/genotyped_data/${i}_all_studies_merged_eur_pruned \
--keep ${path_gwas}pre_imputation/QC/pca_1000gp/list_${j}_ancestry_samples_noDuplicates_perstudy_interstudy_2 \
--maf 0.01 --write-snplist --out ${path_gwas}post_imputation/2022/${i}/genotyped_data/${i}_all_studies_merged_eur_pruned_${j}_maf_pass
done

for i in ${array[@]}
do echo ${i} && \
wc -l ${path_gwas}post_imputation/2022/${i}/genotyped_data/${i}_all_studies_merged_eur_pruned_${j}_maf_pass.snplist && \
wc -l ${path_gwas}post_imputation/2022/${i}/genotyped_data/${i}_all_studies_merged_eur_pruned.bim
done

# illumina370
# 191270 /path/to/ibdgwas/IIBDGC/post_imputation/2022/illumina370/genotyped_data/illumina370_all_studies_merged_eur_pruned_eur_jewish_maf_pass.snplist
# 191680 /path/to/ibdgwas/IIBDGC/post_imputation/2022/illumina370/genotyped_data/illumina370_all_studies_merged_eur_pruned.bim
# 
# gsa
# 237395 /path/to/ibdgwas/IIBDGC/post_imputation/2022/gsa/genotyped_data/gsa_all_studies_merged_eur_pruned_eur_jewish_maf_pass.snplist
# 263372 /path/to/ibdgwas/IIBDGC/post_imputation/2022/gsa/genotyped_data/gsa_all_studies_merged_eur_pruned.bim





########################################################################################################################################
########################################################################################################################################

#######################################
# 2.- MERGE IMPUTED FILES - AUTOSOMAL #
#######################################

##################################
# 2.1 - ONLY ONE STUDY PER ARRAY #
##################################

# 2.1.1 exclude list of variants HWE

humanomniexpress=(australia_omniexome)
array=humanomniexpress

humancoreexome=(all_hce)
array=humancoreexome

humanomni1=(pittsburgh_gsa)
array=humanomni1

quad610=(spain_gsa)
array=quad610

affymetrix500=(gwas1)
array=affymetrix500

# for i in ${affymetrix500[@]}
# for i in ${quad610[@]}
# for i in ${humanomniexpress[@]}
# for i in ${humancoreexome[@]}
# for i in ${humanomni1[@]}
do for chr in {1..22}
do
bsub -J"reg_ibd_test" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -n 2 -q basement \
-e ${path_gwas}post_imputation/2022/log/${i}_${chr}_subset_postimp_stderr \
-o ${path_gwas}post_imputation/2022/log/${i}_${chr}_subset_postimp_stdout \
"/path/to/software/qctool_v2.2.0/./qctool \
-g ${path_gwas}imputed/${i}/2022/eur/chr${chr}.dose.vcf.gz \
-filetype vcf -vcf-genotype-field GP \
-excl-rsids ${path_gwas}imputed/${i}/2022/eur/list_variants_post_imputation_hwe_toremove.tsv \
-os ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.sample \
-og ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.bgen"
done
done

for i in ${humanomni1[@]}
do echo ${i} && for chr in {1..22}
do echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/${i}_${chr}_subset_postimp_stdout | grep "Successfully"
done
done

for i in ${illuminaexome[@]}
do echo ${i} && for chr in {1..22}
do ls -la ${path_gwas}post_imputation/2022/log/${i}_${chr}_subset_postimp_stdout
done
done



########################################################################################################################################

#######################################
# 2.2 - MORE THAN ONE STUDY PER ARRAY #
#######################################


#######################################
# 2.2.1 exclude list of variants HWE

path_gwas=/path/to/ibdgwas/IIBDGC/
  
illumina370=(belgium_inf1_old_gwas belgium_inf2_old_gwas swedish_uc_old_gwas niddk_old_gwas)
array=illumina370

affymetrix6=(german_affy6_old_gwas norway_affy6_old_gwas gwas2)
array=affymetrix6

illuminaexome=(prism_nfe_gwas helmsley_prism_gsa helmsley_xavier_prism_gsa)
array=illuminaexome

# exclude from list farkkila_gsa; palotie_hus_gsa  - both in finngen, plus ccfa_gsa
gsa=(italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa prism_nfe_gsa 
     lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa mccauley_gsa cedars_gsa bernstein_gsa 
     franchimont_gsa franke_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa 
     newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa pekow_share_gsa 
     rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)
array=gsa


MEM=200

# for i in ${illumina370[@]} - SUBMITTED
# for i in ${affymetrix6[@]} - SUBMITTED
# for i in ${illuminaexome[@]} - SUBMITTED
# for i in ${gsa[@]} - RE-SUBMITTED
do
for chr in {1..22}
do
bsub -J"reg_ibd_test" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -n 8 \
-e ${path_gwas}post_imputation/2022/log/${i}_${chr}_subset_postimp_stderr_vcf \
-o ${path_gwas}post_imputation/2022/log/${i}_${chr}_subset_postimp_stdout_vcf \
"/path/to/software/bcftools-1.16/./bcftools filter \
-e 'ID=@/path/to/ibdgwas/IIBDGC/imputed/'${i}'/2022/eur/list_variants_post_imputation_hwe_toremove.tsv' \
${path_gwas}imputed/${i}/2022/eur/chr${chr}.dose.vcf.gz \
-Oz -o ${path_gwas}post_imputation/2022/${array}/imputed_data/${i}_chr${chr}.dose.subset.vcf.gz \
--threads 8"
done
done



# for i in ${illumina370[@]} - SUBMITTED
# for i in ${affymetrix6[@]} - SUBMITTED
# for i in ${illuminaexome[@]} - SUBMITTED
# for i in ${gsa[@]} - SUBMITTED
do
for chr in X
do
bsub -J"reg_ibd_test" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -n 8 \
-e ${path_gwas}post_imputation/2022/log/${i}_${chr}_subset_postimp_stderr_vcf \
-o ${path_gwas}post_imputation/2022/log/${i}_${chr}_subset_postimp_stdout_vcf \
"/path/to/software/bcftools-1.16/./bcftools filter \
-e 'ID=@/path/to/ibdgwas/IIBDGC/imputed/'${i}'/2022/eur/list_variants_post_imputation_hwe_toremove_chrX.tsv' \
${path_gwas}imputed/${i}/2022/eur/chr${chr}.dose.vcf.gz \
-Oz -o ${path_gwas}post_imputation/2022/${array}/imputed_data/${i}_chr${chr}.dose.subset.vcf.gz \
--threads 8"
done
done



# for i in ${illumina370[@]} - COMPLETED
# for i in ${affymetrix6[@]} - COMPLETED
# for i in ${illuminaexome[@]} - COMPLETED
# for i in ${gsa[@]} - COMPLETED
do echo ${i} && for chr in {1..22} X
do echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/${i}_${chr}_subset_postimp_stdout_vcf | grep "Successfully"
done
done

# for i in ${illumina370[@]} - COMPLETED
# for i in ${affymetrix6[@]} - COMPLETED
# for i in ${illuminaexome[@]} - COMPLETED
# for i in ${gsa[@]} - COMPLETED
do echo ${i} && for chr in {1..22} X
do echo ${chr} && ls -la ${path_gwas}post_imputation/2022/log/${i}_${chr}_subset_postimp_stdout_vcf
done
done


# get N SNPs and samples per vcf:
# for i in ${affymetrix6[@]} 
for i in ${gsa[@]}
do 
echo ${i} && \
for chr in {1..22} X
do 
echo ${chr} &&
/path/to/software/bcftools-1.16/./bcftools  query -l \
${path_gwas}post_imputation/2022/${array}/imputed_data/${i}_chr${chr}.dose.subset.vcf.gz | wc -l
done

# german_affy6_old_gwas
# # 2787
# norway_affy6_old_gwas
# # 542
# gwas2
# # 7758


gsa=( lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa 
     newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa pekow_share_gsa 
     rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

# follow up a specific SNP:
# for i in ${affymetrix6[@]} 
for i in ${gsa[@]}
do
echo ${i} && /path/to/software/bcftools-1.16/./bcftools  query -f \
'%ID %CHROM %POS %AF\n' ${path_gwas}post_imputation/2022/${array}/imputed_data/${i}_chr22.dose.subset.vcf.gz \
| grep 'chr22:22914089:T:G'
done

############################
# german_affy6_old_gwas
# chr22:22914089:T:G chr22 22914089 0.24234
# norway_affy6_old_gwas
# chr22:22914089:T:G chr22 22914089 0.23462
# gwas2
# chr22:22914089:T:G chr22 22914089 0.24541

############################
# italy_gsa
# chr22:22914089:T:G chr22 22914089 0.15929
# kiel_austria_sibdcs_gsa
# chr22:22914089:T:G chr22 22914089 0.22838
# netherlands_gsa
# chr22:22914089:T:G chr22 22914089 0.24537
# slovenia_gsa
# chr22:22914089:T:G chr22 22914089 0.25871
# sweden_gsa
# chr22:22914089:T:G chr22 22914089 0.24043
# niddk_broad_gsa
# chr22:22914089:T:G chr22 22914089 0.20175
# niddk_feinstein_gsa
# chr22:22914089:T:G chr22 22914089 0.18817
# basque_gsa
# chr22:22914089:T:G chr22 22914089 0.21445
# prism_nfe_gsa
# chr22:22914089:T:G chr22 22914089 0.19976
# lithuania_gsa
# chr22:22914089:T:G chr22 22914089 0.2517
# belgium_louis_gsa
# chr22:22914089:T:G chr22 22914089 0.21814
# belgium_franchimont_gsa
# chr22:22914089:T:G chr22 22914089 0.20838
# belgium_vermeire_gsa
# chr22:22914089:T:G chr22 22914089 0.21542
# mccauley_gsa
# chr22:22914089:T:G chr22 22914089 0.17627
# cedars_gsa
# chr22:22914089:T:G chr22 22914089 0.19713
# bernstein_gsa
# chr22:22914089:T:G chr22 22914089 0.2413
# franchimont_gsa
# chr22:22914089:T:G chr22 22914089 0.20355
# franke_gsa
# chr22:22914089:T:G chr22 22914089 0.23825
# hyams_protect_gsa
# chr22:22914089:T:G chr22 22914089 0.2163
# lewis_sparc_gsa
# chr22:22914089:T:G chr22 22914089 0.22326
# mccauley_new_gsa
# chr22:22914089:T:G chr22 22914089 0.18834
# mcgovern_gsa
# chr22:22914089:T:G chr22 22914089 0.20316
# moayyedi_imagine_gsa
# chr22:22914089:T:G chr22 22914089 0.23442
# newberry_share_gsa
# chr22:22914089:T:G chr22 22914089 0.22845
# niddk_cho_gsa
# chr22:22914089:T:G chr22 22914089 0.19037
# niddk_duerr_gsa
# chr22:22914089:T:G chr22 22914089 0.21849
# niddk_rioux_gsa
# chr22:22914089:T:G chr22 22914089 0.19727
# niddk_silverberg_gsa
# chr22:22914089:T:G chr22 22914089 0.19941
# pekow_share_gsa
# chr22:22914089:T:G chr22 22914089 0.21264
# rioux_igenomed_gsa
# chr22:22914089:T:G chr22 22914089 0.2214
# sands_msccr_gsa
# chr22:22914089:T:G chr22 22914089 0.16884
# stampfer_gsa
# chr22:22914089:T:G chr22 22914089 0.22434
# vermeire_gsa
# chr22:22914089:T:G chr22 22914089 0.21451
# weersma_gsa
# chr22:22914089:T:G chr22 22914089 0.26036
# xavier_prism_gsa
# chr22:22914089:T:G chr22 22914089 0.20565
# xavier_share_gsa
# chr22:22914089:T:G chr22 22914089 0.21711
############################


#######################################
# 2.2.2 index vcf files

# for i in ${illumina370[@]} # COMPLETED
# for i in ${affymetrix6[@]} # SUBMITTED
# for i in ${illuminaexome[@]} # COMPLETED
# for i in ${gsa[@]} # SUBMITTED
do 
for chr in {1..22} X
do
bsub -J"reg_ibd_test" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -n 4 \
-e ${path_gwas}post_imputation/2022/log/${i}_${chr}_subset_postimp_index_stderr_vcf \
-o ${path_gwas}post_imputation/2022/log/${i}_${chr}_subset_postimp_index_stdout_vcf \
"/path/to/software/bcftools-1.16/./bcftools index \
${path_gwas}post_imputation/2022/${array}/imputed_data/${i}_chr${chr}.dose.subset.vcf.gz"
done
done


# for i in ${illumina370[@]} - COMPLETED
# for i in ${affymetrix6[@]} - COMPLETED
# for i in ${illuminaexome[@]} - COMPLETED
for i in ${gsa[@]}
do echo ${i} && for chr in {1..22} X
do echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/${i}_${chr}_subset_postimp_index_stdout_vcf | grep "Successfully"
done
done


# for i in ${illumina370[@]} - COMPLETED
# for i in ${affymetrix6[@]} - COMPLETED
# for i in ${illuminaexome[@]} - COMPLETED
# for i in ${gsa[@]} - COMPLETED
do echo ${i} && for chr in {1..22} X
do echo ${chr} &&  ls -la ${path_gwas}post_imputation/2022/log/${i}_${chr}_subset_postimp_index_stdout_vcf
done
done

#######################################
# 2.2.3 merge the files

## ILLUMINA 370

illumina370=(belgium_inf1_old_gwas belgium_inf2_old_gwas swedish_uc_old_gwas niddk_old_gwas)
array=illumina370

for chr in {1..22}
do
bsub -J"reg_ibd_test" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -n 8 \
-e ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_merge_stderr_vcf \
-o ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_merge_stdout_vcf \
"/path/to/software/bcftools-1.16/./bcftools merge \
${path_gwas}post_imputation/2022/${array}/imputed_data/${illumina370[0]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${illumina370[1]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${illumina370[2]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${illumina370[3]}_chr${chr}.dose.subset.vcf.gz \
-Oz -o ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.vcf.gz \
-m none \
--threads 8"
done

# swedish_uc_old_gwas has no chrX data

illumina370=(belgium_inf1_old_gwas belgium_inf2_old_gwas niddk_old_gwas)
array=illumina370

bsub -J"reg_ibd_test" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -n 8 \
-e ${path_gwas}post_imputation/2022/log/${array}_X_subset_merge_stderr_vcf \
-o ${path_gwas}post_imputation/2022/log/${array}_X_subset_merge_stdout_vcf \
"/path/to/software/bcftools-1.16/./bcftools merge \
${path_gwas}post_imputation/2022/${array}/imputed_data/${illumina370[0]}_chrX.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${illumina370[1]}_chrX.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${illumina370[2]}_chrX.dose.subset.vcf.gz \
-Oz -o ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chrX.dose.subset.vcf.gz \
-m none --threads 8"


for chr in {1..22} X
do echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_merge_stdout_vcf | grep "Successfully"
done

for chr in {1..22} X
do ls -la ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_merge_stdout_vcf
done

for chr in {1..22} X
do ls -la ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.vcf.gz
done

# convert into bgen

for chr in {1..22} X
do
bsub -J"reg_ibd_test" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -n 2 \
-e ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_vcf_to_bgen_stderr \
-o ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_vcf_to_bgen_stdout \
"/path/to/software/qctool_v2.2.0/./qctool \
-g ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.vcf.gz \
-filetype vcf -vcf-genotype-field GP \
-os ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.sample \
-og ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.bgen"
done

for chr in {1..22} X
do echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_vcf_to_bgen_stdout | grep "Successfully"
done

for chr in {1..22} X
do ls -la ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_vcf_to_bgen_stdout
done


for chr in {1..22} X
do echo ${chr} && ls -la ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.bgen
done


# remove intermediate files

for i in ${illumina370[@]}
do for chr in {1..22} X
do ls ${path_gwas}post_imputation/2022/${array}/imputed_data/${i}_chr${chr}.dose.subset.vcf.gz
done
done

for chr in {1..22} X
do rm ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.vcf.gz
done

#########################

## ILLUMIANEXOME 
illuminaexome=(prism_nfe_gwas helmsley_prism_gsa helmsley_xavier_prism_gsa)
array=illuminaexome

for chr in {1..22}
do
bsub -J"reg_ibd_test" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -n 8 \
-e ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_merge_stderr_vcf \
-o ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_merge_stdout_vcf \
"/path/to/software/bcftools-1.16/./bcftools merge \
${path_gwas}post_imputation/2022/${array}/imputed_data/${illuminaexome[0]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${illuminaexome[1]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${illuminaexome[2]}_chr${chr}.dose.subset.vcf.gz \
-Oz -o ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.vcf.gz \
-m none --threads 8"
done

for chr in {1..22}
do tail -50 ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_merge_stdout_vcf | grep "Successfully"
done

for chr in {1..22}
do ls -la ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_merge_stdout_vcf
done

for chr in {1..22}
do ls -la ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.vcf.gz
done

bsub -J"reg_ibd_test" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -n 8 \
-e ${path_gwas}post_imputation/2022/log/${array}_X_subset_merge_stderr_vcf \
-o ${path_gwas}post_imputation/2022/log/${array}_X_subset_merge_stdout_vcf \
"/path/to/software/bcftools-1.16/./bcftools merge \
${path_gwas}post_imputation/2022/${array}/imputed_data/${illuminaexome[0]}_chrX.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${illuminaexome[1]}_chrX.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${illuminaexome[2]}_chrX.dose.subset.vcf.gz \
-Oz -o ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chrX.dose.subset.vcf.gz \
-m none --threads 8"

tail -50 ${path_gwas}post_imputation/2022/log/${array}_X_subset_merge_stdout_vcf | grep "Successfully"

ls -la ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chrX.dose.subset.vcf.gz

# convert into bgen

for chr in {1..22} X
do
bsub -J"reg_ibd_test" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -n 2 \
-e ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_vcf_to_bgen_stderr \
-o ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_vcf_to_bgen_stdout \
"/path/to/software/qctool_v2.2.0/./qctool \
-g ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.vcf.gz \
-filetype vcf -vcf-genotype-field GP \
-os ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.sample \
-og ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.bgen"
done

for chr in {1..22} X
do ls -la ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_vcf_to_bgen_stdout
done

for chr in {1..22} X
do tail -50 ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_vcf_to_bgen_stdout | grep "Successfully"
done

for chr in {1..22} X
do ls -la ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.bgen
done

# remove intermediate files

for i in ${illuminaexome[@]}
do for chr in {1..22} X
do rm ${path_gwas}post_imputation/2022/${array}/imputed_data/${i}_chr${chr}.dose.subset.vcf.gz
done
done

for chr in {1..22} X
do rm ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.vcf.gz
done
done


#########################

## AFFYMETRIX6 

affymetrix6=(german_affy6_old_gwas norway_affy6_old_gwas gwas2)
array=affymetrix6

for chr in {1..22} X
do
bsub -J"reg_ibd_test" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -n 8 -q long \
-e ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_merge_stderr_vcf \
-o ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_merge_stdout_vcf \
"/path/to/software/bcftools-1.16/./bcftools merge \
${path_gwas}post_imputation/2022/${array}/imputed_data/${affymetrix6[0]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${affymetrix6[1]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${affymetrix6[2]}_chr${chr}.dose.subset.vcf.gz \
-Oz -o ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.vcf.gz \
-m none --threads 8"
done

for chr in {1..22} X
do echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_merge_stdout_vcf | grep "Successfully"
done

for chr in {1..22} X
do ls -la ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.vcf.gz
done

# follow up a specific SNP:
echo ${i} && /path/to/software/bcftools-1.16/./bcftools  query -f \
'%ID %CHROM %POS %AF\n' ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr22.dose.subset.vcf.gz \
| grep 'chr22:22914089:T:G'
# chr22:22914089:T:G chr22 22914089 0.24234


# convert into bgen

for chr in {1..5} 11
do
bsub -J"reg_ibd_test" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -m "modern_hardware" -G ibdgwas -n 2 -q long \
-e ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_vcf_to_bgen_stderr \
-o ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_vcf_to_bgen_stdout \
"/path/to/software/qctool_v2.2.0/./qctool \
-g ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.vcf.gz \
-filetype vcf -vcf-genotype-field GP \
-os ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.sample \
-og ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.bgen"
done

# SUBMITTED


for chr in {1..22} X
do echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_vcf_to_bgen_stdout | grep "Successfully"
done

for chr in {1..22} X
do ls -la ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_vcf_to_bgen_stdout
done

for chr in {1..22} X
do echo ${chr} && ls -la ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.bgen
done


# remove intermediate files

for i in ${affymetrix6[@]}
do for chr in {1..22} X
do rm ${path_gwas}post_imputation/2022/${array}/imputed_data/${i}_chr${chr}.dose.subset.vcf.gz
done
done

for chr in {1..22} X
do rm ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.vcf.gz
done




#########################

## GSA 

# exclude from list farkkila_gsa; palotie_hus_gsa  - both in finngen, plus ccfa_gsa
gsa=(italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa prism_nfe_gsa 
     lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa mccauley_gsa cedars_gsa bernstein_gsa 
     franchimont_gsa franke_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa 
     newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa pekow_share_gsa 
     rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)
array=gsa

for i in ${gsa[@]} 
do echo ${i} && for chr in 22
do ls -la ${path_gwas}post_imputation/2022/${array}/imputed_data/${i}_chr${chr}.dose.subset.vcf.gz
done
done

MEM=300

for chr in {1..22} X
do
bsub -J"reg_ibd_test" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -n 8 -q basement \
-e ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_merge_stderr_vcf \
-o ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_merge_stdout_vcf \
"/path/to/software/bcftools-1.16/./bcftools merge \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[0]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[1]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[2]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[3]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[4]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[5]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[6]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[7]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[8]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[9]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[10]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[11]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[12]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[13]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[14]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[15]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[16]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[17]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[18]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[19]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[20]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[21]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[22]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[23]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[24]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[25]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[26]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[27]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[28]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[29]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[30]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[31]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[32]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[33]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[34]}_chr${chr}.dose.subset.vcf.gz \
${path_gwas}post_imputation/2022/${array}/imputed_data/${gsa[35]}_chr${chr}.dose.subset.vcf.gz \
-Oz -o ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.vcf.gz \
-m none --threads 8"
done

# SUBMITTED


for chr in {1..22} X
do echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_merge_stdout_vcf | grep "Successfully"
done

for chr in {1..22} X
do ls -la ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.vcf.gz
done


# convert into bgen - submitted 30 march

for chr in 2
do
bsub -J"reg_ibd_test" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -n 2 -q basement \
-e ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_vcf_to_bgen_stderr \
-o ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_vcf_to_bgen_stdout \
"/path/to/software/qctool_v2.2.0/./qctool \
-g ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.vcf.gz \
-filetype vcf -vcf-genotype-field GP \
-os ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.sample \
-og ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.bgen"
done


for chr in {1..22} X
do echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/${array}_${chr}_subset_vcf_to_bgen_stdout | grep "Successfully"
done



path_gwas=/path/to/ibdgwas/IIBDGC/
array=(gsa affymetrix6)
chr=22
MEM=200

for i in ${array[@]}
do
bsub -J"snp_stats" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -n 2 -q long \
-e ${path_gwas}post_imputation/2022/log/${i}_${chr}_snp_stats_merged_file_stderr \
-o ${path_gwas}post_imputation/2022/log/${i}_${chr}_snp_stats_merged_file_stdout \
"/path/to/software/qctool_v2.2.0/./qctool \
-g ${path_gwas}post_imputation/2022/${i}/imputed_data/${i}_chr${chr}.dose.subset.vcf.gz \
-snp-stats -osnp ${path_gwas}post_imputation/2022/${i}/imputed_data/${i}_chr${chr}_snp_stats.txt"
done




# remove intermediate files

for i in ${gsa[@]}
do for chr in {1..22} X
do rm ${path_gwas}post_imputation/2022/${array}/imputed_data/${i}_chr${chr}.dose.subset.vcf.gz
done
done


for i in ${gsa[@]}
do for chr in {1..22} X
do rm ${path_gwas}post_imputation/2022/${array}/imputed_data/${i}_chr${chr}.dose.subset.vcf.gz.csi
done
done

for chr in {1..22} X
do rm ${path_gwas}post_imputation/2022/${array}/imputed_data/${array}_chr${chr}.dose.subset.vcf.gz
done



# look into number of merged samples - to see if it is as expected:

arrays=(illumina370 affymetrix6 gsa illuminaexome)

for i in ${arrays[@]}
do echo ${i} && for chr in {1..22} X
do echo ${chr} && tail ${path_gwas}post_imputation/2022/log/${i}_${chr}_subset_vcf_to_bgen_stderr | grep -E "Number of samples|Number of samples"
done
done


##########################
# illumina370
# 1
# Number of samples in input file(s):   5633.
# 2
# Number of samples in input file(s):   5633.
# 3
# Number of samples in input file(s):   5633.
# 4
# Number of samples in input file(s):   5633.
# 5
# Number of samples in input file(s):   5633.
# 6
# Number of samples in input file(s):   5633.
# 7
# Number of samples in input file(s):   5633.
# 8
# Number of samples in input file(s):   5633.
# 9
# Number of samples in input file(s):   5633.
# 10
# Number of samples in input file(s):   5633.
# 11
# Number of samples in input file(s):   5633.
# 12
# Number of samples in input file(s):   5633.
# 13
# Number of samples in input file(s):   5633.
# 14
# Number of samples in input file(s):   5633.
# 15
# Number of samples in input file(s):   5633.
# 16
# Number of samples in input file(s):   5633.
# 17
# Number of samples in input file(s):   5633.
# 18
# Number of samples in input file(s):   5633.
# 19
# Number of samples in input file(s):   5633.
# 20
# Number of samples in input file(s):   5633.
# 21
# Number of samples in input file(s):   5633.
# 22
# Number of samples in input file(s):   5633.
# X
# Number of samples in input file(s):   4378.
# affymetrix6
# 1
# Number of samples in input file(s):   11087.
# 2
# Number of samples in input file(s):   11087.
# 3
# Number of samples in input file(s):   11087.
# 4
# Number of samples in input file(s):   11087.
# 5
# Number of samples in input file(s):   11087.
# 6
# Number of samples in input file(s):   11087.
# 7
# Number of samples in input file(s):   11087.
# 8
# Number of samples in input file(s):   11087.
# 9
# Number of samples in input file(s):   11087.
# 10
# Number of samples in input file(s):   11087.
# 11
# Number of samples in input file(s):   11087.
# 12
# Number of samples in input file(s):   11087.
# 13
# Number of samples in input file(s):   11087.
# 14
# Number of samples in input file(s):   11087.
# 15
# Number of samples in input file(s):   11087.
# 16
# Number of samples in input file(s):   11087.
# 17
# Number of samples in input file(s):   11087.
# 18
# Number of samples in input file(s):   11087.
# 19
# Number of samples in input file(s):   11087.
# 20
# Number of samples in input file(s):   11087.
# 21
# Number of samples in input file(s):   11087.
# 22
# Number of samples in input file(s):   11087.
# X
# Number of samples in input file(s):   11087.
# gsa
# 1
# Number of samples in input file(s):   78481.
# 2
# Number of samples in input file(s):   78481.
# 3
# Number of samples in input file(s):   78481.
# 4
# Number of samples in input file(s):   78481.
# 5
# Number of samples in input file(s):   78481.
# 6
# Number of samples in input file(s):   78481.
# 7
# Number of samples in input file(s):   78481.
# 8
# Number of samples in input file(s):   78481.
# 9
# Number of samples in input file(s):   78481.
# 10
# Number of samples in input file(s):   78481.
# 11
# Number of samples in input file(s):   78481.
# 12
# Number of samples in input file(s):   78481.
# 13
# Number of samples in input file(s):   78481.
# 14
# Number of samples in input file(s):   78481.
# 15
# Number of samples in input file(s):   78481.
# 16
# Number of samples in input file(s):   78481.
# 17
# Number of samples in input file(s):   78481.
# 18
# Number of samples in input file(s):   78481.
# 19
# Number of samples in input file(s):   78481.
# 20
# Number of samples in input file(s):   78481.
# 21
# Number of samples in input file(s):   78481.
# 22
# Number of samples in input file(s):   78481.
# X
# Number of samples in input file(s):   78480.
# illuminaexome
# 1
# Number of samples in input file(s):   2557.
# 2
# Number of samples in input file(s):   2557.
# 3
# Number of samples in input file(s):   2557.
# 4
# Number of samples in input file(s):   2557.
# 5
# Number of samples in input file(s):   2557.
# 6
# Number of samples in input file(s):   2557.
# 7
# Number of samples in input file(s):   2557.
# 8
# Number of samples in input file(s):   2557.
# 9
# Number of samples in input file(s):   2557.
# 10
# Number of samples in input file(s):   2557.
# 11
# Number of samples in input file(s):   2557.
# 12
# Number of samples in input file(s):   2557.
# 13
# Number of samples in input file(s):   2557.
# 14
# Number of samples in input file(s):   2557.
# 15
# Number of samples in input file(s):   2557.
# 16
# Number of samples in input file(s):   2557.
# 17
# Number of samples in input file(s):   2557.
# 18
# Number of samples in input file(s):   2557.
# 19
# Number of samples in input file(s):   2557.
# 20
# Number of samples in input file(s):   2557.
# 21
# Number of samples in input file(s):   2557.
# 22
# Number of samples in input file(s):   2557.
# X
# Number of samples in input file(s):   2557.








