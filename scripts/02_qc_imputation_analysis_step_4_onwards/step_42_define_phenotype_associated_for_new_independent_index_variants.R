# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# # # FOR EACH ONE OF THE ASSOCIATED variants, DEFINE PHENOTYPE


# # # 1.- EXTRACT VARIANTS FROM LIST OF INDEX VARIANTS DEFINED IN STEP 38:

# # # singularity exec iibdgc_postprocess_10_singularity.sif

# # path_gwas="/path/to/ibdgwas/IIBDGC/"
# # ph=ibd
# # array=(humancoreexome gsa)

# # # list of credible sets by Ruize:
# # /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/fine_mapping/tier1_EUR_array_data_fine_mapping_results_06_01_2025_list_of_variants.csv
# # /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/fine_mapping/tier1_EUR_array_data_fine_mapping_results_06_01_2025.xlsx

# # cut -d',' -f10  ${path_gwas}post_imputation/2022/analysis/metaanalysis/fine_mapping/tier1_EUR_array_data_fine_mapping_results_06_01_2025_list_of_variants.csv > \
# # ${path_gwas}post_imputation/2022/analysis/metaanalysis/fine_mapping/tier1_EUR_array_data_fine_mapping_results_06_01_2025_MarkerNames

# # wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/fine_mapping/tier1_EUR_array_data_fine_mapping_results_06_01_2025_MarkerNames
# # # 7222

# # # list of index variants from COJO:
# # zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal.tsv.gz | cut -f1 > \
# # ${path_gwas}post_imputation/2022/analysis/metaanalysis/fine_mapping/list_tier_2_unsupervised_conditional_frequency_with_old_index_MarkerNames
# # wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/fine_mapping/list_tier_2_unsupervised_conditional_frequency_with_old_index_MarkerNames
# # # 1375

# # # list of old variants:
# # cut -f3 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/forward_regression/list_known_index_gwas_fm_2023_with_iibdgc_tier1_eas_summary_stats.tsv > \
# # ${path_gwas}post_imputation/2022/analysis/metaanalysis/fine_mapping/list_old_index_variants_MarkerNames
# # wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/fine_mapping/list_old_index_variants_MarkerNames
# # # 391

# # cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/fine_mapping/tier1_EUR_array_data_fine_mapping_results_06_01_2025_MarkerNames \
# # ${path_gwas}post_imputation/2022/analysis/metaanalysis/fine_mapping/list_tier_2_unsupervised_conditional_frequency_with_old_index_MarkerNames \
# # ${path_gwas}post_imputation/2022/analysis/metaanalysis/fine_mapping/list_old_index_variants_MarkerNames | sort | uniq > \
# # ${path_gwas}post_imputation/2022/analysis/metaanalysis/fine_mapping/tier1_EUR_array_data_fine_mapping_tier_2_unsupervised_conditional_old_index_markerNames
# # wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/fine_mapping/tier1_EUR_array_data_fine_mapping_tier_2_unsupervised_conditional_old_index_markerNames
# # # 8123 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/fine_mapping/tier1_EUR_array_data_fine_mapping_tier_2_unsupervised_conditional_old_index_markerNames

# # paste <(cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/fine_mapping/tier1_EUR_array_data_fine_mapping_tier_2_unsupervised_conditional_old_index_markerNames) <(cut -d':' -f4 ${path_gwas}post_imputation/2022/analysis/metaanalysis/fine_mapping/tier1_EUR_array_data_fine_mapping_tier_2_unsupervised_conditional_old_index_markerNames| cat) > \
# # ${path_gwas}post_imputation/2022/analysis/metaanalysis/fine_mapping/tier1_EUR_array_data_fine_mapping_tier_2_unsupervised_conditional_old_index_markerNames_AltAllele
# # wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/fine_mapping/tier1_EUR_array_data_fine_mapping_tier_2_unsupervised_conditional_old_index_markerNames_AltAllele
# # # 8123 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/fine_mapping/tier1_EUR_array_data_fine_mapping_tier_2_unsupervised_conditional_old_index_markerNames_AltAllele


# # # retrieve variant frequency:

# # MEM=500

# # for i in ${array[@]}
# # do 
# # for chr in {1..22} X
# # do
# # bsub -J"bgen_bed" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
# # -o ${path_gwas}post_imputation/2022/log/${i}_${chr}_freq_new_and_old_index_gwas_variants_plus_credible_sets_stdout \
# # -e ${path_gwas}post_imputation/2022/log/${i}_${chr}_freq_new_and_old_index_gwas_variants_plus_credible_sets_stderr \
# # "plink2 \
# # --bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${i}_chr${chr}_subset_included_in_${ph}_analysis \
# # --threads 4 \
# # --extract ${path_gwas}post_imputation/2022/analysis/metaanalysis/fine_mapping/tier1_EUR_array_data_fine_mapping_tier_2_unsupervised_conditional_old_index_markerNames_AltAllele \
# # --freq 'counts' --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${i}_chr${chr}_subset_eur_tier_1_fm_tier_2_unsupervised_conditional_old"
# # done
# # done



# # for i in ${array[@]}
# # do 
# # echo ${i} && for chr in {1..22} X
# # do
# # echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/${i}_${chr}_freq_new_and_old_index_gwas_variants_plus_credible_sets_stdout | grep "Successfully"
# # done
# # done


# # ### Check out frequency (counts):

# # MEM=2000
# # bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q normal R 


# # library(data.table)


# # path_gwas<-"/path/to/ibdgwas/IIBDGC/"

# # # keep only samples two largest arrays, including small ones and adjust by array does not allow to converge the model:
# # array<-c("humancoreexome","gsa")

# # rm(dat)
# # for (chr in c(1:22)) {

# #   for (i in 1:length(array)) {

# #     file_raw<-paste(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/",array[i],"_chr",chr,"_subset_eur_tier_1_fm_tier_2_unsupervised_conditional_old.acount",sep="")
    
# #     if(file.exists(file_raw)) {

# #       tmp<-read.table(file_raw,head=F,check.names=F)
# #       tmp$min_freq<-pmin(tmp$V5,tmp$V6)/rowSums(tmp[,c(5:6)])
# #       tmp<-tmp[,c(2,7)]
# #       colnames(tmp)[2]<-paste(colnames(tmp)[2],array[i],sep="_")
      
# #       if(i==1){
# #         tmp1<-tmp
# #       } else {
# #         tmp1<-merge(tmp1,tmp,by="V2",all=T)
# #       }
# #     }
# #     rm(tmp)
# #   }

# #   if(file.exists(file_raw)) {
# #     if(!exists("dat")){
# #       dat<-tmp1
# #     }else{
# #       dat<-rbind(dat,tmp1)
# #     }
# #     rm(tmp1)
# #   } else {
# #     next
# #   }
  
# # }


# # min(dat$min_freq_humancoreexome,na.rm=T)
# # # [1] 0.00026826
# # min(dat$min_freq_gsa,na.rm=T)
# # # [1] 0.000778788

# # variants<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/fine_mapping/tier1_EUR_array_data_fine_mapping_tier_2_unsupervised_conditional_old_index_markerNames_AltAllele"),head=F)
# # dim(variants)
# # # [1] 8122     1

# # variants<-variants[!duplicated(variants),]
# # dim(variants)
# # # [1] 8122    1

# # dim(dat)
# # # [1] 8108    3

# # variants$V1[which(!variants$V1 %in% dat$V2)]
# # #  [1] "chr10:101436625:AATAGATAGATAGATAG:AATAGATAGATAGATAGATAG"
# # #  [2] "chr10:98808605:C:T"                                     
# # #  [3] "chr11:118273990:A:G"                                    
# # #  [4] "chr1:151281565:T:A"                                     
# # #  [5] "chr11:73152652:A:G"                                     
# # #  [6] "chr12:110259525:G:A"                                    
# # #  [7] "chr1:28131939:C:T"                                      
# # #  [8] "chr14:105963206:C:A"                                    
# # #  [9] "chr14:105968520:G:A"                                    
# # # [10] "chr21:44193942:T:C"                                     
# # # [11] "chr6:32247045:T:G"                                      
# # # [12] "chr6:34550494-G-A"                                      
# # # [13] "chr7:5436753:C:T"   

# # # keep only present in boht

# # dat<-dat[which(!is.na(dat$min_freq_humancoreexome) & !is.na(dat$min_freq_gsa)),]
# # dim(dat)
# # # [1] 8070    3
# # variants$V1[which(!variants$V1 %in% dat$V2)]
# #  [1] "chr10:101436625:AATAGATAGATAGATAG:AATAGATAGATAGATAGATAG"
# #  [2] "chr10:98808605:C:T"                                     
# #  [3] "chr11:118273990:A:G"                                    
# #  [4] "chr1:151281565:T:A"                                     
# #  [5] "chr11:73152652:A:G"                                     
# #  [6] "chr12:110259525:G:A"                                    
# #  [7] "chr1:28131939:C:T"                                      
# #  [8] "chr14:105963206:C:A"                                    
# #  [9] "chr14:105968520:G:A"                                    
# # [10] "chr21:44193942:T:C"                                     
# # [11] "chr6:32247045:T:G"                                      
# # [12] "chr6:34550494-G-A"                                      
# # [13] "chr7:5436753:C:T"                                       
# # [14] "MarkerName"                                             
# # > dat<-dat[which(!is.na(dat$min_freq_humancoreexome) & !is.na(dat$min_freq_gsa)),]
# # > dim(dat)
# # [1] 8070    3
# # > variants$V1[which(!variants$V1 %in% dat$V2)]
# # #  [1] "chr10:101436625:AATAGATAGATAGATAG:AATAGATAGATAGATAGATAG"
# # #  [2] "chr10:105644358:TA:T"                                   
# # #  [3] "chr10:125821742:T:G"                                    
# # #  [4] "chr10:60750152:A:G"                                     
# # #  [5] "chr10:98808605:C:T"                                     
# # #  [6] "chr11:118273990:A:G"                                    
# # #  [7] "chr1:151281565:T:A"                                     
# # #  [8] "chr11:73152652:A:G"                                     
# # #  [9] "chr1:19814567:G:A"                                      
# # # [10] "chr12:110259525:G:A"                                    
# # # [11] "chr1:28131939:C:T"                                      
# # # [12] "chr14:105686797:A:G"                                    
# # # [13] "chr14:105702719:A:G"                                    
# # # [14] "chr14:105721972:C:T"                                    
# # # [15] "chr14:105722906:G:T"                                    
# # # [16] "chr14:105731743:G:A"                                    
# # # [17] "chr14:105737956:C:T"                                    
# # # [18] "chr14:105742782:T:C"                                    
# # # [19] "chr14:105749143:G:A"                                    
# # # [20] "chr14:105751077:G:T"                                    
# # # [21] "chr14:105753554:C:T"                                    
# # # [22] "chr14:105757817:A:G"                                    
# # # [23] "chr14:105760146:T:A"                                    
# # # [24] "chr14:105761095:T:C"                                    
# # # [25] "chr14:105761758:A:G"                                    
# # # [26] "chr14:105769791:T:A"                                    
# # # [27] "chr14:105772740:G:A"                                    
# # # [28] "chr14:105777786:C:A"                                    
# # # [29] "chr14:105778578:A:G"                                    
# # # [30] "chr14:105779242:C:A"                                    
# # # [31] "chr14:105779492:C:T"                                    
# # # [32] "chr14:105781212:A:G"                                    
# # # [33] "chr14:105781644:C:T"                                    
# # # [34] "chr14:105782062:G:A"                                    
# # # [35] "chr14:105783193:A:C"                                    
# # # [36] "chr14:105783627:T:C"                                    
# # # [37] "chr14:105785012:G:T"                                    
# # # [38] "chr14:105795065:A:T"                                    
# # # [39] "chr14:105963206:C:A"                                    
# # # [40] "chr14:105968520:G:A"                                    
# # # [41] "chr15:34731239:A:C"                                     
# # # [42] "chr15:44802938:T:C"                                     
# # # [43] "chr16:27384341:C:CT"                                    
# # # [44] "chr16:50716899:A:G"                                     
# # # [45] "chr21:44193942:T:C"                                     
# # # [46] "chr4:42396369:G:T"                                      
# # # [47] "chr4:42635348:ATATC:A"                                  
# # # [48] "chr5:4422049:C:T"                                       
# # # [49] "chr6:32247045:T:G"                                      
# # # [50] "chr6:34550494-G-A"                                      
# # # [51] "chr7:5436753:C:T"


# # # exclude as well those with very low freq:
# # dat[which(dat$min_freq_humancoreexome<0.001 | dat$min_freq_gsa<0.001),"V2"]
# # #  [1] "chr1:66914276:C:G"   "chr2:18792705:C:CT"  "chr7:4827597:G:C"   
# # #  [4] "chr7:4881854:T:C"    "chr7:4981470:C:G"    "chr7:5139807:G:A"   
# # #  [7] "chr11:117998788:C:T" "chr16:50476064:A:G"  "chr16:50500888:G:T" 
# # # [10] "chr16:50537867:C:T"  "chr16:50588216:A:G"  "chr16:50593427:A:G" 
# # # [13] "chr16:50699481:T:A"  "chr16:50787956:A:G"

# # dat[which(dat$V2=="chr16:50537867:C:T"),]
# # #                      V2 min_freq_humancoreexome min_freq_gsa
# # # 6578 chr16:50537867:C:T             0.000656535  0.004124066


# # dat<-dat[which(dat$min_freq_humancoreexome>=0.001 & dat$min_freq_gsa>=0.001),]
# # dim(dat)
# # # [1] 8056    3

# # variants<-variants[which(variants$V1 %in% dat$V2),]
# # dim(variants)
# # # [1] 8056    2

# # dat$alt<-gsub("chr[0-9]{1,2}:[0-9]*:[A-Z]*:","",dat$V2)

# # fwrite(dat[,c("V2","alt")],paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/tier1_EUR_array_data_fine_mapping_tier_2_unsupervised_conditional_old_index_markerNames_AltAllele"),
# # col.names=T,row.names=F,quote=F,sep="\t")

# # fwrite(variants[,1,drop=F],paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/tier1_EUR_array_data_fine_mapping_tier_2_unsupervised_conditional_old_index_markerNames"),
# # col.names=T,row.names=F,quote=F,sep="\t")

# # q("no")

# # ###########################################


# # MEM=900

# # for i in ${array[@]}
# # do 
# # for chr in {1..22}
# # do
# # bsub -J"bgen_bed" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G humgen-priority -n 4 \
# # -o ${path_gwas}post_imputation/2022/log/${i}_${chr}_bed_to_raw_stdout \
# # -e ${path_gwas}post_imputation/2022/log/${i}_${chr}_bed_to_raw_stderr \
# # "plink2 \
# # --bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${i}_chr${chr}_subset_included_in_${ph}_analysis \
# # --threads 4 \
# # --extract ${path_gwas}post_imputation/2022/analysis/metaanalysis/forward_regression/tier1_EUR_array_data_fine_mapping_tier_2_unsupervised_conditional_old_index_markerNames \
# # --export-allele ${path_gwas}post_imputation/2022/analysis/metaanalysis/forward_regression/tier1_EUR_array_data_fine_mapping_tier_2_unsupervised_conditional_old_index_markerNames_AltAllele \
# # --export A --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${i}_chr${chr}_tier1_EUR_array_data_fine_mapping_tier_2_unsupervised_conditional_old_index"
# # done
# # done

# # # completed


# # for i in ${array[@]}
# # do 
# # echo ${i} && for chr in {1..22} X
# # do
# # echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/${i}_${chr}_bed_to_raw_stdout | grep "Successfully"
# # done
# # done

# # for i in ${array[@]}
# # do 
# # for chr in {1..22}
# # do
# # ls -la ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${i}_chr${chr}_tier1_EUR_array_data_fine_mapping_tier_2_unsupervised_conditional_old_index.raw
# # done
# # done

# # for i in ${array[@]}
# # do 
# # for chr in {1..22}
# # do
# # echo ${chr} && tail -100 ${path_gwas}post_imputation/2022/log/${i}_${chr}_bed_to_raw_stdout | grep "variants remaining"
# # done
# # done


# # # Define the phenotype using a multinomial model:


# # MEM=15000

# # for chr in {1..22}
# # do
# # bsub -J"ph" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 16 -q long \
# # -e ${path_gwas}post_imputation/2022/log/def_associated_ph_index_credible_set_${chr}_stderr \
# # -o ${path_gwas}post_imputation/2022/log/def_associated_ph_index_credible_set_${chr}_stdout \
# # "Rscript ~/git/IIBDGC_GWAS/scripts/other/define_associated_phenotype_index_variants_mc16.R ${chr} > \
# # ${path_gwas}post_imputation/2022/log/define_associated_phenotype_index_variants_mc16_${chr}.Rout"
# # done


# # for chr in {1..22}
# # do
# # echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/def_associated_ph_index_credible_set_${chr}_stdout | grep -E "Successfully|Exited"
# # done



# # ########################################################################
# # # combine all summary stats to share with Rui

# # # MEM=5000
# # # bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R


# # library(data.table)

# # path_gwas<-"/path/to/ibdgwas/IIBDGC/"

# # files<-list.files(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/"))

# # files<-files[grep("_tier1_EUR_array_data_fine_mapping_tier_2_unsupervised_conditional_old_index_subphenotype_models_adjusted_by_array_mc16",files)]

# # rm(dat)

# # for (i in 1:length(files)) {

# #   tmp<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/",files[i]),head=T)
# #   tmp$MarkerName<-gsub("_[A-Z]*$","",tmp$var_id)

# #   if(i==1) {
# #     dat<-tmp
# #   }else{
# #     dat<-rbind(dat,tmp)
# #   }
# #   rm(tmp)

# # }

# # dim(dat)
# # # [1] 8055    11

# # dat<-dat[!duplicated(dat),]
# # dat[which(duplicated(dat$MarkerName)),]
# # # 0

# # dim(dat)
# # # [1] 8055    11

# # fwrite(dat,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/allchr_tier1_EUR_array_data_fine_mapping_tier_2_unsupervised_conditional_old_index_subphenotype_models_adjusted_by_array_mc16.tsv.gz"),
# # col.names=T,row.names=F,quote=F,sep="\t")

# # q("no")



# ###########################################
# # === ACTIVE CODE (steps below are not commented out) ===
# ###########################################

# # singularity exec iibdgc_postprocess_10_singularity.sif


# # RE RUN THE MODELS WITH PCS JUST FOR THE FINAL LIST OF VARIANTS

# MEM=3000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q yesterday R

# library(data.table)
# rm(list=ls())

# path_gwas<-"/path/to/ibdgwas/IIBDGC/"


# all<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model.tsv.gz",sep=""))
# all<-as.data.frame(all)
# dim(all)
# # [1] 689 163

# variants<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/fine_mapping/tier1_EUR_array_data_fine_mapping_tier_2_unsupervised_conditional_old_index_markerNames_AltAllele"),head=F)
# dim(variants)

# variants<-variants[which(variants$V1 %in% all$MarkerName),]
# dim(variants)
# # [1] 689   2

# fwrite(variants,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_final_join_model_markerNames_AltAllele"),
# col.names=T,row.names=F,quote=F,sep="\t")

# fwrite(variants[,1,drop=F],paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_final_join_model_markerNames"),
# col.names=T,row.names=F,quote=F,sep="\t")

# q("no")

# ##########################################

# MEM=900


# path_gwas="/path/to/ibdgwas/IIBDGC/"
# ph=ibd
# array=(humancoreexome gsa)

# for i in ${array[@]}
# do 
# for chr in {1..22}
# do
# bsub -J"bgen_bed" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -n 4 \
# -o ${path_gwas}post_imputation/2022/log/${i}_${chr}_bed_to_raw_stdout \
# -e ${path_gwas}post_imputation/2022/log/${i}_${chr}_bed_to_raw_stderr \
# "plink2 \
# --bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${i}_chr${chr}_subset_included_in_${ph}_analysis \
# --threads 4 \
# --extract ${path_gwas}post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_final_join_model_markerNames \
# --export-allele ${path_gwas}post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_final_join_model_markerNames_AltAllele \
# --export A --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${i}_chr${chr}_list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_with_fm_signals_with_multiancestry_signals"
# done
# done




# for chr in {1..22}
# do 
# echo chr${chr} && for i in ${array[@]}
# do 
# head -1 ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${i}_chr${chr}_list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_with_fm_signals_with_multiancestry_signals.raw | wc -w
# done
# done



# for chr in {1..22}
# do 
# echo chr${chr} && for i in ${array[@]}
# do 
# ls -la ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${i}_chr${chr}_list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_with_fm_signals_with_multiancestry_signals.raw 
# done
# done



# # chr1
# # 90
# # 91
# # chr2
# # 66
# # 66
# # chr3
# # 34
# # 34
# # chr4
# # 30
# # 32
# # chr5
# # 47
# # 47
# # chr6
# # 40
# # 40
# # chr7
# # 63
# # 63
# # chr8
# # 31
# # 31
# # chr9
# # 35
# # 35
# # chr10
# # 37
# # 40
# # chr11
# # 40
# # 40
# # chr12
# # 29
# # 29
# # chr13
# # 20
# # 20
# # chr14
# # 23
# # 24
# # chr15
# # 23
# # 24
# # chr16
# # 46
# # 47
# # chr17
# # 38
# # 38
# # chr18
# # 16
# # 16
# # chr19
# # 34
# # 34
# # chr20
# # 30
# # 30
# # chr21
# # 17
# # 17
# # chr22
# # 21
# # 21

# MEM=15000
# for chr in {1..22}
# do
# bsub -J"ph" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 16 -q normal \
# -e ${path_gwas}post_imputation/2022/log/def_associated_ph_index_credible_set_with_PCs_${chr}_stderr \
# -o ${path_gwas}post_imputation/2022/log/def_associated_ph_index_credible_set_with_PCs_${chr}_stdout \
# "Rscript ~/git/IIBDGC_GWAS/scripts/other/define_associated_phenotype_index_variants_mc16_with_PCs.R ${chr} > \
# ${path_gwas}post_imputation/2022/log/define_associated_phenotype_index_variants_mc16_with_PCs_${chr}.Rout"
# done

# for chr in {1..22}
# do
# echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/def_associated_ph_index_credible_set_with_PCs_${chr}_stdout | grep -E "Successfully|Exited"
# done


# for chr in {1..22}
# do
# ls -la ${path_gwas}post_imputation/2022/log/def_associated_ph_index_credible_set_with_PCs_${chr}_stdout
# done

# ######################


# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=2000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q yesterday R

library(data.table)

rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

gwas_dat<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_no_subset_final_2.tsv.gz",sep=""))
gwas_dat<-as.data.frame(gwas_dat)
dim(gwas_dat)
# [1] 619 173


# ADD THE ASSOCIATED PHENOTYPE
ph_files <- paste0(path_gwas,
  "post_imputation/2022/analysis/metaanalysis/chr", 1:22,
  "_tier1_EUR_array_data_fine_mapping_tier_2_unsupervised_conditional_old_index_subphenotype_models_adjusted_by_array_by_pcs_mc16_with_PCs.tsv.gz")
ph <- rbindlist(lapply(ph_files, fread))

ph$MarkerName<-gsub("_[A-Z]*$","",ph$var_id)
dim(ph)
# [1] 686  11

dim(gwas_dat[which(gwas_dat$MarkerName %in% ph$MarkerName),])
# [1] 618 211

gwas_dat$MarkerName[which(!gwas_dat$MarkerName %in% ph$MarkerName)]
# [1] "chr1:66934185:C:T"

gwas_dat<-merge(gwas_dat,ph,by="MarkerName",all.x=T)

table(gwas_dat$phenotype,useNA="ifany")
            #  CD   IBD_saturated IBD_unsaturated              UC            <NA>
            # 132             114             277              95               1


# fill in variants where the model could not converge, or those not significant in EUR tier 1

list_ids_1<-gwas_dat$MarkerName[which(is.na(gwas_dat$phenotype))]
list_ids_2<-gwas_dat$MarkerName[which((gwas_dat$'P-value_ibd_eur_tier_1'>1E-5) & (gwas_dat$'P-value_cd_eur_tier_1'>1E-5) & (gwas_dat$'P-value_uc_eur_tier_1'>1E-5))]
length(list_ids_2)
# [1] 296

list_ids<-unique(c(list_ids_1,list_ids_2))

p_cols_eur   <- c("P-value_uc_eur_tier_2",        "P-value_cd_eur_tier_2",        "P-value_ibd_eur_tier_2")
p_cols_multi <- c("P-value_uc_eur_tier_2_eas_sas", "P-value_cd_eur_tier_2_eas_sas", "P-value_ibd_eur_tier_2_eas_sas")
ph_labels    <- c("UC", "CD", "IBD_unsaturated")

assign_phenotype_fallback <- function(row) {
  p_vals <- suppressWarnings(as.numeric(row[p_cols_eur]))
  if (all(is.na(p_vals))) {
    p_vals <- suppressWarnings(as.numeric(row[p_cols_multi]))
  }
  if (all(is.na(p_vals))) return(NA_character_)
  ph_labels[which.min(p_vals)]
}

idx <- gwas_dat$MarkerName %in% list_ids
gwas_dat$phenotype[idx] <- apply(gwas_dat[idx, ], 1, assign_phenotype_fallback)

table(gwas_dat$phenotype,useNA="ifany")
            #  CD   IBD_saturated IBD_unsaturated              UC 
            # 128             104             269             118 
            
dim(gwas_dat)
# [1] 619 177

fwrite(gwas_dat,paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_pheno.tsv.gz",sep=""),
col.names=T,row.names=F,quote=F,sep="\t")

q("no")
