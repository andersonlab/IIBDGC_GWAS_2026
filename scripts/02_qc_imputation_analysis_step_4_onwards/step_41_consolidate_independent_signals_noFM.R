# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# ############################################################################################################################
# ############################################################################################################################

# # two different versions of FM resutls stored here:
# # /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/fine_mapping/

# # singularity exec iibdgc_postprocess_10_singularity.sif

# MEM=5000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

# library(data.table)
# rm(list=ls())

# path_gwas<-"/path/to/ibdgwas/IIBDGC/"


# all<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal.tsv.gz",sep=""))
# all<-as.data.frame(all)
# dim(all)
# # [1] 1374  160

# length(names(table(all$updated_region)))
# # [1] 420


# all$class[which(all$updated_region=="2_229716690_231482974")]<-"old"


# list_ids<-all$MarkerName
# list_ids<-list_ids[!duplicated(list_ids)]
# length(list_ids)
# # [1] 1374

# write.table(list_ids,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_variants_for_ld_with_old_data.tsv",sep=""),col.names=T,row.names=F,quote=F,sep="\t")


# q("no")

# #############################################################################

# path_gwas=/path/to/ibdgwas/IIBDGC/
# ph=ibd

# # note this will creat hardcalls, and those will be used for ld estimation, all variants good imputation, it should not be biased...

# MEM=2000
# for chr in {1..22} 
# do
# bsub -J"cojo_i" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
# -o ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_ld1_stdout \
# -e ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_ld1_stderr \
# "/path/to/software/username/./plink2 \
# --bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_subset_included_in_${ph}_analysis \
# --extract ${path_gwas}post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_variants_for_ld_with_old_data.tsv \
# --make-bed --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_list_tier_2_unsupervised_conditional_variants_for_ld_with_old_data"
# done

# for chr in {1..22} 
# do
# echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_ld1_stdout | grep "Successfully"
# done

# for chr in {1..22} 
# do
# bsub -J"cojo_i" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
# -o ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_ld2_stdout \
# -e ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_ld2_stderr \
# "/path/to/software/username/plink_linux_x86_64_20181202/./plink \
# --bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_list_tier_2_unsupervised_conditional_variants_for_ld_with_old_data \
# --r2 --ld-window-kb '1000000000' --ld-window '1000000000' --ld-window-r2 '0' --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_list_tier_2_unsupervised_conditional_variants_for_ld_with_old_data"
# done

# for chr in {1..22} 
# do
# echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_ld2_stdout | grep "Successfully"
# done


#################################################################################################################################################

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=5000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"


# Helper: load plink LD files for all autosomes, filter by R2 threshold
load_ld <- function(path_gwas, r2_min = 0) {
  ld_list <- vector("list", 22)
  for (chr in 1:22) {
    ld_tmp <- fread(paste0(path_gwas,
      "post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr", chr,
      "_list_tier_2_unsupervised_conditional_variants_for_ld_with_old_data.ld"))
    if (r2_min > 0) ld_tmp <- ld_tmp[R2 >= r2_min]
    ld_list[[chr]] <- ld_tmp
  }
  rbindlist(ld_list)
}

# Helper: load COJO joint model results across all autosomes and phenotypes
load_cojo <- function(path_gwas, suffix, cols = c("SNP", "bJ", "bJ_se", "pJ", "LD_r")) {
  pheno <- c("ibd", "cd", "uc")
  dat <- NULL
  for (ph in pheno) {
    dat_list <- vector("list", 22)
    for (chr in 1:22) {
      dat_list[[chr]] <- fread(paste0(path_gwas,
        "post_imputation/2022/analysis/conditional_analysis/final_independent_model/",
        chr, "_", ph,
        "_independent_index_variants_eur_tier_2_ld_allarrays_conditioning_on_final_independent_model_0.1_no_subset",
        suffix, ".jma.cojo"))
    }
    dat_tmp <- rbindlist(dat_list)[, ..cols]
    sel_cols <- cols[cols != "SNP"]
    setnames(dat_tmp, sel_cols, paste(sel_cols, ph, sep = "_"))
    print(dim(dat_tmp))
    if (is.null(dat)) {
      dat <- dat_tmp
    } else {
      dat <- merge(dat, dat_tmp, by = "SNP", all = TRUE)
    }
  }
  dat
}

all <- as.data.frame(fread(paste0(path_gwas, "post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal.tsv.gz")))
dim(all)
# [1] 1374  160

all$class[which(all$updated_region=="2_229716690_231482974")]<-"old"

length(unique(all$updated_region))
# [1] 420

reg<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_regions_420_cojo_supervised_plus_unsupervised.tsv.gz"))

# fm file includes regions later excluded (post FM) - update accordingly:
fm_reg<-fread("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/fine_mapping/fine_mapping_region_reclassification.csv")
dim(fm_reg)
fm_reg$region_label[which(fm_reg$region_label=="2_229716690_230796312")]<-"2_229716690_231482974"
table(fm_reg$reason_noCS)
  #                                .              classA_no_significant 
  #                              170                                  1 
  #                    classB_purity classC1_nosignificant_after_filter 
  #                                1                                231 
  # classC2_significant_after_filter 
  #                               16 


all$region_fm_2025<-NA
all$region_fm_2025[which(all$updated_region %in% fm_reg$region_label[which(fm_reg$reason_noCS %in% c("classC1_nosignificant_after_filter","classA_no_significant"))])]<-"N_filtered_out_preFM"
all$region_fm_2025[which(all$updated_region %in% fm_reg$region_label[which(fm_reg$reason_noCS %in% c("classB_purity","classC2_significant_after_filter"))])]<-"N_filtered_out_postFM"
all$region_fm_2025[which(all$updated_region %in% fm_reg$region_label[which(fm_reg$reason_noCS %in% c("."))])]<-"Y"

table(all$region_fm_2025,useNA="ifany")

table(all$updated_region[which(is.na(all$region_fm_2025))],all$class_signal[which(is.na(all$region_fm_2025))])
  #                        in_ld_new_cojo_unsupervised new_cojo_unsupervised
  # 12_8422106_10257536                              1                     2
  # 14_106197560_107197560                           0                     1
  # 16_69127624_70127624                             0                     1
  # 18_21986192_22986192                             0                     1
  # 4_41896369_42896369                              1                     5
  # 9_93754276_94754276                              1                     1

all$region_fm_2025[which(is.na(all$region_fm_2025))]<-"N_filtered_out_preFM"

table(all$region_fm_2025,useNA="ifany")
# N_filtered_out_postFM  N_filtered_out_preFM                     Y 
#                    65                   438                   871  

########################################################
# retain any signal from new_cojo_unsupervised

length(unique(all$updated_region))
# [1] 420

table(all$class_signal[which(all$class!="old")])
#                             new_cojo_unsupervised   # report as new
#                                               150 
# new_cojo_unsupervised_in_ld_new_cojo_unsupervised   # do not report
#                                                40 
#              new_cojo_unsupervised_in_ld_with_old   # report as known
#                                                 3 
#             new_signal_new_region_cojo_supervised  # do not report
#                                                24 


table(all$class_signal[which(all$class=="old")])
#                new_cojo_supervised_in_ld_with_old   # do not report
#                                                19 
#                             new_cojo_unsupervised   # report as new
#                                               284 
# new_cojo_unsupervised_in_ld_new_cojo_unsupervised   # do not report
#                                                97 
#              new_cojo_unsupervised_in_ld_with_old   # report as known
#                                               367 
#                                               old   # do not report
#                                               317 
#                         old|new_cojo_unsupervised   # report as known
#                                                73 


# SUBSET OF SIGNALS 'DO NOT REPORT'
all_tmp<-all[which(!all$class_signal %in% c("new_cojo_unsupervised","old|new_cojo_unsupervised","new_cojo_unsupervised_in_ld_with_old")),]
dim(all_tmp)
# [1] 497 161


# SUBSET OF SIGNALS TO REPORT
all<-all[which(all$class_signal %in% c("new_cojo_unsupervised","old|new_cojo_unsupervised","new_cojo_unsupervised_in_ld_with_old")),]

table(all$class_signal,useNA="ifany")
          #      new_cojo_unsupervised new_cojo_unsupervised_in_ld_with_old 
          #                        434                                  370 
          #  old|new_cojo_unsupervised 
          #                         73 
all$class_signal_final<-NA

all$class_signal_final[which(all$class_signal %in% c("old|new_cojo_unsupervised","new_cojo_unsupervised_in_ld_with_old"))]<-"new_cojo_unsupervised_known_signal"
all$class_signal_final[which(all$class_signal %in% c("new_cojo_unsupervised"))]<-"new_cojo_unsupervised_new_signal"
table(all$class_signal_final)
# new_cojo_unsupervised_known_signal   new_cojo_unsupervised_new_signal 
#                                443                                434

table(all$class_signal,all$class)
  #                                      new old
  # new_cojo_unsupervised                150 284
  # new_cojo_unsupervised_in_ld_with_old   3 367
  # old|new_cojo_unsupervised              0  73


table(all$class_signal_final,all$class)                           
  #                                    new old
  # new_cojo_unsupervised_known_signal   3 440
  # new_cojo_unsupervised_new_signal   150 284



###############################################################################################
# table to report final signals from our analyses:


table(all$class_signal_final,useNA="ifany")
# new_cojo_unsupervised_known_signal   new_cojo_unsupervised_new_signal 
#                                443                                434 


table(all$class,all$region_fm_2025,useNA="ifany")
  #     N_filtered_out_postFM N_filtered_out_preFM   Y
  # new                     6                  125  22
  # old                    35                  120 569

length(names(table(all$updated_region)))
# [1] 352

length(names(table(all$updated_region[which(all$region_fm_2025=="Y")])))
# [1] 165
length(names(table(all$updated_region[which(all$region_fm_2025 %in% c("N_filtered_out_postFM","N_filtered_out_preFM"))])))
# [1] 187

length(names(table(all$updated_region[which(all$region_fm_2025 %in% c("N_filtered_out_postFM"))])))
# [1] 17
length(names(table(all$updated_region[which(all$region_fm_2025 %in% c("N_filtered_out_preFM"))])))
# [1] 170

dim(reg[which(!reg$updated_region %in% all$updated_region)])
# [1] 68   2

table(fm_reg$reason_noCS[!fm_reg$region_label %in% all$updated_region])
#                                  .              classA_no_significant 
#                                  5                                  1 
# classC1_nosignificant_after_filter 
#                                 67 

dim(all[which(all$region_fm_2025=="Y"),])
# [1] 591  213

dim(all[which(all$region_fm_2025!="Y"),])
# [1] 286 213


all$chr<-gsub("chr","",all$chr)

fwrite(all,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_no_overlap.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")

all_keep<-all

table(all_keep$class_signal)
          #      new_cojo_unsupervised new_cojo_unsupervised_in_ld_with_old 
          #                        434                                  370 
          #  old|new_cojo_unsupervised 
          #                         73


############################################
# annotate the remaining regions

dim(reg[which(reg$updated_region %in% all_keep$updated_region),])
# [1] 352   

reg_tmp<-reg$updated_region[which(!reg$updated_region %in% all_keep$updated_region)]
length(reg_tmp)
# [1] 68

length(table(all_tmp$updated_region[which(all_tmp$updated_region %in% reg_tmp)]))
# [1] 68


all<-all_tmp[which(all_tmp$updated_region %in% reg_tmp),]

table(all$class_signal)
all$class_signal[which(all$class_signal=="new_signal_new_region_cojo_supervised")]<-"new_cojo_supervised"
table(all$class_signal)
              #  new_cojo_supervised new_cojo_supervised_in_ld_with_old 
              #                   24                                  1 
              #                  old 
              #                   44 


all$"P-value_cd_eur_tier_1"<-as.numeric(all$"P-value_cd_eur_tier_1")
all$"P-value_uc_eur_tier_1"<-as.numeric(all$"P-value_uc_eur_tier_1")
all$"P-value_ibd_eur_tier_1"<-as.numeric(all$"P-value_ibd_eur_tier_1")


all$"P-value_cd_eur_tier_2"<-as.numeric(all$"P-value_cd_eur_tier_2")
all$"P-value_uc_eur_tier_2"<-as.numeric(all$"P-value_uc_eur_tier_2")
all$"P-value_ibd_eur_tier_2"<-as.numeric(all$"P-value_ibd_eur_tier_2")


all$"P-value_cd_eur_tier_2_eas_sas"<-as.numeric(all$"P-value_cd_eur_tier_2_eas_sas")
all$"P-value_uc_eur_tier_2_eas_sas"<-as.numeric(all$"P-value_uc_eur_tier_2_eas_sas")
all$"P-value_ibd_eur_tier_2_eas_sas"<-as.numeric(all$"P-value_ibd_eur_tier_2_eas_sas")

all$"P-value_association_cd_eur_tier_2_eas_sas_mrmega"<-as.numeric(all$"P-value_association_cd_eur_tier_2_eas_sas_mrmega")
all$"P-value_association_uc_eur_tier_2_eas_sas_mrmega"<-as.numeric(all$"P-value_association_uc_eur_tier_2_eas_sas_mrmega")
all$"P-value_association_ibd_eur_tier_2_eas_sas_mrmega"<-as.numeric(all$"P-value_association_ibd_eur_tier_2_eas_sas_mrmega")


all$genomewide_significance_eur_tier_1<-"N"
all$genomewide_significance_eur_tier_1[which(all$"P-value_cd_eur_tier_1"<5E-8 | all$"P-value_uc_eur_tier_1"<5E-8 | all$"P-value_ibd_eur_tier_1"<5E-8)]<-"Y"
table(all$genomewide_significance_eur_tier_1)
#  N  Y 
# 66  3  

all$genomewide_significance_eur_tier_2<-"N"
all$genomewide_significance_eur_tier_2[which(all$"P-value_cd_eur_tier_2"<5E-8 | all$"P-value_uc_eur_tier_2"<5E-8 | all$"P-value_ibd_eur_tier_2"<5E-8)]<-"Y"
table(all$genomewide_significance_eur_tier_2)
#  N  Y 
# 65  4 

all$genomewide_significance_eur_eas_sas_tier_2<-"N"
all$genomewide_significance_eur_eas_sas_tier_2[which(all$"P-value_cd_eur_tier_2_eas_sas"<5E-8 | all$"P-value_uc_eur_tier_2_eas_sas"<5E-8 | all$"P-value_ibd_eur_tier_2_eas_sas"<5E-8 |
all$"P-value_association_cd_eur_tier_2_eas_sas_mrmega"<5E-8 | all$"P-value_association_uc_eur_tier_2_eas_sas_mrmega"<5E-8 | all$"P-value_association_ibd_eur_tier_2_eas_sas_mrmega"<5E-8)]<-"Y"
table(all$genomewide_significance_eur_eas_sas_tier_2)
#  N  Y 
# 42 27 

table(all$genomewide_significance_eur_tier_2,all$genomewide_significance_eur_tier_1,useNA="ifany")
  #    N  Y
  # N 63  2
  # Y  3  1


table(all$genomewide_significance_eur_tier_2,all$genomewide_significance_eur_eas_sas_tier_2,useNA="ifany")
  #    N  Y
  # N 42 23
  # Y  0  4

### create one flag per region:

all$region_genomewide_significance_eur_tier_1<-"N"
all$region_genomewide_significance_eur_tier_2<-"N"
all$region_genomewide_significance_eur_eas_sas_tier_2<-"N"

regions<-all$updated_region
regions<-regions[!duplicated(regions)]
length(regions)
# [1] 68

for (i in 1:length(regions)) {

  tmp<-all[which(all$updated_region==regions[i]),]

  if (nrow(tmp[which(tmp$genomewide_significance_eur_tier_1=="Y"),])>0) {
    all$region_genomewide_significance_eur_tier_1[which(all$updated_region==regions[i])]<-"Y"
  } 

  if (nrow(tmp[which(tmp$genomewide_significance_eur_tier_2=="Y"),])>0) {
    all$region_genomewide_significance_eur_tier_2[which(all$updated_region==regions[i])]<-"Y"
  }

  if (nrow(tmp[which(tmp$genomewide_significance_eur_eas_sas_tier_2=="Y"),])>0) {
    all$region_genomewide_significance_eur_eas_sas_tier_2[which(all$updated_region==regions[i])]<-"Y"
  }

}

# regions with significant variants in eur tier 1
length(names(table(all$updated_region[which(all$genomewide_significance_eur_tier_1=="Y")])))
# [1] 3

# regions with significant variants in eur tier 2
length(names(table(all$updated_region[which(all$genomewide_significance_eur_tier_2=="Y")])))
# [1] 4

# regions with significant variants in eur_eas_sas tier 2
length(names(table(all$updated_region[which(all$region_genomewide_significance_eur_eas_sas_tier_2=="Y")])))
# [1] 26

all$class_signal_final<-NA

fwrite(all,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_new_and_old_index_gwas_variants_regions_no_conditional_significant_eur_tier2_no_fm_imput.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")



######################################################
# recapture those that reach genome wide significance in multiancestry data

all_recap<-all[which(all$genomewide_significance_eur_eas_sas_tier_2=="Y"),colnames(all_keep)]
table(all_recap$class_signal)
              #  new_cojo_supervised new_cojo_supervised_in_ld_with_old 
              #                   19                                  1 
              #                  old 
              #                    7 

# remove new_cojo_supervised_in_ld_with_old
all_recap<-all_recap[which(all_recap$class_signal!="new_cojo_supervised_in_ld_with_old"),]

all_recap$class_signal_final[which(all_recap$class_signal=="old")]<-"new_cojo_supervised_gw_significant_multiancestry_known_signal"
all_recap$class_signal_final[which(all_recap$class_signal=="new_cojo_supervised")]<-"new_cojo_supervised_gw_significant_multiancestry_new_signal"
table(all_recap$class_signal_final)
# new_cojo_supervised_gw_significant_multiancestry_known_signal 
#                                                             7 
#   new_cojo_supervised_gw_significant_multiancestry_new_signal 
#                                                            19

# variant not present in EUR
all_recap<-all_recap[which(!all_recap$MarkerName %in% "chr12:110259525:G:A"),]

table(all_recap$class_signal_final)
# new_cojo_supervised_gw_significant_multiancestry_known_signal 
#                                                             6 
#   new_cojo_supervised_gw_significant_multiancestry_new_signal 
#                                                            19 

ld <- load_ld(path_gwas, r2_min = 0.1)

ld<-ld[which( (ld$SNP_A %in% all_recap$MarkerName) | (ld$SNP_B %in% all_recap$MarkerName)),]
ld
#      CHR_A      BP_A              SNP_A CHR_B      BP_B              SNP_B
# 3312     5 130440872 chr5:130440872:C:T     5 132143203 chr5:132143203:T:C
# 3313     5 130440872 chr5:130440872:C:T     5 132290918 chr5:132290918:C:T
# 3314     5 130440872 chr5:130440872:C:T     5 132336964 chr5:132336964:T:C
# 3317     5 130440872 chr5:130440872:C:T     5 132435113 chr5:132435113:C:T
# 3318     5 130440872 chr5:130440872:C:T     5 132461230 chr5:132461230:G:A
#            R2
# 3312 0.178712
# 3313 0.110882
# 3314 0.349184
# 3317 0.312514
# 3318 0.111288

# chr5:130440872:C:T in LD with other

dim(all_recap[which(all_recap$MarkerName %in% ld$SNP_A),])
# [1]   1 162 
dim(all_recap[which(all_recap$MarkerName %in% ld$SNP_B),])
# [1]   0 162

all_recap[which(all_recap$MarkerName %in% ld$SNP_A),c("updated_region","class","class_signal")]
#             updated_region class class_signal
# 1014 5_129940872_130940873   old          old

dim(all_keep[which(all_keep$MarkerName %in% ld$SNP_A),])
# [1]   0 162 
dim(all_keep[which(all_keep$MarkerName %in% ld$SNP_B),])
# [1]   4 162 - why those 4????

all_keep[which(all_keep$MarkerName %in% ld$SNP_B),c("updated_region","class")]
#             updated_region class
# 1015 5_131190525_132974928   old
# 1017 5_131190525_132974928   old
# 1020 5_131190525_132974928   old
# 1021 5_131190525_132974928   old

all_keep[which(all_keep$MarkerName %in% ld$SNP_B),c("updated_region","class_signal")]
#             updated_region                         class_signal
# 1015 5_131190525_132974928 new_cojo_unsupervised_in_ld_with_old
# 1017 5_131190525_132974928 new_cojo_unsupervised_in_ld_with_old
# 1020 5_131190525_132974928            old|new_cojo_unsupervised
# 1021 5_131190525_132974928 new_cojo_unsupervised_in_ld_with_old

132143203-130440872
# [1] 1702331 > 1Mb away


all_recap<-rbind(all_recap,all_keep)

table(all_recap$class_signal)
#                  new_cojo_supervised                new_cojo_unsupervised 
#                                   19                                  434 
# new_cojo_unsupervised_in_ld_with_old                                  old 
#                                  370                                    6 
#            old|new_cojo_unsupervised 
#                                   73 


# DOUBLE CHECK THAT THE FINAL LIST IS ROBUST:


# SEVEN of them to exclude: no significant at all in EUR tier1 only signif when tier 2 is included:

all_recap$'P-value_ibd_eur_tier_1'<-as.numeric(all_recap$'P-value_ibd_eur_tier_1')
all_recap$'P-value_cd_eur_tier_1'<-as.numeric(all_recap$'P-value_cd_eur_tier_1')
all_recap$'P-value_uc_eur_tier_1'<-as.numeric(all_recap$'P-value_uc_eur_tier_1')

all_recap$'P-value_ibd_eur_tier_2'<-as.numeric(all_recap$'P-value_ibd_eur_tier_2')
all_recap$'P-value_cd_eur_tier_2'<-as.numeric(all_recap$'P-value_cd_eur_tier_2')
all_recap$'P-value_uc_eur_tier_2'<-as.numeric(all_recap$'P-value_uc_eur_tier_2')


all_recap$rate_Neff_ibd_eur_tier_2<-as.numeric(all_recap$rate_Neff_ibd_eur_tier_2)
all_recap$rate_Neff_cd_eur_tier_2<-as.numeric(all_recap$rate_Neff_cd_eur_tier_2)
all_recap$rate_Neff_uc_eur_tier_2<-as.numeric(all_recap$rate_Neff_uc_eur_tier_2)



all_recap[which( 
(all_recap$'P-value_ibd_eur_tier_1' > 0.01 & all_recap$'P-value_ibd_eur_tier_2' < 5E-8 & all_recap$rate_Neff_ibd_eur_tier_2<0.5) | 
(all_recap$'P-value_cd_eur_tier_1' > 0.01 & all_recap$'P-value_cd_eur_tier_2' < 5E-8 & all_recap$rate_Neff_cd_eur_tier_2<0.5) |
(all_recap$'P-value_uc_eur_tier_1' > 0.01 & all_recap$'P-value_uc_eur_tier_2' < 5E-8 & all_recap$rate_Neff_uc_eur_tier_2<0.5) ),
c('MarkerName','P-value_ibd_eur_tier_1','P-value_ibd_eur_tier_2','P-value_cd_eur_tier_1','P-value_cd_eur_tier_2','P-value_uc_eur_tier_1','P-value_uc_eur_tier_2',
"rate_Neff_ibd_eur_tier_2","rate_Neff_cd_eur_tier_2","rate_Neff_uc_eur_tier_2")]
#               MarkerName P-value_ibd_eur_tier_1 P-value_ibd_eur_tier_2
# 188 chr10:105644358:TA:T                0.65780              2.887e-09
# 193  chr10:125821742:T:G                0.03577              1.094e-15
# 221   chr10:60750152:A:G                0.07403              5.847e-18
# 438   chr15:34731239:A:C                0.27510              2.751e-01
# 994    chr4:42396369:G:T                0.63690              2.645e-11
#     P-value_cd_eur_tier_1 P-value_cd_eur_tier_2 P-value_uc_eur_tier_1
# 188               0.23250             2.244e-06               0.35930
# 193               0.04845             1.148e-12               0.06247
# 221               0.16390             8.196e-13               0.18570
# 438               0.11490             3.080e-10               0.58850
# 994               0.56470             1.809e-05               0.16220
#     P-value_uc_eur_tier_2 rate_Neff_ibd_eur_tier_2 rate_Neff_cd_eur_tier_2
# 188             2.858e-05                0.4659048               0.6721715
# 193             1.980e-10                0.4253783               0.6095054
# 221             5.362e-14                0.4091426               0.5928558
# 438             6.038e-13                0.2263522               0.6193028
# 994             2.250e-07                0.4111789               0.5982854
#     rate_Neff_uc_eur_tier_2
# 188               0.4470857
# 193               0.4463048
# 221               0.4163770
# 438               0.4496562
# 994               0.4168224

to_exclude1<-all_recap$MarkerName[which( 
(all_recap$'P-value_ibd_eur_tier_1' > 0.01 & all_recap$'P-value_ibd_eur_tier_2' < 5E-8 & all_recap$rate_Neff_ibd_eur_tier_2<0.5) | 
(all_recap$'P-value_cd_eur_tier_1' > 0.01 & all_recap$'P-value_cd_eur_tier_2' < 5E-8 & all_recap$rate_Neff_cd_eur_tier_2<0.5) |
(all_recap$'P-value_uc_eur_tier_1' > 0.01 & all_recap$'P-value_uc_eur_tier_2' < 5E-8 & all_recap$rate_Neff_uc_eur_tier_2<0.5) )]

to_exclude2<-all_recap$MarkerName[which( 
  (is.na(all_recap$'P-value_ibd_eur_tier_1') & all_recap$'P-value_ibd_eur_tier_2' < 5E-8 & all_recap$rate_Neff_ibd_eur_tier_2<0.5) | 
(is.na(all_recap$'P-value_cd_eur_tier_1') & all_recap$'P-value_cd_eur_tier_2' < 5E-8 & all_recap$rate_Neff_cd_eur_tier_2<0.5) |
(is.na(all_recap$'P-value_uc_eur_tier_1') & all_recap$'P-value_uc_eur_tier_2' < 5E-8 & all_recap$rate_Neff_uc_eur_tier_2<0.5) )]


# "chr14:105968520:G:A" excluded from all tier 1 cohorts for UC, imputation <0.3 in gsa
to_exclude3<-c("chr2:18792705:C:CT","chr9:19972835:AG:A","chr4:42635348:ATATC:A","chr14:105963206:C:A","chr14:105968520:G:A")

to_exclude<-c(to_exclude1,to_exclude2,to_exclude3)
length(to_exclude)
# [1] 10

dim(all_recap)
# [1] 902 162

all_recap<-all_recap[which(!all_recap$MarkerName %in% to_exclude),]
dim(all_recap)
# [1] 892 162


table(all_recap$class_signal_final,all_recap$class,useNA="ifany") 
  #                                                               new old
  # new_cojo_supervised_gw_significant_multiancestry_known_signal   0   6  # known signal
  # new_cojo_supervised_gw_significant_multiancestry_new_signal    18   0  # new 
  # new_cojo_unsupervised_known_signal                              3 440  # known signal
  # new_cojo_unsupervised_new_signal                              144 281  # new 



# double check that no signal is in strong LD with each other:

ld <- load_ld(path_gwas, r2_min = 0.1)


ld<-ld[which( (ld$SNP_A %in% all_recap$MarkerName) & (ld$SNP_B %in% all_recap$MarkerName)),]
dim(ld)
# [1] 142   7 # at r2 = 0.6, chr 20 issue with colinearity
# [1] 281   7


ldtest<-c(ld$SNP_A,ld$SNP_B)
ldtest<-ldtest[!duplicated(ldtest)]
length(ldtest)
# [1] 370

table(all_recap$class_signal[which(all_recap$MarkerName %in% ldtest)])
              #  new_cojo_unsupervised new_cojo_unsupervised_in_ld_with_old 
              #                     30                                  290 
              #                    old            old|new_cojo_unsupervised 
              #                      1                                   49 

all_recap$exclude<-NA

for (i in 1:length(ldtest)) {

  print(i)

  ld_tmp<-ld[which( (ld$SNP_A %in% ldtest[i]) | (ld$SNP_B %in% ldtest[i])),]
  print(ld_tmp$R2)

  # retain the one with the smallest conditional pvalue:
  df<-all_recap[which(all_recap$MarkerName %in% c(ld_tmp$SNP_A,ld_tmp$SNP_B)),c("MarkerName","Pvalue_cond_cd","Pvalue_cond_uc","Pvalue_cond_ibd")]
  df$min<-apply(df[2:4], 1, min, na.rm = TRUE)

  to_retain<-df$MarkerName[which(df$min==min(df$min,na.rm=T))]
  
  # exclude any other variant
  print(all_recap$exclude[which(all_recap$MarkerName %in% df$MarkerName)])
  all_recap$exclude[which(all_recap$MarkerName %in% df$MarkerName[which(df$MarkerName!=to_retain)])]<-1
  all_recap$exclude[which(all_recap$MarkerName %in% df$MarkerName[which(df$MarkerName==to_retain)])]<-0
  print(all_recap$exclude[which(all_recap$MarkerName %in% df$MarkerName)])

}

table(all_recap$class_signal_final,all_recap$exclude)
  #                                                                 0   1
  # new_cojo_supervised_gw_significant_multiancestry_known_signal   0   1
  # new_cojo_supervised_gw_significant_multiancestry_new_signal     0   0
  # new_cojo_unsupervised_known_signal                            136 203
  # new_cojo_unsupervised_new_signal                               21   9

table(all_recap$class_signal_final)
# new_cojo_supervised_gw_significant_multiancestry_known_signal 
#                                                             6 
#   new_cojo_supervised_gw_significant_multiancestry_new_signal 
#                                                            18 
#                            new_cojo_unsupervised_known_signal 
#                                                           443 
#                              new_cojo_unsupervised_new_signal 
#                                                           425 

dim(all_recap)
# [1] 892 163

all_recap<-all_recap[which(is.na(all_recap$exclude) | all_recap$exclude==0),]
dim(all_recap)
# [1] 679 163

table(all_recap$class_signal_final)
# new_cojo_supervised_gw_significant_multiancestry_known_signal 
#                                                             5 
#   new_cojo_supervised_gw_significant_multiancestry_new_signal 
#                                                            18 
#                            new_cojo_unsupervised_known_signal 
#                                                           240 
#                              new_cojo_unsupervised_new_signal 
#                                                           416 

table(all_recap$class_signal)
#                  new_cojo_supervised                new_cojo_unsupervised 
#                                   18                                  416 
# new_cojo_unsupervised_in_ld_with_old                                  old 
#                                  193                                    5 
#            old|new_cojo_unsupervised 
#                                   47 

dim(all_recap)
# [1] 679 163


## exclude the variants tagging the TNRC18 variant 7-5397122-C-T variant

# rs61197566 - chr7:5139807:G:A (best LD proxi in 1000GP finnish)
finngen<-c("chr7:5420919:C:T","chr7:4468400:T:G","chr7:4827597:G:C","chr7:4881854:T:C","chr7:4981470:C:G","chr7:5139807:G:A","chr7:5147785:G:A","chr7:5266889:G:T","chr7:5285792:A:T","chr7:5319942:C:T","chr7:5333452:G:A","chr7:5433979:G:C","chr7:5453006:C:T","chr7:5611768:C:T","chr7:6161112:CCAAT:C","chr7:6352123:G:C")

# retain only chr7:5433979:G:C - more significant
finngen<-finngen[which(!finngen %in% "chr7:5433979:G:C")]
all_recap<-all_recap[which(!all_recap$MarkerName %in% finngen),]

dim(all_recap)
# [1] 664 163

fwrite(all_recap,paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_no_tags.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")

### save the list of variants to run a join model:

table(all_recap$chr)
#  1 10 11 12 13 14 15 16 17 18 19  2 20 21 22  3  4  5  6  7  8  9 
# 85 31 34 23 14 20 17 41 32 10 28 59 24 11 15 28 24 41 34 57 25 28 


write.table(all_recap[,c("MarkerName")],paste(path_gwas,"post_imputation/2022/analysis/conditional_analysis/final_independent_model/final_independent_model_no_tags.snplist",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")


q("no")


####################################################
# estimate the variant effec in the joint model - but retain all variants, no subset by Neff>=0.5 and hetPval >=1E-15 


echo "SNP A1 A2 freq b se p N" \
> ${path_gwas}post_imputation/2022/analysis/conditional_analysis/header.ma

MEM=800
path_gwas=/path/to/ibdgwas/IIBDGC/
pheno=(ibd cd uc)

for chr in {1..22}
do
for ph in ${pheno[@]}
do 
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 4 \
-e ${path_gwas}post_imputation/log/${chr}_${ph}_cojo_sumstats_eur_tier2_stderr \
-o ${path_gwas}post_imputation/log/${chr}_${ph}_cojo_sumstats_eur_tier2_stdout \
"zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/${ph}/${chr}_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz \
 | awk -F ' ' '\$23 >= 0.2 { print }' | awk -F '\t' '{ print \$1,\$4,\$3,\$16,\$5,\$6,\$7,\$15}' > \
${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/${chr}_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format_no_subset.ma"
done
done

for chr in {1..22}
do 
echo $chr && for ph in ${pheno[@]}
do
echo $ph && tail -50 ${path_gwas}post_imputation/log/${chr}_${ph}_cojo_sumstats_eur_tier2_stdout | grep -E "Successfully|Exited"
done
done

for ph in ${pheno[@]}
do
ls -la ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/*_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format_no_subset.ma | wc -l
done

rm ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allchr_*_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format_no_subset.ma

for ph in ${pheno[@]}
do
wc -l ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/*_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format_no_subset.ma
done
# 18646697 total
# 18861653 total
# 18832762 total

for ph in ${pheno[@]}
do
cat ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/*_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format_no_subset.ma \
| grep -v "MarkerName" \
> ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allchr_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format_no_subset_tmp.ma
done

for ph in ${pheno[@]}
do
ls -la ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allchr_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format_no_subset_tmp.ma
done

for ph in ${pheno[@]}
do
rm ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/*_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format_no_subset.ma
done


for ph in ${pheno[@]}
do 
cat ${path_gwas}post_imputation/2022/analysis/conditional_analysis/header.ma  ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allchr_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format_no_subset_tmp.ma \
> ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allchr_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format_no_subset.ma
done

for ph in ${pheno[@]}
do 
rm ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allchr_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format_no_subset_tmp.ma
done


for ph in ${pheno[@]}
do echo ${ph} && wc -l ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allchr_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format_no_subset.ma
done

# ibd
# 18646676 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/allchr_ibd_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format_no_subset.ma
# cd
# 18861632 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/allchr_cd_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format_no_subset.ma
# uc
# 18832741 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/allchr_uc_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format_no_subset.ma



path_gwas=/path/to/ibdgwas/IIBDGC/

pheno=(ibd cd uc)
MEM=2500

for ph in ${pheno[@]}
do
for chr in {1..22}
do
bsub -J"cojo_i" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" \
-o ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_conditioning_final_model_ld_allarrays_no_subset_stdout \
-e ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_conditioning_final_model_ld_allarrays_no_subset_stderr \
"/path/to/software/username/gcta-1.94.1-linux-kernel-3-x86_64/./gcta64 \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_subset_included_in_${ph}_analysis \
--chr ${chr} \
--extract ${path_gwas}post_imputation/2022/analysis/conditional_analysis/final_independent_model/final_independent_model_no_tags.snplist \
--cojo-file ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allchr_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format_no_subset.ma \
--cojo-joint \
--out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/final_independent_model/${chr}_${ph}_independent_index_variants_eur_tier_2_ld_allarrays_conditioning_on_final_independent_model_0.1_no_subset"
done
done


for ph in ${pheno[@]}
do
echo ${ph} && for chr in {1..22}
do 
echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_conditioning_final_model_ld_allarrays_no_subset_stdout | grep -E "Successfully|Exited"
done
done




################################################################

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=2000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
library(ggplot2)

rm(list=ls())

# Helper: load plink LD files for all autosomes, filter by R2 threshold
load_ld <- function(path_gwas, r2_min = 0) {
  ld_list <- vector("list", 22)
  for (chr in 1:22) {
    ld_tmp <- fread(paste0(path_gwas,
      "post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr", chr,
      "_list_tier_2_unsupervised_conditional_variants_for_ld_with_old_data.ld"))
    if (r2_min > 0) ld_tmp <- ld_tmp[R2 >= r2_min]
    ld_list[[chr]] <- ld_tmp
  }
  rbindlist(ld_list)
}

# Helper: load COJO joint model results across all autosomes and phenotypes
load_cojo <- function(path_gwas, suffix, cols = c("SNP", "bJ", "bJ_se", "pJ", "LD_r")) {
  pheno <- c("ibd", "cd", "uc")
  dat <- NULL
  for (ph in pheno) {
    dat_list <- vector("list", 22)
    for (chr in 1:22) {
      dat_list[[chr]] <- fread(paste0(path_gwas,
        "post_imputation/2022/analysis/conditional_analysis/final_independent_model/",
        chr, "_", ph,
        "_independent_index_variants_eur_tier_2_ld_allarrays_conditioning_on_final_independent_model_0.1_no_subset",
        suffix, ".jma.cojo"))
    }
    dat_tmp <- rbindlist(dat_list)[, ..cols]
    sel_cols <- cols[cols != "SNP"]
    setnames(dat_tmp, sel_cols, paste(sel_cols, ph, sep = "_"))
    print(dim(dat_tmp))
    if (is.null(dat)) {
      dat <- dat_tmp
    } else {
      dat <- merge(dat, dat_tmp, by = "SNP", all = TRUE)
    }
  }
  dat
}


path_gwas<-"/path/to/ibdgwas/IIBDGC/"
pheno<-c("ibd","cd","uc")

all_recap<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_no_tags.tsv.gz",sep=""))
dim(all_recap)
# [1] 664 163

dat <- load_cojo(path_gwas, suffix = "")

all_recap[which(!all_recap$MarkerName %in% dat$SNP),]
# 0

all_recap<-merge(all_recap,dat,by.x="MarkerName",by.y="SNP")
dim(all_recap)
# [1] 664 175

ibd<-all_recap[,c("MarkerName","P-value_ibd_eur_tier_2","pJ_ibd","class_signal")]
colnames(ibd)<-c("MarkerName","p","pJ","class_signal")
ibd$pheno<-"ibd"
cd<-all_recap[,c("MarkerName","P-value_cd_eur_tier_2","pJ_cd","class_signal")]
colnames(cd)<-c("MarkerName","p","pJ","class_signal")
cd$pheno<-"cd"
uc<-all_recap[,c("MarkerName","P-value_uc_eur_tier_2","pJ_uc","class_signal")]
colnames(uc)<-c("MarkerName","p","pJ","class_signal")
uc$pheno<-"uc"

all<-rbind(ibd,cd,uc)
table(all$pheno)
#  cd ibd  uc 
# 664 664 664

all$p<-as.numeric(all$p)
head(all)

all<-all[which(all$p<=5E-8 & !is.na(all$pJ)),]
table(all$pheno)
#  cd ibd  uc 
# 333 472 333 

# ggplot(all, aes(x=-log10(pJ), y=-log10(p),color=all$class_signal)) +
#   geom_point() + xlim(0,10) + ylim(0,10)

tmp<-all[which(all$pJ>=5E-8),c("MarkerName","class_signal")]
tmp<-tmp[!duplicated(tmp),]
dim(tmp)
# [1] 83

table(tmp$class_signal)
#                  new_cojo_supervised                new_cojo_unsupervised 
#                                    1                                   62 
# new_cojo_unsupervised_in_ld_with_old            old|new_cojo_unsupervised 
#                                   17                                    3 

# out of these retain those which are significant in join model in at least one  of the models, ibd, cd, uc:

tmp$retain<-0
for (i in 1:nrow(tmp)) {

  tmp1<-all[which(all$MarkerName==tmp$MarkerName[i]),]
  tmp1<-tmp1[which(tmp1$pJ<=5E-8),]

  if (nrow(tmp1)>0) {
    tmp$retain[i]<-1
  }

  rm(tmp1)
}

 table(tmp$retain)
#  0  1 
# 27 56

 table(tmp$class_signal,tmp$retain)
  #                                       0  1
  # new_cojo_supervised                   1  0
  # new_cojo_unsupervised                21 41
  # new_cojo_unsupervised_in_ld_with_old  4 13
  # old|new_cojo_unsupervised             1  2


all_recap<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_no_tags.tsv.gz",sep=""))
dim(all_recap)
# [1] 664 163

all_recap<-all_recap[which(!all_recap$MarkerName %in% tmp$MarkerName[which(tmp$retain==0)])]
dim(all_recap)
# [1] 637 163


fwrite(all_recap,paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_no_subset.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")

### save the list of variants to run a join model:

table(all_recap$chr)
#  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 
# 80 58 27 26 37 34 54 25 29 32 31 23 12 19 17 41 30 10 28 24  9 15 

write.table(all_recap[,c("MarkerName")],paste(path_gwas,"post_imputation/2022/analysis/conditional_analysis/final_independent_model/final_independent_model_final_join_model_no_subset.snplist",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")


q("no")

###################################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/

pheno=(ibd cd uc)
MEM=2500

for ph in ${pheno[@]}
do
for chr in {1..22}
do
bsub -J"cojo_i" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" \
-o ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_conditioning_final_model_ld_allarrays_stdout \
-e ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_conditioning_final_model_ld_allarrays_stderr \
"/path/to/software/username/gcta-1.94.1-linux-kernel-3-x86_64/./gcta64 \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_subset_included_in_${ph}_analysis \
--chr ${chr} \
--extract ${path_gwas}post_imputation/2022/analysis/conditional_analysis/final_independent_model/final_independent_model_final_join_model_no_subset.snplist \
--cojo-file ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allchr_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format_no_subset.ma \
--cojo-joint \
--out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/final_independent_model/${chr}_${ph}_independent_index_variants_eur_tier_2_ld_allarrays_conditioning_on_final_independent_model_0.1_no_subset_final_join_model"
done
done


for ph in ${pheno[@]}
do
echo ${ph} && for chr in {1..22}
do 
echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_conditioning_final_model_ld_allarrays_stdout | grep -E "Successfully|Exited"
done
done


################################################################

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=2000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
library(ggplot2)

rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"
pheno<-c("ibd","cd","uc")

# Helper: load plink LD files for all autosomes, filter by R2 threshold
load_ld <- function(path_gwas, r2_min = 0) {
  ld_list <- vector("list", 22)
  for (chr in 1:22) {
    ld_tmp <- fread(paste0(path_gwas,
      "post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr", chr,
      "_list_tier_2_unsupervised_conditional_variants_for_ld_with_old_data.ld"))
    if (r2_min > 0) ld_tmp <- ld_tmp[R2 >= r2_min]
    ld_list[[chr]] <- ld_tmp
  }
  rbindlist(ld_list)
}

# Helper: load COJO joint model results across all autosomes and phenotypes
load_cojo <- function(path_gwas, suffix, cols = c("SNP", "bJ", "bJ_se", "pJ", "LD_r")) {
  pheno <- c("ibd", "cd", "uc")
  dat <- NULL
  for (ph in pheno) {
    dat_list <- vector("list", 22)
    for (chr in 1:22) {
      dat_list[[chr]] <- fread(paste0(path_gwas,
        "post_imputation/2022/analysis/conditional_analysis/final_independent_model/",
        chr, "_", ph,
        "_independent_index_variants_eur_tier_2_ld_allarrays_conditioning_on_final_independent_model_0.1_no_subset",
        suffix, ".jma.cojo"))
    }
    dat_tmp <- rbindlist(dat_list)[, ..cols]
    sel_cols <- cols[cols != "SNP"]
    setnames(dat_tmp, sel_cols, paste(sel_cols, ph, sep = "_"))
    print(dim(dat_tmp))
    if (is.null(dat)) {
      dat <- dat_tmp
    } else {
      dat <- merge(dat, dat_tmp, by = "SNP", all = TRUE)
    }
  }
  dat
}

all_recap<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_no_subset.tsv.gz",sep=""))
dim(all_recap)
# [1] 637 163

dat <- load_cojo(path_gwas, suffix = "_final_join_model")

all_recap[which(!all_recap$MarkerName %in% dat$SNP),]
# 0

all_recap<-merge(all_recap,dat,by.x="MarkerName",by.y="SNP")
dim(all_recap)
# [1] 637 175

ibd<-all_recap[,c("MarkerName","P-value_ibd_eur_tier_2","pJ_ibd","class_signal")]
colnames(ibd)<-c("MarkerName","p","pJ","class_signal")
ibd$pheno<-"ibd"
cd<-all_recap[,c("MarkerName","P-value_cd_eur_tier_2","pJ_cd","class_signal")]
colnames(cd)<-c("MarkerName","p","pJ","class_signal")
cd$pheno<-"cd"
uc<-all_recap[,c("MarkerName","P-value_uc_eur_tier_2","pJ_uc","class_signal")]
colnames(uc)<-c("MarkerName","p","pJ","class_signal")
uc$pheno<-"uc"

all<-rbind(ibd,cd,uc)
table(all$pheno)
#  cd ibd  uc 
# 637 637 637

all$p<-as.numeric(all$p)
head(all)

all<-all[which(all$p<=5E-8),]
table(all$pheno)
#  cd ibd  uc 
# 319 451 323 


# chr4:42413003:C:T

# ggplot(all, aes(x=-log10(pJ), y=-log10(p),color=all$class_signal)) +
#   geom_point() + xlim(0,10) + ylim(0,10)

tmp<-all[which(all$pJ>5E-8),c("MarkerName","class_signal")]
tmp<-tmp[!duplicated(tmp),]
dim(tmp)
# [1] 52

table(tmp$class_signal)
          #      new_cojo_unsupervised new_cojo_unsupervised_in_ld_with_old 
          #                         38                                   12 
          #  old|new_cojo_unsupervised 
          #                          2 

# out of these retain those which are significant in join model in at least one  of the models, ibd, cd, uc:

tmp$retain<-0
for (i in 1:nrow(tmp)) {

  tmp1<-all[which(all$MarkerName==tmp$MarkerName[i]),]
  tmp1<-tmp1[which(tmp1$pJ<=5E-8),]

  if (nrow(tmp1)>0) {
    tmp$retain[i]<-1
  }

  rm(tmp1)
}

 table(tmp$retain)
#  0  1 
#  1 52

 table(tmp$class_signal,tmp$retain)
  #                                       0  1
  # new_cojo_unsupervised                 1 37
  # new_cojo_unsupervised_in_ld_with_old  0 12
  # old|new_cojo_unsupervised             0  2

all_recap<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_no_subset.tsv.gz",sep=""))
dim(all_recap)
# [1] 637 163

all_recap<-all_recap[which(!all_recap$MarkerName %in% tmp$MarkerName[which(tmp$retain==0)]),]
dim(all_recap)
# [1] 636 163


table(all_recap$class_signal_final)
# new_cojo_supervised_gw_significant_multiancestry_known_signal 
#                                                             5 # known
#   new_cojo_supervised_gw_significant_multiancestry_new_signal 
#                                                            16  # new
#                            new_cojo_unsupervised_known_signal 
#                                                           235 # known
#                              new_cojo_unsupervised_new_signal 
#                                                           397 # new



## add the conditional resutls:

dat <- load_cojo(path_gwas, suffix = "_final_join_model", cols = c("SNP", "bJ", "bJ_se", "pJ"))

# get max LD friend between variants:

ld <- load_ld(path_gwas, r2_min = 0)

ld<-ld[which(ld$SNP_A %in% dat$SNP & ld$SNP_B %in% dat$SNP),]

dim(dat[which(!dat$SNP %in% ld$SNP_A | !dat$SNP %in% ld$SNP_A),])
# [1] 22 12
dat$SNP[which(!dat$SNP %in% ld$SNP_A | !dat$SNP %in% ld$SNP_A)]
#  [1] "chr1:234599210:C:T"   "chr10:125794381:G:A"  "chr11:128523008:T:TC"
#  [4] "chr12:124548792:T:G"  "chr13:99433020:G:A"   "chr14:106697560:G:T" 
#  [7] "chr15:100948377:T:C"  "chr16:86154861:T:C"   "chr17:83136753:C:G"  
# [10] "chr18:79490197:G:C"   "chr19:58185431:A:G"   "chr2:241816853:C:T"  
# [13] "chr20:64080106:C:T"   "chr21:44095466:T:G"   "chr22:50532837:T:C"  
# [16] "chr3:196688484:G:C"   "chr4:155676868:C:T"   "chr5:177372991:T:C"  
# [19] "chr6:170103419:C:G"   "chr7:150971426:T:C"   "chr8:144103704:G:A"  
# [22] "chr9:136920601:G:A"  

dat$best_tag_variant<-NA
dat$best_tag_variant_r2<-NA


for (i in 1:nrow(dat)) {

  tmp<-ld[which(ld$SNP_A==dat$SNP[i] | ld$SNP_B==dat$SNP[i]),]
  tmp<-tmp[which(tmp$R2==max(tmp$R2)),]

  dat$best_tag_variant[i]<-c(tmp$SNP_A,tmp$SNP_B)[!c(tmp$SNP_A,tmp$SNP_B) %in% dat$SNP[i]]
  dat$best_tag_variant_r2[i]<-tmp$R2

}

summary(dat$best_tag_variant_r2)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 2.311e-05 1.571e-04 6.498e-04 9.366e-03 6.741e-03 2.384e-01  

dat[which(dat$best_tag_variant_r2>0.1),]
#                    SNP     bJ_ibd  bJ_se_ibd      pJ_ibd      bJ_cd  bJ_se_cd
# 303 chr19:33253966:T:A  0.0636130 0.00760967 6.29859e-17  0.0782363 0.0108303
# 304 chr19:33263138:G:A  0.0861012 0.00895274 6.75828e-22  0.0920288 0.0127287
# 393 chr20:63639284:T:G -0.0715871 0.00911230 3.96277e-15 -0.0577149 0.0125838
# 396 chr20:63697746:G:A  0.0821402 0.00856142 8.45431e-22  0.0957216 0.0124412
#           pJ_cd      bJ_uc   bJ_se_uc       pJ_uc   best_tag_variant
# 303 5.05465e-13  0.0585851 0.00937551 4.13780e-10 chr19:33263138:G:A
# 304 4.82839e-13  0.0864925 0.01104980 4.97659e-15 chr19:33253966:T:A
# 393 4.50888e-06 -0.0872593 0.01117260 5.71359e-15 chr20:63697746:G:A
# 396 1.42666e-14  0.0784149 0.01051610 8.87600e-14 chr20:63639284:T:G
#     best_tag_variant_r2
# 303            0.228762
# 304            0.228762
# 393            0.238448
# 396            0.238448



tmp<-dat[which(dat$pJ_ibd>5E-8 & dat$pJ_cd>5E-8 & dat$pJ_uc>5E-8),c("SNP","pJ_ibd","pJ_cd","pJ_uc","best_tag_variant_r2")]
dim(tmp)
# [1] 27  5

tmp[order(tmp$SNP),]
#                      SNP      pJ_ibd       pJ_cd       pJ_uc
# 12     chr1:15949011:T:C 4.04591e-07 4.64865e-03 1.94707e-06
# 55  chr1:46217400:ATGT:A 7.55688e-08 5.94847e-06 5.88846e-05
# 117  chr11:118195867:T:C 2.19618e-07 7.85791e-08 1.13342e-03
# 147  chr12:124548792:T:G 3.47085e-03 2.11551e-05 7.15152e-01
# 149   chr12:30654261:C:A 5.25007e-08 1.85982e-03 3.58163e-07
# 154   chr12:56076060:G:A 2.03398e-03 1.65631e-02 2.60514e-03
# 181   chr14:24146226:T:C 5.83790e-03 1.00982e-05 2.91450e-01
# 183   chr14:39135131:G:A 8.21045e-08 2.19788e-04 1.91317e-06
# 194  chr15:100948377:T:C 2.56859e-03 3.67520e-07 8.43412e-01
# 216  chr16:27384341:C:CT 8.56845e-05 7.43152e-01 7.48829e-06
# 225   chr16:50588216:A:G 1.69834e-03 7.62063e-08 8.84952e-01
# 239   chr16:75463237:C:T 1.61174e-05 1.98367e-02 1.01607e-05
# 318   chr2:100146324:T:G 4.49780e-03 3.13011e-01 1.41155e-07
# 353   chr2:230982974:G:A 6.75747e-07 1.17952e-02 9.25351e-08
# 386   chr20:46120612:T:C 7.65669e-03 9.51449e-07 6.00478e-01
# 397   chr20:63826135:G:C 8.56388e-04 1.45696e-02 3.08907e-03
# 403   chr21:35350264:T:C 2.40203e-06 1.50806e-06 2.04873e-03
# 409   chr22:22909672:C:T 8.68376e-06 1.79395e-02 2.03228e-06
# 465    chr4:42413003:C:T 2.19465e-02 1.92145e-01 2.99649e-04
# 469    chr4:54212027:A:T 3.82528e-02 3.36221e-07 6.05474e-01
# 470     chr4:7039650:A:G 7.03325e-04 6.34889e-08 4.79968e-01
# 485   chr5:151097526:C:T 7.05263e-07 2.03084e-02 5.47036e-08
# 509  chr5:88820075:TTA:T 2.12594e-07 1.95641e-03 4.89515e-06
# 542    chr6:87705790:C:T 2.91399e-06 3.02301e-06 4.20525e-04
# 575    chr7:48997332:T:A 2.82100e-03 1.72030e-06 5.35835e-01
# 596      chr7:931488:G:A 4.77209e-05 9.92104e-08 6.35706e-02
# 640    chr9:29702044:G:A 7.51087e-08 6.19060e-05 1.66139e-06
# 648    chr9:91166134:C:T 3.41506e-05 1.09655e-01 3.30313e-07

# which ones were brought by Multiancestry analyses:
table(all_recap$class_signal_final)
# new_cojo_supervised_gw_significant_multiancestry_known_signal  
#                                                             5 - genome wide only in multiancestry
#   new_cojo_supervised_gw_significant_multiancestry_new_signal 
#                                                            17 - genome wide only in multiancestry
#                            new_cojo_unsupervised_known_signal 
#                                                           234 
#                              new_cojo_unsupervised_new_signal 
#                                                           393 

17+5
# [1] 22

all_recap$'P-value_ibd_eur_tier_2_eas_sas'<-as.numeric(all_recap$'P-value_ibd_eur_tier_2_eas_sas')
all_recap$'P-value_cd_eur_tier_2_eas_sas'<-as.numeric(all_recap$'P-value_cd_eur_tier_2_eas_sas')
all_recap$'P-value_uc_eur_tier_2_eas_sas'<-as.numeric(all_recap$'P-value_uc_eur_tier_2_eas_sas')

# exclude them all:
tmp[which(tmp$SNP %in% all_recap$MarkerName[which(all_recap$'P-value_ibd_eur_tier_2_eas_sas'>5E-8 & all_recap$'P-value_cd_eur_tier_2_eas_sas'>5E-8 & all_recap$'P-value_uc_eur_tier_2_eas_sas'>5E-8)]),]
#                      SNP      pJ_ibd       pJ_cd       pJ_uc
# 55  chr1:46217400:ATGT:A 7.55688e-08 5.94847e-06 5.88846e-05
# 117  chr11:118195867:T:C 2.19618e-07 7.85791e-08 1.13342e-03
# 397   chr20:63826135:G:C 8.56388e-04 1.45696e-02 3.08907e-03
# 465    chr4:42413003:C:T 2.19465e-02 1.92145e-01 2.99649e-04
# 485   chr5:151097526:C:T 7.05263e-07 2.03084e-02 5.47036e-08
#     best_tag_variant_r2
# 55          0.000409506
# 117         0.001711690
# 397         0.027757500
# 465         0.000292806
# 485         0.002041210

# retain the remaining - double check with the team
tmp[which(!tmp$SNP %in% all_recap$MarkerName[which(all_recap$'P-value_ibd_eur_tier_2_eas_sas'>5E-8 & all_recap$'P-value_cd_eur_tier_2_eas_sas'>5E-8 & all_recap$'P-value_uc_eur_tier_2_eas_sas'>5E-8)]),]
#                     SNP      pJ_ibd       pJ_cd       pJ_uc best_tag_variant_r2
# 12    chr1:15949011:T:C 4.04591e-07 4.64865e-03 1.94707e-06         3.75368e-04
# 147 chr12:124548792:T:G 3.47085e-03 2.11551e-05 7.15152e-01         1.33722e-04
# 149  chr12:30654261:C:A 5.25007e-08 1.85982e-03 3.58163e-07         1.59699e-04
# 154  chr12:56076060:G:A 2.03398e-03 1.65631e-02 2.60514e-03         1.06157e-04
# 181  chr14:24146226:T:C 5.83790e-03 1.00982e-05 2.91450e-01         3.83167e-05
# 183  chr14:39135131:G:A 8.21045e-08 2.19788e-04 1.91317e-06         4.93349e-04
# 194 chr15:100948377:T:C 2.56859e-03 3.67520e-07 8.43412e-01         2.59266e-05
# 216 chr16:27384341:C:CT 8.56845e-05 7.43152e-01 7.48829e-06         1.54965e-04
# 225  chr16:50588216:A:G 1.69834e-03 7.62063e-08 8.84952e-01         3.00030e-02
# 239  chr16:75463237:C:T 1.61174e-05 1.98367e-02 1.01607e-05         1.05405e-04
# 318  chr2:100146324:T:G 4.49780e-03 3.13011e-01 1.41155e-07         1.11875e-03
# 353  chr2:230982974:G:A 6.75747e-07 1.17952e-02 9.25351e-08         1.90575e-03
# 386  chr20:46120612:T:C 7.65669e-03 9.51449e-07 6.00478e-01         1.40137e-04
# 403  chr21:35350264:T:C 2.40203e-06 1.50806e-06 2.04873e-03         1.48533e-04
# 409  chr22:22909672:C:T 8.68376e-06 1.79395e-02 2.03228e-06         1.63124e-04
# 469   chr4:54212027:A:T 3.82528e-02 3.36221e-07 6.05474e-01         1.04314e-04
# 470    chr4:7039650:A:G 7.03325e-04 6.34889e-08 4.79968e-01         9.57797e-05
# 509 chr5:88820075:TTA:T 2.12594e-07 1.95641e-03 4.89515e-06         1.63096e-04
# 542   chr6:87705790:C:T 2.91399e-06 3.02301e-06 4.20525e-04         1.11050e-04
# 572   chr7:48997332:T:A 2.82100e-03 1.72030e-06 5.35835e-01         4.31107e-04
# 628   chr9:29702044:G:A 7.51087e-08 6.19060e-05 1.66139e-06         8.69377e-05
# 636   chr9:91166134:C:T 3.41506e-05 1.09655e-01 3.30313e-07         2.17638e-04


dim(all_recap)
# [1] 636 163

all_recap<-all_recap[which(!all_recap$MarkerName %in% tmp$SNP[which(tmp$SNP %in% all_recap$MarkerName[which(all_recap$'P-value_ibd_eur_tier_2_eas_sas'>5E-8 & all_recap$'P-value_cd_eur_tier_2_eas_sas'>5E-8 & all_recap$'P-value_uc_eur_tier_2_eas_sas'>5E-8)])]),]

dim(all_recap)
# [1] 631 163


# list of variants with het-pval lower than 0.05/631 - and
pval_threshold<-0.05/nrow(all_recap)

all_recap$HetPVal_uc_eur_tier_2<-as.numeric(all_recap$HetPVal_uc_eur_tier_2)
all_recap$HetPVal_cd_eur_tier_2<-as.numeric(all_recap$HetPVal_cd_eur_tier_2)
all_recap$HetPVal_ibd_eur_tier_2<-as.numeric(all_recap$HetPVal_ibd_eur_tier_2)


tmp<-all_recap[which(all_recap$HetPVal_uc_eur_tier_2<pval_threshold | all_recap$HetPVal_cd_eur_tier_2<pval_threshold | all_recap$HetPVal_ibd_eur_tier_2<pval_threshold),]
dim(tmp)
# [1]  29 163

# after regional plot inspection, exclude those where the only one study was contributing to the summary statistics:
vec<-c("chr11:85454082:A:G","chr2:10840215:C:T","chr4:42459885:G:A","chr17:66833937:C:T","chr9:29702044:G:A","chr10:49517763:C:T","chr8:15592628:A:C","chr2:18822378:C:A","chr1:53476954:G:A","chr5:77636984:G:A","chr4:42689519:A:G","chr12:14388669:C:T")
dim(tmp[which(tmp$MarkerName %in% vec),])
# [1]  12 166
length(vec)
# [1] 12

all_recap<-all_recap[which(!all_recap$MarkerName %in% vec),]
dim(all_recap)
# [1] 619 163

### save the list of variants to run a join model:

table(all_recap$chr)
#  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 
# 77 55 27 21 35 34 40 24 27 30 29 22 12 17 16 39 29 10 28 23  9 15 


fwrite(all_recap,paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_no_subset_final.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")


write.table(all_recap[,c("MarkerName")],paste(path_gwas,"post_imputation/2022/analysis/conditional_analysis/final_independent_model/final_independent_model_final_join_model_no_subset_final.snplist",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")

q("no")

################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/

pheno=(ibd cd uc)
MEM=2500

for ph in ${pheno[@]}
do
for chr in {1..22}
do
bsub -J"cojo_i" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" \
-o ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_conditioning_final_model_ld_allarrays_stdout \
-e ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_conditioning_final_model_ld_allarrays_stderr \
"/path/to/software/username/gcta-1.94.1-linux-kernel-3-x86_64/./gcta64 \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_subset_included_in_${ph}_analysis \
--chr ${chr} \
--extract ${path_gwas}post_imputation/2022/analysis/conditional_analysis/final_independent_model/final_independent_model_final_join_model_no_subset_final.snplist \
--cojo-file ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allchr_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format_no_subset.ma \
--cojo-joint \
--out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/final_independent_model/${chr}_${ph}_independent_index_variants_eur_tier_2_ld_allarrays_conditioning_on_final_independent_model_0.1_no_subset_final_join_model_final"
done
done


for ph in ${pheno[@]}
do
echo ${ph} && for chr in {1..22}
do 
echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_conditioning_final_model_ld_allarrays_stdout | grep -E "Successfully|Exited"
done
done


################################################################

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=15000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
library(ggplot2)

rm(list=ls())

# Helper: load plink LD files for all autosomes, filter by R2 threshold
load_ld <- function(path_gwas, r2_min = 0) {
  ld_list <- vector("list", 22)
  for (chr in 1:22) {
    ld_tmp <- fread(paste0(path_gwas,
      "post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr", chr,
      "_list_tier_2_unsupervised_conditional_variants_for_ld_with_old_data.ld"))
    if (r2_min > 0) ld_tmp <- ld_tmp[R2 >= r2_min]
    ld_list[[chr]] <- ld_tmp
  }
  rbindlist(ld_list)
}

# Helper: load COJO joint model results across all autosomes and phenotypes
# type: "final"       -> final_independent_model dir, conditioning_on_final_independent_model_0.1_no_subset stem
#       "unsupervised"-> eur dir, no_conditioning stem
load_cojo <- function(path_gwas, suffix, cols = c("SNP", "bJ", "bJ_se", "pJ", "LD_r"),
                      type = "final") {
  pheno <- c("ibd", "cd", "uc")
  dat <- NULL
  for (ph in pheno) {
    dat_list <- vector("list", 22)
    for (chr in 1:22) {
      base <- switch(type,
        final = paste0(path_gwas,
          "post_imputation/2022/analysis/conditional_analysis/final_independent_model/",
          chr, "_", ph,
          "_independent_index_variants_eur_tier_2_ld_allarrays_conditioning_on_final_independent_model_0.1_no_subset"),
        unsupervised = paste0(path_gwas,
          "post_imputation/2022/analysis/conditional_analysis/eur/",
          chr, "_", ph,
          "_independent_index_variants_eur_tier_2_ld_allarrays_no_conditioning"),
        stop("Unknown type: ", type)
      )
      dat_list[[chr]] <- fread(paste0(base, suffix, ".jma.cojo"))
    }
    dat_tmp <- rbindlist(dat_list)[, ..cols]
    sel_cols <- cols[cols != "SNP"]
    setnames(dat_tmp, sel_cols, paste(sel_cols, ph, sep = "_"))
    print(dim(dat_tmp))
    if (is.null(dat)) {
      dat <- dat_tmp
    } else {
      dat <- merge(dat, dat_tmp, by = "SNP", all = TRUE)
    }
  }
  dat
}

# Helper: load COJO supervised conditional results (conditioning_on_known_ibd_index, .cma.cojo),
# subsetting on the fly to specific index variants to avoid loading all variants.
# index_variants: a character vector of SNP IDs to retain (subset applied per chromosome as loaded).
load_cojo_supervised <- function(path_gwas, suffix = "", cols = c("SNP", "bC", "bC_se", "pC"),
                                 index_variants = NULL) {
  pheno <- c("ibd", "cd", "uc")
  dat <- NULL
  for (ph in pheno) {
    dat_list <- vector("list", 22)
    for (chr in 1:22) {
      tmp <- fread(paste0(path_gwas,
        "post_imputation/2022/analysis/conditional_analysis/eur/",
        chr, "_", ph,
        "_independent_index_variants_eur_tier_2_ld_allarrays_conditioning_on_known_ibd_index",
        suffix, ".cma.cojo"))
      if (!is.null(index_variants) && length(index_variants) > 0) {
        tmp <- tmp[SNP %in% index_variants]
      }
      dat_list[[chr]] <- tmp
    }
    dat_tmp <- rbindlist(dat_list)[, ..cols]
    sel_cols <- cols[cols != "SNP"]
    setnames(dat_tmp, sel_cols, paste(sel_cols, ph, sep = "_"))
    print(dim(dat_tmp))
    if (is.null(dat)) {
      dat <- dat_tmp
    } else {
      dat <- merge(dat, dat_tmp, by = "SNP", all = TRUE)
    }
  }
  dat
}


path_gwas<-"/path/to/ibdgwas/IIBDGC/"
pheno<-c("ibd","cd","uc")

all_recap<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_no_subset_final.tsv.gz",sep=""))
dim(all_recap)
# [1] 619 163

all_recap$'P-value_ibd_eur_tier_2_eas_sas'<-as.numeric(all_recap$'P-value_ibd_eur_tier_2_eas_sas')
all_recap$'P-value_cd_eur_tier_2_eas_sas'<-as.numeric(all_recap$'P-value_cd_eur_tier_2_eas_sas')
all_recap$'P-value_uc_eur_tier_2_eas_sas'<-as.numeric(all_recap$'P-value_uc_eur_tier_2_eas_sas')

dat <- load_cojo(path_gwas, suffix = "_final_join_model_final", cols = c("SNP", "bJ", "bJ_se", "pJ"))

all<-merge(all_recap,dat,by.x="MarkerName",by.y="SNP")

all[which((all$'P-value_ibd_eur_tier_2_eas_sas'>5E-8 & all$'P-value_cd_eur_tier_2_eas_sas'>5E-8 & all$'P-value_uc_eur_tier_2_eas_sas'>5E-8) & (all$pJ_ibd>5E-8 & all$pJ_cd>5E-8 & all$pJ_uc>5E-8)),c("MarkerName")]
# 0

dim(all[which((all$pJ_ibd>5E-8 & all$pJ_cd>5E-8 & all$pJ_uc>5E-8)),c("MarkerName")])
# [1] 21  1

all[which((all$'P-value_ibd_eur_tier_2_eas_sas'<5E-8 | all$'P-value_cd_eur_tier_2_eas_sas'<5E-8 | all$'P-value_uc_eur_tier_2_eas_sas'<5E-8) & (all$pJ_ibd>5E-8 & all$pJ_cd>5E-8 & all$pJ_uc>5E-8)),c("MarkerName",
"pJ_ibd","P-value_ibd_eur_tier_2","P-value_ibd_eur_tier_2_eas_sas",
"pJ_cd","P-value_cd_eur_tier_2","P-value_cd_eur_tier_2_eas_sas",
"pJ_uc","P-value_uc_eur_tier_2","P-value_uc_eur_tier_2_eas_sas")]
#              MarkerName      pJ_ibd P-value_ibd_eur_tier_2
#                  <char>       <num>                  <num>
#  1: chr12:124548792:T:G 3.47085e-03              3.103e-03
#  2:  chr12:30654261:C:A 5.25007e-08              1.088e-07
#  3:  chr12:56076060:G:A 2.03398e-03              1.812e-03
#  4:  chr14:24146226:T:C 5.83790e-03              5.607e-03
#  5:  chr14:39135131:G:A 8.21045e-08              1.113e-07
#  6: chr15:100948377:T:C 2.56859e-03              2.586e-03
#  7: chr16:27384341:C:CT 8.56845e-05              6.842e-05
#  8:  chr16:75463237:C:T 1.61174e-05              1.180e-05
#  9:   chr1:15949011:T:C 4.04591e-07              1.888e-07
# 10:  chr20:46120612:T:C 7.65666e-03              4.833e-03
# 11:  chr21:35350264:T:C 2.40203e-06              1.564e-06
# 12:  chr22:22909672:C:T 8.68376e-06              6.747e-06
# 13:  chr2:100146324:T:G 4.49780e-03              4.328e-03
# 14:  chr2:230982974:G:A 6.75747e-07              7.221e-08
# 15:   chr4:54212027:A:T 3.82574e-02              3.671e-02
# 16:    chr4:7039650:A:G 7.03325e-04              6.845e-04
# 17: chr5:88820075:TTA:T 2.12594e-07              1.294e-07
# 18:   chr6:87705790:C:T 2.91399e-06              3.347e-06
# 19:   chr7:48997332:T:A 2.82100e-03              4.622e-03
# 20:   chr9:29702044:G:A 7.51087e-08              7.209e-08
# 21:   chr9:91166134:C:T 3.41506e-05              4.804e-05
#              MarkerName      pJ_ibd P-value_ibd_eur_tier_2
#     P-value_ibd_eur_tier_2_eas_sas       pJ_cd P-value_cd_eur_tier_2
#                              <num>       <num>                 <num>
#  1:                      1.218e-04 2.11551e-05             2.009e-05
#  2:                      2.319e-08 1.85982e-03             3.319e-03
#  3:                      2.682e-08 1.65631e-02             1.603e-02
#  4:                      1.680e-04 1.00982e-05             9.261e-06
#  5:                      1.388e-08 2.19788e-04             2.453e-04
#  6:                      1.327e-03 3.67520e-07             3.377e-07
#  7:                      3.313e-08 7.43152e-01             7.020e-01
#  8:                      3.422e-09 1.98367e-02             2.047e-02
#  9:                      3.228e-09 4.64865e-03             4.142e-03
# 10:                      5.301e-06 9.51455e-07             6.546e-07
# 11:                      2.085e-08 1.50806e-06             1.230e-06
# 12:                      2.003e-07 1.79395e-02             2.091e-02
# 13:                      2.772e-02 3.13011e-01             2.718e-01
# 14:                      1.424e-08 1.17952e-02             2.462e-03
# 15:                      1.047e-02 3.36280e-07             2.327e-07
# 16:                      1.147e-04 6.34889e-08             5.746e-08
# 17:                      9.292e-11 1.95641e-03             1.331e-03
# 18:                      2.851e-08 3.02301e-06             3.674e-06
# 19:                      7.841e-04 1.72030e-06             2.327e-06
# 20:                      3.968e-09 6.19060e-05             6.373e-05
# 21:                      7.189e-06 1.09655e-01             1.262e-01
#     P-value_ibd_eur_tier_2_eas_sas       pJ_cd P-value_cd_eur_tier_2
#     P-value_cd_eur_tier_2_eas_sas       pJ_uc P-value_uc_eur_tier_2
#                             <num>       <num>                 <num>
#  1:                     1.375e-08 7.15152e-01             7.256e-01
#  2:                     2.574e-03 3.58163e-07             4.632e-07
#  3:                     8.668e-08 2.60514e-03             2.268e-03
#  4:                     2.318e-08 2.91450e-01             2.920e-01
#  5:                     4.699e-05 1.91317e-06             1.921e-06
#  6:                     8.323e-10 8.43412e-01             8.480e-01
#  7:                     1.425e-01 7.48829e-06             6.094e-06
#  8:                     7.341e-04 1.01607e-05             1.071e-05
#  9:                     6.245e-04 1.94707e-06             5.030e-07
# 10:                     1.157e-11 6.00483e-01             7.734e-01
# 11:                     2.948e-07 2.04873e-03             2.294e-03
# 12:                     1.670e-02 2.03228e-06             1.812e-06
# 13:                     1.433e-02 1.41155e-07             9.889e-08
# 14:                     5.083e-04 9.25351e-08             7.353e-08
# 15:                     1.069e-09 6.05326e-01             6.070e-01
# 16:                     2.731e-09 4.79968e-01             4.777e-01
# 17:                     5.006e-05 4.89515e-06             3.885e-06
# 18:                     5.763e-09 4.20525e-04             5.094e-04
# 19:                     3.606e-08 5.35835e-01             6.148e-01
# 20:                     3.646e-05 1.66139e-06             1.452e-06
# 21:                     2.912e-01 3.30313e-07             3.453e-07
#     P-value_cd_eur_tier_2_eas_sas       pJ_uc P-value_uc_eur_tier_2
#     P-value_uc_eur_tier_2_eas_sas
#                             <num>
#  1:                     6.091e-01
#  2:                     7.413e-08
#  3:                     4.368e-05
#  4:                     1.467e-01
#  5:                     6.583e-07
#  6:                     7.090e-01
#  7:                     1.898e-10
#  8:                     2.986e-09
#  9:                     1.099e-08
# 10:                     4.221e-01
# 11:                     9.271e-05
# 12:                     8.085e-09
# 13:                     8.098e-09
# 14:                     3.769e-08
# 15:                     5.484e-01
# 16:                     3.779e-01
# 17:                     8.163e-09
# 18:                     7.120e-05
# 19:                     5.213e-01
# 20:                     8.058e-08
# 21:                     2.796e-10
#     P-value_uc_eur_tier_2_eas_sas

## double check those were variants captured by cojo supervised:

all$conditional_nogws<-0
all$conditional_nogws[which((all$'P-value_ibd_eur_tier_2_eas_sas'<5E-8 | all$'P-value_cd_eur_tier_2_eas_sas'<5E-8 | all$'P-value_uc_eur_tier_2_eas_sas'<5E-8) & (all$pJ_ibd>5E-8 & all$pJ_cd>5E-8 & all$pJ_uc>5E-8))]<-1

table(all$conditional_nogws)
#   0   1 
# 598  21 

recaptured<-all[which(all$conditional_nogws==1),]

## add the conditional resutls, supervised
dat <- load_cojo_supervised(path_gwas, suffix="", cols = c("SNP", "bC", "bC_se", "pC"),
                            index_variants = recaptured$MarkerName)


tmp<-merge(tmp,dat,by.x="MarkerName",by.y="SNP",all.x=T)




## add the conditional resutls, unsupervised
dat <- load_cojo(path_gwas, suffix="", type="unsupervised", cols = c("SNP", "bJ", "bJ_se", "pJ"))
colnames(dat)[2:ncol(dat)]<-paste0(colnames(dat)[2:ncol(dat)],"_supervised")

tmp<-merge(tmp,dat,by.x="MarkerName",by.y="SNP",all.x=T)


## manually inspect some:
vec<-c("chr1:50950322:G:A","chr1:53476954:G:A")
all[which(all$MarkerName %in% vec),c("MarkerName","pJ_ibd","pJ_cd","pJ_uc")]
#           MarkerName      pJ_ibd       pJ_cd       pJ_uc
#               <char>       <num>       <num>       <num>
# 1: chr1:50950322:G:A 3.24883e-09 1.27289e-05 2.24804e-05
# 2: chr1:53476954:G:A 3.36587e-09 2.76861e-11 3.28067e-06

vec<-c("chr1:67240217:G:A","chr1:67240275:G:A")
all[which(all$MarkerName %in% vec),c("MarkerName","pJ_ibd","pJ_cd","pJ_uc")]
#           MarkerName       pJ_ibd        pJ_cd       pJ_uc
#               <char>        <num>        <num>       <num>
# 1: chr1:67240217:G:A  5.27055e-08  9.17962e-09 2.50891e-04
# 2: chr1:67240275:G:A 7.20404e-203 9.13123e-174 3.08493e-95

vec<-c("chr7:5266889:G:T","chr7:5285792:A:T")
all[which(all$MarkerName %in% vec),c("MarkerName","pJ_ibd","pJ_cd","pJ_uc")]
#          MarkerName      pJ_ibd       pJ_cd       pJ_uc
#              <char>       <num>       <num>       <num>
# 1: chr7:5266889:G:T 4.81504e-09 1.77842e-02 5.33997e-08
# 2: chr7:5285792:A:T 4.72431e-18 2.44879e-05 4.09263e-11


vec<-c("chr7:5333452:G:A","chr7:5433979:G:C")
all[which(all$MarkerName %in% vec),c("MarkerName","pJ_ibd","pJ_cd","pJ_uc")]
#          MarkerName      pJ_ibd       pJ_cd       pJ_uc
#              <char>       <num>       <num>       <num>
# 1: chr7:5333452:G:A 1.46033e-08 1.38056e-02 8.12691e-08
# 2: chr7:5433979:G:C 1.45666e-39 1.48681e-13 2.27495e-25


vec<-c("chr1:67152482:A:T","chr1:67240275:G:A")
all[which(all$MarkerName %in% vec),c("MarkerName","pJ_ibd","pJ_cd","pJ_uc")]
#           MarkerName       pJ_ibd        pJ_cd       pJ_uc
#               <char>        <num>        <num>       <num>
# 1: chr1:67240275:G:A 7.20404e-203 9.13123e-174 3.08493e-95


dim(all_recap)
# [1] 619 163

dim(all)
# [1] 619 172

table(all$class_signal_final)
# new_cojo_supervised_gw_significant_multiancestry_known_signal 
#                                                             5 
#   new_cojo_supervised_gw_significant_multiancestry_new_signal 
#                                                            16 
#                            new_cojo_unsupervised_known_signal 
#                                                           234 
#                              new_cojo_unsupervised_new_signal 
#                                                           377 

# still need to add the exonic annotation to define the final classification of known/new signal

fwrite(all,paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_no_subset_final_2.tsv.gz",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")


# compare with previous list:
old<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model.tsv.gz",sep=""))
dim(old)

dim(all[which(all$MarkerName %in% old$MarkerName),c("MarkerName")])
# [1] 630   1
dim(all[which(!all$MarkerName %in% old$MarkerName),c("MarkerName")])
# [1] 1 1

all[which(!all$MarkerName %in% old$MarkerName),c("MarkerName","updated_region")]
#            MarkerName       updated_region
#                <char>               <char>
# 1: chr10:45014562:G:A 10_44514562_45514563


old[which(old$updated_region %in% c("10_44514562_45514563")),"MarkerName"]
# region excluded before

################################################################




################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/

pheno=(ibd cd uc)
MEM=55000

for ph in ${pheno[@]}
do
for chr in {1..22}
do
bsub -J"cojo_i" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" \
-o ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_conditioning_final_model_ld_allarrays_stdout \
-e ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_conditioning_final_model_ld_allarrays_stderr \
"/path/to/software/username/gcta-1.94.1-linux-kernel-3-x86_64/./gcta64 \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_subset_included_in_${ph}_analysis \
--chr ${chr} \
--cojo-cond ${path_gwas}post_imputation/2022/analysis/conditional_analysis/final_independent_model/final_independent_model_final_join_model_no_subset_final.snplist \
--cojo-file ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allchr_${ph}_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_cojo_format_no_subset.ma \
--out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/final_independent_model/${chr}_${ph}_eur_tier_2_conditioning_on_final_independent_model_619"
done
done