# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
############################################################################################################################
############################################################################################################################

# different versions of resutls stored here:
# /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/exome/

# latest:

## PREVIOUSLY USED: 
# STb_Aug16_Aug27_scores_update_region.csv
# IBDexome_ASHG2025_shared_tier1_tier2 
# wes_2023_S3.csv
# sequencing_s4_81_20260424.csv


# 1.- ADD THE NEW wes RESULTS FROM Ruifei:

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=2000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
rm(list=ls())

path_gwas <- "/path/to/ibdgwas/IIBDGC/"
path_final <- paste0(path_gwas, "post_imputation/2022/analysis/final_tables/")
path_ld    <- paste0(path_gwas, "post_imputation/2022/analysis/conditional_analysis/eur/")
ld_prefix  <- "allarrays_chr%d_list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_for_ld_with_exome"
out_file   <- paste0(path_final, "list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants.tsv.gz")

all <- fread(paste0(path_final, "list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_pheno.tsv.gz"))
all <- as.data.frame(all)
dim(all)
# [1]  619 179


# exo<-fread("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/exome/IBDexome_ASHG2025_shared_tier1_tier2.txt")
# exo<-as.data.frame(exo)
# dim(exo)
# # 109
# class(exo$accept)

# var_tier2<-c("chr1:6469122:TTCC:T","chr1:6471529:C:T","chr5:172770787:C:T","chr5:172769749:A:G","chr6:111591867:G:A","chr6:111592059:C:T","chr7:99434906:G:A","chr7:99563862:G:A","chr8:81753241:TGTGA:T","chr8:81758536:G:A","chr22:31139886:G:C","chr22:31137810:G:C","chr5:132340627:C:T","chr2:241816853:C:T","chr1:2556714:A:G","chr1:155133751:A:T","chr1:160882036:A:T","chr16:28502082:A:G","chr12:40346884:A:G","chr11:61008737:C:T","chr20:63569223:G:A","chr19:1127616:G:C","chr19:1127721:G:A","chr10:73913343:T:C","chr4:3448188:A:G","chr2:97725153:C:T","chr15:67165360:A:G","chr2:213147743:T:C","chr6:159041392:C:T","chr5:40959520:C:A","chr9:5126343:G:A")
# exo$tier_exome<-"tier_1"
# exo$tier_exome[which(exo$ID %in% var_tier2)]<-"tier_2"

# table(exo$tier_exome)
# # tier_1  tier2 
# #     78     31 

# exo<-exo[which( (exo$tier_exome=="tier_1") | (exo$tier_exome=="tier_2" & exo$overall_ld_decision=="pass" & exo$accept>0.5 )),]
# dim(exo)
# # [1] 90 31

# table(exo$tier_exome)
# # tier_1 tier_2 
# #     78     12

exo <- fread("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/exome/sequencing_s4_83_20260521.csv")
exo <- as.data.frame(exo)
exo <- exo[exo$chr != "", 1:16]
exo <- exo[!is.na(exo$chr), ]
dim(exo)
# [1] 81 16

dim(exo[which(exo$ID %in% all$MarkerName),])
# [1] 36 16

exo[which(exo$ID %in% all$MarkerName),"Gene_symbol"]
#  [1] "IL23R"    "IL23R"    "IL23R"    "PTPN22"   "SLAMF8"   "FCGR2A"  
#  [7] "ADAM17"   "GCKR"     "ARHGAP25" "IFIH1"    "SLC39A8"  "DUSP1"   
# [13] "RCC1L"    "DOK2"     "SHARPIN"  "CARD9"    "NXPE1"    "IL10RA"  
# [19] "PTGER2"   "ADCY7"    "NOD2"     "NOD2"     "NOD2"     "NOD2"    
# [25] "NOD2"     "CARMIL2"  "PLCG2"    "CCR7"     "TYK2"     "TYK2"    
# [31] "TYK2"     "FUT2"     "MFNG"     "TSPO" 

# TRAF3IP2 JAK2 lost from prvious version - wes_2023_S3.csv

#  [1] "IL23R"    "IL23R"    "IL23R"    "PTPN22"   "SLAMF8"   "FCGR2A"  
#  [7] "ADAM17"   "GCKR"     "ARHGAP25" "IFIH1"    "SLC39A8"  "DUSP1"   
# [13] "RCC1L"    "DOK2"     "SHARPIN"  "CARD9"    "NXPE1"    "IL10RA"  
# [19] "PTGER2"   "ADCY7"    "NOD2"     "NOD2"     "NOD2"     "NOD2"    
# [25] "NOD2"     "CARMIL2"  "PLCG2"    "CCR7"     "TYK2"     "TYK2"    
# [31] "TYK2"     "FUT2"     "MFNG"     "TSPO"     "TRAF3IP2" "NOD2"    

# TRAF3IP2 recovered in this last version

exo[which(!exo$ID %in% all$MarkerName),"Gene_symbol"]
#  [1] "PLEKHG5"   "MAP3K6"    "PTAFR"     "ATG4C"     "BTBD8"     "IFIH1"    
#  [7] "IFIH1"     "GPR35"     "MST1"      "MST1"      "HGFAC"     "HGFAC"    
# [13] "MAN2B2"    "ORC3"      "TRAF3IP2"  "ZC3H12D"   "PTCD1"     "MUC12"    
# [19] "CFTR"      "CFTR"      "DOK2"      "CHMP4C"    "RIPK2"     "SPATA31F1"
# [25] "CARD9"     "VPS37C"    "RELA"      "ARHGDIB"   "NFE2"      "CMKLR1"   
# [31] "LACC1"     "SEL1L"     "GPR65"     "NOD2"      "NOD2"      "NOD2"     
# [37] "NOD2"      "NOD2"      "NOD2"      "NOD2"      "ELMO3"     "PLCG2"    
# [43] "ARHGAP45"  "ZNF101"    "DMWD"      "TNFRSF6B"  "TNFRSF6B"  "PLA2G3"   
# [49] "CSF2RB"    "PLEKHG5"   "DUSP1"     "ZNF655"    "CHMP4C"    "PLA2G3"   
# [55] "NOD2" 


# BORROW THE STATUS (NOVEL, KNOWN AS REPORTED BY EXOME PAPER BEFORE)
exo_tmp <- fread("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/exome/STb_Aug16_Aug27_scores_update_region.csv")
exo_tmp$variant_status

exo <- merge(exo, exo_tmp[, c("ID", "variant_status")], by = "ID", all.x = TRUE)
table(exo$variant_status,useNA="ifany")
  #  Novel_variant Reported_variant 
  #              36              55 


list_ids<-c(all$MarkerName,exo$ID)
length(list_ids)
# [1] 710
list_ids<-list_ids[!duplicated(list_ids)]
length(list_ids)
# [1] 674

fwrite(data.table(list_ids), paste0(path_gwas, "post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_for_ld_with_exome_data.tsv"),
       col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# # #############################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
ph=ibd

# note this will creat hardcalls, and those will be used for ld estimation, all variants good imputation, it should not be biased...

MEM=2000
for chr in {1..22} 
do
bsub -J"cojo_i" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-o ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_ld1_final_model_stdout \
-e ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_ld1_final_model_stderr \
"/path/to/software/username/./plink2 \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_subset_included_in_${ph}_analysis \
--extract ${path_gwas}post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_for_ld_with_exome_data.tsv \
--make-bed --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_for_ld_with_exome"
done

for chr in {1..22} 
do
echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_ld1_final_model_stdout | grep "Successfully"
done

for chr in {1..22} 
do
bsub -J"cojo_i" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-o ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_ld2_final_model_stdout \
-e ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_ld2_final_model_stderr \
"/path/to/software/username/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_for_ld_with_exome \
--r2 --ld-window-kb '1000000000' --ld-window '1000000000' --ld-window-r2 '0' --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_for_ld_with_exome"
done

for chr in {1..22} 
do
echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_ld2_final_model_stdout | grep "Successfully"
done

######################################################################
# add those where the lead is the coding variant:

exo_same_index<-exo[which(exo$ID %in% all$MarkerName),]
dim(exo_same_index)
# [1] 36 17

all$exonic_variant_in_ld <- NA_character_
all$exonic_variant_in_ld_gene_aac <- NA_character_
all$exonic_variant_in_ld_status <- NA_character_
all$exonic_variant_in_ld_associated_phenotype <- NA_character_

idx <- match(exo_same_index$ID, all$MarkerName)
valid <- !is.na(idx)
all$exonic_variant_in_ld[idx[valid]]                     <- paste(exo_same_index$ID[valid], exo_same_index$ID[valid], 1, sep = "|")
all$exonic_variant_in_ld_gene_aac[idx[valid]]            <- paste(exo_same_index$Gene_symbol[valid], exo_same_index$AAC[valid], sep = ":")
all$exonic_variant_in_ld_status[idx[valid]]              <- exo_same_index$variant_status[valid]
all$exonic_variant_in_ld_associated_phenotype[idx[valid]] <- exo_same_index$Most_Associated_Trait[valid]

all[!is.na(all$exonic_variant_in_ld),c("exonic_variant_in_ld_gene_aac","exonic_variant_in_ld_status")]

dim(all[which(!is.na(all$exonic_variant_in_ld)),])
# [1]  36 219


###########################################################
# Recover additional variants (if they are in LD)

ld <- rbindlist(lapply(1:22, function(chr) {
  f <- sprintf("%s%s.ld", path_ld, sprintf(ld_prefix, chr))
  tmp <- fread(f)
  tmp[tmp$R2 >= 0.6, ]
}))

dim(ld)
# [1] 24  7

ld1<-ld[which(ld$SNP_A %in% exo$ID),c("SNP_A","SNP_B","R2")]
ld2<-ld[which(ld$SNP_B %in% exo$ID),c("SNP_B","SNP_A","R2")]

colnames(ld1)<-c("SNP_exome","SNP","R2")
colnames(ld2)<-c("SNP_exome","SNP","R2")

ld<-rbind(ld1,ld2)
dim(ld)
# [1] 31  3

ld<-ld[which(ld$SNP %in% all$MarkerName),]
dim(ld)
# [1] 20    3

exo_sub <- exo[match(ld$SNP_exome, exo$ID), ]
idx <- match(ld$SNP, all$MarkerName)
valid <- !is.na(idx)
all$exonic_variant_in_ld[idx[valid]]                      <- paste(ld$SNP_exome[valid], ld$SNP[valid], ld$R2[valid], sep = "|")
all$exonic_variant_in_ld_gene_aac[idx[valid]]             <- paste(exo_sub$Gene_symbol[valid], exo_sub$AAC[valid], sep = ":")
all$exonic_variant_in_ld_status[idx[valid]]               <- exo_sub$variant_status[valid]
all$exonic_variant_in_ld_associated_phenotype[idx[valid]]  <- exo_sub$Most_Associated_Trait[valid]

all[!is.na(all$exonic_variant_in_ld),c("exonic_variant_in_ld_gene_aac","exonic_variant_in_ld_status")]


dim(all[which(!is.na(all$exonic_variant_in_ld)),])
# [1]  51 191

all$class_signal_final_exome<-all$class_signal_final

### RECLASSIFY THE SIGNALS BASED ON THE EXONIC DATA

table(all$class_signal_final,all$exonic_variant_in_ld_status)
  #                                                               Novel_variant
  # new_cojo_supervised_gw_significant_multiancestry_known_signal             0
  # new_cojo_supervised_gw_significant_multiancestry_new_signal               1
  # new_cojo_unsupervised_known_signal                                        2
  # new_cojo_unsupervised_new_signal                                          6
                                                               
  #                                                               Reported_variant
  # new_cojo_supervised_gw_significant_multiancestry_known_signal                0
  # new_cojo_supervised_gw_significant_multiancestry_new_signal                  0
  # new_cojo_unsupervised_known_signal                                          30
  # new_cojo_unsupervised_new_signal                                            12

all[which( (all$class_signal_final %in% c("new_cojo_supervised_gw_significant_multiancestry_new_signal") & 
all$exonic_variant_in_ld_status=="Novel_variant")),
c("exonic_variant_in_ld_gene_aac","exonic_variant_in_ld","exonic_variant_in_ld_status")]
#     exonic_variant_in_ld_gene_aac                         exonic_variant_in_ld
# 546                        ORC3:. chr6:87657917:A:G|chr6:87705790:C:T|0.832395
#     exonic_variant_in_ld_status
# 546               Novel_variant

all$class_signal_final_exome[which(all$exonic_variant_in_ld_gene_aac %in% c("ORC3:."))]<-
paste0(all$class_signal_final[which(all$exonic_variant_in_ld_gene_aac %in% c("ORC3:."))],"_new_exonic_variant")



# RETAIN AS KNOWN

# retain as known GWAS signal

all[which( (all$class_signal_final %in% c("new_cojo_unsupervised_known_signal") & 
all$exonic_variant_in_ld_status=="Novel_variant")),
c("exonic_variant_in_ld_gene_aac","exonic_variant_in_ld","exonic_variant_in_ld_status")]
#     exonic_variant_in_ld_gene_aac                         exonic_variant_in_ld
# 62                  PLEKHG5:A414T   chr1:6471529:C:T|chr1:6459964:C:A|0.624958
# 626                   RIPK2:I259T chr8:89772751:T:C|chr8:89834754:G:A|0.658084
#     exonic_variant_in_ld_status
# 62                Novel_variant
# 626               Novel_variant

# retain as novel GWAS signal/

all[which( (all$class_signal_final %in% c("new_cojo_unsupervised_new_signal") & 
all$exonic_variant_in_ld_status=="Novel_variant")),
c("exonic_variant_in_ld_gene_aac","exonic_variant_in_ld","exonic_variant_in_ld_status")]
#     exonic_variant_in_ld_gene_aac
# 187                   PTGER2:C83G
# 340                       IFIH1:.
# 498                    DUSP1:A56T
# 531                 ZC3H12D:K106R
# 549                   PTCD1:R113W
# 616                    DOK2:L138S
#                               exonic_variant_in_ld exonic_variant_in_ld_status
# 187        chr14:52314795:T:G|chr14:52314795:T:G|1               Novel_variant
# 340 chr2:162279995:C:G|chr2:162280432:A:G|0.995039               Novel_variant
# 498        chr5:172770787:C:T|chr5:172770787:C:T|1               Novel_variant
# 531 chr6:149461959:T:C|chr6:149413099:A:G|0.938823               Novel_variant
# 549  chr7:99434906:G:A|chr7:100021222:G:A|0.657185               Novel_variant
# 616   chr8:21911921:A:G|chr8:21853230:T:G|0.889599               Novel_variant

all$class_signal_final_exome[which(all$exonic_variant_in_ld_gene_aac %in% c("PTGER2:C83G","IFIH1:.","DUSP1:A56T","ZC3H12D:K106R","PTCD1:R113W","DOK2:L138S"))]<-
paste0(all$class_signal_final[which(all$exonic_variant_in_ld_gene_aac %in% c("PTGER2:C83G","IFIH1:.","DUSP1:A56T","ZC3H12D:K106R","PTCD1:R113W","DOK2:L138S"))],"_new_exonic_variant")



all[which( (all$class_signal_final %in% c("new_cojo_unsupervised_new_signal") 
& all$exonic_variant_in_ld_status=="Reported_variant")),
c("exonic_variant_in_ld_gene_aac","exonic_variant_in_ld","exonic_variant_in_ld_status")]
#     exonic_variant_in_ld_gene_aac                      exonic_variant_in_ld
# 113                  IL10RA:P295L chr11:117998788:C:T|chr11:117998788:C:T|1- CD_22seq paper
# 222                    NOD2:R703C   chr16:50712018:C:T|chr16:50712018:C:T|1 CD_22seq
# 228                 CARMIL2:V181M   chr16:67646903:G:A|chr16:67646903:G:A|1 - NEW WES/GWAS
# 252                      CCR7:M7V   chr17:40558934:T:C|chr17:40558934:T:C|1 - CD_22seq
# 284                    TYK2:A928V   chr19:10354167:G:A|chr19:10354167:G:A|1 - CD_22seq
# 285                    TYK2:I684S   chr19:10359299:A:C|chr19:10359299:A:C|1 - CD_22seq
# 361                ARHGAP25:E254K     chr2:68813372:G:A|chr2:68813372:G:A|1 - NEW WES/GWAS
# 363                      ADAM17:.       chr2:9521321:A:G|chr2:9521321:A:G|1 - NEW WES/GWAS
# 403           MFNG:MFNG:Arg302Cys   chr22:37470025:G:A|chr22:37470025:G:A|1 - NEW WES/GWAS
# 407           TSPO:TSPO:Thr147Ala   chr22:43162920:A:G|chr22:43162920:A:G|1 - NEW WES/GWAS
# 565                    RCC1L:G30R     chr7:75073650:C:T|chr7:75073650:C:T|1 - CD_22seq
# 582                    DOK2:P274L     chr8:21909729:G:A|chr8:21909729:G:A|1 - CD_22seq

#     exonic_variant_in_ld_gene_aac                      exonic_variant_in_ld
# 117                  IL10RA:P295L chr11:117998788:C:T|chr11:117998788:C:T|1 - CD_22seq paper
# 232                    NOD2:R703C   chr16:50712018:C:T|chr16:50712018:C:T|1 - CD_22seq
# 239                 CARMIL2:V181M   chr16:67646903:G:A|chr16:67646903:G:A|1 - NEW WES/GWAS
# 263                      CCR7:M7V   chr17:40558934:T:C|chr17:40558934:T:C|1 - CD_22seq
# 296                    TYK2:A928V   chr19:10354167:G:A|chr19:10354167:G:A|1 - CD_22seq
# 297                    TYK2:I684S   chr19:10359299:A:C|chr19:10359299:A:C|1 - known
# 375                ARHGAP25:E254K     chr2:68813372:G:A|chr2:68813372:G:A|1 - NEW WES/GWAS
# 377                      ADAM17:.       chr2:9521321:A:G|chr2:9521321:A:G|1 - NEW WES/GWAS
# 418           MFNG:MFNG:Arg302Cys   chr22:37470025:G:A|chr22:37470025:G:A|1 - NEW WES/GWAS
# 422    TSPO;TTLL12:TSPO:Thr147Ala   chr22:43162920:A:G|chr22:43162920:A:G|1 - NEW WES/GWAS
# 599                    RCC1L:G30R     chr7:75073650:C:T|chr7:75073650:C:T|1 - NEW WES/GWAS
# 617                    DOK2:P274L     chr8:21909729:G:A|chr8:21909729:G:A|1 - CD_22seq

all$class_signal_final_exome[which(all$exonic_variant_in_ld_gene_aac %in% c("IL10RA:P295L","TYK2:I684S","DOK2:P274L","CCR7:M7V","NOD2:R703C","TYK2:A928V","NOD2:.","GPR35:T108M"))]<-
paste0(all$class_signal_final[which(all$exonic_variant_in_ld_gene_aac %in% c("IL10RA:P295L","TYK2:I684S","DOK2:P274L","CCR7:M7V","NOD2:R703C","TYK2:A928V","NOD2:.","GPR35:T108M"))],"_known_exonic_variant")

all$class_signal_final_exome[which(all$exonic_variant_in_ld_gene_aac %in% c("ARHGAP25:E254K","ADAM17:.","MFNG:MFNG:Arg302Cys","CARMIL2:V181M","IL27:L119P","NEU4:.","TSPO;TTLL12:TSPO:Thr147Ala","RCC1L:G30R"))]<-
paste0(all$class_signal_final[which(all$exonic_variant_in_ld_gene_aac %in% c("ARHGAP25:E254K","ADAM17:.","MFNG:MFNG:Arg302Cys","CARMIL2:V181M","IL27:L119P","NEU4:.","TSPO;TTLL12:TSPO:Thr147Ala","RCC1L:G30R"))],"_new_exonic_variant")


all[which( (all$class_signal_final %in% c("new_cojo_unsupervised_new_signal") 
& all$exonic_variant_in_ld_status=="Reported_variant")),
c("exonic_variant_in_ld_gene_aac","exonic_variant_in_ld","exonic_variant_in_ld_status","class_signal_final_exome")]




table(all$class_signal_final_exome,useNA="ifany")
#                  new_cojo_supervised_gw_significant_multiancestry_known_signal 
#                                                                              5 
#                    new_cojo_supervised_gw_significant_multiancestry_new_signal 
#                                                                             14 
# new_cojo_supervised_gw_significant_multiancestry_new_signal_new_exonic_variant 
#                                                                              1 
#                                             new_cojo_unsupervised_known_signal 
#                                                                            233 
#                        new_cojo_unsupervised_known_signal_known_exonic_variant 
#                                                                              1 
#                                               new_cojo_unsupervised_new_signal 
#                                                                            350 
#                          new_cojo_unsupervised_new_signal_known_exonic_variant 
#                                                                              6 
#                            new_cojo_unsupervised_new_signal_new_exonic_variant 
#                                                                             9 

all[which(all$class_signal_final_exome=="new_cojo_unsupervised_new_signal_known_exonic_variant"),c("exonic_variant_in_ld","exonic_variant_in_ld_gene_aac","class_signal_final_exome")]
#                          exonic_variant_in_ld exonic_variant_in_ld_gene_aac
# 117 chr11:117998788:C:T|chr11:117998788:C:T|1                  IL10RA:P295L
# 232   chr16:50712018:C:T|chr16:50712018:C:T|1                    NOD2:R703C
# 263   chr17:40558934:T:C|chr17:40558934:T:C|1                      CCR7:M7V
# 296   chr19:10354167:G:A|chr19:10354167:G:A|1                    TYK2:A928V
# 285   chr19:10359299:A:C|chr19:10359299:A:C|1                    TYK2:I684S
# 617     chr8:21909729:G:A|chr8:21909729:G:A|1                    DOK2:P274L
#                                  class_signal_final_exome
# 117 new_cojo_unsupervised_new_signal_known_exonic_variant
# 232 new_cojo_unsupervised_new_signal_known_exonic_variant
# 263 new_cojo_unsupervised_new_signal_known_exonic_variant
# 296 new_cojo_unsupervised_new_signal_known_exonic_variant
# 617 new_cojo_unsupervised_new_signal_known_exonic_variant

dim(all[which(all$class_signal_final_exome %in% c("new_cojo_unsupervised_new_signal_new_exonic_variant","new_cojo_supervised_gw_significant_multiancestry_new_signal_new_exonic_variant")),c("exonic_variant_in_ld","exonic_variant_in_ld_gene_aac","class_signal_final_exome")])
# [1] 10  3

all[which(all$class_signal_final_exome %in% c("new_cojo_unsupervised_new_signal_new_exonic_variant","new_cojo_supervised_gw_significant_multiancestry_new_signal_new_exonic_variant")),c("exonic_variant_in_ld","exonic_variant_in_ld_gene_aac","class_signal_final_exome")]
#     exonic_variant_in_ld_gene_aac
# 178                   PTGER2:C83G
# 228                 CARMIL2:V181M
# 327                       IFIH1:.
# 361                ARHGAP25:E254K
# 363                      ADAM17:.
# 403           MFNG:MFNG:Arg302Cys
# 511                 ZC3H12D:K106R
# 526                        ORC3:.
# 565                    RCC1L:G30R
# 581                    DOK2:L138S

#                                                           class_signal_final_exome
# 185                            new_cojo_unsupervised_new_signal_new_exonic_variant
# 237                            new_cojo_unsupervised_new_signal_new_exonic_variant
# 338                            new_cojo_unsupervised_new_signal_new_exonic_variant
# 373                            new_cojo_unsupervised_new_signal_new_exonic_variant
# 375                            new_cojo_unsupervised_new_signal_new_exonic_variant
# 416                            new_cojo_unsupervised_new_signal_new_exonic_variant
# 420                            new_cojo_unsupervised_new_signal_new_exonic_variant
# 496                            new_cojo_unsupervised_new_signal_new_exonic_variant
# 529                            new_cojo_unsupervised_new_signal_new_exonic_variant
# 544 new_cojo_supervised_gw_significant_multiancestry_new_signal_new_exonic_variant
# 547                            new_cojo_unsupervised_new_signal_new_exonic_variant
# 597                            new_cojo_unsupervised_new_signal_new_exonic_variant
# 614                            new_cojo_unsupervised_new_signal_new_exonic_variant


# using most associated trait as reported by Exome
table(all$exonic_variant_in_ld_associated_phenotype[which(!is.na(all$exonic_variant_in_ld))],all$phenotype[which(!is.na(all$exonic_variant_in_ld))])
  #     CD IBD_saturated IBD_unsaturated UC
  # CD  12             7               1  0
  # IBD  4             5              15  3
  # UC   0             4               0  1

table(all$class_signal_final_exome,all$class_signal_final,useNA="ifany")



table(all$class_signal_final_exome,useNA="ifany")
#                  new_cojo_supervised_gw_significant_multiancestry_known_signal 
#                                                                              5 
#                    new_cojo_supervised_gw_significant_multiancestry_new_signal  # new
#                                                                             14 
# new_cojo_supervised_gw_significant_multiancestry_new_signal_new_exonic_variant  # new
#                                                                              1 
#                                             new_cojo_unsupervised_known_signal 
#                                                                            233 
#                        new_cojo_unsupervised_known_signal_known_exonic_variant 
#                                                                              1  
#                                               new_cojo_unsupervised_new_signal # new
#                                                                            350 
#                          new_cojo_unsupervised_new_signal_known_exonic_variant 
#                                                                              6 
#                            new_cojo_unsupervised_new_signal_new_exonic_variant # new
#                                                                             9 


# novel
14+1+350+9
# [1] 374
# known:
5+233+1+6
# [1] 245



exonic<-all[which(!is.na(all$exonic_variant_in_ld)),c("MarkerName","phenotype",
"exonic_variant_in_ld","exonic_variant_in_ld_gene_aac" ,"exonic_variant_in_ld_status" ,"class_signal_final_exome")]

dim(exonic)
# [1] 51  5

# fwrite(exonic, out_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


dim(all[grep("new_exonic",all$class_signal_final_exome),])
# 10

fwrite(all, out_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

q("no")

# plot these results using:
# /path/to/user/git/IIBDGC_GWAS/scripts/other/plot_effect_sizes_CD_vs_UC_known_new_regions_consensus_list.R
# /path/to/user/git/IIBDGC_GWAS/scripts/other/plot_figure_1_version_5.R

# Run before starting R (in shell, not here):
# singularity exec iibdgc_postprocess_10_singularity.sif
MEM=6000


bsub -J"plot1" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-o ${path_gwas}post_imputation/2022/log/plot_figure_1_stdout \
-e ${path_gwas}post_imputation/2022/log/plot_figure_1_stderr \
" Rscript /path/to/user/git/IIBDGC_GWAS/scripts/other/plot_figure_1_version_5.R > \
${path_gwas}post_imputation/2022/log/plot_figure_1_version_5.Rout"

