# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
############################################################################################################################
############################################################################################################################


# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=15000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q yesterday R

library(data.table)
library(ggplot2)
library(colorspace)
library(dplyr)
library(tibble)
library(ggrepel)
library(ggpubr)
rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

all<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants_join_conditinal_model.tsv.gz",sep=""))
dim(all)
# [1]  619 205

fm<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/fine_mapping/tier1_EUR_array_data_fine_mapping_results_08_16_2025_list_of_variants.csv",sep=""))
fm<-as.data.frame(fm)
dim(fm)
# 7221

table(fm$lead)
# FALSE  TRUE 
#  6953   268 

fm<-fm[which(fm$lead==TRUE),]
dim(fm)
# [1] 268  51

fm$region_signal<-paste(fm$region_label,fm$CS_index,sep="_")

# see ~/git/IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_40_evaluate_fine_mapping_results.R
for (chr in 1:22) {
  ld.tmp<-read.table(paste0(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr",chr,"_list_tier_2_unsupervised_conditional_variants_for_ld_with_fm_with_old_data.ld"),head=T)
  ld.tmp<-ld.tmp[which(ld.tmp$R2>=0.1),]

  if(chr==1) {
    ld<-ld.tmp
  } else {
    ld<-rbind(ld,ld.tmp)
  }
  rm(ld.tmp)
}
dim(ld)
# [1] 1305    7

### map the FM lead variants to cojo leads:
fm<-fm[,c("MarkerName","region_signal","CS_index","alpha_lead","lead")]

# if lead is the same:
fm1<-fm[which(fm$MarkerName %in% all$MarkerName),]
dim(fm1)
# [1] 83  3
fm1$ld_cojo_lead_fm_lead<-paste(fm1$MarkerName,"1",sep="|")
fm1$ld_cojo_lead_fm_lead_r2<-1

# link the rest by LD:
fm2<-fm[which(!fm$MarkerName %in% all$MarkerName),]
dim(fm2)
# [1] 185   3

fm2$ld_cojo_lead_fm_lead<-NA
fm2$ld_cojo_lead_fm_lead_r2<-NA

for (i in c(1:nrow(fm))) {

    # get LD:
    tmp<-ld[which( (ld$SNP_A %in% fm2$MarkerName[i]) | (ld$SNP_B %in% fm2$MarkerName[i])),]
    tmp<-tmp[which( (tmp$SNP_A %in% all$MarkerName) | (tmp$SNP_B %in% all$MarkerName)),]

    tmp$ids<-paste(as.character(tmp$SNP_A),as.character(tmp$SNP_B),tmp$R2,sep="|")

    if(nrow(tmp)>0) {
        fm2$ld_cojo_lead_fm_lead[i]<-paste(tmp$ids,collapse=";")
        fm2$ld_cojo_lead_fm_lead_r2[i]<-paste(tmp$R2,collapse=";")
    }
    rm(tmp)
    
}

fm<-rbind(fm1,fm2)
# [1] 268  51

dim(fm[which(is.na(fm$ld_cojo_lead_fm_lead)),])
# [1] 19  4

fm_all<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/fine_mapping/tier1_EUR_array_data_fine_mapping_results_08_16_2025_list_of_variants.csv",sep=""))
fm_all<-as.data.frame(fm_all)
dim(fm_all)
# 7221


head(fm_all)

fm_all$region_signal<-paste(fm_all$region_label,fm_all$CS_index,sep="_")

### map the FM lead variants to cojo leads:
fm_all<-fm_all[,c("MarkerName","region_signal","CS_index","alpha_lead","lead")]


fm$ld_cojo_lead_credible_set<-NA
fm$ld_cojo_lead_credible_set_r2<-NA

for (i in 1:nrow(fm)) {

    tmp<-fm_all[which(fm_all$region_signal==fm$region_signal[i]),]

    tmp_ld<-ld[which(ld$SNP_A %in% tmp$MarkerName | ld$SNP_B %in% tmp$MarkerName),]
    tmp_ld<-tmp_ld[which( (tmp_ld$SNP_A %in% all$MarkerName) | (tmp_ld$SNP_B %in% all$MarkerName)),]

    if(nrow(tmp_ld)>0) {
        # retain only one
        tmp_ld<-tmp_ld[which(tmp_ld$R2==max(tmp_ld$R2)),]
        tmp_ld$ids<-paste(as.character(tmp_ld$SNP_A),as.character(tmp_ld$SNP_B),tmp_ld$R2,sep="|")
        fm$ld_cojo_lead_credible_set[i]<-paste(tmp_ld$ids,collapse=";")
        fm$ld_cojo_lead_credible_set_r2[i]<-max(tmp_ld$R2)
    }

    rm(tmp,tmp_ld)

}

# overwrite those with the same lead
fm$ld_cojo_lead_credible_set[which(fm$MarkerName %in% all$MarkerName)]<-paste(fm$MarkerName[which(fm$MarkerName %in% all$MarkerName)],"1",sep="|")
fm$ld_cojo_lead_credible_set_r2[which(fm$MarkerName %in% all$MarkerName)]<-1




dim(fm[which(is.na(fm$ld_cojo_lead_fm_lead)),])
# [1] 19  7


table(fm$CS_index[which(is.na(fm$ld_cojo_lead_fm_lead))])
# 1 2 3 4 7 9 
# 9 4 2 2 1 1

fm[which(is.na(fm$ld_cojo_lead_fm_lead) & fm$CS_index==1),]
#                MarkerName            region_signal CS_index alpha_lead lead
# 791  chr1:212839129:T:TTA  1_212213971_213339130_1        1 0.69783623 TRUE - no region
# 794    chr1:226627077:G:C  1_226075302_227400891_1        1 0.23686940 TRUE - best variant chr1:226754854:T:C
# 2361     chr5:4422049:C:T      5_3922049_4922050_1        1 1.00000000 TRUE - no region
# 4872  chr11:125059389:C:T 11_124471689_125471689_1        1 0.04424202 TRUE - no region
# 5579  chr14:105753554:C:T 14_104456243_106231744_1        1 0.07565553 TRUE - chr14:104956890:A:G, chr14:105795065:A:T,chr14:105968520:G:A
# 5604   chr15:44802938:T:C   15_44302938_45302939_1        1 0.99875021 TRUE - no region
# 5862   chr16:31926275:A:C   16_29971173_31965789_1        1 1.00000000 TRUE - chr16:30471173:T:C, chr16:30760683:GAA:G, chr16:31357553:T:C
# 6511   chr18:79395695:T:G   18_78960616_80020064_1        1 0.31167180 TRUE - chr18:79490197:G:C
# 6633    chr20:6084192:T:C     20_5613242_6613243_1        1 0.95642433 TRUE - no region

tmp<-fm[which(is.na(fm$ld_cojo_lead_fm_lead) & fm$CS_index==1),]


fm_all<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/fine_mapping/tier1_EUR_array_data_fine_mapping_results_08_16_2025_list_of_variants.csv",sep=""))
fm_all<-as.data.frame(fm_all)
fm_all[which(fm_all$MarkerName %in% tmp$MarkerName),]

# GET CONDITIONAL SUMMARY STATS FOR THEM:

## add the conditional resutls:

pheno<-c("ibd","cd","uc")
for (ph in pheno) {

  for (chr in 1:22) {

    tmp<-read.table(paste0(path_gwas,"post_imputation/2022/analysis/conditional_analysis/final_independent_model/",chr,"_",ph,"_independent_index_variants_eur_tier_2_ld_allarrays_conditioning_on_final_independent_model_0.1_other_variants.cma.cojo"),head=T)
    tmp<-tmp[which(tmp$SNP %in% fm_all$MarkerName),c("SNP","bC","bC_se","pC")]

    if(chr==1) {
      dat_tmp<-tmp
    } else {
      dat_tmp<-rbind(dat_tmp,tmp)
    }
    rm(tmp)
  }

  colnames(dat_tmp)[2:ncol(dat_tmp)]<-paste(colnames(dat_tmp)[2:ncol(dat_tmp)],ph,sep="_")
  
  print(dim(dat_tmp))

  if(ph=="ibd") {
    dat<-dat_tmp
  } else {
    dat<-merge(dat,dat_tmp,by="SNP",all=T)
  }
  rm(dat_tmp)
}


head(dat)
dat<-dat[which(dat$SNP %in% fm$MarkerName),]

fm<-merge(fm,dat,by.x="MarkerName",by.y="SNP",all.x=T)
fm[which(is.na(fm$ld_cojo_lead_fm_lead) & (fm$pC_ibd<5E-8 | fm$pC_cd<5E-8 | fm$pC_uc<5E-8)),]
#               MarkerName           region_signal CS_index alpha_lead lead
# 25  chr1:212839129:T:TTA 1_212213971_213339130_1        1 0.69783623 TRUE
# 47    chr10:62617743:G:A  10_60240494_64101958_4        4 0.17661424 TRUE
# 98    chr16:50537867:C:T  16_49599118_51605296_9        9 0.48880951 TRUE
# 103   chr16:50712383:A:C  16_49599118_51605296_4        4 0.99998835 TRUE
# 215   chr5:39999103:T:TA   5_39390457_41455460_3        3 0.04942866 TRUE
#     ld_cojo_lead_fm_lead ld_cojo_lead_fm_lead_r2 ld_cojo_lead_credible_set
# 25                  <NA>                    <NA>                      <NA>
# 47                  <NA>                    <NA>                      <NA>
# 98                  <NA>                    <NA>                      <NA>
# 103                 <NA>                    <NA>                      <NA>
# 215                 <NA>                    <NA>                      <NA>
#     ld_cojo_lead_credible_set_r2     bC_ibd  bC_se_ibd      pC_ibd      bC_cd
# 25                            NA -0.0361360 0.00780008 3.60799e-06 -0.0580920
# 47                            NA  0.1036990 0.01000030 3.40968e-25  0.0855829
# 98                            NA  0.5398450 0.11840400 5.13119e-06  0.7372300
# 103                           NA  0.0760630 0.02290010 8.95287e-04  0.3271070
# 215                           NA  0.0418294 0.00760029 3.72000e-08  0.0500667
#      bC_se_cd       pC_cd      bC_uc   bC_se_uc       pC_uc
# 25  0.0103001 1.70124e-08 -0.0200715 0.00970002 3.85252e-02
# 47  0.0132001 8.96084e-11  0.1116970 0.01230040 1.07795e-19
# 98  0.1254040 4.13153e-09  0.0547504 0.16480000 7.39720e-01
# 103 0.0325010 7.92767e-24 -0.1016240 0.02910020 4.79033e-04
# 215 0.0110004 5.33000e-06  0.0325470 0.00940013 5.35393e-04


fwrite(fm,paste0(path_gwas,"post_imputation/2022/analysis/final_tables/fine_mapping_captured_by_cojo_631_indep.txt.gz"),
col.names=T,row.names=F,quote=F,sep="\t")

q("no")



# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=4000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q yesterday R

library(data.table)
library(ggplot2)
library(colorspace)
library(dplyr)
library(tibble)
library(ggrepel)
library(ggpubr)
rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"


all<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model.tsv.gz",sep=""))
dim(all)
# [1] 653 163

###################
# link to effector genes:

fm<-fread(paste0(path_gwas,"post_imputation/2022/analysis/final_tables/fine_mapping_captured_by_cojo_653_indep.txt.gz"))
fm$MarkerName_indep_cojo<-NA

for (i in c(1:nrow(fm))) {

  ids<-unlist(strsplit(fm$ld_cojo_lead_fm_lead[i],"\\|"))
  ids<-unlist(strsplit(ids,"\\;"))
  

  if (length(ids)>0) {
    fm$MarkerName_indep_cojo[i]<-paste(all$MarkerName[which(all$MarkerName %in% ids)],collapse="|")
    rm(ids)
  }

  # if (length(ids)>3) {
  #   print(i)
  # }


}



list_genes<-fread(paste0(path_gwas,"post_imputation/2022/analysis/final_tables/list_effector_genes_coloc_mr_exonic_with_coloc_rank1_noFM_with_otar_evidence_with_ligand_receptor_with_protein_complexes_with_monogenic_with_gps.tsv.gz"))
list_genes<-as.data.frame(list_genes)

head(fm)
fm<-as.data.frame(fm)
fm$effector_gene<-NA
fm$effector_gene_mr_exonic_rank1<-NA

for (i in 1:nrow(fm)) {

  ids<-unlist(strsplit(fm$MarkerName_indep_cojo[i],"\\|"))
  ids<-unlist(strsplit(ids,"\\|"))

  tmp<-list_genes[grep(paste(ids,collapse="|"),list_genes$independent_index_variants),]
  fm$effector_gene[i]<-paste(tmp$list_genes,collapse="|")

  tmp<-tmp[which(tmp$mr==1 | tmp$exonic==1 | tmp$coloc_rank1==1),]
  fm$effector_gene_mr_exonic_rank1[i]<-paste(tmp$list_genes,collapse="|")
  

}

fwrite(fm,paste0(path_gwas,"post_imputation/2022/analysis/final_tables/fine_mapping_captured_by_cojo_633_indep.txt.gz"),
col.names=T,row.names=F,quote=F,sep="\t")

# subset the ones highlighted in table 1:
t1<-fread("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/fine_mapping/table_1_finemapping.csv")

t1<-merge(t1,fm[c("MarkerName","effector_gene","effector_gene_mr_exonic_rank1")],by="MarkerName")


fwrite(t1,paste0(path_gwas,"post_imputation/2022/analysis/final_tables/table_1_finemapping_with_effector_genes.tsv.gz"),
col.names=T,row.names=F,quote=F,sep="\t")