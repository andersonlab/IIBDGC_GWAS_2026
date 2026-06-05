# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
## ADD COLOC RESULTS FROM SNAKEMAKE_COLOC PIPELINE
# See ~/git/snakemake_colocalisation

# singularity exec iibdgc_postprocess_10_singularity.sif

##########################################################################
# 1.- EXTRACT LD BETWEEN ALL VARIANTS, eQTL and GWAS leads

MEM=25000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q yesterday R

rm(list=ls())

library(data.table)
library(stringr)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

all<-fread(paste0(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants.tsv.gz"))
dim(all)
# [1] 619 184

#################################################################################################
# map to rename conditions - see /path/to/user/git/snakemake_colocalisation/scripts/create_maps_for_IIBDGC.R

map<-fread("/path/to/project",head=T)

coloc_colnames<-c("condition_name","phenotype_id","snp_id",".row","nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf","qtl_pval","gwas_pval","qtl_lead","gwas_lead","chr","gwas_lead_pos","qtl_lead_pos","gwas_trait")

pheno<-c("cd","uc","ibd")

marker_list<-c()

for (ph in pheno) {

    print(ph)
    ph_upper<-toupper(ph)

    tmp1<-fread(paste0("/path/to/project",ph_upper,".gz"),head=F)
    tmp2<-fread(paste0("/path/to/project",ph_upper,".gz"),head=F)
    # do not include edQTL data
    # tmp3<-fread(paste0("/path/to/project",ph_upper,".gz"),head=F)
    tmp4<-fread(paste0("/path/to/project",ph_upper,".gz"),head=F)
    tmp5<-fread(paste0("/path/to/project",ph,".gz"),head=F,skip=1)
    tmp6<-fread(paste0("/path/to/project",ph_upper,".gz"),head=F)
    tmp7<-fread(paste0("/path/to/project",ph_upper,".gz"),head=F)

    # pQTL_sparc chr columns lack "chr" prefix - add before rbind
    tmp7$V3 <-paste0("chr",tmp7$V3)
    tmp7$V13<-paste0("chr",tmp7$V13)
    tmp7$V14<-paste0("chr",tmp7$V14)

    tmp<-rbind(tmp1,tmp2,tmp4,tmp5,tmp6,tmp7)
    rm(tmp1,tmp2,tmp4,tmp5,tmp6,tmp7)

    colnames(tmp)<-coloc_colnames
    tmp$condition_name<-gsub("-","_",tmp$condition_name)
    tmp<-merge(tmp,map,by.x="condition_name",by.y="id_map")

    print(nrow(tmp))
    tmp$PP.H4.abf<-as.numeric(tmp$PP.H4.abf)
    tmp<-tmp[!is.na(tmp$PP.H4.abf) & tmp$PP.H4.abf>=0.8,]
    print(nrow(tmp))

    print(table(tmp$study_label))

    marker_list<-c(marker_list,
                   gsub("_",":",tmp$gwas_lead),
                   gsub("_",":",tmp$qtl_lead))
    rm(tmp)
}

marker_list<-unique(c(marker_list, all$MarkerName))

length(marker_list)
# [1] 7994


write.table(marker_list,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_MarkerNames",sep=""),
       col.names=F,row.names=F,quote=F,sep="\t")

q("no")

###############################################################

path_gwas="/path/to/ibdgwas/IIBDGC/"
ph="ibd"

MEM=2000

for chr in {1..22}
do
bsub -J"kg" -M"$MEM" -n8 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/log/stderr_ld_coloc_1_${chr} \
-o ${path_gwas}pre_imputation/QC/log/stdout_ld_coloc_1_${chr} \
"/path/to/software/username/./plink2 \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_subset_included_in_ibd_analysis \
--extract ${path_gwas}post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_MarkerNames \
--make-bed --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/chr${chr}_list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_for_ld"
done


for chr in {1..22}
do
echo ${chr} && tail -50 ${path_gwas}pre_imputation/QC/log/stdout_ld_coloc_1_${chr} | grep "Successfully"
done

for chr in {1..22}
do
echo ${chr} && ls -la ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/chr${chr}_list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_for_ld.bim
done

for chr in {1..22}
do
bsub -J"kg" -M"$MEM" -n8 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/log/stderr_ld_coloc_2_${chr} \
-o ${path_gwas}pre_imputation/QC/log/stdout_ld_coloc_2_${chr} \
"/path/to/software/username/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/chr${chr}_list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_for_ld \
--r2 --ld-window-kb '1000000000' --ld-window '1000000000' --ld-window-r2 '0' --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/chr${chr}_list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_for_ld"
done

for chr in {1..22}
do
echo ${chr} && tail -50 ${path_gwas}pre_imputation/QC/log/stdout_ld_coloc_2_${chr} | grep "Successfully"
done

for chr in {1..22}
do
echo ${chr} && wc -l ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/chr${chr}_list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_for_ld.ld
done


###############################################################

# 2.-  SUBSET THE FILES BY LD LEADS:

# retain only:
# - colocs where IIBDGC index variant and eqtl lead have r2>=0.6
# - colocs where GWAS lead and eqtl lead have r2>=0.6; and IIBDGC index variant and GWAS lead OR eqtl lead have r2>=0.6

MEM=3000

path_gwas="/path/to/ibdgwas/IIBDGC/"

bsub -J"ph" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}post_imputation/2022/log/subset_files_by_ld_leads_allpheno_stderr \
-o ${path_gwas}post_imputation/2022/log/subset_files_by_ld_leads_allpheno_stdout \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/coloc_subset_results_by_ld_leads_noFM.R ${ph} > \
${path_gwas}post_imputation/2022/log/coloc_subset_results_by_ld_leads_noFM.Rout"

wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_filtered_by_gwas_r2_0.6_allpheno_noFM.tsv
# 10421 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_filtered_by_gwas_r2_0.6_allpheno_noFM.tsv

###############################################################

# 3.- ADD BETA per coloc:

MEM=2000

bsub -J"ph" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}post_imputation/2022/log/add_beta_allpheno_stderr \
-o ${path_gwas}post_imputation/2022/log/add_beta_allpheno_stdout \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/coloc_subset_add_eqtl_beta_noFM.R > \
${path_gwas}post_imputation/2022/log/coloc_subset_add_eqtl_beta_noFM.Rout"

wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_filtered_by_gwas_r2_0.6_allpheno_noFM_with_beta.tsv
# 10421 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_filtered_by_gwas_r2_0.6_allpheno_noFM_with_beta.tsv

###############################################################

# 4 - CREATE MATRICES PER PHENOTYPE

pheno=(cd uc ibd)
path_gwas="/path/to/ibdgwas/IIBDGC/"

MEM=1500

for ph in ${pheno[@]}
do
bsub -J"ph" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}post_imputation/2022/log/create_enrichment_zscore_matrices_noFM_${ph}_stderr \
-o ${path_gwas}post_imputation/2022/log/create_enrichment_zscore_matrices_noFM_${ph}_stdout \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/coloc_subset_create_weight_matrices_noFM.R ${ph} > \
${path_gwas}post_imputation/2022/log/coloc_subset_create_weight_matrices_noFM_${ph}.Rout"
done


for ph in ${pheno[@]}
do
echo ${ph} && \
tail -50  ${path_gwas}post_imputation/2022/log/create_enrichment_zscore_matrices_noFM_${ph}_stdout | grep -E "Exited|Successfully"
done


for ph in ${pheno[@]}
do
echo ${ph} && \
tail -50 ${path_gwas}post_imputation/2022/log/coloc_subset_create_weight_matrices_noFM_${ph}.Rout
done 


for ph in ${pheno[@]}
do
ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_matrices/matrix_${ph}_zscore_genes_tissue_noFM.tsv.gz
done

for ph in ${pheno[@]}
do
zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/coloc_results/coloc_matrices/matrix_${ph}_zscore_genes_tissue_noFM.tsv.gz | wc -l
done

# 652
# 652
# 652
