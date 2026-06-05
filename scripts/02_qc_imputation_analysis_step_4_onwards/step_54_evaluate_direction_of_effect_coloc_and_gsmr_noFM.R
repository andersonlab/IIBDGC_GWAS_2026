# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
## ADD COLOC RESULTS FROM SNAKEMAKE_COLOC PIPELINE
# See ~/git/snakemake_colocalisation

# singularity exec iibdgc_postprocess_10_singularity.sif

############################################################################################
# 1.- EXTRACT ALL VARIANTS IN LD FOR EACH ONE OF OUR FINAL LIST OF INDEX VARIANTS (AT R2 0.6)

path_gwas="/path/to/ibdgwas/IIBDGC/"
MEM=4500

variants=($(zcat ${path_gwas}post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants.tsv.gz | cut -f1))
chromosomes=($(zcat ${path_gwas}post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants.tsv.gz | cut -f13))

rm  ~/tmp.txt
for i in {1..653} 
do
echo ${i} >> ~/tmp.txt && echo ${variants[i]} >> ~/tmp.txt && echo ${chromosomes[i]} >> ~/tmp.txt
done


for i in {1..653} 
do
bsub -J"ld" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-o ${path_gwas}post_imputation/2022/log/allarray_chr${chromosomes[i]}_${variants[i]}_final_index_ld_stdout \
-e ${path_gwas}post_imputation/2022/log/allarray_chr${chromosomes[i]}_${variants[i]}_final_index_ld_stderr \
"/path/to/software/username/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chromosomes[i]}_subset_included_in_ibd_analysis \
--show-tags <(echo ${variants[i]}) \
--tag-r2 0.6 \
--tag-kb 2000 \
--out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chromosomes[i]}_${variants[i]}_list_653_index_ld_r2_0.6"
done



rm ~/tmp.txt
for i in {1..653} 
do
echo ${i} >> ~/tmp.txt && echo ${variants[i]} >> ~/tmp.txt && echo ${chromosomes[i]} >> ~/tmp.txt && tail -50  ${path_gwas}post_imputation/2022/log/allarray_chr${chromosomes[i]}_${variants[i]}_final_index_ld_stdout | grep -E "Successfully|Exited" >> ~/tmp.txt
done



##################################################################################################################################
# 2.- FOR ALL VARIANTS IN LD WITH LEAD AT R2 = 0.6, ESTIMATE THE DIRECTION OF EFFECT BY USING THE MEAN OF THE WALD TEST

MEM=3800
pheno=(ibd cd uc)

for chr in {1..22}
do 
for ph in ${pheno[@]}
do
bsub -J"ph" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}post_imputation/2022/log/add_direction_effect_${chr}_${ph}_stderr \
-o ${path_gwas}post_imputation/2022/log/add_direction_effect_${chr}_${ph}_stdout \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/evaluate_direction_effect_QTL_risk_variant_noFM.R ${chr} ${ph} > \
${path_gwas}post_imputation/2022/log/evaluate_direction_effect_QTL_risk_variant_noFM_${chr}_${ph}.Rout"
done
done

# CONTINUE HERE

# 1_ibd
# 1_cd
# 1_uc
# 5_ibd
# 11_ibd
# 11_cd


for chr in {1..22}
do 
echo ${chr} && for ph in ${pheno[@]}
do
echo ${ph} && tail -50 ${path_gwas}post_imputation/2022/log/add_direction_effect_${chr}_${ph}_stdout | grep -E "Successfully|Exited"
done
done

for ph in ${pheno[@]}
do
ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/coloc_results/tmp_coloc_filtered_by_gwas_r2_0.6_allpheno_with_beta_with_wald_ratio_*_${ph}.tsv
done
 
for ph in ${pheno[@]}
do
echo ${ph} && wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/coloc_results/tmp_coloc_filtered_by_gwas_r2_0.6_allpheno_with_beta_with_wald_ratio_*_${ph}.tsv
done

# ibd
#    6991 total
# cd
#    6152 total
# uc
#    5399 total

##################################################################################################################################
# 3.- PLOT THE RESULTS, AND ADD DIRECION OF EFFECT PER GENE VIA COLOC

~/git/IIBDGC_GWAS/scripts/other/plot_direction_effect_QTL_risk_variant.R



##################################################################################################################################
##################################################################################################################################

# run the analysis, but this time just including the GWAS lead variant (if present in the eQTL study)

##################################################################################################################################
# 4.- HOW MANY COLOCS WE RETAIN IF WE JUST USE THE GWAS LEAD, AND HOW THE ESTIMATES WILL COMPARE TO USING THE MEAN OF VARIANTS IN LD


# MEM=4500
# pheno=(ibd cd uc)

# for chr in {1..2}
# do 
# for ph in ${pheno[@]}
# do
# bsub -J"ph" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
# -e ${path_gwas}post_imputation/2022/log/add_direction_effect_keeping_only_gwas_lead_${chr}_${ph}_stderr \
# -o ${path_gwas}post_imputation/2022/log/add_direction_effect_keeping_only_gwas_lead_${chr}_${ph}_stdout \
# "Rscript ~/git/IIBDGC_GWAS/scripts/other/evaluate_direction_effect_QTL_risk_variant_keeping_only_gwas_lead.R ${chr} ${ph} > \
# ${path_gwas}post_imputation/2022/log/evaluate_direction_effect_QTL_risk_variant_keeping_only_gwas_lead_${chr}_${ph}.Rout"
# done
# done


# for chr in {1..22}
# do 
# echo ${chr} && for ph in ${pheno[@]}
# do
# echo ${ph} && tail -50 ${path_gwas}post_imputation/2022/log/add_direction_effect_keeping_only_gwas_lead_${chr}_${ph}_stdout | grep -E "Successfully"
# done
# done

# for ph in ${pheno[@]}
# do
# ls -la ${path_gwas}post_imputation/2022/analysis/metaanalysis/coloc_results/tmp_coloc_filtered_by_gwas_r2_0.6_allpheno_with_beta_with_wald_ratio_*_${ph}_only_gwas_lead.tsv | wc -l
# done
 
# for ph in ${pheno[@]}
# do
# echo ${ph} && wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/coloc_results/tmp_coloc_filtered_by_gwas_r2_0.6_allpheno_with_beta_with_wald_ratio_*_${ph}_only_gwas_lead.tsv
# done

# # ibd
# #    7050 total
# # cd
# #    6211 total
# # uc
# #    5636 total


##################################################################################################################################
# 4.- Extract the direction of effect from MR data (plots already created)


# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=2000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(stringr)

rm(list=ls())
path_gwas<-"/path/to/ibdgwas/IIBDGC/"


mr<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_with_fm_signals_no_overlap_with_mr_results.tsv.gz",sep=""),
head=T)

mr_genes<-c(mr$mr_gene_ibd,mr$mr_gene_cd,mr$mr_gene_uc)
mr_genes<-mr_genes[!mr_genes %in% c("")]
mr_genes<-unlist(strsplit(mr_genes, "\\|"))
mr_genes<-as.data.frame(mr_genes[!mr_genes %in% c("")])
colnames(mr_genes)<-"list_genes"
mr_genes<-mr_genes[!duplicated(mr_genes),,drop=F]

reg<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_regions_420_cojo_supervised_plus_unsupervised_with_gene_names_with_GSMR_results.tsv.gz",sep=""))

mr_genes$direction_effect_gsmr<-NA

for (i in 1:nrow(mr_genes)) {
    
    tmp1<-reg[grep(mr_genes$list_genes[i],reg$pqtl_mr_cd),"pqtl_mr_cd"]
    tmp2<-reg[grep(mr_genes$list_genes[i],reg$pqtl_mr_uc),"pqtl_mr_uc"]
    tmp3<-reg[grep(mr_genes$list_genes[i],reg$pqtl_mr_ibd),"pqtl_mr_ibd"]

    tmp<-c(tmp1$pqtl_mr_cd,tmp2$pqtl_mr_uc,tmp3$pqtl_mr_ibd)
    tmp<-unlist(strsplit(tmp, "\\|"))

    tmp<-tmp[grep(mr_genes$list_genes[i],tmp)]

    betas<-as.numeric(matrix(unlist(strsplit(tmp, ",")),ncol=2)[3,])

    mr_genes$direction_effect_gsmr[i]<-max(betas,na.rm=T)


}

## combine all sources with the list of genes:

# LOAD THE LIST OF EFFECTOR GENES :
list_genes<-fread(paste0(path_gwas,"post_imputation/2022/analysis/final_tables/list_effector_genes_coloc_mr_exonic_with_coloc_rank1_noFM_with_otar_evidence_with_coloc_direction.tsv.gz"))
list_genes_tmp<-as.data.frame(list_genes)
list_genes_tmp<-list_genes_tmp[,c("list_genes","negative","positive")]

# LOAD THE LIST OF EFFECTOR GENES :
list_genes<-fread(paste0(path_gwas,"post_imputation/2022/analysis/final_tables/list_effector_genes_coloc_mr_exonic_with_coloc_rank1_noFM_with_otar_evidence.tsv.gz"))
list_genes<-as.data.frame(list_genes)

list_genes<-merge(list_genes,list_genes_tmp,by="list_genes")


list_genes<-merge(list_genes,mr_genes,by="list_genes",all.x=T)

# col<-fread("/path/to/user/git/IIBDGC_GWAS/plots/paper_tables/Source_data_coloc_results_beta_wald.tsv.gz",head=T)
# list_genes<-merge(list_genes,col,by="list_genes",all.x=T)
colnames(list_genes)

list_genes$beta_gsmr<-list_genes$direction_effect_gsmr

list_genes$direction_effect_gsmr<-NA
list_genes$direction_effect_gsmr[which(list_genes$beta_gsmr>0)]<-"same"
list_genes$direction_effect_gsmr[which(list_genes$beta_gsmr<0)]<-"opposite"


list_genes$direction_effect_colocalization<-NA
list_genes$direction_effect_colocalization[which(list_genes$negative>0)]<-"opposite"
list_genes$direction_effect_colocalization[which(list_genes$positive>0)]<-"same"
list_genes$direction_effect_colocalization[which(list_genes$positive>0 & list_genes$negative>0)]<-"opposite;same"

table(list_genes$direction_effect_gsmr)
# opposite     same 
#        2       14

table(list_genes$direction_effect_colocalization)
    #  opposite opposite;same          same 
    #       294            78           253 

fwrite(list_genes,paste0(path_gwas,"post_imputation/2022/analysis/final_tables/list_effector_genes_coloc_mr_exonic_with_coloc_rank1_noFM_with_otar_evidence_with_coloc_mr_direction.tsv.gz"),
col.names=T,row.names=F,quote=F,sep="\t")

q("no")
