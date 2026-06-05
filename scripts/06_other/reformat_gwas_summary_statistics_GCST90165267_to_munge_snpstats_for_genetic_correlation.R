# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
######################################################################################################################################################################
######################################################################################################################################################################
### MANUALLY ADD - # download clonal hematopoyesis:


path_gwas="/path/to/ibdgwas/IIBDGC/"
gwas_id="GCST90165267"

cd ${path_gwas}resources/gwas_summary_statistics/

# download clonal hematopoyesis:
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90165001-GCST90166000/GCST90165267/harmonised/GCST90165267.h.tsv.gz

# singularity exec iibdgc_postprocess_10_singularity.sif
# singularity exec polyfun_2_singularity.sif


# harmonise nomenclature:

# MEM=45000
# bsub -Is -M"$MEM" -R"select[model==Intel_Platinum && mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group R \

library(data.table)
library(R.utils)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"
gwas_id<-"GCST90165267"

df<-fread(paste0(path_gwas,"resources/gwas_summary_statistics/",gwas_id,".h.tsv.gz"),)
dim(df)
#[1] 56831521 6

df<-df[which(df$effect_allele_frequency>=0.001 & df$effect_allele_frequency<=0.999),]
dim(df)
# [1] 14115013       22


df$SNP<-paste("chr",df$chromosome,":",df$base_pair_location,":",df$other_allele,":",df$effect_allele,sep="")
df<-df[which(df$effect_allele!=df$other_allele),]
dim(df)
#[1] 14115013       10 

# Successfully harmonised variant
# 95.52% ( 56831521 of 59494575 ) sites successfully harmonised.
# hm_code	Number	Percentage	Explanation
# 10	49017552	82.39%	Forward strand; Correct orientation; Already harmonised
# 11	17193	0.03%	Forward strand; Flipped orientation; Requires harmonisation
# 12	7	0.00%	Reverse strand; Correct orientation; Already harmonised
# 13	2	0.00%	Reverse strand; Flipped orientation; Requires harmonisation
# 5	7796767	13.11%	Palindromic; Assume forward strand; Correct orientation; Already harmonised

df$beta<-log(df$odds_ratio)

df$BETA<-NA
df$A1FREQ<-NA

df$BETA[which(df$hm_code %in% c(10,12))]<-df$beta[which(df$hm_code %in% c(10,12))]
df$A1FREQ[which(df$hm_code %in% c(10,12))]<-df$effect_allele_frequency[which(df$hm_code %in% c(10,12))]

df$BETA[which(df$hm_code %in% c(11,13))]<-df$beta[which(df$hm_code %in% c(11,13))]*(-1)
df$A1FREQ[which(df$hm_code %in% c(11,13))]<-1-(df$effect_allele_frequency[which(df$hm_code %in% c(11,13))])

# exclude the unclear:
df<-df[which(!df$hm_code %in% c(5,6)),]
dim(df)
# [1] 6663283      15

df$SE<-(log(df$ci_upper)-log(df$ci_lower))/3.919928

df$N<-368526
df$INFO<-1

summary(df$BETA)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -36.50842  -0.01854  -0.00003   0.00010   0.01830  28.11599

summary(df$effect_allele_frequency)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.00100 0.00396 0.03107 0.16928 0.23674 0.99900

summary(df$A1FREQ)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.00100 0.00396 0.03107 0.16931 0.23679 0.99900


# length(df$beta_ed[which(df$other_allele==df$A1 & df$effect_allele==df$A0)])

df<-df[,c("SNP","chromosome","base_pair_location","other_allele","effect_allele","BETA","SE","p_value","INFO","N","A1FREQ")]

colnames(df)<-c("SNP","CHR","BP","A0","A1","BETA","SE","PVALUE","INFO","N","A1FREQ")
dim(df)
# [1] 12055364       11

# remove duplicated ids:
dups<-df$SNP[duplicated(df$SNP)]
length(dups)
# [1] 1665
df<-df[!df$SNP %in% dups,]
dim(df)
# [1] 12052034       11

fwrite(df,paste(path_gwas,"resources/gwas_summary_statistics/",gwas_id,"_edited.tsv.gz",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

q("no")

#########################################


MEM=8200

for i in 0
do
bsub -J"subset_gwas" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-o ${path_gwas}post_imputation/2022/log/subset_gwas_variants_${gwas_id[i]}_stdout \
-e ${path_gwas}post_imputation/2022/log/subset_gwas_variants_${gwas_id[i]}_stderr \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/subset_gwas_variants_for_genetic_correlation.R ${gwas_id[i]} > \
${path_gwas}post_imputation/2022/log/subset_gwas_variants_for_genetic_correlation_${gwas_id[i]}.Rout"
done

for i in  0
do
echo ${i} && tail -50 ${path_gwas}post_imputation/2022/log/subset_gwas_variants_${gwas_id[i]}_stdout  | grep -E "Successfully|Exit"
done

for i in 0
do
echo ${i} && ls -la ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_edited_2.tsv.gz
done

for i in 0
do
cat ${path_gwas}post_imputation/2022/log/subset_gwas_variants_for_genetic_correlation_${gwas_id[i]}.Rout
done

#######################################################################################################
# reformat to munged file

MEM=7000

for i in 0
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${gwas_id[i]}_ldsc_sumstats_only_SNPs_polyfun_mungen_format_stderr \
-o ${path_gwas}post_imputation/log/${gwas_id[i]}_ldsc_sumstats_only_SNPs_polyfun_mungen_format_stdout \
"munge_polyfun_sumstats.py \
--sumstats ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_edited_2.tsv.gz \
--min-info 0 \
--min-maf 0 \
--remove-strand-ambig \
--out ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_only_SNPs_sumstats_munged.parquet"
done

for i in 0
do
echo ${i} && tail -50 ${path_gwas}post_imputation/log/${gwas_id[i]}_ldsc_sumstats_only_SNPs_polyfun_mungen_format_stdout  | grep -E "Successfully|Exit"
done

for i in 0
do
echo ${i} && ls -la ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_only_SNPs_sumstats_munged.parquet
done


# SAVE PARQUET FILE AS FLAT TEXT FILE AS INPUT FOR LDSC/LDSC.PY

module unload HGI/softpack/groups/team152/iibdgc_postprocess/10
module unload HGI/softpack/groups/team152/polyfun-2/1

MEM=3000
for i in 0
do
bsub -J"subset_gwas" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-o ${path_gwas}post_imputation/2022/log/munge_parquet_to_sumstats_${gwas_id[i]}_stdout \
-e ${path_gwas}post_imputation/2022/log/munge_parquet_to_sumstats_${gwas_id[i]}_stderr \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/reformat_munge_parquet_to_sumtats_for_genetic_correlation.R ${gwas_id[i]} > \
${path_gwas}post_imputation/2022/log/reformat_munge_parquet_to_sumtats_for_genetic_correlation_${gwas_id[i]}.Rout"
done

for i in 0
do
echo ${i} && tail -50 ${path_gwas}post_imputation/2022/log/munge_parquet_to_sumstats_${gwas_id[i]}_stdout | grep -E "Successfully|Exit"
done

for i in 0
do
echo ${i} && ls -la ${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_only_SNPs_sumstats_munged.sumstats
done


# COMPLETED - JUNE 25
