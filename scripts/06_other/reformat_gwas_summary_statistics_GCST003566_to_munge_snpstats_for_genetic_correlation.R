# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
######################################################################################################################################################################
######################################################################################################################################################################
### MANUALLY ADD - partially harmonised by GWAS catalog:

# 5. Successfully harmonised variants

# 97.78% ( 7790531 of 7967643 ) sites successfully harmonised.
# hm_code	Number	Percentage	Explanation
# 10	1372090	17.22%	Forward strand; Correct orientation; Already harmonised
# 11	5285453	66.34%	Forward strand; Flipped orientation; Requires harmonisation
# 12	1488	0.02%	Reverse strand; Correct orientation; Already harmonised
# 13	4252	0.05%	Reverse strand; Flipped orientation; Requires harmonisation
# 5	233936	2.94%	Palindromic; Assume forward strand; Correct orientation; Already harmonised
# 6	893312	11.21%	Palindromic; Assume forward strand; Flipped orientation; Requires harmonisation

# path_gwas="/path/to/ibdgwas/IIBDGC/"
path_gwas="/path/to/ibdgwas/IIBDGC/"

cd /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/

# download clonal hematopoyesis:
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003566/harmonised/GCST003566.h.tsv.gz


# singularity exec iibdgc_postprocess_10_singularity.sif
# singularity exec polyfun_2_singularity.sif

# harmonise nomenclature based on the instructructions provided by GWAS catalog

# MEM=45000
# bsub -Is -M"$MEM" -R"select[model==Intel_Platinum && mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group R 

library(data.table)
library(R.utils)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

df<-fread(paste(path_gwas,"resources/gwas_summary_statistics/GCST003566.h.tsv.gz",sep=""),)
dim(df)
#[1] 61937846 6

df<-df[which(df$effect_allele_frequency>=0.001 & df$effect_allele_frequency<=0.999),]
dim(df)
# [1] 7790531       22


df$SNP<-paste("chr",df$chromosome,":",df$base_pair_location,":",df$other_allele,":",df$effect_allele,sep="")
df<-df[which(df$effect_allele!=df$other_allele),]
dim(df)
#[1] 15033621       10 

# hm_code	Number	Percentage	Explanation
# 10	1372090	17.22%	Forward strand; Correct orientation; Already harmonised
# 11	5285453	66.34%	Forward strand; Flipped orientation; Requires harmonisation
# 12	1488	0.02%	Reverse strand; Correct orientation; Already harmonised
# 13	4252	0.05%	Reverse strand; Flipped orientation; Requires harmonisation
# 5	233936	2.94%	Palindromic; Assume forward strand; Correct orientation; Already harmonised
# 6	893312	11.21%	Palindromic; Assume forward strand; Flipped orientation; Requires harmonisation

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

summary(df$BETA)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -1.0888740 -0.0281953 -0.0005142 -0.0002690  0.0272304  1.2572400 

summary(df$effect_allele_frequency)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.009501 0.565906 0.819736 0.725024 0.942434 0.990525

summary(df$A1FREQ) # MAF distribution?
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.009475 0.053930 0.152302 0.188392 0.307683 0.505662 

df$SE<-df$standard_error

df$N<-15283+6226
df$INFO<-1


df<-df[,c("SNP","chromosome","base_pair_location","other_allele","effect_allele","BETA","SE","p_value","INFO","N","A1FREQ")]

colnames(df)<-c("SNP","CHR","BP","A0","A1","BETA","SE","PVALUE","INFO","N","A1FREQ")


# remove duplicated ids:
dups<-df$SNP[duplicated(df$SNP)]
length(dups)
# [1] 29
df<-df[!df$SNP %in% dups,]
dim(df)
# [1] 6663225       11

fwrite(df,paste(path_gwas,"resources/gwas_summary_statistics/GCST003566_edited.tsv.gz",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

q("no")

#########################################


gwas_id=GCST003566

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

gwas_id=GCST003566
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
