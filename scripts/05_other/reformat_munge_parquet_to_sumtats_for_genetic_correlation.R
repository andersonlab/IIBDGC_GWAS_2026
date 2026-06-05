# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#


# MEM=3000
# for i in 0
# do
# bsub -J"subset_gwas" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
# -o ${path_gwas}post_imputation/2022/log/munge_parquet_to_sumstats_${gwas_id[i]}_stdout \
# -e ${path_gwas}post_imputation/2022/log/munge_parquet_to_sumstats_${gwas_id[i]}_stderr \
# "Rscript ~/git/IIBDGC_GWAS/scripts/other/reformat_munge_parquet_to_sumtats_for_genetic_correlation.R ${gwas_id[i]} > \
# ${path_gwas}post_imputation/2022/log/reformat_munge_parquet_to_sumtats_for_genetic_correlation_${gwas_id[i]}.Rout"
# done

# MEM=3000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q yesterday R

library(data.table)
library(arrow)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

args = commandArgs(trailingOnly=TRUE)
study_id<-args[1]

# to test 
# study_id<-"GCST90165267"

print(study_id)

tmp<-read_parquet(paste(path_gwas,"resources/gwas_summary_statistics/",study_id,"_only_SNPs_sumstats_munged.parquet",sep=""))
tmp<-as.data.frame(tmp)

print(head(tmp))
print(nrow(tmp))

tmp<-tmp[,c("SNP","CHR","BP","A2","A1","N","MAF","Z")]

fwrite(tmp,paste(path_gwas,"resources/gwas_summary_statistics/",study_id,"_only_SNPs_sumstats_munged.sumstats",sep=""),
col.names=T,row.names=F,quote=F,sep="\t")

q("no")