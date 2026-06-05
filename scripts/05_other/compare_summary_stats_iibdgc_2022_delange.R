# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# get summary stats from Liu Liu Gao 2017 for EUR (= delange + FinnGen) with our latest IIBDGC summary stats:


# singularity exec iibdgc_postprocess_10_singularity.sif

path_gwas="/path/to/ibdgwas/IIBDGC/"


MEM=15000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q normal R 

library(data.table)

path_gwas="/path/to/ibdgwas/IIBDGC/"

file<-list.files("/path/to/project")
file<-file[which(!file %in% "README")]
pheno<-gsub("_.*","",file)

rm(dat)

for (i in 1:length(file)) {

    print(pheno[i])
    dat<-fread(paste0("/path/to/project",file[i]))
    dat$pos<-gsub("[0-9]{1,2}:","",dat$MarkerName)
    dat$pos<-as.numeric(gsub("_[A-Z]*_[A-Z]*","",dat$pos))
    dat$end<-dat$pos+1
    dat<-dat[order(dat$pos,decreasing=F),]

    dat$chr<-as.numeric(gsub(":.*","",dat$MarkerName))
    dat<-dat[order(dat$chr,decreasing=F),]
    dat$chr<-paste0("chr",dat$chr)

    dat$chisq<-(dat$Effect/dat$StdErr)^2

    dat<-dat[,c("chr","pos","end","MarkerName","P.value","chisq")]

    dat$pos<-format(dat$pos, scientific=FALSE)
    dat$end<-format(dat$end, scientific=FALSE)

    fwrite(dat,paste0(path_gwas,"summary_files/summary_stats_",pheno[i],"_delange_sorted_noheader_chr_pos_b37_pval"),col.names=F,row.names=F,sep="\t")
    rm(dat)

}

q("no")

##############################

pheno=(cd ibd uc)
MEM=200

for ph in ${pheno[@]}
do
bsub -J"lift" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-o ${path_gwas}post_imputation/2022/log/liftover_${ph}_delange_stdout \
-e ${path_gwas}post_imputation/2022/log/liftover_${ph}_delange_stderr \
"${path_gwas}previous_qced_b38/liftover/liftOver \
${path_gwas}summary_files/summary_stats_${ph}_delange_sorted_noheader_chr_pos_b37_pval \
${path_gwas}previous_qced_b38/liftover/hg19ToHg38.over.chain.gz \
${path_gwas}summary_files/summary_stats_${ph}_delange_sorted_noheader_chr_pos_b37_pval_lifted_hg38 \
${path_gwas}summary_files/summary_stats_${ph}_delange_sorted_noheader_chr_pos_b37_pval_lifted_no_lifted_hg38"
done

for ph in ${pheno[@]}
do
echo ${ph} && cat ${path_gwas}summary_files/summary_stats_${ph}_delange_sorted_noheader_chr_pos_b37_pval_lifted_no_lifted_hg38 | sed '/^#/d' | wc -l
done
# cd
# 1180
# ibd
# 1193
# uc
# 1189


for ph in ${pheno[@]}
do
wc -l ${path_gwas}summary_files/summary_stats_${ph}_delange_sorted_noheader_chr_pos_b37_pval_lifted_hg38
done
# 9569607 /path/to/ibdgwas/IIBDGC/summary_files/summary_stats_cd_delange_sorted_noheader_chr_pos_b37_pval_lifted_hg38
# 9734253 /path/to/ibdgwas/IIBDGC/summary_files/summary_stats_ibd_delange_sorted_noheader_chr_pos_b37_pval_lifted_hg38
# 9586827 /path/to/ibdgwas/IIBDGC/summary_files/summary_stats_uc_delange_sorted_noheader_chr_pos_b37_pval_lifted_hg38

# compress data:
for ph in ${pheno[@]}
do
gzip ${path_gwas}summary_files/summary_stats_${ph}_delange_sorted_noheader_chr_pos_b37_pval_lifted_hg38
done


#########################

pheno=(cd ibd uc)
MEM=15000

release=(eur_eas_sas_tier_2 eur_tier_2 eur_tier_1)
i=2

# for i in {0..2}
# do
for ph in ${pheno[@]}
do
bsub -J"lift" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-o ${path_gwas}post_imputation/2022/log/plot_${ph}_delange_vs_iibdgc_${release[i]}_pval_stdout \
-e ${path_gwas}post_imputation/2022/log/plot_${ph}_delange_vs_iibdgc_${release[i]}_pval_stderr \
"Rscript /path/to/user/git/IIBDGC_GWAS/scripts/other/plot_comparison_iibdgc_2022_delange.R ${ph} ${release[i]} > \
${path_gwas}post_imputation/2022/log/plot_comparison_iibdgc_2022_delange_${ph}_${release[i]}.Rout"
done
# done


#########################

pheno=(cd ibd uc)
MEM=15000

release=(eur_eas_sas_tier_2 eur_tier_2 eur_tier_1)
i=1
echo ${release[i]}

# for i in {0..2}
# do
for ph in ${pheno[@]}
do
bsub -J"lift" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-o ${path_gwas}post_imputation/2022/log/plot_${ph}_delange_vs_iibdgc_${release[i]}_pval_stdout \
-e ${path_gwas}post_imputation/2022/log/plot_${ph}_delange_vs_iibdgc_${release[i]}_pval_stderr \
"Rscript /path/to/user/git/IIBDGC_GWAS/scripts/other/plot_comparison_iibdgc_2022_delange.R ${ph} ${release[i]} > \
${path_gwas}post_imputation/2022/log/plot_comparison_iibdgc_2022_delange_${ph}_${release[i]}.Rout"
done
# done


