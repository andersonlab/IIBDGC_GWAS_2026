# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# See: https://github.com/josefin-werme/LAVA

# local genetic correlation will allow us later to pinpoint shared effector genes between traits

# singularity exec iibdgc_postprocess_10_singularity.sif
path_gwas="/path/to/ibdgwas/IIBDGC/"

######################################################################################################
## do not restrict to the list of significant global rg

/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/genetic_correlation/list_traits_significantly_correlated_with_ibd_cd_uc

gwas_id=($(cat /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/genetic_correlation/list_traits_significantly_correlated_with_ibd_cd_uc)) 
echo ${#gwas_id[@]}

for i in {0..86}
do
echo ${gwas_id[i]}
done


######################################################################################
# 1 - CREATE AN INFO FILE REQUIRED TO HANDLE MULTIPLE PHENOTYPES IN THE ANALYSIS:

# Input info file, used for convenient processing of multiple phenotypes. Requires the columns:

# phenotype: phenotype IDs
# cases: number of cases (set to NA for continuous phenotypes)
# controls: number of controls (set to NA for continuous phenotypes)
# prevalence (optional): the population prevalence of binary phenotypes
# this is only relevant if you want an estimate of the local population h2 for binary phenotypes. Estimates of the observed local sample h2 are still provided
# filename: paths and file names to the relevant summary statistics

MEM=800
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q normal R 

library(data.table)
library(R.utils)
library(stringr)

path_gwas="/path/to/ibdgwas/IIBDGC/"

ids<-fread("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/genetic_correlation/list_traits_significantly_correlated_with_ibd_cd_uc",head=F)
colnames(ids)<-"phenotype"

# retrieve N cases and controls:
all<-fread(paste(path_gwas,"resources/gwas_summary_statistics/gwas-catalog-v1.0.3.1-studies-r2025-05-13.tsv.gz",sep=""),quote="")
all<-all[which(all$'STUDY ACCESSION' %in% ids$phenotype),]


# save sample size:
all$N_tmp<-gsub("[A-z]*","",all$'INITIAL SAMPLE SIZE')
all$N_tmp<-gsub(",","",all$N_tmp)

tmp<-str_split(all$N_tmp," ")


all$cases<-NA
all$controls<-NA
all$N<-NA

for(i in c(1:nrow(all))) {

    all$N[i]<-sum(as.numeric(tmp[[i]]),na.rm=T)

    all$cases[i]<-sum(as.numeric(tmp[[i]][1]),as.numeric(tmp[[i]][9]),na.rm=T)
    all$controls[i]<-sum(as.numeric(tmp[[i]][5]),as.numeric(tmp[[i]][14]),na.rm=T)

}

# all cuantitative:
all[which(all$controls==0),"MAPPED_TRAIT"]

# binomial:
all[which(all$controls!=0),"MAPPED_TRAIT"]

# set to missing Ncases and Ncontrols for cuantitative traits:
all$cases[which(all$controls==0)]<-NA
all$controls[which(all$controls==0)]<-NA

all<-all[,c("STUDY ACCESSION","MAPPED_TRAIT","cases","controls")]


ids[which(!ids$phenotype %in% all$'STUDY ACCESSION'),]
#    phenotype
#       <char>
# 1:  Hyper_He
# 2:   IRF_FSC
# 3:    MicroR
# 4:   RBC_FSC
# 5:    RBC_He
# 6:   RBC_SSC
# 7:   RET_FSC
# 8:    RET_He

# all interval, all cuantitative traits:
ids<-merge(ids,all[,c("STUDY ACCESSION","cases","controls")],by.x="phenotype",by.y="STUDY ACCESSION",all.x=T)

# add IBD phenotypes - for EUR_tier2 use Neff, and 50% split in cases and controls:
ibd<-as.data.frame(matrix(c("ibd","cd","uc","120193","59461","72549","120193","59461","72549"),nrow=3,ncol=3))
colnames(ibd)<-colnames(ids)
ids<-rbind(ids,ibd)

ids$filename<-paste0(path_gwas,"resources/gwas_summary_statistics/",ids$phenotype,"_only_SNPs_sumstats_munged.sumstats")
ids$filename[which(ids$phenotype %in% c("ibd","cd","uc"))]<-paste0(path_gwas,"resources/gwas_summary_statistics/",ids$phenotype[which(ids$phenotype %in% c("ibd","cd","uc"))],"_eur_tier1_list_variants_only_SNPs_rate_0.5_nohla_sumstats_munged.sumstats")


write.table(ids,paste0(path_gwas,"post_imputation/2022/analysis/local_genetic_correlation/input_files/input_ibd_cd_uc_local_genetic_correlation_list_traits_significantly_correlated.info.txt"),
col.names=T,row.names=F,quote=F,sep="\t")

q("no")

######################################################################################
# 2 - CREATE A REGIONS FILE:

MEM=800
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q normal R 

library(data.table)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

reg<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_new_and_old_index_gwas_variants_conditional_on_known_index_gwas_fm_2023_with_iibdgc_tier1_tier2_eas_sas_summary_stats_plus_unsupervised_conditional.tsv",sep=""))

reg<-reg[which(reg$updated_region!=""),c("updated_region"),]
reg<-reg[!duplicated(reg),]

write.table(reg,paste0(path_gwas,"post_imputation/2022/analysis/local_genetic_correlation/input_files/regions.txt"),
col.names=F,row.names=F,quote=F,sep="\t")


#####################################


path_gwas="/path/to/ibdgwas/IIBDGC/"

regions=($(cat /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/local_genetic_correlation/input_files/regions.txt)) 
echo ${#regions[@]}
# 421

for i in {0..420}
do
echo ${i} && echo ${regions[i]}
done


MEM=50000

for i in {0..420}
for i in 5 108 225 230 231 239 261 318 338 372 415
do
bsub -J"lava" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q week \
-o ${path_gwas}post_imputation/2022/log/lava_${regions[i]}_stdout \
-e ${path_gwas}post_imputation/2022/log/lava_${regions[i]}_stderr \
"Rscript /path/to/user/git/IIBDGC_GWAS/scripts/other/run_lava.R ${regions[i]} > \
${path_gwas}post_imputation/2022/log/run_lava_${regions[i]}.Rout"
done

# COMPLETED

for i in {0..420}
do
echo ${i} && tail -50  ${path_gwas}post_imputation/2022/log/lava_${regions[i]}_stdout | grep -E "after reaching LSF memory usage limit"
done


for i in {0..420}
do
echo ${i} && tail -50  ${path_gwas}post_imputation/2022/log/lava_${regions[i]}_stdout | grep -E "job killed after reaching LSF run time limit."
done

for i in {0..420}
do
echo ${i} && tail -50  ${path_gwas}post_imputation/2022/log/lava_${regions[i]}_stdout | grep -E "Successfully|Exited"
done

for i in {0..420}
do
echo ${i} && cat ${path_gwas}post_imputation/2022/log/run_lava_${regions[i]}.Rout
done.


## plot results and create data release:
~/git/IIBDGC_GWAS/scripts/other/plot_local_genetic_correlation_lava.R

 