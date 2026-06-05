# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
############################################################################################################################################################
# RUN THE SAME ANALYSIS JUST INCLUDING HAPMAP3 VARIANTS:

# create the reference file to subset and compare:

# singularity exec iibdgc_postprocess_10_singularity.sif
path_gwas="/path/to/ibdgwas/IIBDGC/"

MEM=4000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q yesterday R

library(data.table)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

list_hapmap_variants<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/w_hm3.snplist",sep=""))

for (chr in c(1:22)) {

    tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/GRCh38/plink_files/1000G.EUR.hg38.",chr,".bim",sep=""))
    tmp<-tmp[which(tmp$V2 %in% list_hapmap_variants$SNP),]

    if(chr==1) {
        all<-tmp
    } else {
        all<-rbind(all,tmp)
    }
    rm(tmp)
}

# double check that Alele notation in bim does not reflect frequency
dim(all)
# [1] 1189841       6

all$SNP<-paste("chr",all$V1,":",all$V4,":",all$V6,":",all$V5,sep="")
all$X<-paste("chr",all$V1,":",all$V4,sep="")

list_snps<-fread(paste(path_gwas,"resources/gwas_summary_statistics/subset_gwas_variants_for_genetic_correlation.snplist",sep=""),head=T)
list_snps$X<-gsub(":[A-Z]{1}:[A-Z]{1}","",list_snps$SNP)

dim(all[which(all$X %in% list_snps$X),])
# [1] 1185395       8

dim(all[which(all$SNP %in% list_snps$SNP),])
# [1] 833577      7

# subset the original list to the 1185395 variants present in both
dim(list_snps)
# [1] 14601777        4
list_snps<-list_snps[which(list_snps$X %in% all$X),]
dim(list_snps)
# [1] 1185418       4

fwrite(list_snps[,1:3],paste(path_gwas,"resources/gwas_summary_statistics/subset_gwas_variants_hm3_for_genetic_correlation.snplist",sep=""),
col.names=T,row.names=F,quote=F,sep="\t")

rm(all)

# subset the summary statistics so they only those SNPs

path_gwas<-"/path/to/ibdgwas/IIBDGC/"
pheno<-c("ibd","cd","uc")


for (ph in pheno) {

    print(ph)
    tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/allchr_",ph,"_eur_tier1_list_variants_rate_0.5_with_pval.tsv.gz",sep=""),head=T)
    tmp<-tmp[which(tmp$SNP %in% list_snps$SNP),]

    tmp<-tmp[,c("SNP","CHR","BP","ALLELE2","ALLELE1","BETA","SE","P-value","INFO","N","A1FREQ")]

    fwrite(tmp,paste(path_gwas,"resources/gwas_summary_statistics/allchr_",ph,"_eur_tier1_list_variants_only_SNPs_hm3_rate_0.5_with_pval.tsv.gz",sep=""),
    col.names=T,row.names=F,quote=F,sep="\t")
    rm(tmp)

}

q("no")


###################################################
# REFORMAT IBD SUMMARY STATISTICS - rg from polyfun/ldsc does not work use a combination of munge_polyfun_sumstats -> parquet to sumstats -> ldsc -rg (from ldsc)

# module unload HGI/softpack/groups/team152/iibdgc_postprocess/10
path_gwas="/path/to/ibdgwas/IIBDGC/"

# version for polyfun/munge
# singularity exec polyfun_2_singularity.sif

MEM=7000

pheno=(ibd cd uc)

for ph in ${pheno[@]}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${ph}_allchr_ldsc_sumstats_only_SNPs_hm3_nohla_polyfun_mungen_format_stderr \
-o ${path_gwas}post_imputation/log/${ph}_allchr_ldsc_sumstats_only_SNPs_hm3_nohla_polyfun_mungen_format_stdout \
"munge_polyfun_sumstats.py \
--sumstats ${path_gwas}resources/gwas_summary_statistics/allchr_${ph}_eur_tier1_list_variants_only_SNPs_hm3_rate_0.5_with_pval.tsv.gz \
--min-info 0 \
--min-maf 0 \
--remove-strand-ambig \
--out ${path_gwas}resources/gwas_summary_statistics/${ph}_eur_tier1_list_variants_only_SNPs_hm3_rate_0.5_nohla_sumstats_munged.parquet"
done

less ${path_gwas}post_imputation/log/ibd_allchr_ldsc_sumstats_only_SNPs_hm3_nohla_polyfun_mungen_format_stdout 
# [INFO]  Reading sumstats file...
# [INFO]  Done in 3.71 seconds
# [INFO]  1185283 SNPs are in the sumstats file
# [INFO]  Removing 4 SNPs with strand ambiguity
# [INFO]  Removing 1371 HLA SNPs
# [INFO]  1183908 SNPs with sumstats remained after all filtering stages
# [INFO]  Saving munged sumstats of 1183908 SNPs to /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/ibd_eur_tier1_list_var
# iants_only_SNPs_hm3_rate_0.5_nohla_sumstats_munged.parquet
# [INFO]  Done

less ${path_gwas}post_imputation/log/cd_allchr_ldsc_sumstats_only_SNPs_hm3_nohla_polyfun_mungen_format_stdout 
# [INFO]  Reading sumstats file...
# [INFO]  Done in 3.75 seconds
# [INFO]  1185272 SNPs are in the sumstats file
# [INFO]  Removing 4 SNPs with strand ambiguity
# [INFO]  Removing 1371 HLA SNPs
# [INFO]  1183897 SNPs with sumstats remained after all filtering stages
# [INFO]  Saving munged sumstats of 1183897 SNPs to /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/cd_eur_tier1_list_variants_only_SNPs_hm3_rate_0.5_nohla_sumstats_munged.parquet
# [INFO]  Done

less ${path_gwas}post_imputation/log/uc_allchr_ldsc_sumstats_only_SNPs_hm3_nohla_polyfun_mungen_format_stdout 
# [INFO]  Reading sumstats file...
# [INFO]  Done in 2.43 seconds
# [INFO]  1185217 SNPs are in the sumstats file
# [INFO]  Removing 3 SNPs with strand ambiguity
# [INFO]  Removing 1371 HLA SNPs
# [INFO]  1183843 SNPs with sumstats remained after all filtering stages
# [INFO]  Saving munged sumstats of 1183843 SNPs to /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/uc_eur_tier1_list_vari
# ants_only_SNPs_hm3_rate_0.5_nohla_sumstats_munged.parquet
# [INFO]  Done


module unload HGI/softpack/groups/team152/iibdgc_postprocess/10

MEM=6000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q yesterday R

library(arrow)
library(data.table)

pheno<-c("ibd","cd","uc")

for (ph in pheno) {

    print(ph)

    df<-read_parquet(paste("/path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/",ph,"_eur_tier1_list_variants_only_SNPs_hm3_rate_0.5_nohla_sumstats_munged.parquet",sep=""))
    print(head(df))

    print(dim(df))

    fwrite(df[,2:9],paste("/path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/",ph,"_eur_tier1_list_variants_only_SNPs_hm3_rate_0.5_nohla_sumstats_munged.sumstats",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    rm(df)

}

# [1] "ibd"
#   index                SNP CHR        BP A2 A1      N        MAF          Z
# 1     0  chr1:15494905:C:A   1  15494905  C  A 102561 0.29822650  0.6637168
# 2     1 chr1:117591017:C:T   1 117591017  C  T 102561 0.08735635  0.7580645
# 3     2  chr1:89264267:C:T   1  89264267  C  T 102561 0.43524176 -0.5045872
# 4     3 chr1:115378803:T:G   1 115378803  T  G 102561 0.38727522  3.1944444
# 5     4 chr1:225400227:G:A   1 225400227  G  A 102561 0.48124888  0.8761905
# 6     5   chr1:2613201:G:A   1   2613201  G  A 102561 0.31201597 -3.1071429
# [1] 1183908       9
# [1] "cd"
#   index                SNP CHR        BP A2 A1     N        MAF           Z
# 1     0  chr1:15494905:C:A   1  15494905  C  A 72834 0.29660994 -0.16058394
# 2     1 chr1:117591017:C:T   1 117591017  C  T 72834 0.08791637  0.01345291
# 3     2  chr1:89264267:C:T   1  89264267  C  T 72834 0.43443696 -0.12781955
# 4     3 chr1:115378803:T:G   1 115378803  T  G 72834 0.38715478  3.64885496
# 5     4 chr1:225400227:G:A   1 225400227  G  A 72834 0.48107217  0.68503937
# 6     5   chr1:2613201:G:A   1   2613201  G  A 72834 0.31705557  0.78518519
# [1] 1183897       9
# [1] "uc"
#   index                SNP CHR        BP A2 A1     N        MAF          Z
# 1     0 chr1:117591017:C:T   1 117591017  C  T 64380 0.08648040  0.9957265
# 2     1 chr1:225400227:G:A   1 225400227  G  A 64380 0.47744037  0.9848485
# 3     2   chr1:2613201:G:A   1   2613201  G  A 64380 0.31042243 -5.3829787
# 4     3 chr1:114717147:T:C   1 114717147  T  C 64380 0.07204753 -0.2371542
# 5     4 chr1:246967275:C:A   1 246967275  C  A 64380 0.33354908  1.4927536
# 6     5 chr1:241229625:A:C   1 241229625  A  C 64380 0.14164343  1.6719577
# [1] 1183843       9

q("no")



##########################################################################################################################################################
##########################################################################################################################################################

files=($(cat /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/list_path_summary_stats_gwas_catalog)) 
echo ${#files[@]}
# 269

gwas_id=($(cat /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/list_path_summary_stats_gwas_catalog | \
sed 's/\/lustre\/scratch123\/hgi\/projects\/ibdgwas\/IIBDGC\/resources\/gwas_summary_statistics\///g' | \
sed 's/\/lustre\/scratch125\/humgen\/resources\/GWAScatalog\/.*\///g' | sed  's/.h.tsv.gz//g')) 
echo ${#gwas_id[@]}
# 269


module unload HGI/softpack/groups/team152/polyfun-2/1
# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=7500

pheno=(ibd cd uc)

i=259
echo ${gwas_id[i]}
# GCST004030 - PSC Ji study

#######################################################################
# test using 1000GP weigths - rename variants to rsid

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=25000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q yesterday R

library(data.table)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"
pheno<-c("ibd","cd","uc")

list_hapmap_variants<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/w_hm3.snplist",sep=""))

# list of IDs on UKB ldsc
ukbldsc<-fread("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/UKBB.ALL.ldscore/UKBB.EUR.rsid.l2.ldscore.gz")

for (chr in c(1:22)) {

    tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/GRCh38/plink_files/1000G.EUR.hg38.",chr,".bim",sep=""))
    tmp<-tmp[which(tmp$V2 %in% list_hapmap_variants$SNP),]

    if(chr==1) {
        all<-tmp
    } else {
        all<-rbind(all,tmp)
    }
    rm(tmp)
    
}

all$X<-paste("chr",all$V1,":",all$V4,sep="")

# update PSC data - check that the alleles match

psc<-fread(paste(path_gwas,"resources/gwas_summary_statistics/GCST004030_only_SNPs_sumstats_munged.sumstats",sep=""))
psc$X<-gsub(":[A-Z]{1}:[A-Z]{1}","",psc$SNP)


for (ph in pheno) {

    print(ph)

    df<-fread(paste("/path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/",ph,"_eur_tier1_list_variants_only_SNPs_hm3_rate_0.5_nohla_sumstats_munged.sumstats",sep=""),)
    print(head(df))

    df$X<-gsub(":[A-Z]{1}:[A-Z]{1}","",df$SNP)
    df<-merge(df,all[,c("V2","X")],by="X")

    df<-df[which(df$SNP %in% psc$SNP),]
    psc<-psc[which(psc$SNP %in% df$SNP),]

    df1<-df[,c(10,3:9)]
    colnames(df1)[1]<-"SNP"

    fwrite(df1,paste("/path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/",ph,"_eur_tier1_list_variants_only_SNPs_hm3_rsids_rate_0.5_nohla_sumstats_munged.sumstats",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    rm(df1)

    df1<-df[,c(2,3:9)]

    fwrite(df1,paste("/path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/",ph,"_eur_tier1_list_variants_only_SNPs_hm3_chrposid_rate_0.5_nohla_sumstats_munged.sumstats",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    
    df1$SNP<-gsub("chr","",df1$SNP)
    
    fwrite(df1,paste("/path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/",ph,"_eur_tier1_list_variants_only_SNPs_hm3_chrposid_ukbid_rate_0.5_nohla_sumstats_munged.sumstats",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    
    rm(df1)

    # rm(df)

}

dim(df)
dim(psc)

psc<-merge(psc,all[,c("V2","X")],by="X")
table(psc$V2==df$SNP)
table(psc$A2==df$A2)
table(psc$A1==df$A1)

psc1<-psc[,c(10,3:9)]
colnames(psc1)[1]<-"SNP"

fwrite(psc1,paste(path_gwas,"resources/gwas_summary_statistics/GCST004030_only_SNPs_hm3_rsids_sumstats_munged.sumstats",sep=""),
col.names=T,row.names=F,quote=F,sep="\t")
rm(psc1)

psc1<-psc[,c(2,3:9)]

fwrite(psc1,paste(path_gwas,"resources/gwas_summary_statistics/GCST004030_only_SNPs_hm3_chrposid_sumstats_munged.sumstats",sep=""),
col.names=T,row.names=F,quote=F,sep="\t")



q("no")

################################

j=0
i=259

bsub -J"ldsc_pf" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${pheno[j]}_${gwas_id[i]}_ldsc_only_SNPs_hm3_hm3_weights_pf_mungen_format_sumstats_ldsc_script_stderr \
-o ${path_gwas}post_imputation/log/${pheno[j]}_${gwas_id[i]}_ldsc_only_SNPs_hm3_hm3_weights_pf_mungen_format_sumstats_ldsc_script_stdout \
"ldsc.py \
--rg ${path_gwas}resources/gwas_summary_statistics/${pheno[j]}_eur_tier1_list_variants_only_SNPs_hm3_rsids_rate_0.5_nohla_sumstats_munged.sumstats,${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_only_SNPs_hm3_rsids_sumstats_munged.sumstats \
--w-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/GRCh38/weights/weights.hm3_noMHC. \
--ref-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/GRCh38/weights/weights.hm3_noMHC. \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/${pheno[j]}_${gwas_id[i]}_hm3_hm3_weights_gc"


bsub -J"ldsc_pf" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${pheno[j]}_${gwas_id[i]}_ldsc_only_SNPs_hm3_pf_mungen_format_sumstats_ldsc_script_stderr \
-o ${path_gwas}post_imputation/log/${pheno[j]}_${gwas_id[i]}_ldsc_only_SNPs_hm3_pf_mungen_format_sumstats_ldsc_script_stdout \
"ldsc.py \
--rg ${path_gwas}resources/gwas_summary_statistics/${pheno[j]}_eur_tier1_list_variants_only_SNPs_hm3_chrposid_rate_0.5_nohla_sumstats_munged.sumstats,${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_only_SNPs_hm3_chrposid_sumstats_munged.sumstats \
--no-check-alleles \
--w-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/regression_weight_files/weights_only_SNPs.june2024. \
--ref-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/regression_weight_files/weights_only_SNPs.june2024. \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/${pheno[j]}_${gwas_id[i]}_hm3_gc"

# ukbb ld 
# downloaded from : https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.ALL.ldscore.tar.gz




#######################################################################
# test using UKB weigths - rename variants to rsid

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=25000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q yesterday R

library(data.table)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"
pheno<-c("ibd","cd","uc")

list_hapmap_variants<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/w_hm3.snplist",sep=""))

# list of IDs on UKB ldsc
ukbldsc<-fread("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/UKBB.ALL.ldscore/UKBB.EUR.rsid.l2.ldscore.gz")
ukbvar<-fread("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/UKBB.ALL.ldscore/full_variant_qc_metrics.txt.bgz")
ukbvar<-ukbvar[which(ukbvar$rsid %in% ukbldsc$SNP),]
ukbvar$ID_tmp<-paste(ukbvar$rsid,ukbvar$ref,ukbvar$alt,sep=":")

dim(ukbldsc)
dim(ukbvar)

# convert psc and ibd into ukb variants:

psc1<-fread(paste(path_gwas,"resources/gwas_summary_statistics/GCST004030_only_SNPs_hm3_rsids_sumstats_munged.sumstats",sep=""))
psc1$ID_tmp<-paste(psc1$SNP,psc1$A2,psc1$A1,sep=":")
psc1$ID_tmp2<-paste(psc1$SNP,psc1$A1,psc1$A2,sep=":")

dim(psc1)
dim(psc1[which(psc1$ID_tmp %in% ukbvar$ID_tmp),])
# [1] 944844      9
dim(psc1[which(psc1$ID_tmp2 %in% ukbvar$ID_tmp),])
# [1] 1864   10

dim(psc_tmp1)
# [1] 944844     10
psc_tmp1<-psc1[which(psc1$ID_tmp %in% ukbvar$ID_tmp),]
psc_tmp1<-merge(psc_tmp1,ukbvar[,c("pos","ref","alt","ID_tmp")],by="ID_tmp")
dim(psc_tmp1)
# [1] 944844     13

# save again rsid+b37 as well as rsid+b38 to run both analysis togeher and compare ldsc ukb vs ldsc 1000gp with same N variants

psc_b38<-psc_tmp1[,c("SNP","CHR","BP","A2","A1","N","MAF","Z")]
psc_b37<-psc_tmp1[,c("SNP","CHR","pos","A2","A1","N","MAF","Z")]
colnames(psc_b37)<-colnames(psc_b38)

fwrite(psc_b37,paste(path_gwas,"resources/gwas_summary_statistics/GCST004030_only_SNPs_hm3_rsid_b37_sumstats_munged.sumstats",sep=""),
col.names=T,row.names=F,quote=F,sep="\t")

fwrite(psc_b38,paste(path_gwas,"resources/gwas_summary_statistics/GCST004030_only_SNPs_hm3_rsid_b38_sumstats_munged.sumstats",sep=""),
col.names=T,row.names=F,quote=F,sep="\t")

# save chrpos to compare with our own weights:
psc_b38$SNP<-paste("chr",psc_b38$CHR,":",psc_b38$BP,":",psc_b38$A2,":",psc_b38$A1,sep="")
head(psc_b38)
fwrite(psc_b38,paste(path_gwas,"resources/gwas_summary_statistics/GCST004030_only_SNPs_hm3_chrposid_b38_sumstats_munged.sumstats",sep=""),
col.names=T,row.names=F,quote=F,sep="\t")

# do the same for the IBD phenotypes

for (ph in pheno) {

    print(ph)

    df<-fread(paste(path_gwas,"resources/gwas_summary_statistics/",ph,"_eur_tier1_list_variants_only_SNPs_hm3_rsids_rate_0.5_nohla_sumstats_munged.sumstats",sep=""),)
    print(head(df))

    df$ID_tmp<-paste(df$SNP,df$A2,df$A1,sep=":")
    dim(df) 
    dim(df[which(df$ID_tmp %in% ukbvar$ID_tmp),])
    # [1] 944844      9

    df_tmp1<-df[which(df$ID_tmp %in% ukbvar$ID_tmp),]
    df_tmp1<-merge(df_tmp1,ukbvar[,c("pos","ref","alt","ID_tmp")],by="ID_tmp")
    dim(df_tmp1)
    # [1] 944844     13

    # save again rsid+b37 as well as rsid+b38 to run both analysis togeher and compare ldsc ukb vs ldsc 1000gp with same N variants

    df_b38<-df_tmp1[,c("SNP","CHR","BP","A2","A1","N","MAF","Z")]
    df_b37<-df_tmp1[,c("SNP","CHR","pos","A2","A1","N","MAF","Z")]
    colnames(df_b37)<-colnames(df_b38)

    fwrite(df_b37,paste(path_gwas,"resources/gwas_summary_statistics/",ph,"_eur_tier1_list_variants_only_SNPs_hm3_rsid_b37_rsids_rate_0.5_nohla_sumstats_munged.sumstats",sep=""),
    col.names=T,row.names=F,quote=F,sep="\t")

    fwrite(df_b38,paste(path_gwas,"resources/gwas_summary_statistics/",ph,"_eur_tier1_list_variants_only_SNPs_hm3_rsid_b38_rsids_rate_0.5_nohla_sumstats_munged.sumstats",sep=""),
    col.names=T,row.names=F,quote=F,sep="\t")

    # save chrpos to compare with our own weights:
    df_b38$SNP<-paste("chr",df_b38$CHR,":",df_b38$BP,":",df_b38$A2,":",df_b38$A1,sep="")
    head(df_b38)
    fwrite(df_b38,paste(path_gwas,"resources/gwas_summary_statistics/",ph,"_eur_tier1_list_variants_only_SNPs_hm3_chrposid_b38_chrposid_rate_0.5_nohla_sumstats_munged.sumstats",sep=""),
    col.names=T,row.names=F,quote=F,sep="\t")

    rm(df,df_tmp1,df_b37,df_b38)

}

q("no")

###################

path_gwas="/path/to/ibdgwas/IIBDGC/"

i=259
pheno=(ibd cd uc)

for j in {0..2}
do
bsub -J"ldsc_pf" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${pheno[j]}_${gwas_id[i]}_ldsc_only_SNPs_hm3_rsid_b38_hm3_weights_pf_mungen_format_sumstats_ldsc_script_stderr \
-o ${path_gwas}post_imputation/log/${pheno[j]}_${gwas_id[i]}_ldsc_only_SNPs_hm3_rsid_b38_hm3_weights_pf_mungen_format_sumstats_ldsc_script_stdout \
"ldsc.py \
--rg ${path_gwas}resources/gwas_summary_statistics/${pheno[j]}_eur_tier1_list_variants_only_SNPs_hm3_rsid_b38_rsids_rate_0.5_nohla_sumstats_munged.sumstats,${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_only_SNPs_hm3_rsid_b38_sumstats_munged.sumstats \
--w-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/GRCh38/weights/weights.hm3_noMHC. \
--ref-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/GRCh38/weights/weights.hm3_noMHC. \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/${pheno[j]}_${gwas_id[i]}_hm3_rsid_b38_hm3_weights_gc"
done

for j in {0..2}
do
bsub -J"ldsc_pf" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${pheno[j]}_${gwas_id[i]}_ldsc_only_SNPs_hm3_chrposid_b38_gsa_weights_pf_mungen_format_sumstats_ldsc_script_stderr \
-o ${path_gwas}post_imputation/log/${pheno[j]}_${gwas_id[i]}_ldsc_only_SNPs_hm3_chrposid_b38_gsa_weights_pf_mungen_format_sumstats_ldsc_script_stdout \
"ldsc.py \
--rg ${path_gwas}resources/gwas_summary_statistics/${pheno[j]}_eur_tier1_list_variants_only_SNPs_hm3_chrposid_b38_chrposid_rate_0.5_nohla_sumstats_munged.sumstats,${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_only_SNPs_hm3_chrposid_b38_sumstats_munged.sumstats \
--w-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/regression_weight_files/weights_only_SNPs.june2024. \
--ref-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/regression_weight_files/weights_only_SNPs.june2024. \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/${pheno[j]}_${gwas_id[i]}_hm3_chrposid_b38_gsa_weights_gc"
done

for j in {0..2}
do
bsub -J"ldsc_pf" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${pheno[j]}_${gwas_id[i]}_ldsc_only_SNPs_hm3_rsid_b37_hm3_weights_pf_mungen_format_sumstats_ldsc_script_stderr \
-o ${path_gwas}post_imputation/log/${pheno[j]}_${gwas_id[i]}_ldsc_only_SNPs_hm3_rsid_b37_hm3_weights_pf_mungen_format_sumstats_ldsc_script_stdout \
"ldsc.py \
--rg ${path_gwas}resources/gwas_summary_statistics/${pheno[j]}_eur_tier1_list_variants_only_SNPs_hm3_rsid_b37_rsids_rate_0.5_nohla_sumstats_munged.sumstats,${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_only_SNPs_hm3_rsid_b37_sumstats_munged.sumstats \
--w-ld ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--ref-ld ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/UKBB.ALL.ldscore/UKBB.EUR.rsid \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/${pheno[j]}_${gwas_id[i]}_hm3_rsid_b37_ukb_weights_gc"
done



##################################################################################################################################################################################
## CREATE .BGEN FILES FROM GSA USING SAME FILES (and N samples) USED FOR REGENIE STEP2, BUT CONVERTED AND SUBSET BY SNPs INCLUDED IN FINAL ANALYSIS:

# singularity exec iibdgc_postprocess_10_singularity.sif
path_gwas="/path/to/ibdgwas/IIBDGC/"


MEM=22000

for chr in {1..22} X
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 4 \
-e ${path_gwas}post_imputation/log/${chr}_ldsc_h3_ldsc_plink_stderr \
-o ${path_gwas}post_imputation/log/${chr}_ldsc_h3_ldsc_plink_stdout \
"plink2 \
--bgen ${path_gwas}post_imputation/2022/gsa/imputed_data/gsa_chr${chr}.dose.subset.bgen 'ref-first' \
--sample ${path_gwas}post_imputation/2022/gsa/imputed_data/gsa_chr${chr}.dose.subset.sample \
--keep ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/gsa_list_samples_included_in_ibd_analysis \
--extract <(cut -f1 ${path_gwas}resources/gwas_summary_statistics/ibd_eur_tier1_list_variants_only_SNPs_hm3_chrposid_b38_chrposid_rate_0.5_nohla_sumstats_munged.sumstats) \
--threads 4 \
--make-bed --out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_june2024/gsa_chr${chr}_hm3_SNPs_cd_uc_ibd_metaanalysis"
done

for chr in {1..22} X
do
echo ${chr} && tail -30 ${path_gwas}post_imputation/log/${chr}_ldsc_h3_ldsc_plink_stdout  | grep -E "Successfully|Exit"
done

wc -l ${path_gwas}resources/gwas_summary_statistics/ibd_eur_tier1_list_variants_only_SNPs_hm3_chrposid_b38_chrposid_rate_0.5_nohla_sumstats_munged.sumstats
# 944845 - with header

wc -l ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_june2024/gsa_chr*_hm3_SNPs_cd_uc_ibd_metaanalysis.bim
# 944844


### estimate regression weights, no partitioned:
# singularity exec polyfun_2_singularity.sif

MEM=35000
for chr in {1..22} X
do
bsub -J"ldsc_w" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q long \
-e ${path_gwas}post_imputation/log/${chr}_ldsc_h3_w_stderr \
-o ${path_gwas}post_imputation/log/${chr}_ldsc_h3_w_stdout \
"ldsc.py \
--bfile ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/plink_files_version_june2024/gsa_chr${chr}_hm3_SNPs_cd_uc_ibd_metaanalysis \
--l2 \
--ld-wind-kb 1000 \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/regression_weight_files/weights_only_hm3_SNPs.june2024.${chr}"
done


####### run analysis with those weights


### estimate regression weights, no partitioned:
module unload HGI/softpack/groups/team152/polyfun-2/1
# singularity exec iibdgc_postprocess_10_singularity.sif

path_gwas="/path/to/ibdgwas/IIBDGC/"

MEM=7000

i=259
pheno=(ibd cd uc)

for j in {0..2}
do
bsub -J"ldsc_pf" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${pheno[j]}_${gwas_id[i]}_ldsc_only_SNPs_hm3_chrposid_b38_gsa_weights_estimated_with_only_hm3_SNPs_pf_mungen_format_sumstats_ldsc_script_stderr \
-o ${path_gwas}post_imputation/log/${pheno[j]}_${gwas_id[i]}_ldsc_only_SNPs_hm3_chrposid_b38_gsa_weights_estimated_with_only_hm3_SNPs_pf_mungen_format_sumstats_ldsc_script_stdout \
"ldsc.py \
--rg ${path_gwas}resources/gwas_summary_statistics/${pheno[j]}_eur_tier1_list_variants_only_SNPs_hm3_chrposid_b38_chrposid_rate_0.5_nohla_sumstats_munged.sumstats,${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_only_SNPs_hm3_chrposid_b38_sumstats_munged.sumstats \
--w-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/regression_weight_files/weights_only_hm3_SNPs.june2024. \
--ref-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/regression_weight_files/weights_only_hm3_SNPs.june2024. \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/${pheno[j]}_${gwas_id[i]}_hm3_chrposid_b38_gsa_weights_estimated_with_only_hm3_SNPs_gc"
done





# compare ldscore from GSA, 1000GP and UKB

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=15000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G ibdgwas R

library(data.table)
library(ggplot2)
library(ggpubr)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"
pheno<-c("ibd","cd","uc")

# use the same set of variants as tested above:

vars<-fread(paste(path_gwas,"resources/gwas_summary_statistics/GCST004030_only_SNPs_hm3_rsid_b38_sumstats_munged.sumstats",sep=""))
head(vars)
#           SNP   CHR        BP     A2     A1     N   MAF          Z
#        <char> <int>     <int> <char> <char> <int> <num>      <num>
# 1:  rs1000000    12 126406434      G      A 14890  0.22 -0.4682528
# 2: rs10000010     4  21617051      T      C 14890  0.49 -0.8167857
# 3: rs10000023     4  94812755      G      T 14890  0.41  0.0666000
# 4:  rs1000002     3 183917980      C      T 14890  0.49 -0.8850931
# 5: rs10000033     4 138678744      T      C 14890  0.46  0.6273933
# 6: rs10000037     4  38922709      G      A 14890  0.25  0.7726351

vars$chrposid<-paste("chr",vars$CHR,":",vars$BP,":",vars$A2,":",vars$A1,sep="")
vars<-as.data.frame(vars)

# UKB ldsc
ukbldsc<-fread("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/UKBB.ALL.ldscore/UKBB.EUR.rsid.l2.ldscore.gz")
ukbldsc<-ukbldsc[which(ukbldsc$SNP %in% vars$SNP),]
dim(ukbldsc)
ukbldsc<-as.data.frame(ukbldsc)

# GSA ldsc
for (chr in c(1:22)) {
    tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/regression_weight_files/weights_only_SNPs.june2024.",chr,".l2.ldscore.gz",sep=""))
    tmp<-tmp[which(tmp$SNP %in% vars$chrposid),]

    if(chr==1) {
        gsaldsc<-tmp
    } else {
        gsaldsc<-rbind(gsaldsc,tmp)
    }

}
dim(gsaldsc)
# [1] 944844      4
gsaldsc<-as.data.frame(gsaldsc)
colnames(gsaldsc)[2]<-"chrposid"
gsaldsc<-merge(gsaldsc,vars[,c("SNP","chrposid")],by="chrposid",all.x=T)


# GSA ldsc - with only hm3 variants
for (chr in c(1:22)) {
    tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/regression_weight_files/weights_only_hm3_SNPs.june2024.",chr,".l2.ldscore.gz",sep=""))
    tmp<-tmp[which(tmp$SNP %in% vars$chrposid),]

    if(chr==1) {
        gsaldschm3<-tmp
    } else {
        gsaldschm3<-rbind(gsaldschm3,tmp)
    }

}
dim(gsaldschm3)
# [1] 944844      4
gsaldschm3<-as.data.frame(gsaldschm3)
colnames(gsaldschm3)[2]<-"chrposid"
gsaldschm3<-merge(gsaldschm3,vars[,c("SNP","chrposid")],by="chrposid",all.x=T)


# 1000GP ldsc
for (chr in c(1:22)) {
    tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/test/GRCh38/weights/weights.hm3_noMHC.",chr,".l2.ldscore.gz",sep=""))
    tmp<-tmp[which(tmp$SNP %in% vars$SNP),]

    if(chr==1) {
        oneldsc<-tmp
    } else {
        oneldsc<-rbind(oneldsc,tmp)
    }

}
dim(oneldsc)
oneldsc<-as.data.frame(oneldsc)

colnames(ukbldsc)[3:4]<-paste(colnames(ukbldsc)[3:4],"_ukb",sep="")
colnames(oneldsc)[3:4]<-paste(colnames(oneldsc)[3:4],"_1000gp_hm3",sep="")
colnames(gsaldsc)[3:4]<-paste(colnames(gsaldsc)[3:4],"_gsa",sep="")
colnames(gsaldschm3)[3:4]<-paste(colnames(gsaldschm3)[3:4],"_gsa_hm3",sep="")

all<-merge(ukbldsc[,c("SNP","L2_ukb")],oneldsc[,c("SNP","L2_1000gp_hm3")],by="SNP")
all<-merge(all,gsaldschm3[,c("SNP","L2_gsa_hm3")],by="SNP")
all<-merge(all,gsaldsc[,c("SNP","L2_gsa","chrposid")],by="SNP")


cor.test(all$L2_ukb, all$L2_1000gp, method=c("pearson"))
cor.test(all$L2_ukb, all$L2_gsa, method=c("pearson"))
cor.test(all$L2_1000gp, all$L2_gsa, method=c("pearson"))
cor.test(all$L2_gsa_hm3, all$L2_gsa, method=c("pearson"))
cor.test(all$L2_gsa_hm3, all$L2_1000gp, method=c("pearson"))


max(all$L2_1000gp,all$L2_gsa,all$L2_ukb)
# 1] 3781.3


### compare ldsc:

p1<-ggplot(all, aes(y=L2_gsa,x=L2_ukb)) + xlim(0,4000) + ylim(0,4000) + 
  geom_point() + stat_cor(method="pearson",label.x = 2500,label.y = 0) +
  geom_smooth(method=lm)

p2<-ggplot(all, aes(y=L2_1000gp_hm3, x=L2_ukb)) + xlim(0,4000) +
  geom_point() + stat_cor(method="pearson",label.x = 2500,label.y = 0) +
  geom_smooth(method=lm)

p3<-ggplot(all, aes(y=L2_1000gp_hm3, x=L2_gsa)) + xlim(0,4000) +
  geom_point() + stat_cor(method="pearson",label.x = 2500,label.y = 0) +
  geom_smooth(method=lm)

p4<-ggplot(all, aes(y=L2_1000gp_hm3, x=L2_gsa_hm3)) + 
  geom_point() + stat_cor(method="pearson",label.x = 150,label.y = 0) +
  geom_smooth(method=lm)

px<-ggplot() + theme_void()


p<-ggarrange(p1,px,px,p2,p3,p4,ncol=3,nrow=2)

ggsave(
  "~/git/IIBDGC_GWAS/plots/paper_figures/Supplementary_Figure_ldsc_reference_comparison.pdf",
  p,
  width = 150,
  height = 90,
  dpi = 400,
  units = c("mm"),
  limitsize = T,scale=4
)


### check out regions where ldsc is not consistent:

all[which(all$L2_ukb>2000 & all$L2_1000gp_hm3,100),]

all[which(all$L2_ukb>2000 & all$L2_1000gp_hm3<100),]
all[which(all$L2_ukb>2000 & all$L2_1000gp_hm3>100),]



all[which(all$L2_gsa>2000 & all$L2_1000gp_hm3<100),]
all[which(all$L2_gsa>2000 & all$L2_1000gp_hm3>100),]

################################


###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################

# ESTIMATE RG FOR ALL METABOLIC BIOMARKERS TO SUPPORT THE FACT THAT IE CHOLESTEROL ARE NOW INDEPENDENT:

files=($(cat /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/list_path_summary_stats_cholesterol)) 
echo ${#files[@]}
# 77


gwas_id=($(cat /path/to/ibdgwas/IIBDGC/resources/gwas_summary_statistics/list_path_summary_stats_cholesterol | \
sed 's/\/lustre\/scratch123\/hgi\/projects\/ibdgwas\/IIBDGC\/resources\/gwas_summary_statistics\///g' | \
sed 's/\/lustre\/scratch125\/humgen\/resources\/GWAScatalog\/.*\///g' | sed  's/.h.tsv.gz//g')) 
echo ${#gwas_id[@]}
# 77

module unload HGI/softpack/groups/team152/polyfun-2/1
# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=7300

for j in {0..76}
do 
for i in {0..76}
do 
bsub -J"ldsc_pf" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${gwas_id[j]}_${gwas_id[i]}_ldsc_only_SNPs_pf_mungen_format_sumstats_ldsc_script_stderr \
-o ${path_gwas}post_imputation/log/${gwas_id[j]}_${gwas_id[i]}_ldsc_only_SNPs_pf_mungen_format_sumstats_ldsc_script_stdout \
"ldsc.py \
--rg ${path_gwas}resources/gwas_summary_statistics/${gwas_id[j]}_only_SNPs_sumstats_munged.sumstats,${path_gwas}resources/gwas_summary_statistics/${gwas_id[i]}_only_SNPs_sumstats_munged.sumstats \
--no-check-alleles \
--w-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/regression_weight_files/weights_only_SNPs.june2024. \
--ref-ld-chr ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/regression_weight_files/weights_only_SNPs.june2024. \
--out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/genetic_correlation/${gwas_id[j]}_${gwas_id[i]}_metabolic_biomarkers_test_gc"
done
done

##########################################################################################################################################################

