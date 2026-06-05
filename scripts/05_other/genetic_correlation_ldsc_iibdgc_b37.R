# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
#genetic_correlation_ldsc_iibdgc.R

# note that this is a tmp script to generate some enrichmnet analysis for ECCO'24, the final version will have new weights and LDScores 
# genereated for all variants included in the analysis, so all can contribute to estimate the priors for FM

# ldsc files were donwloaded:
/path/to/ibdgwas/IIBDGC/post_imputation/analysis/metaanalysis/annotation/ldscore/test/

# hapmap3 SNPs
wget https://zenodo.org/records/8292725/files/hm3_no_MHC.list.txt?download=1

# weights
wget https://zenodo.org/records/8292725/files/1000G_Phase3_weights_hm3_no_MHC.tgz?download=1
tar -xvzf 1000G_Phase3_weights_hm3_no_MHC.tgz

# frequencies:
wget https://zenodo.org/records/8292725/files/1000G_Phase3_frq.tgz?download=1
tar -xvzf 1000G_Phase3_frq.tgz

# ldscores:
wget https://zenodo.org/records/8292725/files/1000G_Phase3_ldscores.tgz?download=1
tar -xvzf 1000G_Phase3_ldscores.tgz

# there rationale to use v1.2 here https://github.com/RHReynolds/LDSCforRyten instead of new releases 
wget wget https://zenodo.org/records/8292725/files/1000G_Phase3_baseline_v1.2_ldscores.tgz?download=1
tar -xvzf 1000G_Phase3_baselineLD_v1.2_ldscores.tgz

# plink files:
wget https://zenodo.org/records/7768714/files/1000G_Phase3_plinkfiles.tgz?download=1
tar -xvzf 1000G_Phase3_plinkfiles.tgz

wget https://ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/eur_w_ld_chr/w_hm3.snplist



########################################################################################################
# 1 - liftover summary statistics to b37

# combine all variants from the three analyses into one unique bed (UCSC format):

# MEM=4000
# bsub -Is -m "modern_hardware" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -n 4 R \

library(data.table)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"
pheno<-c("cd","ibd","uc")

for (chr in c(1:22)) {

    print(chr)

    for (i in 1:length(pheno)) {
        tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr,"_",pheno[i],"_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""),head=T)
        tmp<-tmp[,c("Position_b38","MarkerName")]
        tmp$chr<-paste("chr",chr,sep="")
        tmp$end<-as.numeric(as.character(tmp$Position_b38))+1
        tmp<-tmp[,c("chr","Position_b38","end","MarkerName")]

        if(i==1) {
            dat<-tmp
        } else {
            dat<-rbind(dat,tmp)
            dat<-dat[!duplicated(dat)]
        }

    }

    dat<-dat[order(dat$Position_b38,decreasing=F),]
    dat$Position_b38<-format(dat$Position_b38, scientific=F)
    dat$end<-format(dat$end, scientific=F)
    fwrite(dat,paste(path_gwas,"post_imputation/analysis/metaanalysis/annotation/ldscore/test/iibdgc_b37/",chr,"_meta_eur_tier_2_b38.bed",sep=""),
    col.names=F,row.names=F,quote=F,sep="\t")
    rm(dat)
}

q("no")


path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=200

for chr in {1..22}
do 
bsub -J"lift" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}post_imputation/log/stderr_liftover_to_b37_${chr} \
-o ${path_gwas}post_imputation/log/stdout_liftover_to_b37_${chr} \
"${path_gwas}previous_qced_b38/liftover/liftOver \
${path_gwas}post_imputation/analysis/metaanalysis/annotation/ldscore/test/iibdgc_b37/${chr}_meta_eur_tier_2_b38.bed \
${path_gwas}previous_qced_b38/liftover/hg38ToHg19.over.chain.gz \
${path_gwas}post_imputation/analysis/metaanalysis/annotation/ldscore/test/iibdgc_b37/${chr}_meta_eur_tier_2_lifted_hg37 \
${path_gwas}post_imputation/analysis/metaanalysis/annotation/ldscore/test/iibdgc_b37/${chr}_meta_eur_tier_2_no_lifted_hg37"
done


########################################################################################################
# 2 - get the alternative allele frequency for the variants included in ldsc:

MEM=600
for chr in {1..22}
do
bsub -J"ldsc" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -m 'modern_hardware'  -n 4 \
-e /path/to/ibdgwas/IIBDGC/post_imputation/analysis/metaanalysis/log/${chr}_plink_freq_stderr \
-o /path/to/ibdgwas/IIBDGC/post_imputation/analysis/metaanalysis/log/${chr}_plink_frea_stdout \
"/path/to/software/./plink2 \
--bfile ${path_gwas}post_imputation/analysis/metaanalysis/annotation/ldscore/test/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
--freq --out  ${path_gwas}post_imputation/analysis/metaanalysis/annotation/ldscore/test/1000G_Phase3_frq/1000G.EUR.QC.${chr}"
done


########################################################################################################
# 3 - combine IIBDGC eur_tier2 summary stats with liftover position; intersect with the variants in ldsc (compare allele freq)

# A note on file formats
# PolyFun uses input files that are very similar to the input files of S-LDSC. The main differences are:
# The .annot files must contain two additional columns called A1,A2 which encode the identifies of the effect and alternative allele 
# (i.e., the sign of effect sizes is with respect to A1)
# The .l2.ldscore files may contain the additional columns A1,A2. We strongly encourage including these columns.
# Polyfun supports files in .parquet format in addition to .gzip/.bzip2 formats. Parquet files can be loaded substantially faster than 
# alternative formats, at the cost of slightly larger file sizes.


# MEM=15000
# bsub -Is -m "modern_hardware" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -n 2 R \

library(data.table)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

args = commandArgs(trailingOnly=TRUE)

# phenotypes<-c("ibd","cd","uc")
phenotypes<-c("cd","uc")

# keep only variants in the plink files supplied by ldscore - as tmp solution before building our weight files
for (chr in c(1:22)) {

    tmp1<-fread(paste("/path/to/ibdgwas/IIBDGC/post_imputation/analysis/metaanalysis/annotation/ldscore/test/1000G_EUR_Phase3_plink/1000G.EUR.QC.",chr,".bim",sep=""),head=F)
    tmp2<-fread(paste("/path/to/project",chr,".afreq",sep=""),head=T)

    tmp<-merge(tmp1,tmp2,by.x="V2",by.y="ID",)
    rm(tmp1,tmp2)

    if(chr==1) {
        bim<-tmp
    } else {
        bim<-rbind(bim,tmp)
    }
    rm(tmp)
}

bim$MarkerName<-paste("chr",bim$'#CHROM',":",bim$V4,":",bim$REF,":",bim$ALT,sep="")
bim$MarkerName_2<-paste("chr",bim$'#CHROM',":",bim$V4,":",bim$ALT,":",bim$REF,sep="")
# table(bim$V1)
#      1      2      3      4      5      6      7      8      9     10     11 
# 779354 839590 706350 729645 633015 664016 589569 549971 438106 510501 493922 
#     12     13     14     15     16     17     18     19     20     21     22 
    # 480110 366200 324698 287001 316981 269222 285156 232363 221626 138712 141123

ids<-c(bim$MarkerName,bim$MarkerName_2)

phenotypes<-c("cd","uc")

for (pheno in phenotypes) {

    print(pheno)

    # # keep variants in hapmap3 list:
    # ref<-fread("/path/to/ibdgwas/IIBDGC/post_imputation/analysis/metaanalysis/annotation/ldscore/test/w_hm3.snplist")

    # dim(ref)
    # [1] 1217311       3 



    # load eur_tier_2 dataset ~ LD as reference used here:

    for (chr in c(1:22)) {

        tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno,"/",chr,"_",pheno,"_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""),head=T)

        # map to b37
        tmp2<-fread(paste(path_gwas,"post_imputation/analysis/metaanalysis/annotation/ldscore/test/iibdgc_b37/",chr,"_meta_eur_tier_2_lifted_hg37",sep=""),head=F)
        colnames(tmp2)[4]<-"MarkerName"
        colnames(tmp2)[2]<-"position_b37"

        tmp<-merge(tmp,tmp2[,c("MarkerName","position_b37")])
        rm(tmp2)

        tmp$MarkerName_b37<-paste("chr",chr,":",tmp$position_b37,":",tmp$A1,":",tmp$A2,sep="")
        tmp<-tmp[which(tmp$MarkerName_b37 %in% ids),]
    
        if(chr==1) {
            dat<-tmp
        } else {
            dat<-rbind(dat,tmp)
        }
        rm(tmp)
    }

    dat$chr<-gsub(":.*","",dat$MarkerName_b37)
    dat$chr<-gsub("chr","",dat$chr)

    dat<-dat[,c("chr","MarkerName_b37","position_b37","A2","A1","N","N_CASES","N_CONTROLS","avgA2FREQ","BETA","SE","P-value","INFO")]
    dat<-as.data.frame(dat)

    # samples with same ref and alt alleles:
    dat1<-merge(dat,bim[,c("MarkerName","REF","ALT","ALT_FREQS","V2")],by.x="MarkerName_b37",by.y="MarkerName")
    dim(dat1)

    dat2<-merge(dat,bim[,c("MarkerName_2","REF","ALT","ALT_FREQS","V2")],by.x="MarkerName_b37",by.y="MarkerName_2")
    dim(dat2)

    # flip beta and frequency:
    dat2$BETA_ed<-as.numeric(dat2$BETA)*-1
    dat2$avgA2FREQ_ed<-1-as.numeric(dat2$avgA2FREQ)
    dat2<-dat2[,c("MarkerName_b37","chr","position_b37","A2","A1","N","N_CASES","N_CONTROLS","avgA2FREQ_ed","BETA_ed","SE","P-value","INFO","REF","ALT","ALT_FREQS","V2")]
    colnames(dat2)<-colnames(dat1)

    dat<-rbind(dat1,dat2)
    rm(dat1,dat2)

    # keep only those whose freq matches the expected freq:

    dat$value<-((dat$avgA2FREQ-dat$ALT_FREQS)^2)/
        ((dat$avgA2FREQ+dat$ALT_FREQS)*
            (2-dat$avgA2FREQ-dat$ALT_FREQS))
    summary(dat$value)

    dat<-dat[which(dat$value<0.125),]

    # exclude dup SNPs introduced during liftover:
    dups<-dat$V2[which(duplicated(dat$V2))]
    dat<-dat[which(!dat$V2 %in% dups),]

    dat$X<-paste(dat$chr,dat$position_b37,sep="_")
    dups<-dat$X[which(duplicated(dat$X))]
    dat<-dat[which(!dat$X %in% dups),]

    # re0run keeeping the SNP id as in ref files
    dat<-dat[,c("V2","chr","position_b37","REF","ALT","avgA2FREQ","BETA","SE","P-value","N_CASES","N_CONTROLS","INFO")]
    colnames(dat)<-c("SNP","CHR","BP","A2","A1","A1FREQ","BETA","SE","pval","N_CASES","N_CONTROLS","INFO")

    dat$CHR<-as.numeric(dat$CHR)
    dat$BP<-as.numeric(dat$BP)

    dat<-dat[order(dat$CHR,dat$BP,decreasing=F),]

    # dat<-fread(paste(path_gwas,"post_imputation/analysis/metaanalysis/annotation/ldscore/test/iibdgc_b37/",pheno,"_summary_stats_eur_tier_2.txt.gz",sep=""),head=T)

    dat$N<-dat$N_CASES+dat$N_CONTROLS
    nmax<-max(dat$N)

    # retain only variants with >40% of samples contributing
    dat<-dat[which(dat$N>=(nmax*0.4)),]

    dat<-dat[,c("SNP","CHR","BP","A2","A1","A1FREQ","BETA","SE","pval","N_CASES","N_CONTROLS","INFO")]

    dat$BP<-format(dat$BP, scientific=F)

    fwrite(dat,paste(path_gwas,"post_imputation/analysis/metaanalysis/annotation/ldscore/test/iibdgc_b37/",pheno,"_summary_stats_eur_tier_2.txt.gz",sep="")
    ,col.names=T,row.names=F,quote=F,sep="\t")

    # double check that the variants are also included in the weights files:
    ref<-fread(paste(path_gwas,"post_imputation/analysis/metaanalysis/annotation/ldscore/test/w_hm3.snplist",sep=""),head=T)
    dim(ref)
    # [1] 1217311       3
    dim(dat[which(dat$SNP %in% ref$SNP),])
    # [1] 1186122      13


    # create a bed file with all the variants retained

    dat<-dat[,c("CHR","BP","SNP")]
    dat$CHR<-paste("chr",dat$CHR,sep="")
    dat$end<-as.numeric(dat$BP)+1
    colnames(dat)<-c("CHR","start","SNP","end")
    dat<-dat[,c("CHR","start","end","SNP")]

    dat$start<-format(dat$start, scientific=F)
    dat$end<-format(dat$end, scientific=F)

    dat<-dat[!duplicated(dat),]

    fwrite(dat,paste(path_gwas,"post_imputation/analysis/metaanalysis/annotation/ldscore/test/iibdgc_b37/",pheno,"_summary_stats_eur_tier_2.bed",sep=""),
    col.names=F,row.names=F,sep="\t")
    rm(dat)

}

q("no")

#########################

conda activate polyfun

# to run cd and uc with rsids
path_gwas=/path/to/ibdgwas/IIBDGC/

MEM=5500

pheno=(ibd cd uc)

for ph in ${pheno[@]}
do
bsub -J"munge" -m "modern_hardware" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/analysis/metaanalysis/log/munge_polyfun_stderr \
-o ${path_gwas}post_imputation/analysis/metaanalysis/log/munge_polyfun_stdout \
"python /software/hgi/installs/polyfun/./munge_polyfun_sumstats.py \
--sumstats ${path_gwas}post_imputation/analysis/metaanalysis/annotation/ldscore/test/iibdgc_b37/${ph}_summary_stats_eur_tier_2.txt.gz \
--min-info 0.4 \
--min-maf 0.001 \
--out ${path_gwas}post_imputation/analysis/metaanalysis/annotation/ldscore/test/iibdgc_b37/${ph}_summary_stats_eur_tier_2_munged.parquet"
done

# [INFO]  Reading sumstats file...
# [INFO]  Done in 20.08 seconds
# [INFO]  9421936 SNPs are in the sumstats file
# [INFO]  Removing 47226 HLA SNPs
# [INFO]  9374710 SNPs with sumstats remained after all filtering stages
# [INFO]  Computing the effective sample size for case-control data...
# [INFO]  Saving munged sumstats of 9374710 SNPs to /path/to/ibdgwas/IIBDGC/post_imputation/analysis/metaanalysis/annotation/ldscore/test/iibdgc_b
# 37/uc_summary_stats_eur_tier_2_munged.parquet
# [INFO]  Done

# [INFO]  Reading sumstats file...
# [INFO]  Done in 32.50 seconds
# [INFO]  9450456 SNPs are in the sumstats file
# [INFO]  Removing 47339 HLA SNPs
# [INFO]  9403117 SNPs with sumstats remained after all filtering stages
# [INFO]  Computing the effective sample size for case-control data...
# [INFO]  Saving munged sumstats of 9403117 SNPs to /path/to/ibdgwas/IIBDGC/post_imputation/analysis/metaanalysis/annotation/ldscore/test/iibdgc_b
# 37/ibd_summary_stats_eur_tier_2_munged.parquet
# [INFO]  Done

# [INFO]  Reading sumstats file...
# [INFO]  Done in 35.37 seconds
# [INFO]  9446728 SNPs are in the sumstats file
# [INFO]  Removing 47331 HLA SNPs
# [INFO]  9399397 SNPs with sumstats remained after all filtering stages
# [INFO]  Computing the effective sample size for case-control data...
# [INFO]  Saving munged sumstats of 9399397 SNPs to /path/to/ibdgwas/IIBDGC/post_imputation/analysis/metaanalysis/annotation/ldscore/test/iibdgc_b
# 37/cd_summary_stats_eur_tier_2_munged.parquet
# [INFO]  Done


# run partitioned heritability - test files

MEM=4500

for ph in ${pheno[@]}
do
bsub -J"munge" -m "modern_hardware" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/analysis/metaanalysis/log/${ph}_baseline_partitioned_polyfun_stderr \
-o ${path_gwas}post_imputation/analysis/metaanalysis/log/${ph}_baseline_partitioned_polyfun_stdout \
"python /software/hgi/installs/polyfun/./ldsc.py \
--h2 ${path_gwas}post_imputation/analysis/metaanalysis/annotation/ldscore/test/iibdgc_b37/${ph}_summary_stats_eur_tier_2_munged.parquet \
--w-ld-chr /path/to/ibdgwas/IIBDGC/post_imputation/analysis/metaanalysis/annotation/ldscore/test/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--ref-ld-chr /path/to/ibdgwas/IIBDGC/post_imputation/analysis/metaanalysis/annotation/ldscore/test/baseline_v1.2/baseline. \
--frqfile-chr /path/to/ibdgwas/IIBDGC/post_imputation/analysis/metaanalysis/annotation/ldscore/test/1000G_Phase3_frq/1000G.EUR.QC. \
--overlap-annot \
--thin-annot \
--print-coefficients \
--print-delete-vals \
--out ${path_gwas}post_imputation/analysis/metaanalysis/annotation/ldscore/test/iibdgc_b37/${ph}_summary_stats_eur_tier_2_munged_baseline"
done

# IBD
# subset by maf >0.001 info >0.4 and % samples >0.4
# Total Observed scale h2: 0.1525 (0.0089)
# Lambda GC: 1.0197e-06
# Mean Chi^2: 1.4654e-06
# Intercept: 1.0784 (0.0247)
# Ratio: NA (mean chi^2 < 1)

# subset by maf >0.001 info >0.4
# Total Observed scale h2: 0.1521 (0.0088)
# Lambda GC: 1.0156e-06
# Mean Chi^2: 1.4599e-06
# Intercept: 1.0798 (0.0242)
# Ratio: NA (mean chi^2 < 1)

# CD
# subset by maf >0.001 info >0.4 and % samples >0.4
# Total Observed scale h2: 0.2236 (0.0136)
# Lambda GC: 9.6134e-07
# Mean Chi^2: 1.3063e-06
# Intercept: 1.0467 (0.0218)
# Ratio: NA (mean chi^2 < 1)

# UC
# subset by maf >0.001 info >0.4 and % samples >0.4
# Total Observed scale h2: 0.2213 (0.0118)
# Lambda GC: 9.7124e-07
# Mean Chi^2: 1.3037e-06
# Intercept: 1.0684 (0.0186)
# Ratio: NA (mean chi^2 < 1)


########################################################################################################
# 4 - Get enrichment:


cd /path/to/project
files_atac=($(ls | grep allchr_ATAC_counts))

echo ${#files_atac[@]}
# 45
  
files_atac=($(echo ${files_atac[*]} | sed 's/\.bed\.gz//g' | sed 's/allchr_//g' | sed 's/\-/\_/g'))

for i in ${files_atac[@]}
do
echo ${i}
done


MEM=5500

conda activate polyfun

MEM=5500
pheno=(cd uc)
# baseline_v1.2
for ph in ${pheno[@]}
do
for i in ${files_atac[@]}
do
bsub -J"munge" -m "modern_hardware" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/analysis/metaanalysis/log/${ph}_baseline_partitioned_v1.2_${i}_polyfun_stderr \
-o ${path_gwas}post_imputation/analysis/metaanalysis/log/${ph}_baseline_partitioned_v1.2_${i}_polyfun_stdout \
"python /software/hgi/installs/polyfun/./ldsc.py \
--h2 ${path_gwas}post_imputation/analysis/metaanalysis/annotation/ldscore/test/iibdgc_b37/${ph}_summary_stats_eur_tier_2_munged.parquet \
--w-ld-chr /path/to/ibdgwas/IIBDGC/post_imputation/analysis/metaanalysis/annotation/ldscore/test/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--ref-ld-chr /path/to/ibdgwas/IIBDGC/post_imputation/analysis/metaanalysis/annotation/ldscore/test/baseline_v1.2/baseline.,/path/to/project \
--frqfile-chr /path/to/ibdgwas/IIBDGC/post_imputation/analysis/metaanalysis/annotation/ldscore/test/1000G_Phase3_frq/1000G.EUR.QC. \
--overlap-annot \
--print-coefficients \
--print-delete-vals \
--out ${path_gwas}post_imputation/analysis/metaanalysis/annotation/ldscore/test/iibdgc_b37/enrichment/${ph}_summary_stats_munged_baseline_v1.2_${i}"
done
done

###

cd /path/to/project
files_cre=($(ls | grep allchr_cCRE_Accessibility))

echo ${#files_cre[@]}
# 54
  
files_cre=($(echo ${files_cre[*]} | sed 's/\.bed\.gz//g' | sed 's/allchr_//g' | sed 's/\-/\_/g'))

for i in ${files_cre[@]}
do
echo ${i}
done


# baseline_v1.2
for ph in ${pheno[@]}
do
for i in ${files_cre[@]}
do
bsub -J"munge" -m "modern_hardware" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e /path/to/project \
-o /path/to/project \
"python /software/hgi/installs/polyfun/./ldsc.py \
--h2 ${path_gwas}post_imputation/analysis/metaanalysis/annotation/ldscore/test/iibdgc_b37/${ph}_summary_stats_eur_tier_2_munged.parquet \
--w-ld-chr /path/to/ibdgwas/IIBDGC/post_imputation/analysis/metaanalysis/annotation/ldscore/test/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--ref-ld-chr /path/to/ibdgwas/IIBDGC/post_imputation/analysis/metaanalysis/annotation/ldscore/test/baseline_v1.2/baseline.,/path/to/project \
--frqfile-chr /path/to/ibdgwas/IIBDGC/post_imputation/analysis/metaanalysis/annotation/ldscore/test/1000G_Phase3_frq/1000G.EUR.QC. \
--overlap-annot \
--print-coefficients \
--print-delete-vals \
--out ${path_gwas}post_imputation/analysis/metaanalysis/annotation/ldscore/test/iibdgc_b37/enrichment/${ph}_summary_stats_munged_baseline_v1.2_${i}"
done
done

##########################################################################################################################################

# MEM=8000
# bsub -Is -m "modern_hardware" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -n 2  /software/bin/R-4.3.1 \

library(data.table)
library(ggplot2)
library(ggpubr)
library(colorspace)
library(rcartocolor)
library(qvalue)

set.seed(10) 

files<-list.files("/path/to/ibdgwas/IIBDGC/post_imputation/analysis/metaanalysis/annotation/ldscore/test/iibdgc_b37/enrichment")
files<-files[grep(".results",files)]
files_ref<-files[grep("_v1.2",files)]

pheno<-c("cd","uc")

for (ph in c(1:length(pheno))) {

    print(pheno[ph])
    
    files<-files_ref[grep(pheno[ph],files_ref)]
    
    ##########
    # zhang

    files_cre<-files[grep("cCRE",files)]
    length(files_cre)
    # [1] 54

    #######################
    # ATAC-SEQ

    data<-gsub("_b37_psc_summary_stats.results","",files_cre)
    data<-gsub(paste(pheno[ph],"_summary_stats_munged_baseline_v1.2_cCRE_Accessibility_",sep=""),"cCRE_",data)

    for (i in 1:length(files_cre)) {

        tmp<-fread(paste("/path/to/ibdgwas/IIBDGC/post_imputation/analysis/metaanalysis/annotation/ldscore/test/iibdgc_b37/enrichment/",files_cre[i],sep=""),head=T)
        tmp$Enrichment_p[which(is.na(tmp$Enrichment_p))]<-1
        tmp$qvalue<-qvalue(tmp$Enrichment_p)$qvalue
        tmp<-tmp[which(tmp$Category=="L2_1",)]
        tmp$Category<-data[i]

        if(i==1) {
            enrich<-tmp
        } else {
            enrich<-rbind(enrich,tmp)
        }
    }

    enrich$cell_type<-gsub("cCRE_","",enrich$Category)
    enrich$cell_type<-gsub("_b37_psc_summary_stats.results$","",enrich$cell_type)

    enrich$cell<-enrich$cell_type
    enrich$cell<-gsub("\\.[0-9]{1}$","",enrich$cell)
    enrich$cell<-gsub("[0-9]","",enrich$cell)
    enrich$cell[which(enrich$cell=="Adi")]<-"Adipocyte"
    enrich$cell[which(enrich$cell=="Mes")]<-"Mesothelial cell"
    enrich$cell[which(enrich$cell=="Grn")]<-"Granulosa cell"
    enrich$cell[which(enrich$cell=="Smm")]<-"Smooth muscle"
    enrich$cell[which(enrich$cell=="Vsm")]<-"Vascular smooth muscle"
    enrich$cell[which(enrich$cell=="Esm")]<-"Esophageal smooth muscle"
    enrich$cell[which(enrich$cell=="Mfb")]<-"Myofibroblast"
    enrich$cell[which(enrich$cell=="Swn")]<-"Schwann cell"
    enrich$cell[which(enrich$cell=="Fib")]<-"Fibroblast"
    enrich$cell[which(enrich$cell=="Stl")]<-"Myosatellite cell"
    enrich$cell[which(enrich$cell=="Skm")]<-"Skeletal myocyte"
    enrich$cell[which(enrich$cell=="Cam")]<-"Cardiomyocyte"
    enrich$cell[which(enrich$cell=="Mac")]<-"Macrophage"
    enrich$cell[which(enrich$cell=="Tly")]<-"T-cell"
    enrich$cell[which(enrich$cell=="Bly")]<-"B-cell"
    enrich$cell[which(enrich$cell=="Mst")]<-"Mast cell"
    enrich$cell[which(enrich$cell=="End")]<-"Endothelial cell"
    enrich$cell[which(enrich$cell=="Msc")]<-"Miscellaneous stromal cell"
    enrich$cell[which(enrich$cell=="Pal")]<-"Pheumocyte"
    enrich$cell[which(enrich$cell=="Krt")]<-"Keratinocyte"
    enrich$cell[which(enrich$cell=="Eep")]<-"Esophageal epithelial cell"
    enrich$cell[which(enrich$cell=="Enc")]<-"Enterocyte"
    enrich$cell[which(enrich$cell=="Gbl")]<-"Goblet cell"
    enrich$cell[which(enrich$cell=="Fol")]<-"Follicular cell of thyroid"
    enrich$cell[which(enrich$cell=="Lue")]<-"Luminal epithelial cell"
    enrich$cell[which(enrich$cell=="Epc")]<-"Epithelial cells"
    enrich$cell[which(enrich$cell=="Bas")]<-"Basal cell"
    enrich$cell[which(enrich$cell=="Agb")]<-"Airway goblet cell"
    enrich$cell[which(enrich$cell=="Dut")]<-"Ductal cell"
    enrich$cell[which(enrich$cell=="Acn")]<-"Acinar cell of pancreas"
    enrich$cell[which(enrich$cell=="Prt")]<-"Gastric parietal cell"
    enrich$cell[which(enrich$cell=="Gcf")]<-"Gastric chief cell"
    enrich$cell[which(enrich$cell=="Adc")]<-"Cortical cell of adrenal gland"
    enrich$cell[which(enrich$cell=="Nec")]<-"Neuroendocrine cell"
    enrich$cell[which(enrich$cell=="Hpc")]<-"Hepatocyte"

    enrich$cell<-as.factor(enrich$cell)
    levels(enrich$cell)
    #  [1] "Acinar cell of pancreas"        "Adipocyte"                     
    #  [3] "Airway goblet cell"             "Basal cell"                    
    #  [5] "B-cell"                         "Cardiomyocyte"                 
    #  [7] "Cortical cell of adrenal gland" "Ductal cell"                   
    #  [9] "Endothelial cell"               "Enterocyte"                    
    # [11] "Epithelial cells"               "Esophageal epithelial cell"    
    # [13] "Esophageal smooth muscle"       "Fibroblast"                    
    # [15] "Follicular cell of thyroid"     "Gastric chief cell"            
    # [17] "Gastric parietal cell"          "Goblet cell"                   
    # [19] "Granulosa cell"                 "Hepatocyte"                    
    # [21] "Keratinocyte"                   "Luminal epithelial cell"       
    # [23] "Macrophage"                     "Mast cell"                     
    # [25] "Mesothelial cell"               "Miscellaneous stromal cell"    
    # [27] "Myofibroblast"                  "Myosatellite cell"             
    # [29] "Neuroendocrine cell"            "Pheumocyte"                    
    # [31] "Schwann cell"                   "Skeletal myocyte"              
    # [33] "Smooth muscle"                  "T-cell"                        
    # [35] "Vascular smooth muscle"



    enrich$cell<-factor(enrich$cell,levels=c("B-cell","T-cell","Macrophage","Monocyte/Granulocyte","NK","Acinar cell of pancreas","Adipocyte","Airway goblet cell","Basal cell","Cardiomyocyte",
    "Cortical cell of adrenal gland","Ductal cell","Endothelial cell","Enterocyte","Epithelial cells",
    "Esophageal epithelial cell","Esophageal smooth muscle","Fibroblast","Follicular cell of thyroid",
    "Gastric chief cell","Gastric parietal cell","Goblet cell","Granulosa cell","Hepatocyte","Keratinocyte",
    "Luminal epithelial cell","Mast cell","Mesothelial cell","Miscellaneous stromal cell",
    "Myofibroblast","Myosatellite cell","Neuroendocrine cell","Pheumocyte","Schwann cell","Skeletal myocyte",
    "Smooth muscle","Vascular smooth muscle"))

    enrich$class<-"cis CRE"
    enrich_cre<-enrich
    dim(enrich_cre)
    # [1] 54 14


    #######################
    # ATAC-SEQ


    files_atac<-files[grep("ATAC",files)]

    data<-gsub("_merged_samples_b37_psc_summary_stats.results","",files_atac)
    data<-gsub(paste(pheno[ph],"_summary_stats_munged_baseline_v1.2_ATAC_counts_",sep=""),"ATAC_",data)

    ##########
    # immune cell atlas:

    for (i in 1:length(files_atac)) {

        tmp<-fread(paste("/path/to/ibdgwas/IIBDGC/post_imputation/analysis/metaanalysis/annotation/ldscore/test/iibdgc_b37/enrichment/",files_atac[i],sep=""),head=T)
        tmp$Enrichment_p[which(is.na(tmp$Enrichment_p))]<-1
        tmp$qvalue<-qvalue(tmp$Enrichment_p)$qvalue
        tmp<-tmp[which(tmp$Category=="L2_1",)]
        tmp$Category<-data[i]

        if(i==1) {
            enrich<-tmp
        } else {
            enrich<-rbind(enrich,tmp)
        }
    }

    enrich$cell_type<-gsub("ATAC_","",enrich$Category)
    enrich$cell_type<-gsub("_U$","",enrich$cell_type)
    enrich$cell_type<-gsub("_S$","",enrich$cell_type)

    enrich$cell_type<-factor(enrich$cell_type,levels = (c("pDCs","Myeloid_DCs","Monocytes",
                                                    "Immature_NK","Mature_NK","Memory_NK",
                                                    "Gamma_delta_T",
                                                    "Effector_CD4pos_T",
                                                    "Naive_Teffs","Memory_Teffs","Th1_precursors","Th2_precursors","Th17_precursors","Follicular_T_Helper",
                                                    "Naive_Tregs","Memory_Tregs",
                                                    "CD8pos_T","Naive_CD8_T","Central_memory_CD8pos_T","Effector_memory_CD8pos_T",
                                                    "Regulatory_T",
                                                    "Bulk_B","Naive_B","Mem_B","Plasmablasts")))



    enrich$cell<-"T-cell"
    enrich$cell[enrich$cell_type %in% c("Mature_NK","Mature_NK","Memory_NK","Immature_NK")]<-"NK"
    enrich$cell[enrich$cell_type %in% c("Mem_B","Plasmablasts","Naive_B","Bulk_B")]<-"B-cell"
    enrich$cell[enrich$cell_type %in% c("Monocytes","Myeloid_DCs","pDCs")]<-"Monocyte/Granulocyte"

    enrich$cell<-factor(enrich$cell,levels=c("B-cell","T-cell","Macrophage","Monocyte/Granulocyte","NK","Acinar cell of pancreas","Adipocyte","Airway goblet cell","Basal cell","Cardiomyocyte",
    "Cortical cell of adrenal gland","Ductal cell","Endothelial cell","Enterocyte","Epithelial cells",
    "Esophageal epithelial cell","Esophageal smooth muscle","Fibroblast","Follicular cell of thyroid",
    "Gastric chief cell","Gastric parietal cell","Goblet cell","Granulosa cell","Hepatocyte","Keratinocyte",
    "Luminal epithelial cell","Mast cell","Mesothelial cell","Miscellaneous stromal cell",
    "Myofibroblast","Myosatellite cell","Neuroendocrine cell","Pheumocyte","Schwann cell","Skeletal myocyte",
    "Smooth muscle","Vascular smooth muscle"))


    enrich$class<-"ATAC-seq"
    enrich_atac<-enrich
    dim(enrich_atac)
    # [1] 45 14

    enrich<-rbind(enrich_atac,enrich_cre)

    col<-carto_pal(12, "Safe")
    col_ligth_1<-sample(lighten(col, 0.90, space = "HCL"))
    col_ligth_2<-sample(lighten(col, 0.60, space = "HCL"))
    col_ligth_3<-sample(lighten(col[1], 0.30, space = "HCL"))

    colo<-c(col,col_ligth_1,col_ligth_2,col_ligth_3,"darkgrey")

    # add fdr 
    # enrich$qvalue<-qvalue(enrich$Enrichment_p)$qvalues

    bonferroni_qvalue<-0.05/nrow(enrich)

    # p1<-ggplot(enrich[which(enrich$class=="ATAC-seq"),], aes(x=Category, y=Enrichment,fill=cell,label = ifelse(qvalue < bonferroni_qvalue, "*", ""))) + 
    # geom_bar(stat="identity", color="black", position=position_dodge()) +  geom_text(vjust = 0,hjust=1.3) +
    # scale_fill_manual(limits = levels(enrich$cell),values=c(colo)) +
    # geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.2,
    #                 position=position_dodge(.9)) + theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
    #                 axis.text.x=element_text(angle=45,hjust=1),plot.margin = unit(c(1,2,1,1), "cm"))

    # p2<-ggplot(enrich[which(enrich$class=="cis CRE"),], aes(x=Category, y=Enrichment,fill=cell,label = ifelse(qvalue < bonferroni_qvalue, "*", ""))) + 
    # geom_bar(stat="identity", color="black", position=position_dodge()) +  geom_text(vjust = 0,hjust=1.3) +
    # scale_fill_manual(limits = levels(enrich$cell),values=c(colo)) +
    # geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.2,
    #                 position=position_dodge(.9)) + theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
    #                 axis.text.x=element_text(angle=45,hjust=1),legend.position="none",plot.margin = unit(c(1,2,1,2), "cm"))

    
    # p3<-ggarrange(p1,p2,heights=c(2,1.5),nrow=2,common.legend=T,legend="bottom")
    # p3

    # ggsave(
    # paste("~/git/IIBDGC_GWAS/plots/",pheno[ph],"_enrichment_analyses.png",sep=""),
    # p3,
    # width = 10,
    # height = 12,
    # dpi = 600,
    # units = c("in"),
    # limitsize = FALSE
    # )


    # create a index by level of enrichemnt:
    tmp<-enrich[,c("Enrichment","Category")]
    tmp<-tmp[order(tmp$Enrichment,decreasing=T),]


    enrich$Category<-as.factor(enrich$Category)
    enrich$Category<-factor(enrich$Category,levels=tmp$Category)


    p1<-ggplot(enrich[which(enrich$class=="ATAC-seq"),], aes(x=Category, y=Enrichment,fill=cell,label = ifelse(qvalue < bonferroni_qvalue, "*", ""))) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +  geom_text(vjust = 0,hjust=0.9,size=7) +
    scale_fill_manual(limits = levels(enrich$cell),values=c(colo)) +
    geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.2,
                    position=position_dodge(.9)) + theme(axis.title.x = element_blank(),
                    axis.text.x=element_text(angle=45,hjust=1),plot.margin = unit(c(1,1,1,4), "cm")) + labs(y="Enrichment") + 
                    ggtitle(paste("Enrichement",toupper(pheno[ph]),"variants in ATAC-seq regions"))


    assign(paste("p1",pheno[ph],sep="_"),p1)

    p2<-ggplot(enrich[which(enrich$class=="cis CRE"),], aes(x=Category, y=Enrichment,fill=cell,label = ifelse(qvalue < bonferroni_qvalue, "*", ""))) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +  geom_text(vjust = 0,hjust=0.9,size=7) +
    scale_fill_manual(limits = levels(enrich$cell),values=c(colo)) +
    geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.2,
                    position=position_dodge(.9)) + theme(axis.title.x = element_blank(),
                    axis.text.x=element_text(angle=45,hjust=1),legend.position="none",plot.margin = unit(c(1,1,1,4), "cm")) + labs(y="Enrichment") + 
                    ggtitle(paste("Enrichement",toupper(pheno[ph]),"variants in cCRE regions"))

    assign(paste("p2",pheno[ph],sep="_"),p2)

    p3<-ggarrange(p1,p2,heights=c(2,1.5),nrow=2,common.legend=T,legend="bottom")
    p3

    ggsave(
    paste("~/git/IIBDGC_GWAS/plots/",pheno[ph],"_enrichment_analyses_sorted.png",sep=""),
    p3,
    width = 12,
    height = 12,
    dpi = 600,
    units = c("in"),
    limitsize = FALSE
    )

    assign(paste("enrich",pheno[ph],sep="_"),enrich)

    rm(enrich,enrich_cre,enrich_atac,files_cre,files_atac)

}

p1<-ggarrange(p1_cd,p1_uc,nrow=2,common.legend=T,legend="right")

ggsave(
    paste("~/git/IIBDGC_GWAS/plots/ATACseq_immune_cell_atlas_enrichment_analyses_sorted.png",sep=""),
    p1,
    width = 15,
    height = 12,
    dpi = 600,
    units = c("in"),
    limitsize = FALSE
    )

p2<-ggarrange(p2_cd,p2_uc,nrow=2,common.legend=T,legend="right")

ggsave(
    paste("~/git/IIBDGC_GWAS/plots/cCRE_immune_cell_atlas_enrichment_analyses_sorted.png",sep=""),
    p2,
    width = 15,
    height = 12,
    dpi = 600,
    units = c("in"),
    limitsize = FALSE
    )


# change plot so that only the associated are significant:

### CD

enrich_cd$cell<-factor(enrich_cd$cell,levels=c("B-cell","T-cell","Macrophage","Monocyte/Granulocyte","NK","Acinar cell of pancreas","Adipocyte","Airway goblet cell","Basal cell","Cardiomyocyte",
    "Cortical cell of adrenal gland","Ductal cell","Endothelial cell","Enterocyte","Epithelial cells",
    "Esophageal epithelial cell","Esophageal smooth muscle","Fibroblast","Follicular cell of thyroid",
    "Gastric chief cell","Gastric parietal cell","Goblet cell","Granulosa cell","Hepatocyte","Keratinocyte",
    "Luminal epithelial cell","Mast cell","Mesothelial cell","Miscellaneous stromal cell",
    "Myofibroblast","Myosatellite cell","Neuroendocrine cell","Pheumocyte","Schwann cell","Skeletal myocyte",
    "Smooth muscle","Vascular smooth muscle","n.s."))

enrich_cd$cell_significant<-enrich_cd$cell
enrich_cd$cell_significant[which(enrich_cd$qvalue > bonferroni_qvalue)]<-"n.s."

table(as.character(enrich_cd$cell_significant))
        #       B-cell           Fibroblast           Macrophage 
        #            8                    1                    3 
        #    Mast cell Monocyte/Granulocyte                   NK 
        #            1                    4                    4 
        #         n.s.        Smooth muscle               T-cell 
        #           45                    1                   32 


# enrich_cd$cell_significant<-factor(enrich_cd$cell_significant,
# levels=c("B-cell","T-cell","Macrophage","Monocyte/Granulocyte","NK",
# "Mast cell","Fibroblast","Smooth muscle","n.s."))

### UC

enrich_uc$cell<-factor(enrich_uc$cell,levels=c("B-cell","T-cell","Macrophage","Monocyte/Granulocyte","NK","Acinar cell of pancreas","Adipocyte","Airway goblet cell","Basal cell","Cardiomyocyte",
    "Cortical cell of adrenal gland","Ductal cell","Endothelial cell","Enterocyte","Epithelial cells",
    "Esophageal epithelial cell","Esophageal smooth muscle","Fibroblast","Follicular cell of thyroid",
    "Gastric chief cell","Gastric parietal cell","Goblet cell","Granulosa cell","Hepatocyte","Keratinocyte",
    "Luminal epithelial cell","Mast cell","Mesothelial cell","Miscellaneous stromal cell",
    "Myofibroblast","Myosatellite cell","Neuroendocrine cell","Pheumocyte","Schwann cell","Skeletal myocyte",
    "Smooth muscle","Vascular smooth muscle","n.s."))

enrich_uc$cell_significant<-enrich_uc$cell
enrich_uc$cell_significant[which(enrich_uc$qvalue > bonferroni_qvalue)]<-"n.s."

table(as.character(enrich_uc$cell_significant))
#               B-cell           Enterocyte           Fibroblast 
#                    8                    3                    1 
#          Goblet cell           Macrophage            Mast cell 
#                    2                    3                    1 
# Monocyte/Granulocyte                   NK                 n.s. 
#                    4                    4                   40 
#        Smooth muscle               T-cell 
#                    1                   32 

# enrich_uc$cell_significant<-factor(enrich_uc$cell_significant,
# levels=c("B-cell","T-cell","Macrophage","Monocyte/Granulocyte","NK",
# "Mast cell","Fibroblast","Smooth muscle","Enterocyte","Goblet cell","n.s."))

ph=1

p1_cd<-ggplot(enrich_cd[which(enrich_cd$class=="ATAC-seq"),], aes(x=Category, y=Enrichment,fill=cell_significant)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +  
    # geom_text(vjust = 0,hjust=0.9,size=7) +
    scale_fill_manual(breaks = levels(enrich_cd$cell_significant),values=c(colo)) +
    geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.2,
                    position=position_dodge(.9)) + theme(axis.title.x = element_blank(),
                    axis.text.x=element_text(angle=45,hjust=1),plot.margin = unit(c(1,1,1,4), "cm"),
                    legend.text = element_text(size=12),legend.title = element_text(size=12)) + labs(y="Enrichment") + 
                    ggtitle(paste("Enrichement",toupper(pheno[ph]),"variants in ATAC-seq regions")) + labs(fill = "Significantly enriched\ncell group")

p2_cd<-ggplot(enrich_cd[which(enrich_cd$class=="cis CRE"),], aes(x=Category, y=Enrichment,fill=cell_significant)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +  
    # geom_text(vjust = 0,hjust=0.9,size=7) +
    scale_fill_manual(breaks = levels(enrich_cd$cell_significant),values=c(colo)) +
    geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.2,
                    position=position_dodge(.9)) + theme(axis.title.x = element_blank(),
                    axis.text.x=element_text(angle=45,hjust=1),plot.margin = unit(c(1,1,1,4), "cm"),
                    legend.text = element_text(size=12),legend.title = element_text(size=12)) + labs(y="Enrichment") + 
                    ggtitle(paste("Enrichement",toupper(pheno[ph]),"variants in cCRE regions")) + labs(fill = "Significantly enriched\ncell group")


ph=2

p1_uc<-ggplot(enrich_uc[which(enrich_uc$class=="ATAC-seq"),], aes(x=Category, y=Enrichment,fill=cell_significant)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +  
    # geom_text(vjust = 0,hjust=0.9,size=7) +
    scale_fill_manual(breaks = levels(enrich_uc$cell_significant),values=c(colo)) +
    geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.2,
                    position=position_dodge(.9)) + theme(axis.title.x = element_blank(),
                    axis.text.x=element_text(angle=45,hjust=1),plot.margin = unit(c(1,1,1,4), "cm"),
                    legend.text = element_text(size=12),legend.title = element_text(size=12)) + labs(y="Enrichment") + 
                    ggtitle(paste("Enrichement",toupper(pheno[ph]),"variants in ATAC-seq regions")) + labs(fill = "Significantly enriched\ncell group")

p2_uc<-ggplot(enrich_uc[which(enrich_uc$class=="cis CRE"),], aes(x=Category, y=Enrichment,fill=cell_significant)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +  
    # geom_text(vjust = 0,hjust=0.9,size=7) +
    scale_fill_manual(breaks = levels(enrich_uc$cell_significant),values=c(colo)) +
    geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.2,
                    position=position_dodge(.9)) + theme(axis.title.x = element_blank(),
                    axis.text.x=element_text(angle=45,hjust=1),plot.margin = unit(c(1,1,1,4), "cm"),
                    legend.text = element_text(size=12),legend.title = element_text(size=12)) + labs(y="Enrichment") + 
                    ggtitle(paste("Enrichement",toupper(pheno[ph]),"variants in cCRE regions")) + labs(fill = "Significantly enriched\ncell group")



p1<-ggarrange(p1_cd,p1_uc,nrow=2,common.legend=T,legend="right")

ggsave(
    paste("~/git/IIBDGC_GWAS/plots/ATACseq_immune_cell_atlas_enrichment_analyses_sorted.png",sep=""),
    p1,
    width = 15,
    height = 12,
    dpi = 600,
    units = c("in"),
    limitsize = FALSE
    )

p2<-ggarrange(p2_cd,p2_uc,nrow=2,common.legend=F,legend="right")

ggsave(
    paste("~/git/IIBDGC_GWAS/plots/cCRE_immune_cell_atlas_enrichment_analyses_sorted.png",sep=""),
    p2,
    width = 15,
    height = 12,
    dpi = 600,
    units = c("in"),
    limitsize = FALSE
    )



enrich_uc$cell_significant<-as.factor(as.character(enrich_uc$cell_significant))
enrich_cd$cell_significant<-as.factor(as.character(enrich_cd$cell_significant))

enrich_uc$cell_significant<-factor(enrich_uc$cell_significant,levels=c("B-cell","T-cell","Macrophage","Monocyte/Granulocyte","Mast cell","NK","Fibroblast","Enterocyte","Goblet cell","Smooth muscle"))
enrich_cd$cell_significant<-factor(enrich_cd$cell_significant,levels=c("B-cell","T-cell","Macrophage","Monocyte/Granulocyte","Mast cell","NK","Fibroblast","Enterocyte","Goblet cell","Smooth muscle"))

ph=1
p2_cd_signif<-ggplot(enrich_cd[which(enrich_cd$class=="cis CRE" & (enrich_cd$qvalue <= bonferroni_qvalue)),], aes(x=Category, y=Enrichment,fill=cell_significant)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +  
    # geom_text(vjust = 0,hjust=0.9,size=7) +
    scale_fill_manual(breaks = levels(enrich_cd$cell_significant),values=c(colo)) +
    geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.2,
                    position=position_dodge(.9)) + theme(axis.title.x = element_blank(),
                    axis.text.x=element_text(angle=45,hjust=1),plot.margin = unit(c(1,1,1,4), "cm"),
                    legend.text = element_text(size=12),legend.title = element_text(size=12)) + labs(y="Enrichment") + 
                    ggtitle(paste("Enrichement",toupper(pheno[ph]),"variants in cCRE regions")) + labs(fill = "Significantly enriched\ncell group")


ph=2
p2_uc_signif<-ggplot(enrich_uc[which(enrich_uc$class=="cis CRE" & (enrich_uc$qvalue < bonferroni_qvalue)),], aes(x=Category, y=Enrichment,fill=cell_significant)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +  
    # geom_text(vjust = 0,hjust=0.9,size=7) +
    scale_fill_manual(breaks = levels(enrich_uc$cell_significant),values=c(colo)) +
    geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.2,
                    position=position_dodge(.9)) + theme(axis.title.x = element_blank(),
                    axis.text.x=element_text(angle=45,hjust=1),plot.margin = unit(c(1,1,1,4), "cm"),
                    legend.text = element_text(size=12),legend.title = element_text(size=12)) + labs(y="Enrichment") + 
                    ggtitle(paste("Enrichement",toupper(pheno[ph]),"variants in cCRE regions")) + labs(fill = "Significantly enriched\ncell group")



p2<-ggarrange(p2_cd_signif,p2_uc_signif,nrow=2,widths=c(0.8,1),common.legend=F,legend="right")

ggsave(
    paste("~/git/IIBDGC_GWAS/plots/cCRE_immune_cell_atlas_enrichment_analyses_sorted_cd_signifcant.png",sep=""),
    p2_cd_signif,
    width = 11,
    height = 8,
    dpi = 600,
    units = c("in"),
    limitsize = FALSE
    )


ggsave(
    paste("~/git/IIBDGC_GWAS/plots/cCRE_immune_cell_atlas_enrichment_analyses_sorted_uc_signifcant.png",sep=""),
    p2_uc_signif,
    width = 16,
    height = 8,
    dpi = 600,
    units = c("in"),
    limitsize = FALSE
    )
