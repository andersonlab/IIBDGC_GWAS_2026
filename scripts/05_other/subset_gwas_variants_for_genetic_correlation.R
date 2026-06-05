# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# double check overlap of sample ids:

# # singularity exec iibdgc_postprocess_10_singularity.sif

# MEM=10000
# bsub -J"gsmr_pqtl" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 2 \
# -o ${path_gwas}post_imputation/2022/log/subset_gwas_variants_${gwas_id[i]}_stdout \
# -e ${path_gwas}post_imputation/2022/log/subset_gwas_variants_${gwas_id[i]}_stderr \
# "/software/R-4.3.1/bin/Rscript ~/git/IIBDGC_GWAS/scripts/other/subset_gwas_variants_for_genetic_correlation.R ${gwas_id[i]} > \
# ${path_gwas}post_imputation/2022/log/subset_gwas_variants_for_genetic_correlation_${gwas_id[i]}.Rout"

# MEM=20000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
library(R.utils)
library(ggplot2)

rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

args = commandArgs(trailingOnly=TRUE)
study_id<-args[1]

# for testing purposes
# study_id<-"GCST90480502"
# study_id<-"GCST90002316"
# study_id<-"GCST90301994"
# study_id<-"GCST004030"

print(study_id)

# list of good quality SNPS (info>=0.9 and Neff>0.5 in EUR tier 2) in IIBDGC
ibd<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/allchr_ibd_cd_uc_eur_tier2_list_variants_rate_0.5_info_0.9_only_SNPs.tsv.gz",sep=""),head=T)

ot<-fread(paste(path_gwas,"resources/gwas_summary_statistics/",study_id,"_edited.tsv.gz",sep=""),head=T)
dim(ot)

# correct instances where there is systematically a large number of missmatches
if (nrow(ot[which(ot$SNP %in% ibd$SNP),])<=(nrow(ot)*0.5)) {

      ot1<-ot[which(ot$SNP %in% ibd$SNP),]
      ot$SNP2<-paste0("chr",ot$CHR,":",ot$BP,":",ot$A1,":",ot$A0)

      ot2<-ot[which(ot$SNP2 %in% ibd$SNP),]
      ot2<-ot2[which(!ot2$SNP2 %in% ot1$SNP),]
      rm(ot)

      ot2$A1FREQ2<-(1-(ot2$A1FREQ))
      ot2$BETA2<-(-1*(ot2$BETA))

      ot2<-ot2[,c("SNP2","CHR","BP","A1","A0","A1FREQ2","BETA2","SE","PVALUE","N")]
      colnames(ot2)<-c("SNP","CHR","BP","A0","A1","A1FREQ","BETA","SE","PVALUE","N")

      ot1<-ot1[,c("SNP","CHR","BP","A0","A1","A1FREQ","BETA","SE","PVALUE","N")]

      ot<-rbind(ot1,ot2)
      rm(ot1,ot2)

}

nrow(ot[which(ot$SNP %in% ibd$SNP),])

# ibd$X<-paste(ibd$CHR,ibd$BP,sep="_")
# ot$X<-paste(ot$CHR,ot$BP,sep="_")
# dim(ibd[which(ibd$X %in% ot$X),])
# head(ot[which(ot$X %in% ibd$X),])

all<-merge(ibd,ot,by="SNP")
rm(ibd)
# plot(all$A1FREQ.x,all$A1FREQ.y)

# retain only the list of SNPs pre selected in:
# IIBDGC_GWAS/scripts/qc_imputation_analysis_step_4_onwards/step_45_compute_genetic_correlation.R
list_snps<-fread(paste(path_gwas,"resources/gwas_summary_statistics/subset_allchr_ibd_cd_uc_eur_tier2_variants_for_genetic_correlation.snplist",sep=""),head=T)
all<-all[which(all$SNP %in% list_snps$SNP),]
rm(list_snps)

all$value<-((all$A1FREQ.x-all$A1FREQ.y)^2)/
      ((all$A1FREQ.x+all$A1FREQ.y)*
         (2-all$A1FREQ.x-all$A1FREQ.y))


all$test<-"no"
all$test[which(all$value>=0.05)]<-"yes"

print(table(all$test))
#      no     yes 
# 9175844    7722 

table(all$CHR,all$test)

all<-all[which(all$test=="no"),]
table(all$ALLELE2==all$A0)
table(all$ALLELE1==all$A1)


ot<-ot[which(ot$SNP %in% all$SNP),]
print(dim(ot))

# if no allelle frequency - use freq from IBD:
if (nrow(ot[which(is.na(ot$A1FREQ)),])==nrow(ot)) {
      all<-all[match(ot$SNP,all$SNP),]
      print(table(all$SNP==ot$SNP))
      ot$A1FREQ<-all$A1FREQ.x
}

rm(all)

# final check to assess that no duplicates remain in the results 
ids_dup<-ot$SNP[duplicated(ot$SNP)]
ot<-ot[which(!ot$SNP %in% ids_dup),]

ot$INFO<-1

# rename columns to same format as IIBDGC input files:
ot<-ot[,c("SNP","CHR","BP","A0","A1","BETA","SE","PVALUE","INFO","N","A1FREQ")]
colnames(ot)<-c("SNP","CHR","BP","ALLELE2","ALLELE1","BETA","SE","P-value","INFO","N","A1FREQ")

fwrite(ot,paste(path_gwas,"resources/gwas_summary_statistics/",study_id,"_edited_2.tsv.gz",sep=""),
col.names=T,row.names=F,quote=F,sep="\t")

q("no")
