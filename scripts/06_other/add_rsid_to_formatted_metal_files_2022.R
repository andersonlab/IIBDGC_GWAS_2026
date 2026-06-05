# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# add rsids to formatted metal output files:

library(data.table)
library(ggpubr)
library(ggplot2)
library(ggExtra)
library(stringr)
library(textclean)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

args <- commandArgs()
chr<-args[6]

print(chr)

pheno<-c("ibd","cd","uc")

ref<-fread(paste("grep -v '^#' ",path_gwas,"post_imputation/2022/analysis/stage_2_summary_statistics/iibdgc/chr",chr,"_eur_all_list_chr_position_alleles_with_rsid_dbsnp154",sep=""),
           head=F)

ref$ID<-paste("chr",ref$V1,":",ref$V2,":",ref$V4,":",ref$V5,sep="")
dim(ref[duplicated(ref$ID),])

# there are some duplicated chr:pos:ref:alt, keep both

dup<-ref[duplicated(ref$ID),"ID"]

ref$dbsnp154<-ref$V3

for (i in 1:nrow(dup)) {
  
  tmp1<-ref[which(ref$ID %in% dup[i]),]
  
  ref$dbsnp154[which(ref$ID %in% dup[i])]<-paste(tmp1$V3,collapse="|")
  
  rm(tmp1)
}

ref<-ref[,c("ID","dbsnp154")]
ref<-ref[!duplicated(ref),]


for (i in 1:length(pheno)) {
   
  dat<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr,"_",pheno[i],"_meta_eur_tier_1_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz",sep=""),
             head=T)

  dat<-merge(dat,ref,by.x="MarkerName",by.y="ID",all.x=T,sort=F)
  dat<-dat[,c("MarkerName","dbsnp154","Position_b38","A1","A2","BETA","SE","P-value","Direction_ed",
              "HetISq","HetChiSq","HetDf","HetPVal","INFO","N","avgA2FREQ","N_CASES","avgA2FREQ_CASES","N_CONTROLS",
              "avgA2FREQ_CONTROLS","Neff","rate_total_sample","rate_Neff")]
  
  fwrite(dat,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr,"_",pheno[i],"_meta_eur_tier_1_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff_with_rsid_dbsnp154.txt.gz",sep=""),
              col.names=T,row.names=F,quote=F,sep="\t")
  
}


q("no")


# # merge all regenie output files, and create one unique file per array and phenotype
# pheno=(ibd cd uc)
# 
# for ph in ${pheno[@]}
# do 
# for chr in {1..22}
# do
# mv ${path_gwas}post_imputation/analysis/metaanalysis/${ph}/${chr}_meta_noGC_PCs_firthse_no_illumina550_1_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA1freq.txt \
# ${path_gwas}post_imputation/analysis/metaanalysis/${ph}/${chr}_${ph}_meta_noGC_PCs_firthse_no_illumina550_1_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt 
# done
# done





