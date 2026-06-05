# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=15000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G ibdgwas R


library(data.table)
rm(list=ls())
path_gwas<-"/path/to/ibdgwas/IIBDGC/"

var<-fread("/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_eur_tier_1_metaanalysis_rate_0.5_sorted.bed")


pheno<-c("ibd")

for (i in 1:length(pheno)) {
  
  print(pheno[i])
  
  for (chr in 1:22) {
    
    tmp1<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr,"_",pheno[i],"_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""),head=T)
    tmp1<-tmp1[which(tmp1$MarkerName %in% var$V4),]
    
    if(chr==1) {
      tmp<-tmp1
    } else {
      tmp<-rbind(tmp,tmp1)
    }
    
    rm(tmp1)
    
  }
  
  if(i==1) {
    tmp<-tmp[,c("MarkerName","Position_b38","A1","A2","BETA","SE","P-value","Direction_ed","HetISq",
                "HetChiSq","HetDf","HetPVal","rate_total_sample","avgA2FREQ","avgA2FREQ_CASES","avgA2FREQ_CONTROLS")]
    colnames(tmp)[2:4]<-paste(colnames(tmp)[2:4],"eur_tier_2",sep="_")
    colnames(tmp)[5:ncol(tmp)]<-paste(colnames(tmp)[5:ncol(tmp)],pheno[i],"eur_tier_2",sep="_")
    all<-tmp

  } else {
    tmp<-tmp[,c("MarkerName","BETA","SE","P-value","Direction_ed","HetISq",
                "HetChiSq","HetDf","HetPVal","rate_total_sample","avgA2FREQ","avgA2FREQ_CASES","avgA2FREQ_CONTROLS")]
    colnames(tmp)[2:ncol(tmp)]<-paste(colnames(tmp)[2:ncol(tmp)],pheno[i],"eur_tier_2",sep="_")
    all<-merge(all,tmp,by="MarkerName",sort=F)
  }
  rm(tmp)
}

dim(all[which(all$rate_total_sample_ibd_eur_tier_2>=0.5),])
[1] 12225856       16
> 12225856/dim(var)
[1] 7.795384e-01 3.056464e+06
> dim(var)
[1] 15683456        4