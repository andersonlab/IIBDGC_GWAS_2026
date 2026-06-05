# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# MEM=15000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q normal R 

library(data.table)
library(ggplot2)

path_gwas="/path/to/ibdgwas/IIBDGC/"


args = commandArgs(trailingOnly=TRUE)

pheno<-args[1]
release<-args[2]
locus<-args[3]

# pheno<-c("cd","uc","ibd")
# pheno<-c("cd")
# release<-"eur_eas_sas_tier_2"
# release<-"eur_tier_2"
# release<-"eur_tier_1"
# locus<-"5_55648856_56648857"

print(release)
print(pheno)

chr<-gsub("_.*","",locus)
start_locus<-gsub("^[0-9]{1,2}_","",locus)
end_locus<-as.numeric(gsub("^[0-9]*_","",start_locus))
start_locus<-as.numeric(gsub("_[0-9]*","",start_locus))

for (i in 1:length(pheno)) {
  
  print(pheno[i])
  
  # delange:
  old<-fread(paste("/path/to/ibdgwas/IIBDGC/summary_files/summary_stats_",pheno[i],"_delange_sorted_noheader_chr_pos_b37_pval_lifted_hg38.gz",sep=""),head=F)
  colnames(old)<-c("chr","position_b38","end","ID_b37","pvalue")
  old<-old[,c("chr","position_b38","ID_b37","pvalue")]

  old$ref<-gsub("[0-9]{1,2}:[0-9]*_","",old$ID_b37)
  old$alt<-gsub(".*_","",old$ref)
  old$ref<-gsub("_.*","",old$ref)

  old$chr<-gsub("_.*","",old$chr)

  old$MarkerName_1<-paste(old$chr,old$position_b38,old$ref,old$alt,sep=":")
  old$MarkerName_2<-paste(old$chr,old$position_b38,old$alt,old$ref,sep=":")

  old<-as.data.frame(old)
  old<-old[which(old$chr==paste0("chr",chr) & old$position_b38>=start_locus & old$position_b38<=end_locus),]


    print(chr)

    if(release=="eur_tier_1") {
        tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr,"_",pheno[i],"_meta_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""))
    } else {
        tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr,"_",pheno[i],"_meta_",release,"_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""))
    }

    tmp1<-merge(old,tmp,by.x="MarkerName_1",by.y="MarkerName")
    tmp2<-merge(old,tmp,by.x="MarkerName_2",by.y="MarkerName")

    new<-rbind(tmp1,tmp2)
    rm(tmp1,tmp2)

  }

    colnames(new)[14]<-"pvalue_new"
   
    print(dim(new))
    # [1] 9569439      29
    print(dim(old))
    # [1] 9569607       9

    rm(old)
    
    new$pvalue<-as.numeric(new$pvalue)
    new$pvalue_new<-as.numeric(new$pvalue_new)

    new<-new[which(new$pvalue_new!=0),]

    maxlim<-ceiling(max(-log10(new$pvalue),-log10(new$pvalue_new)))

}

all<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_new_and_old_index_gwas_variants_conditional_on_known_index_gwas_fm_2023_with_iibdgc_tier1_tier2_eas_sas_summary_stats_plus_unsupervised_conditional.tsv",sep=""))
all<-as.data.frame(all)
all<-all[which(all$updated_region==locus),]

new$index_var<-"no"
new$index_var[which( (new$MarkerName_1 %in% all$MarkerName) | (new$MarkerName_2 %in% all$MarkerName))]<-"Yes"

if (locus=='5_55648856_56648857') {
    new$index_var[which( (new$MarkerName_1=="chr5:56158115:A:T") | (new$MarkerName_2=="chr5:56158115:A:T"))]<-"Yes"
}


p1<-ggplot(new[which(new$avgA2FREQ_CASES>=0.001 & new$avgA2FREQ_CASES<=0.999),], aes(x=-log10(pvalue), y=(-log10(pvalue_new)), color=index_var)) + 
              geom_point(shape=20) + xlim(0,maxlim) + ylim(0,maxlim) + 
              geom_abline(intercept = 0, slope = 1, color="red",linetype="dashed", size=0.5) +
              ggtitle(paste("P-value variants MAF 0.001",pheno[i])) + geom_smooth(method='lm', color="darkgrey",linetype="dashed", size=0.8) + scale_color_manual(values=c("black","red"))

ggsave(paste0("~/git/IIBDGC_GWAS/plots/metaanalysis/pvalues_delange_vs_iibdgc_",release,"_",pheno[i],"_maf_0.01_delange_pval_0.1_",locus,".pdf",sep=""),
    p1,
    width = 10,
    height = 10,
    dpi = 200,
    units = c("in"),
    limitsize = T
)

q("no")