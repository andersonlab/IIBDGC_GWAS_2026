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

# pheno<-c("cd","uc","ibd")
# pheno<-c("ibd")
# release<-"eur_eas_sas_tier_2"
# release<-"eur_tier_2"
# release<-"eur_tier_1"

print(release)
print(pheno)


for (i in 1:length(pheno)) {
  
  print(pheno[i])
  
  # delange:
  old<-fread(paste("/path/to/ibdgwas/IIBDGC/summary_files/summary_stats_",pheno[i],"_liu_eur_gwas_immunochip_sorted_noheader_chr_pos_b37_pval_lifted_hg38.gz",sep=""),head=F)
  old<-old[,1:5]
  colnames(old)<-c("chr","position_b38","end","ID_b37","pvalue")
  old<-old[,c("chr","position_b38","ID_b37","pvalue")]


  # add chisq
  old2<-fread(paste0(path_gwas,"summary_files/summary_stats_",pheno[i],"_liu_eur_gwas_immunochip_sorted_noheader_chr_pos_b37_pval"),head=F)
  colnames(old2)<-c("chr","position_b37","end","ID_b37","pvalue","chisq")

  old<-merge(old,old2[,c("ID_b37","chisq")],by="ID_b37")


  old$ref<-gsub("[0-9]{1,2}:[0-9]*_","",old$ID_b37)
  old$alt<-gsub(".*_","",old$ref)
  old$ref<-gsub("_.*","",old$ref)

  old$chr<-gsub("_.*","",old$chr)

  old$MarkerName_1<-paste(old$chr,old$position_b38,old$ref,old$alt,sep=":")
  old$MarkerName_2<-paste(old$chr,old$position_b38,old$alt,old$ref,sep=":")


    table(old$chr)
    #   chr1  chr10  chr11  chr12  chr13  chr14  chr15  chr16  chr17  chr18  chr19 
    # 736648 493122 475085 464218 349248 311343 273654 296661 250138 276141 208255 
    #   chr2  chr20  chr21  chr22   chr3   chr4   chr5   chr6   chr7   chr8   chr9 
    # 803493 216245 130367 127800 682870 714006 616314 629374 560822 536922 416880 

  rm(new)
  for (chr in c(1:22)) {

    print(chr)

    if(release=="eur_tier_1") {
        tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr,"_",pheno[i],"_meta_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""))
    } else {
        tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr,"_",pheno[i],"_meta_",release,"_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""))
    }

    tmp1<-merge(old,tmp,by.x="MarkerName_1",by.y="MarkerName")
    tmp2<-merge(old,tmp,by.x="MarkerName_2",by.y="MarkerName")

    tmp<-rbind(tmp1,tmp2)
    rm(tmp1,tmp2)

    print(nrow(tmp))

    if (chr==1) {
        new<-tmp
    } else {
        new<-rbind(new,tmp)
    }
    rm(tmp)
  }

    colnames(new)[15]<-"pvalue_new"
    new$chisq_new<-(new$BETA/new$SE)^2
   
    print(dim(new))
    # [1] 69466      29
    print(dim(old))
    # [1] 126081       9

    rm(old)
    
    new$pvalue<-as.numeric(new$pvalue)
    new$pvalue_new<-as.numeric(new$pvalue_new)

    new<-new[which(new$pvalue_new!=0),]

    maxlim<-ceiling(max(-log10(new$pvalue),-log10(new$pvalue_new)))

}

new<-as.data.frame(new)

p1<-ggplot(new[which(new$avgA2FREQ_CASES>=0.01 & new$avgA2FREQ_CASES<=0.99 & new$pvalue<0.01),], aes(x=-log10(pvalue), y=(-log10(pvalue_new)))) + 
              geom_point(shape=20) + xlim(0,maxlim) + ylim(0,maxlim) + 
              geom_abline(intercept = 0, slope = 1, color="red",linetype="dashed", size=0.5) +
              ggtitle(paste(toupper(pheno[i]),"P-value variants MAF 0.01, liu pval < 0.01",pheno[i])) + geom_smooth(method='lm', color="darkgrey",linetype="dashed", size=0.8)

ggsave(paste0("~/git/IIBDGC_GWAS/plots/metaanalysis/pvalues_liu_eur_gwas_immunochip_vs_iibdgc_",release,"_",pheno[i],"_maf_0.01_liu_pval_0.1.pdf",sep=""),
    p1,
    width = 10,
    height = 10,
    dpi = 200,
    units = c("in"),
    limitsize = T
)

maxlim<-max(new$chisq,new$chisq_new)

new_significant<-new[which(new$pvalue_new<5E-8),]
mod<-lm(new_significant$chisq_new ~ new_significant$chisq)
coef(mod)[2]
# new_significant$chisq 
#             0.9143975

p1<-ggplot(new_significant, aes(x=chisq, y=chisq_new)) + 
              geom_point(shape=20) + xlim(0,maxlim) + ylim(0,maxlim) + 
              geom_abline(intercept = 0, slope = 1, color="red",linetype="dashed", size=0.5) + 
              ggtitle(paste(toupper(pheno[i]),"Chisq signif variants new iibdgc",release,pheno[i],"vs liu eur gwas immunochip\nslope chisq_new ~ chisq_old",coef(mod)[2])) + geom_smooth(method='lm', color="darkgrey",linetype="dashed", size=0.8)

ggsave(paste0("~/git/IIBDGC_GWAS/plots/metaanalysis/chisq_liu_eur_gwas_immunochip_vs_iibdgc_",release,"_",pheno[i],"_signif.pdf",sep=""),
    p1,
    width = 10,
    height = 10,
    dpi = 200,
    units = c("in"),
    limitsize = T
)

q("no")