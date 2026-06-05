# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# How to submit
# pheno=(ibd cd uc)
# MEM=45000

# for ph in ${pheno[@]}
# do
# bsub -J"format" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 2 \
# -e ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_compare_metal_mrmega_multiancestr_stderr \
# -o ${path_gwas}post_imputation/2022/analysis/metaanalysis/log/${ph}_compare_metal_mrmega_multiancestr_stdout \
# "/software/R-4.3.1/bin/Rscript ~/git/IIBDGC_GWAS/scripts/other/compare_metal_mrmega_multiancestry.R ${ph} \
# > ${path_gwas}scripts/logs/compare_metal_mrmega_multiancestry_${ph}.Rout"
# done

# MEM=45000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1

library(data.table)
library(ggplot2)
library(ggpubr)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

args <- commandArgs()
pheno<-args[6]

print(pheno)

  for (chr in c(1:22)) {

    metal_tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno,"/",chr,"_",pheno,"_meta_eur_eas_sas_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""),head=T)
    metal_tmp<-metal_tmp[,c("MarkerName","BETA","SE","P-value","HetDf","HetPVal","rate_total_sample")]
    metal_tmp<-metal_tmp[which(metal_tmp$"P-value"<=0.01),]

    if (chr==1) {
        metal<-metal_tmp
        rm(metal_tmp)
    } else {
        metal<-rbind(metal,metal_tmp)
        rm(metal_tmp)
    }
  }

  mrmega<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/mrmega/output_files/",pheno,"/mrmega_eur_eas_sas_tier_2_",pheno,"_allchr_pc1.result",sep=""))
  all<-merge(metal,mrmega,by="MarkerName",all.x=T)
  
  all$'P-value_association'<-as.numeric(all$'P-value_association')
  all$'P-value'<-as.numeric(all$'P-value')

  all<-as.data.frame(all[which(is.na(all$Comments)),])
 # remove site

  all$pvalue_ancestry_het<-pchisq(all$chisq_ancestry_het, 1, lower.tail = F)


  p1<-ggplot(all, aes(y= -log10(all$'P-value_association'), x= -log10(all$"P-value"),color=-log10(all$HetPVal))) + geom_point() + geom_smooth(method = "lm", se = FALSE)
  p2<-ggplot(all, aes(y= -log10(all$'P-value_association'), x= -log10(all$'P-value'),color=-log10(abs(all$pvalue_ancestry_het)))) + geom_point() + geom_smooth(method = "lm", se = FALSE)

  pdf(paste("~/git/IIBDGC_GWAS/plots/",pheno,"_pvalue_metal_vs_mrmega.pdf",sep=""))
  ggarrange(p1,p2,legend="bottom")
  dev.off()

q("no")