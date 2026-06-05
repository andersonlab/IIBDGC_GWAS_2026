# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
  
cd /path/to/ibdgwas/IIBDGC/post_imputation/analysis/stage_2_summary_statistics/danish/
tar -xvzf Denmark_IBD_CON_20210816.tar.gz


# CHR: chromosome
# POS: genome position
# SNPID: variant ID
# Allele1: Ref allele
# Allele2: Alt allele
# AC_Allele2: allele count of Alt allele
# AF_Allele2: allele frequency of Alt allele
# N: sample size
# BETA: effect size
# SE: standard error of BETA
# Tstat: score statistic
# p.value: p value with SPA applied
# p.value.NA: p value when SPA is not applied
# Is.SPA.converge: whether SPA has converged
# varT: estimated variance of score statistic with sample relatedness incorporated
# varTstar: variance of score statistic without sample relatedness incorporated


# CREATE FILES TO META-ANALYSE:


MEM=8000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1

library(data.table)
path<-"/path/to/ibdgwas/IIBDGC/"

pheno_1<-c("CD","IBD","UC")
pheno_2<-c("cd","ibd","uc")

# list of variants used in our analysis:

for (i in 1:length(pheno_1)) {
  
  print(pheno_1[i])
  dat<-fread(paste(path,"post_imputation/analysis/stage_2_summary_statistics/danish/Denmark_IBD_CON_20210816/",pheno_1[i],"_CON.SAIGE.stats.rsq04.MAF01perc.txt.gz",sep=""),head=T)
  
  # split into chromosomes to meta-analyse with other European data:
  for (chr in 1:22) {
    tmp<-dat[which(dat$CHR==paste("chr",chr,sep="")),]
    tmp$test<-NA
    tmp$CHISQ.Y1<-NA
    tmp$LOG10P.Y1<-(-log10(tmp$p.value))
    tmp<-tmp[,c("CHR","POS","SNPID","Allele1","Allele2","AF_Allele2","imputationInfo","test","BETA","SE","CHISQ.Y1","LOG10P.Y1")]
    colnames(tmp)<-c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","TEST","BETA.Y1","SE.Y1","CHISQ.Y1","LOG10P.Y1")
    
    write.table(tmp,paste(path,"post_imputation/analysis/stage_2_summary_statistics/danish/",pheno_2[i],"/chr",chr,"_danish_gsa_",pheno_2[i],"_eur_sex_10PCs_saige_spa",sep=""),
                col.names=T,row.names=F,quote=F)
  }
}

q("no")


# continue with metal_metaanalysis.R

