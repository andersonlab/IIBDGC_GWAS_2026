# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# # singularity exec iibdgc_postprocess_10_singularity.sif
# # singularity exec iibdgc_postprocess_10_singularity.sif
# path_gwas="/path/to/ibdgwas/IIBDGC/"

# to test
# MEM=15000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q normal R 

library(LCV)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

args = commandArgs(trailingOnly=TRUE)

dataset<-args[2]
pheno<-args[1]

# # to test model
# pheno<-"cd"
# dataset<-"GCST90029070"

print(dataset)
print(pheno)

#Load trait 1 data and calculate Zs
d1<-read.table(paste(path_gwas,"resources/gwas_summary_statistics/",pheno,"_eur_tier1_list_variants_only_SNPs_rate_0.5_nohla_sumstats_munged.sumstats",sep=""),header=TRUE,sep="\t",stringsAsFactors = FALSE)

#Load trait 2 data and calculate Zs
d2<-read.table(paste(path_gwas,"resources/gwas_summary_statistics/",dataset,"_only_SNPs_sumstats_munged.sumstats",sep=""),header=TRUE,sep="\t",stringsAsFactors = FALSE)

list_snsp<-intersect(d1$SNP,d2$SNP)

#Load LD scores
for (chr in c(1:22)) {
    tmp<-read.table(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/regression_weight_files/weights_only_SNPs.april2025.",chr,".l2.ldscore.gz",sep=""),header=TRUE,sep='\t',stringsAsFactors=FALSE)
    tmp<-tmp[which(tmp$SNP %in% list_snsp),]
    if(chr==1) {
        d3<-tmp
    }else{
        d3<-rbind(d3,tmp)
    }
}

#Merge
data<-merge(d3,d1,by="SNP")
data<-merge(data,d2,by="SNP")

#Sort by position 
data<-data[order(data[,"CHR"],data[,"BP"]),]

#Flip sign of one z-score if opposite alleles-shouldn't occur with UKB data
#If not using munged data, will have to check that alleles match-not just whether they're opposite A1/A2
mismatch<-which(data$A1.x!=data$A1.y,arr.ind=TRUE)
data[mismatch,]$Z.y = data[mismatch,]$Z.y*-1
data[mismatch,]$A1.y = data[mismatch,]$A1.x
data[mismatch,]$A2.y = data[mismatch,]$A2.x


#Run LCV-need to setwd to directory containing LCV package - patch before HGI sets up as module
source("/path/to/software/username/LCV-master/R/RunLCV.R")
setwd("/path/to/software/username/LCV-master/R/")

LCV = RunLCV(data$L2,data$Z.x,data$Z.y)

# LCV1 = RunLCV(data$L2,data$Z.y,data$Z.x)
# sprintf("Estimated posterior gcp=%.2f(%.2f), log10(p)=%.1f; estimated rho=%.2f(%.2f)",LCV$gcp.pm, LCV$gcp.pse, log(LCV$pval.gcpzero.2tailed)/log(10), LCV$rho.est, LCV$rho.err)

output<-as.data.frame(t(unlist(LCV)))
output$d1<-pheno
output$d2<-dataset

file_out<-paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_april2025/genetic_correlation/lcv/",pheno,"_",dataset,"_lcv_output.tsv",sep="")
write.table(output,file_out,col.names=T,row.names=F,quote=F,sep="\t")

system(paste("gzip ",file_out))

q("no")


# p<-ggplot(data,aes(x=Z.x, y=Z.y))+
#  geom_point()
# pdf("~/git/IIBDGC_GWAS/plots/paper_figures/test.pdf")
# p
# dev.off()