# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

MEM=18000
bsub -Is -m "modern_hardware" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -q yesterday -G your_hpc_group /software/bin/R-4.3.1


library(data.table)
library(ggpubr)
library(ggExtra)
library(gridExtra)
library(grid)
library(data.table)
library(qqman)
library(colorspace)

args <- commandArgs()
pheno<-args[6]

print(pheno)

path<-"/path/to/ibdgwas/IIBDGC/"

# array<-c("illumina370","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")

array<-c("illumina370","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","illuminaexome")

for (chr in 1:22) {
  
  tmp<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/",pheno,"/chr",chr,"_eur_all_",pheno,"_meta_noGC_PCs_firthse_1.txt",sep=""),head=T)
  tmp<-tmp[which(tmp$HetPVal<=1E-5 & tmp$TotalSampleSize>=max(tmp$TotalSampleSize)*0.8),]
  
  if(chr==1) {
    all<-tmp
  } else {
    all<-rbind(all,tmp)
  }
  rm(tmp)
}

all$pos<-gsub("chr[0-9]{1,2}:","",all$MarkerName)
all$pos<-as.numeric(gsub(":.*","",all$pos))
all$chr<-gsub(":.*","",all$MarkerName)
all$chr<-as.numeric(gsub("chr","",all$chr))


# plot this as manhattan plot


pdf(paste(path,"post_imputation/2022/analysis/metaanalysis/plots/",pheno,"_meta_HetPval_noGC_PCs_firthse.pdf",sep=""),width=14, height=5)
manhattan(all, genomewideline = FALSE, suggestiveline = FALSE, main = toupper(pheno),
          chr = "chr", bp = "pos", p = "HetPVal", snp = "MarkerName")
dev.off()

# extract if one study is driving the het:

for (j in 1:length(array)) {
  
  tmp1<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/",pheno,"/allchr_eur_all_",pheno,"_meta_noGC_PCs_firthse_no_",array[j],"_1.txt",sep=""),head=T)
  tmp1<-tmp1[which(tmp1$MarkerName %in% all$MarkerName),]
  
  tmp1$pos<-gsub("chr[0-9]{1,2}:","",tmp1$MarkerName)
  tmp1$pos<-as.numeric(gsub(":.*","",tmp1$pos))
  tmp1$chr<-gsub(":.*","",tmp1$MarkerName)
  tmp1$chr<-as.numeric(gsub("chr","",tmp1$chr))
  
  
  pdf(paste(path,"post_imputation/2022/analysis/metaanalysis/plots/",pheno,"_meta_HetPval_noGC_PCs_firthse_no",array[j],".pdf",sep=""),width=14, height=5)
  manhattan(tmp1, genomewideline = FALSE, suggestiveline = FALSE, main = paste(toupper(pheno),"HetPval no",array[j]),
            chr = "chr", bp = "pos", p = "HetPVal", snp = "MarkerName")
  dev.off()
  
  
  tmp1<-tmp1[,c("MarkerName","P-value","HetPVal")]
  colnames(tmp1)[2:ncol(tmp1)]<-paste(c("pvalue","HetPVal"),"no",array[j],sep="_")
                                      
  
  tmp2<-fread(paste(path,"post_imputation/2022/analysis/regenie/",array[j],"/",pheno,"/allchr_",array[j],"_eur_all_step2_",pheno,"_eur_sex_PCs_firthse_",pheno,".regenie",sep=""),head=T)
  tmp2<-tmp2[which(tmp2$ID %in% all$MarkerName),]
  
  tmp2<-tmp2[,c("ID","INFO","A1FREQ_CONTROLS","BETA","LOG10P")]
  colnames(tmp2)[2:ncol(tmp2)]<-paste(c("INFO","A1FREQ_CONTROLS","BETA","LOG10P"),array[j],sep="_")

  tmp<-merge(tmp1,tmp2,by.x="MarkerName",by.y="ID")
    
  all<-merge(all,tmp,by="MarkerName")
  
  rm(tmp,tmp1,tmp2)
  
}





library(data.table)
library(qqman)
library(colorspace)


  ylim_max<-1E-270


