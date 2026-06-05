# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

###########################################################################################

MEM=4000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -m "modern_hardware" -G ibdgwas /software/bin/R-4.3.1

library(data.table)
library(ggpubr)
library(ggplot2)
library(ggExtra)

path<-"/path/to/ibdgwas/IIBDGC/"


array<-c("illumina370","illumina550","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")
sample_size_ibd<-c(3589,7646,8147,1235,4581,21196,1998,3343,46774,630)
sample_size_cd<-c(2165,7073,0,0,4581,15708,1694,2566,34983,406)
sample_size_uc<-c(2160,6515,8147,1235,0,14955,1657,2225,24482,419)

sample_size<-matrix(ncol=4,nrow=length(array))
sample_size<-as.data.frame(sample_size)
colnames(sample_size)<-c("array","n_ibd","n_cd","n_uc")
sample_size$array<-array
sample_size$n_ibd<-sample_size_ibd
sample_size$n_cd<-sample_size_cd
sample_size$n_uc<-sample_size_uc


pheno<-c("ibd","cd","uc")

for (j in 1:length(pheno)) {
  
  for (chr in c(1:22)) {
    
    tmp<-fread(paste(path,"post_imputation/analysis/metaanalysis/",pheno[j],"/",chr,"_meta_noGC_PCs_firthse_1.txt",sep=""),head=T)
  
    for (jj in 1:length(array)) {
      
      file_tmp11<-paste(path,"post_imputation/analysis/regenie/",array[jj],"/",pheno[j],"/chr",chr,"_",array[jj],"_step2_",pheno[j],"_eur_sex_PCs_firthse.regenie",sep="")
      
      if(file.exists(file_tmp11)) {
        
        tmp11<-fread(file_tmp11,head=T)
        tmp11$sample<-sample_size[jj,j+1]
        tmp11<-tmp11[,c("ID","sample")]
        colnames(tmp11)[2]<-paste(colnames(tmp11)[2],array[jj],sep="_")
        
        if(jj==1) {
          tmp1<-tmp11
        } else {
          tmp1<-merge(tmp1,tmp11,by="ID",all=T)
        }
        
        rm(tmp11)
      }
      
    }
    
    tmp1$total_sample<-rowSums(tmp1[,2:ncol(tmp1)],na.rm=T)
    tmp1$rate_total_sample<-tmp1$total_sample/sum(sample_size[,j+1])
    
    tmp<-merge(tmp,tmp1[,c("ID","rate_total_sample")],by.x="MarkerName",by.y="ID")
    
    if (chr==1) {
      tmp_final<-tmp
    } else {
      tmp_final<-rbind(tmp_final,tmp)
    }
    rm(tmp)
  }
  assign(pheno[j],tmp_final)
  rm(tmp_final)
}

# > dim(ibd)
# [1] 21050046       12
# > dim(cd)
# [1] 20280921       12
# > dim(uc)
# [1] 20638832       12

fwrite(ibd,paste(path,"post_imputation/analysis/metaanalysis/ibd/allchr_meta_noGC_PCs_firthse_1_with_per_sample_rate.txt",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")
fwrite(cd,paste(path,"post_imputation/analysis/metaanalysis/cd/allchr_meta_noGC_PCs_firthse_1_with_per_sample_rate.txt",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")
fwrite(uc,paste(path,"post_imputation/analysis/metaanalysis/uc/allchr_meta_noGC_PCs_firthse_1_with_per_sample_rate.txt",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")

cd<-fread(paste(path,"post_imputation/analysis/metaanalysis/cd/allchr_meta_noGC_PCs_firthse_1_with_per_sample_rate.txt",sep=""),head=T)
uc<-fread(paste(path,"post_imputation/analysis/metaanalysis/uc/allchr_meta_noGC_PCs_firthse_1_with_per_sample_rate.txt",sep=""),head=T)


ibd_tab<-table(cut(ibd$rate_total_sample,  breaks = seq(from=0.0,to=1,by=0.01)))
ibd_tab<-as.data.frame(ibd_tab)
ibd_tab$pheno<-"ibd"
ibd_tab<-ibd_tab[order(ibd_tab$Var1,decreasing=T),]
ibd_tab$cum_sum<-cumsum(ibd_tab$Freq)
ibd_tab$cum_rate_var<-ibd_tab$cum_sum/sum(ibd_tab$Freq)

cd_tab<-table(cut(cd$rate_total_sample, breaks = seq(from=0.0,to=1,by=0.01)))
cd_tab<-as.data.frame(cd_tab)
cd_tab$pheno<-"cd"
cd_tab<-cd_tab[order(cd_tab$Var1,decreasing=T),]
cd_tab$cum_sum<-cumsum(cd_tab$Freq)
cd_tab$cum_rate_var<-cd_tab$cum_sum/sum(cd_tab$Freq)

uc_tab<-table(cut(uc$rate_total_sample, breaks = seq(from=0.0,to=1,by=0.01)))
uc_tab<-as.data.frame(uc_tab)
uc_tab$pheno<-"uc"
uc_tab<-uc_tab[order(uc_tab$Var1,decreasing=T),]
uc_tab$cum_sum<-cumsum(uc_tab$Freq)
uc_tab$cum_rate_var<-uc_tab$cum_sum/sum(uc_tab$Freq)

all<-rbind(ibd_tab,cd_tab,uc_tab)
all$Var1<-factor(all$Var1,levels = rev(levels(all$Var1)))

all<-all[which(all$Var1 %in% c("(0.99,1]","(0.98,0.99]","(0.97,0.98]","(0.96,0.97]","(0.95,0.96]"
                               ,"(0.94,0.95]","(0.93,0.94]","(0.92,0.93]","(0.91,0.92]","(0.9,0.91]","(0.89,0.9]"
                               ,"(0.88,0.89]","(0.87,0.88]","(0.86,0.87]","(0.85,0.86]","(0.84,0.85]","(0.83,0.84]"
                               ,"(0.82,0.83]","(0.81,0.82]","(0.8,0.81]","(0.79,0.8]","(0.78,0.79]","(0.77,0.78]","(0.76,0.77]","(0.75,0.76]")),]


pdf(paste("~/tmp_plots/mean_per_sample_rate_per_study.pdf",sep=""),height = 5,width = 18) 
ggplot(all, aes(x=Var1, y=cum_sum,color=pheno)) + geom_point() + geom_line() + ylim(9000000,10500000) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(labels = scales::comma) +
  xlab("Rate samples with data") + ylab("N Variants") + scale_color_manual(values=c("#3182bd","#a38d2e","#de2d26")) + facet_wrap(~pheno,ncol=3)
dev.off()




   