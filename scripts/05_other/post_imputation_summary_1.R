# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# get imputation summary per chromosome

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(dplyr)
library(stringr)
library(viridis)
library(ggExtra)

rm(list=ls())

path<-"/path/to/ibdgwas/IIBDGC/"
args<-commandArgs()

study<-args[6]
print(study)

i<-args[7]
print(i)


# ancestry<-c("eur","noneur")

ancestry<-c("eur")

for (j in 1:length(ancestry)) {
  
  print(ancestry[j])
  
  if(i==23) {
    tmp<-fread(paste(path,"imputed/",study,"/2022/",ancestry[j],"/chrX.info.gz",sep=""),head=T)
  } else {
    tmp<-fread(paste(path,"imputed/",study,"/2022/",ancestry[j],"/chr",i,".info.gz",sep=""),head=T)
  }

  tmp<-as.data.frame(tmp)
  
  tmp<-tmp[which(tmp$Genotyped!="Typed_Only"),]
  tmp$Rsq<-as.numeric(tmp$Rsq)
    
  tmp$chr<-i
  tmp$chr<-as.factor(as.character(tmp$chr))
  
  tmp$Rsq_2<-as.factor(as.character(round(tmp$Rsq,digits=1)))
  tmp$MAF_2<-as.factor(as.character(round(tmp$MAF,digits=1)))
    
  ##################################################
  
  xx<-str_split(tmp$SNP, ":",simplify=T)
  xx<-as.data.frame(xx)
  
  tmp<-cbind(tmp,xx$V2)
  rm(xx)
  colnames(tmp)[ncol(tmp)]<-"Position"
  tmp$Position<-as.numeric(as.character(tmp$Position))
  
  tmp$Rsq<-as.numeric(tmp$Rsq)
  tmp$EmpR<-as.numeric(tmp$EmpR)
  tmp$EmpRsq<-as.numeric(tmp$EmpRsq)
  
  tmp1<-tmp[which(tmp$MAF>0.001),]

  ##################################################
  # plot rainfall imputation rsqr and EmpRsq:
  
  xlimits<-c(min(tmp$Position),max(tmp$Position))
  
  p2 <- ggplot(tmp1, aes(x = Position, y = Rsq,color=MAF)) + ylim(0,1) +
    geom_point(alpha = 0.2,shape=20,size=0.5) +
    geom_rug(alpha = 0.1,colour="black",sides = "b") + theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_blank(),legend.position="bottom",plot.margin = margin(2,0,0,0, "pt"),
          plot.title = element_text(margin = margin(t = 10, b = 5))) +
    scale_x_continuous(labels = function(x) format(x, scientific = FALSE),limits=xlimits) + ggtitle(paste(study," Chr",i," all Variants MAF >= 0.001",sep="")) + scale_color_viridis(option = "D")

  tmp<-tmp[which(tmp$Genotyped=="Genotyped"),]
  
  p1 <- ggplot(tmp, aes(x = Position, y = EmpRsq,color=MAF)) + ylim(0,1) +
    geom_point(alpha = 0.3,shape=20,size=1.5) +
    geom_rug(alpha = 0.1,colour="black",sides = "b") + theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_blank(),legend.position="bottom",plot.margin = margin(2,0,0,0, "pt"),
          plot.title = element_text(margin = margin(t = 10, b = 5))) +
    scale_x_continuous(labels = function(x) format(x, scientific = FALSE),limits=xlimits) + ggtitle(paste(study," Chr",i," Genotyped Variants",sep="")) + scale_color_viridis(option = "D")
  
  pdf(paste(path,"imputed/",study,"/2022/",ancestry[j],"/plots/",study,"_",ancestry[j],"_rainfalls_EmpRsq_Rsqr_chr",i,".pdf",sep=""),width = 14, height = 12,onefile=FALSE)
  print(ggarrange(ggMarginal(p1, type="density",margins = 'y',colour="grey",fill = "grey"),ggMarginal(p2, type="density",margins = 'y',colour="grey",fill = "grey"),nrow=2))
  dev.off()
  
  rm(p1,p2)
  
  print("")
  print("")
  
}

q("no")

