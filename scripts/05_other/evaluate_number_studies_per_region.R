# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

### /software/bin/R-4.3.1

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(stringr)

path<-"/path/to/ibdgwas/IIBDGC/"

all<-read.table(paste(path,"post_imputation/analysis/metaanalysis/summary_results/list_genome_wide_significant_variants_variants_forward_regression_known_signals_plus_old_no_overlaps.bed",sep=""),head=F,sep="\t",check.names=F)
all$class_region<-NA
all$class_region[grep("_",all$V4)]<-"new"
table(all$class_region)
# new 
# 79
# 2 new merged into 2 old in between
all$class_region[grep(":",all$V4)]<-"old"
table(all$class_region)
# new old 
# 40 212 

all$class_region[grep(",",all$V4)]<-"old_extended"
table(all$class_region)
# new          old old_extended 
# 40          173           39 

keep<-fread(paste(path,"post_imputation/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis",sep=""),head=F)
# note that this file includes also variants with high heterogeneity

colnames(all)[1:3]<-c("chr","start","end")


dat<-matrix(ncol=5,nrow=nrow(all))
dat<-as.data.frame(dat)
colnames(dat)<-c("region","N","N_intersection","N_hetfilter","N_interset_hetfilter")
dat$chr<-all$chr

pheno<-c("cd","uc")

for(j in 1:length(pheno)) {
  
  print(pheno[j])
  
  for(chr in 1:22) {
    
    print(chr)
    
    tmp<-fread(paste(path,"post_imputation/analysis/metaanalysis/",pheno[j],"/",chr,"_meta_noGC_PCs_firthse_1_formatted_A2_effect_and_alternative.txt",sep=""))
    
    tmp$chr<-paste("chr",chr,sep="")
    tmp$pos<-gsub("chr[0-9]{1,2}:","",tmp$MarkerName)
    tmp$pos<-as.numeric(gsub(":[A-Z]*:[A-Z]*","",tmp$pos))
    tmp$pheno<-pheno[j]
    
    all_tmp<-all[which(all$chr==chr),]
    dat_tmp<-dat[which(dat$chr==chr),]
    
    for (ii in 1:nrow(all_tmp)) {
      
      dat_tmp$region[ii]<-paste(chr,all_tmp$start[ii],all_tmp$end[ii],sep="_")
      
      tmp1<-tmp[which(tmp$chr==paste("chr",all_tmp$chr[ii],sep="") & tmp$pos>=all_tmp$start[ii] & tmp$pos<=all_tmp$end[ii] ),]
      
      dat_tmp$N[ii]<-nrow(tmp1)
      
      dat_tmp$N_hetfilter[ii]<-nrow(tmp1[which(tmp1$HetPVal>0.001),])
      
      dat_tmp$N_intersection[ii]<-nrow(tmp1[which( (tmp1$pheno=="ibd" & tmp1$HetDf==9) | (tmp1$pheno=="cd" & tmp1$HetDf==7) | (tmp1$pheno=="uc" & tmp1$HetDf==8)),])
      
      tmp1<-tmp1[which( (tmp1$pheno=="ibd" & tmp1$HetDf==9) | (tmp1$pheno=="cd" & tmp1$HetDf==7) | (tmp1$pheno=="uc" & tmp1$HetDf==8)),]
      tmp1<-tmp1[which(tmp1$HetPVal>0.001),]
      
      dat_tmp$N_interset_hetfilter[ii]<-nrow(tmp1)
      
    }
    
    if(chr==1) {
      dat_final<-dat_tmp
    } else {
      dat_final<-rbind(dat_final,dat_tmp)
    }
    
  }
  
  assign(pheno[j],dat_final)
  
}

cd1<-cd
uc1<-uc
rm(cd,uc)


#### WHAT IF WE KEEP VARIANS IN N STUDIES > N TOTAL STUDIES -1, BEING SURE those are not GSA or hce

for (j in 1:length(pheno)) {
  
  print(pheno[j])
  
  for(chr in 1:22) {
    
    print(chr)
    
    tmp<-fread(paste(path,"post_imputation/analysis/metaanalysis/",pheno[j],"/",chr,"_meta_noGC_PCs_firthse_1_formatted_A2_effect_and_alternative.txt",sep=""))
    
    tmp$chr<-paste("chr",chr,sep="")
    tmp$pos<-gsub("chr[0-9]{1,2}:","",tmp$MarkerName)
    tmp$pos<-as.numeric(gsub(":[A-Z]*:[A-Z]*","",tmp$pos))
    tmp$pheno<-pheno[j]
    
    all_tmp<-all[which(all$chr==chr),]
    dat_tmp<-dat[which(dat$chr==chr),]
    
    gsa<-read.table(paste(path,"post_imputation/analysis/regenie/gsa/",pheno[j],"/chr",chr,"_gsa_step2_",pheno[j],"_eur_sex_PCs_firthse.regenie",sep=""),head=T)
    humancoreexome<-read.table(paste(path,"post_imputation/analysis/regenie/humancoreexome/",pheno[j],"/chr",chr,"_humancoreexome_step2_",pheno[j],"_eur_sex_PCs_firthse.regenie",sep=""),head=T)
    
    for (ii in 1:nrow(all_tmp)) {
      
      dat_tmp$region[ii]<-paste(chr,all_tmp$start[ii],all_tmp$end[ii],sep="_")
      tmp1<-tmp[which(tmp$chr==paste("chr",all_tmp$chr[ii],sep="") & tmp$pos>=all_tmp$start[ii] & tmp$pos<=all_tmp$end[ii] ),]
      
      if(pheno[j]=="cd") {
        hetdf<-7
      }
      if(pheno[j]=="uc") {
        hetdf<-8
      }
      
      # to exclude any variant in < all-1 studies
      list_exclude1<-tmp1[which(tmp1$HetDf<(hetdf-1)),"MarkerName"]
      
      # but also any variant not in gsa or Humancoreexome
      list_exclude2<-tmp1[which(!tmp1$MarkerName %in% gsa$ID),"MarkerName"]
      
      # but also any variant not in gsa or Humancoreexome
      list_exclude3<-tmp1[which(!tmp1$MarkerName %in% humancoreexome$ID),"MarkerName"]
      
      list_exclude<-rbind(list_exclude1,list_exclude2,list_exclude3)
      list_exclude<-list_exclude[!duplicated(list_exclude$MarkerName),]
      
      tmp1<-tmp1[which(!tmp1$MarkerName %in% list_exclude$MarkerName),]
      
      tmp1<-tmp1[which(tmp1$HetPVal>0.001),]
      
      dat_tmp$N_all_minus_one_study_hetfilter[ii]<-nrow(tmp1)
      
    }
    
    if(chr==1) {
      dat_final<-dat_tmp
    } else {
      dat_final<-rbind(dat_final,dat_tmp)
    }
    
  }
  
  assign(pheno[j],dat_final)
  
}

cd<-cd[,c(1,7)]
uc<-uc[,c(1,7)]

colnames(cd1)[2:ncol(cd1)]<-paste("cd",colnames(cd1)[2:ncol(cd1)],sep="_")
colnames(uc1)[2:ncol(uc1)]<-paste("uc",colnames(uc1)[2:ncol(uc1)],sep="_")
colnames(cd)[2:ncol(cd)]<-paste("cd",colnames(cd)[2:ncol(cd)],sep="_")
colnames(uc)[2:ncol(uc)]<-paste("uc",colnames(uc)[2:ncol(uc)],sep="_")

all<-merge(cd1,cd,by="region",all.x=T,sort=F)
all<-merge(all,uc1,by="region",all.x=T,sort=F)
all<-merge(all,uc,by="region",all.x=T,sort=F)


write.table(all,paste(path,"post_imputation/analysis/metaanalysis/summary_results/summary_number_variants_per_regions_252_using_differenent_filters.tsv",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")

############################################################################################################################################################
############################################################################################################################################################

##### plot:

### /software/bin/R-4.3.1

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(stringr)

path<-"/path/to/ibdgwas/IIBDGC/"
all<-read.table(paste(path,"post_imputation/analysis/metaanalysis/summary_results/summary_number_variants_per_regions_252_using_differenent_filters.tsv",sep=""),head=T)


# all<-all[order(all$cd_N,decreasing=F),]
# 
# p1<-ggplot(all, aes(x=region, y=cd_N)) +
#   geom_point() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p2<-ggplot(all, aes(x=region, y=uc_N)) +
#   geom_point() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# pdf(paste(path,"post_imputation/analysis/metaanalysis/plots/number_variants_per_region",sep=""),width = 25, height = 10)
# print(ggarrange(p1,p2,nrow=2))
# dev.off()
# 
# system(paste("cp ",path,"post_imputation/analysis/metaanalysis/plots/number_variants_per_region ~/tmp_plots/",sep=""))

#########

library(reshape)
library(khroma)

bright <- colour("bright")
bright(5)
# blue       red     green    yellow      cyan 
# "#4477AA" "#EE6677" "#228833" "#CCBB44" "#66CCEE" 

cd<-all[,c(1:5,7)]
colnames(cd)<-c("region","N","N_intersection","N_hetfilter","N_intersect_hetfilter","N_all_minus_one_study_hetfilter")
uc<-all[,c(1,8:11,13)]
colnames(uc)<-c("region","N","N_intersection","N_hetfilter","N_intersect_hetfilter","N_all_minus_one_study_hetfilter")

# cd$percent_intersection <-cd$cd_N_intersection/cd$cd_N
# cd$percent_hetfilter <-cd$cd_N_hetfilter/cd$cd_N
# cd$percent_interset_hetfilter <-cd$cd_N_interset_hetfilter/cd$cd_N
# cd$percent_all_minus_one_study_hetfilter <-cd$cd_N_all_minus_one_study_hetfilter/cd$cd_N
# cd<-cd[,c(1,7:10)]

cd <- melt(cd, id=c("region"))
cd$pheno<-"CD"
uc <- melt(uc, id=c("region"))
uc$pheno<-"UC"

all1<-rbind(cd,uc)
bright(5)
# blue       red     green    yellow      cyan 
# "#4477AA" "#EE6677" "#228833" "#CCBB44" "#66CCEE"

p1<-ggplot(all1[which(all1$variable %in% c("N") ),], aes(x=region, y=value, color=variable)) + 
  geom_point() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(~pheno,nrow=2) + 
  scale_color_manual(values=c("#4477AA")) + theme(legend.position="top")

pdf(paste(path,"post_imputation/analysis/metaanalysis/plots/number_variants_per_region_all_N",sep=""),width = 28, height = 10)
p1
dev.off()
system(paste("cp ",path,"post_imputation/analysis/metaanalysis/plots/number_variants_per_region* ~/tmp_plots/",sep=""))


p1<-ggplot(all1[which(all1$variable %in% c("N","N_hetfilter") ),], aes(x=region, y=value, color=variable)) + 
  geom_point() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(~pheno,nrow=2) + 
  scale_color_manual(values=c("#4477AA","#228833")) + theme(legend.position="top")

pdf(paste(path,"post_imputation/analysis/metaanalysis/plots/number_variants_per_region_all_N_het",sep=""),width = 28, height = 10)
p1
dev.off()
system(paste("cp ",path,"post_imputation/analysis/metaanalysis/plots/number_variants_per_region* ~/tmp_plots/",sep=""))


#####################

p1<-ggplot(all1[which(all1$variable %in% c("N","N_hetfilter","N_intersect_hetfilter") ),], aes(x=region, y=value, color=variable)) + 
  geom_point() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(~pheno,nrow=2) + 
  scale_color_manual(values=c("#4477AA","#228833","#EE6677")) + theme(legend.position="top")

pdf(paste(path,"post_imputation/analysis/metaanalysis/plots/number_variants_per_region_all_N_het_intersect",sep=""),width = 28, height = 10)
p1
dev.off()
system(paste("cp ",path,"post_imputation/analysis/metaanalysis/plots/number_variants_per_region* ~/tmp_plots/",sep=""))

ids<-all1[which(all1$value<=2500),"region"]
ids<-ids[!duplicated(ids)]


#####################

p1<-ggplot(all1[which(all1$variable %in% c("N","N_hetfilter","N_intersect_hetfilter") & all1$region %in% ids),], aes(x=region, y=value, color=variable)) + 
  geom_point() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(~pheno,nrow=2) + 
  scale_color_manual(values=c("#4477AA","#228833","#EE6677")) + theme(legend.position="top")

pdf(paste(path,"post_imputation/analysis/metaanalysis/plots/number_variants_per_region_all_N_het_intersect_zoom",sep=""),width = 5, height = 10)
p1
dev.off()
system(paste("cp ",path,"post_imputation/analysis/metaanalysis/plots/number_variants_per_region* ~/tmp_plots/",sep=""))


#####################

p1<-ggplot(all1[which(all1$variable %in% c("N","N_hetfilter","N_intersect_hetfilter","N_all_minus_one_study_hetfilter") ),], aes(x=region, y=value, color=variable)) + 
  geom_point() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(~pheno,nrow=2) + 
  scale_color_manual(values=c("#4477AA","#228833","#EE6677","#66CCEE")) + theme(legend.position="top")
pdf(paste(path,"post_imputation/analysis/metaanalysis/plots/number_variants_per_region_all_N_het_intersect_all_minus_one",sep=""),width = 28, height = 10)
p1
dev.off()
system(paste("cp ",path,"post_imputation/analysis/metaanalysis/plots/number_variants_per_region* ~/tmp_plots/",sep=""))

#####################

p1<-ggplot(all1[which(all1$variable %in% c("N","N_hetfilter","N_intersect_hetfilter","N_all_minus_one_study_hetfilter") & all1$region %in% ids),], aes(x=region, y=value, color=variable)) + 
  geom_point() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(~pheno,nrow=2) + 
  scale_color_manual(values=c("#4477AA","#228833","#EE6677","#66CCEE")) + theme(legend.position="top")

pdf(paste(path,"post_imputation/analysis/metaanalysis/plots/number_variants_per_region_all_N_het_intersect_all_minus_one_zoom",sep=""),width = 5, height = 10)
p1
dev.off()
system(paste("cp ",path,"post_imputation/analysis/metaanalysis/plots/number_variants_per_region* ~/tmp_plots/",sep=""))

#####################

p1<-ggplot(all1, aes(x=region, y=value, color=variable)) + 
  geom_point() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(~pheno,nrow=2) 

pdf(paste(path,"post_imputation/analysis/metaanalysis/plots/number_variants_per_region_diff_thresholds",sep=""),width = 28, height = 10)
p1
dev.off()


all1<-all1[grep("^22_",all1$region),]

p1<-ggplot(all1[which(all1$variable %in% c("N") & all1$pheno=="CD" ),], aes(x=region, y=value, color=variable)) + 
  geom_point() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(~pheno,nrow=2) + 
  scale_color_manual(values=c("#4477AA","#228833","#EE6677")) + theme(legend.position="top") + ylim(0,max(all1$value))

pdf(paste(path,"post_imputation/analysis/metaanalysis/plots/number_variants_per_region_all_N_zoom_chr22",sep=""),width = 5, height = 6)
p1
dev.off()
system(paste("cp ",path,"post_imputation/analysis/metaanalysis/plots/number_variants_per_region* ~/tmp_plots/",sep=""))



#####################################################################################################################################################################
#####################################################################################################################################################################

#### /software/bin/R-4.3.1

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(stringr)

path<-"/path/to/ibdgwas/IIBDGC/"

all<-read.table(paste(path,"post_imputation/analysis/metaanalysis/summary_results/list_genome_wide_significant_variants_variants_forward_regression_known_signals_plus_old_no_overlaps.bed",sep=""),head=F,sep="\t",check.names=F)
all$class_region<-NA
all$class_region[grep("_",all$V4)]<-"new"
table(all$class_region)
# new 
# 79
# 2 new merged into 2 old in between
all$class_region[grep(":",all$V4)]<-"old"
table(all$class_region)
# new old 
# 40 212 

all$class_region[grep(",",all$V4)]<-"old_extended"
table(all$class_region)
# new          old old_extended 
# 40          173           39 

keep<-fread(paste(path,"post_imputation/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis",sep=""),head=F)
# note that this file includes also variants with high heterogeneity

colnames(all)[1:3]<-c("chr","start","end")


dat<-matrix(ncol=18,nrow=nrow(all))
dat<-as.data.frame(dat)

colnames(dat)<-c("region","N_1_study","N_2_study","N_3_study","N_4_study","N_5_study","N_6_study","N_7_study","N_8_study","chr",
                 "maf_avg_1_study","maf_avg_2_study","maf_avg_3_study","maf_avg_4_study","maf_avg_5_study","maf_avg_6_study","maf_avg_7_study","maf_avg_8_study")

dat$chr<-all$chr


array<-c("illumina370","illumina550","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")

pheno<-c("cd","uc")
for(j in 1:length(pheno)) {
  
  print(pheno[j])
  
  chr=22
  
  print(chr)
  
  tmp<-fread(paste(path,"post_imputation/analysis/metaanalysis/",pheno[j],"/",chr,"_meta_noGC_PCs_firthse_1_formatted_A2_effect_and_alternative.txt",sep=""))
  
  
  
  for (jj in 1:length(array)) {
    
    file_tmp11<-paste(path,"post_imputation/analysis/regenie/",array[jj],"/",pheno[j],"/chr",chr,"_",array[jj],"_step2_",pheno[j],"_eur_sex_PCs_firthse.regenie",sep="")
    
    if(file.exists(file_tmp11)) {
      
      tmp11<-read.table(file_tmp11,head=T)
      tmp11$MAF<-pmin(tmp11$A1FREQ,1-tmp11$A1FREQ)
      tmp11<-tmp11[,c("ID","MAF")]
      colnames(tmp11)[2]<-paste(colnames(tmp11)[2],array[jj],sep="_")
      
      if(jj==1) {
        tmp1<-tmp11
      } else {
        tmp1<-merge(tmp1,tmp11,by="ID",all=T)
      }
      
      rm(tmp11)
    }
    
  }
  
  tmp1$MAF<-NA
  tmp1$MAF<-rowMeans(tmp1[,2:(ncol(tmp1)-1)],na.rm=T)
  
  tmp<-merge(tmp,tmp1[,c("ID","MAF")],by.x="MarkerName",by.y="ID",all.x=T)
  summary(tmp$MAF)
  
  tmp$chr<-paste("chr",chr,sep="")
  tmp$pos<-gsub("chr[0-9]{1,2}:","",tmp$MarkerName)
  tmp$pos<-as.numeric(gsub(":[A-Z]*:[A-Z]*","",tmp$pos))
  tmp$pheno<-pheno[j]
  
  all_tmp<-all[which(all$chr==chr),]
  dat_tmp<-dat[which(dat$chr==chr),]
  
  
  
  for (ii in 1:nrow(all_tmp)) {
    
    tmp1<-tmp[which(tmp$chr==paste("chr",all_tmp$chr[ii],sep="") & tmp$pos>=all_tmp$start[ii] & tmp$pos<=all_tmp$end[ii] ),]
    
    dat_tmp$region[ii]<-paste(chr,all_tmp$start[ii],all_tmp$end[ii],sep="_")
    
    for (iii in 1:length(table(tmp1$HetDf)) ) {
      dat_tmp[ii,1+iii]<-table(tmp1$HetDf)[iii]
    }
    dat_tmp$maf_avg_1_study[ii]<-mean(tmp1$MAF[which(tmp1$HetDf==0)])
    dat_tmp$maf_avg_2_study[ii]<-mean(tmp1$MAF[which(tmp1$HetDf==1)])
    dat_tmp$maf_avg_3_study[ii]<-mean(tmp1$MAF[which(tmp1$HetDf==2)])
    dat_tmp$maf_avg_4_study[ii]<-mean(tmp1$MAF[which(tmp1$HetDf==3)])
    dat_tmp$maf_avg_5_study[ii]<-mean(tmp1$MAF[which(tmp1$HetDf==4)])
    dat_tmp$maf_avg_6_study[ii]<-mean(tmp1$MAF[which(tmp1$HetDf==5)])
    dat_tmp$maf_avg_7_study[ii]<-mean(tmp1$MAF[which(tmp1$HetDf==6)])
    dat_tmp$maf_avg_8_study[ii]<-mean(tmp1$MAF[which(tmp1$HetDf==7)])
  }
  
  assign(pheno[j],dat_tmp)
  
}

cd1<-cd[,c(1,11:18)]
colnames(cd1)[2:ncol(cd1)]<-seq(1:8)
cd1<-melt(cd1, id=c("region"))

cd2<-cd[,c(1:9)]
colnames(cd2)[2:ncol(cd2)]<-seq(1:8)
cd2<-melt(cd2, id=c("region"))

p1<-ggplot(cd2, aes(x=variable, y=value, group=region, colour =region)) + geom_line() + geom_point() + 
  ylab("N Variants") + xlab("Number studies") + facet_wrap(~pheno,nrow=2)
p2<-ggplot(cd1, aes(x=variable, y=value, group=region, colour =region)) + geom_line() + geom_point() + 
  ylab("Mean MAF") + xlab("Number studies")

pdf(paste(path,"post_imputation/analysis/metaanalysis/plots/number_variants_per_region_number_studies",sep=""),width = 10, height = 10)
ggarrange(p1,p2,nrow=2,common.legend = T)
dev.off()
system(paste("cp ",path,"post_imputation/analysis/metaanalysis/plots/number_variants_per_region* ~/tmp_plots/",sep=""))


#####################################################################################################################################################################
#####################################################################################################################################################################



###########################################################################################

MEM=4000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -m "modern_hardware" -G ibdgwas /software/bin/R-4.3.1

library(data.table)
library(ggpubr)
library(ggplot2)
library(ggExtra)

path<-"/path/to/ibdgwas/IIBDGC/"

pheno<-c("ibd")

# get MAF distribution:


array<-c("illumina370","illumina550","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")


dat<-matrix(nrow=22,ncol=10)
dat<-as.data.frame(dat)

for (j in 1:length(pheno)) {
  
  for (chr in c(1:22)) {
    
    tmp<-fread(paste(path,"post_imputation/analysis/metaanalysis/",pheno[j],"/",chr,"_meta_noGC_PCs_firthse_1.txt",sep=""),head=T)
    
    # tmp$pval_range<-NA
    # tmp$pval_range[which(tmp$HetPVal>=0.1)]<-"more_0.1"
    # tmp$pval_range[which(tmp$HetPVal<0.1 & tmp$HetPVa>=0.01)]<-"1e1_1e2"
    # tmp$pval_range[which(tmp$HetPVal<0.01 & tmp$HetPVa>=0.001)]<-"1e2_1e3"
    # tmp$pval_range[which(tmp$HetPVal<0.001 & tmp$HetPVa>=0.0001)]<-"1e3_1e4"
    # tmp$pval_range[which(tmp$HetPVal<0.0001 & tmp$HetPVa>=0.00001)]<-"1e4_1e5"
    # tmp$pval_range[which(tmp$HetPVal<0.00001 & tmp$HetPVa>=0.000001)]<-"1e5_1e6"
    # tmp$pval_range[which(tmp$HetPVal<0.000001 & tmp$HetPVa>=0.0000001)]<-"1e6_1e7"
    # tmp$pval_range[which(tmp$HetPVal<0.0000001 & tmp$HetPVa>=0.00000001)]<-"1e7_1e8"
    # tmp$pval_range[which(tmp$HetPVal<0.00000001)]<-"less_1e8"
    
    for (jj in 1:length(array)) {
      
      file_tmp11<-paste(path,"post_imputation/analysis/regenie/",array[jj],"/",pheno[j],"/chr",chr,"_",array[jj],"_step2_",pheno[j],"_eur_sex_PCs_firthse.regenie",sep="")
      
      if(file.exists(file_tmp11)) {
        
        tmp11<-fread(file_tmp11,head=T)
        tmp11$MAF<-pmin(tmp11$A1FREQ,1-tmp11$A1FREQ)
        tmp11<-tmp11[,c("ID","MAF")]
        colnames(tmp11)[2]<-paste(colnames(tmp11)[2],array[jj],sep="_")
        
        if(jj==1) {
          tmp1<-tmp11
        } else {
          tmp1<-merge(tmp1,tmp11,by="ID",all=T)
        }
        
        rm(tmp11)
      }
      
    }
    
    tmp1$MAF<-NA
    tmp1$MAF<-rowMeans(tmp1[,2:(ncol(tmp1)-1)],na.rm=T)
    
    tmp<-merge(tmp,tmp1[,c("ID","MAF")],by.x="MarkerName",by.y="ID",all.x=T)
    rm(tmp1)
    
    # dat<-tapply(tmp$MAF, tmp$HetDf, mean)
    # dat<-as.data.frame(t(dat))
    # dat$chr<-chr
    # dat$N_var<-nrow(tmp)
    
    if(chr==1){
      tmp_final<-tmp
    } else {
      tmp_final<-rbind(tmp,tmp_final)
    }
    
    rm(tmp)
  }
  
  assign(pheno[j],tmp_final)
  rm(tmp_final)
}

fwrite(ibd,paste(path,"post_imputation/analysis/metaanalysis/ibd/allchr_meta_noGC_PCs_firthse_1_with_MAF.txt",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")
fwrite(cd,paste(path,"post_imputation/analysis/metaanalysis/cd/allchr_meta_noGC_PCs_firthse_1_with_MAF.txt",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")
fwrite(uc,paste(path,"post_imputation/analysis/metaanalysis/uc/allchr_meta_noGC_PCs_firthse_1_with_MAF.txt",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")

# combine with per sample rate files and save (see file evaluate_percentage_total_sample_size

ibd1<-fread(paste(path,"post_imputation/analysis/metaanalysis/ibd/allchr_meta_noGC_PCs_firthse_1_with_MAF.txt",sep=""),head=T)
ibd2<-fread(paste(path,"post_imputation/analysis/metaanalysis/ibd/allchr_meta_noGC_PCs_firthse_1_with_per_sample_rate.txt",sep=""),head=T)
table(nrow(ibd1)==nrow(ibd2))
table(nrow(ibd)==nrow(ibd1))
ibd<-merge(ibd1,ibd2[,c("MarkerName","rate_total_sample")],by="MarkerName")
table(nrow(ibd)==nrow(ibd1))
fwrite(ibd,paste(path,"post_imputation/analysis/metaanalysis/ibd/allchr_meta_noGC_PCs_firthse_1_with_MAF_with_per_sample_rate.txt",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")
rm(ibd1,ibd2)

cd1<-fread(paste(path,"post_imputation/analysis/metaanalysis/cd/allchr_meta_noGC_PCs_firthse_1_with_MAF.txt",sep=""),head=T)
cd2<-fread(paste(path,"post_imputation/analysis/metaanalysis/cd/allchr_meta_noGC_PCs_firthse_1_with_per_sample_rate.txt",sep=""),head=T)
table(nrow(cd1)==nrow(cd2))
cd<-merge(cd1,cd2[,c("MarkerName","rate_total_sample")],by="MarkerName")
table(nrow(cd)==nrow(cd1))
fwrite(cd,paste(path,"post_imputation/analysis/metaanalysis/cd/allchr_meta_noGC_PCs_firthse_1_with_MAF_with_per_sample_rate.txt",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")
rm(cd1,cd2)

uc1<-fread(paste(path,"post_imputation/analysis/metaanalysis/uc/allchr_meta_noGC_PCs_firthse_1_with_MAF.txt",sep=""),head=T)
uc2<-fread(paste(path,"post_imputation/analysis/metaanalysis/uc/allchr_meta_noGC_PCs_firthse_1_with_per_sample_rate.txt",sep=""),head=T)
table(uc1$MarkerName==uc2$MarkerName)
table(nrow(uc1)==nrow(uc2))
uc<-merge(uc1,uc2[,c("MarkerName","rate_total_sample")],by="MarkerName")
table(nrow(uc)==nrow(uc1))
fwrite(uc,paste(path,"post_imputation/analysis/metaanalysis/uc/allchr_meta_noGC_PCs_firthse_1_with_MAF_with_per_sample_rate.txt",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")
rm(uc1,uc2)


# rm intermediate files

cd1<-tapply(cd$MAF, cd$HetDf, mean)
cd2<-tapply(cd$MAF, cd$HetDf, sd)

dat1<-rbind(cd1,cd2)
dat1<-as.data.frame(t(dat1))
dat1$pheno<-"CD"
colnames(dat1)[1:2]<-c("MAF_mean","MAF_sd")
dat1$N_studies<-seq(1:nrow(dat1))

uc1<-tapply(uc$MAF, uc$HetDf, mean)
uc2<-tapply(uc$MAF, uc$HetDf, sd)

dat2<-rbind(uc1,uc2)
dat2<-as.data.frame(t(dat2))
dat2$pheno<-"UC"
colnames(dat2)[1:2]<-c("MAF_mean","MAF_sd")
dat2$N_studies<-seq(1:nrow(dat2))

all<-rbind(dat1,dat2)

ggplot(df2, aes(x=dose, y=len, group=supp, color=supp)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2,
                position=position_dodge(0.05))



all<-rbind(cd1,uc1)

all$min<-all$MAF_mean-all$MAF_sd
all$min[which(all$min<0)]<-0

all$max<-all$MAF_mean+all$MAF_sd

pdf(paste("~/tmp_plots/mean_maf_per_study.pdf",sep=""),height = 5,width = 12) 
ggplot(all, aes(x=N_studies, y=MAF_mean,color=pheno)) + geom_point() + geom_line() + ylim(-0.001,0.4) +
  scale_x_continuous(breaks =seq(1:9)) +
  xlab("Number Array groups with data") + ylab("Mean MAF") + scale_color_manual(values=c("#3182bd","#de2d26")) + facet_wrap(~pheno,ncol=2) +
  geom_errorbar(aes(ymin=min, ymax=max), width=.2, position=position_dodge(0.05))
dev.off()


# write.table(cd,paste(path,"post_imputation/analysis/metaanalysis/summary_results/maf_per_study_included_bin_per_chr_cd_analysis.tsv",sep=""),
# col.names=T,row.names=F,quote=F,sep="\t")

# write.table(uc,paste(path,"post_imputation/analysis/metaanalysis/summary_results/maf_per_study_included_bin_per_chr_uc_analysis.tsv",sep=""),
#             col.names=T,row.names=F,quote=F,sep="\t")
# 
#     
# cd<-read.table(paste(path,"post_imputation/analysis/metaanalysis/summary_results/maf_per_study_included_bin_per_chr_cd_analysis.tsv",sep=""),head=T)
# uc<-read.table(paste(path,"post_imputation/analysis/metaanalysis/summary_results/maf_per_study_included_bin_per_chr_uc_analysis.tsv",sep=""),head=T)

# dat<-matrix(ncol=2,nrow=8)
# dat<-as.data.frame(dat)
# colnames(dat)<-c("N_studies","mean_MAF")
# dat$N_studies<-seq(1:nrow(dat))
# 
# for(i in 1:nrow(dat)) {
#   
#   tmp<-cd[,c(i,ncol(cd))]
#   colnames(tmp)[1]<-"MAF"
#   
#   dat$mean_MAF[i]<-sum(tmp$MAF*tmp$N_var)/sum(tmp$N_var)
#   
# }
# 
# dat$pheno<-"CD"
# cd1<-dat
# 
# 
# dat<-matrix(ncol=2,nrow=9)
# dat<-as.data.frame(dat)
# colnames(dat)<-c("N_studies","mean_MAF")
# dat$N_studies<-seq(1:nrow(dat))
# 
# for(i in 1:nrow(dat)) {
#   
#   tmp<-uc[,c(i,ncol(uc))]
#   colnames(tmp)[1]<-"MAF"
#   
#   dat$mean_MAF[i]<-sum(tmp$MAF*tmp$N_var)/sum(tmp$N_var)
#   
# }
# 
# dat$pheno<-"UC"
# uc1<-dat

# create bins of MAF:



                              