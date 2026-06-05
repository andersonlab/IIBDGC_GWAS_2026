# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# compare effect sizes between GSA and hce, just CD/UC because difference in N cd/UC could drive diff in ibd

MEM=4000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G ibdgwas /software/bin/R-4.3.1

library(data.table)
library(ggpubr)
path<-"/path/to/ibdgwas/IIBDGC/"
pheno<-c("cd","uc")

rm(all1)
for (i in 1:length(pheno)){

  print(toupper(pheno[i]))
  
  for (chr in 1:22) {
    print(chr)
    
    tmp2<-fread(paste(path,"post_imputation/analysis/regenie/gsa/",pheno[i],"/chr",chr,"_gsa_step2_",pheno[i],"_eur_sex_PCs_firthse.regenie",sep=""),head=T)
    colnames(tmp2)[4:ncol(tmp2)]<-paste(colnames(tmp2)[4:ncol(tmp2)],"gsa",sep="_")
    
    tmp3<-fread(paste(path,"post_imputation/analysis/regenie/humancoreexome/",pheno[i],"/chr",chr,"_humancoreexome_step2_",pheno[i],"_eur_sex_PCs_firthse.regenie",sep=""),head=T)
    colnames(tmp3)[4:ncol(tmp3)]<-paste(colnames(tmp3)[4:ncol(tmp3)],"hce",sep="_")
    
    all<-merge(tmp2[,c("ID","ALLELE0_gsa","ALLELE1_gsa","A1FREQ_gsa","INFO_gsa","BETA.Y1_gsa","LOG10P.Y1_gsa")],
               tmp3[,c("ID","ALLELE0_hce","ALLELE1_hce","A1FREQ_hce","INFO_hce","BETA.Y1_hce","LOG10P.Y1_hce")],by="ID")

    # exclude anything non significant in both studies:
    all<-all[which(all$LOG10P.Y1_gsa >= -log10(0.001) & all$LOG10P.Y1_hce >= -log10(0.001)),]
    all$pheno<-pheno[i]
    
    print(nrow(all))
    print(cor(all$LOG10P.Y1_gsa,all$LOG10P.Y1_hce))
    print(cor(all$BETA.Y1_gsa,all$BETA.Y1_hce))
    
    
    tmp<-read.table(paste(path,"post_imputation/analysis/metaanalysis/",pheno[i],"/",chr,"_meta_noGC_PCs_firthse_1_formatted_A2_effect_and_alternative.txt",sep=""),head=T)
    all<-merge(tmp,all,by.x="MarkerName",by.y="ID",all.y=T,order=F)
    rm(tmp)
    
    if(!exists("all1")){
      all1<-all
    }else{
      all1<-rbind(all1,all)
    }
    rm(all)
  }
  
  assign(pheno[i],all1)
  rm(all1)
}

write.table(cd,paste(path,"post_imputation/analysis/metaanalysis/summary_results/comparison_heterogeneity_cd_gsa_humancoreexome.tsv",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")
write.table(uc,paste(path,"post_imputation/analysis/metaanalysis/summary_results/comparison_heterogeneity_uc_gsa_humancoreexome.tsv",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")

q("no")

#############

### /software/bin/R-4.3.1

library(data.table)
library(ggpubr)
path<-"/path/to/ibdgwas/IIBDGC/"

cd<-read.table(paste(path,"post_imputation/analysis/metaanalysis/summary_results/comparison_heterogeneity_cd_gsa_humancoreexome.tsv",sep=""),head=T)
uc<-read.table(paste(path,"post_imputation/analysis/metaanalysis/summary_results/comparison_heterogeneity_uc_gsa_humancoreexome.tsv",sep=""),head=T)

nrow(cd)
# [1] 12754
nrow(uc)
# [1] 8439


print(cor(cd$LOG10P.Y1_gsa,cd$LOG10P.Y1_hce))
# [1] 0.8104906
print(cor(cd$BETA.Y1_gsa,cd$BETA.Y1_hce))
# [1] 0.9602837

print(cor(uc$LOG10P.Y1_gsa,uc$LOG10P.Y1_hce))
# [1] 0.8170538
print(cor(uc$BETA.Y1_gsa,uc$BETA.Y1_hce))
# [1] 0.9482833


summary(uc$HetDf,useNA="ifany")
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6       8       8       8       8       8 
summary(cd$HetDf,useNA="ifany")
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.000   7.000   7.000   6.993   7.000   7.000 

######################

cd$Heterogeneticy_effec<-NA
cd$Heterogeneticy_effec[which(cd$HetPVal>0.001)]<-">0.001"
cd$Heterogeneticy_effec[which(cd$HetPVal<0.001)]<-"<0.001"
cd$Heterogeneticy_effec[which(cd$HetPVal<0.0001)]<-"<0.0001"
cd$Heterogeneticy_effec[which(cd$HetPVal<0.00001)]<-"<0.00001"
cd$Heterogeneticy_effec[which(cd$HetPVal<0.000001)]<-"<0.000001"
  
p1<-ggscatter(data=cd, 'BETA.Y1_gsa',"BETA.Y1_hce",
              xlim=c(min(cd$BETA.Y1_gsa,cd$BETA.Y1_hce),max(cd$BETA.Y1_gsa,cd$BETA.Y1_hce)),ylim=c(min(cd$BETA.Y1_gsa,cd$BETA.Y1_hce),max(cd$BETA.Y1_gsa,cd$BETA.Y1_hce)),
              # add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "Regenie GSA CD Beta", merge=T,
              ylab = "Regenie GWAS3 CD Beta",color="Heterogeneticy_effec") + scale_color_manual(values=rev(c("#4575b4", "#fee090", "#fdae61","#f46d43","#d73027")))
p2<-ggscatter(data=cd, 'LOG10P.Y1_gsa', "LOG10P.Y1_hce",
              xlim=c(0,max(cd$LOG10P.Y1_gsa,cd$LOG10P.Y1_hce)),ylim=c(0,max(cd$LOG10P.Y1_gsa,cd$LOG10P.Y1_hce)),
              # add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "Regenie GSA CD log10(P)", ylab = "Regenie GWAS3 CD log10(P)",color="Heterogeneticy_effec")+ scale_color_manual(values=rev(c("#4575b4", "#fee090", "#fdae61","#f46d43","#d73027")))

p<-ggarrange(p1,p2,ncol=2,common.legend = T)
pdf(paste("~/tmp_plots/allchr_CD_effect_size_pvalue_gsa_vs_hce.pdf",sep=""),width = 14)
print(annotate_figure(p,top = text_grob("CD p-value < 0.001 in GSA and GWAS3")))
dev.off()

###

p1<-ggscatter(data=cd[which(cd$HetPVal<0.001),], 'BETA.Y1_gsa',"BETA.Y1_hce",
              xlim=c(min(cd$BETA.Y1_gsa,cd$BETA.Y1_hce),max(cd$BETA.Y1_gsa,cd$BETA.Y1_hce)),ylim=c(min(cd$BETA.Y1_gsa,cd$BETA.Y1_hce),max(cd$BETA.Y1_gsa,cd$BETA.Y1_hce)),
              # add = "reg.line", conf.int = TRUE, 
              # cor.coef = TRUE, cor.method = "pearson",
              xlab = "Regenie GSA CD Beta", merge=T,
              ylab = "Regenie GWAS3 CD Beta",color="Heterogeneticy_effec") + 
  scale_color_manual(values=rev(c("#4575b4", "#fee090", "#fdae61","#f46d43","#d73027")))

p2<-ggscatter(data=cd[which(cd$HetPVal<0.001),], 'LOG10P.Y1_gsa', "LOG10P.Y1_hce",
              xlim=c(0,max(cd$LOG10P.Y1_gsa,cd$LOG10P.Y1_hce)),ylim=c(0,max(cd$LOG10P.Y1_gsa,cd$LOG10P.Y1_hce)),
              # add = "reg.line", conf.int = TRUE, 
              # cor.coef = TRUE, cor.method = "pearson",
              xlab = "Regenie GSA CD log10(P)", ylab = "Regenie GWAS3 cd log10(P)",color="Heterogeneticy_effec") + 
  scale_color_manual(values=rev(c("#4575b4", "#fee090", "#fdae61","#f46d43","#d73027")))

p<-ggarrange(p1,p2,ncol=2,common.legend = T)
pdf(paste("~/tmp_plots/allchr_CD_effect_size_pvalue_gsa_vs_hce_hetPval_filter.pdf",sep=""),width = 14)
print(annotate_figure(p,top = text_grob("CD p-value < 0.001 in GSA and GWAS3")))
dev.off()


dim(cd[which(cd$HetPVal<0.001),])
# [1] 240  25

dim(cd[which(cd$HetPVal<0.0001),])
# [1] 168  25

cd[which( cd$HetPVal<0.0000001),c("MarkerName","Direction_ed","HetPVal",
                                  "BETA.Y1_gsa","LOG10P.Y1_gsa","A1FREQ_gsa","INFO_gsa",
                                  "BETA.Y1_hce","LOG10P.Y1_hce","A1FREQ_hce","INFO_hce")]
# variants across 4 chr, get one example per region:

cd[which(cd$MarkerName %in% c("chr3:49588749:T:G","chr5:96897759:C:G","chr10:99514301:T:C","chr22:39312328:T:C")),c("MarkerName","Direction_ed","HetPVal",
                                  "BETA.Y1_gsa","LOG10P.Y1_gsa","A1FREQ_gsa","INFO_gsa",
                                  "BETA.Y1_hce","LOG10P.Y1_hce","A1FREQ_hce","INFO_hce")]

# MarkerName Direction_ed   HetPVal BETA.Y1_gsa LOG10P.Y1_gsa
# 2660   chr3:49588749:T:G     -------- 1.175e-12  -0.0686843       3.89966
# 5168   chr5:96897759:C:G     +------- 1.485e-24  -0.0806001       4.82479
# 10066 chr10:99514301:T:C     +------- 2.917e-08  -0.1273200      11.73970
# 12681 chr22:39312328:T:C     ---+---- 3.831e-08  -0.0866397       5.60291
# A1FREQ_gsa INFO_gsa BETA.Y1_hce LOG10P.Y1_gsa.1 A1FREQ_hce INFO_hce
# 2660    0.503528 1.005990   -0.136785         3.89966   0.520367 0.989985
# 5168    0.364858 0.980572   -0.112663         4.82479   0.347316 0.960733
# 10066   0.493363 0.975136   -0.119209        11.73970   0.500973 0.965241
# 12681   0.593435 0.992789   -0.120049         5.60291   0.583044 0.934130

######################

uc$Heterogeneticy_effec<-NA
uc$Heterogeneticy_effec[which(uc$HetPVal>0.001)]<-">0.001"
uc$Heterogeneticy_effec[which(uc$HetPVal<0.001)]<-"<0.001"
uc$Heterogeneticy_effec[which(uc$HetPVal<0.0001)]<-"<0.0001"
uc$Heterogeneticy_effec[which(uc$HetPVal<0.00001)]<-"<0.00001"
uc$Heterogeneticy_effec[which(uc$HetPVal<0.000001)]<-"<0.000001"

levels(as.factor(uc$Heterogeneticy_effec))

p1<-ggscatter(data=uc, 'BETA.Y1_gsa',"BETA.Y1_hce",
              xlim=c(min(uc$BETA.Y1_gsa,uc$BETA.Y1_hce),max(uc$BETA.Y1_gsa,uc$BETA.Y1_hce)),ylim=c(min(uc$BETA.Y1_gsa,uc$BETA.Y1_hce),max(uc$BETA.Y1_gsa,uc$BETA.Y1_hce)),
              # add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "Regenie GSA UC Beta", merge=T,
              ylab = "Regenie GWAS3 UC Beta",color="Heterogeneticy_effec") + scale_color_manual(values=rev(c("#4575b4", "#fee090", "#fdae61","#f46d43","#d73027")))

p2<-ggscatter(data=uc, 'LOG10P.Y1_gsa', "LOG10P.Y1_hce",
              xlim=c(0,max(uc$LOG10P.Y1_gsa,uc$LOG10P.Y1_hce)),ylim=c(0,max(uc$LOG10P.Y1_gsa,uc$LOG10P.Y1_hce)),
              # add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "Regenie GSA UC log10(P)", ylab = "Regenie GWAS3 UC log10(P)",color="Heterogeneticy_effec") + 
              scale_color_manual(values=rev(c("#4575b4", "#fee090", "#fdae61","#f46d43","#d73027")))

p<-ggarrange(p1,p2,ncol=2,common.legend = T)
pdf(paste("~/tmp_plots/allchr_UC_effect_size_pvalue_gsa_vs_hce.pdf",sep=""),width = 14)
print(annotate_figure(p,top = text_grob("UC p-value < 0.001 in GSA and GWAS3")))
dev.off()

####

p1<-ggscatter(data=uc[which(uc$HetPVal<0.001),], 'BETA.Y1_gsa',"BETA.Y1_hce",
              xlim=c(min(uc$BETA.Y1_gsa,uc$BETA.Y1_hce),max(uc$BETA.Y1_gsa,uc$BETA.Y1_hce)),ylim=c(min(uc$BETA.Y1_gsa,uc$BETA.Y1_hce),max(uc$BETA.Y1_gsa,uc$BETA.Y1_hce)),
              # add = "reg.line", conf.int = TRUE, 
              # cor.coef = TRUE, cor.method = "pearson",
              xlab = "Regenie GSA UC Beta", merge=T,
              ylab = "Regenie GWAS3 UC Beta",color="Heterogeneticy_effec") + 
  geom_hline(yintercept=0, linetype="dashed",color = "grey") + geom_vhline(xintercept=0, linetype="dashed",color = "grey") + 
  scale_color_manual(values=rev(c("#4575b4", "#fee090", "#fdae61","#f46d43","#d73027")))

p2<-ggscatter(data=uc[which(uc$HetPVal<0.001),], 'LOG10P.Y1_gsa', "LOG10P.Y1_hce",
              xlim=c(0,max(uc$LOG10P.Y1_gsa,uc$LOG10P.Y1_hce)),ylim=c(0,max(uc$LOG10P.Y1_gsa,uc$LOG10P.Y1_hce)),
              # add = "reg.line", conf.int = TRUE, 
              # cor.coef = TRUE, cor.method = "pearson",
              xlab = "Regenie GSA UC log10(P)", ylab = "Regenie GWAS3 UC log10(P)",color="Heterogeneticy_effec") + 
  scale_color_manual(values=rev(c("#4575b4", "#fee090", "#fdae61","#f46d43","#d73027")))

p<-ggarrange(p1,p2,ncol=2,common.legend = T)
pdf(paste("~/tmp_plots/allchr_UC_effect_size_pvalue_gsa_vs_hce_hetPval_filter.pdf",sep=""),width = 14)
print(annotate_figure(p,top = text_grob("UC p-value < 0.001 in GSA and GWAS3")))
dev.off()

dim(uc[which(uc$HetPVal<0.001),])
# [1] 738  25

dim(uc[which(uc$HetPVal<0.0001),])
# [1] 340  25


uc[which( ((uc$BETA.Y1_gsa>0 & uc$BETA.Y1_hce<0) | (uc$BETA.Y1_gsa<0 & uc$BETA.Y1_hce>0)) & uc$Heterogeneticy_effec==">0.001"),
   c("MarkerName","Direction_ed","HetPVal",
     "BETA.Y1_gsa","LOG10P.Y1_gsa","A1FREQ_gsa","INFO_gsa",
     "BETA.Y1_hce","LOG10P.Y1_hce","A1FREQ_hce","INFO_hce")]

# MarkerName Direction_ed  HetPVal BETA.Y1_gsa LOG10P.Y1_gsa
# 6470   chr6:56452312:G:A    +++-+---- 0.001247  -0.0918517       3.13260
# 8244 chr22:36840586:T:TA    -+++---+- 0.001333  -0.1003570       3.09158
# A1FREQ_gsa INFO_gsa BETA.Y1_hce LOG10P.Y1_hce A1FREQ_hce INFO_hce
# 6470   0.162428 0.983883    0.133358       3.35819   0.159980 0.921039
# 8244   0.146078 0.934620    0.137572       3.06503   0.149173 0.849635



######################


dat$chr<-gsub(":[0-9]*:[A-Z]*:[A-Z]*","",dat$MarkerName)
dat$position<-gsub(":[A-Z]*:[A-Z]*","",dat$MarkerName)
dat$position<-gsub("chr[0-9]{1,2}:","",dat$position)

table(dat$chr)
# 1 12  2 22  3  6  7 
# 23  1  1  1  3  1  1

dat<-as.data.frame(dat)

# extract data from info file first:

studies<-c("all_hce","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"mccauley_gsa","ccfa_gsa","cedars_gsa")

i=1
for (chr in c(1,2,3,6,7,12,22)) {
  tmp<-fread(paste("/path/to/ibdgwas/IIBDGC/imputed/",studies[i],"/eur/chr",chr,".info.gz",sep=""),head=T)
  tmp<-tmp[which(tmp$SNP %in% dat$MarkerName),]
  if(chr==1){
    tmp1<-tmp
  }else{
    tmp1<-rbind(tmp1,tmp)
  }
}

colnames(tmp1)[2:ncol(tmp)]<-paste(colnames(tmp1)[2:ncol(tmp)],studies[i],sep="_")
dat<-merge(dat,tmp1,by.x="MarkerName",by.y="SNP",all.x=T)
rm(tmp,tmp1)

for (i in 2:length(studies)) {
  
  for (chr in c(1,2,3,6,7,12,22)) {
    tmp<-fread(paste("/path/to/ibdgwas/IIBDGC/imputed/",studies[i],"/eur/chr",chr,".info.gz",sep=""),head=T)
    tmp<-tmp[which(tmp$SNP %in% dat$MarkerName),]
    if(chr==1){
      tmp1<-tmp
    }else{
      tmp1<-rbind(tmp1,tmp)
    }
  }
  
  dat$study<-studies[i]
  dat_tmp<-merge(dat,tmp1,by.x="MarkerName",by.y="SNP",all.x=T)
  
  if(i==2){
    dat_2<-dat_tmp
  }else{
    dat_2<-rbind(dat_2,dat_tmp)
  }
  
  
}



p1<-ggplot(dat_2, aes(x = ALT_Frq, y = ALT_Frq_all_hce)) +
  geom_point(aes(color = factor(study))) + labs(x = "Alt Freq GSA",y = "Alt Freq UK-GWAS3") +
  xlim(0, 1) + ylim(0, 1)
p2<-ggplot(dat_2, aes(x = Rsq, y = Rsq_all_hce)) +
  geom_point(aes(color = factor(study))) + labs(x = "Imputation Rsq GSA",y = "Imputation Rsq UK-GWAS3")+
  xlim(0.45, 1) + ylim(0.45, 1)


pdf(paste("~/tmp_plots/allchr_",pheno[i],"_frequency_imputationRsq_gsa_vs_hce_31_variants.pdf",sep=""),width = 12)
ggarrange(p1,p2,common.legend = T,legend="bottom")
dev.off()



# are these variants in known regions:
dat<-dat[,1:16]

noov<-read.table(paste(path,"post_imputation/analysis/metaanalysis/summary_results/list_genome_wide_significant_variants_variants_forward_regression_known_signals_plus_old_no_overlaps.bed",sep=""),head=F,sep="\t",check.names=F)
noov$class_region<-NA
noov$class_region[grep("_",noov$V4)]<-"new"
table(noov$class_region)
# new 
# 79
# 2 new merged into 2 old in between
noov$class_region[grep(":",noov$V4)]<-"old"
table(noov$class_region)
# new old 
# 40 212 

noov$class_region[grep(",",noov$V4)]<-"old_extended"
table(noov$class_region)
# new          old old_extended 
# 40          173           39 

noov$new_region_ids<-paste(noov$V1,":",noov$V2,"-",noov$V3,sep="")


dat$chr<-gsub(":[0-9]*:[A-Z]*:[A-Z]*","",dat$MarkerName)
dat$position<-gsub(":[A-Z]*:[A-Z]*","",dat$MarkerName)
dat$position<-as.numeric(gsub("chr[0-9]{1,2}:","",dat$position))
dat$region<-NA
dat$class_region<-NA

noov$chr<-paste("chr",noov$V1,sep="")


# compare freq data from impute output file and the post-analysis data (use ibd results, larger N)

for (chr in c(1,2,3,6,7,12,22)) {
  
  tmp2<-fread(paste(path,"post_imputation/analysis/regenie/gsa/ibd/chr",chr,"_gsa_step2_ibd_eur_sex_PCs_firthse.regenie",sep=""),head=T)
  colnames(tmp2)[4:ncol(tmp2)]<-paste(colnames(tmp2)[4:ncol(tmp2)],"gsa",sep="_")
  
  tmp3<-fread(paste(path,"post_imputation/analysis/regenie/humancoreexome/ibd/chr",chr,"_humancoreexome_step2_ibd_eur_sex_PCs_firthse.regenie",sep=""),head=T)
  colnames(tmp3)[4:ncol(tmp3)]<-paste(colnames(tmp3)[4:ncol(tmp3)],"hce",sep="_")
  
  tmp<-merge(tmp2[,c("ID","ALLELE0_gsa","ALLELE1_gsa","A1FREQ_gsa","INFO_gsa")],tmp3[,c("ID","ALLELE0_hce","ALLELE1_hce","A1FREQ_hce","INFO_hce")],by="ID")
  
  tmp<-tmp[which(tmp$ID %in% dat$MarkerName),]
  
  if(chr==1){
    tmp1<-tmp
  }else{
    tmp1<-rbind(tmp1,tmp)
  }
}


table(tmp1$ALLELE0_gsa==tmp1$ALLELE0_hce)
# TRUE 
# 31 
table(tmp1$ALLELE1_gsa==tmp1$ALLELE1_hce)
# TRUE 
# 31

tmp1<-as.data.frame(tmp1)

p1<-ggplot(tmp1, aes(x = A1FREQ_gsa, y = A1FREQ_hce)) +
  geom_point() + labs(x = "Alt Freq GSA",y = "Alt Freq UK-GWAS3") +
  xlim(0, 1) + ylim(0, 1)
p2<-ggplot(tmp1, aes(x = INFO_gsa, y = INFO_hce)) +
  geom_point() + labs(x = "Imputation Rsq GSA",y = "Imputation Rsq UK-GWAS3")+
  xlim(0.45, 1.01) + ylim(0.45, 1.01)


pdf(paste("~/tmp_plots/allchr_",pheno[i],"_frequency_imputationRsq_gsa_vs_hce_31_variants_regenie data.pdf",sep=""),width = 12,height = 6)
ggarrange(p1,p2,common.legend = T,legend="bottom")
dev.off()

##########################################################################################################################################################
##########################################################################################################################################################

# get heterogenetity distribution:

# CD/UC

MEM=4000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G ibdgwas /software/bin/R-4.3.1

library(data.table)
library(ggpubr)
library(ggplot2)
library(ggExtra)

path<-"/path/to/ibdgwas/IIBDGC/"

pheno<-c("cd","uc","ibd")
for (i in 1:length(pheno)) {
  for (chr in c(1:22)) {
    
    tmp<-fread(paste(path,"post_imputation/analysis/metaanalysis/",pheno[i],"/",chr,"_meta_noGC_PCs_firthse_1.txt",sep=""),head=T)
    tmp$pval_range<-NA
    tmp$pval_range[which(tmp$HetPVal>=0.1)]<-"more_0.1"
    tmp$pval_range[which(tmp$HetPVal<0.1 & tmp$HetPVa>=0.01)]<-"1e1_1e2"
    tmp$pval_range[which(tmp$HetPVal<0.01 & tmp$HetPVa>=0.001)]<-"1e2_1e3"
    tmp$pval_range[which(tmp$HetPVal<0.001 & tmp$HetPVa>=0.0001)]<-"1e3_1e4"
    tmp$pval_range[which(tmp$HetPVal<0.0001 & tmp$HetPVa>=0.00001)]<-"1e4_1e5"
    tmp$pval_range[which(tmp$HetPVal<0.00001 & tmp$HetPVa>=0.000001)]<-"1e5_1e6"
    tmp$pval_range[which(tmp$HetPVal<0.000001 & tmp$HetPVa>=0.0000001)]<-"1e6_1e7"
    tmp$pval_range[which(tmp$HetPVal<0.0000001 & tmp$HetPVa>=0.00000001)]<-"1e7_1e8"
    tmp$pval_range[which(tmp$HetPVal<0.00000001)]<-"less_1e8"
    
    if(chr==1) {
      x<-table(tmp$HetDf,tmp$pval_range)
    } else {
      x1<-table(tmp$HetDf,tmp$pval_range)
      x<-x+x1
    }
    assign(pheno[i],x)
  }
  
}


cd<-as.data.frame(cd)
cd$Var1<-as.numeric(as.character(cd$Var1))
cd$N_studies<-cd$Var1+1


uc<-as.data.frame(uc)
uc$Var1<-as.numeric(as.character(uc$Var1))
uc$N_studies<-uc$Var1+1

ibd<-as.data.frame(ibd)
ibd$Var1<-as.numeric(as.character(ibd$Var1))
ibd$N_studies<-ibd$Var1+1

cd$pheno<-"cd"
colnames(cd)[3]<-"N_variants"
colnames(cd)[2]<-"HetPvalue"

uc$pheno<-"uc"
colnames(uc)[3]<-"N_variants"
colnames(uc)[2]<-"HetPvalue"

ibd$pheno<-"ibd"
colnames(ibd)[3]<-"N_variants"
colnames(ibd)[2]<-"HetPvalue"

all<-rbind(cd,uc,ibd)

all$N_studies<-as.factor(as.character(all$N_studies))
all$HetPvalue<-factor(all$HetPvalue,levels=rev(c("more_0.1","1e1_1e2","1e2_1e3","1e3_1e4","1e4_1e5","1e5_1e6","1e6_1e7","1e7_1e8","less_1e8")))


# number of variants per hetPval bin

p1<-ggplot(data=all, aes(x=HetPvalue, y=N_variants, fill=pheno)) +
  geom_bar(stat="identity") +
  ylab("Number Variants") + xlab("HetPval") + facet_wrap(~pheno,ncol=3) + 
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("#3182bd","#de2d26","#ffeda0"))

pdf(paste("~/tmp_plots/barplot_cd_uc_ibd_number_variants_per_heterogeneity_pval.pdf",sep=""),height = 5,width = 15)
p1
dev.off()

all<-rbind(cd,uc)
all$N_studies<-as.factor(as.character(all$N_studies))
all$HetPvalue<-factor(all$HetPvalue,levels=rev(c("more_0.1","1e1_1e2","1e2_1e3","1e3_1e4","1e4_1e5","1e5_1e6","1e6_1e7","1e7_1e8","less_1e8")))


p1<-ggplot(data=all, aes(x=HetPvalue, y=N_variants, fill=pheno)) +
  geom_bar(stat="identity") +
  ylab("Number Variants") + xlab("HetPval") + facet_wrap(~pheno,ncol=2)+ 
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("#3182bd","#de2d26"))

pdf(paste("~/tmp_plots/barplot_cd_uc_number_variants_per_heterogeneity_pval.pdf",sep=""),height = 5,width = 12)
p1
dev.off()


# number of variants per N studies bin

all$N_studies<-as.numeric(as.character(all$N_studies))

p1<-ggplot(data=all, aes(x=N_studies, y=N_variants, fill=pheno)) +
  geom_bar(stat="identity") +
  scale_x_continuous(breaks =seq(1:9)) +
  ylab("Number Variants") + xlab("Number Array groups with data") + facet_wrap(~pheno,ncol=2) + 
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("#3182bd","#de2d26"))

pdf(paste("~/tmp_plots/barplot_cd_uc_number_variants_per_studies_available.pdf",sep=""),height = 5,width = 12)
p1
dev.off()


# pdf(paste("~/tmp_plots/matrix_cd_heterogeneity_pval_per_studies_available.pdf",sep=""),height = 8,width = 9)
# ggplot(all[which(all$pheno=="cd"),], aes(HetPvalue,N_studies,  fill = N_variants)) + 
#   geom_tile() + geom_text(aes(HetPvalue,N_studies, label = N_variants), color = "black", size = 4) + 
#   scale_fill_gradient(high = "#3182bd", low="white",) + 
#   theme(legend.position = "none",panel.background = element_blank())+ scale_color_manual(values=c("#3182bd","#de2d26"))
# dev.off()
# 
# pdf(paste("~/tmp_plots/matrix_uc_heterogeneity_pval_per_studies_available.pdf",sep=""),height = 9,width = 9)
# ggplot(all[which(all$pheno=="uc"),], aes(HetPvalue,N_studies,  fill = N_variants)) +
#   geom_tile() + geom_text(aes(HetPvalue,N_studies, label = N_variants), color = "black", size = 4) + 
#   scale_fill_gradient(high = "#de2d26", low="white") +
#   theme(legend.position = "none",panel.background = element_blank())
# dev.off()


# Number dropped by n studies:

head(cd)

cd_cum_nstudies<-matrix(nrow=length(levels(as.factor(cd$Var1))),ncol=2)
cd_cum_nstudies<-as.data.frame(cd_cum_nstudies)
colnames(cd_cum_nstudies)<-c("N_studies","N_variants")
cd_cum_nstudies$N_studies<-levels(as.factor(cd$N_studies))

for (i in 1:nrow(cd_cum_nstudies)) {
  cd_cum_nstudies$N_variants[i]<-sum(cd$N_variants[which(cd$N_studies==cd_cum_nstudies$N_studies[i])])
}

cd_cum_nstudies$N_variants_excluded<-cumsum(cd_cum_nstudies$N_variants)
cd_cum_nstudies$percent_variants_excluded<-cd_cum_nstudies$N_variants_excluded/sum(cd_cum_nstudies$N_variants)


head(uc)

uc_cum_nstudies<-matrix(nrow=length(levels(as.factor(uc$Var1))),ncol=2)
uc_cum_nstudies<-as.data.frame(uc_cum_nstudies)
colnames(uc_cum_nstudies)<-c("N_studies","N_variants")
uc_cum_nstudies$N_studies<-levels(as.factor(uc$N_studies))

for (i in 1:nrow(uc_cum_nstudies)) {
  uc_cum_nstudies$N_variants[i]<-sum(uc$N_variants[which(uc$N_studies==uc_cum_nstudies$N_studies[i])])
}

uc_cum_nstudies$N_variants_excluded<-cumsum(uc_cum_nstudies$N_variants)
uc_cum_nstudies$percent_variants_excluded<-uc_cum_nstudies$N_variants_excluded/sum(uc_cum_nstudies$N_variants)


cd_cum_nstudies$pheno<-"CD"
uc_cum_nstudies$pheno<-"UC"

all<-rbind(cd_cum_nstudies,uc_cum_nstudies)
all$N_studies<-as.numeric(as.character(all$N_studies))

pdf(paste("~/tmp_plots/number_variants_dropped_per_available_study.pdf",sep=""),height = 5,width = 12) 
ggplot(all, aes(x=N_studies, y=percent_variants_excluded,color=pheno)) + geom_point() + geom_line() + ylim(0,1) +
  scale_x_continuous(breaks =seq(1:9)) +
  xlab("Number Array groups with data") + ylab("cummulative rate dropped variants") + scale_color_manual(values=c("#3182bd","#de2d26")) + facet_wrap(~pheno,ncol=2)
dev.off()



# Number dropped by het-pvalue:

head(cd)

cd_cum_nstudies<-matrix(nrow=length(levels(as.factor(cd$HetPvalue))),ncol=2)
cd_cum_nstudies<-as.data.frame(cd_cum_nstudies)
colnames(cd_cum_nstudies)<-c("HetPvalue","N_variants")
cd_cum_nstudies$HetPvalue<-levels(as.factor(cd$HetPvalue))

for (i in 1:nrow(cd_cum_nstudies)) {
  cd_cum_nstudies$N_variants[i]<-sum(cd$N_variants[which(cd$HetPvalue==cd_cum_nstudies$HetPvalue[i])])
}
cd_cum_nstudies$HetPvalue<-as.factor(cd_cum_nstudies$HetPvalue)
cd_cum_nstudies$HetPvalue<-factor(cd_cum_nstudies$HetPvalue,levels=rev(c("more_0.1","1e1_1e2","1e2_1e3","1e3_1e4","1e4_1e5","1e5_1e6","1e6_1e7","1e7_1e8","less_1e8")))
cd_cum_nstudies<-cd_cum_nstudies[order(cd_cum_nstudies$HetPvalue),]

cd_cum_nstudies$N_variants_excluded<-cumsum(cd_cum_nstudies$N_variants)
cd_cum_nstudies$percent_variants_excluded<-cd_cum_nstudies$N_variants_excluded/sum(cd_cum_nstudies$N_variants)

head(uc)

uc_cum_nstudies<-matrix(nrow=length(levels(as.factor(uc$HetPvalue))),ncol=2)
uc_cum_nstudies<-as.data.frame(uc_cum_nstudies)
colnames(uc_cum_nstudies)<-c("HetPvalue","N_variants")
uc_cum_nstudies$HetPvalue<-levels(as.factor(uc$HetPvalue))

for (i in 1:nrow(uc_cum_nstudies)) {
  uc_cum_nstudies$N_variants[i]<-sum(uc$N_variants[which(uc$HetPvalue==uc_cum_nstudies$HetPvalue[i])])
}

uc_cum_nstudies$HetPvalue<-as.factor(uc_cum_nstudies$HetPvalue)
uc_cum_nstudies$HetPvalue<-factor(uc_cum_nstudies$HetPvalue,levels=rev(c("more_0.1","1e1_1e2","1e2_1e3","1e3_1e4","1e4_1e5","1e5_1e6","1e6_1e7","1e7_1e8","less_1e8")))
uc_cum_nstudies<-uc_cum_nstudies[order(uc_cum_nstudies$HetPvalue),]

uc_cum_nstudies$N_variants_excluded<-cumsum(uc_cum_nstudies$N_variants)
uc_cum_nstudies$percent_variants_excluded<-uc_cum_nstudies$N_variants_excluded/sum(uc_cum_nstudies$N_variants)

cd_cum_nstudies$pheno<-"CD"
uc_cum_nstudies$pheno<-"UC"

all<-rbind(cd_cum_nstudies,uc_cum_nstudies)

pdf(paste("~/tmp_plots/number_variants_dropped_per_HetPvalue.pdf",sep=""),height = 5,width = 12) 
ggplot(all, aes(x=HetPvalue, y=percent_variants_excluded,color=pheno)) + geom_point() + geom_line() + ylim(0,0.2)  + 
  scale_x_discrete(labels=c("<1E-8","<1E-7","<1E-6","<1E-5","<1E-4","<1E-3","<1E-2","<1E-1",">1E-1")) +
  xlab("HetPvalue") + ylab("cummulative rate dropped variants") + scale_color_manual(values=c("#3182bd","#de2d26")) + facet_wrap(~pheno,ncol=2)
dev.off()

###########################################################################################################################################################################
###########################################################################################################################################################################

MEM=4000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -m "modern_hardware" -G ibdgwas /software/bin/R-4.3.1

library(data.table)
library(ggpubr)
library(ggplot2)
library(ggExtra)

path<-"/path/to/ibdgwas/IIBDGC/"

ibd<-fread(paste(path,"post_imputation/analysis/metaanalysis/ibd/allchr_meta_noGC_PCs_firthse_1_with_MAF_with_per_sample_rate.txt",sep=""),head=T)
cd<-fread(paste(path,"post_imputation/analysis/metaanalysis/cd/allchr_meta_noGC_PCs_firthse_1_with_MAF_with_per_sample_rate.txt",sep=""),head=T)
uc<-fread(paste(path,"post_imputation/analysis/metaanalysis/uc/allchr_meta_noGC_PCs_firthse_1_with_MAF_with_per_sample_rate.txt",sep=""),head=T)

summary(ibd$HetPVal)
min(ibd$HetPVal)
# [1] 2.812e-169

ibd_tab<-table(cut(ibd$HetPVal,  breaks = c(0,1E-8,1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1,1E-0)))
ibd_tab<-as.data.frame(ibd_tab)
ibd_tab$pheno<-"IBD"
cd_tab<-table(cut(cd$HetPVal,  breaks = c(0,1E-8,1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1,1E-0)))
cd_tab<-as.data.frame(cd_tab)
cd_tab$pheno<-"CD"
uc_tab<-table(cut(uc$HetPVal,  breaks = c(0,1E-8,1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1,1E-0)))
uc_tab<-as.data.frame(uc_tab)
uc_tab$pheno<-"UC"

all<-rbind(ibd_tab,cd_tab,uc_tab)

p1<-ggplot(data=all, aes(x=Var1, y=Freq, fill=pheno)) +
  geom_bar(stat="identity") +
  ylab("Number Variants") + xlab("HetPval") + facet_wrap(~pheno,ncol=3)+ 
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values=c("#3182bd","#a38d2e","#de2d26")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

pdf(paste("~/tmp_plots/barplot_cd_uc_number_variants_per_heterogeneity_pval.pdf",sep=""),height = 5,width = 18)
p1
dev.off()


ibd_tab$cum_var<-cumsum(ibd_tab$Freq)
ibd_tab$percent_variants_excluded<-ibd_tab$cum_var/sum(ibd_tab$Freq)

cd_tab$cum_var<-cumsum(cd_tab$Freq)
cd_tab$percent_variants_excluded<-cd_tab$cum_var/sum(cd_tab$Freq)

uc_tab$cum_var<-cumsum(uc_tab$Freq)
uc_tab$percent_variants_excluded<-uc_tab$cum_var/sum(uc_tab$Freq)

all<-rbind(ibd_tab,cd_tab,uc_tab)


pdf(paste("~/tmp_plots/number_variants_dropped_per_HetPvalue.pdf",sep=""),height = 5,width = 18) 
ggplot(all, aes(x=Var1, y=percent_variants_excluded,color=pheno)) + geom_point() + geom_line() + ylim(0,0.2)  + 
  xlab("HetPvalue") + ylab("Cumulative rate dropped variants") + scale_color_manual(values=c("#3182bd","#a38d2e","#de2d26")) + facet_wrap(~pheno,ncol=3) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
dev.off()

#########################################


cd$HetPval_bin<-cut(cd$HetPVal,breaks = c(0,5E-8,1E-8,1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1,1E-0))
cd$Pval_bin<-cut(cd$'P-value',breaks = c(0,5E-8,1E-8,1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1,1E-0))

uc$HetPval_bin<-cut(uc$HetPVal,breaks = c(0,5E-8,1E-8,1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1,1E-0))
uc$Pval_bin<-cut(uc$'P-value',breaks = c(0,5E-8,1E-8,1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1,1E-0))

ibd$HetPval_bin<-cut(ibd$HetPVal,breaks = c(0,5E-8,1E-8,1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1,1E-0))
ibd$Pval_bin<-cut(ibd$'P-value',breaks = c(0,5E-8,1E-8,1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1,1E-0))

cd_matrix<-table(cd$HetPval_bin,cd$Pval_bin)
ibd_matrix<-table(ibd$HetPval_bin,ibd$Pval_bin)
uc_matrix<-table(uc$HetPval_bin,uc$Pval_bin)


library(reshape2)
library(ggplot2)

cd_matrix_1<-melt(cd_matrix)
colnames(cd_matrix_1)[1:2]<-c("HetPval_bin","Pval_bin")
cd_matrix_1$pheno<-"CD"
cd_matrix_1$rate<-cd_matrix_1$value/sum(cd_matrix_1$value)


uc_matrix_1<-melt(uc_matrix)
colnames(uc_matrix_1)[1:2]<-c("HetPval_bin","Pval_bin")
uc_matrix_1$pheno<-"UC"
uc_matrix_1$rate<-uc_matrix_1$value/sum(uc_matrix_1$value)


ibd_matrix_1<-melt(ibd_matrix)
colnames(ibd_matrix_1)[1:2]<-c("HetPval_bin","Pval_bin")
ibd_matrix_1$pheno<-"IBD"
ibd_matrix_1$rate<-ibd_matrix_1$value/sum(ibd_matrix_1$value)


# all<-rbind(cd_matrix_1,uc_matrix_1,ibd_matrix_1)
# 
# p1<-ggplot(cd_matrix_1, aes(HetPval_bin,Pval_bin,  fill = rate)) +
#   geom_tile() + geom_text(aes(HetPval_bin,Pval_bin, label = rate), color = "black", size = 4) +
#   scale_fill_gradient(high = "#3182bd", low="white",) +
#   theme(legend.position = "none",panel.background = element_blank()) 
# 
# p2<-ggplot(ibd_matrix_1, aes(HetPval_bin,Pval_bin,  fill = rate)) +
#   geom_tile() + geom_text(aes(HetPval_bin,Pval_bin, label = rate), color = "black", size = 4) +
#   scale_fill_gradient(high = "#a38d2e", low="white",) +
#   theme(legend.position = "none",panel.background = element_blank()) 
# 
# p3<-ggplot(uc_matrix_1, aes(HetPval_bin,Pval_bin,  fill = rate)) +
#   geom_tile() + geom_text(aes(HetPval_bin,Pval_bin, label = rate), color = "black", size = 4) +
#   scale_fill_gradient(high = "#de2d26", low="white",) +
#   theme(legend.position = "none",panel.background = element_blank())
# 
# pdf(paste("~/tmp_plots/matrix_cd_heterogeneity_pval_effect_pval.pdf",sep=""),height = 5,width = 15)
# ggarrange(p1,p2,p3,ncol=3)
# dev.off()

##############################################################


MEM=4000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -m "modern_hardware" -G ibdgwas /software/bin/R-4.3.1

library(data.table)
library(ggpubr)
path<-"/path/to/ibdgwas/IIBDGC/"

cd<-read.table(paste(path,"post_imputation/analysis/metaanalysis/summary_results/comparison_heterogeneity_cd_gsa_humancoreexome.tsv",sep=""),head=T)
uc<-read.table(paste(path,"post_imputation/analysis/metaanalysis/summary_results/comparison_heterogeneity_uc_gsa_humancoreexome.tsv",sep=""),head=T)

nrow(cd)
# [1] 12754
nrow(uc)
# [1] 8439

######################

uc$Heterogeneticy_effec<-cut(uc$HetPVal,breaks = c(0,1E-7,1E-6,1E-5,1E-4,1E-3,1E-0))

p1<-ggscatter(data=uc, 'BETA.Y1_gsa',"BETA.Y1_hce",
              xlim=c(min(cd$BETA.Y1_gsa,cd$BETA.Y1_hce),max(cd$BETA.Y1_gsa,cd$BETA.Y1_hce)),ylim=c(min(cd$BETA.Y1_gsa,cd$BETA.Y1_hce),max(cd$BETA.Y1_gsa,cd$BETA.Y1_hce)),
              # add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "Regenie GSA CD Beta", merge=T,
              ylab = "Regenie GWAS3 CD Beta",color="Heterogeneticy_effec") + 
  scale_color_manual(values=rev(c("#4575b4","#ffffb2","#fecc5c","#fd8d3c","#f03b20")))
# scale_color_manual(values=rev(c("#4575b4","#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026")))

p2<-ggscatter(data=uc, 'LOG10P.Y1_gsa', "LOG10P.Y1_hce",
              xlim=c(0,max(cd$LOG10P.Y1_gsa,cd$LOG10P.Y1_hce)),ylim=c(0,max(cd$LOG10P.Y1_gsa,cd$LOG10P.Y1_hce)),
              # add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "Regenie GSA CD log10(P)", ylab = "Regenie GWAS3 CD log10(P)",color="Heterogeneticy_effec") + 
  scale_color_manual(values=rev(c("#4575b4","#ffffb2","#fecc5c","#fd8d3c","#f03b20")))
# scale_color_manual(values=rev(c("#4575b4","#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026")))
p<-ggarrange(p1,p2,ncol=2,common.legend = T)
dev.off()

pdf(paste("~/tmp_plots/allchr_UC_effect_size_pvalue_gsa_vs_hce.pdf",sep=""),width = 14)
print(annotate_figure(p,top = text_grob("UC p-value < 0.001 in GSA and GWAS3")))
dev.off()

uc[which(uc$BETA.Y1_gsa<0 & uc$BETA.Y1_hce>0 & uc$HetPVal>0.001),]

####


p1<-ggscatter(data=uc[which(uc$HetPVal<0.001),], 'BETA.Y1_gsa',"BETA.Y1_hce",
              xlim=c(min(cd$BETA.Y1_gsa,cd$BETA.Y1_hce),max(cd$BETA.Y1_gsa,cd$BETA.Y1_hce)),ylim=c(min(cd$BETA.Y1_gsa,cd$BETA.Y1_hce),max(cd$BETA.Y1_gsa,cd$BETA.Y1_hce)),
              # add = "reg.line", conf.int = TRUE, 
              # cor.coef = TRUE, cor.method = "pearson",
              xlab = "Regenie GSA CD Beta", merge=T,
              ylab = "Regenie GWAS3 CD Beta",color="Heterogeneticy_effec") + 
  scale_color_manual(values=rev(c("#4575b4","#ffffb2","#fecc5c","#fd8d3c","#f03b20")))
# scale_color_manual(values=rev(c("#4575b4","#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026")))

p2<-ggscatter(data=uc[which(uc$HetPVal<0.001),], 'LOG10P.Y1_gsa', "LOG10P.Y1_hce",
              xlim=c(0,max(cd$LOG10P.Y1_gsa,cd$LOG10P.Y1_hce)),ylim=c(0,max(cd$LOG10P.Y1_gsa,cd$LOG10P.Y1_hce)),
              # add = "reg.line", conf.int = TRUE, 
              # cor.coef = TRUE, cor.method = "pearson",
              xlab = "Regenie GSA CD log10(P)", ylab = "Regenie GWAS3 CD log10(P)",color="Heterogeneticy_effec") + 
  scale_color_manual(values=rev(c("#4575b4","#ffffb2","#fecc5c","#fd8d3c","#f03b20")))
# scale_color_manual(values=rev(c("#4575b4","#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026")))
p<-ggarrange(p1,p2,ncol=2,common.legend = T)
dev.off()

              
pdf(paste("~/tmp_plots/allchr_UC_effect_size_pvalue_gsa_vs_hce_hetPval_filter.pdf",sep=""),width = 14)
print(annotate_figure(p,top = text_grob("UC p-value < 0.001 in GSA and GWAS3")))
dev.off()

dim(uc[which(uc$HetPVal<1E-3),])
# [1] 738  25

dim(uc[which(uc$HetPVal<1E-4),])
# [1] 340  25

dim(uc[which(uc$HetPVal<1E-5),])
# [1] 45 25

uc[which( ((uc$BETA.Y1_gsa>0 & uc$BETA.Y1_hce<0) | (uc$BETA.Y1_gsa<0 & uc$BETA.Y1_hce>0)) & uc$Heterogeneticy_effec=="(0.001,1]"),
   c("MarkerName","Direction_ed","HetPVal",
     "BETA.Y1_gsa","LOG10P.Y1_gsa","A1FREQ_gsa","INFO_gsa",
     "BETA.Y1_hce","LOG10P.Y1_hce","A1FREQ_hce","INFO_hce")]


uc[which(uc$MarkerName %in% c("chr6:32279268:T:G","chr10:99524480:T:G","chr11:114453250:A:G")),c("MarkerName","Direction_ed","BETA","P.value","HetPVal",
                                                                                               "BETA.Y1_gsa","LOG10P.Y1_gsa","A1FREQ_gsa","INFO_gsa",
                                                                                               "BETA.Y1_hce","LOG10P.Y1_hce","A1FREQ_hce","INFO_hce")]
# MarkerName Direction_ed    BETA    P.value   HetPVal BETA.Y1_gsa
# 2860   chr6:32279268:T:G    +++++++++  0.9893 2.720e-131 8.868e-06    1.167680
# 7445  chr10:99524480:T:G    --------- -0.1499  1.798e-28 5.738e-06   -0.129345
# 7510 chr11:114453250:A:G    -+------+ -0.0957  8.580e-12 6.412e-09   -0.130352
# LOG10P.Y1_gsa A1FREQ_gsa INFO_gsa BETA.Y1_hce LOG10P.Y1_hce A1FREQ_hce
# 2860      55.63330  0.0201216 0.962635   0.8595040      36.62520   0.037829
# 7445       9.71900  0.4943140 0.992683  -0.1384480       6.40234   0.500335
# 7510       9.22708  0.3716490 1.001730  -0.0968758       3.13879   0.355447
# INFO_hce
# 2860 0.978284
# 7445 0.994067
# 7510 0.987128

array<-c("illumina370","illumina550","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")

pheno<-c("cd","uc")
j=2

for (chr in c(6,10,11)) {
  
  print(chr)
  
  tmp<-fread(paste(path,"post_imputation/analysis/metaanalysis/",pheno[j],"/",chr,"_meta_noGC_PCs_firthse_1_formatted_A2_effect_and_alternative.txt",sep=""))
  tmp<-tmp[which(tmp$MarkerName %in% c("chr6:32279268:T:G","chr10:99524480:T:G","chr11:114453250:A:G")),]
  
  for (jj in 1:length(array)) {
    
    file_tmp11<-paste(path,"post_imputation/analysis/regenie/",array[jj],"/",pheno[j],"/chr",chr,"_",array[jj],"_step2_",pheno[j],"_eur_sex_PCs_firthse.regenie",sep="")
    
    if(file.exists(file_tmp11)) {
      
      tmp11<-fread(file_tmp11,head=T)
      tmp11<-tmp11[which(tmp11$ID %in%c("chr6:32279268:T:G","chr10:99524480:T:G","chr11:114453250:A:G")),]
      tmp11<-tmp11[,c("ID","A1FREQ","INFO","BETA.Y1","LOG10P.Y1")]
      colnames(tmp11)[2:ncol(tmp11)]<-paste(colnames(tmp11)[2:ncol(tmp11)],array[jj],sep="_")
      
      if(jj==1) {
        tmp1<-tmp11
      } else {
        tmp1<-merge(tmp1,tmp11,by="ID",all=T)
      }
      
      rm(tmp11)
    }
  }
  
  tmp<-merge(tmp,tmp1,by.x="MarkerName",by.y="ID")
  
  if(chr==6) {
    tmp_final<-tmp
  }else {
    tmp_final<-rbind(tmp_final,tmp)
  }
}

> tmp_final
# MarkerName A1 A2    BETA     SE    P-value Direction_ed HetISq
# 1:   chr6:32279268:T:G  T  G  0.9893 0.0406 2.720e-131    +++++++++   78.7
# 2:  chr10:99524480:T:G  T  G -0.1499 0.0135  1.798e-28    ---------   79.3
# 3: chr11:114453250:A:G  A  G -0.0957 0.0140  8.580e-12    -+------+   85.2
# HetChiSq HetDf   HetPVal A1FREQ_illumina370 INFO_illumina370
# 1:   37.615     8 8.868e-06          0.0189447         0.949577
# 2:   38.638     8 5.738e-06          0.4756480         0.952060
# 3:   54.165     8 6.412e-09          0.3835060         0.976007
# BETA.Y1_illumina370 LOG10P.Y1_illumina370 A1FREQ_illumina550
# 1:          0.67579200             1.8192000          0.0211946
# 2:         -0.00978619             0.0465522          0.4887210
# 3:         -0.09328910             0.6360630          0.3741310
# INFO_illumina550 BETA.Y1_illumina550 LOG10P.Y1_illumina550
# 1:         0.975905            0.809968               4.02141
# 2:         0.932998           -0.496402              11.98120
# 3:         0.970576            0.369972               7.28258
# A1FREQ_affymetrix6 INFO_affymetrix6 BETA.Y1_affymetrix6
# 1:          0.0326951         1.003780            0.768708
# 2:          0.4954950         0.980827           -0.153737
# 3:          0.3613120         0.991913           -0.123028
# LOG10P.Y1_affymetrix6 A1FREQ_humanomniexpress INFO_humanomniexpress
# 1:              14.61230               0.0536631              0.969786
# 2:               4.96552               0.4902830              1.013390
# 3:               3.15522               0.3676490              0.965495
# BETA.Y1_humanomniexpress LOG10P.Y1_humanomniexpress A1FREQ_humancoreexome
# 1:                1.9407400                  16.356900              0.037829
# 2:               -0.2603440                   2.638640              0.500335
# 3:                0.0385604                   0.172993              0.355447
# INFO_humancoreexome BETA.Y1_humancoreexome LOG10P.Y1_humancoreexome
# 1:            0.978284              0.8595040                 36.62520
# 2:            0.994067             -0.1384480                  6.40234
# 3:            0.987128             -0.0968758                  3.13879
# A1FREQ_humanomni1 INFO_humanomni1 BETA.Y1_humanomni1 LOG10P.Y1_humanomni1
# 1:         0.0191345        0.938783          1.2003100            4.1517500
# 2:         0.5178030        1.018670         -0.0119055            0.0459513
# 3:         0.3671590        0.988869         -0.1508140            0.8978700
# A1FREQ_quad610 INFO_quad610 BETA.Y1_quad610 LOG10P.Y1_quad610 A1FREQ_gsa
# 1:      0.0301917     0.956209        1.320320          11.06060  0.0201216
# 2:      0.4542160     0.942088       -0.154784           1.59401  0.4943140
# 3:      0.3762140     1.013690       -0.116665           1.04155  0.3716490
# INFO_gsa BETA.Y1_gsa LOG10P.Y1_gsa A1FREQ_illuminaexome INFO_illuminaexome
# 1: 0.962635    1.167680      55.63330            0.0323471           0.959036
# 2: 0.992683   -0.129345       9.71900            0.4827730           0.906165
# 3: 1.001730   -0.130352       9.22708            0.3547850           1.061950
# BETA.Y1_illuminaexome LOG10P.Y1_illuminaexome
# 1:              0.656918                0.785173
# 2:             -0.552094                2.905340
# 3:             -0.282172                1.105370


########################################
# run meta-analysis excluding illumina 550

array<-c("illumina370","illumina550","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")

pheno<-c("cd","uc")
j=2

for (chr in c(6,10,11)) {
  
  print(chr)
  
  tmp<-fread(paste(path,"post_imputation/analysis/metaanalysis/",pheno[j],"/",chr,"_meta_noGC_PCs_firthse_no_ilu550_1.txt",sep=""))
  tmp<-tmp[which(tmp$MarkerName %in% c("chr6:32279268:T:G","chr10:99524480:T:G","chr11:114453250:A:G")),]
  
  for (jj in 1:length(array)) {
    
    file_tmp11<-paste(path,"post_imputation/analysis/regenie/",array[jj],"/",pheno[j],"/chr",chr,"_",array[jj],"_step2_",pheno[j],"_eur_sex_PCs_firthse.regenie",sep="")
    
    if(file.exists(file_tmp11)) {
      
      tmp11<-fread(file_tmp11,head=T)
      tmp11<-tmp11[which(tmp11$ID %in%c("chr6:32279268:T:G","chr10:99524480:T:G","chr11:114453250:A:G")),]
      tmp11<-tmp11[,c("ID","A1FREQ","INFO","BETA.Y1","LOG10P.Y1")]
      colnames(tmp11)[2:ncol(tmp11)]<-paste(colnames(tmp11)[2:ncol(tmp11)],array[jj],sep="_")
      
      if(jj==1) {
        tmp1<-tmp11
      } else {
        tmp1<-merge(tmp1,tmp11,by="ID",all=T)
      }
      
      rm(tmp11)
    }
  }
  
  tmp<-merge(tmp,tmp1,by.x="MarkerName",by.y="ID")
  
  if(chr==6) {
    tmp_final<-tmp
  }else {
    tmp_final<-rbind(tmp_final,tmp)
  }
}

tmp_final
# MarkerName Allele1 Allele2  Effect StdErr    P-value Direction
# 1:   chr6:32279268:T:G       t       g  0.9964 0.0414 3.780e-128  ++++++++
#   2:  chr10:99524480:T:G       t       g -0.1363 0.0138  5.544e-23  --------
#   3: chr11:114453250:A:G       a       g -0.1164 0.0143  4.421e-16  -------+
#   HetISq HetChiSq HetDf   HetPVal A1FREQ_illumina370 INFO_illumina370
# 1:   81.0   36.839     7 5.032e-06          0.0189447         0.949577
# 2:   45.9   12.931     7 7.380e-02          0.4756480         0.952060
# 3:    0.0    5.120     7 6.453e-01          0.3835060         0.976007
# BETA.Y1_illumina370 LOG10P.Y1_illumina370 A1FREQ_illumina550
# 1:          0.67579200             1.8192000          0.0211946
# 2:         -0.00978619             0.0465522          0.4887210
# 3:         -0.09328910             0.6360630          0.3741310
# INFO_illumina550 BETA.Y1_illumina550 LOG10P.Y1_illumina550
# 1:         0.975905            0.809968               4.02141
# 2:         0.932998           -0.496402              11.98120
# 3:         0.970576            0.369972               7.28258
# A1FREQ_affymetrix6 INFO_affymetrix6 BETA.Y1_affymetrix6
# 1:          0.0326951         1.003780            0.768708
# 2:          0.4954950         0.980827           -0.153737
# 3:          0.3613120         0.991913           -0.123028
# LOG10P.Y1_affymetrix6 A1FREQ_humanomniexpress INFO_humanomniexpress
# 1:              14.61230               0.0536631              0.969786
# 2:               4.96552               0.4902830              1.013390
# 3:               3.15522               0.3676490              0.965495
# BETA.Y1_humanomniexpress LOG10P.Y1_humanomniexpress A1FREQ_humancoreexome
# 1:                1.9407400                  16.356900              0.037829
# 2:               -0.2603440                   2.638640              0.500335
# 3:                0.0385604                   0.172993              0.355447
# INFO_humancoreexome BETA.Y1_humancoreexome LOG10P.Y1_humancoreexome
# 1:            0.978284              0.8595040                 36.62520
# 2:            0.994067             -0.1384480                  6.40234
# 3:            0.987128             -0.0968758                  3.13879
# A1FREQ_humanomni1 INFO_humanomni1 BETA.Y1_humanomni1 LOG10P.Y1_humanomni1
# 1:         0.0191345        0.938783          1.2003100            4.1517500
# 2:         0.5178030        1.018670         -0.0119055            0.0459513
# 3:         0.3671590        0.988869         -0.1508140            0.8978700
# A1FREQ_quad610 INFO_quad610 BETA.Y1_quad610 LOG10P.Y1_quad610 A1FREQ_gsa
# 1:      0.0301917     0.956209        1.320320          11.06060  0.0201216
# 2:      0.4542160     0.942088       -0.154784           1.59401  0.4943140
# 3:      0.3762140     1.013690       -0.116665           1.04155  0.3716490
# INFO_gsa BETA.Y1_gsa LOG10P.Y1_gsa A1FREQ_illuminaexome INFO_illuminaexome
# 1: 0.962635    1.167680      55.63330            0.0323471           0.959036
# 2: 0.992683   -0.129345       9.71900            0.4827730           0.906165
# 3: 1.001730   -0.130352       9.22708            0.3547850           1.061950
# BETA.Y1_illuminaexome LOG10P.Y1_illuminaexome
# 1:              0.656918                0.785173
# 2:             -0.552094                2.905340
# 3:             -0.282172                1.105370


###############

cd$Heterogeneticy_effec<-cut(cd$HetPVal,breaks = c(0,1E-7,1E-6,1E-5,1E-4,1E-3,1E-0))

p1<-ggscatter(data=cd, 'BETA.Y1_gsa',"BETA.Y1_hce",
              xlim=c(min(cd$BETA.Y1_gsa,cd$BETA.Y1_hce),max(cd$BETA.Y1_gsa,cd$BETA.Y1_hce)),ylim=c(min(cd$BETA.Y1_gsa,cd$BETA.Y1_hce),max(cd$BETA.Y1_gsa,cd$BETA.Y1_hce)),
              # add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "Regenie GSA CD Beta", merge=T,
              ylab = "Regenie GWAS3 CD Beta",color="Heterogeneticy_effec") + 
  scale_color_manual(values=rev(c("#4575b4","#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026")))

p2<-ggscatter(data=cd, 'LOG10P.Y1_gsa', "LOG10P.Y1_hce",
              xlim=c(0,max(cd$LOG10P.Y1_gsa,cd$LOG10P.Y1_hce)),ylim=c(0,max(cd$LOG10P.Y1_gsa,cd$LOG10P.Y1_hce)),
              # add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "Regenie GSA CD log10(P)", ylab = "Regenie GWAS3 CD log10(P)",color="Heterogeneticy_effec") + 
  scale_color_manual(values=rev(c("#4575b4","#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026")))
p<-ggarrange(p1,p2,ncol=2,common.legend = T)
dev.off()

pdf(paste("~/tmp_plots/allchr_CD_effect_size_pvalue_gsa_vs_hce.pdf",sep=""),width = 14)
print(annotate_figure(p,top = text_grob("CD p-value < 0.001 in GSA and GWAS3")))
dev.off()

###

p1<-ggscatter(data=cd[which(cd$HetPVal<0.001),], 'BETA.Y1_gsa',"BETA.Y1_hce",
              xlim=c(min(cd$BETA.Y1_gsa,cd$BETA.Y1_hce),max(cd$BETA.Y1_gsa,cd$BETA.Y1_hce)),ylim=c(min(cd$BETA.Y1_gsa,cd$BETA.Y1_hce),max(cd$BETA.Y1_gsa,cd$BETA.Y1_hce)),
              # add = "reg.line", conf.int = TRUE, 
              # cor.coef = TRUE, cor.method = "pearson",
              xlab = "Regenie GSA CD Beta", merge=T,
              ylab = "Regenie GWAS3 CD Beta",color="Heterogeneticy_effec") + 
  scale_color_manual(values=rev(c("#4575b4","#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026")))

p2<-ggscatter(data=cd[which(cd$HetPVal<0.001),], 'LOG10P.Y1_gsa', "LOG10P.Y1_hce",
              xlim=c(0,max(cd$LOG10P.Y1_gsa,cd$LOG10P.Y1_hce)),ylim=c(0,max(cd$LOG10P.Y1_gsa,cd$LOG10P.Y1_hce)),
              # add = "reg.line", conf.int = TRUE, 
              # cor.coef = TRUE, cor.method = "pearson",
              xlab = "Regenie GSA CD log10(P)", ylab = "Regenie GWAS3 CD log10(P)",color="Heterogeneticy_effec") + 
  scale_color_manual(values=rev(c("#4575b4","#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026")))
p<-ggarrange(p1,p2,ncol=2,common.legend = T)
dev.off()

p<-ggarrange(p1,p2,ncol=2,common.legend = T)
pdf(paste("~/tmp_plots/allchr_CD_effect_size_pvalue_gsa_vs_hce_hetPval_filter.pdf",sep=""),width = 14)
print(annotate_figure(p,top = text_grob("CD p-value < 0.001 in GSA and GWAS3")))
dev.off()


dim(cd[which(cd$HetPVal<0.001),])
# [1] 240  25

dim(cd[which(cd$HetPVal<0.0001),])
# [1] 168  25

dim(cd[which(cd$HetPVal<1E-5),])
# [1] 93 25

cd[which( (cd$BETA.Y1_gsa>0 & cd$BETA.Y1_hce<0) | (cd$BETA.Y1_gsa<0 & cd$BETA.Y1_hce>0) ),]


cd[which( cd$HetPVal<1E-5),c("MarkerName")]
# variants across 4 chr, get one example per region:

cd[which(cd$MarkerName %in% c("chr3:49588749:T:G","chr5:96897759:C:G","chr10:99514301:T:C","chr22:39312328:T:C")),c("MarkerName","Direction_ed","BETA","P.value","HetPVal",
                                                                                                                    "BETA.Y1_gsa","LOG10P.Y1_gsa","A1FREQ_gsa","INFO_gsa",
                                                                                                                    "BETA.Y1_hce","LOG10P.Y1_hce","A1FREQ_hce","INFO_hce")]
# MarkerName Direction_ed    BETA   P.value   HetPVal BETA.Y1_gsa
# 2660   chr3:49588749:T:G     -------- -0.1267 6.197e-23 1.175e-12  -0.0686843
# 5168   chr5:96897759:C:G     +------- -0.1210 2.908e-19 1.485e-24  -0.0806001
# 10066 chr10:99514301:T:C     +------- -0.1440 8.763e-29 2.917e-08  -0.1273200
# 12681 chr22:39312328:T:C     ---+---- -0.1099 7.862e-17 3.831e-08  -0.0866397
# LOG10P.Y1_gsa A1FREQ_gsa INFO_gsa BETA.Y1_hce LOG10P.Y1_hce A1FREQ_hce
# 2660        3.89966   0.503528 1.005990   -0.136785       6.99464   0.520367
# 5168        4.82479   0.364858 0.980572   -0.112663       4.41225   0.347316
# 10066      11.73970   0.493363 0.975136   -0.119209       5.36467   0.500973
# 12681       5.60291   0.593435 0.992789   -0.120049       5.13840   0.583044
# INFO_hce
# 2660  0.989985
# 5168  0.960733
# 10066 0.965241
# 12681 0.934130


array<-c("illumina370","illumina550","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")

pheno<-c("cd","uc")
j=1

for (chr in c(3,5,10,22)) {
  
  print(chr)
  
  tmp<-fread(paste(path,"post_imputation/analysis/metaanalysis/",pheno[j],"/",chr,"_meta_noGC_PCs_firthse_1_formatted_A2_effect_and_alternative.txt",sep=""))
  tmp<-tmp[which(tmp$MarkerName %in% c("chr3:49588749:T:G","chr5:96897759:C:G","chr10:99514301:T:C","chr22:39312328:T:C")),]
  
  for (jj in 1:length(array)) {
    
    file_tmp11<-paste(path,"post_imputation/analysis/regenie/",array[jj],"/",pheno[j],"/chr",chr,"_",array[jj],"_step2_",pheno[j],"_eur_sex_PCs_firthse.regenie",sep="")
    
    if(file.exists(file_tmp11)) {
      
      tmp11<-fread(file_tmp11,head=T)
      tmp11<-tmp11[which(tmp11$ID %in% c("chr3:49588749:T:G","chr5:96897759:C:G","chr10:99514301:T:C","chr22:39312328:T:C")),]
      tmp11<-tmp11[,c("ID","A1FREQ","INFO","BETA.Y1","LOG10P.Y1")]
      colnames(tmp11)[2:ncol(tmp11)]<-paste(colnames(tmp11)[2:ncol(tmp11)],array[jj],sep="_")
      
      if(jj==1) {
        tmp1<-tmp11
      } else {
        tmp1<-merge(tmp1,tmp11,by="ID",all=T)
      }
      
      rm(tmp11)
    }
  }
  
  tmp<-merge(tmp,tmp1,by.x="MarkerName",by.y="ID")
  
  if(chr==3) {
    tmp_final<-tmp
  }else {
    tmp_final<-rbind(tmp_final,tmp)
  }
}



# > tmp_final
# MarkerName A1 A2    BETA     SE   P-value Direction_ed HetISq
# 1:  chr3:49588749:T:G  T  G -0.1267 0.0128 6.197e-23     --------   90.1
# 2:  chr5:96897759:C:G  C  G -0.1210 0.0135 2.908e-19     +-------   94.5
# 3: chr10:99514301:T:C  T  C -0.1440 0.0129 8.763e-29     +-------   85.5
# 4: chr22:39312328:T:C  T  C -0.1099 0.0132 7.862e-17     ---+----   85.4
# HetChiSq HetDf   HetPVal A1FREQ_illumina370 INFO_illumina370
# 1:   70.492     7 1.175e-12           0.486282         1.020980
# 2:  128.212     7 1.485e-24           0.389467         0.937127
# 3:   48.443     7 2.917e-08           0.476490         0.944886
# 4:   47.837     7 3.831e-08           0.619102         1.005340
# BETA.Y1_illumina370 LOG10P.Y1_illumina370 A1FREQ_illumina550
# 1:          -0.1368630              1.183260           0.503104
# 2:           0.0873753              0.575021           0.345108
# 3:           0.0365760              0.199943           0.491817
# 4:          -0.1122600              0.827504           0.596665
# INFO_illumina550 BETA.Y1_illumina550 LOG10P.Y1_illumina550
# 1:         0.883229           -0.535864               23.1107
# 2:         0.930349           -0.716976               36.6416
# 3:         0.919529           -0.444532               16.9820
# 4:         0.968137           -0.427005               16.3860
# A1FREQ_affymetrix500 INFO_affymetrix500 BETA.Y1_affymetrix500
# 1:             0.511273           1.002110            -0.1592290
# 2:             0.337494           0.958660            -0.0750376
# 3:             0.508788           0.986499            -0.2133530
# 4:             0.581717           0.977043            -0.0330527
# LOG10P.Y1_affymetrix500 A1FREQ_humancoreexome INFO_humancoreexome
# 1:                3.223690              0.520367            0.989985
# 2:                0.882928              0.347316            0.960733
# 3:                5.318530              0.500973            0.965241
# 4:                0.312238              0.583044            0.934130
# BETA.Y1_humancoreexome LOG10P.Y1_humancoreexome A1FREQ_humanomni1
# 1:              -0.136785                  6.99464          0.512021
# 2:              -0.112663                  4.41225          0.360700
# 3:              -0.119209                  5.36467          0.521547
# 4:              -0.120049                  5.13840          0.602090
# INFO_humanomni1 BETA.Y1_humanomni1 LOG10P.Y1_humanomni1 A1FREQ_quad610
# 1:        1.017600        -0.12504700            0.7815840       0.495334
# 2:        0.999423        -0.02954170            0.1225650       0.359949
# 3:        1.026530        -0.05435870            0.2627680       0.461556
# 4:        0.990131         0.00934319            0.0359121       0.605079
# INFO_quad610 BETA.Y1_quad610 LOG10P.Y1_quad610 A1FREQ_gsa INFO_gsa
# 1:     0.988379      -0.1489230         1.8623200   0.503528 1.005990
# 2:     0.966239      -0.0987577         0.9250580   0.364858 0.980572
# 3:     0.951528      -0.0592817         0.4736470   0.493363 0.975136
# 4:     0.955922      -0.0118347         0.0704414   0.593435 0.992789
# BETA.Y1_gsa LOG10P.Y1_gsa A1FREQ_illuminaexome INFO_illuminaexome
# 1:  -0.0686843       3.89966             0.479031           0.933071
# 2:  -0.0806001       4.82479             0.363639           0.866038
# 3:  -0.1273200      11.73970             0.501118           0.910533
# 4:  -0.0866397       5.60291             0.606689           0.943233
# BETA.Y1_illuminaexome LOG10P.Y1_illuminaexome
# 1:            -0.0671534                0.157833
# 2:            -0.3659780                1.291500
# 3:            -0.4094520                1.773630
# 4:            -0.2213200                0.702555



########################################
# run meta-analysis excluding illumina 550

array<-c("illumina370","illumina550","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")

pheno<-c("cd","uc")
j=1

for (chr in c(3,5,10,22)) {
  
  print(chr)
  
  tmp<-fread(paste(path,"post_imputation/analysis/metaanalysis/",pheno[j],"/",chr,"_meta_noGC_PCs_firthse_no_ilu550_1.txt",sep=""))
  tmp<-tmp[which(tmp$MarkerName %in% c("chr3:49588749:T:G","chr5:96897759:C:G","chr10:99514301:T:C","chr22:39312328:T:C")),]
  
  for (jj in 1:length(array)) {
    
    file_tmp11<-paste(path,"post_imputation/analysis/regenie/",array[jj],"/",pheno[j],"/chr",chr,"_",array[jj],"_step2_",pheno[j],"_eur_sex_PCs_firthse.regenie",sep="")
    
    if(file.exists(file_tmp11)) {
      
      tmp11<-fread(file_tmp11,head=T)
      tmp11<-tmp11[which(tmp11$ID %in% c("chr3:49588749:T:G","chr5:96897759:C:G","chr10:99514301:T:C","chr22:39312328:T:C")),]
      tmp11<-tmp11[,c("ID","A1FREQ","INFO","BETA.Y1","LOG10P.Y1")]
      colnames(tmp11)[2:ncol(tmp11)]<-paste(colnames(tmp11)[2:ncol(tmp11)],array[jj],sep="_")
      
      if(jj==1) {
        tmp1<-tmp11
      } else {
        tmp1<-merge(tmp1,tmp11,by="ID",all=T)
      }
      
      rm(tmp11)
    }
  }
  
  tmp<-merge(tmp,tmp1,by.x="MarkerName",by.y="ID")
  
  if(chr==3) {
    tmp_final<-tmp
  }else {
    tmp_final<-rbind(tmp_final,tmp)
  }
}


# MarkerName Allele1 Allele2  Effect StdErr   P-value Direction HetISq
# 1:  chr3:49588749:T:G       t       g -0.1014 0.0132 1.899e-14   -------   22.5
# 2:  chr5:96897759:C:G       c       g -0.0845 0.0139 1.204e-09   +------   30.0
# 3: chr10:99514301:T:C       t       c -0.1241 0.0134 1.567e-20   +------   52.6
# 4: chr22:39312328:T:C       t       c -0.0870 0.0137 1.927e-10   --+----    0.3
# HetChiSq HetDf HetPVal A1FREQ_illumina370 INFO_illumina370
# 1:    7.742     6 0.25770           0.486282         1.020980
# 2:    8.565     6 0.19950           0.389467         0.937127
# 3:   12.660     6 0.04877           0.476490         0.944886
# 4:    6.019     6 0.42110           0.619102         1.005340
# BETA.Y1_illumina370 LOG10P.Y1_illumina370 A1FREQ_illumina550
# 1:          -0.1368630              1.183260           0.503104
# 2:           0.0873753              0.575021           0.345108
# 3:           0.0365760              0.199943           0.491817
# 4:          -0.1122600              0.827504           0.596665
# INFO_illumina550 BETA.Y1_illumina550 LOG10P.Y1_illumina550
# 1:         0.883229           -0.535864               23.1107
# 2:         0.930349           -0.716976               36.6416
# 3:         0.919529           -0.444532               16.9820
# 4:         0.968137           -0.427005               16.3860
# A1FREQ_affymetrix500 INFO_affymetrix500 BETA.Y1_affymetrix500
# 1:             0.511273           1.002110            -0.1592290
# 2:             0.337494           0.958660            -0.0750376
# 3:             0.508788           0.986499            -0.2133530
# 4:             0.581717           0.977043            -0.0330527
# LOG10P.Y1_affymetrix500 A1FREQ_humancoreexome INFO_humancoreexome
# 1:                3.223690              0.520367            0.989985
# 2:                0.882928              0.347316            0.960733
# 3:                5.318530              0.500973            0.965241
# 4:                0.312238              0.583044            0.934130
# BETA.Y1_humancoreexome LOG10P.Y1_humancoreexome A1FREQ_humanomni1
# 1:              -0.136785                  6.99464          0.512021
# 2:              -0.112663                  4.41225          0.360700
# 3:              -0.119209                  5.36467          0.521547
# 4:              -0.120049                  5.13840          0.602090
# INFO_humanomni1 BETA.Y1_humanomni1 LOG10P.Y1_humanomni1 A1FREQ_quad610
# 1:        1.017600        -0.12504700            0.7815840       0.495334
# 2:        0.999423        -0.02954170            0.1225650       0.359949
# 3:        1.026530        -0.05435870            0.2627680       0.461556
# 4:        0.990131         0.00934319            0.0359121       0.605079
# INFO_quad610 BETA.Y1_quad610 LOG10P.Y1_quad610 A1FREQ_gsa INFO_gsa
# 1:     0.988379      -0.1489230         1.8623200   0.503528 1.005990
# 2:     0.966239      -0.0987577         0.9250580   0.364858 0.980572
# 3:     0.951528      -0.0592817         0.4736470   0.493363 0.975136
# 4:     0.955922      -0.0118347         0.0704414   0.593435 0.992789
# BETA.Y1_gsa LOG10P.Y1_gsa A1FREQ_illuminaexome INFO_illuminaexome
# 1:  -0.0686843       3.89966             0.479031           0.933071
# 2:  -0.0806001       4.82479             0.363639           0.866038
# 3:  -0.1273200      11.73970             0.501118           0.910533
# 4:  -0.0866397       5.60291             0.606689           0.943233
# BETA.Y1_illuminaexome LOG10P.Y1_illuminaexome
# 1:            -0.0671534                0.157833
# 2:            -0.3659780                1.291500
# 3:            -0.4094520                1.773630
# 4:            -0.2213200                0.702555


############################################################################################################################################################
############################################################################################################################################################


############# compare pvalue and het pvalue with and without studies - RUN NOW AS INDEP SCRIPT 'sensitivity_test_meta_analysis_results.R'

pheno=(ibd cd uc)
MEM="4000"
for ph in ${pheno[@]}
do
bsub -J"senst" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/analysis/metaanalysis/log/stderr_sensitivity_test_${ph} \
-o ${path_gwas}post_imputation/analysis/metaanalysis/log/stdout_def_sensitivity_test_${ph} \
"/software/R-4.3.1/bin/Rscript ${path_gwas}scripts/sensitivity_test_meta_analysis_results.R ${ph} > \
${path_gwas}scripts/sensitivity_test_meta_analysis_results_${ph}.Rout"
done
# Job <600788> is submitted to default queue <normal>.
# Job <600789> is submitted to default queue <normal>.
# Job <600790> is submitted to default queue <normal>.

for ph in ${pheno[@]}
do
tail -200  ${path_gwas}post_imputation/analysis/metaanalysis/log/stdout_def_sensitivity_test_${ph} | grep "Successfully completed."
done

pheno=(ibd cd uc)
MEM="4000"
for ph in ${pheno[@]}
do
bsub -J"senst" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/analysis/metaanalysis/log/stderr_sensitivity_test_${ph}_no_illumina550 \
-o ${path_gwas}post_imputation/analysis/metaanalysis/log/stdout_def_sensitivity_test_${ph}_no_illumina550 \
"/software/R-4.3.1/bin/Rscript ${path_gwas}scripts/sensitivity_test_meta_analysis_results_no_illumina550.R ${ph} > \
${path_gwas}scripts/sensitivity_test_meta_analysis_results_no_illumina550_${ph}.Rout"
done
# Job <582185> is submitted to default queue <normal>.
# Job <582186> is submitted to default queue <normal>.
# Job <582187> is submitted to default queue <normal>.




# MEM=4000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -m "modern_hardware" -G ibdgwas /software/bin/R-4.3.1
# 
# library(data.table)
# library(ggpubr)
# library(ggExtra)
# library(gridExtra)
# library(grid)
# 
# path<-"/path/to/ibdgwas/IIBDGC/"
# 
# array<-c("illumina370","illumina550","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")
# # pheno<-c("cd","uc","ibd")
# 
# pheno<-c("cd")
# 
# for (j in 1:length(pheno)) {
# 
#   for (chr in 1:22) {
#     
#     print(chr)
#     
#     tmp1<-fread(paste(path,"post_imputation/analysis/metaanalysis/",pheno[j],"/",chr,"_meta_noGC_PCs_firthse_1.txt",sep=""))
#     
#     for (i in 1:length(array)) {
#       
#       file_tmp2<-paste(path,"post_imputation/analysis/metaanalysis/",pheno[j],"/",chr,"_meta_noGC_PCs_firthse_no_",array[i],"_1.txt",sep="")
#       
#       if(file.exists(file_tmp2)) {
#         tmp2<-fread(file_tmp2,head=T)
#         tmp2<-tmp2[,c("MarkerName","Effect","P-value","HetPVal")]
#         colnames(tmp2)<-c("MarkerName","Effect","PVal","HetPVal")
#         colnames(tmp2)[2:ncol(tmp2)]<-paste(colnames(tmp2)[2:ncol(tmp2)],"_no_",array[i],sep="")
#         
#         if(i==1) {
#           tmp2_f<-tmp2
#           rm(tmp2,file_tmp2)
#         } else {
#           tmp2_f<-merge(tmp2_f,tmp2,by="MarkerName")
#           rm(tmp2,file_tmp2)
#         }
#       }
#       
#     }
#     
#     tmp<-merge(tmp1,tmp2_f,by="MarkerName")
#     rm(tmp1,tmp2)
#     
#     if(chr==1) {
#       tmp_final<-tmp
#     } else {
#       tmp_final<-rbind(tmp_final,tmp)
#     }
#     rm(tmp,tmp1)
#   }
# 
#   # assign(pheno[j],tmp_final)
#   # rm(tmp_final)
#   
#   # keep only variants with data from all studies
#   
#   tmp_final<-tmp_final[which(tmp_final$HetDf==max(tmp_final$HetDf)),]
#   colnames(tmp_final)[6]<-"PVal"
#   
#   # trim the data to keep just variants with effect PVal<0.1
# 
#   tmp_final<-tmp_final[which(tmp_final$PVal<0.01),]
# 
#   tmp_final<-as.data.frame(tmp_final)
#  
# 
#   if (pheno[j]=="cd") {
#     array_pheno<-c("illumina370","illumina550","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")
#   }
#   
#   if (pheno[j]=="uc") {
#     array_pheno<-c("illumina370","illumina550","affymetrix6","humanomniexpress","humancoreexome","humanomni1","quad610","gsa","illuminaexome")
#   }
#   
#   if (pheno[j]=="ibd") {
#     array_pheno<-c("illumina370","illumina550","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")
#   }
#   
#   ######################
#   ####### plot HetPval
#   
#   for (i in 1:length(array_pheno)) {
#     ids<-c("HetPVal",paste("HetPVal_no_",array_pheno[i],sep=""))
#     tmp<-tmp_final[,ids]
#     tmp$logHetPVal<- -log10(tmp$HetPVal)
#     tmp$logHetPVal_no_study<- -log10(tmp[,2])
#     tmp<-tmp[,3:4]
#     
#     uplim<-max(tmp)
#     
#     p1<-ggplot(tmp, aes(logHetPVal,logHetPVal_no_study)) +
#       geom_point() + ylim(0,uplim) + xlim(0,uplim) + geom_abline(intercept = 0, slope= 1, linetype="dashed", color = "red") +
#       ggtitle(array_pheno[i]) + xlab(" -log10(Heterogeneity Pval)\nall studies") + ylab(paste(" -log10(Heterogeneity Pval)\nall studies but ",array_pheno[i],sep=""))
#     
#     assign(paste("plot_HetPVal_",pheno[j],"_",array_pheno[i],sep=""),p1)
#   }
#   
#   l = mget(ls()[grep("plot_HetPVal_",ls())])
#   ggsave(paste("~/tmp_plots/",pheno[j],"_HetPval_sensitivity_test.pdf",sep=""), arrangeGrob(grobs = l), width = 14, height = 14)
#   
#   #####################
#   # PLOT BETA
#   
#   for (i in 1:length(array_pheno)) {
#     ids<-c("Effect",paste("Effect_no_",array_pheno[i],sep=""))
#     tmp<-tmp_final[,ids]
#     colnames(tmp)[2]<-"Effect_no_study"
#     
#     uplim<-max(tmp)
#     lowlim<-min(tmp)
#     
#     p1<-ggplot(tmp, aes(Effect,Effect_no_study)) +
#       geom_point() + ylim(lowlim,uplim) + xlim(lowlim,uplim) + geom_abline(intercept = 0, slope= 1, linetype="dashed", color = "red") +
#       ggtitle(array_pheno[i]) + xlab("Beta\nall studies") + ylab(paste("Beta\nall studies but ",array_pheno[i],sep=""))
#     
#     assign(paste("plot_beta_",pheno[j],"_",array_pheno[i],sep=""),p1)
#   }
#   
#   l = mget(ls()[grep("plot_beta_",ls())])
#   ggsave(paste("~/tmp_plots/",pheno[j],"_Beta_sensitivity_test.pdf",sep=""), arrangeGrob(grobs = l), width = 14, height = 14)
#   
#   #####################
#   # PLOT as well PVal
# 
#   for (i in 1:length(array_pheno)) {
#     ids<-c("PVal",paste("PVal_no_",array_pheno[i],sep=""))
#     tmp<-tmp_final[,ids]
#     tmp$logPVal<- -log10(tmp$PVal)
#     tmp$logPVal_no_study<- -log10(tmp[,2])
#     tmp<-tmp[,3:4]
#     
#     uplim<-max(tmp)
#     
#     p1<-ggplot(tmp, aes(logPVal,logPVal_no_study)) +
#       geom_point() + ylim(0,uplim) + xlim(0,uplim) + geom_abline(intercept = 0, slope= 1, linetype="dashed", color = "red") +
#       ggtitle(array_pheno[i]) + xlab(" -log10(PVal)\nall studies") + ylab(paste(" -log10(PVal)\nall studies but ",array_pheno[i],sep=""))
# 
#     assign(paste("plot_PVal_",pheno[j],"_",array_pheno[i],sep=""),p1)
#   }
#   
#   l = mget(ls()[grep("plot_PVal_",ls())])
#   ggsave(paste("~/tmp_plots/",pheno[j],"_PVal_sensitivity_test.pdf",sep=""),arrangeGrob(grobs = l), width = 14, height = 14)
#   
# }
# 
# q("no")

# MEM="4000"
# bsub -J"ph" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 16 \
# -e ${path_gwas}post_imputation/analysis/metaanalysis/log/stderr_def_associated_ph_known \
# -o ${path_gwas}post_imputation/analysis/metaanalysis/log/stdout_def_associated_ph_known \
# "/software/R-4.3.1/bin/R CMD BATCH ${path_gwas}scripts/define_associated_phenotype_known_variants_mc16.R"
# # Job <189389> is submitted to default queue <normal>


# # keep only variants with data from all studies
# cd<-cd[which(cd$HetDf==max(cd$HetDf)),]
# colnames(cd)[6]<-"PVal"
# uc<-uc[which(uc$HetDf==max(uc$HetDf)),]
# colnames(uc)[6]<-"PVal"
# ibd<-ibd[which(ibd$HetDf==max(ibd$HetDf)),]
# colnames(ibd)[6]<-"PVal"

# dim(cd)
# # [1] 9942255      27
# dim(uc)
# # [1] 9878930      29
# dim(ibd)
# # [1] 9878930      29
# 
# dim(cd[which(cd$PVal<0.01),])
# # [1] 347132     27
# dim(uc[which(uc$PVal<0.01),])
# # [1] 300886     29
# dim(ibd[which(ibd$PVal<0.01),])
# # [1] 372733     31

# # trim the data to keep just variants with effect PVal<0.1
# 
# cd<-cd[which(cd$PVal<0.01),]
# uc<-uc[which(uc$PVal<0.01),]
# ibd<-ibd[which(ibd$PVal<0.01),]
# 
# cd<-as.data.frame(cd)
# uc<-as.data.frame(uc)
# ibd<-as.data.frame(ibd)


######################
####### plot HetPval

########### CD

array_cd<-c("illumina370","illumina550","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")

for (i in 1:length(array_cd)) {
  ids<-c("HetPVal",paste("HetPVal_no_",array_cd[i],sep=""))
  tmp<-cd[,ids]
  tmp$logHetPVal<- -log10(tmp$HetPVal)
  tmp$logHetPVal_no_study<- -log10(tmp[,2])
  tmp<-tmp[,3:4]
  
  max(tmp)
  uplim<-max(tmp)
  
  p1<-ggplot(tmp, aes(logHetPVal,logHetPVal_no_study)) +
    geom_point() + ylim(0,uplim) + xlim(0,uplim) + geom_abline(intercept = 0, slope= 1, linetype="dashed", color = "red") +
    ggtitle(array_cd[i]) + xlab(" -log10(Heterogeneity Pval)\nall studies") + ylab(paste(" -log10(Heterogeneity Pval)\nall studies but ",array_cd[i],sep=""))
  # p1<-ggMarginal(tmp)
  
  assign(paste("plot_cd_",array_cd[i],sep=""),p1)
}

l = mget(ls()[grep("plot_cd",ls())])
# ggsave("~/tmp_plots/CD_HetPval_sensitibity_test.pdf", arrangeGrob(grobs = l), width = 14, height = 14)
ggsave("~/tmp_plots/CD_HetPval_sensitibity_test_new.pdf", arrangeGrob(grobs = l), width = 14, height = 14)


########### UC

array_uc<-c("illumina370","illumina550","affymetrix6","humanomniexpress","humancoreexome","humanomni1","quad610","gsa","illuminaexome")

for (i in 1:length(array_uc)) {
  ids<-c("HetPVal",paste("HetPVal_no_",array_uc[i],sep=""))
  tmp<-uc[,ids]
  tmp$logHetPVal<- -log10(tmp$HetPVal)
  tmp$logHetPVal_no_study<- -log10(tmp[,2])
  tmp<-tmp[,3:4]
  
  max(tmp)
  uplim<-max(tmp)
  
  p1<-ggplot(tmp, aes(logHetPVal,logHetPVal_no_study)) +
    geom_point() + ylim(0,uplim) + xlim(0,uplim) + geom_abline(intercept = 0, slope= 1, linetype="dashed", color = "red") +
    ggtitle(array_uc[i]) + xlab(" -log10(Heterogeneity Pval)\nall studies") + ylab(paste(" -log10(Heterogeneity Pval)\nall studies but ",array_uc[i],sep=""))
  # p1<-ggMarginal(tmp)
  
  assign(paste("plot_uc_",array_uc[i],sep=""),p1)
}

l = mget(ls()[grep("plot_uc",ls())])
# ggsave("~/tmp_plots/UC_HetPval_sensitibity_test.pdf", arrangeGrob(grobs = l), width = 14, height = 14)
ggsave("~/tmp_plots/UC_HetPval_sensitibity_test_new.pdf", arrangeGrob(grobs = l), width = 14, height = 14)


########### IBD

array_ibd<-c("illumina370","illumina550","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")

for (i in 1:length(array_ibd)) {
  ids<-c("HetPVal",paste("HetPVal_no_",array_ibd[i],sep=""))
  tmp<-ibd[,ids]
  tmp$logHetPVal<- -log10(tmp$HetPVal)
  tmp$logHetPVal_no_study<- -log10(tmp[,2])
  tmp<-tmp[,3:4]
  
  max(tmp)
  uplim<-max(tmp)
  
  p1<-ggplot(tmp, aes(logHetPVal,logHetPVal_no_study)) +
    geom_point() + ylim(0,uplim) + xlim(0,uplim) + geom_abline(intercept = 0, slope= 1, linetype="dashed", color = "red") +
    ggtitle(array_ibd[i]) + xlab(" -log10(Heterogeneity Pval)\nall studies") + ylab(paste(" -log10(Heterogeneity Pval)\nall studies but ",array_ibd[i],sep=""))
  # p1<-ggMarginal(tmp)
  
  assign(paste("plot_ibd_",array_ibd[i],sep=""),p1)
}

l = mget(ls()[grep("plot_ibd",ls())])
# ggsave("~/tmp_plots/IBD_HetPval_sensitivity_test.pdf", arrangeGrob(grobs = l), width = 14, height = 14)
ggsave("~/tmp_plots/IBD_HetPval_sensitivity_test_new.pdf", arrangeGrob(grobs = l), width = 14, height = 14)

# NOTE: NEW HETPVAL DATA CREATED WITH ILLUMINA550 DATA EXCLUDING PCA-JEWISH SAMPLES

#####################
# PLOT BETA

########### CD

array_cd<-c("illumina370","illumina550","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")

for (i in 1:length(array_cd)) {
  ids<-c("Effect",paste("Effect_no_",array_cd[i],sep=""))
  tmp<-cd[,ids]
  colnames(tmp)[2]<-"Effect_no_study"
  
  max(tmp)
  uplim<-max(tmp)
  lowlim<-min(tmp)
  
  p1<-ggplot(tmp, aes(Effect,Effect_no_study)) +
    geom_point() + ylim(lowlim,uplim) + xlim(lowlim,uplim) + geom_abline(intercept = 0, slope= 1, linetype="dashed", color = "red") +
    ggtitle(array_cd[i]) + xlab("Beta\nall studies") + ylab(paste("Beta\nall studies but ",array_cd[i],sep=""))
  # p1<-ggMarginal(tmp)
  
  assign(paste("plot_cd_",array_cd[i],sep=""),p1)
}

l = mget(ls()[grep("plot_cd",ls())])
ggsave("~/tmp_plots/CD_Beta_sensitibity_test.pdf", arrangeGrob(grobs = l), width = 14, height = 14)


########### UC

array_uc<-c("illumina370","illumina550","affymetrix6","humanomniexpress","humancoreexome","humanomni1","quad610","gsa","illuminaexome")

for (i in 1:length(array_uc)) {
  ids<-c("Effect",paste("Effect_no_",array_uc[i],sep=""))
  tmp<-uc[,ids]
  colnames(tmp)[2]<-"Effect_no_study"
  
  max(tmp)
  uplim<-max(tmp)
  lowlim<-min(tmp)
  
  p1<-ggplot(tmp, aes(Effect,Effect_no_study)) +
    geom_point() + ylim(lowlim,uplim) + xlim(lowlim,uplim) + geom_abline(intercept = 0, slope= 1, linetype="dashed", color = "red") +
    ggtitle(array_uc[i]) + xlab("Beta\nall studies") + ylab(paste("Beta\nall studies but ",array_uc[i],sep=""))
  # p1<-ggMarginal(tmp)
  
  assign(paste("plot_uc_",array_uc[i],sep=""),p1)
}

l = mget(ls()[grep("plot_uc",ls())])
ggsave("~/tmp_plots/UC_Beta_sensitibity_test.pdf", arrangeGrob(grobs = l), width = 14, height = 14)

########### IBD

array_ibd<-c("illumina370","illumina550","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")

for (i in 1:length(array_ibd)) {
  ids<-c("Effect",paste("Effect_no_",array_ibd[i],sep=""))
  tmp<-ibd[,ids]
  colnames(tmp)[2]<-"Effect_no_study"
  
  max(tmp)
  uplim<-max(tmp)
  lowlim<-min(tmp)
  
  p1<-ggplot(tmp, aes(Effect,Effect_no_study)) +
    geom_point() + ylim(lowlim,uplim) + xlim(lowlim,uplim) + geom_abline(intercept = 0, slope= 1, linetype="dashed", color = "red") +
    ggtitle(array_ibd[i]) + xlab("Beta\nall studies") + ylab(paste("Beta\nall studies but ",array_ibd[i],sep=""))
  # p1<-ggMarginal(tmp)
  
  assign(paste("plot_ibd_",array_ibd[i],sep=""),p1)
}

l = mget(ls()[grep("plot_ibd",ls())])
ggsave("~/tmp_plots/IBD_Beta_sensitivity_test.pdf", arrangeGrob(grobs = l), width = 14, height = 14)


############ PLOT as well PVal

########### CD

array_cd<-c("illumina370","illumina550","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")


for (i in 1:length(array_cd)) {
  ids<-c("PVal",paste("PVal_no_",array_cd[i],sep=""))
  tmp<-cd[,ids]
  tmp$logPVal<- -log10(tmp$PVal)
  tmp$logPVal_no_study<- -log10(tmp[,2])
  tmp<-tmp[,3:4]
  
  max(tmp)
  uplim<-max(tmp)
  
  p1<-ggplot(tmp, aes(logPVal,logPVal_no_study)) +
    geom_point() + ylim(0,uplim) + xlim(0,uplim) + geom_abline(intercept = 0, slope= 1, linetype="dashed", color = "red") +
    ggtitle(array_cd[i]) + xlab(" -log10(PVal)\nall studies") + ylab(paste(" -log10(PVal)\nall studies but ",array_cd[i],sep=""))
  # p1<-ggMarginal(tmp)
  
  assign(paste("plot_cd_PVal_",array_cd[i],sep=""),p1)
}

l = mget(ls()[grep("plot_cd_PVal",ls())])
ggsave("~/tmp_plots/CD_PVal_sensitibity_test.pdf", arrangeGrob(grobs = l), width = 14, height = 14)


########### UC

array_uc<-c("illumina370","illumina550","affymetrix6","humanomniexpress","humancoreexome","humanomni1","quad610","gsa","illuminaexome")

for (i in 1:length(array_uc)) {
  ids<-c("PVal",paste("PVal_no_",array_uc[i],sep=""))
  tmp<-uc[,ids]
  tmp$logPVal<- -log10(tmp$PVal)
  tmp$logPVal_no_study<- -log10(tmp[,2])
  tmp<-tmp[,3:4]
  
  max(tmp)
  uplim<-max(tmp)
  
  p1<-ggplot(tmp, aes(logPVal,logPVal_no_study)) +
    geom_point() + ylim(0,uplim) + xlim(0,uplim) + geom_abline(intercept = 0, slope= 1, linetype="dashed", color = "red") +
    ggtitle(array_uc[i]) + xlab(" -log10(PVal)\nall studies") + ylab(paste(" -log10(PVal)\nall studies but ",array_uc[i],sep=""))
  # p1<-ggMarginal(tmp)
  
  assign(paste("plot_uc_PVal_",array_uc[i],sep=""),p1)
}

l = mget(ls()[grep("plot_uc_PVal",ls())])
ggsave("~/tmp_plots/UC_PVal_sensitibity_test.pdf", arrangeGrob(grobs = l), width = 14, height = 14)

########### IBD

array_ibd<-c("illumina370","illumina550","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")

for (i in 1:length(array_ibd)) {
  ids<-c("PVal",paste("PVal_no_",array_ibd[i],sep=""))
  tmp<-ibd[,ids]
  tmp$logPVal<- -log10(tmp$PVal)
  tmp$logPVal_no_study<- -log10(tmp[,2])
  tmp<-tmp[,3:4]
  
  max(tmp)
  uplim<-max(tmp)
  
  p1<-ggplot(tmp, aes(logPVal,logPVal_no_study)) +
    geom_point() + ylim(0,uplim) + xlim(0,uplim) + geom_abline(intercept = 0, slope= 1, linetype="dashed", color = "red") +
    ggtitle(array_ibd[i]) + xlab(" -log10(PVal)\nall studies") + ylab(paste(" -log10(PVal)\nall studies but ",array_ibd[i],sep=""))
  # p1<-ggMarginal(tmp)
  
  assign(paste("plot_ibd_PVal_",array_ibd[i],sep=""),p1)
}

l = mget(ls()[grep("plot_ibd_PVal",ls())])
ggsave("~/tmp_plots/IBD_PVal_sensitivity_test.pdf", arrangeGrob(grobs = l), width = 14, height = 14)



vars1<-ibd[which(ibd$HetPVal<1E-10 & ibd$HetPVal_no_illumina550),"MarkerName"]
vars2<-cd[which(cd$HetPVal<1E-10 & cd$HetPVal_no_illumina550),"MarkerName"]
vars3<-cd[which(uc$HetPVal<1E-10 & uc$HetPVal_no_illumina550),"MarkerName"]

length(vars1)
# [1] 3552
length(vars2)
# [1] 2042
length(vars3)
# [1] 729

vars<-c(vars1,vars2,vars3)
vars<-vars[!duplicated(vars)]
length(vars)
# [1] 4463




###################

# double check MAF:

/software/bin/R-4.3.1

library(data.table)
library(ggpubr)
library(ggExtra)

path<-"/path/to/ibdgwas/IIBDGC/"

for (chr in 1:22) {
  tmp1<-fread(paste(path,"post_imputation/analysis/regenie/gsa/ibd/chr",chr,"_gsa_step2_ibd_eur_sex_PCs_firthse.regenie",sep=""),head=T)
  colnames(tmp1)[4:ncol(tmp1)]<-paste(colnames(tmp1)[4:ncol(tmp1)],"_gsa",sep="")
  tmp2<-fread(paste(path,"post_imputation/analysis/regenie/illumina550/ibd/chr",chr,"_illumina550_step2_ibd_eur_sex_PCs_firthse.regenie",sep=""),head=T)
  colnames(tmp2)[4:ncol(tmp2)]<-paste(colnames(tmp2)[4:ncol(tmp2)],"_illumina550",sep="")
  
  tmp<-merge(tmp1[,c("ID","A1FREQ_gsa","INFO_gsa","BETA.Y1_gsa","LOG10P.Y1_gsa")],tmp2[,c("ID","A1FREQ_illumina550","INFO_illumina550","BETA.Y1_illumina550","LOG10P.Y1_illumina550")],by="ID")
  rm(tmp1,tmp2)
  
  if(chr==1) {
    tmp_final<-tmp
  } else {
    tmp_final<-rbind(tmp_final,tmp)
  }
  rm(tmp)
}

tmp_final$high_het<-0
tmp_final$high_het[which(tmp_final$ID %in% vars)]<-1

tmp1<-tmp_final[which(tmp_final$ID %in% vars),]

uplim<-max(tmp1[,c("A1FREQ_gsa","A1FREQ_illumina550")])
p1<-ggplot(tmp1, aes(A1FREQ_gsa,A1FREQ_illumina550)) +
  geom_point() + ylim(0,uplim) + xlim(0,uplim) + geom_abline(intercept = 0, slope= 1, linetype="dashed", color = "red") + 
  ggtitle("Frequency\nHigh Het Variants")
uplim<-max(tmp1[,c("INFO_gsa","INFO_illumina550")])
p2<-ggplot(tmp1, aes(INFO_gsa,INFO_illumina550)) +
  geom_point() + ylim(0,uplim) + xlim(0,uplim) + geom_abline(intercept = 0, slope= 1, linetype="dashed", color = "red") + 
  ggtitle("Imputation info\nHigh Het Variants")
uplim<-max(tmp1[,c("BETA.Y1_gsa","BETA.Y1_illumina550")])
p3<-ggplot(tmp1, aes(BETA.Y1_gsa,BETA.Y1_illumina550)) +
  geom_point() + ylim(0,uplim) + xlim(0,uplim) + geom_abline(intercept = 0, slope= 1, linetype="dashed", color = "red") + 
  ggtitle("Beta\nHigh Het Variants")
uplim<-max(tmp1[,c("LOG10P.Y1_gsa","LOG10P.Y1_illumina550")])
p4<-ggplot(tmp1, aes(LOG10P.Y1_gsa,LOG10P.Y1_illumina550)) +
  geom_point() + ylim(0,uplim) + xlim(0,uplim) + geom_abline(intercept = 0, slope= 1, linetype="dashed", color = "red") + 
  ggtitle("Log10P\nHigh Het Variants")


l = mget(ls()[grep("^p[0-9]{1}",ls())])
ggsave("~/tmp_plots/IBD_4K_variants_HetPVal_lower_1E10.pdf", arrangeGrob(grobs = l), width = 14, height = 14)

###########################

head(tmp1[which(tmp1$BETA.Y1_illumina550>3 & tmp1$A1FREQ_gsa>0.01),])
# ID A1FREQ_gsa INFO_gsa BETA.Y1_gsa LOG10P.Y1_gsa
# 1: chr13:29006927:G:A  0.0141523 0.989452   0.0225546     0.1341040
# 2: chr13:29006962:C:A  0.0140573 0.989665   0.0165484     0.0947617
# A1FREQ_illumina550 INFO_illumina550 BETA.Y1_illumina550
# 1:          0.0370368         0.907174             3.06989
# 2:          0.0369189         0.908447             3.06395
# LOG10P.Y1_illumina550 high_het
# 1:               170.085        1
# 2:               169.250        1

# 
# Sequence name	Change
# GRCh37.p13 chr 13	NC_000013.10:g.29581064G>A
# GRCh38.p12 chr 13	NC_000013.11:g.29006927G>A

grep "13:29581064" ibd_build37_59957_20161107.txt 








