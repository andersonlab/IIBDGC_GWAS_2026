# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# UK-IIBGBD:

# See how variants were defined in:

# create_list_common_variants_among_cohorts.R

##################################################################
# 1.- SUBSET LIST OF COMMON VARIANTS FROM 1000GP FILES AND MERGE #
##################################################################

path=/path/to/ibdgwas/IIBDGC/

for chr in {1..22}; do /path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path}pre_imputation/QC/1000gp/1000GP_chr${chr}_b37.bed \
--bim ${path}pre_imputation/QC/1000gp/1000GP_chr${chr}_b37_edited.bim \
--fam ${path}pre_imputation/QC/1000gp/1000GP_chr${chr}_b37.fam \
--extract ${path}pre_imputation/QC/relatedness/list_common_variants_among_cohorts_iibdgc \
--make-bed --out ${path}pre_imputation/QC/1000gp/1000GP_chr${chr}_b37_subset_iibdgc;done



#####################
# 2.- COMBINE FILES #
#####################

#### R

# no chrx or chry
dat<-matrix(ncol=3,nrow=21)
dat<-as.data.frame(dat)
for (i in 1:nrow(dat)){
  dat[i,1]<-paste(path,"pre_imputation/QC/1000gp/1000GP_chr",i+1,"_b37_subset_iibdgc.bed",sep="")
  dat[i,2]<-paste(path,"pre_imputation/QC/1000gp/1000GP_chr",i+1,"_b37_subset_iibdgc.bim",sep="")
  dat[i,3]<-paste(path,"pre_imputation/QC/1000gp/1000GP_chr",i+1,"_b37_subset_iibdgc.fam",sep="")
}

write.table(dat,paste(path,"pre_imputation/QC/1000gp/list_1000GP_files_subset_iibdgc_tomerge.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#######

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/1000gp/1000GP_chr1_b37_subset_iibdgc \
--merge-list ${path}pre_imputation/QC/1000gp/list_1000GP_files_subset_iibdgc_tomerge.txt \
--allow-no-sex \
--make-bed --out ${path}pre_imputation/QC/1000gp/1000GP_all_b37_subset_iibdgc
# 4036 variants and 2504 people pass filters and QC.
# Note: No phenotypes present.

  

#######################
# 3.- MERGE ALL FILES #
#######################

#### R

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"chop_old_gwas","german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_cd_old_gwas")

length(cohorts)
# [1] 27

dat<-matrix(ncol=3,nrow=length(cohorts))
dat<-as.data.frame(dat)

for (i in 1:nrow(dat)){
  dat[i,1]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset.bed",sep="")
  dat[i,2]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset.bim",sep="")
  dat[i,3]<-paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset.fam",sep="")
}

write.table(dat,paste(path,"pre_imputation/QC/pca_1000gp/list_cohorts_tomerge_iibdgc.txt",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

#############################

# hpc-server

/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path}pre_imputation/QC/1000gp/1000GP_all_b37_subset_iibdgc \
--allow-no-sex \
--merge-list ${path}pre_imputation/QC/pca_1000gp/list_cohorts_tomerge_iibdgc.txt \
--make-bed --out ${path}pre_imputation/QC/pca_1000gp/1000GP_iibdgc_merged
# 4036 variants and 102641 people pass filters and QC.
# Among remaining phenotypes, 56391 are cases and 43746 are controls.  (2504
#                                                                       phenotypes are missing.)

###############
# 4.- RUN PCA #
###############


path=/path/to/ibdgwas/IIBDGC/
MEM=5000 # next time only this

bsub -J"bcftools" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal -n 10 \
-e ${path}pre_imputation/QC/pca_1000gp/logs/stderr_pca_iibdgc_1000gp \
-o ${path}pre_imputation/QC/pca_1000gp/logs/stdout_pca_iibdgc_1000gp \
"/path/to/software/./plink2  \
--bfile ${path}pre_imputation/QC/pca_1000gp/1000GP_iibdgc_merged \
--pca approx --out ${path}pre_imputation/QC/pca_1000gp/1000GP_iibdgc_merged_pca \
--threads 8 --allow-no-sex --memory $MEM"
# Job <216531> is submitted to queue <normal>.


# using what is recommented:
# The 'approx' modifier causes the standard deterministic computation to be replaced with the randomized algorithm originally implemented for 
# Galinsky KJ, Bhatia G, Loh PR, Georgiev S, Mukherjee S, Patterson NJ, Price AL (2016) Fast Principal-Component Analysis Reveals 
# Convergent Evolution of ADH1B in Europe and East Asia. This can be a good idea when you have >5000 samples, and is almost required once you have >50000.

# CPU time :                                   66.40 sec.
# Max Memory :                                 6047 MB
# Average Memory :                             4217.25 MB
# Total Requested Memory :                     50000.00 MB
# Delta Memory :                               43953.00 MB
# Max Swap :                                   -
# Max Processes :                              3
# Max Threads :                                19
# Run time :                                   28 sec.
# Turnaround time :                            20 sec.


###### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(viridis)

path<-"/path/to/ibdgwas/IIBDGC/"


# no german_illu
cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"chop_old_gwas","german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_cd_old_gwas")

length(cohorts)
# [1] 27

# needs to be QCed:
# "swedish_uc_old_gwas"
# "niddk_uc_old_gwas"
# "cedars"

for (i in 1:length(cohorts)){
  a<-read.table(paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset.fam",sep=""),head=F)
  a$cohort<-cohorts[i]
  if(i==1){
    fam<-a
  }else{
    fam<-rbind(a,fam)
  }
}

table(fam$cohort)

samples1000<-read.table("/path/to/project",head=T)
table(samples1000$pop,samples1000$super_pop)
#     AFR AMR EAS EUR SAS
# ACB  96   0   0   0   0
# ASW  61   0   0   0   0
# BEB   0   0   0   0  86
# CDX   0   0  93   0   0
# CEU   0   0   0  99   0
# CHB   0   0 103   0   0
# CHS   0   0 105   0   0
# CLM   0  94   0   0   0
# ESN  99   0   0   0   0
# FIN   0   0   0  99   0
# GBR   0   0   0  91   0
# GIH   0   0   0   0 103
# GWD 113   0   0   0   0
# IBS   0   0   0 107   0
# ITU   0   0   0   0 102
# JPT   0   0 104   0   0
# KHV   0   0  99   0   0
# LWK  99   0   0   0   0
# MSL  85   0   0   0   0
# MXL   0  64   0   0   0
# PEL   0  85   0   0   0
# PJL   0   0   0   0  96
# PUR   0 104   0   0   0
# STU   0   0   0   0 102
# TSI   0   0   0 107   0
# YRI 108   0   0   0   0

fam<-fam[,c(1,7)]
colnames(fam)[1]<-c("sample")
fam$super_pop<-"IIBDGC"
fam$pop<-NA

samples1000$cohort<-"1000GP"
samples1000<-samples1000[,c("sample","cohort","super_pop","pop")]
fam<-rbind(fam,samples1000)
dim(fam)
# [1] 102641      3


pca<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/1000GP_iibdgc_merged_pca.eigenvec",sep=""),head=F)
dim(pca)
# [1] 102641     12
colnames(pca)<-c("FID","IID",paste("PC",seq(1:10),sep=""))

pca<-merge(pca,fam[,c("sample","cohort","super_pop","pop")],by.x="FID",by.y="sample",all.x=T)


table(pca$super_pop)
# AFR    AMR    EAS    EUR IIBDGC    SAS 
# 661    347    504    503 100137    489 

# TOTAL VARIANCE EXPLAINED BY EACH PC:

eigenval<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/1000GP_iibdgc_merged_pca.eigenval",sep=""),head=F)

eigenval$var_exp<-NA
for (i in 1:nrow(eigenval)){
  eigenval$var_exp[i]<-(eigenval$V1[i] / sum(eigenval$V1))*100
}

eigenval
#          V1   var_exp
# PC1  994.4520 43.743286
# PC2  444.6950 19.560945
# PC3  250.3040 11.010204
# PC4  114.7790  5.048822
# PC5   94.0627  4.137567
# PC6   89.5257  3.937996
# PC7   81.2687  3.574793
# PC8   72.5730  3.192292
# PC9   69.2292  3.045208
# PC10  62.4927  2.748887

# PC1 EUR - non EUR
# PC2 AFR  - non AFR
# PC3 north EUR -South EUR
# PC4 EAS - SAS
# PC5 Finish - non-Finish
# PC6: AMR - non-AMR
# PC7 - IIBDGC outliers (from more than 1 cohort)
# PC8 transition EUR - EAS - SAS/AFR - AMR

table(pca[which(pca$PC7< -0.05),"cohort"])
# all_hce belgium_vermeire_gsa       prism_nfe_gwas 
# 11                    2                    6 


cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"chop_old_gwas","german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_cd_old_gwas")

pca$study<-NA

pca$study[which(pca$cohort=="australia_omniexome")]<-"Australia"
pca$study[which(pca$cohort %in% c("gwas1","gwas2","all_hce"))]<-"UK"
pca$study[which(pca$cohort %in% c("pittsburgh_gsa"))]<-"Pittsburgh"
pca$study[which(pca$cohort %in% c("spain_gsa","basque_gsa"))]<-"Spain"
pca$study[which(pca$cohort %in% c("kiel_austria_sibdcs_gsa","german_affy6_old_gwas"))]<-"Germany"                 
pca$study[which(pca$cohort %in% c("italy_gsa"))]<-"Italy"   
pca$study[which(pca$cohort %in% c("netherlands_gsa"))]<-"Netherlands"                     
pca$study[which(pca$cohort %in% c("slovenia_gsa"))]<-"Slovenia"                     
pca$study[which(pca$cohort %in% c("sweden_gsa"))]<-"Sweden"                    
pca$study[which(pca$cohort %in% c("niddk_broad_gsa","niddk_feinstein_gsa","niddk_cd_old_gwas"))]<-"NIDDK"              
pca$study[which(pca$cohort %in% c("finland_illugwas"))]<-"Finland"
pca$study[which(pca$cohort %in% c("chop_old_gwas"))]<-"CHOP"
pca$study[which(pca$cohort %in% c("norway_affy6_old_gwas"))]<-"Norway"
pca$study[which(pca$cohort %in% c("belgium_inf1_old_gwas","belgium_inf2_old_gwas","belgium_louis_gsa","belgium_franchimont_gsa"
                                  ,"belgium_vermeire_gsa"))]<-"Belgium"
pca$study[which(pca$cohort %in% c("lithuania_gsa"))]<-"Lithuania"
pca$study[which(pca$cohort %in% c("prism_nfe_gwas","prism_nfe_gsa"))]<-"PRISM"




####################
# 4.- PLOT RESULTS #
####################

######################################################
# 4.1 PLOT 1000GP ALONE, COLOUR BY SUPER POPULATION:

p0<-qplot(pca$PC1,pca$PC2, data = pca, colour = super_pop)
pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2.pdf",sep=""),width=10,height=7)
ggarrange(p0,legend=c("right"))
dev.off()

p0<-qplot(pca$PC2[which(pca$cohort=="1000GP")],pca$PC3[which(pca$cohort=="1000GP")], data = pca, colour = super_pop)
pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC2_PC3.pdf",sep=""),width=10,height=7)
ggarrange(p0,legend=c("right"))
dev.off()

p0<-qplot(pca$PC2[which(pca$super_pop=="EUR")],pca$PC3[which(pca$super_pop=="EUR")], data = pca[which(pca$super_pop=="EUR"),], colour = pop)
pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_only_PC2_PC3.pdf",sep=""),width=10,height=7)
ggarrange(p0,legend=c("right"))
dev.off()


p0<-qplot(pca$PC3,pca$PC4, data = pca, colour = super_pop)
pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC3_PC4.pdf",sep=""),width=10,height=7)
ggarrange(p0,legend=c("right"))
dev.off()

p0<-qplot(pca$PC4,pca$PC5, data = pca, colour = super_pop)
pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC4_PC5.pdf",sep=""),width=10,height=7)
ggarrange(p0,legend=c("right"))
dev.off()

p0<-qplot(pca$PC4[which(pca$super_pop=="EUR")],pca$PC5[which(pca$super_pop=="EUR")], data = pca[which(pca$super_pop=="EUR"),], colour = pop)
pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_only_PC4_PC5.pdf",sep=""),width=10,height=7)
ggarrange(p0,legend=c("right"))
dev.off()


p0<-qplot(pca$PC5[which(pca$cohort=="1000GP")],pca$PC6[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop)
pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_only_PC5_PC6.pdf",sep=""),width=10,height=7)
ggarrange(p0,legend=c("right"))
dev.off()

p0<-qplot(pca$PC6,pca$PC7, data = pca, colour = super_pop)
pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC6_PC7.pdf",sep=""),width=10,height=7)
ggarrange(p0,legend=c("right"))
dev.off()

p0<-qplot(pca$PC7,pca$PC8, data = pca, colour = super_pop)
pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC7_PC8.pdf",sep=""),width=10,height=7)
ggarrange(p0,legend=c("right"))
dev.off()

p0<-qplot(pca$PC8,pca$PC9, data = pca, colour = super_pop)
pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC8_PC9.pdf",sep=""),width=10,height=7)
ggarrange(p0,legend=c("right"))
dev.off()

p0<-qplot(pca$PC8[which(pca$cohort=="1000GP")],pca$PC9[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop)
pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_only_PC8_PC9.pdf",sep=""),width=10,height=7)
ggarrange(p0,legend=c("right"))
dev.off()


p0<-qplot(pca$PC9,pca$PC10, data = pca, colour = super_pop)
pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC9_PC10.pdf",sep=""),width=10,height=7)
ggarrange(p0,legend=c("right"))
dev.off()

p0<-qplot(pca$PC9[which(pca$cohort=="1000GP")],pca$PC10[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop)
pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_only_PC9_PC10.pdf",sep=""),width=10,height=7)
ggarrange(p0,legend=c("right"))
dev.off()

# creat plots PC1 to PC6, 1000GP only
pca$pop<-as.factor(pca$pop)
pca$pop<-factor(pca$pop, levels=c("CEU","FIN","GBR","IBS","TSI",
                                  "CDX","CHB","CHS","JPT","KHV",
                                  "BEB","GIH","ITU","PJL","STU",
                                  "ACB","ASW","ESN","GWD","LWK","MSL","YRI",
                                  "CLM","MXL","PEL","PUR"))

# EUR = "CEU","FIN","GBR","IBS","TSI"
# EAS = "CDX","CHB","CHS","JPT","KHV"
# SAS =  "BEB","GIH","ITU","PJL","STU"
# AFR = "ACB","ASW","ESN","GWD","LWK","MSL","YRI"
# AMF = "CLM","MXL","PEL","PUR"

pca$super_pop<-as.factor(pca$super_pop)
pca$super_pop<-factor(pca$super_pop, levels=c("EUR","AFR","EAS","SAS","AMR","IIBDGC"))

  
pna<-qplot()+theme(
  panel.background = element_rect(fill = "transparent") # bg of the panel
  , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
  , panel.grid.major = element_blank() # get rid of major grid
  , panel.grid.minor = element_blank() # get rid of minor grid
  , legend.background = element_rect(fill = "transparent") # get rid of legend bg
  , axis.title=element_blank()
  , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
)

p11<-qplot(pca$PC1[which(pca$cohort=="1000GP")],pca$PC1[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p12<-qplot(pca$PC2[which(pca$cohort=="1000GP")],pca$PC1[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p13<-qplot(pca$PC3[which(pca$cohort=="1000GP")],pca$PC1[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p14<-qplot(pca$PC4[which(pca$cohort=="1000GP")],pca$PC1[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p15<-qplot(pca$PC5[which(pca$cohort=="1000GP")],pca$PC1[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p16<-qplot(pca$PC6[which(pca$cohort=="1000GP")],pca$PC1[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")

p21<-qplot(pca$PC1[which(pca$cohort=="1000GP")],pca$PC2[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p22<-qplot(pca$PC2[which(pca$cohort=="1000GP")],pca$PC2[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p23<-qplot(pca$PC3[which(pca$cohort=="1000GP")],pca$PC2[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p24<-qplot(pca$PC4[which(pca$cohort=="1000GP")],pca$PC2[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p25<-qplot(pca$PC5[which(pca$cohort=="1000GP")],pca$PC2[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p26<-qplot(pca$PC6[which(pca$cohort=="1000GP")],pca$PC2[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")

p31<-qplot(pca$PC1[which(pca$cohort=="1000GP")],pca$PC3[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p32<-qplot(pca$PC2[which(pca$cohort=="1000GP")],pca$PC3[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p33<-qplot(pca$PC3[which(pca$cohort=="1000GP")],pca$PC3[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p34<-qplot(pca$PC4[which(pca$cohort=="1000GP")],pca$PC3[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p35<-qplot(pca$PC5[which(pca$cohort=="1000GP")],pca$PC3[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p36<-qplot(pca$PC6[which(pca$cohort=="1000GP")],pca$PC3[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")

p41<-qplot(pca$PC1[which(pca$cohort=="1000GP")],pca$PC4[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p42<-qplot(pca$PC2[which(pca$cohort=="1000GP")],pca$PC4[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p43<-qplot(pca$PC3[which(pca$cohort=="1000GP")],pca$PC4[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p44<-qplot(pca$PC4[which(pca$cohort=="1000GP")],pca$PC4[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p45<-qplot(pca$PC5[which(pca$cohort=="1000GP")],pca$PC4[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p46<-qplot(pca$PC6[which(pca$cohort=="1000GP")],pca$PC4[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")

p51<-qplot(pca$PC1[which(pca$cohort=="1000GP")],pca$PC5[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p52<-qplot(pca$PC2[which(pca$cohort=="1000GP")],pca$PC5[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p53<-qplot(pca$PC3[which(pca$cohort=="1000GP")],pca$PC5[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p54<-qplot(pca$PC4[which(pca$cohort=="1000GP")],pca$PC5[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p55<-qplot(pca$PC5[which(pca$cohort=="1000GP")],pca$PC5[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p56<-qplot(pca$PC6[which(pca$cohort=="1000GP")],pca$PC5[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")

p61<-qplot(pca$PC1[which(pca$cohort=="1000GP")],pca$PC6[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p62<-qplot(pca$PC2[which(pca$cohort=="1000GP")],pca$PC6[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p63<-qplot(pca$PC3[which(pca$cohort=="1000GP")],pca$PC6[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p64<-qplot(pca$PC4[which(pca$cohort=="1000GP")],pca$PC6[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p65<-qplot(pca$PC5[which(pca$cohort=="1000GP")],pca$PC6[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")
p66<-qplot(pca$PC6[which(pca$cohort=="1000GP")],pca$PC6[which(pca$cohort=="1000GP")], data = pca[which(pca$cohort=="1000GP"),], colour = super_pop) + scale_color_viridis(discrete = TRUE, option = "D")


r1<-ggarrange(pna,p12,p13,p14,p15,p16,ncol=6,legend=c("none"))
r2<-ggarrange(p21,pna,p23,p24,p25,p26,ncol=6,legend=c("none"))
r3<-ggarrange(p31,p32,pna,p34,p35,p36,ncol=6,legend=c("none"))
r4<-ggarrange(p41,p42,p43,pna,p45,p36,ncol=6,legend=c("none"))
r5<-ggarrange(p51,p52,p53,p54,pna,p56,ncol=6,legend=c("none"))
r6<-ggarrange(p61,p62,p63,p64,p65,pna,ncol=6,common.legend = TRUE,legend=c("bottom"))
dev.off()


pdf(paste(path,"pre_imputation/QC/pca_1000gp/pca_PC1_PC2_PC3_PC5_PC6_3kvariants.pdf",sep=""),width=30,height=31)
ggarrange(r1,r2,r3,r4,r5,r6,nrow=6,common.legend = TRUE,legend=c("bottom"),heights = c(1,1,1,1,1,1.2))
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/pca_PC1_PC2_PC3_PC5_PC6_3kvariants.pdf ~/",sep=""))



#######################################################################
# 4.1 PLOT 1000GP ALONE, AND SEPARATELY , COLOUR BY SUPER POPULATION:


# PC1 EUR - non EUR
# PC2 AFR  - non AFR
# PC3 north EUR -South EUR
# PC4 EAS - SAS
# PC5 Finish - non-Finish
# PC6: AMR - non-AMR
# PC7 - IIBDGC outliers (from more than 1 cohort)
# PC8 transition EUR - EAS - SAS/AFR - AMR


[1]                          
[6]                        
[11]                        
[16]   


pca$study<-as.factor(pca$study)
pca$study<-factor(pca$study, levels=c("Spain","Italy","Slovenia","Belgium","Germany","UK","Netherlands","Lithuania","Norway","Sweden","Finland",
                                    "CHOP","NIDDK","Pittsburgh","PRISM","Australia"))

# PC1 - PC2

ylims<-c(min(pca$PC1),max(pca$PC1))
xlims<-c(min(pca$PC2),max(pca$PC2))

p0<-ggplot(pca[which(pca$cohort=="1000GP"),],aes(PC2,PC1)) +
  geom_point(aes(color = super_pop)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims)

p1<-ggplot(pca[which(pca$cohort!="1000GP"),],aes(PC2,PC1)) +
  geom_point(aes(color = study)) + 
  xlim(xlims) + ylim(ylims)


pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2.pdf",sep=""),width=14,height=7)
ggarrange(p0,p1,legend=c("right"),widths = c(1,1.1))
dev.off()

for (i in 1:length(levels(pca$study))) {
  
  print(i)
  p1<-ggplot(pca[which(pca$study==levels(pca$study)[i]),],aes(PC2,PC1,color=study)) +
    geom_point() + 
    xlim(xlims) + ylim(ylims)
  
  pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC2_",levels(pca$study)[i],".pdf",sep=""),width=14,height=7)
  print(ggarrange(p0,p1,legend=c("right"),widths = c(1,1.1)))
  dev.off()
  
}

# PC1 - PC3

ylims<-c(min(pca$PC1),max(pca$PC1))
xlims<-c(min(pca$PC3),max(pca$PC3))

p0<-ggplot(pca[which(pca$super_pop=="1000GP"),],aes(PC3,PC1)) +
  geom_point(aes(color = pop)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims)

p1<-ggplot(pca[which(pca$cohort!="1000GP"),],aes(PC3,PC1)) +
  geom_point(aes(color = study)) + 
  xlim(xlims) + ylim(ylims)


pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC3.pdf",sep=""),width=14,height=7)
ggarrange(p0,p1,legend=c("right"),widths = c(1,1))
dev.off()

#####

ylims<-c(min(pca$PC1),max(pca$PC1))
xlims<-c(min(pca$PC3),max(pca$PC3))

p0<-ggplot(pca[which(pca$super_pop=="EUR"),],aes(PC3,PC1)) +
  geom_point(aes(color = pop)) + scale_color_viridis(discrete = TRUE, option = "D") + 
  xlim(xlims) + ylim(ylims)

p1<-ggplot(pca[which(pca$cohort!="1000GP"),],aes(PC3,PC1)) +
  geom_point(aes(color = study)) + 
  xlim(xlims) + ylim(ylims)

pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_EUR_PC1_PC3.pdf",sep=""),width=14,height=7)
ggarrange(p0,p1,legend=c("right"),widths = c(1,1))
dev.off()

for (i in 1:length(levels(pca$study))) {
  
  print(i)
  p1<-ggplot(pca[which(pca$study==levels(pca$study)[i]),],aes(PC3,PC1,color=study)) +
    geom_point() + 
    xlim(xlims) + ylim(ylims)
  
  pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_1000GP_PC1_PC3_",levels(pca$study)[i],".pdf",sep=""),width=14,height=7)
  print(ggarrange(p0,p1,legend=c("right"),widths = c(1,1.1)))
  dev.off()
  
}




# OLD
# 
# 
# 
# peur1<-qplot(pca$PC3[which(pca$super_pop %in% c("EUR","GWAS1","GWAS2","GWAS3","New_Wave","Slovenia","Sweden","Netherlands"))],
#             pca$PC6[which(pca$super_pop %in% c("EUR","GWAS1","GWAS2","GWAS3","New_Wave","Slovenia","Sweden","Netherlands"))],
#             data = pca[which(pca$super_pop %in% c("EUR","GWAS1","GWAS2","GWAS3","New_Wave","Slovenia","Sweden","Netherlands")),], colour = pop) + 
#   scale_color_manual(values=c("FIN"="#662506","CEU"="#cc4c02","GBR"="#fe9929","IBS"="#fee391","TSI"="#fff7bc",
#                               "GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04"),na.value="darkgrey",name="pop") + 
#   xlim(min(pca$PC3),max(pca$PC3)) + 
#   ylim(min(pca$PC6),max(pca$PC6)) + 
#   xlab("PC3") + 
#   ylab("PC6")
# 
# peur2<-qplot(pca$PC3[which(pca$super_pop %in% c("EUR"))],
#             pca$PC6[which(pca$super_pop %in% c("EUR"))],
#             data = pca[which(pca$super_pop %in% c("EUR")),], colour = pop) + 
#   scale_color_manual(values=c("FIN"="#662506","CEU"="#cc4c02","GBR"="#fe9929","IBS"="#fee391","TSI"="#fff7bc"
#                               ),na.value="darkgrey",name="pop") + 
#   xlim(min(pca$PC3),max(pca$PC3)) + 
#   ylim(min(pca$PC6),max(pca$PC6)) + 
#   xlab("PC3") + 
#   ylab("PC6")
# 
# pdf(paste(path,"pre_imputation/QC/pca_1000gp/pca_gwas1_gwas1_gwas3_new_wave_gsa_1000GP_PC3_PC6.pdf",sep=""),width=12,height=8)
# ggarrange(peur1,peur2,ncol=2)
# dev.off()
# 
# peur1<-qplot(pca$PC1[which(pca$super_pop %in% c("EUR","GWAS1","GWAS2","GWAS3","New_Wave","Slovenia","Sweden","Netherlands"))],
#              pca$PC6[which(pca$super_pop %in% c("EUR","GWAS1","GWAS2","GWAS3","New_Wave","Slovenia","Sweden","Netherlands"))],
#              data = pca[which(pca$super_pop %in% c("EUR","GWAS1","GWAS2","GWAS3","New_Wave","Slovenia","Sweden","Netherlands")),], colour = pop) + 
#   scale_color_manual(values=c("FIN"="#662506","CEU"="#cc4c02","GBR"="#fe9929","IBS"="#fee391","TSI"="#fff7bc",
#                               "GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04"),na.value="darkgrey",name="pop") + 
#   xlim(min(pca$PC1),max(pca$PC1)) + 
#   ylim(min(pca$PC6),max(pca$PC6)) + 
#   xlab("PC1") + 
#   ylab("PC6")
# 
# peur2<-qplot(pca$PC1[which(pca$super_pop %in% c("EUR"))],
#              pca$PC6[which(pca$super_pop %in% c("EUR"))],
#              data = pca[which(pca$super_pop %in% c("EUR")),], colour = pop) + 
#   scale_color_manual(values=c("FIN"="#662506","CEU"="#cc4c02","GBR"="#fe9929","IBS"="#fee391","TSI"="#fff7bc"
#   ),na.value="darkgrey",name="pop") + 
#   xlim(min(pca$PC1),max(pca$PC1)) + 
#   ylim(min(pca$PC6),max(pca$PC6)) + 
#   xlab("PC1") + 
#   ylab("PC6")
# 
# pdf(paste(path,"pre_imputation/QC/pca_1000gp/pca_gwas1_gwas1_gwas3_new_wave_gsa_1000GP_PC1_PC6.pdf",sep=""),width=12,height=8)
# ggarrange(peur1,peur2,ncol=2)
# dev.off()
# 
# 
# #V3 thershold:
# table(pca[which(pca$V3< 0.0025),"super_pop"]) # OK
# #GWAS1    GWAS2    GWAS3 New_Wave      EUR      AFR      EAS      SAS 
# #4676     7778    20885     6299      503        0        0        0
# 
# table(pca[which( pca$V4< 0.0035 & pca$V4> -0.0035),"super_pop"]) # OK
# #GWAS1    GWAS2    GWAS3 New_Wave      EUR      AFR      EAS      SAS 
# # 4679     7778    20990     6342      503        2        0        0 
# 
# table(pca[which( pca$V5< 0.005 & pca$V5> -0.005),"super_pop"]) # OK
# #GWAS1    GWAS2    GWAS3 New_Wave      EUR      AFR      EAS      SAS 
# #4649     7778    20896     6326      503      650        0        0
# 
# table(pca[which(pca$V3< 0.0025 & pca$V4< 0.005 & pca$V4> -0.005 & pca$V5< 0.005 & pca$V5> -0.005),"super_pop"])
# #GWAS1    GWAS2    GWAS3 New_Wave      EUR      AFR      EAS      SAS 
# # 4648     7778    20707     6254      503        0        0        0 
# 
# 
# #########################
# # now plot just 1000GP
# 
# pcaall<-pca
# pca<-pca[which(is.na(pca$study)),]
# 
# pna<-qplot()+theme(
#   panel.background = element_rect(fill = "transparent") # bg of the panel
#   , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
#   , panel.grid.major = element_blank() # get rid of major grid
#   , panel.grid.minor = element_blank() # get rid of minor grid
#   , legend.background = element_rect(fill = "transparent") # get rid of legend bg
#   , axis.title=element_blank()
#   , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
# )
# 
# 
# p12<-qplot(pca$PC2,pca$PC1, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# 
# p13<-qplot(pca$PC3,pca$PC1, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p14<-qplot(pca$PC4,pca$PC1, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p15<-qplot(pca$PC5,pca$PC1, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p16<-qplot(pca$PC6,pca$PC1, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p17<-qplot(pca$PC7,pca$PC1, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p18<-qplot(pca$PC8,pca$PC1, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p19<-qplot(pca$PC9,pca$PC1, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# 
# 
# pc1_pcall<-ggarrange(p12,p13,p14,p15,p16,p17,p18,p19,nrow=2,ncol=4,common.legend = TRUE,legend=c("bottom"))
# 
# pdf(paste(path,"pre_imputation/QC/pca_1000gp/pca_gwas1_gwas1_gwas3_new_wave_gsa_1000GP_PC1_vs_other.pdf",sep=""),width=12,height=8)
# pc1_pcall
# dev.off()
# 
# # UP to PC6
# 
# p21<-qplot(pca$PC1,pca$PC2, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p23<-qplot(pca$PC3,pca$PC2, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p24<-qplot(pca$PC4,pca$PC2, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p25<-qplot(pca$PC5,pca$PC2, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p26<-qplot(pca$PC6,pca$PC2, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# 
# 
# p31<-qplot(pca$PC1,pca$PC3, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p32<-qplot(pca$PC2,pca$PC3, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p34<-qplot(pca$PC4,pca$PC3, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p35<-qplot(pca$PC5,pca$PC3, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p36<-qplot(pca$PC6,pca$PC3, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# 
# p41<-qplot(pca$PC1,pca$PC4, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p42<-qplot(pca$PC2,pca$PC4, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p43<-qplot(pca$PC3,pca$PC4, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p45<-qplot(pca$PC5,pca$PC4, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p46<-qplot(pca$PC6,pca$PC4, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# 
# 
# p51<-qplot(pca$PC1,pca$PC5, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p52<-qplot(pca$PC2,pca$PC5, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p53<-qplot(pca$PC3,pca$PC5, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p54<-qplot(pca$PC4,pca$PC5, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p56<-qplot(pca$PC6,pca$PC5, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# 
# p61<-qplot(pca$PC1,pca$PC6, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p62<-qplot(pca$PC2,pca$PC6, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p63<-qplot(pca$PC3,pca$PC6, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p64<-qplot(pca$PC4,pca$PC6, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# p65<-qplot(pca$PC5,pca$PC6, data = pca, colour = super_pop) + 
#   # geom_vline(xintercept=0.0025, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=-0.0035, linetype="dashed",color = "darkgrey") + 
#   # geom_hline(yintercept=0.0035, linetype="dashed",color = "darkgrey") + 
#   scale_color_manual(values=c("GWAS1"="#a6cee3","GWAS2"="#4497cf","GWAS3"="#1f78b4","New_Wave"="#0c69a8",
#                               "Slovenia"="#b2df8a","Sweden"="#33a02c","Netherlands"="#0b7d04",
#                               "EUR"="#e31a1c","AFR"="#cab2d6","AMR"="#fb9a99","EAS"="#fdbf6f","SAS"="#ff7f00","AMR"="#fb9a99"),na.value="darkgrey",name="super_pop")
# 
# 
# 
# 
# r1<-ggarrange(pna,p12,p13,p14,p15,p16,ncol=6,legend=c("none"))
# r2<-ggarrange(p21,pna,p23,p24,p25,p26,ncol=6,legend=c("none"))
# r3<-ggarrange(p31,p32,pna,p34,p35,p36,ncol=6,legend=c("none"))
# r4<-ggarrange(p41,p42,p43,pna,p45,p36,ncol=6,legend=c("none"))
# r5<-ggarrange(p51,p52,p53,p54,pna,p56,ncol=6,legend=c("none"))
# r6<-ggarrange(p61,p62,p63,p64,p65,pna,ncol=6,common.legend = TRUE,legend=c("bottom"))
# dev.off()
# 
# pdf(paste(path,"pre_imputation/QC/pca_1000gp/pca_1000GP_PC1_PC2_PC3_PC5_PC6.pdf",sep=""),width=18,height=19)
# ggarrange(r1,r2,r3,r4,r5,r6,nrow=6,common.legend = TRUE,legend=c("bottom"),heights = c(1,1,1,1,1,1.2))
# dev.off()
# 
# 
# ############
# # EUROPEAN #
# ############
# 
# eur_samples<-pca[which(pca$V3< 0.0025 & pca$V4< 0.005 & pca$V4> -0.005 & pca$V5< 0.005 & pca$V5> -0.005),c("V1","V2","super_pop")]
# eur_samples<-eur_samples[which(eur_samples$super_pop %in% c("GWAS1","GWAS2","GWAS3","New_Wave")),]
# table(eur_samples$super_pop)
# #GWAS1    GWAS2    GWAS3 New_Wave      EUR      AFR      EAS      SAS 
# # 4648     7778    20707     6254        0        0        0        0 
# colnames(eur_samples)[1:2]<-c("FID","IID")
# write.table(eur_samples[,1:2],paste(path,"pre_imputation/QC/pca_1000gp/list_european_samples",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
# 
# 
# #######
# # SAS #
# #######
# 
# sas_samples<-pca[which(pca$V4> 0.0049  & pca$V5> 0.02),c("V1","V2","super_pop")]
# table(sas_samples$super_pop)
# #GWAS1    GWAS2    GWAS3 New_Wave      EUR      AFR      EAS      SAS 
# #    0        0      406       96        0        0        0      488 # captures all SAS
# sas_samples<-sas_samples[which(sas_samples$super_pop %in% c("GWAS1","GWAS2","GWAS3","New_Wave")),]
# colnames(sas_samples)[1:2]<-c("FID","IID")
# write.table(sas_samples[,1:2],paste(path,"pre_imputation/QC/pca_1000gp/list_south_asian_samples",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
# 
# 
# #######
# # AFR #
# #######
# 
# afr_samples<-pca[which(pca$V4< -0.005  & pca$V5< 0.005 & pca$V3> 0.0125),c("V1","V2","super_pop")]
# table(afr_samples$super_pop)
# #GWAS1    GWAS2    GWAS3 New_Wave      EUR      AFR      EAS      SAS 
# #    0        0       85       35        0      655        0        0  # loosing 6 (AFR= 661)
# afr_samples<-afr_samples[which(afr_samples$super_pop %in% c("GWAS1","GWAS2","GWAS3","New_Wave")),]
# colnames(afr_samples)[1:2]<-c("FID","IID")
# write.table(afr_samples[,1:2],paste(path,"pre_imputation/QC/pca_1000gp/list_african_samples",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
# 
# 
# samples1000_ids<-samples1000[,c(1,1)]
# colnames(samples1000_ids)[1:2]<-c("FID","IID")
# write.table(samples1000_ids,paste(path,"pre_imputation/QC/pca_1000gp/list_1000gp_samples",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
# 
# 
# # EXPLORE PROJECTIONS!!!!!
# 
# # 
# # Yes, it's likely that something like PLINK 1.9's --pca-clusters/--pca-cluster-names will eventually make it into 2.0.  With that said, PCA projection is actually already supported, the workflow is just a bit more convoluted for now.
# # 
# # Step 1: Export allele frequencies and PCA variant weights from your reference dataset.  E.g.
# # plink2 --bfile hapmap --freq --pca var-wts --out pca_hapmap
# # 
# # Step 2: Use --score to compute the necessary dot products with the variant weights.  E.g.
# # plink2 --bfile mydata --read-freq pca_hapmap.afreq --score pca_hapmap.eigenvec.var 2 3 header-read no-mean-imputation variance-normalize --score-col-nums 5-14 --out pca_proj_mydata
# # 
# # These PCs will be scaled a bit differently from pca_hapmap.eigenvec (you need to multiply or divide the PCs by sqrt(eigenvalue) to put them on the same scale).  One way to make the PCs directly comparable is to also run step 2 on the original dataset.
# # 
# 
# 
# http://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_20200409.zip