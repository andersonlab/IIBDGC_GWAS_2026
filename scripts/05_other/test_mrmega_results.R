# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# compare results using different PCs for chr22 - CD

# MEM=4000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1

library(data.table)
library(ggplot2)
library(ggpubr)

# path_gwas<-"/path/to/ibdgwas/IIBDGC/"

# pc1<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/mrmega/output_files/cd/mrmega_eur_eas_sas_tier_2_cd_chr22_pc1.result",sep=""))
# pc2<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/mrmega/output_files/cd/mrmega_eur_eas_sas_tier_2_cd_chr22_pc2.result",sep=""))
# pc3<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/mrmega/output_files/cd/mrmega_eur_eas_sas_tier_2_cd_chr22_pc3.result",sep=""))
# pc4<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/mrmega/output_files/cd/mrmega_eur_eas_sas_tier_2_cd_chr22_pc4.result",sep=""))
# pc5<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/mrmega/output_files/cd/mrmega_eur_eas_sas_tier_2_cd_chr22_pc5.result",sep=""))


# for (i in c(1:5)) {
#     print(paste("pc",i,sep=""))

#     tmp<-get(paste("pc",i,sep=""))
#     print(paste("N variants in pc",i,nrow(tmp)))
#     print(paste("N variants with mrmega results in pc",i,nrow(tmp[!is.na(tmp$beta_0),])))
#     print(table(tmp$Ncohort[is.na(tmp$beta_0)]))

#     rm(tmp)

# }

# #[1] "pc1"
# #[1] "N variants in pc 1 597136"
# #[1] "N variants with mrmega results in pc 1 222801"
# #     1      2      3 
# #289809  59257  25269 

# #[1] "pc2"
# #[1] "N variants in pc 2 597136"
# #[1] "N variants with mrmega results in pc 2 208331"
#      1      2      3      4 
# #289809  59257  25269  14470 

# #[1] "pc3"
# #[1] "N variants in pc 3 597136"
# #[1] "N variants with mrmega results in pc 3 198630"
# #     1      2      3      4      5 
# #289809  59257  25269  14470   9701 

# #[1] "pc4"
# #[1] "N variants in pc 4 597136"
# #[1] "N variants with mrmega results in pc 4 191551"
# #     1      2      3      4      5      6 
# #289809  59257  25269  14470   9701   7079 

# #[1] "pc5"
# #[1] "N variants in pc 5 597136"
# #[1] "N variants with mrmega results in pc 5 185111"
# #     1      2      3      4      5      6      7 
# #289809  59257  25269  14470   9701   7079   6440 

# # N variants with an output driven by PC input - only results when studies contributing to snp > pc+2


# # how resutls compare for the variants in all outputs?

# pc5<-pc5[!is.na(pc5$beta_0),]
# colnames(pc5)[2:ncol(pc5)]<-paste(colnames(pc5)[2:ncol(pc5)],"_pc5",sep="")
# pc4<-pc4[which(pc4$MarkerName %in% pc5$MarkerName),]
# colnames(pc4)[2:ncol(pc4)]<-paste(colnames(pc4)[2:ncol(pc4)],"_pc4",sep="")
# pc3<-pc3[which(pc3$MarkerName %in% pc5$MarkerName),]
# colnames(pc3)[2:ncol(pc3)]<-paste(colnames(pc3)[2:ncol(pc3)],"_pc3",sep="")
# pc2<-pc2[which(pc2$MarkerName %in% pc5$MarkerName),]
# colnames(pc2)[2:ncol(pc2)]<-paste(colnames(pc2)[2:ncol(pc2)],"_pc2",sep="")
# pc1<-pc1[which(pc1$MarkerName %in% pc5$MarkerName),]
# colnames(pc1)[2:ncol(pc1)]<-paste(colnames(pc1)[2:ncol(pc1)],"_pc1",sep="")

# pc<-merge(pc1,pc2,by="MarkerName")
# pc<-merge(pc,pc3,by="MarkerName")
# pc<-merge(pc,pc4,by="MarkerName")
# pc<-merge(pc,pc5,by="MarkerName")


# p0<-ggplot(pc, aes(x=-log10(pc$'P-value_association_pc2'), y=-log10(pc$'P-value_association_pc1'))) +
#   geom_point() + geom_abline(intercept = 0, slope = 1,color="blue")

# p1<-ggplot(pc, aes(x=-log10(pc$'P-value_association_pc3'), y=-log10(pc$'P-value_association_pc1'))) +
#   geom_point() + geom_abline(intercept = 0, slope = 1,color="blue")

# p2<-ggplot(pc, aes(x=-log10(pc$'P-value_association_pc4'), y=-log10(pc$'P-value_association_pc1'))) +
#   geom_point() + geom_abline(intercept = 0, slope = 1,color="blue")

# p3<-ggplot(pc, aes(x=-log10(pc$'P-value_association_pc5'), y=-log10(pc$'P-value_association_pc1'))) +
#   geom_point() + geom_abline(intercept = 0, slope = 1,color="blue")

# p<-ggarrange(p0,p1,p2,p3,ncol=4)


# pdf("~/git/IIBDGC_GWAS/mrmega_pvalue_association.pdf",width = 20, height = 5)
#  annotate_figure(p, top = text_grob("P-value association"))
# dev.off()

# p0<-ggplot(pc, aes(x=-log10(pc$'P-value_ancestry_het_pc2'), y=-log10(pc$'P-value_ancestry_het_pc1'))) +
#   geom_point() + geom_abline(intercept = 0, slope = 1,color="blue")

# p1<-ggplot(pc, aes(x=-log10(pc$'P-value_ancestry_het_pc3'), y=-log10(pc$'P-value_ancestry_het_pc1'))) +
#   geom_point() + geom_abline(intercept = 0, slope = 1,color="blue")

# p2<-ggplot(pc, aes(x=-log10(pc$'P-value_ancestry_het_pc4'), y=-log10(pc$'P-value_ancestry_het_pc1'))) +
#   geom_point() + geom_abline(intercept = 0, slope = 1,color="blue")

# p3<-ggplot(pc, aes(x=-log10(pc$'P-value_ancestry_het_pc5'), y=-log10(pc$'P-value_ancestry_het_pc1'))) +
#   geom_point() + geom_abline(intercept = 0, slope = 1,color="blue")

# p<-ggarrange(p0,p1,p2,p3,ncol=4)

# pdf("~/git/IIBDGC_GWAS/mrmega_pvalue_ancestry_het.pdf",width = 15, height = 5)
#  annotate_figure(p, top = text_grob("P-value ancestry Het"))
# dev.off()

# p0<-ggplot(pc, aes(x=-log10(pc$'P-value_residual_het_pc2'), y=-log10(pc$'P-value_residual_het_pc1'))) +
#   geom_point() + geom_abline(intercept = 0, slope = 1,color="blue")

# p1<-ggplot(pc, aes(x=-log10(pc$'P-value_residual_het_pc3'), y=-log10(pc$'P-value_residual_het_pc1'))) +
#   geom_point() + geom_abline(intercept = 0, slope = 1,color="blue")

# p2<-ggplot(pc, aes(x=-log10(pc$'P-value_residual_het_pc4'), y=-log10(pc$'P-value_residual_het_pc1'))) +
#   geom_point() + geom_abline(intercept = 0, slope = 1,color="blue")

# p3<-ggplot(pc, aes(x=-log10(pc$'P-value_residual_het_pc5'), y=-log10(pc$'P-value_residual_het_pc1'))) +
#   geom_point() + geom_abline(intercept = 0, slope = 1,color="blue")

# p<-ggarrange(p0,p1,p2,p3,ncol=4)

# pdf("~/git/IIBDGC_GWAS/mrmega_pvalue_residual_het.pdf",width = 15, height = 5)
#  annotate_figure(p, top = text_grob("P-value residual Het"))
# dev.off()


# ## mds all chr, plus
# vec<-c("0","-0.017557","0.00369821","-0.00216237",
# "1","-0.018128","0.0039589","-0.0011246",
# "2","-0.0178759","0.00596411","0.00487938",
# "3","-0.0185202","0.00568787","0.00804003",
# "4","0.166458","0.0201149","-0.00358016",
# "5","-0.0059465","0.00924409","0.039498",
# "6","-0.0155228","0.00140931","-0.00419584",
# "7","-0.0170316","0.00353529","-0.00227943",
# "8","-0.0166462","0.00174916","-0.00291361",
# "9","-0.0145416","-0.00122709","-0.012414",
# "10","-0.0148536","-0.000347214","-0.00940286",
# "11","-0.017846","0.00372435","-0.000889448",
# "12","0.0424821","-0.0632909","0.0050225",
# "13","-0.01653","0.00170933","-0.0175688",
# "14","-0.0179404","0.00406966","-0.000908843")


# mds<-as.data.frame(matrix(vec,ncol=4,byrow=T))
# colnames(mds)<-c("n","pc1","pc2","pc3")
# mds$studies<-c("affymetrix500_eur","affymetrix6_eur","danish_eur","decode_eur",
# "iibdgc_eas","finngen_eur","gsa_eur","humancoreexome_eur","humanomni1_eur","illumina370_eur","illuminaexome_eur","interval_ibdbioresource_eur","interval_ibdbioresource_sas","quad610_eur","ukbb_eur")

# mds$pc1<-as.numeric(as.character(mds$pc1))
# mds$pc2<-as.numeric(as.character(mds$pc2))
# mds$pc3<-as.numeric(as.character(mds$pc3))

# pdf("~/git/IIBDGC_GWAS/plots/mrmega_mds_pc1_pc2.pdf",width = 8, height = 5)
#  ggplot(mds,aes(x=pc1,y=pc2,color=studies)) + geom_point()
# dev.off()

# pdf("~/git/IIBDGC_GWAS/plots/mrmega_mds_pc1_pc3.pdf",width = 8, height = 5)
#  ggplot(mds,aes(x=pc1,y=pc3,color=studies)) + geom_point()
# dev.off()

# based on results above, two pcs are enough to separate studies by broad continental ancestry, test now per phenotype

# MEM=4000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1


library(data.table)
library(ggplot2)
library(ggpubr)

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

# cd<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/mrmega/output_files/cd/mrmega_eur_eas_sas_tier_2_cd_allchr_pc2.result",sep=""))
# uc<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/mrmega/output_files/cd/mrmega_eur_eas_sas_tier_2_cd_allchr_pc2.result",sep=""))
# ibd<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/mrmega/output_files/cd/mrmega_eur_eas_sas_tier_2_cd_allchr_pc2.result",sep=""))


## mds all chr, plus
vec_cd<-c("0","-0.017557","0.00369821",
"1","-0.018128","0.0039589",
"2","-0.0178759","0.00596411",
"3","-0.0185202","0.00568787",
"4","0.166458","0.0201149",
"5","-0.0059465","0.00924409",
"6","-0.0155228","0.00140931",
"7","-0.0170316","0.00353529",
"8","-0.0166462","0.00174916",
"9","-0.0145416","-0.00122709",
"10","-0.0148536","-0.000347214",
"11","-0.017846","0.00372435",
"12","0.0424821","-0.0632909",
"13","-0.01653","0.00170933",
"14","-0.0179404","0.00406966")

vec_ibd<-c("0","-0.0198379","0.00495123",
"1","-0.0200489","0.00539814",
"2","-0.0203912","0.00686363",
"3","-0.0209915","0.00668613",
"4","0.160467"," 0.0415474",
"5","-0.00890562","0.0107364",
"6","-0.0174339","0.00289503",
"7","-0.0192345","0.00465156",
"8","-0.018761","0.00356705",
"9","-0.0195854","0.00480767",
"10","-0.0168369","0.00337431",
"11","-0.016591","0.00117784",
"12","-0.0202293","0.00525734",
"13","0.0464062","-0.0541089",
"14","-0.0186626","0.00389909",
"15","-0.0203002","0.005383",
"16","0.0509369","-0.0570869")

vec_uc<-c("0","-0.0165979","0.00347956",
"1","-0.0165988","0.00352493",
"2","-0.016863","0.00535587",
"3","-0.0174592","0.00504758",
"4","0.16694","0.0200669",
"5","-0.00521471","0.00849603",
"6","-0.0144005","0.00112054",
"7","-0.0159037","0.00299708",
"8","-0.01539"," 0.0013554",
"9","-0.0161813","0.00282403",
"10","-0.0128778","0.00182792",
"11","-0.0135593","-0.000935049",
"12","-0.0167677","0.00322311",
"13","0.0430927","-0.0634233",
"14","-0.0153809","0.00152972",
"15","-0.0168378","0.00350966")

pheno<-c("cd","uc","ibd")

for (ph in pheno) {

  vec<-get(paste("vec_",ph,sep=""))

  if (ph=="cd") {
    mds<-as.data.frame(matrix(vec,ncol=3,byrow=T))
    mds$V2<-as.numeric(as.character(mds$V2))
    mds$V3<-as.numeric(as.character(mds$V3))

    colnames(mds)<-c("n",paste("pc",seq(1:2),"_",ph,sep=""))
    mds$studies<-c("affymetrix500_eur","affymetrix6_eur","danish_eur","decode_eur","iibdgc_eas","finngen_eur","gsa_eur","humancoreexome_eur","humanomni1_eur","illumina370_eur","illuminaexome_eur","interval_ibdbioresource_eur","interval_ibdbioresource_sas","quad610_eur","ukbb_eur")

  } else {

    vec<-get(paste("vec_",ph,sep=""))

    tmp<-as.data.frame(matrix(vec,ncol=3,byrow=T))
    tmp$V2<-as.numeric(as.character(tmp$V2))
    tmp$V3<-as.numeric(as.character(tmp$V3))

    colnames(tmp)<-c("n",paste("pc",seq(1:2),"_",ph,sep=""))

    if(ph=="uc") {
      tmp$studies<-c("affymetrix500_eur","affymetrix6_eur","danish_eur","decode_eur","iibdgc_eas","finngen_eur","gsa_eur","humancoreexome_eur","humanomni1_eur","humanomniexpress_eur","illumina370_eur","illuminaexome_eur","interval_ibdbioresource_eur","interval_ibdbioresource_sas","quad610_eur","ukbb_eur")
    } else {
      tmp$studies<-c("affymetrix500_eur","affymetrix6_eur","danish_eur","decode_eur","iibdgc_eas","finngen_eur","gsa_eur","humancoreexome_eur","humanomni1_eur","humanomniexpress_eur","illumina370_eur","illuminaexome_eur","interval_ibdbioresource_eur","interval_ibdbioresource_sas","quad610_eur","ukbb_eur","ukbb_sas")
    }
    mds<-merge(mds,tmp[,c(2:ncol(tmp))],by="studies",all=T)
  }
  

}

# add global ancestry flag:
mds$global_ancestry<-NA
mds$global_ancestry[grep("_eur",mds$studies)]<-"EUR_nf"
mds$global_ancestry[grep("finngen",mds$studies)]<-"EUR_finish"
mds$global_ancestry[grep("_sas",mds$studies)]<-"SAS"
mds$global_ancestry[grep("_eas",mds$studies)]<-"EAS"

p1<-ggplot(mds,aes(x=pc1_cd,y=pc2_cd,color=global_ancestry)) + geom_point() + ggtitle("CD")
p2<-ggplot(mds,aes(x=pc1_uc,y=pc2_uc,color=global_ancestry)) + geom_point() + ggtitle("UC")
p3<-ggplot(mds,aes(x=pc1_ibd,y=pc2_ibd,color=global_ancestry)) + geom_point() + ggtitle("IBD")

pdf("~/git/IIBDGC_GWAS/plots/mrmega_mds_pc1_pc2_CD_UC_IBD.pdf",width = 15, height = 5,onefile=FALSE)
 ggarrange(p1,p2,p3,ncol=3,common.legend=T,legend=c("top"))
dev.off()



q("no")
