# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# compare different population frequency sources:

path_gwas=/path/to/ibdgwas/IIBDGC/

studies=(australia_omniexome all_hce pittsburgh_gsa italy_gsa niddk_broad_gsa niddk_feinstein_gsa belgium_louis_gsa spain_gsa)

for i in ${studies[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01_nomissidiscrep_perstudy_nodupsample_perstudy \
--allow-no-sex \
--extract ${path_gwas}pre_imputation/QC/relatedness/list_common_variants_among_cohorts_iibdgc \
--make-bed \
--out ${path_gwas}pre_imputation/QC/relatedness/${i}_subset
done

# ## TOPMED
# awk 'NR==FNR{vals[\$2];next} (\$1) in vals' ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim \
# ${path_gwas}resources/gnomad/gnomad_freq_edited > ${path_gwas}pre_imputation/QC/${i}/${i}_TOPMed_variants
# 
# # ALFA
# awk 'NR==FNR{vals[$2];next} ($6) in vals' ${path_gwas}/pre_imputation/QC/${i}/niddk_broad_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim \
# ${path_gwas}resources/alfa/freq_all_continental_populations > ${path_gwas}pre_imputation/QC/${i}/${i}_alfa_variants
# 
# # GNOMAD:
# awk 'NR==FNR{vals[\$2];next} (\$1) in vals' ${path_gwas}/pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim \
# ${path_gwas}resources/gnomad/gnomad_freq_edited > ${path_gwas}pre_imputation/QC/${i}/${i}_gnomad_variants
# 

awk 'NR==FNR{vals[$2];next} ($6) in vals' ${path_gwas}/pre_imputation/QC/spain_gsa/spain_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed.bim \
/path/to/project > ${path_gwas}pre_imputation/QC/spain_gsa/spain_gsa_alfa_variants
  
# after evaluating A/T C/G alleles:

/software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"


cbPalette <- c("#999999", "#E69F00", "#56B4E9")
cohorts<-c("australia_omniexome","all_hce","pittsburgh_gsa","italy_gsa","niddk_broad_gsa","niddk_feinstein_gsa","belgium_louis_gsa","spain_gsa")


# for (xx in 1:length(cohorts)) {

for (xx in 8) {
  
  file_anc<-paste(path,"pre_imputation/QC/",cohorts[xx],"/",cohorts[xx],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_eur_ctr.frq.counts",sep="")
  vfreq<-read.table(file_anc,head=T)
  vfreq[,paste("eur","freq_alt",sep="_")]<-vfreq$C1/(vfreq$C1+vfreq$C2)
  
  print(cohorts[xx])
  print(nrow(vfreq))

  # ALFA
  alfa<-read.table(paste(path,"pre_imputation/QC/",cohorts[xx],"/",cohorts[xx],"_alfa_variants",sep=""),head=F)
  colnames(alfa)<-c("CHROM","POS","ID","REF","ALT","SNP","CHR_POS","EUR_freq_Alt","AFR_freq_Alt","EAS_freq_Alt","SAS_freq_Alt","AMR_freq_Alt")
  alfa<-alfa[which(alfa$SNP %in% vfreq$SNP),]
  print("ALFA")
  print(nrow(alfa))
  print(nrow(alfa)/nrow(vfreq))
  
  # Gnomad
  gnomad<-read.table(paste(path,"pre_imputation/QC/",cohorts[xx],"/",cohorts[xx],"_gnomad_variants",sep=""),head=F)
  colnames(gnomad)<-c("SNP","CHROM","POS","REF","ALT","AF","AF_nfe","AF_afr","AF_amr","AF_eas","AF_sas","AF_asj")
  gnomad$AF_nfe<-as.numeric(as.character(gnomad$AF_nfe))
  gnomad<-gnomad[which(gnomad$SNP %in% vfreq$SNP),]
  print("Gnomad")
  print(nrow(gnomad))
  print(nrow(gnomad)/nrow(vfreq))
  
  # TOPMEd
  tm<-read.table(paste(path,"pre_imputation/QC/",cohorts[xx],"/",cohorts[xx],"_TOPMed_variants",sep=""),head=F)
  tm<-tm[,c(9,8)]
  colnames(tm)<-c("SNP","topmed_freq_alt")
  tm<-tm[which(tm$SNP %in% vfreq$SNP),]
  print("TOPMEd")
  print(nrow(tm))
  print(nrow(tm)/nrow(vfreq))
  

  all<-merge(vfreq[,c(2,8)],alfa,by="SNP",sort=F)
  all<-merge(all,tm,by="SNP",all.x=T,sort=F)
  all<-merge(all,gnomad[,c("SNP","AF","AF_nfe","AF_afr","AF_amr","AF_eas","AF_sas","AF_asj")],by="SNP",all.x=T,sort=F)


  p1<-ggplot(all, aes(y=eur_freq_alt, x=EUR_freq_Alt)) +
    geom_point() + ylab(paste(cohorts[xx],"FRQ Alt")) + xlab(paste("ALFA EUR FRQ Alt\nN shared variants ",nrow(alfa),sep="")) + scale_colour_manual(values=cbPalette) +
    ggtitle(paste(cohorts[xx],"EUR",sep=" "))

  p2<-ggplot(all, aes(y=eur_freq_alt, x=topmed_freq_alt)) +
    geom_point() + ylab(paste(cohorts[xx],"FRQ Alt")) + xlab(paste("TOPMed FRQ Alt\nN shared variants ",nrow(tm),sep="")) + scale_colour_manual(values=cbPalette) +
    ggtitle(paste(cohorts[xx],"EUR",sep=" "))

  p3<-ggplot(all, aes(y=eur_freq_alt, x=AF_nfe)) +
    geom_point() + ylab(paste(cohorts[xx],"FRQ Alt")) + xlab(paste("Gnomad EUR non Finish FRQ Alt\nN shared variants ",nrow(gnomad),sep="")) + scale_colour_manual(values=cbPalette) +
    ggtitle(paste(cohorts[xx],"EUR",sep=" "))

  p<-ggarrange(p1,p2,p3,ncol=3,nrow=2,common.legend=T)

  pdf(paste(path,"pre_imputation/QC/",cohorts[xx],"/",cohorts[xx],"_plot_maf_list_variants_toexclude_maf_differs_ALFA_TOPMed_Gnomad_eur.pdf",sep=""),width = 15, height = 10)
  print(p)
  dev.off()
  system(paste("cp ",path,"pre_imputation/QC/",cohorts[xx],"/",cohorts[xx],"_plot_maf_list_variants_toexclude_maf_differs_ALFA_TOPMed_Gnomad_eur.pdf ~/tmp_plots/",sep=""))

}




# [1] "australia_omniexome"
# [1] 727259
# [1] "ALFA"
# [1] 726884
# [1] 0.9994844
# [1] "Gnomad"
# [1] 727204
# [1] 0.9999244
# [1] "TOPMEd"
# [1] 727259
# [1] 1


# [1] "all_hce"
# [1] 407285
# [1] "ALFA"
# [1] 407159
# [1] 0.9996906
# [1] "Gnomad"
# [1] 404867
# [1] 0.9940631
# [1] "TOPMEd"
# [1] 407285
# [1] 1


# [1] "pittsburgh_gsa"
# [1] 873827
# [1] "ALFA"
# [1] 873798
# [1] 0.9999668
# [1] "Gnomad"
# [1] 873746
# [1] 0.9999073
# [1] "TOPMEd"
# [1] 873827
# [1] 1

# [1] "italy_gsa"
# [1] 544961
# [1] "ALFA"
# [1] 544353
# [1] 0.9988843
# [1] "Gnomad"
# [1] 544886
# [1] 0.9998624
# [1] "TOPMEd"
# [1] 544961
# [1] 1

# [1] "niddk_broad_gsa"
# [1] 570569
# [1] "ALFA"
# [1] 570096
# [1] 0.999171
# [1] "Gnomad"
# [1] 570362
# [1] 0.9996372
# [1] "TOPMEd"
# [1] 570569
# [1] 1

# [1] "niddk_feinstein_gsa"
# [1] 581841
# [1] "ALFA"
# [1] 581380
# [1] 0.9992077
# [1] "Gnomad"
# [1] 581598
# [1] 0.9995824
# [1] "TOPMEd"
# [1] 581841
# [1] 1

# [1] "belgium_louis_gsa"
# [1] 543463
# [1] "ALFA"
# [1] 543019
# [1] 0.999183
# [1] "Gnomad"
# [1] 543276
# [1] 0.9996559
# [1] "TOPMEd"
# [1] 543463
# [1] 1

# [1] "spain_gsa"
# [1] 568656
# [1] "ALFA"
# [1] 568624
# [1] 0.9999437
# [1] "Gnomad"
# [1] 568645
# [1] 0.9999807
# [1] "TOPMEd"
# [1] 568656
# [1] 1

#######################################################################################################################################################################################
#### variants excluded by Kyle:

/software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

## BROAD:

br<-read.table("/path/to/ibdgwas/IIBDGC/from_kyle/Broad_flagged67262.txt",head=F)
dim(br)
# [1] 67262     1



bim<-read.table(paste(path,"pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_hg19.bim",sep=""),head=F)
dim(br[which(br$V1 %in% bim$V2),,drop=F])
# 67262 53401     2
dim(br[which(!br$V1 %in% bim$V2),,drop=F])
# [1] 0   2


bim<-read.table(paste(path,"pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_hg19_noind_posstr.bim",sep=""),sep="\t",head=F)
ids<-read.table(paste(path,"pre_imputation/QC/niddk_broad_gsa/list_variants_niddk_broad_gsa_hg19_posstrandaligned",sep=""),sep=" ",head=F,skip=33)
colnames(ids)[2]<-"ids_hg19"
bim.1<-cbind(bim,ids[,"ids_hg19",drop=F])

dim(br[which(br$V1 %in% bim.1$V2),,drop=F])
# [1] 56124     1
dim(br[which(!br$V1 %in% bim.1$V2),,drop=F])
# [1] 11138     1                            # already excluded in steps before, 

br<-merge(br,bim.1[,c("V2","ids_hg19")],by.x="V1",by.y="V2",all.x=T,sort=F)

exclusion<-read.table(paste(path,"pre_imputation/QC/niddk_broad_gsa/list_indel_var_exclude",sep=""),head=F)
br$excluded<-NA
br$excluded[which(br$V1 %in% exclusion$V1)]<-1
table(br$excluded)
# 1 
# 11138

bim<-read.table(paste(path,"pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned.bim",sep=""),sep="\t",head=F)
ids<-read.table(paste(path,"pre_imputation/QC/niddk_broad_gsa/list_variants_niddk_broad_gsa_hg38_posstrandaligned",sep=""),sep=" ",head=F,skip=34)
colnames(ids)[2]<-"ids_b38"
bim.1<-cbind(bim,ids[,"ids_b38",drop=F])

dim(br[which(br$ids_hg19 %in% bim.1$V2),,drop=F])
# [1] 28366     1
dim(br[which(!br$ids_hg19 %in% bim.1$V2),,drop=F])
# [1] 38896     1                            # already excluded in steps before,

br<-merge(br,bim.1[,c("V2","ids_b38")],by.x="ids_hg19",by.y="V2",all.x=T,sort=F)
dim(br[which(!is.na(br$ids_b38)),])
# [1] 28366     4

bim_final<-read.table(paste(path,"pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2.bim",sep=""),sep="")
bim_final$to_exclude<-0
bim_final$to_exclude[which(bim_final$V2 %in% br$ids_b38)]<-1
#      0      1 
# 544875  25694 

67262-25694
# [1] 41568

cbPalette <- c("#999999", "#E69F00", "#56B4E9")
cohorts<-c("niddk_broad_gsa")

for (xx in 1) {
  
  file_anc<-paste(path,"pre_imputation/QC/",cohorts[xx],"/",cohorts[xx],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_eur_ctr.frq.counts",sep="")
  vfreq<-read.table(file_anc,head=T)
  vfreq[,paste("eur","freq_alt",sep="_")]<-vfreq$C1/(vfreq$C1+vfreq$C2)
  
  print(cohorts[xx])
  print(nrow(vfreq))
  
  # source
  
  gnomad<-read.table(paste(path,"pre_imputation/QC/",cohorts[xx],"/",cohorts[xx],"_gnomad_variants",sep=""),head=F)
  colnames(gnomad)<-c("SNP","CHROM","POS","REF","ALT","AF","AF_nfe","AF_afr","AF_amr","AF_eas","AF_sas","AF_asj")
  gnomad$AF_nfe<-as.numeric(as.character(gnomad$AF_nfe))
  gnomad<-gnomad[which(gnomad$SNP %in% vfreq$SNP),]
  print("Gnomad")
  print(nrow(gnomad))
  print(nrow(gnomad)/nrow(vfreq))
  
  all<-merge(vfreq[,c(2,8)],gnomad[,c("SNP","AF","AF_nfe","AF_afr","AF_amr","AF_eas","AF_sas","AF_asj")],by="SNP",sort=F)
  all$toexclude<-"0"
  all$SNP<-as.character(all$SNP)
  vec_to_exclude<-as.character(bim_final$V2[which(bim_final$to_exclude=="1")])
  all$toexclude[which(all$SNP %in% vec_to_exclude)]<-"1"
  
  all$value<-((all$eur_freq_alt-all$AF_nfe)^2)/((all$eur_freq_alt+all$AF_nfe)*(2-all$eur_freq_alt-all$AF_nfe))
  value<-0.025
  all$col[which(all$value<=value)]<-paste("<=",value,sep="")
  all$col[which(all$value>value)]<-paste(">",value,sep="")
  all$col[which(is.na(all$value))]<-"monomorphic"
  
  print(table(all$col,all$toexclude))
  #                  0      1
  # <=0.025     531134  24828
  # >0.025         765    189
  # monomorphic  12785    661
  
  tmp1<-all[which(all$col=="<=0.025"),]
  tmp2<-all[which(all$col==">0.025"),]
  
  p1<-ggplot(tmp1,aes(y=eur_freq_alt, x=AF_nfe)) +
    geom_point(aes(colour = toexclude)) + ylab(paste(cohorts[xx],"FRQ Alt")) + xlab(paste("Gnomad EUR non Finish FRQ Alt\nN shared variants ",nrow(gnomad),sep="")) + scale_colour_manual(values=cbPalette) +
    ggtitle(paste(cohorts[xx],"Keep based on allele freq",sep=" "))
  p2<-ggplot(tmp2,aes(y=eur_freq_alt, x=AF_nfe)) +
    geom_point(aes(colour = toexclude)) + ylab(paste(cohorts[xx],"FRQ Alt")) + xlab(paste("Gnomad EUR non Finish FRQ Alt\nN shared variants ",nrow(gnomad),sep="")) + scale_colour_manual(values=cbPalette) +
    ggtitle(paste(cohorts[xx],"Exclude based on allele freq",sep=" "))
  
  
  p<-ggarrange(p1,p2,ncol=2,nrow=1,common.legend=T)
  
  pdf(paste(path,"pre_imputation/QC/",cohorts[xx],"/",cohorts[xx],"_plot_maf_list_variants_toexclude_maf_differs_ALFA_TOPMed_Gnomad_eur_exclusionKyleTalin.pdf",sep=""),width = 10, height = 5)
  print(p)
  dev.off()
  system(paste("cp ",path,"pre_imputation/QC/",cohorts[xx],"/",cohorts[xx],"_plot_maf_list_variants_toexclude_maf_differs_ALFA_TOPMed_Gnomad_eur_exclusionKyleTalin.pdf ~/tmp_plots/",sep=""))
  
  

  
  all$SNP2<-sub(":","_",all$SNP)
  all$SNP3<-sub("_[A-Z]*_[A-Z]*$","",all$SNP)
  
  tmp3<-all[which(all$SNP2 %in% r2$SNP2),]
  dim(tmp3)
  # [1] 66 13

  tmp3<-all[which(all$SNP3 %in% r2$SNP3),]
  dim(tmp3)
  # [1] 66 14
  
}

# VARIANTS LOW R2
r2<-read.csv("/path/to/ibdgwas/IIBDGC/from_kyle/LowER2snpsNIDDK.csv",head=T)
r2$SNP2<-gsub(":","_",r2$SNP)
r2$SNP3<-gsub(":[A-Z]*:[A-Z]*$","",r2$SNP)
r2$chr<-gsub(":.*","",r2$SNP3)
r2$position<-gsub(".*:","",r2$SNP3)
r2$chr<-as.numeric(r2$chr)
r2$position<-as.numeric(r2$position)
r2<-r2[order(r2$chr,r2$position,decreasing=F),]


bim_final$SNP2<-sub(":","_",bim_final$V2)
bim_final$SNP3<-sub("_[A-Z]*_[A-Z]*$","",bim_final$V2)

dim(r2[which(r2$SNP2 %in% bim_final$SNP2),])
# [1] 66  8

# file post liftover to b38
bim$SNP2<-sub(":","_",bim$V2)
bim$SNP3<-sub("_[A-Z]*_[A-Z]*$","",bim$V2)
dim(r2[which(r2$SNP2 %in% bim$SNP2),])
# [1] 4150    8


topmed<-read.table(paste(path,"pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_TOPMed_variants_ids",sep=""))
topmed$SNP2<-sub(":","_",topmed$V1)
dim(r2[which(r2$SNP2 %in% topmed$SNP2),])
# [1] 66  8



