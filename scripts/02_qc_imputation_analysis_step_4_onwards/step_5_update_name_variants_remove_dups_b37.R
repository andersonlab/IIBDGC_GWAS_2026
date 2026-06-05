# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#

# step_5_update_name_variants_remov_b37

################################################################################
# 5.- UPDATE NAME VARIANTS TO CHR:POSITION_REF_ALT; REMOVE DUPLICATED VARIANTS #
################################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
  
studies=(australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)
  
for i in ${studies[@]}
do
zcat ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_posstrandaligned_2.vcf.gz | cut -f '1-5' | awk '{print $3,$1":"$2"_"$4"_"$5}' \
> ${path_gwas}pre_imputation/QC/${i}/list_variants_${i}_hg19_posstrandaligned
done



##############################
### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

# niddk_old_gwas all_hce - add in later step
path<-"/path/to/ibdgwas/IIBDGC/"
studies2<-c("bernstein_gsa","farkkila_gsa","franchimont_gsa","franke_gsa",
           "hyams_protect_gsa","lewis_sparc_gsa","mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa",
           "niddk_cho_gsa","niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa","xavier_share_gsa",
           "australia_omniexome",
           "basque_gsa","belgium_franchimont_gsa","belgium_louis_gsa","belgium_vermeire_gsa","ccfa_gsa","cedars_gsa",
           "italy_gsa","kiel_austria_sibdcs_gsa","lithuania_gsa","mccauley_gsa","netherlands_gsa",
           "niddk_broad_gsa","niddk_feinstein_gsa","prism_nfe_gsa","slovenia_gsa","sweden_gsa",
           "helmsley_prism_gsa","helmsley_xavier_prism_gsa","prism_nfe_gwas","finland_illugwas",
           "cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas",
           "chop_old_gwas","german_affy6_old_gwas","norway_affy6_old_gwas",
           "belgium_inf1_old_gwas","belgium_inf2_old_gwas","swedish_uc_old_gwas",
           "pittsburgh_gsa","spain_gsa","gwas1","gwas2")

for (i in 52:length(studies)) {
  
  print(studies[i])
  bim<-read.table(paste(path,"pre_imputation/QC/",studies[i],"/",studies[i],"_hg19_noind_posstr.bim",sep=""),sep="\t",head=F)
  
  ids<-read.table(paste(path,"pre_imputation/QC/",studies[i],"/list_variants_",studies[i],"_hg19_posstrandaligned",sep=""),sep=" ",head=F
                  ,skip=30)
  ids<-ids[which(ids$V2!=""),]
  
  # if(studies[i] %in% non_gsa) {
  #   ids<-read.table(paste(path,"pre_imputation/QC/",studies[i],"/list_variants_",studies[i],"_hg19_posstrandaligned",sep=""),sep=" ",head=F
  #                   ,skip=32)
  #   ids<-ids[which(ids$V2!=""),]
  # } else {
  #   ids<-read.table(paste(path,"pre_imputation/QC/",studies[i],"/list_variants_",studies[i],"_hg19_posstrandaligned",sep=""),sep=" ",head=F
  #                   ,skip=33)
  #   ids<-ids[which(ids$V2!=""),]
  # }
  
  varmiss<-read.table(paste(path,"pre_imputation/QC/",studies[i],"/",studies[i],"_hg19_noind_posstr.lmiss",sep=""),head=T)
  colnames(ids)[2]<-"ids"
  bim$V2<-as.character(bim$V2)
  ids$V1<-as.character(ids$V1)
  print(table(bim$V2==ids$V1))
  
  #the major allele is set to A2 by default by Plink, keep ids with real ref/alt as in vcf using the ids file
  bim.1<-cbind(bim,ids[,"ids",drop=F])
  
  bim.1<-merge(bim.1,varmiss[,c("SNP","F_MISS")],by.x="V2",by.y="SNP",sort=F)
  print(table(bim.1$V2==bim$V2))
  # 
  bim.1$ids<-as.character(bim.1$ids)
  
  # identify duplicated variants (same chr position ref and alt)
  dups<-bim.1[which(duplicated(bim.1$ids)),"ids"]
  
  print(paste("N variants duplicated:",length(dups)))
  
  # bim.1[which(duplicated(bim.1$V2)),]
  
  if (length(dups)>0) {
    for (ii in 1:length(dups)){
      
      tmp<-bim.1[which(bim.1$ids %in% dups[ii]),]
      keep<-tmp[which(tmp$F_MISS==min(tmp$F_MISS)),]
      
      if (nrow(keep)>1){
        keep<-keep[1,]
      }
      
      exclude<-tmp[which(!tmp$V2 %in% keep$V2),]
      bim.1$ids[which(bim.1$V2 %in% exclude$V2)]<-paste(bim.1$ids[which(bim.1$V2 %in% exclude$V2)],"_rm",sep="")
      
    }
  }
 
  print(paste("N IDs more than 2 probes:",nrow(bim.1[which(duplicated(bim.1$ids)),])))
  
  duplicated_variants<-bim.1[grep("_rm",bim.1$ids),"ids"]
  print(paste("N duplicated variants to exclude:",length(duplicated_variants)))
  
  write.table(duplicated_variants,paste(path,"pre_imputation/QC/",studies[i],"/list_duplicated_var_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  
  print(table(bim.1$V2==bim$V2))
  
  # print(table(bim.1$V1))
  #
  # NOTE chrx pseudoautosomal now in chr23, recode this
  
  # https://www.ncbi.nlm.nih.gov/grc/human
  # Name Chr Start Stop
  # PAR#1	X	60,001	2,699,520
  # PAR#2	X	154,931,044	155,260,560
  # PAR#1	Y	10,001	2,649,520
  # PAR#2	Y	59,034,050	59,363,566
  
  print(dim(bim.1[which( (bim.1$V1==23) & ((bim.1$V4>=60001 & bim$V4<=2699520) | (bim.1$V4>=154931044 & bim$V4<=155260560)) ),]))
  print(dim(bim.1[which( (bim.1$V1==23) & ((bim.1$V4>=60001 & bim$V4<=2699520) | (bim.1$V4>=154931044 & bim$V4<=155260560)) & (!bim.1$ids %in% duplicated_variants) ),]))
  print(dim(bim.1[which( (bim.1$V1==24) & ((bim.1$V4>=10001 & bim$V4<=2649520) | (bim.1$V4>=59034050 & bim$V4<=59363566)) ),]))
  
  write.table(bim.1[,c(2,7,3:6)],paste(path,"pre_imputation/QC/",studies[i],"/",studies[i],"_hg19_noind_posstr_edited.bim",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  
  rm(list=ls()[!ls() %in% c("i","studies","path","non_gsa")])
}

q("no")


##############################

# [1] "bernstein_gsa"
# 
# TRUE 
# 565758 
# 
# TRUE 
# 565758 
# [1] "N variants duplicated: 290"
# [1] "N IDs more than 2 probes: 4"
# [1] "N duplicated variants to exclude: 291"
# 
# TRUE 
# 565758 
# [1] 452   8
# [1] 451   8
# [1] 0 8


# [1] "farkkila_gsa"
# 
# TRUE 
# 528311 
# 
# TRUE 
# 528311 
# [1] "N variants duplicated: 271"
# [1] "N IDs more than 2 probes: 3"
# [1] "N duplicated variants to exclude: 272"
# 
# TRUE 
# 528311 
# [1] 403   8
# [1] 402   8
# [1] 0 8


# [1] "franchimont_gsa"
# 
# TRUE 
# 671132 
# 
# TRUE 
# 671132 
# [1] "N variants duplicated: 451"
# [1] "N IDs more than 2 probes: 5"
# [1] "N duplicated variants to exclude: 452"
# 
# TRUE 
# 671132 
# [1] 515   8
# [1] 514   8
# [1] 0 8


# [1] "franke_gsa"
# 
# TRUE 
# 571118 
# 
# TRUE 
# 571118 
# [1] "N variants duplicated: 290"
# [1] "N IDs more than 2 probes: 4"
# [1] "N duplicated variants to exclude: 291"
# 
# TRUE 
# 571118 
# [1] 447   8
# [1] 446   8
# [1] 0 8


# [1] "hyams_protect_gsa"
# 
# TRUE 
# 580928 
# 
# TRUE 
# 580928 
# [1] "N variants duplicated: 297"
# [1] "N IDs more than 2 probes: 4"
# [1] "N duplicated variants to exclude: 298"
# 
# TRUE 
# 580928 
# [1] 483   8
# [1] 482   8
# [1] 0 8


# [1] "lewis_sparc_gsa"
# 
# TRUE 
# 619895 
# 
# TRUE 
# 619895 
# [1] "N variants duplicated: 331"
# [1] "N IDs more than 2 probes: 5"
# [1] "N duplicated variants to exclude: 332"
# 
# TRUE 
# 619895 
# [1] 493   8
# [1] 492   8
# [1] 0 8


# [1] "mccauley_new_gsa"
# 
# TRUE 
# 591351 
# 
# TRUE 
# 591351 
# [1] "N variants duplicated: 310"
# [1] "N IDs more than 2 probes: 4"
# [1] "N duplicated variants to exclude: 311"
# 
# TRUE 
# 591351 
# [1] 484   8
# [1] 483   8
# [1] 0 8

# [1] "mcgovern_gsa"
# 
# TRUE 
# 652046 
# 
# TRUE 
# 652046 
# [1] "N variants duplicated: 370"
# [1] "N IDs more than 2 probes: 5"
# [1] "N duplicated variants to exclude: 371"
# 
# TRUE 
# 652046 
# [1] 500   8
# [1] 499   8
# [1] 0 8

# [1] "moayyedi_imagine_gsa"
# 
# TRUE 
# 615521 
# 
# TRUE 
# 615521 
# [1] "N variants duplicated: 323"
# [1] "N IDs more than 2 probes: 5"
# [1] "N duplicated variants to exclude: 324"
# 
# TRUE 
# 615521 
# [1] 490   8
# [1] 489   8
# [1] 0 8

# [1] "newberry_share_gsa"
# 
# TRUE 
# 592708 
# 
# TRUE 
# 592708 
# [1] "N variants duplicated: 312"
# [1] "N IDs more than 2 probes: 4"
# [1] "N duplicated variants to exclude: 313"
# 
# TRUE 
# 592708 
# [1] 486   8
# [1] 485   8
# [1] 0 8

# [1] "niddk_cho_gsa"
# 
# TRUE 
# 608759 
# 
# TRUE 
# 608759 
# [1] "N variants duplicated: 315"
# [1] "N IDs more than 2 probes: 4"
# [1] "N duplicated variants to exclude: 316"
# 
# TRUE 
# 608759 
# [1] 486   8
# [1] 485   8
# [1] 0 8

# [1] "niddk_duerr_gsa"
# 
# TRUE 
# 596287 
# 
# TRUE 
# 596287 
# [1] "N variants duplicated: 319"
# [1] "N IDs more than 2 probes: 5"
# [1] "N duplicated variants to exclude: 320"
# 
# TRUE 
# 596287 
# [1] 488   8
# [1] 487   8
# [1] 0 8

# [1] "niddk_rioux_gsa"
# 
# TRUE 
# 589690 
# 
# TRUE 
# 589690 
# [1] "N variants duplicated: 306"
# [1] "N IDs more than 2 probes: 5"
# [1] "N duplicated variants to exclude: 307"
# 
# TRUE 
# 589690 
# [1] 481   8
# [1] 480   8
# [1] 0 8

# [1] "niddk_silverberg_gsa"
# 
# TRUE 
# 626861 
# 
# TRUE 
# 626861 
# [1] "N variants duplicated: 337"
# [1] "N IDs more than 2 probes: 4"
# [1] "N duplicated variants to exclude: 338"
# 
# TRUE 
# 626861 
# [1] 492   8
# [1] 491   8
# [1] 0 8

# [1] "palotie_hus_gsa"
# 
# TRUE 
# 565416 
# 
# TRUE 
# 565416 
# [1] "N variants duplicated: 294"
# [1] "N IDs more than 2 probes: 4"
# [1] "N duplicated variants to exclude: 295"
# 
# TRUE 
# 565416 
# [1] 454   8
# [1] 453   8
# [1] 0 8

# [1] "pekow_share_gsa"
# 
# TRUE 
# 581240 
# 
# TRUE 
# 581240 
# [1] "N variants duplicated: 301"
# [1] "N IDs more than 2 probes: 5"
# [1] "N duplicated variants to exclude: 302"
# 
# TRUE 
# 581240 
# [1] 482   8
# [1] 481   8
# [1] 0 8

# [1] "rioux_igenomed_gsa"
# 
# TRUE 
# 551803 
# 
# TRUE 
# 551803 
# [1] "N variants duplicated: 282"
# [1] "N IDs more than 2 probes: 3"
# [1] "N duplicated variants to exclude: 283"
# 
# TRUE 
# 551803 
# [1] 447   8
# [1] 446   8
# [1] 0 8

# [1] "sands_msccr_gsa"
# 
# TRUE 
# 604196 
# 
# TRUE 
# 604196 
# [1] "N variants duplicated: 309"
# [1] "N IDs more than 2 probes: 5"
# [1] "N duplicated variants to exclude: 310"
# 
# TRUE 
# 604196 
# [1] 486   8
# [1] 485   8
# [1] 0 8

# [1] "stampfer_gsa"
# 
# TRUE 
# 672634 
# 
# TRUE 
# 672634 
# [1] "N variants duplicated: 465"
# [1] "N IDs more than 2 probes: 5"
# [1] "N duplicated variants to exclude: 466"
# 
# TRUE 
# 672634 
# [1] 513   8
# [1] 512   8
# [1] 0 8

# [1] "vermeire_gsa"
# 
# TRUE 
# 675377 
# 
# TRUE 
# 675377 
# [1] "N variants duplicated: 461"
# [1] "N IDs more than 2 probes: 5"
# [1] "N duplicated variants to exclude: 462"
# 
# TRUE 
# 675377 
# [1] 516   8
# [1] 515   8
# [1] 0 8

# [1] "weersma_gsa"
# 
# TRUE 
# 573894 
# 
# TRUE 
# 573894 
# [1] "N variants duplicated: 294"
# [1] "N IDs more than 2 probes: 5"
# [1] "N duplicated variants to exclude: 295"
# 
# TRUE 
# 573894 
# [1] 464   8
# [1] 463   8
# [1] 0 8

# [1] "xavier_prism_gsa"
# 
# TRUE 
# 584364 
# 
# TRUE 
# 584364 
# [1] "N variants duplicated: 305"
# [1] "N IDs more than 2 probes: 5"
# [1] "N duplicated variants to exclude: 306"
# 
# TRUE 
# 584364 
# [1] 482   8
# [1] 481   8
# [1] 0 8

# [1] "xavier_share_gsa"
# 
# TRUE 
# 584949 
# 
# TRUE 
# 584949 
# [1] "N variants duplicated: 303"
# [1] "N IDs more than 2 probes: 4"
# [1] "N duplicated variants to exclude: 304"
# 
# TRUE 
# 584949 
# [1] 482   8
# [1] 481   8
# [1] 0 8

# [1] "australia_omniexome"
# 
# TRUE 
# 807107
# 
# TRUE 
# 807107 
# [1] "N variants duplicated: 18843"
# [1] "N IDs more than 2 probes: 89"
# [1] "N duplicated variants to exclude: 18877"
# 
# TRUE 
# 807107
# [1] 472   8
# [1] 467   8
# [1] 0 8

# [1] "basque_gsa"
# 
# TRUE 
# 579360 
# 
# TRUE 
# 579360 
# [1] "N variants duplicated: 1112"
# [1] "N IDs more than 2 probes: 304"
# [1] "N duplicated variants to exclude: 1327"
# 
# TRUE 
# 579360 
# [1] 669   8
# [1] 649   8
# [1] 0 8

# [1] "belgium_franchimont_gsa"
# 
# TRUE 
# 585266 
# 
# TRUE 
# 585266 
# [1] "N variants duplicated: 321"
# [1] "N IDs more than 2 probes: 5"
# [1] "N duplicated variants to exclude: 322"
# 
# TRUE 
# 585266 
# [1] 477   8
# [1] 476   8
# [1] 0 8

# [1] "belgium_louis_gsa"
# 
# TRUE 
# 605806 
# 
# TRUE 
# 605806 
# [1] "N variants duplicated: 349"
# [1] "N IDs more than 2 probes: 5"
# [1] "N duplicated variants to exclude: 350"
# 
# TRUE 
# 605806 
# [1] 472   8
# [1] 471   8
# [1] 0 8

# [1] "belgium_vermeire_gsa"
# 
# TRUE 
# 674795 
# 
# TRUE 
# 674795 
# [1] "N variants duplicated: 459"
# [1] "N IDs more than 2 probes: 5"
# [1] "N duplicated variants to exclude: 460"
# 
# TRUE 
# 674795 
# [1] 516   8
# [1] 515   8
# [1] 0 8

# [1] "ccfa_gsa"
# 
# TRUE 
# 611354 
# 
# TRUE 
# 611354 
# [1] "N variants duplicated: 322"
# [1] "N IDs more than 2 probes: 5"
# [1] "N duplicated variants to exclude: 323"
# 
# TRUE 
# 611354 
# [1] 490   8
# [1] 489   8
# [1] 0 8

# [1] "cedars_gsa"
# 
# TRUE 
# 669422 
# 
# TRUE 
# 669422 
# [1] "N variants duplicated: 410"
# [1] "N IDs more than 2 probes: 5"
# [1] "N duplicated variants to exclude: 411"
# 
# TRUE 
# 669422 
# [1] 506   8
# [1] 505   8
# [1] 0 8

# [1] "italy_gsa"
# 
# TRUE 
# 586354 
# 
# TRUE 
# 586354 
# [1] "N variants duplicated: 1111"
# [1] "N IDs more than 2 probes: 299"
# [1] "N duplicated variants to exclude: 1325"
# 
# TRUE 
# 586354 
# [1] 666   8
# [1] 645   8
# [1] 0 8

# [1] "kiel_austria_sibdcs_gsa"
# 
# TRUE 
# 663705 
# 
# TRUE 
# 663705 
# [1] "N variants duplicated: 1435"
# [1] "N IDs more than 2 probes: 342"
# [1] "N duplicated variants to exclude: 1656"
# 
# TRUE 
# 663705 
# [1] 719   8
# [1] 699   8
# [1] 0 8

# [1] "lithuania_gsa"
# 
# TRUE 
# 580098 
# 
# TRUE 
# 580098 
# [1] "N variants duplicated: 306"
# [1] "N IDs more than 2 probes: 5"
# [1] "N duplicated variants to exclude: 307"
# 
# TRUE 
# 580098 
# [1] 443   8
# [1] 442   8
# [1] 0 8

# [1] "mccauley_gsa"
# 
# TRUE 
# 571173 
# 
# TRUE 
# 571173 
# [1] "N variants duplicated: 294"
# [1] "N IDs more than 2 probes: 4"
# [1] "N duplicated variants to exclude: 295"
# 
# TRUE 
# 571173 
# [1] 472   8
# [1] 471   8
# [1] 0 8

# [1] "netherlands_gsa"
# 
# TRUE 
# 684826 
# 
# TRUE 
# 684826 
# [1] "N variants duplicated: 489"
# [1] "N IDs more than 2 probes: 8"
# [1] "N duplicated variants to exclude: 491"
# 
# TRUE 
# 684826 
# [1] 527   8
# [1] 526   8
# [1] 0 8

# [1] "niddk_broad_gsa"
# 
# TRUE 
# 622147 
# 
# TRUE 
# 622147 
# [1] "N variants duplicated: 335"
# [1] "N IDs more than 2 probes: 7"
# [1] "N duplicated variants to exclude: 336"
# 
# TRUE 
# 622147 
# [1] 501   8
# [1] 500   8
# [1] 0 8

# [1] "niddk_feinstein_gsa"
# 
# TRUE 
# 645538 
# 
# TRUE 
# 645538 
# [1] "N variants duplicated: 690"
# [1] "N IDs more than 2 probes: 8"
# [1] "N duplicated variants to exclude: 692"
# 
# TRUE 
# 645538 
# [1] 497   8
# [1] 496   8
# [1] 0 8

# [1] "prism_nfe_gsa"
# 
# TRUE 
# 573760 
# 
# TRUE 
# 573760 
# [1] "N variants duplicated: 294"
# [1] "N IDs more than 2 probes: 4"
# [1] "N duplicated variants to exclude: 295"
# 
# TRUE 
# 573760 
# [1] 472   8
# [1] 471   8
# [1] 0 8

# [1] "slovenia_gsa"
# 
# TRUE 
# 557690 
# 
# TRUE 
# 557690 
# [1] "N variants duplicated: 303"
# [1] "N IDs more than 2 probes: 6"
# [1] "N duplicated variants to exclude: 305"
# 
# TRUE 
# 557690 
# [1] 423   8
# [1] 422   8
# [1] 0 8

# [1] "sweden_gsa"
# 
# TRUE 
# 575192 
# 
# TRUE 
# 575192 
# [1] "N variants duplicated: 306"
# [1] "N IDs more than 2 probes: 7"
# [1] "N duplicated variants to exclude: 308"
# 
# TRUE 
# 575192 
# [1] 465   8
# [1] 464   8
# [1] 0 8

# [1] "helmsley_prism_gsa"
# 
# TRUE 
# 243983 
# 
# TRUE 
# 243983 
# [1] "N variants duplicated: 0"
# [1] "N IDs more than 2 probes: 0"
# [1] "N duplicated variants to exclude: 0"
# 
# TRUE 
# 243983 
# [1] 0 8
# [1] 0 8
# [1] 0 8

# [1] "helmsley_xavier_prism_gsa"
# 
# TRUE 
# 243952 
# 
# TRUE 
# 243952 
# [1] "N variants duplicated: 0"
# [1] "N IDs more than 2 probes: 0"
# [1] "N duplicated variants to exclude: 0"
# 
# TRUE 
# 243952 
# [1] 0 8
# [1] 0 8
# [1] 0 8

# [1] "prism_nfe_gwas"
# 
# TRUE 
# 243983 
# 
# TRUE 
# 243983 
# [1] "N variants duplicated: 0"
# [1] "N IDs more than 2 probes: 0"
# [1] "N duplicated variants to exclude: 0"
# 
# TRUE 
# 243983 
# [1] 0 8
# [1] 0 8
# [1] 0 8

# [1] "finland_illugwas"
# 
# TRUE 
# 239244 
# 
# TRUE 
# 239244 
# [1] "N variants duplicated: 0"
# [1] "N IDs more than 2 probes: 0"
# [1] "N duplicated variants to exclude: 0"
# 
# TRUE 
# 239244 
# [1] 0 8
# [1] 0 8
# [1] 0 8

# [1] "cedars_370k_old_gwas"
# 
# TRUE 
# 339975 
# 
# TRUE 
# 339975 
# [1] "N variants duplicated: 1"
# [1] "N IDs more than 2 probes: 0"
# [1] "N duplicated variants to exclude: 1"
# 
# TRUE 
# 339975 
# [1] 25  8
# [1] 25  8
# [1] 0 8

# [1] "cedars_610k_old_gwas"
# 
# TRUE 
# 586376 
# 
# TRUE 
# 586376 
# [1] "N variants duplicated: 2"
# [1] "N IDs more than 2 probes: 0"
# [1] "N duplicated variants to exclude: 2"
# 
# TRUE 
# 586376 
# [1] 39  8
# [1] 39  8
# [1] 0 8

# [1] "cedars_omni_old_gwas"
# 
# TRUE 
# 721654 
# 
# TRUE 
# 721654 
# [1] "N variants duplicated: 3"
# [1] "N IDs more than 2 probes: 0"
# [1] "N duplicated variants to exclude: 3"
# 
# TRUE 
# 721654 
# [1] 434   8
# [1] 434   8
# [1] 0 8

# [1] "chop_old_gwas"
# 
# TRUE 
# 535308 
# 
# TRUE 
# 535308 
# [1] "N variants duplicated: 0"
# [1] "N IDs more than 2 probes: 0"
# [1] "N duplicated variants to exclude: 0"
# 
# TRUE 
# 535308 
# [1] 0 8
# [1] 0 8
# [1] 0 8

# [1] "german_affy6_old_gwas"
# 
# TRUE 
# 893281 
# 
# TRUE 
# 893281 
# [1] "N variants duplicated: 2918"
# [1] "N IDs more than 2 probes: 0"
# [1] "N duplicated variants to exclude: 2918"
# 
# TRUE 
# 893281 
# [1] 364   8
# [1] 361   8
# [1] 0 8

# [1] "norway_affy6_old_gwas"
# 
# TRUE 
# 846280 
# 
# TRUE 
# 846280 
# [1] "N variants duplicated: 13"
# [1] "N IDs more than 2 probes: 0"
# [1] "N duplicated variants to exclude: 13"
# 
# TRUE 
# 846280 
# [1] 350   8
# [1] 350   8
# [1] 0 8

# [1] "belgium_inf1_old_gwas"
# 
# TRUE 
# 302700 
# 
# TRUE 
# 302700 
# [1] "N variants duplicated: 0"
# [1] "N IDs more than 2 probes: 0"
# [1] "N duplicated variants to exclude: 0"
# 
# TRUE 
# 302700 
# [1] 0 8
# [1] 0 8
# [1] 0 8

# [1] "belgium_inf2_old_gwas"
# 
# TRUE 
# 290730 
# 
# TRUE 
# 290730 
# [1] "N variants duplicated: 0"
# [1] "N IDs more than 2 probes: 0"
# [1] "N duplicated variants to exclude: 0"
# 
# TRUE 
# 290730 
# [1] 0 8
# [1] 0 8
# [1] 0 8
# 
# [1] "swedish_uc_old_gwas"
# 
# TRUE 
# 297004 
# 
# TRUE 
# 297004 
# [1] "N variants duplicated: 0"
# [1] "N IDs more than 2 probes: 0"
# [1] "N duplicated variants to exclude: 0"
# 
# TRUE 
# 297004 
# [1] 0 8
# [1] 0 8
# [1] 0 8

# [1] "pittsburgh_gsa"
# 
# TRUE 
# 945292 
# 
# TRUE 
# 945292 
# [1] "N variants duplicated: 30"
# [1] "N IDs more than 2 probes: 0"
# [1] "N duplicated variants to exclude: 30"
# 
# TRUE 
# 945292 
# [1] 374   8
# [1] 374   8
# [1] 0 8

# [1] "spain_gsa"
# 
# TRUE 
# 585730 
# 
# TRUE 
# 585730 
# [1] "N variants duplicated: 2"
# [1] "N IDs more than 2 probes: 0"
# [1] "N duplicated variants to exclude: 2"
# 
# TRUE 
# 585730 
# [1] 41  8
# [1] 41  8
# [1] 0 8

# [1] "gwas1"
# 
# TRUE 
# 459425 
# 
# TRUE 
# 459425 
# [1] "N variants duplicated: 1"
# [1] "N IDs more than 2 probes: 0"
# [1] "N duplicated variants to exclude: 1"
# 
# TRUE 
# 459425 
# [1] 0 8
# [1] 0 8
# [1] 0 8
# 
# [1] "gwas2"
# 
# TRUE 
# 799263 
# 
# TRUE 
# 799263 
# [1] "N variants duplicated: 2937"
# 
# [1] "N IDs more than 2 probes: 10"
# [1] "N duplicated variants to exclude: 2937"
# 
# TRUE 
# 799263 
# [1] 429   8
# [1] 414   8
# [1] 0 8

##############################

for i in ${studies[@]}
do
echo ${i} && wc -l ${path_gwas}pre_imputation/QC/${i}/list_duplicated_var_exclude
done

##############################

# australia_omniexome
# 18877 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/australia_omniexome/list_duplicated_var_exclude
# gwas1
# 1 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas1/list_duplicated_var_exclude
# gwas2
# 2937 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas2/list_duplicated_var_exclude
# pittsburgh_gsa
# 30 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pittsburgh_gsa/list_duplicated_var_exclude
# spain_gsa
# 2 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/spain_gsa/list_duplicated_var_exclude
# italy_gsa
# 1325 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/italy_gsa/list_duplicated_var_exclude
# kiel_austria_sibdcs_gsa
# 1656 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/kiel_austria_sibdcs_gsa/list_duplicated_var_exclude
# netherlands_gsa
# 491 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/netherlands_gsa/list_duplicated_var_exclude
# slovenia_gsa
# 305 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/slovenia_gsa/list_duplicated_var_exclude
# sweden_gsa
# 308 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sweden_gsa/list_duplicated_var_exclude
# niddk_broad_gsa
# 336 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_broad_gsa/list_duplicated_var_exclude
# niddk_feinstein_gsa
# 692 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_feinstein_gsa/list_duplicated_var_exclude
# basque_gsa
# 1327 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/basque_gsa/list_duplicated_var_exclude
# lithuania_gsa
# 307 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lithuania_gsa/list_duplicated_var_exclude
# belgium_louis_gsa
# 350 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_louis_gsa/list_duplicated_var_exclude
# belgium_franchimont_gsa
# 322 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_franchimont_gsa/list_duplicated_var_exclude
# belgium_vermeire_gsa
# 460 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_vermeire_gsa/list_duplicated_var_exclude
# prism_nfe_gsa
# 295 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gsa/list_duplicated_var_exclude
# prism_nfe_gwas
# 0 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gwas/list_duplicated_var_exclude
# finland_illugwas
# 0 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/finland_illugwas/list_duplicated_var_exclude
# german_affy6_old_gwas
# 2918 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/german_affy6_old_gwas/list_duplicated_var_exclude
# norway_affy6_old_gwas
# 13 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/norway_affy6_old_gwas/list_duplicated_var_exclude
# belgium_inf1_old_gwas
# 0 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf1_old_gwas/list_duplicated_var_exclude
# belgium_inf2_old_gwas
# 0 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf2_old_gwas/list_duplicated_var_exclude
# cedars_370k_old_gwas
# 1 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_370k_old_gwas/list_duplicated_var_exclude
# cedars_610k_old_gwas
# 2 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_610k_old_gwas/list_duplicated_var_exclude
# cedars_omni_old_gwas
# 3 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_omni_old_gwas/list_duplicated_var_exclude
# swedish_uc_old_gwas
# 0 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/list_duplicated_var_exclude
# mccauley_gsa
# 295 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_gsa/list_duplicated_var_exclude
# ccfa_gsa
# 323 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/ccfa_gsa/list_duplicated_var_exclude
# cedars_gsa
# 411 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_gsa/list_duplicated_var_exclude
# bernstein_gsa
# 291 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/bernstein_gsa/list_duplicated_var_exclude
# farkkila_gsa
# 272 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/farkkila_gsa/list_duplicated_var_exclude
# franchimont_gsa
# 452 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franchimont_gsa/list_duplicated_var_exclude
# franke_gsa
# 291 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franke_gsa/list_duplicated_var_exclude
# helmsley_prism_gsa
# 0 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_prism_gsa/list_duplicated_var_exclude
# helmsley_xavier_prism_gsa
# 0 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_xavier_prism_gsa/list_duplicated_var_exclude
# hyams_protect_gsa
# 298 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/hyams_protect_gsa/list_duplicated_var_exclude
# lewis_sparc_gsa
# 332 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lewis_sparc_gsa/list_duplicated_var_exclude
# mccauley_new_gsa
# 311 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_new_gsa/list_duplicated_var_exclude
# mcgovern_gsa
# 371 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mcgovern_gsa/list_duplicated_var_exclude
# moayyedi_imagine_gsa
# 324 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/moayyedi_imagine_gsa/list_duplicated_var_exclude
# newberry_share_gsa
# 313 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/newberry_share_gsa/list_duplicated_var_exclude
# niddk_cho_gsa
# 316 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_cho_gsa/list_duplicated_var_exclude
# niddk_duerr_gsa
# 320 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_duerr_gsa/list_duplicated_var_exclude
# niddk_rioux_gsa
# 307 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_rioux_gsa/list_duplicated_var_exclude
# niddk_silverberg_gsa
# 338 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_silverberg_gsa/list_duplicated_var_exclude
# palotie_hus_gsa
# 295 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/palotie_hus_gsa/list_duplicated_var_exclude
# pekow_share_gsa
# 302 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pekow_share_gsa/list_duplicated_var_exclude
# rioux_igenomed_gsa
# 283 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/rioux_igenomed_gsa/list_duplicated_var_exclude
# sands_msccr_gsa
# 310 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sands_msccr_gsa/list_duplicated_var_exclude
# stampfer_gsa
# 466 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/stampfer_gsa/list_duplicated_var_exclude
# vermeire_gsa
# 462 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/vermeire_gsa/list_duplicated_var_exclude
# weersma_gsa
# 295 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/weersma_gsa/list_duplicated_var_exclude
# xavier_prism_gsa
# 306 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_prism_gsa/list_duplicated_var_exclude
# xavier_share_gsa
# 304 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_share_gsa/list_duplicated_var_exclude

##############################

# NOTE: use pre-vcf fam file which keeps the ca/ctr info:


path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500

for i in ${studies[@]}
do
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_5_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_5_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bed ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr.bed \
--bim ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_edited.bim \
--fam ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind.fam \
--exclude ${path_gwas}pre_imputation/QC/${i}/list_duplicated_var_exclude \
--split-x 'b37' no-fail \
--make-bed --out ${path_gwas}/pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup"
done
# Job <809962..809987> is submitted to queue <normal>.

for i in ${studies[@]}
do
echo ${i} && tail -70 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_5_${i} | grep -E "completed|pass filters|chromosome codes changed"
done

## some need to be repeated 

##############################
# australia_omniexome
# 788230 variants and 1310 people pass filters and QC.
# --split-x: 467 chromosome codes changed.
# Successfully completed.
# gwas1
# 459424 variants and 4684 people pass filters and QC.
# Successfully completed.
# gwas2
# 796326 variants and 7778 people pass filters and QC.
# --split-x: 414 chromosome codes changed.
# Successfully completed.
# pittsburgh_gsa
# 945262 variants and 2781 people pass filters and QC.
# --split-x: 374 chromosome codes changed.
# Successfully completed.
# spain_gsa
# 585728 variants and 3443 people pass filters and QC.
# --split-x: 41 chromosome codes changed.
# Successfully completed.
# italy_gsa
# 585029 variants and 1016 people pass filters and QC.
# --split-x: 645 chromosome codes changed.
# Successfully completed.
# kiel_austria_sibdcs_gsa
# 662049 variants and 14660 people pass filters and QC.
# --split-x: 699 chromosome codes changed.
# Successfully completed.
# netherlands_gsa
# 684335 variants and 4705 people pass filters and QC.
# --split-x: 526 chromosome codes changed.
# Successfully completed.
# slovenia_gsa
# 557385 variants and 270 people pass filters and QC.
# --split-x: 422 chromosome codes changed.
# Successfully completed.
# sweden_gsa
# 574884 variants and 1414 people pass filters and QC.
# --split-x: 464 chromosome codes changed.
# Successfully completed.
# niddk_broad_gsa
# 621811 variants and 5515 people pass filters and QC.
# --split-x: 500 chromosome codes changed.
# Successfully completed.
# niddk_feinstein_gsa
#644846 variants and 8301 people pass filters and QC.
# --split-x: 496 chromosome codes changed.
# Successfully completed.
# basque_gsa
# 578033 variants and 1516 people pass filters and QC.
# --split-x: 649 chromosome codes changed.
# Successfully completed.
# lithuania_gsa
# 579791 variants and 2277 people pass filters and QC.
# --split-x: 442 chromosome codes changed.
# Successfully completed.
# belgium_louis_gsa
# 605456 variants and 1531 people pass filters and QC.
# --split-x: 471 chromosome codes changed.
# Successfully completed.
# belgium_franchimont_gsa
# 584944 variants and 1532 people pass filters and QC.
# --split-x: 476 chromosome codes changed.
# Successfully completed.
# belgium_vermeire_gsa
# 674335 variants and 4014 people pass filters and QC.
# --split-x: 515 chromosome codes changed.
# Successfully completed.
# prism_nfe_gsa
# 573465 variants and 466 people pass filters and QC.
# --split-x: 471 chromosome codes changed.
# Successfully completed.
# prism_nfe_gwas
# 243983 variants and 862 people pass filters and QC.
# Successfully completed.
# finland_illugwas
# 239244 variants and 457 people pass filters and QC.
# Successfully completed.
# german_affy6_old_gwas
# 890363 variants and 2827 people pass filters and QC.
# --split-x: 361 chromosome codes changed.
# Successfully completed.
# norway_affy6_old_gwas
# 846267 variants and 550 people pass filters and QC.
# --split-x: 350 chromosome codes changed.
# Successfully completed.
# belgium_inf1_old_gwas
# 302700 variants and 1417 people pass filters and QC.
# Successfully completed.
# belgium_inf2_old_gwas
# 290730 variants and 272 people pass filters and QC.
# Successfully completed.
# cedars_370k_old_gwas
# 339974 variants and 608 people pass filters and QC.
# --split-x: 25 chromosome codes changed.
# Successfully completed.
# cedars_610k_old_gwas
# 586374 variants and 893 people pass filters and QC.
# --split-x: 39 chromosome codes changed.
# Successfully completed.
# cedars_omni_old_gwas
# 721651 variants and 1226 people pass filters and QC.
# --split-x: 434 chromosome codes changed.
# Successfully completed.
# swedish_uc_old_gwas
# 297004 variants and 1264 people pass filters and QC.
# Successfully completed.
# mccauley_gsa
# 570878 variants and 788 people pass filters and QC.
# --split-x: 471 chromosome codes changed.
# Successfully completed.
# ccfa_gsa
# 611031 variants and 2188 people pass filters and QC.
# --split-x: 489 chromosome codes changed.
# Successfully completed.
# cedars_gsa
# 669011 variants and 3096 people pass filters and QC.
# --split-x: 505 chromosome codes changed.
# Successfully completed.
# bernstein_gsa
# 565467 variants and 514 people pass filters and QC.
# --split-x: 451 chromosome codes changed.
# Successfully completed.
# farkkila_gsa
# 528039 variants and 68 people pass filters and QC.
# --split-x: 402 chromosome codes changed.
# Successfully completed.
# franchimont_gsa
# 670680 variants and 2789 people pass filters and QC.
# --split-x: 514 chromosome codes changed.
# Successfully completed.
# franke_gsa
# 570827 variants and 885 people pass filters and QC.
# --split-x: 446 chromosome codes changed.
# Successfully completed.
# helmsley_prism_gsa
# 243983 variants and 780 people pass filters and QC.
# Successfully completed.
# helmsley_xavier_prism_gsa
# 243952 variants and 1297 people pass filters and QC.
# Successfully completed.
# hyams_protect_gsa
# 580630 variants and 418 people pass filters and QC.
# --split-x: 482 chromosome codes changed.
# Successfully completed.
# lewis_sparc_gsa
# 619563 variants and 2859 people pass filters and QC.
# --split-x: 492 chromosome codes changed.
# Successfully completed.
# mccauley_new_gsa
# 591040 variants and 1628 people pass filters and QC.
# --split-x: 483 chromosome codes changed.
# Successfully completed.
# mcgovern_gsa
# 651675 variants and 6050 people pass filters and QC.
# --split-x: 499 chromosome codes changed.
# Successfully completed.
# moayyedi_imagine_gsa
# 615197 variants and 1128 people pass filters and QC.
# --split-x: 489 chromosome codes changed.
# Successfully completed.
# newberry_share_gsa
# 592395 variants and 865 people pass filters and QC.
# --split-x: 485 chromosome codes changed.
# Successfully completed.
# niddk_cho_gsa
# 608443 variants and 1764 people pass filters and QC.
# --split-x: 485 chromosome codes changed.
# Successfully completed.
# niddk_duerr_gsa
# 595967 variants and 1944 people pass filters and QC.
# --split-x: 487 chromosome codes changed.
# Successfully completed.
# niddk_rioux_gsa
# 589383 variants and 919 people pass filters and QC.
# --split-x: 480 chromosome codes changed.
# Successfully completed.
# niddk_silverberg_gsa
# 626523 variants and 2371 people pass filters and QC.
# --split-x: 491 chromosome codes changed.
# Successfully completed.
# palotie_hus_gsa
# 565121 variants and 878 people pass filters and QC.
# --split-x: 453 chromosome codes changed.
# Successfully completed.
# pekow_share_gsa
# 580938 variants and 634 people pass filters and QC.
# --split-x: 481 chromosome codes changed.
# Successfully completed.
# rioux_igenomed_gsa
# 551520 variants and 182 people pass filters and QC.
# --split-x: 446 chromosome codes changed.
# Successfully completed.
# sands_msccr_gsa
# 603886 variants and 1430 people pass filters and QC.
# --split-x: 485 chromosome codes changed.
# Successfully completed.
# stampfer_gsa
# 672168 variants and 1477 people pass filters and QC.
# --split-x: 512 chromosome codes changed.
# Successfully completed.
# vermeire_gsa
# 674915 variants and 4713 people pass filters and QC.
# --split-x: 515 chromosome codes changed.
# Successfully completed.
# weersma_gsa
# 573599 variants and 709 people pass filters and QC.
# --split-x: 463 chromosome codes changed.
# Successfully completed.
# xavier_prism_gsa
# 584058 variants and 692 people pass filters and QC.
# --split-x: 481 chromosome codes changed.
# Successfully completed.
# xavier_share_gsa
# 584645 variants and 696 people pass filters and QC.
# --split-x: 481 chromosome codes changed.
# Successfully completed.

##############################


##############################################################################################################################################
