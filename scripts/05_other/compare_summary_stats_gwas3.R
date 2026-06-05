# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# # get gwas3 only summary stats:

path_gwas=/path/to/ibdgwas/IIBDGC/


### IBD
cd /path/to/project
cat *.assoc | awk ' $29 > 0.01 ' | cut -d " " -f 1-2,4-6,42 > ${path_gwas}/summary_files/summary_stats_ibd_gwas3_old

wc -l ${path_gwas}/summary_files/summary_stats_ibd_gwas3_old
# 10126380 ${path_gwas}/summary_files/summary_stats_ibd_gwas3_old

sort -k1n -k3n ${path_gwas}/summary_files/summary_stats_ibd_gwas3_old > ${path_gwas}/summary_files/summary_stats_ibd_gwas3_old_sorted
sed 1,22d ${path_gwas}/summary_files/summary_stats_ibd_gwas3_old_sorted > ${path_gwas}/summary_files/summary_stats_ibd_gwas3_old_sorted_noheader
sed -i s/^/chr/ ${path_gwas}summary_files/summary_stats_ibd_gwas3_old_sorted_noheader

awk -F' ' '{OFS=" ";print $1,$3,$3+1,$2,$6}' ${path_gwas}summary_files/summary_stats_ibd_gwas3_old_sorted_noheader >\
${path_gwas}summary_files/summary_stats_ibd_gwas3_old_sorted_noheader_chr_pos_b37


${path_gwas}previous_qced_b38/liftover/liftOver \
${path_gwas}summary_files/summary_stats_ibd_gwas3_old_sorted_noheader_chr_pos_b37 \
${path_gwas}previous_qced_b38/liftover/hg19ToHg38.over.chain.gz \
${path_gwas}summary_files/summary_stats_ibd_gwas3_old_lifted_hg38 \
${path_gwas}summary_files/summary_stats_ibd_gwas3_old_lifted_no_lifted_hg38

wc -l ${path_gwas}summary_files/summary_stats_ibd_gwas3_old_lifted_hg38
# 10121899 ${path_gwas}summary_files/summary_stats_ibd_gwas3_old_lifted_hg38

awk -F' ' '{OFS="";print $1,"\t",$2,"\t","chr",$4,"\t",$5}' ${path_gwas}summary_files/summary_stats_ibd_gwas3_old_lifted_hg38 > \
${path_gwas}summary_files/summary_stats_ibd_gwas3_old_lifted_hg38_ed


### CD
cd /path/to/project
cat *.assoc | awk ' $29 > 0.01 ' | cut -d " " -f 1-2,4-6,42 > ${path_gwas}/summary_files/summary_stats_cd_gwas3_old

wc -l ${path_gwas}/summary_files/summary_stats_cd_gwas3_old
# 10126834 ${path_gwas}/summary_files/summary_stats_cd_gwas3_old

sort -k1n -k3n ${path_gwas}/summary_files/summary_stats_cd_gwas3_old > ${path_gwas}/summary_files/summary_stats_cd_gwas3_old_sorted
sed 1,22d ${path_gwas}summary_files/summary_stats_cd_gwas3_old_sorted > ${path_gwas}summary_files/summary_stats_cd_gwas3_old_sorted_noheader
sed -i s/^/chr/ ${path_gwas}summary_files/summary_stats_cd_gwas3_old_sorted_noheader

awk -F' ' '{OFS=" ";print $1,$3,$3+1,$2,$6}' ${path_gwas}summary_files/summary_stats_cd_gwas3_old_sorted_noheader >\
${path_gwas}summary_files/summary_stats_cd_gwas3_old_sorted_noheader_chr_pos_b37


${path_gwas}previous_qced_b38/liftover/liftOver \
${path_gwas}summary_files/summary_stats_cd_gwas3_old_sorted_noheader_chr_pos_b37 \
${path_gwas}previous_qced_b38/liftover/hg19ToHg38.over.chain.gz \
${path_gwas}summary_files/summary_stats_cd_gwas3_old_lifted_hg38 \
${path_gwas}summary_files/summary_stats_cd_gwas3_old_lifted_no_lifted_hg38

wc -l ${path_gwas}summary_files/summary_stats_cd_gwas3_old_lifted_hg38
# 10122351 ${path_gwas}summary_files/summary_stats_cd_gwas3_old_lifted_hg38

awk -F' ' '{OFS="";print $1,"\t",$2,"\t","chr",$4,"\t",$5}' ${path_gwas}summary_files/summary_stats_cd_gwas3_old_lifted_hg38 > \
${path_gwas}summary_files/summary_stats_cd_gwas3_old_lifted_hg38_ed



### UC
cd /path/to/project
cat *.assoc | awk ' $29 > 0.01 ' | cut -d " " -f 1-2,4-6,42 > ${path_gwas}/summary_files/summary_stats_uc_gwas3_old

wc -l ${path_gwas}/summary_files/summary_stats_uc_gwas3_old
# 10124873 ${path_gwas}/summary_files/summary_stats_uc_gwas3_old

sort -k1n -k3n ${path_gwas}/summary_files/summary_stats_uc_gwas3_old > ${path_gwas}/summary_files/summary_stats_uc_gwas3_old_sorted
sed 1,22d ${path_gwas}summary_files/summary_stats_uc_gwas3_old_sorted > ${path_gwas}summary_files/summary_stats_uc_gwas3_old_sorted_noheader
sed -i s/^/chr/ ${path_gwas}summary_files/summary_stats_uc_gwas3_old_sorted_noheader

awk -F' ' '{OFS=" ";print $1,$3,$3+1,$2,$6}' ${path_gwas}summary_files/summary_stats_uc_gwas3_old_sorted_noheader >\
${path_gwas}summary_files/summary_stats_uc_gwas3_old_sorted_noheader_chr_pos_b37


${path_gwas}previous_qced_b38/liftover/liftOver \
${path_gwas}summary_files/summary_stats_uc_gwas3_old_sorted_noheader_chr_pos_b37 \
${path_gwas}previous_qced_b38/liftover/hg19ToHg38.over.chain.gz \
${path_gwas}summary_files/summary_stats_uc_gwas3_old_lifted_hg38 \
${path_gwas}summary_files/summary_stats_uc_gwas3_old_lifted_no_lifted_hg38

wc -l ${path_gwas}summary_files/summary_stats_uc_gwas3_old_lifted_hg38
# 10120394 ${path_gwas}summary_files/summary_stats_uc_gwas3_old_lifted_hg38

awk -F' ' '{OFS="";print $1,"\t",$2,"\t","chr",$4,"\t",$5}' ${path_gwas}summary_files/summary_stats_uc_gwas3_old_lifted_hg38 > \
${path_gwas}summary_files/summary_stats_uc_gwas3_old_lifted_hg38_ed


MEM=15000
bsub -Is -m "modern_hardware" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

pheno<-c("cd","uc","ibd")

for (i in 1:length(pheno)) {
# for (i in 1:2) {
  
  print(pheno[i])
  
  old<-fread(paste("/path/to/ibdgwas/IIBDGC/summary_files/summary_stats_",pheno[i],"_gwas3_old_lifted_hg38_ed.gz",sep=""),head=F)
  colnames(old)<-c("chr","position_b38","ID_b37","pvalue")
  old<-old[,c("chr","position_b38","ID_b37")]
  old<-old[!duplicated(old),]
  old$X<-paste(gsub("chr","",old$chr),old$position_b38,sep="_")
  old$ID_b37<-gsub("chr","",old$ID_b37)

  # old GWAS3 summary stats
  for (chr in c(1:22)) {
    
    tmp<-fread(paste("/path/to/project",pheno[i],"/",chr,".assoc",sep=""),head=T)
    tmp<-tmp[,c("rsid","info","all_maf","cases_total","controls_total","frequentist_add_beta_1","frequentist_add_se_1","frequentist_add_pvalue")]
    
    tmp<-subset(tmp, rsid %in% old$ID_b37)

     
    if(chr==1) {
      old_tmp<-tmp
    } else {
      old_tmp<-rbind(old_tmp,tmp)
    }
    rm(tmp)
  }
  
  # dim(old_tmp)
  # [1] 10121940        8
  # dim(old)
  # [1] 10120913        4
  
  old<-merge(old,old_tmp,by.x="ID_b37",by.y="rsid")
  
  for (chr in c(1:22)) {
    
    tmp<-fread(paste(path,"post_imputation/analysis/regenie/humancoreexome/",pheno[i],"/chr",chr,"_humancoreexome_step2_",pheno[i],"_eur_sex_PCs_firthse.regenie",sep=""),head=T)
    tmp$X<-paste(tmp$CHROM,tmp$GENPOS,sep="_")
    
    tmp<-subset(tmp, X %in% old$X)
    
    # tmp<-merge(tmp,old,by="X")
    
    if(chr==1) {
      new<-tmp
    } else {
      new<-rbind(new,tmp)
    }
    rm(tmp)
  }
  
  all<-merge(old,new,by="X")
  
  print(dim(all))
  print(dim(old))
  print(dim(new))
  
  rm(old,new)

  fwrite(all,paste("/path/to/ibdgwas/IIBDGC/summary_files/summary_stats_",pheno[i],"_gwas3_old_lifted_hg38_ed_plus_new.tsv.gz",sep=""),
         col.names=T,row.names=F,quote=F,sep="\t")
  
  
}

q("no")

#### plot comparisons:

MEM=8000
bsub -Is -m "modern_hardware" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

pheno<-c("cd","uc","ibd")

# for (i in 1:length(pheno)) {
for (i in 2:3) {
  
  all<-fread(paste("/path/to/ibdgwas/IIBDGC/summary_files/summary_stats_",pheno[i],"_gwas3_old_lifted_hg38_ed_plus_new.tsv.gz",sep=""),head=T)
  print(head(all))

  ## remove some unlikely results - very significant before - no signif now:
  print(dim(all[which(all$frequentist_add_pvalue<1E-10 & all$LOG10P.Y1<6),]))
  
  # summary(all$all_maf[which(all$frequentist_add_pvalue<1E-10 & all$LOG10P.Y1<3)]) - not all low freq
  
  large_discrep<-all[which(all$frequentist_add_pvalue<1E-10 & all$LOG10P.Y1<6),]
  write.table(large_discrep,paste("~/tmp_plots/large_discrepancies_pvalues_gwas3_previous_current_zoom_PCs_firthse_",pheno[i],".tsv",sep=""),
              col.names=T,row.names=F,quote=F,sep="\t")


  # plot only variants p-value<0.05 and maf >=0.01
  tmp<-all[which(all$frequentist_add_pvalue<=0.05 & all$all_maf>=0.01),]
  tmp<-tmp[order(tmp$frequentist_add_pvalue),]
  
  pdf(paste("~/tmp_plots/pvalues_gwas3_previous_current_zoom_PCs_firthse_",pheno[i],"_maf0.01.pdf",sep=""),width=7, height=7)
  print(ggplot(tmp, aes(x=-log10(frequentist_add_pvalue), y=(LOG10P.Y1))) + geom_point(shape=20) + xlim(0,121) + ylim(0,121) + 
    geom_abline(intercept = 0, slope = 1, color="red",linetype="dashed", size=1) + ggtitle(pheno[i]))
  dev.off()
  
  rm(tmp,all)
}

q("no")

#####################

MEM=4000
bsub -Is -m "modern_hardware" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1


library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)


path<-"/path/to/ibdgwas/IIBDGC/"
pheno<-c("cd","uc","ibd")

for (i in 1:length(pheno)) {
  tmp<-read.table(paste("~/tmp_plots/large_discrepancies_pvalues_gwas3_previous_current_zoom_PCs_firthse_",pheno[i],".tsv",sep=""),head=T)
  tmp$pheno<-pheno[i]
  
  if(i==1) {
    all<-tmp
  }else{
    all<-rbind(all,tmp)
  }
}

table(all$pheno)
# cd ibd  uc 
# 90 118 104

# combine into regions:

regions<-all[,c("ID_b37"),drop=F]

regions$chr<-as.numeric(gsub(":.*","",regions$ID_b37))
regions$position<-gsub(".*:","",regions$ID_b37)
regions$position<-gsub("_.*","",regions$position)
regions<-regions[order(regions$chr,all$position,decreasing=F),]

regions$start<-as.numeric(regions$position)
regions$end<-as.numeric(regions$position)+1

regions<-regions[,c("chr","start","end","ID_b37")]
regions<-regions[!duplicated(regions),]

write.table(regions,paste(path,"post_imputation/analysis/metaanalysis/summary_results/list_large_discrepancies_pvalues_gwas3_previous_current_PCs_firthse.bed",sep=""),col.names=T,row.names=F,quote=F,sep="\t")


# get temporary new region to look into, by merging into one variants located <5000000 bp away:
system(paste("/path/to/software/bedtools merge -i ",path,"post_imputation/analysis/metaanalysis/summary_results/list_large_discrepancies_pvalues_gwas3_previous_current_PCs_firthse.bed -d 500000 > ",
             path,"post_imputation/analysis/metaanalysis/summary_results/list_regions_large_discrepancies_pvalues_gwas3_previous_current_PCs_firthse.bed",sep=""))


reg<-read.table(paste(path,"post_imputation/analysis/metaanalysis/summary_results/list_regions_large_discrepancies_pvalues_gwas3_previous_current_PCs_firthse.bed",sep=""),head=F)
dim(reg)
# [1] 26  3

reg$region<-paste(reg$V1,reg$V2,reg$V3,sep="_")

all$position_b37<-gsub(".*:","",all$ID_b37)
all$position_b37<-gsub("_.*","",all$position_b37)
  
all$reg<-NA

for (i in 1:nrow(reg)) {
  all$reg[which( all$chr==paste("chr",reg$V1[i],sep="") & 
                   (all$position_b37>=reg$V2[i] & all$position_b37<=reg$V3[i]))]<-reg$region[i]
}

table(all$reg,all$pheno)
#                        cd ibd uc
# 1_111778824_111780908  15  15 15
# 11_56137709_56227005    0   1  3
# 1_156633358_156640679   0   2  2
# 1_169511555_169511556   0   0  1
# 1_180226284_180227224   6   6  6
# 12_104296808_104304982  5   5  4
# 12_38397036_38397037    0   1  0
# 13_52658068_52658069    0   1  1
# 14_104550688_104619120 42  11 40
# 16_2205821_2223001      4   5  4
# 17_6744513_6748798      0   3  4
# 18_2314077_2314078      1   1  0
# 18_53012155_53134557    7  27  5
# 19_1186316_1186317      1   0  0
# 2_173423066_173424408   2   0  0
# 2_7129248_7129249       0   1  1
# 4_100205819_100205830   0   2  2
# 4_84165738_84206005     0   2  0
# 5_110998611_111010074   0   4  1
# 5_40394661_40451191     1   1  0
# 6_113705570_113710638   5   7  0
# 6_32577563_32652306     0  11  7
# 6_57498623_57537585     0   4  0
# 7_35987379_36001018     0   4  5
# 8_21982997_21986851     1   2  1
# 9_135542694_135542776   0   2  2

q("no")

###########################################################################################################
# extract sum stats from GWAS index variants in regions <500Kb from these regions:

MEM=8000
bsub -Is -m "modern_hardware" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"
pheno<-c("cd","uc","ibd")

ids<-c("chr1:169549811:C:T","chr5:40410482:T:C","chr19:1123653:G:A")


for (i in 1:length(pheno)) {
  tmp<-fread(paste("/path/to/ibdgwas/IIBDGC/summary_files/summary_stats_",pheno[i],"_gwas3_old_lifted_hg38_ed_plus_new.tsv.gz",sep=""),head=T)
  tmp<-tmp[which(tmp$ID %in% ids),]

  tmp$pheno<-pheno[i]
  
  if(i==1){
    all<-tmp
  }else {
    all<-rbind(all,tmp)
  } 
}

write.table(all,paste("~/tmp_plots/nearby_gwas_hits_large_discrepancies_pvalues_gwas3_previous_current_zoom_PCs_firthse_",pheno[i],".tsv",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")


q("no")


###########################################################################################################
# Plot pvalue from GWAS index variants

MEM=10000
bsub -Is -m "modern_hardware" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"
pheno<-c("cd","uc","ibd")

ids<-c("chr19:1124032:G:A","chr1:1312114:T:C","chr1:2569899:C:T","chr1:7962137:G:T","chr1:19845367:G:A","chr1:22375738:T:C",
       "chr1:62583922:C:G","chr1:67242007:G:A","chr1:70529879:G:A","chr1:78157942:C:T","chr1:92088726:G:A","chr1:101000498:C:T",
       "chr1:119908567:T:C","chr1:151829204:G:A","chr1:155260340:C:T","chr1:159830120:G:A","chr1:160887174:A:G","chr1:161509955:A:G",
       "chr1:169549811:C:T","chr1:172884320:C:T","chr1:186906327:G:T","chr1:191590226:A:T","chr1:197732149:G:C","chr1:198629533:G:A",
       "chr1:200132792:A:G","chr1:200908434:C:A","chr1:206766559:G:A","chr1:209797265:G:A","chr2:5523876:G:T","chr2:24874775:A:G",
       "chr2:27508073:T:C","chr2:28391927:C:T","chr2:43579779:C:T","chr2:60977721:A:G","chr2:62325186:A:G","chr2:65440138:C:A",
       "chr2:102322576:C:T","chr2:144734815:C:G","chr2:159937497:A:G","chr2:162254026:A:G","chr2:181443625:A:G","chr2:186711651:A:AT",
       "chr2:191066738:A:C","chr2:198032171:A:G","chr2:198696033:A:G","chr2:203727298:G:A","chr2:218286495:A:C","chr2:227795396:C:T",
       "chr2:227805760:C:G","chr2:230232414:T:C","chr2:233252802:G:A","chr2:240634984:G:A","chr2:241545286:C:T","chr2:241797926:T:C",
       "chr3:18725912:C:T","chr3:46415921:T:C","chr3:49684099:G:A","chr3:53028645:T:C","chr3:53099133:C:G","chr3:71126344:C:T",
       "chr3:101304904:A:G","chr3:101850882:A:C","chr3:141386728:A:G","chr3:188683372:C:T","chr4:3442776:G:A","chr4:26130739:T:C",
       "chr4:38323415:T:C","chr4:38586832:G:A","chr4:48361228:A:G","chr4:73991991:A:G","chr4:101944147:G:A","chr4:105154341:T:G",
       "chr4:122430276:A:T","chr5:583327:G:A","chr5:10695414:T:C","chr5:35876172:A:G","chr5:38867630:C:T","chr5:40410482:T:C",
       "chr5:56143024:C:T","chr5:72398072:C:T","chr5:73244023:C:A","chr5:96917099:G:A","chr5:102611094:A:G","chr5:130681594:T:C",
       "chr5:132435113:C:T","chr5:135107916:G:A","chr5:142133639:A:T","chr5:150226431:C:T","chr5:150898347:A:G","chr5:159400761:C:A",
       "chr5:172897975:G:A","chr5:173852839:C:T","chr5:177361569:T:C","chr6:382559:G:A","chr6:3433084:T:G","chr6:14719265:G:A",
       "chr6:19780778:C:T","chr6:20812357:G:A","chr6:21430499:G:T","chr6:32644620:A:G","chr6:42039665:C:T","chr6:43828231:G:T",
       "chr6:90263440:C:A","chr6:105987394:G:C","chr6:111526988:A:T","chr6:116447754:A:G","chr6:127134977:A:G","chr6:127956006:A:G",
       "chr6:137685367:G:A","chr6:143577757:A:G","chr6:149255943:A:G","chr6:159069404:C:T","chr6:166960059:T:G","chr7:2750246:A:C",
       "chr7:6505557:A:G","chr7:17403055:G:A","chr7:20537675:G:A","chr7:26852821:G:A","chr7:28151111:ATTT:A","chr7:50264865:C:T",
       "chr7:99161494:C:A","chr7:100717894:A:G","chr7:107839870:T:C","chr7:117252792:T:C","chr7:128933913:G:A","chr7:148523356:A:G",
       "chr7:148738247:A:G","chr8:27370037:A:G","chr8:48216682:G:A","chr8:73095112:G:T","chr8:89863690:T:C","chr8:125522429:A:G",
       "chr8:128554935:T:C","chr8:129611859:G:C","chr9:4981602:C:A","chr9:34736161:A:G","chr9:91166134:C:T","chr9:114783386:C:T",
       "chr9:117713024:A:G","chr9:136371953:G:A","chr10:6052734:C:T","chr10:26890667:T:G","chr10:30439172:T:C","chr10:35177257:T:G",
       "chr10:58153390:C:T","chr10:62685804:A:G","chr10:73913343:T:C","chr10:79300560:C:A","chr10:80494291:A:G","chr10:92677094:T:C",
       "chr10:99524480:T:G","chr10:102472959:G:A","chr10:110426390:C:T","chr10:124750812:G:A","chr10:131373856:G:A","chr11:1852842:G:A",
       "chr11:57435536:C:T","chr11:58641214:G:T","chr11:61008737:C:T","chr11:61796827:G:T","chr11:64382898:T:C","chr11:65889093:C:T",
       "chr11:76588605:C:A","chr11:87414396:A:G","chr11:96290263:G:A","chr11:118883644:C:A","chr11:128511079:C:T","chr12:6381959:G:A",
       "chr12:12504579:A:G","chr12:40346421:T:C","chr12:47814585:T:C","chr12:68114342:G:A","chr12:111569952:C:T","chr12:119709120:G:A",
       "chr13:26957130:T:C","chr13:40439840:T:C","chr13:40983974:A:G","chr13:42343725:T:C","chr13:42478744:G:A","chr13:43883789:A:G",
       "chr13:49021195:A:G","chr13:99304368:T:C","chr14:68743482:G:A","chr14:75275048:C:T","chr14:88006251:C:T","chr15:38606989:T:C",
       "chr15:41271752:A:G","chr15:67150258:C:T","chr15:90629669:C:T","chr16:11279463:C:T","chr16:23853269:T:C","chr16:28517460:G:A",
       "chr16:30471173:T:C","chr16:50301163:C:A","chr16:50712015:C:T","chr16:68557327:A:C","chr16:81883307:A:G","chr16:82833851:C:A",
       "chr16:85980635:T:C","chr17:27542007:A:C","chr17:34266955:G:A","chr17:39756124:C:T","chr17:42375526:A:G","chr17:56803632:T:C",
       "chr17:59886176:A:G","chr17:72646784:A:C","chr17:78741036:A:G","chr18:12818923:G:A","chr18:48868651:A:G","chr18:59212595:T:C",
       "chr18:69863203:G:A","chr18:79460616:C:T","chr19:1124032:G:A","chr19:10402235:G:A","chr19:33240645:G:A","chr19:34165501:C:T",
       "chr19:46346549:G:T","chr19:48702915:C:T","chr19:54871595:G:A","chr20:6113242:T:G","chr20:32137845:A:G","chr20:32788476:T:G",
       "chr20:35211477:T:G","chr20:44436388:A:C","chr20:46113425:C:A","chr20:50338887:T:C","chr20:59249254:A:G","chr20:63697746:G:A",
       "chr21:15445619:G:A","chr21:33404389:A:G","chr21:39093608:G:A","chr22:21568615:G:T","chr22:30097893:G:A","chr22:35333728:G:A",
       "chr22:36862461:T:C","chr22:39263768:C:T","chr22:41471373:C:T","chr22:49997051:A:G")


for (i in 2:length(pheno)) {
  
  tmp<-fread(paste("/path/to/ibdgwas/IIBDGC/summary_files/summary_stats_",pheno[i],"_gwas3_old_lifted_hg38_ed_plus_new.tsv.gz",sep=""),head=T)
  
  
  tmp<-tmp[which(tmp$ID %in% ids),]
  tmp$pvalue<-10^-(tmp$LOG10P.Y1)
  
  dim(tmp)
  dim(tmp[which(tmp$frequentist_add_pvalue<=5E-8)])
  tmp_s<-tmp[which(tmp$frequentist_add_pvalue<=5E-8 | tmp$pvalue<=5E-8),]
  
  ## PVALUES
  
  print(paste("N significant variants",pheno[i],nrow(tmp_s)))
  print(paste("N more significant new results",pheno[i],nrow(tmp_s[which(tmp_s$pvalue<tmp_s$frequentist_add_pvalue),])))
  print(paste("N more significant old results",pheno[i],nrow(tmp_s[which(tmp_s$pvalue>tmp_s$frequentist_add_pvalue),])))
  
  p11<-ggplot(tmp, aes(x=-log10(frequentist_add_pvalue), y=(LOG10P.Y1))) + geom_point(shape=20) + xlim(0,40) + ylim(0,40) + 
              geom_abline(intercept = 0, slope = 1, color="red",linetype="dashed", size=0.5) +
              ggtitle(paste("P-value All GWAS index variants",pheno[i])) + geom_smooth(method='lm', color="darkgrey",linetype="dashed", size=0.8)
             
  
  p12<-ggplot(tmp_s, aes(x=-log10(frequentist_add_pvalue), y=(LOG10P.Y1))) + geom_point(shape=20) + xlim(0,40) + ylim(0,40) + 
              geom_abline(intercept = 0, slope = 1, color="red",linetype="dashed", size=0.5) +
              ggtitle(paste("P-value GWAS3 significant variants",pheno[i])) + geom_smooth(method='lm', color="darkgrey",linetype="dashed", size=0.8)
  
  p1<-ggarrange(p11,p12,ncol=2)
  
  ## BETAS
  
  beta_lim<-max(abs(tmp_s$BETA.Y1),abs(tmp_s$frequentist_add_beta_1))
    
  print(paste("N larger beta new results",pheno[i],nrow(tmp_s[which(abs(tmp_s$BETA.Y1)>abs(tmp_s$frequentist_add_beta_1)),])))
  print(paste("N larger beta old results",pheno[i],nrow(tmp_s[which(abs(tmp_s$BETA.Y1)<abs(tmp_s$frequentist_add_beta_1)),])))
  
  p21<-ggplot(tmp, aes(x=frequentist_add_beta_1, y=(BETA.Y1))) + geom_point(shape=20) + xlim(-beta_lim,beta_lim) + ylim(-beta_lim,beta_lim) +
    geom_abline(intercept = 0, slope = 1, color="red",linetype="dashed", size=0.5) +
    ggtitle(paste("Beta All GWAS index variants",pheno[i])) + geom_smooth(method='lm', color="darkgrey",linetype="dashed", size=0.8)
  
  p22<-ggplot(tmp_s, aes(x=frequentist_add_beta_1, y=(BETA.Y1))) + geom_point(shape=20) + xlim(-beta_lim,beta_lim) + ylim(-beta_lim,beta_lim) + 
    geom_abline(intercept = 0, slope = 1, color="red",linetype="dashed", size=0.5) +
    ggtitle(paste("Beta GWAS3 significant variants",pheno[i])) + geom_smooth(method='lm', color="darkgrey",linetype="dashed", size=0.8)

  p2<-ggarrange(p21,p22,ncol=2)

  pdf(paste("~/tmp_plots/pvalues_beta_significant_known_gwas_index_variants_gwas3_previous_current_zoom_PCs_firthse_",pheno[i],"_maf0.01.pdf",sep=""),width=10, height=10)
  print(ggarrange(p1,p2,nrow=2))
  dev.off()
  
  print(paste("Correlation p-value",cor(tmp$frequentist_add_pvalue, tmp$pvalue)))
  print(paste("Correlation beta",cor(tmp$frequentist_add_beta_1, tmp$BETA.Y1)))
  
  ## compare imputation info
  
  
  if(pheno[i]=="ibd"){
    
    tmp$maf_new<-pmin(tmp$A1FREQ,1-tmp$A1FREQ)
    
    p31<-ggplot(tmp, aes(x=info, y=(INFO))) + geom_point(shape=20) + xlim(0.6,1.1) + ylim(0.6,1.1) +
      geom_abline(intercept = 0, slope = 1, color="red",linetype="dashed", size=0.5) +
      ggtitle(paste("Imputation info All GWAS index variants",pheno[i])) + geom_smooth(method='lm', color="darkgrey",linetype="dashed", size=0.8)
    
    p32<-ggplot(tmp, aes(x=all_maf, y=(maf_new))) + geom_point(shape=20) + xlim(0,0.5) + ylim(0,0.5) +
      geom_abline(intercept = 0, slope = 1, color="red",linetype="dashed", size=0.5) +
      ggtitle(paste("MAF All GWAS index variants",pheno[i])) + geom_smooth(method='lm', color="darkgrey",linetype="dashed", size=0.8)
    
    p3<-ggarrange(p31,p32,ncol=2)
    
    pdf(paste("~/tmp_plots/info_maf_significant_known_gwas_index_variants_gwas3_previous_current_zoom_PCs_firthse_",pheno[i],"_maf0.01.pdf",sep=""),width=10, height=5)
    print(p3)
    dev.off()
    
    
    print(paste("Correlation info",cor(tmp$info, tmp$INFO)))
    print(paste("Correlation maf",cor(tmp$all_maf, tmp$maf_new)))
    
  }
  
}

q("no")



# |--------------------------------------------------|
# |==================================================|
# [1] "N significant variants cd 17"
# [1] "N more significant new results cd 10"
# [1] "N more significant old results cd 7"
# [1] "Correlation p-value 0.7494448148788"
# [1] "Correlation beta 0.97361687573511"
# 
# |--------------------------------------------------|
# |==================================================|
# [1] "N significant variants uc 12"
# [1] "N more significant new results uc 6"
# [1] "N more significant old results uc 6"
# [1] "Correlation p-value 0.763146356862267"
# [1] "Correlation beta 0.956028048639135"
# |--------------------------------------------------|
# |==================================================|
# [1] "N significant variants ibd 23"
# [1] "N more significant new results ibd 16"
# [1] "N more significant old results ibd 7"
# [1] "Correlation p-value 0.840567060090743"
# [1] "Correlation beta 0.97429158009529"


# [1] "Correlation info 0.679702281493886"
# [1] "Correlation maf 0.97828448506909"



###########################################################################################################
# Plot imputation info

MEM=10000
bsub -Is -m "modern_hardware" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(dplyr)
library(stringr)
library(viridis)
library(ggExtra)
library(scales)
library(reshape2)

path<-"/path/to/ibdgwas/IIBDGC/"
pheno<-c("ibd")
i=1

# old GWAS3 summary stats
for (chr in c(1:22)) {
  
  tmp<-fread(paste("/path/to/project",pheno[i],"/",chr,".assoc",sep=""),head=T)
  tmp<-tmp[,c("rsid","info","all_maf","cases_total","controls_total","frequentist_add_beta_1","frequentist_add_se_1","frequentist_add_pvalue")]
  
  if(chr==1) {
    old<-tmp
  } else {
    old<-rbind(old,tmp)
  }
  rm(tmp)
}

dim(old)
# [1] 18954471        8

summary(old$info)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4000  0.6470  0.8196  0.7931  0.9683  1.0000 
summary(old$all_maf)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001001 0.002352 0.013842 0.094026 0.144173 0.500000 


sum_tmp1<-old %>%
  group_by(MAF_range = cut(all_maf, breaks = c(1E-3,2.5E-3,5E-3,7.5E-3,1E-2,5E-2,2.5E-2,1E-1,2E-1,3E-1,4E-1,1E0), 
                           right = F)) %>% 
  summarise(info_mean = mean(info),
            info_n = n())

sum_tmp1
# # A tibble: 11 x 3
# MAF_range      info_mean  info_n
# <fct>              <dbl>   <int>
# 1 [0.001,0.0025)     0.633 4980857
# 2 [0.0025,0.005)     0.683 2221764
# 3 [0.005,0.0075)     0.724 1004180
# 4 [0.0075,0.01)      0.750  621048
# 5 [0.01,0.025)       0.803 1755618
# 6 [0.025,0.05)       0.866 1290961
# 7 [0.05,0.1)         0.908 1433781
# 8 [0.1,0.2)          0.936 1843579
# 9 [0.2,0.3)          0.949 1397900
# 10 [0.3,0.4)          0.952 1236068
# 11 [0.4,1)            0.950 1168715

rm(old)

for (chr in c(1:22)) {
  
  tmp<-fread(paste(path,"post_imputation/analysis/regenie/humancoreexome/",pheno[i],"/chr",chr,"_humancoreexome_step2_",pheno[i],"_eur_sex_PCs_firthse.regenie",sep=""),head=T)

  if(chr==1) {
    new<-tmp
  } else {
    new<-rbind(new,tmp)
  }
  rm(tmp)
}

dim(new)
# [1] 14056619       12

summary(new$INFO)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3417  0.8421  0.9366  0.8932  0.9781  1.0962 

new$all_maf<-pmin(new$A1FREQ,1-new$A1FREQ)
summary(new$all_maf)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0007814 0.0039496 0.0300564 0.1099355 0.1849230 0.5000000 

new<-new[which(new$all_maf>=0.001),]


sum_tmp2<-new %>%
  group_by(MAF_range = cut(all_maf, breaks = c(1E-3,2.5E-3,5E-3,7.5E-3,1E-2,5E-2,2.5E-2,1E-1,2E-1,3E-1,4E-1,1E0), 
                           right = F)) %>% 
  summarise(info_mean = mean(INFO),
            info_n = n())



sum_df<-as.data.frame(sum_tmp1)
sum_df$group<-"1000GP|HRC|UK10K"
sum_df$tmp_maf<-gsub("\\[","",sum_df$MAF_range)
sum_df$tmp_maf<-gsub(")","",sum_df$tmp_maf)

xx<-as.data.frame(str_split(sum_df$tmp_maf,",",simplify=T))
xx$V1<-as.numeric(as.character(xx$V1))
xx$V2<-as.numeric(as.character(xx$V2))

sum_df$MAF<-NA
for(xxx in 1:nrow(xx)) {
  sum_df$MAF[xxx]<-sum(xx$V1[xxx],xx$V2[xxx])/2
}

sum_df_old<-sum_df

sum_df<-as.data.frame(sum_tmp2)
sum_df$group<-"TOPMed"
sum_df$tmp_maf<-gsub("\\[","",sum_df$MAF_range)
sum_df$tmp_maf<-gsub(")","",sum_df$tmp_maf)

xx<-as.data.frame(str_split(sum_df$tmp_maf,",",simplify=T))
xx$V1<-as.numeric(as.character(xx$V1))
xx$V2<-as.numeric(as.character(xx$V2))

sum_df$MAF<-NA
for(xxx in 1:nrow(xx)) {
  sum_df$MAF[xxx]<-sum(xx$V1[xxx],xx$V2[xxx])/2
}

sum_df_new<-sum_df

sum_df<-rbind(sum_df_old,sum_df_new)

pdf("~/imputation_common_variants_old_new_gwas3.pdf",width = 11, height = 6)
ggplot(sum_df, aes(x = -log10(MAF), y = info_mean, group = group, color = group)) + ylim(0,1) + ylab("Info") +
  geom_point() + geom_line() + expand_limits(x=c(-log10(1E-2),-log10(1E0)) ) + scale_x_continuous(name ="MAF",trans='reverse', 
                                                                                                  breaks= -log10(c(1E-2,1E-1,1E0)) ,
                                                                                                  labels=scientific(c(1E-2,1E-1,1E0)))
dev.off()
  



#####################

MEM=4000
bsub -Is -m "modern_hardware" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1


library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

nod2<-c("16:50335074","16:50663477","16:50694011","16:50744688","16:50745926","16:50746199","16:50750810","16:50756540","16:50762771","16:50763778","16:50816078")
tomatch<-paste(nod2,collapse="|")

old<-read.table("/path/to/project",head=T)
nod2_old<-old[grep(tomatch,old$rsid),]
write.table(nod2_old,"~/tmp_plots/old_gwas3_cd_nod2.tsv",col.names=T,row.names=F,quote=F)

old_meta<-read.table("/path/to/project",head=T)
nod2_old_meta<-old_meta[grep(tomatch,old_meta$MarkerName),]
write.table(nod2_old_meta,"~/tmp_plots/old_iibdgc_meta_cd_nod2.tsv",col.names=T,row.names=F,quote=F)









# all[which(all$V4==min(all$V4)),]
# #                X   V1       V2                V3           V4              ID.x
# # 1: chr8_22125484 chr8 22125484 chr8:21982997_A_G 3.30425e-121 chr8:21982997:A:G
# #    CHROM   GENPOS              ID.y ALLELE0 ALLELE1   A1FREQ TEST    BETA.Y1
# # 1:     8 22125484 chr8:22125484:A:G       A       G 0.107638  ADD -0.0108009
# # not among results they report, remove:
# all<-all[!(X=="chr8_22125484")]
# 
# all[which(all$V4==min(all$V4)),]
# # X   V1       V2                V3           V4              ID.x
# # 1: chr4_83284851 chr4 83284851 chr4:84206004_T_A 2.22756e-106 chr4:84206004:T:A
# # CHROM   GENPOS              ID.y ALLELE0 ALLELE1    A1FREQ TEST  BETA.Y1
# # 1:     4 83284851 chr4:83284851:T:A       T       A 0.0175505  ADD 0.095986
# # SE.Y1 CHISQ.Y1 LOG10P.Y1
# # 1: 0.0797239  1.44957  0.640929
# # not among results they report, remove:
# all<-all[!(X=="chr4_83284851")]


# all[which(all$V4==min(all$V4)),]
# X   V1        V2                 V3           V4
# 1: chr1_111237792 chr1 111237792 chr1:111780414_G_A 8.56854e-102
# ID.x CHROM    GENPOS               ID.y ALLELE0 ALLELE1
# 1: chr1:111780414:G:A     1 111237792 chr1:111237792:G:A       G       A
# A1FREQ   INFO TEST   BETA.Y1     SE.Y1 CHISQ.Y1 LOG10P.Y1
# 1: 0.421552 0.9619  ADD 0.0303373 0.0217733  1.94135  0.786424

# dim(all)
# # [1] 8814278      18
# # dim(all[!(V4<1E-10 & LOG10P.Y1<5)])
# # [1] 8814170      18
# 
# all<-all[!(V4<1E-10 & LOG10P.Y1<5)]
# summary(all$LOG10P.Y1[which(all$V4>0.05)])
# 

######## WHICH ARE THE VARIANTS NOT PRESENT IN CURRENT ANALYSIS

path_gwas=/path/to/ibdgwas/IIBDGC/

### IBD
cd /path/to/project
cat *.assoc | cut -d " " -f 1-2,4-6,42 > ${path_gwas}/summary_files/summary_stats_ibd_gwas3_old_all_variants 

wc -l ${path_gwas}/summary_files/summary_stats_ibd_gwas3_old_all_variants 
# 18954493 /path/to/ibdgwas/IIBDGC//summary_files/summary_stats_ibd_gwas3_old_all_variants

sort -k1n -k3n ${path_gwas}/summary_files/summary_stats_ibd_gwas3_old_all_variants > ${path_gwas}/summary_files/summary_stats_ibd_gwas3_old_all_variants_sorted
sed 1,22d ${path_gwas}/summary_files/summary_stats_ibd_gwas3_old_all_variants_sorted > ${path_gwas}/summary_files/summary_stats_ibd_gwas3_old_all_variants_sorted_noheader
sed -i s/^/chr/ ${path_gwas}summary_files/summary_stats_ibd_gwas3_old_all_variants_sorted_noheader

awk -F' ' '{OFS=" ";print $1,$3,$3+1,$2,$6}' ${path_gwas}summary_files/summary_stats_ibd_gwas3_old_all_variants_sorted_noheader >\
${path_gwas}summary_files/summary_stats_ibd_gwas3_old_all_variants_sorted_noheader_chr_pos_b37


${path_gwas}previous_qced_b38/liftover/liftOver \
${path_gwas}summary_files/summary_stats_ibd_gwas3_old_all_variants_sorted_noheader_chr_pos_b37 \
${path_gwas}previous_qced_b38/liftover/hg19ToHg38.over.chain.gz \
${path_gwas}summary_files/summary_stats_ibd_gwas3_old_all_variants_lifted_hg38 \
${path_gwas}summary_files/summary_stats_ibd_gwas3_old_all_variants_lifted_no_lifted_hg38

wc -l ${path_gwas}summary_files/summary_stats_ibd_gwas3_old_all_variants_lifted_hg38
# 18948607 ${path_gwas}summary_files/summary_stats_ibd_gwas3_old_lifted_hg38

awk -F' ' '{OFS="";print $1,"\t",$2,"\t","chr",$4,"\t",$5}' ${path_gwas}summary_files/summary_stats_ibd_gwas3_old_all_variants_lifted_hg38 > \
${path_gwas}summary_files/summary_stats_ibd_gwas3_old_all_variants_lifted_hg38_ed

gzip ${path_gwas}summary_files/summary_stats_ibd_gwas3_old_all_variants_lifted_hg38_ed


######


MEM=30000
bsub -Is -m "modern_hardware" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

pheno<-c("cd","uc","ibd")

i=3

for (i in 1:length(pheno)) {
  # for (i in 1:2) {
  
  print(pheno[i])
  
  old<-fread(paste("/path/to/ibdgwas/IIBDGC/summary_files/summary_stats_",pheno[i],"_gwas3_old_all_variants_lifted_hg38_ed.gz",sep=""),head=F)
  colnames(old)<-c("chr","position_b38","ID_b37","pvalue")
  old<-old[,c("chr","position_b38","ID_b37")]
  old<-old[!duplicated(old),]
  old$X<-paste(gsub("chr","",old$chr),old$position_b38,sep="_")
  old$ID_b37<-gsub("chr","",old$ID_b37)
  
  print(nrow(old))
  # [1] 18947291
  
  # new resuls
  for (chr in c(1:22)) {
    
    tmp<-fread(paste(path,"post_imputation/analysis/regenie/humancoreexome/",pheno[i],"/chr",chr,"_humancoreexome_step2_",pheno[i],"_eur_sex_PCs_firthse.regenie",sep=""),head=T)
    tmp$X<-paste(tmp$CHROM,tmp$GENPOS,sep="_")
    
    # keep all new variants
    # tmp<-subset(tmp, X %in% old$X)
  
    if(chr==1) {
      new<-tmp
    } else {
      new<-rbind(new,tmp)
    }
    rm(tmp)
  }
  
  print(nrow(new))
  # [1] 14056619
  
  
  # old GWAS3 summary stats
  for (chr in c(1:22)) {
    
    tmp<-fread(paste("/path/to/project",pheno[i],"/",chr,".assoc",sep=""),head=T)
    tmp<-tmp[,c("rsid","info","all_maf","cases_total","controls_total","frequentist_add_beta_1","frequentist_add_se_1","frequentist_add_pvalue")]
    
    tmp<-subset(tmp, rsid %in% old$ID_b37)
    
    if(chr==1) {
      old_tmp<-tmp
    } else {
      old_tmp<-rbind(old_tmp,tmp)
    }
    rm(tmp)
  }
  
  print(nrow(old_tmp))
  # [1] 18948607        8
  old_tmp<-old_tmp[!duplicated(old_tmp),]
  print(nrow(old_tmp))
  # [1] 18948605
  
  print(nrow(old))
  # [1] 18947291        4
  
  # some sample IDs duplicated, but different summary stats
  dim(old_tmp[duplicated(old_tmp$rsid),])
  # [1] 1314    8
  
  old<-merge(old,old_tmp,by.x="ID_b37",by.y="rsid")
  dim(old)
  # [1] 18948605       11
  
  rm(old_tmp)
  
  old_no_new<-old[which(!old$X %in% new$X),]
  dim(old_no_new)
  # [1] 5685142      11
  
  old_in_new<-old[which(old$X %in% new$X),] 
  dim(old_in_new)
  # [1] 13263463       11
  
  5685142+13263463
  # [1] 18948605
  
  ###########################################
  ## systematic differences among those:
  
  # how many in new alternative contings:
  
  table(old_no_new$chr)
  table(old_in_new$chr)
  
  old_no_new$chr<-gsub(".*_KI.*","alt",old_no_new$chr)
  old_no_new$chr<-gsub(".*_GL.*","alt",old_no_new$chr)
  old_no_new$chr<-gsub("chrUn_KI270742v1","alt",old_no_new$chr)
  old_no_new$chr<-gsub("altalt","alt",old_no_new$chr)
  old_no_new$chr<-gsub("altrandom","alt",old_no_new$chr)
  table(old_no_new$chr)
  
  189+ 5942+12
  # [1] 6143
  
  old_no_new<-old_no_new[which(!old_no_new$chr %in% c("alt","chrX"))]

  summary(old$info)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 0.4000  0.6470  0.8196  0.7931  0.9683  1.0000 
  
  summary(old_no_new$info)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 0.4000  0.5234  0.6148  0.6437  0.7367  1.0000 
  summary(old_in_new$info)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 0.4000  0.7541  0.9156  0.8572  0.9815  1.0000
  
  summary(old$all_maf)
  # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
  # 0.001001 0.002351 0.013831 0.094001 0.144096 0.500000 
  
  summary(old_no_new$all_maf)
  # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
  # 0.001001 0.001312 0.002140 0.043284 0.007878 0.500000 
  summary(old_in_new$all_maf)
  # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
  # 0.001001 0.005223 0.037939 0.115729 0.197140 0.500000

  old$group<-NA
  old$group[which(old$ID_b37 %in% old_no_new$ID_b37)]<-"not_in_new"
  old$group[which(old$ID_b37 %in% old_in_new$ID_b37)]<-"in_new"
  table(old$group,useNA="ifany")
  #   in_new not_in_new       <NA> 
  # 13263463    5678999       6143
  
  
  # remove NAs - those are the ones in alt cont
  old<-old[which(!is.na(old$group)),]
  rm(old_in_new,old_no_new)
  
  p1<-ggplot(old, aes(x=log10(all_maf), color=group)) +
    geom_density()
  p2<-ggplot(old, aes(x=info, color=group)) +
    geom_density()
  
  p<-ggarrange(p1,p2,nrow=2,common.legend=T)
  
  # plot those 
  pdf(paste("~/tmp_plots/gwas3_previous_variants_in_notin_new_maf_info_distribution_",pheno[i],".pdf",sep=""),width=10, height=10)
  print(p)
  dev.off()
  rm(p1,p2)
  
  p1<-ggplot(old, aes(x=-log10(frequentist_add_pvalue), color=group)) +
    geom_histogram(position="identity", alpha=0.5)
  p2<-ggplot(old[which(old$frequentist_add_pvalue<=5e-8)], aes(x=-log10(frequentist_add_pvalue), color=group)) +
    geom_histogram(position="identity", alpha=0.5) + xlim(0,max(-log10(frequentist_add_pvalue)))
  p<-ggarrange(p1,p2,nrow=2,common.legend=T)
  
  pdf(paste("~/tmp_plots/gwas3_previous_variants_in_notin_new_pval_distribution_",pheno[i],".pdf",sep=""),width=10, height=10)
  print(p)
  dev.off()
  
  # there are gwas significant results in "not_in_new" group - are those at reported regions?
  
  dim(old[which((old$group=="not_in_new" & old$frequentist_add_pvalue<=5e-8)),])
  # [1] 4116   12
  dim(old[which((old$group=="in_new" & old$frequentist_add_pvalue<=5e-8)),])
  # [1] 3182   12
  
  table(old$chr[which((old$group=="not_in_new" & old$frequentist_add_pvalue<=5e-8))])
  # chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 
  # 57   142    86   110    15     7     2    79    65    35    18   138    29 
  # chr21 chr22  chr3  chr4  chr5  chr6  chr7  chr8  chr9 
  # 8     8    60    24    51  3105    26     7    44
  
  table(old$chr[which((old$group=="in_new" & old$frequentist_add_pvalue<=5e-8))])
  # chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 
  # 281   103    34    70     4    51    27    46    85    64    19   103     7 
  # chr21 chr22  chr3  chr4  chr5  chr6  chr7  chr8  chr9 
  # 136    65   137     7   415  1242    36     9   241 
  
  
  # plot info vs r2 for variants in both
   
  info<-merge(new[,c("X","ID","ALLELE0","ALLELE1","INFO")],old[,c("X","ID_b37","info","group")],by="X")
  
  p<-ggplot(info, aes(x=info, y=(INFO))) + xlim(0.3,1) + ylim(0.3,1) +
    geom_abline(intercept = 0, slope = 1, color="red",linetype="dashed", size=0.5) +
    geom_smooth(method='lm', color="darkgrey",linetype="dashed", size=0.8, se=T)
  
  pdf(paste("~/tmp_plots/gwas3_previous_variants_in_info_new_old_distribution.pdf",sep=""))
  print(p)
  dev.off()
  
  rm(new)
  
  # how many are indeed in topmed results but with lower info, maf than threshold:
  for (chr in c(1:22)) {
    tmp<-fread(paste(path,"imputed/all_hce/eur/chr",chr,".info.gz",sep=""),head=T)
    tmp$X<-gsub("chr","",tmp$SNP)
    tmp$X<-gsub(":[A-Z]*:[A-Z]*$","",tmp$X)
    tmp$X<-gsub(":","_",tmp$X)
    
    tmp<-subset(tmp, X %in% old$X)
    
    if(chr==1) {
      topmed<-tmp
    } else {
      topmed<-rbind(topmed,tmp)
    }
    rm(tmp)
  }
  
  dim(old)
  # [1] 18942462       12
  dim(topmed)
  # [1] 15506915       14
  
  info<-merge(topmed[,c("X","SNP","REF(0)","ALT(1)","Rsq","MAF")],old[,c("X","ID_b37","info","all_maf","group","frequentist_add_pvalue")],by="X")
  
  table(info$group,useNA="ifany")
  # in_new   not_in_new 
  # 13497863    2030920 

  info$ref<-gsub("[0-9]{1,2}:[0-9]*_","",info$ID_b37)
  info$alt<-gsub("[A-Z]*_","",info$ref)
  info$ref<-gsub("_[A-Z]*","",info$ref)
  

  # plot info vs r2 - same varians (position and alleles)
  
  info<-as.data.frame(info)
  info2<-info[which( (info$'REF(0)'==info$ref & info$'ALT(1)'==info$alt) |
                       (info$'REF(0)'==info$alt & info$'ALT(1)'==info$ref)),]
  
  
  p<-ggplot(info2, aes(y=Rsq, x=(info))) + xlim(0,1) + ylim(0,1) + geom_point(alpha = 1/10) +
    geom_abline(intercept = 0, slope = 1, color="red",linetype="dashed", size=0.5) +
    geom_smooth(method='lm', color="darkgrey",linetype="dashed", size=0.8, se=T)
  
  pdf(paste("~/tmp_plots/gwas3_previous_variants_in_info_rsq_new_old_TOPMED.pdf",sep=""))
  print(p)
  dev.off()
  
  
  
  
  p1<-ggplot(info, aes(x=log10(MAF), color=group)) +
    geom_density()
  p2<-ggplot(info, aes(x=Rsq, color=group)) +
    geom_density()
  
  p<-ggarrange(p1,p2,nrow=2,common.legend=T)
  
  # plot those 
  pdf(paste("~/tmp_plots/gwas3_previous_variants_in_notin_new_maf_info_distribution_",pheno[i],"_in_TOPMED.pdf",sep=""),width=10, height=10)
  print(p)
  dev.off()
  
  old<-as.data.frame(old)
  tmp<-old[which(!old$X %in% info$X),]
  dim(tmp)
  
  table(old$group[which(old$ID_b37 %in% tmp$ID_b37)])
  # not_in_new 
  # 3707542
  
  old$group[which(old$ID_b37 %in% tmp$ID_b37)]<-"not_in_new_not_topmed"
  table(old$group)
  # in_new              not_in_new not_in_new_not_topmed 
  # 13263463               1971457               3707542 
  
  # not useful plots
  # p1<-ggplot(old, aes(x=-log10(frequentist_add_pvalue), color=group)) +
  #   geom_histogram(position="identity", alpha=0.5)
  # p2<-ggplot(old[which(old$frequentist_add_pvalue<=5e-8),], aes(x=-log10(frequentist_add_pvalue), color=group)) +
  #   geom_histogram(position="identity", alpha=0.5) + xlim(0,max(-log10(old$frequentist_add_pvalue)))
  # p<-ggarrange(p1,p2,nrow=2,common.legend=T)
  # 
  # pdf(paste("~/tmp_plots/gwas3_previous_variants_in_notin_new_pval_distribution_",pheno[i],"_TOPMed.pdf",sep=""),width=10, height=10)
  # print(p)
  # dev.off()
  
  dim(old[which(old$group=="not_in_new_not_topmed" & old$frequentist_add_pvalue<=5E-8),])
  # [1] 3013   12
  dim(old[which(old$group=="not_in_new" & old$frequentist_add_pvalue<=5E-8),])
  # [1] 1103   12
  dim(old[which(old$group=="in_new" & old$frequentist_add_pvalue<=5E-8),])
  # [1] 3182   12
  
  table(old$chr[which(old$group=="not_in_new_not_topmed" & old$frequentist_add_pvalue<=5E-8)])
  # chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 
  # 35    42    54    49     3     7     1    27    22     4     6    30     6 
  # chr21 chr22  chr3  chr4  chr5  chr6  chr7  chr8  chr9 
  # 7     5    38     3    47  2584    16     5    22 
  table(old$chr[which(old$group=="not_in_new" & old$frequentist_add_pvalue<=5E-8)])
  table(old$chr[which(old$group=="in_new" & old$frequentist_add_pvalue<=5E-8)])
  # chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 
  # 281   103    34    70     4    51    27    46    85    64    19   103     7 
  # chr21 chr22  chr3  chr4  chr5  chr6  chr7  chr8  chr9 
  # 136    65   137     7   415  1242    36     9   241
  
  # MHC region: https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38.p13
  # which defines it as the 5mb region chr6:28510120-33480577 in GRCh38 coordinates
  
  old$allele1<-gsub("[0-9]{1,2}:[0-9]*_","",old$ID_b37)
  old$allele2<-gsub("[A-Z]*_","",old$allele1)
  old$allele1<-gsub("_[A-Z]*","",old$allele1)
  
  old$indel<-0
  old$indel[which(nchar(old$allele1)>1 | nchar(old$allele2)>1)]<-1
  
  table(old$indel)
  # 0        1 
  # 17100952  1841510 
  
  table(old$indel,old$group)
  #     in_new not_in_new not_in_new_not_topmed
  # 0 12449029    1810347               2841576
  # 1   814434     161110                865966
  
  dups<-old$X[duplicated(old$X)]
  old$duplicated_position<-0
  old$duplicated_position[which(old$X %in% dups)]<-1
  
  table(old$duplicated_position,old$group)
  #     in_new not_in_new not_in_new_not_topmed
  # 0 13236918    1964618               3665203
  # 1    26545       6839                 42339
  
  
  
  dim(old[which(old$group=="not_in_new_not_topmed" & old$frequentist_add_pvalue<=5E-8 & old$chr=="chr6" & old$position_b38>=28510120 & old$position_b38<=33480577),])
  # [1] 2454   12
  
  fwrite(all,paste("/path/to/ibdgwas/IIBDGC/summary_files/summary_stats_",pheno[i],"_gwas3_old_lifted_hg38_ed_plus_new.tsv.gz",sep=""),
         col.names=T,row.names=F,quote=F,sep="\t")
  
  
}

q("no")



### compare the variants in old GWAS with all TOPMed ones!

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=5000

# format IDs to match TOPMed IDs
bsub -J"TOPMed" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/niddk_cho_gsa/logs/stderr_TOPMed_subset \
-o ${path_gwas}pre_imputation/QC/niddk_cho_gsa/logs/stdout_TOPMed_subset \
"awk 'NR==FNR{vals[\$5];next} (\"chr\"\$10) in vals' ${path_gwas}summary_files/summary_stats_ibd_gwas3_old_all_variants_lifted_hg38_ed_plus_chr:position \
<(zcat ${path_gwas}resources/TOPMed/PASS.Variantsbravo-dbsnp-all_edited_3.vcf.gz) > ${path_gwas}summary_files/summary_stats_ibd_gwas3_old_all_variants_lifted_hg38_ed_plus_chr:position_TOPMed_variants"


wc -l ${path_gwas}summary_files/summary_stats_ibd_gwas3_old_all_variants_lifted_hg38_ed_plus_chr:position_TOPMed_variants
# 17245005 /path/to/ibdgwas/IIBDGC/summary_files/summary_stats_ibd_gwas3_old_all_variants_lifted_hg38_ed_plus_chr:position_TOPMed_variants
wc -l ${path_gwas}summary_files/summary_stats_ibd_gwas3_old_all_variants_lifted_hg38_ed_plus_chr:position
# 18948607 /path/to/ibdgwas/IIBDGC/summary_files/summary_stats_ibd_gwas3_old_all_variants_lifted_hg38_ed_plus_chr:position



MEM=18000
bsub -Is -m "modern_hardware" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

pheno<-c("cd","uc","ibd")

i=3

print(pheno[i])
  
old<-fread(paste("/path/to/ibdgwas/IIBDGC/summary_files/summary_stats_",pheno[i],"_gwas3_old_all_variants_lifted_hg38_ed.gz",sep=""),head=F)
colnames(old)<-c("chr","position_b38","ID_b37","pvalue")
old<-old[,c("chr","position_b38","ID_b37")]
old<-old[!duplicated(old),]
old$X<-paste(gsub("chr","",old$chr),old$position_b38,sep="_")
old$ID_b37<-gsub("chr","",old$ID_b37)

print(nrow(old))
# [1] 18947291


# old GWAS3 summary stats
for (chr in c(1:22)) {
  
  tmp<-fread(paste("/path/to/project",pheno[i],"/",chr,".assoc",sep=""),head=T)
  tmp<-tmp[,c("rsid","info","all_maf","cases_total","controls_total","frequentist_add_beta_1","frequentist_add_se_1","frequentist_add_pvalue")]
  
  tmp<-subset(tmp, rsid %in% old$ID_b37)
  
  if(chr==1) {
    old_tmp<-tmp
  } else {
    old_tmp<-rbind(old_tmp,tmp)
  }
  rm(tmp)
}

print(nrow(old_tmp))
# [1] 18948607        8
old_tmp<-old_tmp[!duplicated(old_tmp),]
print(nrow(old_tmp))
# [1] 18948605

print(nrow(old))
# [1] 18947291        4

# some sample IDs duplicated, but different summary stats
dim(old_tmp[duplicated(old_tmp$rsid),])
# [1] 1314    8

old<-merge(old,old_tmp,by.x="ID_b37",by.y="rsid")
dim(old)
# [1] 18948605       11

rm(old_tmp)

old$chr<-gsub(".*_KI.*","alt",old$chr)
old$chr<-gsub(".*_GL.*","alt",old$chr)
old$chr<-gsub("chrUn_KI270742v1","alt",old$chr)
old$chr<-gsub("altalt","alt",old$chr)
old$chr<-gsub("altrandom","alt",old$chr)
table(old$chr)

189+ 5942+12
# [1] 6143

old<-old[which(!old$chr %in% c("alt","chrX"))]
  
print(nrow(old))
# [1] 18947291

old_topmed<-fread("/path/to/ibdgwas/IIBDGC/summary_files/summary_stats_ibd_gwas3_old_all_variants_lifted_hg38_ed_plus_chr:position_TOPMed_variants",head=F)
old_topmed<-old_topmed[!duplicated(old_topmed),]
print(nrow(old_topmed))
# [1] 17245005

old_topmed$V10<-gsub(":","_",old_topmed$V10)

old$class<-"variants_not_available_in_topmed"
old$class[which(old$X %in% old_topmed$V10)]<-"variants_available_in_topmed"
table(old$class)
# variants_available_in_topmed variants_not_available_in_topmed 
#                     15782630                          3159832 


# new results
for (chr in c(1:22)) {
  
  tmp<-fread(paste(path,"post_imputation/analysis/regenie/humancoreexome/",pheno[i],"/chr",chr,"_humancoreexome_step2_",pheno[i],"_eur_sex_PCs_firthse.regenie",sep=""),head=T)
  tmp$X<-paste(tmp$CHROM,tmp$GENPOS,sep="_")
  
  # keep all new variants
  # tmp<-subset(tmp, X %in% old$X)
  
  if(chr==1) {
    new<-tmp
  } else {
    new<-rbind(new,tmp)
  }
  rm(tmp)
}

print(nrow(new))
# [1] 14056619

table(old$class[which(old$X %in% new$X)])
# variants_available_in_topmed variants_not_available_in_topmed 
#                     12776714                           486749 

table(old$class[which(!old$X %in% new$X),])
# variants_available_in_topmed variants_not_available_in_topmed 
#                      3005916                          2673083


head(old[which(old$X %in% new$X & old$class=="variants_not_available_in_topmed"),])
#                   ID_b37   chr position_b38           X     info    all_maf
# 1:      10:100000625_A_G chr10     98240868 10_98240868 0.997971 0.42965300
# 2:      10:100019542_G_A chr10     98259785 10_98259785 0.763040 0.00105557

new[which(new$X=="10_98240868"),]
# CHROM   GENPOS                 ID ALLELE0 ALLELE1   A1FREQ     INFO TEST
# 1:    10 98240868 chr10:98240868:A:G       A       G 0.425181 0.973252  ADD
# BETA.Y1     SE.Y1 CHISQ.Y1 LOG10P.Y1           X
# 1: 0.012841 0.0217602 0.348237  0.255619 10_98240868

old_topmed[which(old_topmed$V10=="10_98240868"),]
# Empty data.table (0 rows and 10 cols): V1,V2,V3,V4,V5,V6... # check out original files as well later


old$class[which(old$X %in% new$X)]<-"variants_available_in_new_results_filtered_maf_info"
table(old$class)
# variants_available_in_new_results_filtered_maf_info 
# 13263463 
# variants_available_in_topmed 
# 3005916 
# variants_not_available_in_topmed 
# 2673083 

fwrite(old,"/path/to/ibdgwas/IIBDGC/summary_files/summary_stats_ibd_gwas3_old_all_variants_lifted_hg38_ed_plus_chr:position_topmed_availability.tsv.gz",
       col.names=T,row.names=F,quote=F,sep="\t")

rm(old_topmed)

q()

#############################

MEM=18000
bsub -Is -m "modern_hardware" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

pheno<-c("cd","uc","ibd")

i=3

print(pheno[i])

old<-fread("/path/to/ibdgwas/IIBDGC/summary_files/summary_stats_ibd_gwas3_old_all_variants_lifted_hg38_ed_plus_chr:position_topmed_availability.tsv.gz",head=T)
old$tmp<-NA

table(old$class)
# variants_available_in_new_results_filtered_maf_info 
# 13263463 
# variants_available_in_topmed 
# 1141924 
# variants_available_in_topmed_and_in_imputed_files 
# 1971457 
# variants_not_available_in_topmed 
# 2565618 

# how many are indeed in topmed results but with lower info, maf than threshold:
for (chr in c(1:22)) {
  tmp<-fread(paste(path,"imputed/all_hce/eur/chr",chr,".info.gz",sep=""),head=T)
  tmp$X<-gsub("chr","",tmp$SNP)
  tmp$X<-gsub(":[A-Z]*:[A-Z]*$","",tmp$X)
  tmp$X<-gsub(":","_",tmp$X)
  
  old$class[which( (old$X %in% tmp$X) & 
                     (old$class!="variants_available_in_new_results_filtered_maf_info"))]<-"variants_available_in_topmed_and_in_imputed_files"
  rm(tmp)
}


# indel distribution:

old$allele1<-gsub("[0-9]{1,2}:[0-9]*_","",old$ID_b37)
old$allele2<-gsub("[A-Z]*_","",old$allele1)
old$allele1<-gsub("_[A-Z]*","",old$allele1)

old$indel<-0
old$indel[which(nchar(old$allele1)>1 | nchar(old$allele2)>1)]<-1

table(old$indel)
# 0        1 
# 17100952  1841510 

table(old$class,old$indel)


# how many new but not in old paper!




