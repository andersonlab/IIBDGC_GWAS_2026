# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=2000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R \

library(data.table)
library(ggplot2)
library(ggpubr)
library(colorspace)
library(rcartocolor)
library(qvalue)

path_gwas="/path/to/ibdgwas/IIBDGC/"


threshold<-0.05/((103+163+129+22)*3)
threshold
# [1] 3.996803e-05

col<-carto_pal(12, "Safe")
col_ligth_1<-lighten(col, 0.90, space = "HCL")
col_ligth_2<-lighten(col, 0.70, space = "HCL")
col_ligth_3<-lighten(col, 0.50, space = "HCL")
col_ligth_4<-lighten(col, 0.30, space = "HCL")

colo<-c(col[0:11],col_ligth_1[0:11],col_ligth_2[0:11],col_ligth_3[0:11],col_ligth_4[0:11])


files<-list.files(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/output/",sep=""))
files<-files[grep(".results",files)]

# look for the annotation for endotelial cells:

enteroc<-files[grep("entero|Enteroc|Enc|Paneth|Gbl|Goblet|Turf",files)]
enteroc<-enteroc[grep("ibd",enteroc)]
enteroc<-gsub("\\.results","",enteroc)
enteroc<-gsub("ibd\\.baseline\\.","",enteroc)
enteroc<-gsub("\\.files_version_june2024","",enteroc)

list_snps<-c("chr18:63286830:T:C","chr18:63285183:C:T","chr18:63287344:A:T","chr2:111175424:G:T","chr2:111165509:G:GT")
# chr2:111165509:G:GT - LD 0.99
# chr18:63285183:C:T - LD 0.92
#  chr18:63287344:A:T - LD 0.8

chr<-gsub(":.*","",list_snps)
chr<-chr[!duplicated(chr)]

pos<-gsub("chr[0-9]{1,2}:","",list_snps)
pos<-gsub(":.*","",pos)
pos<-as.numeric(pos)

rm(tmp1,tmp,tmp_all)
for (i in c(1:length(enteroc))) {

    for (ii in 1:length(chr)) {

         tmp<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/ldscore/analysis_june2024/no_baseline_files/",enteroc[i],"/",enteroc[i],".sorted.june2024.",chr[ii],".gz",sep=""),head=T)    
        #  tmp<-tmp[which(tmp$BP>=pos[ii]-2000 & tmp$BP<=pos[ii]+2000),]
        tmp<-tmp[which( (tmp$SNP %in% list_snps) | (tmp$BP>=pos[ii]-2000 & tmp$BP<=pos[ii]+2000)),]

        if( ii==1) {
            tmp1<-tmp
        }else{
            tmp1<-rbind(tmp1,tmp)
        }
    }

    if(i==1) {
        tmp_all<-tmp1
    } else {
        tmp_all<-merge(tmp_all,tmp1[,c(1,7)],by="SNP",all.x=T)
    }

}

# tmp_all[which(tmp_all$SNP %in% list_snps),]
#                    SNP    CHR        BP     A1     A2    ...
#                 <char> <char>     <int> <char> <char> <char>
# 1:  chr18:63285183:C:T  chr18  63285183      T      C    ...
# 2:  chr18:63286830:T:C  chr18  63286830      C      T    ...
# 3:  chr18:63287344:A:T  chr18  63287344      T      A    ...
# 4: chr2:111165509:G:GT   chr2 111165509     GT      G    ...
# 5:  chr2:111175424:G:T   chr2 111175424      T      G    ...
#    cCRE_Accessibility_Enc.1 cCRE_Accessibility_Enc.2 cCRE_Accessibility_Enc.3
#                       <num>                    <num>                    <num>
# 1:                0.0000000                0.0000000               0.00000000
# 2:                0.0000000                0.0000000               0.00000000
# 3:                0.2918785                0.2112494               0.01597996
# 4:                0.2250361                0.1419873               0.01597996
# 5:                0.0000000                0.0000000               0.00000000
#    cCRE_Accessibility_Gbl.1 cCRE_Accessibility_Gbl.2
#                       <num>                    <num>
# 1:                0.0000000               0.00000000
# 2:                0.0000000               0.00000000
# 3:                0.1110372               0.09572079
# 4:                0.3129231               0.87889089
# 5:                0.0000000               0.00000000
#    duodenum_epithelial_Best4_Enterocytes_hickey
#                                           <int>
# 1:                                            0
# 2:                                            0
# 3:                                            0
# 4:                                            1
# 5:                                            0
#    duodenum_epithelial_Enterochromaffin_hickey
#                                          <int>
# 1:                                           0
# 2:                                           0
# 3:                                           0
# 4:                                           0
# 5:                                           0
#    duodenum_epithelial_Enterocytes_hickey duodenum_epithelial_Goblet_hickey
#                                     <int>                             <int>
# 1:                                      0                                 0
# 2:                                      0                                 0
# 3:                                      0                                 0
# 4:                                      1                                 1
# 5:                                      0                                 0
#    duodenum_epithelial_Immature_Enterocytes_hickey
#                                              <int>
# 1:                                               0
# 2:                                               0
# 3:                                               0
# 4:                                               1
# 5:                                               0
#    duodenum_epithelial_Immature_Goblet_hickey duodenum_epithelial_Paneth_hickey
#                                         <int>                             <int>
# 1:                                          0                                 0
# 2:                                          0                                 0
# 3:                                          0                                 0
# 4:                                          1                                 0
# 5:                                          0                                 0
#    ileum_epithelial_Best4_Enterocytes_hickey
#                                        <int>
# 1:                                         0
# 2:                                         0
# 3:                                         0
# 4:                                         0
# 5:                                         0
#    ileum_epithelial_Enterochromaffin_hickey ileum_epithelial_Enterocytes_hickey
#                                       <int>                               <int>
# 1:                                        0                                   0
# 2:                                        0                                   0
# 3:                                        0                                   0
# 4:                                        1                                   1
# 5:                                        0                                   0
#    ileum_epithelial_Goblet_hickey ileum_epithelial_Immature_Enterocytes_hickey
#                             <int>                                        <int>
# 1:                              0                                            0
# 2:                              0                                            0
# 3:                              0                                            0
# 4:                              1                                            1
# 5:                              0                                            0
#    ileum_epithelial_Immature_Goblet_hickey ileum_epithelial_Paneth_hickey
#                                      <int>                          <int>
# 1:                                       0                              0
# 2:                                       0                              0
# 3:                                       0                              0
# 4:                                       0                              1
# 5:                                       0                              0
#    jejunum_epithelial_Best4_Enterocytes_hickey
#                                          <int>
# 1:                                           0
# 2:                                           0
# 3:                                           0
# 4:                                           1
# 5:                                           0
#    jejunum_epithelial_Enterochromaffin_hickey
#                                         <int>
# 1:                                          0
# 2:                                          0
# 3:                                          0
# 4:                                          0
# 5:                                          0
#    jejunum_epithelial_Enterocytes_hickey jejunum_epithelial_Goblet_hickey
#                                    <int>                            <int>
# 1:                                     0                                0
# 2:                                     0                                0
# 3:                                     0                                0
# 4:                                     1                                1
# 5:                                     0                                0
#    jejunum_epithelial_Immature_Enterocytes_hickey
#                                             <int>
# 1:                                              0
# 2:                                              0
# 3:                                              0
# 4:                                              1
# 5:                                              0
#    jejunum_epithelial_Immature_Goblet_hickey jejunum_epithelial_Paneth_hickey
#                                        <int>                            <int>
# 1:                                         0                                0
# 2:                                         0                                0
# 3:                                         0                                0
# 4:                                         1                                1
# 5:                                         0                                0
#    non_multiome_duodenum_epithelial_Best4_Enterocytes_hickey
#                                                        <int>
# 1:                                                         0
# 2:                                                         0
# 3:                                                         0
# 4:                                                         0
# 5:                                                         0
#    non_multiome_duodenum_epithelial_Enterocytes_hickey
#                                                  <int>
# 1:                                                   0
# 2:                                                   0
# 3:                                                   0
# 4:                                                   1
# 5:                                                   0
#    non_multiome_duodenum_epithelial_Goblet_hickey
#                                             <int>
# 1:                                              0
# 2:                                              0
# 3:                                              0
# 4:                                              0
# 5:                                              0
#    non_multiome_duodenum_epithelial_Immature_Enterocytes_hickey
#                                                           <int>
# 1:                                                            0
# 2:                                                            0
# 3:                                                            0
# 4:                                                            1
# 5:                                                            0
#    non_multiome_duodenum_epithelial_Immature_Goblet_hickey
#                                                      <int>
# 1:                                                       0
# 2:                                                       0
# 3:                                                       0
# 4:                                                       1
# 5:                                                       0
#    non_multiome_duodenum_epithelial_Paneth_hickey
#                                             <int>
# 1:                                              0
# 2:                                              0
# 3:                                              0
# 4:                                              0
# 5:                                              0
#    non_multiome_ileum_epithelial_Best4_Enterocytes_hickey
#                                                     <int>
# 1:                                                      0
# 2:                                                      0
# 3:                                                      0
# 4:                                                      0
# 5:                                                      0
#    non_multiome_ileum_epithelial_Enterocytes_hickey
#                                               <int>
# 1:                                                0
# 2:                                                0
# 3:                                                0
# 4:                                                1
# 5:                                                0
#    non_multiome_ileum_epithelial_Goblet_hickey
#                                          <int>
# 1:                                           0
# 2:                                           0
# 3:                                           0
# 4:                                           1
# 5:                                           0
#    non_multiome_ileum_epithelial_Immature_Enterocytes_hickey
#                                                        <int>
# 1:                                                         0
# 2:                                                         0
# 3:                                                         0
# 4:                                                         1
# 5:                                                         0
#    non_multiome_ileum_epithelial_Immature_Goblet_hickey
#                                                   <int>
# 1:                                                    0
# 2:                                                    0
# 3:                                                    0
# 4:                                                    1
# 5:                                                    0
#    non_multiome_ileum_epithelial_Paneth_hickey
#                                          <int>
# 1:                                           0
# 2:                                           0
# 3:                                           0
# 4:                                           0
# 5:                                           0
#    non_multiome_jejunum_epithelial_Best4_Enterocytes_hickey
#                                                       <int>
# 1:                                                        0
# 2:                                                        0
# 3:                                                        0
# 4:                                                        1
# 5:                                                        0
#    non_multiome_jejunum_epithelial_Enterocytes_hickey
#                                                 <int>
# 1:                                                  0
# 2:                                                  0
# 3:                                                  0
# 4:                                                  1
# 5:                                                  0
#    non_multiome_jejunum_epithelial_Goblet_hickey
#                                            <int>
# 1:                                             0
# 2:                                             0
# 3:                                             0
# 4:                                             1
# 5:                                             0
#    non_multiome_jejunum_epithelial_Immature_Enterocytes_hickey
#                                                          <int>
# 1:                                                           0
# 2:                                                           0
# 3:                                                           0
# 4:                                                           1
# 5:                                                           0
#    non_multiome_jejunum_epithelial_Immature_Goblet_hickey
#                                                     <int>
# 1:                                                      0
# 2:                                                      0
# 3:                                                      0
# 4:                                                      0
# 5:                                                      0
#    non_multiome_jejunum_epithelial_Paneth_hickey
#                                            <int>
# 1:                                             0
# 2:                                             0
# 3:                                             0
# 4:                                             0
# 5:                                             0