# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
bj# RUN PCA IN EUROPEAN ANCESTRY SAMPLES, ESTIMATE PCS IN NON DUPLICATED EUR SUBSET, AND PROJECT ALL

# TO DEFINE POPULATION ANCESTRY PCS TO INCLUDE IN THE ANALYSIS (ARRAY BY ARRAY)

###########################################################################################################################################################################################

########################################################
# 0.0 DEFINE A SET OF NON-RELATED SAMPLES - 3rd degree #
########################################################

# START WITH EUR POST-QC NON PRUNED DATA, see script step_25

array=(illumina370 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 gsa illuminaexome)
path_gwas=/path/to/ibdgwas/IIBDGC/


for i in ${array[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}post_imputation/2022/${i}/genotyped_data/${i}_all_studies_merged_eur \
--keep ${path_gwas}pre_imputation/QC/pca_1000gp/list_eur_ancestry_samples_noDuplicates_perstudy_interstudy_2 \
--keep-allele-order \
--autosome \
--maf 0.001 \
--make-bed --out ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_noDuplicates
done


wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/*_all_studies_merged_eur_noDuplicates.fam
# 4648 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix500_all_studies_merged_eur_noDuplicates.fam
# 8246 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix6_all_studies_merged_eur_noDuplicates.fam
# 59509 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_noDuplicates.fam
# 21243 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humancoreexome_all_studies_merged_eur_noDuplicates.fam
# 2053 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomni1_all_studies_merged_eur_noDuplicates.fam
# 1238 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomniexpress_all_studies_merged_eur_noDuplicates.fam
# 3452 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina370_all_studies_merged_eur_noDuplicates.fam
# 1782 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illuminaexome_all_studies_merged_eur_noDuplicates.fam
# 3396 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/quad610_all_studies_merged_eur_noDuplicates.fam
# 105567 total

wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/*_all_studies_merged_eur_noDuplicates.bim
# 398750 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix500_all_studies_merged_eur_noDuplicates.bim
# 604280 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix6_all_studies_merged_eur_noDuplicates.bim
# 310226 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_noDuplicates.bim
# 230857 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humancoreexome_all_studies_merged_eur_noDuplicates.bim
# 785998 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomni1_all_studies_merged_eur_noDuplicates.bim
# 632174 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomniexpress_all_studies_merged_eur_noDuplicates.bim
# 261044 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina370_all_studies_merged_eur_noDuplicates.bim
# 221591 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illuminaexome_all_studies_merged_eur_noDuplicates.bim
# 547964 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/quad610_all_studies_merged_eur_noDuplicates.bim
 
# The amount of computer memory required by KING analysis is modest, at ~N ✕ M / 4 (where N is the number of samples and M is the number of SNPs)
# plus a small percentage of overhead cost. E.g., for a dataset consisting of 100,000 samples each genotyped at 1,000,000 SNPs, 
# the required memory size is ~25GB.

# Max estimated with GSA
# (59590*310226/4)/1E+6
# [1] 4621.592




# accuracy up to 3rd- or 4th-degree (depending on array or WGS) for --related and --ibdseg analyses, and up to 2nd-degree for --kinship analysis

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=6500
array=(illumina370 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 gsa illuminaexome)

for i in ${array[@]}
do
bsub -J"kg" -M"$MEM" -n8 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal -n 8 \
-e ${path_gwas}pre_imputation/QC/log/stderr_king_per_array_${i} \
-o ${path_gwas}pre_imputation/QC/log/stdout_king_per_array_${i} \
"/path/to/software/./king \
-b ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_noDuplicates.bed \
--related --degree 3 --cpus 8 --prefix ${path_gwas}pre_imputation/QC/pca_1000gp/king_relatedness_eur_noDuplicates_${i}"
done

for i in ${array[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/log/stdout_king_per_array_${i} | grep "Successfully"
done


#### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(viridis)
library(ggsci)

path<-"/path/to/ibdgwas/IIBDGC/"
array<-c("illumina370","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")

for (i in c(1:length(array)) ) {
  tmp<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/king_relatedness_eur_noDuplicates_",array[i],".kin0",sep=""),head=T)
  tmp$array<-array[i]
  if(i==1){
    all<-tmp
  }else{
    all<-rbind(all,tmp)
  }
}

table(all$array,all$InfType)
#                    2nd   3rd    FS    UN    PO   4th
# affymetrix500        6     9     0     0     0     0
# affymetrix6          4     8     7     0    56     0
# gsa                539   500   795 64516  1388     6
# humancoreexome      93   112   163  8409   876     5
# humanomni1           2    14    16     0     2     0
# humanomniexpress     1     1     1     0     2     0
# illumina370          1     5     2    17     0     0
# illuminaexome        0     7    12     0    15     0
# quad610              1    10     0     1     0     0


# plot Kinship vs PropIBD
p1<-ggplot(all, aes(x = IBS0, y=Kinship, color = array)) +
  geom_point() 
p2<-ggplot(all, aes(x = PropIBD, y=Kinship, color = array)) +
  geom_point() 

p3<-ggplot(all, aes(x = IBS0, y=Kinship, color = InfType)) +
  geom_point() 
p4<-ggplot(all, aes(x = PropIBD, y=Kinship, color = InfType)) +
  geom_point() 

p5<-ggplot(all, aes(x = HetHet, y=PropIBD, color = array)) +
  geom_point() 
p6<-ggplot(all, aes(x = IBS0, y=PropIBD, color = array)) +
  geom_point() 

p7<-ggplot(all, aes(y = IBD1Seg, x=IBD2Seg, color = InfType)) +
  geom_point() 
p8<-ggplot(all, aes(y = IBD2Seg, x=PropIBD, color = array)) +
  geom_point() 

p12<-ggarrange(p1,p2,ncol=2,common.legend = T,legend="bottom")
p34<-ggarrange(p3,p4,ncol=2,common.legend = T,legend="bottom")
p56<-ggarrange(p5,p6,ncol=2,common.legend = T,legend="bottom")
p78<-ggarrange(p7,p8,ncol=2,legend="bottom")

pdf(paste(path,"pre_imputation/QC/pca_1000gp/kinship_vs_propIBD_eur_nonDuplicates.pdf",sep=""),width = 14, height = 28)
ggarrange(p12,p34,p56,p78,nrow=4,common.legend = T,legend="bottom")
dev.off()


pdf(paste(path,"pre_imputation/QC/pca_1000gp/kinship_histogram_eur_nonDuplicates.pdf",sep=""),width = 14, height = 7)
ggplot(all, aes(x=Kinship, color=array,fill=array)) +
  geom_density(alpha=0.4)
dev.off()


rel_ids<-c(as.character(all$FID1),as.character(all$FID2))
length(rel_ids)
table(length(rel_ids)==(nrow(all)*2))
rel_ids<-rel_ids[!duplicated(rel_ids)]
length(rel_ids)
# [1] 53897

rm(data_remove)
for (i in 1:length(rel_ids)) {
  
  tmp1<-all[which(all$FID1==rel_ids[i]),]
  tmp2<-all[which(all$FID2==rel_ids[i]),]
  colnames(tmp2)[1:2]<-colnames(tmp2)[2:1]
  
  tmp<-rbind(tmp1,tmp2)
  
  ids_tmp<-c(as.character(tmp$FID1),as.character(tmp$FID2))
  ids_tmp<-ids_tmp[!duplicated(ids_tmp)]
  
  #keep same order as in rel_ids and in all
  ids_tmp<-ids_tmp[match(rel_ids[which(rel_ids %in% ids_tmp)],ids_tmp)]
  
  data<-as.data.frame(matrix(ncol=1,nrow=length(ids_tmp)))
  data$V1<-ids_tmp
  
  keep_sample<-data$V1[i]
  if(!exists("data_remove")) {
    data_remove<-data[which(!data$V1 %in% keep_sample),,drop=F]
  } else {
    data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
  }
  
}

dim(data_remove)
# [1] 53994     1

data_remove<-data_remove[!duplicated(data_remove$V1),,drop=F]
dim(data_remove)
# [1] 20159     1

for (i in c(1:length(array)) ) {
  tmp<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/",array[i],"_all_studies_merged_eur_noDuplicates.fam",sep=""),head=F)
  tmp$array<-array[i]
  if(i==1){
    all<-tmp
  }else{
    all<-rbind(all,tmp)
  }
}

data_remove<-merge(data_remove,all[,c(1,7)],by.x="V1",by.y="V1",all.x=T)
table(data_remove$array)
# affymetrix500      affymetrix6              gsa   humancoreexome 
#            15               72            16165             3822 
# humanomni1 humanomniexpress      illumina370    illuminaexome 
#         32                5               17               23 
# quad610 
#      12 


data_remove<-data_remove[,c(1,1)]
colnames(data_remove)<-c("FID","IID")
write.table(data_remove,paste(path,"pre_imputation/QC/pca_1000gp/list_eur_related_3rd_degree_samples",sep=""),col.names=T,row.names=F,sep="\t",quote=F)


#######################################
# 1.0 liftover high-LD regions to b37 #
#######################################

# array=(illumina370 illumina550 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 gsa illuminaexome)
# path_gwas=/path/to/ibdgwas/IIBDGC/
# path_gwas=/path/to/ibdgwas/IIBDGC/
# 
# awk 'BEGIN{OFS="\t"}$1="chr"$1' /nfs/team152/carl/protocols/QC-paper/high-LD-regions.txt > ${path_gwas}pre_imputation/QC/pca_1000gp/high_LD_regions_b37.txt
# 
# ${path_gwas}previous_qced_b38/liftover/liftOver \
# ${path_gwas}pre_imputation/QC/pca_1000gp/high_LD_regions_b37.txt \
# ${path_gwas}previous_qced_b38/liftover/hg19ToHg38.over.chain.gz \
# ${path_gwas}pre_imputation/QC/pca_1000gp/high_LD_regions_lifted_hg38 \
# ${path_gwas}pre_imputation/QC/pca_1000gp/high_LD_regions_lifted_no_lifted_hg38

# some of the regions are split in b38, some are partially deleted, thus use resoure from:
https://github.com/meyer-lab-cshl/plinkQC/blob/master/inst/extdata/high-LD-regions-hg38-GRCh38.txt

path_gwas=/path/to/ibdgwas/IIBDGC/
cd ${path_gwas}pre_imputation/QC/pca_1000gp/
emacs high-LD-regions-hg38-GRCh38.txt

# cp  ${path_gwas}pre_imputation/QC/pca_1000gp/high-LD-regions-hg38-GRCh38.txt /path/to/ibdgwas/IIBDGC/resources/long_ld_regions/
# Update MHC based on https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38.p11
# MHC -- chr6 (NC_000006.12):28,510,120-33,480,577
# 28510120-500000 
# 33480577+500000
# 
# > 28510120-500000 
# [1] 28010120
# > 33480577+500000
# [1] 33980577

emacs /path/to/ibdgwas/IIBDGC/resources/long_ld_regions/high-LD-regions-hg38-GRCh38.txt

# previous PC runs in EUR sample highlighted two additional regions:
# chr2 134537609  137013709 44
# chr15 28083619 28268218 45
# chr17 15031076  16561983 46
# chr11 109625482 112115767 47
# chr11 46945864  57200532  48




###########################################################################################################################################################################################

################################################################
# 2.0 EXCLUDE VARIANTS IN HIGH LD REGIONS, KEEP ONLY AUTOSOMAL #
################################################################

# source files, same as used in regenie stage1
array=(illumina370 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 gsa illuminaexome)

for i in ${array[@]}
do
wc -l ${path_gwas}post_imputation/2022/${i}/genotyped_data/${i}_all_studies_merged_eur_pruned.bim
done
# 191680 /path/to/ibdgwas/IIBDGC/post_imputation/2022/illumina370/genotyped_data/illumina370_all_studies_merged_eur_pruned.bim
# 264863 /path/to/ibdgwas/IIBDGC/post_imputation/2022/affymetrix6/genotyped_data/affymetrix6_all_studies_merged_eur_pruned.bim
# 312061 /path/to/ibdgwas/IIBDGC/post_imputation/2022/humanomniexpress/genotyped_data/humanomniexpress_all_studies_merged_eur_pruned.bim
# 192361 /path/to/ibdgwas/IIBDGC/post_imputation/2022/affymetrix500/genotyped_data/affymetrix500_all_studies_merged_eur_pruned.bim
# 161991 /path/to/ibdgwas/IIBDGC/post_imputation/2022/humancoreexome/genotyped_data/humancoreexome_all_studies_merged_eur_pruned.bim
# 352258 /path/to/ibdgwas/IIBDGC/post_imputation/2022/humanomni1/genotyped_data/humanomni1_all_studies_merged_eur_pruned.bim
# 298845 /path/to/ibdgwas/IIBDGC/post_imputation/2022/quad610/genotyped_data/quad610_all_studies_merged_eur_pruned.bim
# 263372 /path/to/ibdgwas/IIBDGC/post_imputation/2022/gsa/genotyped_data/gsa_all_studies_merged_eur_pruned.bim
# 166340 /path/to/ibdgwas/IIBDGC/post_imputation/2022/illuminaexome/genotyped_data/illuminaexome_all_studies_merged_eur_pruned.bim

for i in ${array[@]}
do
wc -l ${path_gwas}post_imputation/2022/${i}/genotyped_data/${i}_all_studies_merged_eur_pruned.fam
done
# 5633 /path/to/ibdgwas/IIBDGC/post_imputation/2022/illumina370/genotyped_data/illumina370_all_studies_merged_eur_pruned.fam
# 11087 /path/to/ibdgwas/IIBDGC/post_imputation/2022/affymetrix6/genotyped_data/affymetrix6_all_studies_merged_eur_pruned.fam
# 1245 /path/to/ibdgwas/IIBDGC/post_imputation/2022/humanomniexpress/genotyped_data/humanomniexpress_all_studies_merged_eur_pruned.fam
# 4652 /path/to/ibdgwas/IIBDGC/post_imputation/2022/affymetrix500/genotyped_data/affymetrix500_all_studies_merged_eur_pruned.fam
# 22588 /path/to/ibdgwas/IIBDGC/post_imputation/2022/humancoreexome/genotyped_data/humancoreexome_all_studies_merged_eur_pruned.fam
# 2701 /path/to/ibdgwas/IIBDGC/post_imputation/2022/humanomni1/genotyped_data/humanomni1_all_studies_merged_eur_pruned.fam
# 3396 /path/to/ibdgwas/IIBDGC/post_imputation/2022/quad610/genotyped_data/quad610_all_studies_merged_eur_pruned.fam
# 81334 /path/to/ibdgwas/IIBDGC/post_imputation/2022/gsa/genotyped_data/gsa_all_studies_merged_eur_pruned.fam
# 2997 /path/to/ibdgwas/IIBDGC/post_imputation/2022/illuminaexome/genotyped_data/illuminaexome_all_studies_merged_eur_pruned.fam


for i in ${array[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}post_imputation/2022/${i}/genotyped_data/${i}_all_studies_merged_eur_pruned \
--autosome \
--allow-no-sex \
--exclude range /path/to/ibdgwas/IIBDGC/resources/long_ld_regions/high-LD-regions-hg38-GRCh38.txt \
--make-bed --out ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_noHighLD
done


###########################################################################################################################################################################################

######################################
# 3.0 EXCLUDE ANY ASSOCIATED VARIANT #
######################################

for i in ${array[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_noHighLD \
--assoc --out ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_noHighLD
done

for i in ${array[@]}
do
awk '{if($9 < 1E-4)print $2}' ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_noHighLD.assoc > \
${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_list_associated_variants_toexclude
done

for i in ${array[@]}
do
wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_list_associated_variants_toexclude
done
# 391 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina370_all_studies_merged_eur_list_associated_variants_toexclude
# 179 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix6_all_studies_merged_eur_list_associated_variants_toexclude
# 43 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomniexpress_all_studies_merged_eur_list_associated_variants_toexclude
# 90 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix500_all_studies_merged_eur_list_associated_variants_toexclude
# 363 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humancoreexome_all_studies_merged_eur_list_associated_variants_toexclude
# 75 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomni1_all_studies_merged_eur_list_associated_variants_toexclude
# 58 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/quad610_all_studies_merged_eur_list_associated_variants_toexclude
# 31890 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_list_associated_variants_toexclude
# 840 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illuminaexome_all_studies_merged_eur_list_associated_variants_toexclude

####

for i in ${array[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_noHighLD \
--exclude ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_list_associated_variants_toexclude \
--make-bed --out ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_noHighLD_noassoclogP4
done

for i in ${array[@]}
do
wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_noHighLD_noassoclogP4.bim
done
# 187994 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina370_all_studies_merged_eur_pruned_noHighLD_noassoclogP4.bim
# 252310 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix6_all_studies_merged_eur_pruned_noHighLD_noassoclogP4.bim
# 298206 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomniexpress_all_studies_merged_eur_pruned_noHighLD_noassoclogP4.bim
# 184769 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix500_all_studies_merged_eur_pruned_noHighLD_noassoclogP4.bim
# 154730 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humancoreexome_all_studies_merged_eur_pruned_noHighLD_noassoclogP4.bim
# 336097 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomni1_all_studies_merged_eur_pruned_noHighLD_noassoclogP4.bim
# 286020 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/quad610_all_studies_merged_eur_pruned_noHighLD_noassoclogP4.bim
# 221691 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_pruned_noHighLD_noassoclogP4.bim
# 158536 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illuminaexome_all_studies_merged_eur_pruned_noHighLD_noassoclogP4.bim

###########################################################################################################################################################################################

###############################################
# 4.0 SPLIT INTO EUR-NONJEWISH AND EUR-JEWISH #
###############################################


###########################################################################################################################################################################################


# exclude finish samples from PC estimation - they wont be included in the analysis:

cat ${path_gwas}pre_imputation/QC/farkkila_gsa/farkkila_gsa_postqc_preimp_eur.fam ${path_gwas}pre_imputation/QC/palotie_hus_gsa/palotie_hus_gsa_postqc_preimp_eur.fam ${path_gwas}pre_imputation/QC/finland_illugwas/finland_illugwas_postqc_preimp_eur.fam >> ${path_gwas}pre_imputation/QC/pca_1000gp/list_finish_samples

# No duplicates or related 
# split into Eur_nonJewish and EUR_nonJewish

cp ${path_gwas}pre_imputation/QC/pca_1000gp/list_eur_ancestry_samples_noDuplicates_perstudy_interstudy_2 ${path_gwas}pre_imputation/QC/pca_1000gp/list_eur_all_ancestry_samples_noDuplicates_perstudy_interstudy_2

ancestry=(eur_nonjewish eur_jewish eur_all)

for i in ${array[@]}
do 
for j in ${ancestry[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_noHighLD_noassoclogP4 \
--keep ${path_gwas}pre_imputation/QC/pca_1000gp/list_${j}_ancestry_samples_noDuplicates_perstudy_interstudy_2 \
--remove ${path_gwas}pre_imputation/QC/pca_1000gp/list_finish_samples \
--make-bed --out ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_${j}_pruned_noDuplicates
done
done

for i in ${array[@]}
do 
for j in ${ancestry[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_${j}_pruned_noDuplicates \
--remove ${path_gwas}pre_imputation/QC/pca_1000gp/list_eur_related_3rd_degree_samples \
--maf 0.001 \
--make-bed --out ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_${j}_pruned_noDuplicates_noRelated3rd
done
done

for i in ${array[@]}
do 
for j in ${ancestry[@]}
do
rm ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_${j}_pruned_noDuplicates.*
done
done

################
# eur_jewish

wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/*_all_studies_merged_eur_jewish_pruned_noDuplicates_noRelated3rd.fam
# 52 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix500_all_studies_merged_eur_jewish_pruned_noDuplicates_noRelated3rd.fam
# 1 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix6_all_studies_merged_eur_jewish_pruned_noDuplicates_noRelated3rd.fam
# 4428 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_jewish_pruned_noDuplicates_noRelated3rd.fam
# 145 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humancoreexome_all_studies_merged_eur_jewish_pruned_noDuplicates_noRelated3rd.fam
# 84 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomni1_all_studies_merged_eur_jewish_pruned_noDuplicates_noRelated3rd.fam
# 11 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomniexpress_all_studies_merged_eur_jewish_pruned_noDuplicates_noRelated3rd.fam
# 551 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina370_all_studies_merged_eur_jewish_pruned_noDuplicates_noRelated3rd.fam
# 166 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illuminaexome_all_studies_merged_eur_jewish_pruned_noDuplicates_noRelated3rd.fam
# 41 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/quad610_all_studies_merged_eur_jewish_pruned_noDuplicates_noRelated3rd.fam
# 5479 total

wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/*_all_studies_merged_eur_jewish_pruned_noDuplicates_noRelated3rd.bim
# 182458 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix500_all_studies_merged_eur_jewish_pruned_noDuplicates_noRelated3rd.bim
# 71809 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix6_all_studies_merged_eur_jewish_pruned_noDuplicates_noRelated3rd.bim
# 221384 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_jewish_pruned_noDuplicates_noRelated3rd.bim
# 154015 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humancoreexome_all_studies_merged_eur_jewish_pruned_noDuplicates_noRelated3rd.bim
# 333959 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomni1_all_studies_merged_eur_jewish_pruned_noDuplicates_noRelated3rd.bim
# 267772 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomniexpress_all_studies_merged_eur_jewish_pruned_noDuplicates_noRelated3rd.bim
# 187979 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina370_all_studies_merged_eur_jewish_pruned_noDuplicates_noRelated3rd.bim
# 158438 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illuminaexome_all_studies_merged_eur_jewish_pruned_noDuplicates_noRelated3rd.bim
# 282880 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/quad610_all_studies_merged_eur_jewish_pruned_noDuplicates_noRelated3rd.bim

################
# eur_nonjewish

wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/*_all_studies_merged_eur_nonjewish_pruned_noDuplicates_noRelated3rd.fam
# 4581 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix500_all_studies_merged_eur_nonjewish_pruned_noDuplicates_noRelated3rd.fam
# 8173 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix6_all_studies_merged_eur_nonjewish_pruned_noDuplicates_noRelated3rd.fam
# 38001 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_nonjewish_pruned_noDuplicates_noRelated3rd.fam
# 17276 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humancoreexome_all_studies_merged_eur_nonjewish_pruned_noDuplicates_noRelated3rd.fam
# 1937 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomni1_all_studies_merged_eur_nonjewish_pruned_noDuplicates_noRelated3rd.fam
# 1222 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomniexpress_all_studies_merged_eur_nonjewish_pruned_noDuplicates_noRelated3rd.fam
# 2888 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina370_all_studies_merged_eur_nonjewish_pruned_noDuplicates_noRelated3rd.fam
# 1175 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illuminaexome_all_studies_merged_eur_nonjewish_pruned_noDuplicates_noRelated3rd.fam
# 3343 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/quad610_all_studies_merged_eur_nonjewish_pruned_noDuplicates_noRelated3rd.fam
# 78596 total

wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/*_all_studies_merged_eur_nonjewish_pruned_noDuplicates_noRelated3rd.bim
# 184769 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix500_all_studies_merged_eur_nonjewish_pruned_noDuplicates_noRelated3rd.bim
# 252310 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix6_all_studies_merged_eur_nonjewish_pruned_noDuplicates_noRelated3rd.bim
# 221691 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_nonjewish_pruned_noDuplicates_noRelated3rd.bim
# 154730 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humancoreexome_all_studies_merged_eur_nonjewish_pruned_noDuplicates_noRelated3rd.bim
# 336097 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomni1_all_studies_merged_eur_nonjewish_pruned_noDuplicates_noRelated3rd.bim
# 298206 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomniexpress_all_studies_merged_eur_nonjewish_pruned_noDuplicates_noRelated3rd.bim
# 187994 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina370_all_studies_merged_eur_nonjewish_pruned_noDuplicates_noRelated3rd.bim
# 158536 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illuminaexome_all_studies_merged_eur_nonjewish_pruned_noDuplicates_noRelated3rd.bim
# 286020 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/quad610_all_studies_merged_eur_nonjewish_pruned_noDuplicates_noRelated3rd.bim

################
# eur_all

wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/*_all_studies_merged_eur_all_pruned_noDuplicates_noRelated3rd.fam
# 4633 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix500_all_studies_merged_eur_all_pruned_noDuplicates_noRelated3rd.fam
# 8174 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix6_all_studies_merged_eur_all_pruned_noDuplicates_noRelated3rd.fam
# 42429 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_all_pruned_noDuplicates_noRelated3rd.fam
# 17421 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humancoreexome_all_studies_merged_eur_all_pruned_noDuplicates_noRelated3rd.fam
# 2021 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomni1_all_studies_merged_eur_all_pruned_noDuplicates_noRelated3rd.fam
# 1233 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomniexpress_all_studies_merged_eur_all_pruned_noDuplicates_noRelated3rd.fam
# 3439 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina370_all_studies_merged_eur_all_pruned_noDuplicates_noRelated3rd.fam
# 1341 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illuminaexome_all_studies_merged_eur_all_pruned_noDuplicates_noRelated3rd.fam
# 3384 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/quad610_all_studies_merged_eur_all_pruned_noDuplicates_noRelated3rd.fam
# 84075 total

wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/*_all_studies_merged_eur_all_pruned_noDuplicates_noRelated3rd.bim
# 184769 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix500_all_studies_merged_eur_all_pruned_noDuplicates_noRelated3rd.bim
# 252310 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix6_all_studies_merged_eur_all_pruned_noDuplicates_noRelated3rd.bim
# 221691 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_all_pruned_noDuplicates_noRelated3rd.bim
# 154730 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humancoreexome_all_studies_merged_eur_all_pruned_noDuplicates_noRelated3rd.bim
# 336097 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomni1_all_studies_merged_eur_all_pruned_noDuplicates_noRelated3rd.bim
# 298206 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomniexpress_all_studies_merged_eur_all_pruned_noDuplicates_noRelated3rd.bim
# 187994 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina370_all_studies_merged_eur_all_pruned_noDuplicates_noRelated3rd.bim
# 158536 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illuminaexome_all_studies_merged_eur_all_pruned_noDuplicates_noRelated3rd.bim
# 286020 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/quad610_all_studies_merged_eur_all_pruned_noDuplicates_noRelated3rd.bim

################################################################################################################################

cp ${path_gwas}pre_imputation/QC/pca_1000gp/list_eur_ancestry_samples_withDuplicates ${path_gwas}pre_imputation/QC/pca_1000gp/list_eur_all_ancestry_samples_withDuplicates

ancestry=(eur_nonjewish eur_jewish)
  
# all samples
for i in ${array[@]}
do 
for j in ${ancestry[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_noHighLD_noassoclogP4 \
--keep ${path_gwas}pre_imputation/QC/pca_1000gp/list_${j}_ancestry_samples_withDuplicates  \
--remove ${path_gwas}pre_imputation/QC/pca_1000gp/list_finish_samples \
--extract ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_${j}_pruned_noDuplicates_noRelated3rd.bim \
--make-bed --out ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_${j}_pruned_withDuplicates
done
done


################
# eur_jewish

wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/*_all_studies_merged_eur_jewish_pruned_withDuplicates.fam
# 52 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix500_all_studies_merged_eur_jewish_pruned_withDuplicates.fam
# 1 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix6_all_studies_merged_eur_jewish_pruned_withDuplicates.fam
# 7962 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_jewish_pruned_withDuplicates.fam
# 162 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humancoreexome_all_studies_merged_eur_jewish_pruned_withDuplicates.fam
# 95 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomni1_all_studies_merged_eur_jewish_pruned_withDuplicates.fam
# 11 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomniexpress_all_studies_merged_eur_jewish_pruned_withDuplicates.fam
# 772 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina370_all_studies_merged_eur_jewish_pruned_withDuplicates.fam
# 322 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illuminaexome_all_studies_merged_eur_jewish_pruned_withDuplicates.fam
# 41 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/quad610_all_studies_merged_eur_jewish_pruned_withDuplicates.fam
# 9418 total


wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/*_all_studies_merged_eur_jewish_pruned_withDuplicates.bim
# 182458 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix500_all_studies_merged_eur_jewish_pruned_withDuplicates.bim
# 71809 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix6_all_studies_merged_eur_jewish_pruned_withDuplicates.bim
# 221384 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_jewish_pruned_withDuplicates.bim
# 154015 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humancoreexome_all_studies_merged_eur_jewish_pruned_withDuplicates.bim
# 333959 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomni1_all_studies_merged_eur_jewish_pruned_withDuplicates.bim
# 267772 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomniexpress_all_studies_merged_eur_jewish_pruned_withDuplicates.bim
# 187979 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina370_all_studies_merged_eur_jewish_pruned_withDuplicates.bim
# 158438 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illuminaexome_all_studies_merged_eur_jewish_pruned_withDuplicates.bim
# 282880 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/quad610_all_studies_merged_eur_jewish_pruned_withDuplicates.bim

################
# eur_nonjewish

wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/*_all_studies_merged_eur_nonjewish_pruned_withDuplicates.fam
# 4600 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix500_all_studies_merged_eur_nonjewish_pruned_withDuplicates.fam
# 11086 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix6_all_studies_merged_eur_nonjewish_pruned_withDuplicates.fam
# 72451 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_nonjewish_pruned_withDuplicates.fam
# 22426 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humancoreexome_all_studies_merged_eur_nonjewish_pruned_withDuplicates.fam
# 2606 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomni1_all_studies_merged_eur_nonjewish_pruned_withDuplicates.fam
# 1234 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomniexpress_all_studies_merged_eur_nonjewish_pruned_withDuplicates.fam
# 4861 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina370_all_studies_merged_eur_nonjewish_pruned_withDuplicates.fam
# 2235 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illuminaexome_all_studies_merged_eur_nonjewish_pruned_withDuplicates.fam
# 3355 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/quad610_all_studies_merged_eur_nonjewish_pruned_withDuplicates.fam
# 124854 total


wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/*_all_studies_merged_eur_nonjewish_pruned_withDuplicates.bim
# 184769 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix500_all_studies_merged_eur_nonjewish_pruned_withDuplicates.bim
# 252310 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix6_all_studies_merged_eur_nonjewish_pruned_withDuplicates.bim
# 221691 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_nonjewish_pruned_withDuplicates.bim
# 154730 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humancoreexome_all_studies_merged_eur_nonjewish_pruned_withDuplicates.bim
# 336097 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomni1_all_studies_merged_eur_nonjewish_pruned_withDuplicates.bim
# 298206 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomniexpress_all_studies_merged_eur_nonjewish_pruned_withDuplicates.bim
# 187994 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina370_all_studies_merged_eur_nonjewish_pruned_withDuplicates.bim
# 158536 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illuminaexome_all_studies_merged_eur_nonjewish_pruned_withDuplicates.bim
# 286020 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/quad610_all_studies_merged_eur_nonjewish_pruned_withDuplicates.bim

################
# eur_all

wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/*_all_studies_merged_eur_all_pruned_withDuplicates.fam
# 4652 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix500_all_studies_merged_eur_all_pruned_withDuplicates.fam
# 11087 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix6_all_studies_merged_eur_all_pruned_withDuplicates.fam
# 80413 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_all_pruned_withDuplicates.fam
# 22588 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humancoreexome_all_studies_merged_eur_all_pruned_withDuplicates.fam
# 2701 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomni1_all_studies_merged_eur_all_pruned_withDuplicates.fam
# 1245 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomniexpress_all_studies_merged_eur_all_pruned_withDuplicates.fam
# 5633 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina370_all_studies_merged_eur_all_pruned_withDuplicates.fam
# 2557 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illuminaexome_all_studies_merged_eur_all_pruned_withDuplicates.fam
# 3396 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/quad610_all_studies_merged_eur_all_pruned_withDuplicates.fam
# 134272 total

wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/*_all_studies_merged_eur_all_pruned_withDuplicates.bim
# 184769 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix500_all_studies_merged_eur_all_pruned_withDuplicates.bim
# 252310 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix6_all_studies_merged_eur_all_pruned_withDuplicates.bim
# 221691 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_all_pruned_withDuplicates.bim
# 154730 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humancoreexome_all_studies_merged_eur_all_pruned_withDuplicates.bim
# 336097 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomni1_all_studies_merged_eur_all_pruned_withDuplicates.bim
# 298206 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomniexpress_all_studies_merged_eur_all_pruned_withDuplicates.bim
# 187994 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina370_all_studies_merged_eur_all_pruned_withDuplicates.bim
# 158536 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illuminaexome_all_studies_merged_eur_all_pruned_withDuplicates.bim
# 286020 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/quad610_all_studies_merged_eur_all_pruned_withDuplicates.bim


###########################################################################################################################################################################################


#######################################################
# 5.0 ESTIMATE PC, EXPORT PCS, AND ALLELE FREQUENCIES #
#######################################################


# The following command exports PCs to project onto, along with the allele frequencies needed to calibrate the 'variance-standardize' operation:
 

array=(illumina370 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 illuminaexome gsa)
path_gwas=/path/to/ibdgwas/IIBDGC/

ancestry=(eur_all eur_jewish)
MEM=30000 

ancestry=(eur_nonjewish)
MEM=20000 

for i in ${array[@]}
do
for j in ${ancestry[@]}
do
bsub -J"bcftools" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal -n 8 \
-e ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stderr_${i}_pca_${j}_nodup_pcs_for_analysis \
-o ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stdout_${i}_pca_${j}_nodup_pcs_for_analysis \
"/path/to/software/./plink2  \
--bfile ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_${j}_pruned_noDuplicates_noRelated3rd \
--freq counts \
--pca '40' approx allele-wts \
--out ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_${j}_pruned_ref_pcs \
--threads 8 --memory $MEM"
done
done


array=(illumina370 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 gsa illuminaexome)
path_gwas=/path/to/ibdgwas/IIBDGC/
  
ancestry=(eur_nonjewish eur_all eur_jewish)


for i in ${array[@]}
do
echo ${i} && for j in ${ancestry[@]}
do
echo ${j} && tail -50 ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stdout_${i}_pca_${j}_nodup_pcs_for_analysis | grep -E "Successfully|Exited with"
done
done

###############
# illumina370
# eur_nonjewish
# Successfully completed.
# eur_jewish
# Successfully completed.
# affymetrix6
# eur_nonjewish
# Successfully completed.
# eur_jewish
# Exited with exit code 13. # ONLY ONE SAMPLE
# humanomniexpress
# eur_nonjewish
# Successfully completed.
# eur_jewish
# Exited with exit code 13. # ONLY ELEVEN SAMPLES
# affymetrix500
# eur_nonjewish
# Successfully completed.
# eur_jewish
# Successfully completed.
# humancoreexome
# eur_nonjewish
# Successfully completed.
# eur_jewish
# Successfully completed.
# humanomni1
# eur_nonjewish
# Successfully completed.
# eur_jewish
# Successfully completed.
# quad610
# eur_nonjewish
# Successfully completed.
# eur_jewish
# Exited with exit code 13. # ONLY FORTY-ONE SAMPLES
# gsa
# eur_nonjewish
# Successfully completed.
# eur_jewish
# Successfully completed.
# illuminaexome
# eur_nonjewish
# Successfully completed.
# eur_jewish
# Successfully completed.
###############

###########################################
## 1.2 - PROJECT ONTO THOSE PCS ALL IIBDGC
       
array=(illumina370 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 gsa illuminaexome)
ancestry=(eur_nonjewish eur_jewish eur_all)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=5000 # next time only this


for i in ${array[@]}
do
for j in ${ancestry[@]}
do
bsub -J"bcftools" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal -n 10 \
-e ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stdout_${i}_pca_${j}_withdup_pcs_for_analysis \
-o ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stdout_${i}_pca_${j}_withdup_pcs_for_analysis \
"/path/to/software/./plink2 \
--bfile ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_${j}_pruned_withDuplicates \
--autosome \
--read-freq ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_ref_pcs.acount \
--score ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_${j}_pruned_ref_pcs.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
--score-col-nums 6-45 \
--out ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_${j}_pruned_withDuplicates_ref_pcs_new_projection"
done
done




for i in ${array[@]}
do
echo ${i} && for j in ${ancestry[@]}
do
echo ${j} && tail -50 ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stdout_${i}_pca_${j}_withdup_pcs_for_analysis | grep -E "Successfully|Exited with"
done
done

##########
# illumina370
# eur_nonjewish
# Successfully completed.
# eur_jewish
# Successfully completed.
# affymetrix6
# eur_nonjewish
# Successfully completed.
# eur_jewish
# Exited with exit code 3. #### 
# humanomniexpress
# eur_nonjewish
# Successfully completed.
# eur_jewish
# Exited with exit code 3. #### 
# affymetrix500
# eur_nonjewish
# Successfully completed.
# eur_jewish
# Successfully completed.
# humancoreexome
# eur_nonjewish
# Successfully completed.
# eur_jewish
# Successfully completed.
# humanomni1
# eur_nonjewish
# Successfully completed.
# eur_jewish
# Successfully completed.
# quad610
# eur_nonjewish
# Successfully completed.
# eur_jewish
# Exited with exit code 3.. #### 
# gsa
# eur_nonjewish
# Successfully completed.
# eur_jewish
# Successfully completed.
# illuminaexome
# eur_nonjewish
# Successfully completed.
# eur_jewish
# Successfully completed.
##########



############################
### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(viridis)
library(ggsci)
library(lmtest)

path<-"/path/to/ibdgwas/IIBDGC/"

array<-c("illumina370","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")

illumina370<-c("belgium_inf1_old_gwas","belgium_inf2_old_gwas","swedish_uc_old_gwas","niddk_old_gwas")
affymetrix6<-c("german_affy6_old_gwas","norway_affy6_old_gwas","gwas2")
humanomniexpress<-c("australia_omniexome")
affymetrix500<-c("gwas1")
humancoreexome<-c("all_hce")
humanomni1<-c("pittsburgh_gsa")
quad610<-c("spain_gsa")
gsa<-c("italy_gsa","kiel_austria_sibdcs_gsa"
       ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
       ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
       ,"mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa","franke_gsa","hyams_protect_gsa",
       "lewis_sparc_gsa","mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
       "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
       "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
       "xavier_share_gsa")
illuminaexome<-c("prism_nfe_gwas","helmsley_prism_gsa","helmsley_xavier_prism_gsa","finland_illugwas")

cohorts<-c(illumina370,affymetrix6,humanomniexpress,affymetrix500,humancoreexome,humanomni1,quad610,gsa,illuminaexome)


length(cohorts)
# [1] 55
# ok no cedars old studies


for (i in 1:length(cohorts)){
  a<-read.table(paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_withDuplicates.fam",sep=""),head=F)
  a$cohort<-cohorts[i]
  if(i==1){
    fam<-a
  }else{
    fam<-rbind(a,fam)
  }
}

dim(fam)
# [1] 144238      7

# add country of origin

fam$country<-NA
fam$country[which(fam$cohort %in% c("australia_omniexome"))]<-"Australia"
fam$country[which(fam$cohort %in% c("gwas1","gwas2","all_hce"))]<-"UK"
fam$country[which(fam$cohort %in% c("spain_gsa","basque_gsa"))]<-"Spain"
fam$country[which(fam$cohort %in% c("kiel_austria_sibdcs_gsa","german_affy6_old_gwas"))]<-"Germany"                 
fam$country[which(fam$cohort %in% c("italy_gsa"))]<-"Italy"   
fam$country[which(fam$cohort %in% c("netherlands_gsa"))]<-"Netherlands"                     
fam$country[which(fam$cohort %in% c("slovenia_gsa"))]<-"Slovenia"                     
fam$country[which(fam$cohort %in% c("sweden_gsa","swedish_uc_old_gwas"))]<-"Sweden"                    
fam$country[which(fam$cohort %in% c("finland_illugwas","palotie_hus_gsa","farkkila_gsa"))]<-"Finland"
fam$country[which(fam$cohort %in% c("norway_affy6_old_gwas"))]<-"Norway"
fam$country[which(fam$cohort %in% c("belgium_inf1_old_gwas","belgium_inf2_old_gwas","belgium_louis_gsa","belgium_franchimont_gsa"
                                    ,"belgium_vermeire_gsa"))]<-"Belgium"
fam$country[which(fam$cohort %in% c("lithuania_gsa"))]<-"Lithuania"
fam$country[which(fam$cohort %in% c("prism_nfe_gwas","prism_nfe_gsa","cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","cedars_gsa",
                                    "niddk_broad_gsa","niddk_feinstein_gsa","niddk_old_gwas","pittsburgh_gsa","mccauley_gsa","ccfa_gsa",
                                    "mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","franchimont_gsa",
                                    "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
                                    "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
                                    "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","pekow_share_gsa",
                                    "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
                                    "xavier_share_gsa"))]<-"USA"


table(fam$country,useNA="ifany")
fam$country<-as.factor(fam$country)
fam$country<-factor(fam$country, levels=c("Spain","Italy","Slovenia","Belgium","Germany","Netherlands","UK","Lithuania",
                                          "Sweden","Norway","Finland","Australia","USA"))

table(fam$cohort,fam$country,useNA="ifany")

# Add jewish ancestry - ## latest udpate phenotype file - sent by Phil on 18/Jan/2022
pheno<-read.csv("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes-827da7e_updated.csv.gz",head=T)

colnames(pheno)[16]<-"self_jewish"
pheno_cd<-fread(paste(path,"pheno/list_cedars_old_gwas_ids_TH.txt",sep=""),head=T,sep="\t")
pheno<-rbind(pheno[,c("sample_id","self_jewish")],pheno_cd[,c("sample_id","self_jewish")])
pheno$self_jewish<-as.character(pheno$self_jewish)
pheno$self_jewish[which(pheno$self_jewish %in% c("","Unknown"))]<-"Unknown"
pheno$self_jewish[which(pheno$self_jewish %in% c("No","Non-Jewish"))]<-"Non-Jewish"
pheno$self_jewish[which(pheno$self_jewish %in% c("Yes","Jewish"))]<-"Jewish"
table(pheno$self_jewish)
# Jewish Non-Jewish    Unknown 
#   7193      41195      74089 


fam$self_jewish<-"Unknown"
fam$self_jewish[which(fam$V1 %in% pheno$sample_id[which(pheno$self_jewish=="Non-Jewish")])]<-"Non-Jewish"
fam$self_jewish[which(fam$V1 %in% pheno$sample_id[which(pheno$self_jewish=="Jewish")])]<-"Jewish"
table(fam$self_jewish)
# Jewish Non-Jewish    Unknown 
# 5779      36506     101953


tmp<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_EUR_withDuplicates_ref_pcs_new_projection_edited.sscore",sep=""),head=T,sep="\t")
fam$pca_jewish<-NA
fam$pca_jewish[which(fam$V1 %in% tmp$FID[which(tmp$pca_jewish=="Non-Jewish")])]<-"Non-Jewish"
fam$pca_jewish[which(fam$V1 %in% tmp$FID[which(tmp$pca_jewish=="Jewish")])]<-"Jewish"
rm(tmp)


# add pheno CD/UC/IBDu/Ctr:
pheno<-read.table(paste(path,"pheno/phenotype_data_for_all_cohorts_Dec2022_analysis.tsv",sep=""),head=T)

fam$pheno<-NA
fam$pheno[which(fam$V1 %in% pheno$FID[which(pheno$ibd==0)])]<-"Control"
fam$pheno[which(fam$V1 %in% pheno$FID[which(pheno$cd==0)])]<-"Control"
fam$pheno[which(fam$V1 %in% pheno$FID[which(pheno$uc==0)])]<-"Control"

fam$pheno[which(fam$V1 %in% pheno$FID[which(pheno$ibd==1)])]<-"IBDu"
fam$pheno[which(fam$V1 %in% pheno$FID[which(pheno$cd==1)])]<-"CD"
fam$pheno[which(fam$V1 %in% pheno$FID[which(pheno$uc==1)])]<-"UC"

fam$pheno<-factor(fam$pheno, levels=c("CD","IBDu","UC","Control"))

# array by array explore the PCS and estimate how many to include in the analysis



######################################################################################################################################################



ancestry<-c("eur_nonjewish","eur_all","eur_jewish")



for (jj in c(1:length("ancestry"))) {

  
  ######################################################################################################################################################
  
  ###############################################
  # 1.- TOTAL VARIANCE EXPLAINED BY EACH PC:
  
  
  if(ancestry[jj]=="eur_jewish") {
    array<-c("illumina370","gsa")
  }else{
    array<-c("illumina370","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")
  }
    
  for (i in 1:length(array)) {
    
    eigenval<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/",array[i],"_all_studies_merged_",ancestry[jj],"_pruned_ref_pcs.eigenval",sep=""),head=F)
    
    eigenval$var_exp<-NA
    for (j in 1:nrow(eigenval)){
      eigenval$var_exp[j]<-(eigenval$V1[j] / sum(eigenval$V1))*100
    }
    
    eigenval$PC<-paste("PC",seq(1:nrow(eigenval)),sep="")
    eigenval$cumsum_var_expl<-cumsum(eigenval$var_exp)
    eigenval$array<-array[i]
    
    
    if(i==1) {
      all_eigen<-eigenval
    } else {
      all_eigen<-rbind(all_eigen,eigenval)
    }
    
  }
  
  all_eigen[which(all_eigen$PC=="PC1"),]
  all_eigen$PC<-factor(all_eigen$PC,levels=paste("PC",seq(1:nrow(eigenval)),sep=""))
  
  p1<-ggplot(all_eigen, aes(x = PC, y=cumsum_var_expl, group = array, color = array)) + ylab("Cummulative Variance explained") +
    geom_point() + geom_line()+ theme(axis.text.x=element_text(angle=45, hjust=1))
  
  pdf(paste(path,"pre_imputation/QC/pca_1000gp/cummulative_variance_explained_40PCs_for_analysis_per_array_",ancestry[jj],".pdf",sep=""),width = 10, height = 6)
  print(p1)
  dev.off()
  
  
  p1<-ggplot(all_eigen, aes(x = PC, y=var_exp, group = array, color = array)) + ylab("Variance explained") +
    geom_point() + geom_line()+ theme(axis.text.x=element_text(angle=45, hjust=1)) 
  # +scale_y_reverse()
  
  pdf(paste(path,"pre_imputation/QC/pca_1000gp/variance_explained_40PCs_for_analysis_per_array_",ancestry[jj],".pdf",sep=""),width = 10, height = 6)
  print(p1)
  dev.off()
  
  paste(path,"pre_imputation/QC/pca_1000gp/variance_explained_40PCs_for_analysis_per_array_",ancestry[jj],".pdf",sep="")
  
  ###############################################
  # 2.- EXPLORE VARIANT WEIGHTS PER PC:
  
  
  # plot the weights per variant - variant at a particular cluster - related with LD regions
  
  
  for (ii in 1:length(array)) {
    
    
    score<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/",array[ii],"_all_studies_merged_",ancestry[jj],"_pruned_ref_pcs.eigenvec.allele",sep=""),head=F)
    colnames(score)<-c("chr","snp","alt","ref","allele",paste("PC",seq(1:40),sep=""))
    
    score$chr<-as.factor(score$chr)
    score$chr<-factor(score$chr,levels=paste(seq(1:22)))
    
    score$pos<-gsub("[0-9]{1,2}:","",score$snp)
    score$pos<-as.numeric(gsub("_.*","",score$pos))
    
    score<-score[order(score$chr,score$pos),]
    
    score$index<-seq(1:nrow(score))
    
    for (i in 1:8){
      print(i)
      tmp<-score[,c("index","chr",paste("PC",i,sep=""))]
      colnames(tmp)<-c("index","chr","PC")
      
      
      pdf(paste(path,"pre_imputation/QC/pca_1000gp/PC",i,"_loadings_",array[ii],"_",ancestry[jj],".pdf",sep=""),width=7,height=4)
      print(ggplot(tmp, aes(index, PC)) + geom_point(aes(colour = chr)) + scale_color_viridis(discrete=T) + 
              ggtitle(paste("PC",i,sep="")))
      dev.off()
      
      rm(tmp)
    }
    
    rm(score)
    
  }
  
  
  
  for (ii in 1:length(array)) {
    
    score<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/",array[ii],"_all_studies_merged_",ancestry[jj],"_pruned_ref_pcs.eigenvec.allele",sep=""),head=F)
    colnames(score)<-c("chr","snp","alt","ref","allele",paste("PC",seq(1:40),sep=""))
    
    score$chr<-as.factor(score$chr)
    score$chr<-factor(score$chr,levels=paste(seq(1:22)))
    
    score$pos<-gsub("[0-9]{1,2}:","",score$snp)
    score$pos<-as.numeric(gsub("_.*","",score$pos))
    
    score<-score[order(score$chr,score$pos),]
    
    score$index<-seq(1:nrow(score))
    
    for (i in 1:10){
      
      print(i)
      
      tmp<-score[,c("index","chr","pos",paste("PC",i,sep=""))]
      colnames(tmp)<-c("index","chr","pos","PC")
      lim<-mean(tmp$PC)+12*sd(tmp$PC)
      tmp1<-tmp[which(tmp$PC>=lim | tmp$PC<= (-1*lim)),]
      
      if(nrow(tmp1)>0) {
        
        x<-as.data.frame(tapply(tmp1$pos, as.character(tmp1$chr), range))
        
        x$chr<-rownames(x)
        x$array<-array[ii]
        x$pc<-i
        x$n_snp<-nrow(tmp1)
        
        if(!exists("ld_regions")) {
          ld_regions<-x
        }else{
          ld_regions<-rbind(ld_regions,x)
        }
        rm(x)
      }
      rm(tmp,tmp1)
    }
    rm(score)
  }
  
  
  ld_regions_tmp<-ld_regions
  
  summary(ld_regions$n_snp)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 118.0   118.0   118.0   124.6   132.0   132.0 
  
  
  ld_regions<-ld_regions[which(ld_regions$n_snp>50),c(1:3,5)]
  ld_regions<-ld_regions[!duplicated(ld_regions),]
  colnames(ld_regions)[1]<-"range"
  
  ld_regions$start<-gsub(",.*","",ld_regions$range)
  ld_regions$start<-gsub("\\c\\(","",ld_regions$start)
  ld_regions$end<-gsub(".*,","",ld_regions$range)
  ld_regions$end<-gsub("\\)","",ld_regions$end)
  ld_regions$bp<-as.numeric(ld_regions$end)-as.numeric(ld_regions$start)
  ld_regions[which(ld_regions$bp>0 & ld_regions$bp<2500000),]
  # 0
  
  
  
  ###############################################
  # 2.- EXPLORE EACH PC FOR:
  
  
  data<-matrix(nrow=length(array),ncol=11)
  data<-as.data.frame(data)
  colnames(data)<-c("array","m0_vs_m1","m1_vs_m2","m2_vs_m3","m3_vs_m4","m4_vs_m5","m5_vs_m6","m6_vs_m7","m7_vs_m8","m8_vs_m9","m9_vs_m10")
  
  for (ii in 1:length(array)) {
    
    print(array[ii])
    pca<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/",array[ii],"_all_studies_merged_",ancestry[jj],"_pruned_withDuplicates_ref_pcs_new_projection.sscore",sep=""),head=F)
    colnames(pca)<-c("FID","IID","PHENO1","ALLELE_CT","NAMED_ALLELE_DOSAGE_SUM",paste("PC",seq(1:40),sep=""))
    pca<-merge(pca,fam[,c("V1","cohort","country","self_jewish","pca_jewish","pheno")],by.x="FID",by.y="V1",all.x=T)
    
    # pca$country<-as.character(pca$country)
    # pca$pheno<-as.character(pca$pheno)
    # pca$self_jewish<-as.character(pca$self_jewish)
    
    # # In Australia cohort only 1 self reported jewish, density needs at least 2 points
    
    if(ancestry[ii]!="eur_jewish") {
      if(ii==3) {
        pca$self_jewish<-"Unknown"
      }
      
      # in affymetrix500 only 1 pca-jewish , density needs at least 2 points
      if(ii==2) {
        pca$pca_jewish<-"Non-Jewish"
      }
      
      
      if (ii==4){
        pca$pheno[which(pca$pheno=="IBDu")]<-NA
      }
    }

    if(ancestry[ii]!="eur_jewish") {
      # only one individual - no density
      pca$cohort[which(pca$cohort=="weersma_gsa")]<-NA
    }
    
    print(table(as.character(pca$cohort)))
    print(table(as.character(pca$country)))
    print(table(pca$self_jewish))
    print(table(pca$pheno))
    
    
    for (j in 1:10) {
      
      # by Cohort:
      sp <- split(pca[,5+j], pca$cohort)
      a <- lapply(seq_along(sp), function(i){
        d <- density(sp[[i]])
        k <- which.max(d$y)
        data.frame(cohort = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
      })
      a <- do.call(rbind, a)
      
      p1<-ggplot(pca, aes(x=pca[,5+j],color=cohort, fill=cohort)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = cohort, vjust = - 0.5)) + 
        ggtitle("Cohort") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.2) + theme(legend.title = element_blank()) + 
        xlab(paste("PC",j,sep=""))
      
      # by Country:
      sp <- split(pca[,5+j], as.character(pca$country))
      a <- lapply(seq_along(sp), function(i){
        d <- density(sp[[i]])
        k <- which.max(d$y)
        data.frame(country = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
      })
      a <- do.call(rbind, a)
      
      p2<-ggplot(pca, aes(x=pca[,5+j],color=country, fill=country)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = country, vjust = - 0.5)) + 
        ggtitle("Country") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.2) + theme(legend.title = element_blank()) + 
        xlab(paste("PC",j,sep=""))
      
      if (ancestry[jj]=="eur_all") {
        
        # by Jewish:
        sp <- split(pca[,5+j], as.character(pca$self_jewish))
        a <- lapply(seq_along(sp), function(i){
          d <- density(sp[[i]])
          k <- which.max(d$y)
          data.frame(self_jewish = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
        })
        a <- do.call(rbind, a)
        
        p3<-ggplot(pca, aes(x=pca[,5+j],color=self_jewish, fill=self_jewish)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = self_jewish, vjust = - 0.5)) +
          ggtitle("self-reported Jewish") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.2) + theme(legend.title = element_blank()) +
          xlab(paste("PC",j,sep="")) + scale_color_jco() + scale_fill_jco()
        
        
        # by Jewish:
        sp <- split(pca[,5+j], pca$pca_jewish)
        a <- lapply(seq_along(sp), function(i){
          d <- density(sp[[i]])
          k <- which.max(d$y)
          data.frame(pca_jewish = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
        })
        a <- do.call(rbind, a)
        
        p4<-ggplot(pca, aes(x=pca[,5+j],color=pca_jewish, fill=pca_jewish)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = pca_jewish, vjust = - 0.5)) +
          ggtitle("Jewish (inferred by previous PCA, all cohorts together)") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.2) + theme(legend.title = element_blank()) +
          xlab(paste("PC",j,sep="")) + scale_color_jco() + scale_fill_jco()
        
        # by Phenotype:
        sp <- split(pca[,5+j], as.character(pca$pheno))
        a <- lapply(seq_along(sp), function(i){
          d <- density(sp[[i]])
          k <- which.max(d$y)
          data.frame(pheno = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
        })
        a <- do.call(rbind, a)
        
        p5<-ggplot(pca, aes(x=pca[,5+j],color=pheno, fill=pheno)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = pheno, vjust = - 0.5)) + 
          ggtitle("Phenotype") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.2) + theme(legend.title = element_blank()) + 
          xlab(paste("PC",j,sep="")) + scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2")
        
        fig<-ggarrange(p1,p2,p3,p4,p5,nrow=5)
        
        
        pdf(paste(path,"pre_imputation/QC/pca_1000gp/histogram_european_ancestry_samples_PC",j,"_per_array_",array[ii],"_",ancestry[jj],".pdf",sep=""),width=14,height=20)
        print(annotate_figure(fig,top = text_grob(array[ii], face = "bold",size=14)))
        dev.off()
        
      } else {
        # by Phenotype:
        sp <- split(pca[,5+j], as.character(pca$pheno))
        a <- lapply(seq_along(sp), function(i){
          d <- density(sp[[i]])
          k <- which.max(d$y)
          data.frame(pheno = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
        })
        a <- do.call(rbind, a)
        
        p5<-ggplot(pca, aes(x=pca[,5+j],color=pheno, fill=pheno)) + geom_density(alpha=0.4) + geom_text(data = a, aes(x = xmax, y = ymax,label = pheno, vjust = - 0.5)) + 
          ggtitle("Phenotype") +  theme(legend.position="bottom") + ylim(0,max(a$ymax)+(max(a$ymax))*0.2) + theme(legend.title = element_blank()) + 
          xlab(paste("PC",j,sep="")) + scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2")
        
        fig<-ggarrange(p1,p2,p5,nrow=3)
        
        pdf(paste(path,"pre_imputation/QC/pca_1000gp/histogram_european_ancestry_samples_PC",j,"_per_array_",array[ii],"_",ancestry[jj],".pdf",sep=""),width=14,height=20)
        print(annotate_figure(fig,top = text_grob(array[ii], face = "bold",size=14)))
        dev.off()
        
      }
      
    }
    
    # create association report with phenotype, only nondup or non related samples:
    
    # see pre-regenie
    fam_analysis<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/",array[ii],"_all_studies_merged_",ancestry[jj],"_pruned_noDuplicates_noRelated3rd.fam",sep=""),head=F)
    
    pca$ibd<-NA
    pca$ibd[which( (pca$pheno %in% c("CD","IBDu","UC")) & (pca$FID %in% fam_analysis$V1))]<-1
    pca$ibd[which( (pca$pheno %in% c("Control")) & (pca$FID %in% fam_analysis$V1))]<-0
    
    
    m0_ibd<-glm(ibd ~ 1, data = pca, family = "binomial")
    m1_ibd<-glm(ibd ~ PC1, data = pca, family = "binomial")
    m2_ibd<-glm(ibd ~ PC1 + PC2, data = pca, family = "binomial")
    m3_ibd<-glm(ibd ~ PC1 + PC2 + PC3, data = pca, family = "binomial")
    m4_ibd<-glm(ibd ~ PC1 + PC2 + PC3 + PC4, data = pca, family = "binomial")
    m5_ibd<-glm(ibd ~ PC1 + PC2 + PC3 + PC4 + PC5, data = pca, family = "binomial")
    m6_ibd<-glm(ibd ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6, data = pca, family = "binomial")
    m7_ibd<-glm(ibd ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7, data = pca, family = "binomial")
    m8_ibd<-glm(ibd ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8, data = pca, family = "binomial")
    m9_ibd<-glm(ibd ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9, data = pca, family = "binomial")
    m10_ibd<-glm(ibd ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = pca, family = "binomial")
    
    data$array[ii]<-array[ii]
    data$m0_vs_m1[ii]<-lrtest(m0_ibd,m1_ibd)[["Pr(>Chisq)"]][2]
    data$m1_vs_m2[ii]<-lrtest(m1_ibd,m2_ibd)[["Pr(>Chisq)"]][2]
    data$m2_vs_m3[ii]<-lrtest(m2_ibd,m3_ibd)[["Pr(>Chisq)"]][2]
    data$m3_vs_m4[ii]<-lrtest(m3_ibd,m4_ibd)[["Pr(>Chisq)"]][2]
    data$m4_vs_m5[ii]<-lrtest(m4_ibd,m5_ibd)[["Pr(>Chisq)"]][2]
    data$m5_vs_m6[ii]<-lrtest(m5_ibd,m6_ibd)[["Pr(>Chisq)"]][2]
    data$m6_vs_m7[ii]<-lrtest(m6_ibd,m7_ibd)[["Pr(>Chisq)"]][2]
    data$m7_vs_m8[ii]<-lrtest(m7_ibd,m8_ibd)[["Pr(>Chisq)"]][2]
    data$m8_vs_m9[ii]<-lrtest(m8_ibd,m9_ibd)[["Pr(>Chisq)"]][2]
    data$m9_vs_m10[ii]<-lrtest(m9_ibd,m10_ibd)[["Pr(>Chisq)"]][2]
    
  }
  
  write.table(data,paste(path,"pre_imputation/QC/pca_1000gp/PC_selection_all_studies_",ancestry[jj],"_perarray_2022",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  
  
  
  # for (ii in 1:length(array)) {
  #   print(array[ii])
  #   print(all_eigen[which(all_eigen$array==array[ii] & all_eigen$var_exp>=3),])
  # }
  
  
  if (ancestry[jj]=="eur_all") {
    tmp<-all_eigen[which( (all_eigen$array=="illumina370" & all_eigen$PC %in% paste("PC",seq(1:2),sep="")) |
                            (all_eigen$array=="affymetrix6" & all_eigen$PC %in% paste("PC",seq(1),sep="")) |
                            (all_eigen$array=="humanomniexpress" & all_eigen$PC %in% paste("PC",seq(1),sep="")) |
                            (all_eigen$array=="affymetrix500" & all_eigen$PC %in% paste("PC",seq(1:3),sep="")) |
                            (all_eigen$array=="humancoreexome" & all_eigen$PC %in% paste("PC",seq(1:2),sep="")) |
                            (all_eigen$array=="humanomni1" & all_eigen$PC %in% paste("PC",seq(1:3),sep="")) |
                            (all_eigen$array=="quad610" & all_eigen$PC %in% paste("PC",seq(1),sep="")) |
                            (all_eigen$array=="gsa" & all_eigen$PC %in% paste("PC",seq(1:6),sep="")) |
                            (all_eigen$array=="illuminaexome" & all_eigen$PC %in% paste("PC",seq(1),sep=""))
    ),]
  }
  
  
  if (ancestry[jj]=="eur_nonjewish") {
    tmp<-all_eigen[which( (all_eigen$array=="illumina370" & all_eigen$PC %in% paste("PC",seq(1:5),sep="")) |
                            (all_eigen$array=="affymetrix6" & all_eigen$PC %in% paste("PC",seq(1),sep="")) |
                            (all_eigen$array=="humanomniexpress" & all_eigen$PC %in% paste("PC",seq(1),sep="")) |
                            (all_eigen$array=="affymetrix500" & all_eigen$PC %in% paste("PC",seq(1:3),sep="")) |
                            (all_eigen$array=="humancoreexome" & all_eigen$PC %in% paste("PC",seq(1:2),sep="")) |
                            (all_eigen$array=="humanomni1" & all_eigen$PC %in% paste("PC",seq(1:2),sep="")) |
                            (all_eigen$array=="quad610" & all_eigen$PC %in% paste("PC",seq(1),sep="")) |
                            (all_eigen$array=="gsa" & all_eigen$PC %in% paste("PC",seq(1:7),sep="")) |
                            (all_eigen$array=="illuminaexome" & all_eigen$PC %in% paste("PC",seq(1:4),sep=""))
    ),]
  }
  
  if (ancestry[jj]=="eur_jewish") {
    tmp<-all_eigen[which( (all_eigen$array=="illumina370" & all_eigen$PC %in% paste("PC",seq(1:2),sep="")) |
                            (all_eigen$array=="gsa" & all_eigen$PC %in% paste("PC",seq(1:2),sep=""))
    ),]
  }
  
  
  p1<-ggplot(tmp, aes(x = PC, y=cumsum_var_expl, group = array, color = array)) + ylab("Cummulative Variance explained") +
    geom_point() + geom_line()+ theme(axis.text.x=element_text(angle=45, hjust=1))
  
  pdf(paste(path,"pre_imputation/QC/pca_1000gp/PC_for_analysis_per_array_",ancestry[jj],"_2022.pdf",sep=""),width = 5, height = 6)
  print(p1)
  dev.off()
  
  
  
}

# ancestry="eur_all"
# ancestry="eur_nonjewish"
# ancestry="eur_jewish"
# 
# # 
#   
  
#######################################################################################################################################
#######################################################################################################################################

# double check Ash definitions:

# GSA combined genotype file - see step 27


path_gwas=/path/to/ibdgwas/IIBDGC/
i=(gsa)
ancestry=(eur_nonjewish eur_jewish)


for j in ${ancestry[@]} 
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}post_imputation/2022/${i}/genotyped_data/${i}_all_studies_merged_eur_pruned \
--keep ${path_gwas}pre_imputation/QC/pca_1000gp/list_${j}_ancestry_samples_noDuplicates_perstudy_interstudy_2 \
--keep-allele-order --allow-no-sex \
--freq counts --out ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned${j}_ancestry_samples_noDuplicates_perstudy_interstudy_2_maf
done



# --keep: 54503 people remaining.
# --freq: Allele frequencies (founders only) written to
# /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_prunedeur_nonjewish_ancestry_samples_noDuplicates_perstudy_interstudy_2_maf.frq.counts

# --keep: 5006 people remaining.
# --freq: Allele frequencies (founders only) written to
# /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_prunedeur_jewish_ancestry_samples_noDuplicates_perstudy_interstudy_2_maf.frq.counts


# GNOMAD:

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=200


bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_b38_alleles_5_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_b38_alleles_5_${i} \
"awk 'NR==FNR{vals[\$2];next} (\$1) in vals' ${path_gwas}post_imputation/2022/${i}/genotyped_data/${i}_all_studies_merged_eur_pruned.bim \
<(zcat ${path_gwas}resources/gnomad/gnomad_freq_edited.gz) > ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_gnomad_variants"


###########

##### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(viridis)
library(ggsci)
library(lmtest)

path<-"/path/to/ibdgwas/IIBDGC/"

array<-c("gsa")
i=1
ancestry<-c("eur_jewish","eur_nonjewish")

gnomad<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/",array[i],"_all_studies_merged_eur_pruned_gnomad_variants",sep=""),head=F)
colnames(gnomad)<-c("SNP","CHROM","POS","REF","ALT","AF","AF_nfe","AF_afr","AF_amr","AF_eas","AF_sas","AF_asj")
gnomad$AF_nfe<-as.numeric(as.character(gnomad$AF_nfe))
gnomad$AF_asj<-as.numeric(as.character(gnomad$AF_asj))

# identify set of variants with highest discrepancies between ASJ and NFE:


gnomad$value<-((gnomad$AF_nfe-gnomad$AF_asj)^2)/
  ((gnomad$AF_nfe+gnomad$AF_asj)*
     (2-gnomad$AF_nfe-gnomad$AF_asj))

dim(gnomad[which(gnomad$value>0.025),])
# [1] 3015  13

keep<-gnomad[which(gnomad$value>0.05),]

value<-0.025


gnomad$col[which(gnomad$value<=value)]<-paste("<=",value,sep="")
gnomad$col[which(gnomad$value>value)]<-paste(">",value,sep="")
gnomad$col[which(is.na(gnomad$value))]<-"monomorphic"


cbPalette <- c("#999999", "#E69F00", "#56B4E9")


p1<-ggplot(gnomad, aes(x=AF_asj, y=AF_nfe)) +
  geom_point(aes(colour = col)) + xlab("ASJ FRQ Alt") + ylab("NFE FRQ Alt") + scale_colour_manual(values=cbPalette) +
  ggtitle("Gnomad nfe vs asj")

pdf(paste(path,"pre_imputation/QC/pca_1000gp/gnomad_maf_nfe_asj.pdf",sep=""),width = 6, height = 5)
print(p1)
dev.off()



for (j in 1:length(ancestry)) {
  
  tmp<-fread(paste(path,"pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_pruned",ancestry[j],"_ancestry_samples_noDuplicates_perstudy_interstudy_2_maf.frq.counts",sep=""),head=T)
  
  tmp<-tmp[which(tmp$SNP %in% gnomad$SNP),]
  tmp<-tmp[match(gnomad$SNP,tmp$SNP)]
  print(table(tmp$SNP==gnomad$SNP))
  
  tmp$freq<-NA
  tmp$freq[which(tmp$A1==gnomad$ALT)]<-tmp$C1[which(tmp$A1==gnomad$ALT)]/(tmp$C1[which(tmp$A1==gnomad$ALT)]+tmp$C2[which(tmp$A1==gnomad$ALT)])
  tmp$freq[which(tmp$A2==gnomad$ALT)]<-1-(tmp$C1[which(tmp$A2==gnomad$ALT)]/(tmp$C1[which(tmp$A2==gnomad$ALT)]+tmp$C2[which(tmp$A2==gnomad$ALT)]))
  
  tmp[,paste(ancestry[j],"freq_alt",sep="_")]<-tmp$freq
  tmp<-tmp[,c(2,9)]
  
  gnomad<-merge(gnomad,tmp,by.x="SNP",by.y="SNP")
  
}


head(gnomad)


p1<-ggplot(gnomad[which(gnomad$value>0.025),], aes(y=AF_nfe, x=eur_nonjewish_freq_alt)) +
  geom_point() + ylab("NFE FRQ Alt") + xlab("IIBDGC EUR non-Jewish") + theme(panel.background = element_blank())

p2<-ggplot(gnomad[which(gnomad$value>0.025),], aes(y=AF_nfe, x=eur_jewish_freq_alt)) +
  geom_point() + ylab("NFE FRQ Alt") + xlab("IIBDGC EUR Jewish") + theme(panel.background = element_blank())

p3<-ggplot(gnomad[which(gnomad$value>0.025),], aes(y=AF_asj, x=eur_nonjewish_freq_alt)) +
  geom_point() + ylab("ASJ FRQ Alt") + xlab("IIBDGC EUR non-Jewish") + theme(panel.background = element_blank())

p4<-ggplot(gnomad[which(gnomad$value>0.025),], aes(y=AF_asj, x=eur_jewish_freq_alt)) +
  geom_point() + ylab("ASJ FRQ Alt") + xlab("IIBDGC EUR Jewish") + theme(panel.background = element_blank())


pdf(paste(path,"pre_imputation/QC/pca_1000gp/gnomad_maf_nfe_asj_iibdgc_eur_jewish_nonjewish.pdf",sep=""),width = 8, height = 7)
ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
dev.off()




gnomad$value2<-((gnomad$eur_nonjewish_freq_alt-gnomad$eur_jewish_freq_alt)^2)/
  ((gnomad$eur_nonjewish_freq_alt+gnomad$eur_jewish_freq_alt)*
     (2-gnomad$eur_nonjewish_freq_alt-gnomad$eur_jewish_freq_alt))




p1<-ggplot(gnomad[which(gnomad$value>0.025),], aes(x=eur_jewish_freq_alt, y=eur_nonjewish_freq_alt)) +
  geom_point() + xlab("IIBDGC EUR Jewish") + ylab("IIBDGC EUR non-Jewish") + scale_colour_manual(values=cbPalette) +
  ggtitle("IIBDGC EUR-nonJewish vs EUR-Jewish")

pdf(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_eur_jewish_nonjewish.pdf",sep=""),width = 5, height = 5)
print(p1)
dev.off()






