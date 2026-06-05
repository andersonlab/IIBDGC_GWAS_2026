# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# RUN PCA IN EUROPEAN ANCESTRY SAMPLES, ESTIMATE PCS IN NON DUPLICATED EUR SUBSET, AND PROJECT ALL

# TO DEFINE POPULATION ANCESTRY PCS TO INCLUDE IN THE ANALYSIS (ARRAY BY ARRAY)

###########################################################################################################################################################################################

########################################################
# 0.0 DEFINE A SET OF NON-RELATED SAMPLES - 3rd degree #
########################################################

# START WITH EUR POST-QC NON PRUNED DATA, see script pre-regenie

array=(illumina370 illumina550 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 gsa illuminaexome)
path_gwas=/path/to/ibdgwas/IIBDGC/
 
for i in ${array[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}post_imputation/${i}/genotyped_data/${i}_all_studies_merged_eur \
--keep ${path_gwas}pre_imputation/QC/pca_1000gp/list_eur_ancestry_samples_noDuplicates \
--autosome \
--maf 0.001 \
--make-bed --out ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_noDuplicates
done

wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/*_all_studies_merged_eur_noDuplicates.fam
# 4637 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix500_all_studies_merged_eur_noDuplicates.fam
# 8238 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix6_all_studies_merged_eur_noDuplicates.fam
# 46860 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_noDuplicates.fam
# 21238 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humancoreexome_all_studies_merged_eur_noDuplicates.fam
# 2010 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomni1_all_studies_merged_eur_noDuplicates.fam
# 1236 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomniexpress_all_studies_merged_eur_noDuplicates.fam
# 3706 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina370_all_studies_merged_eur_noDuplicates.fam
# 7749 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina550_all_studies_merged_eur_noDuplicates.fam
# 635 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illuminaexome_all_studies_merged_eur_noDuplicates.fam
# 3400 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/quad610_all_studies_merged_eur_noDuplicates.fam
# 99709 total

wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/*_all_studies_merged_eur_noDuplicates.bim
# 398870 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix500_all_studies_merged_eur_noDuplicates.bim
# 606323 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix6_all_studies_merged_eur_noDuplicates.bim
# 405483 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_noDuplicates.bim
# 235846 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humancoreexome_all_studies_merged_eur_noDuplicates.bim
# 786196 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomni1_all_studies_merged_eur_noDuplicates.bim
# 632211 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomniexpress_all_studies_merged_eur_noDuplicates.bim
# 261281 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina370_all_studies_merged_eur_noDuplicates.bim
# 473207 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina550_all_studies_merged_eur_noDuplicates.bim
# 227289 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illuminaexome_all_studies_merged_eur_noDuplicates.bim
# 548911 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/quad610_all_studies_merged_eur_noDuplicates.bim

# The amount of computer memory required by KING analysis is modest, at ~N ✕ M / 4 (where N is the number of samples and M is the number of SNPs)
# plus a small percentage of overhead cost. E.g., for a dataset consisting of 100,000 samples each genotyped at 1,000,000 SNPs, 
# the required memory size is ~25GB.

# Max estimated with GSA
# (413894*46860/4)/1E+6
# [1] 4848.768

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=5000

# accuracy up to 3rd- or 4th-degree (depending on array or WGS) for --related and --ibdseg analyses, and up to 2nd-degree for --kinship analysis

array=(illumina370 illumina550 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 gsa illuminaexome)

for i in ${array[@]}
do
bsub -J"kg" -M"$MEM" -n8 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal -n 8 \
-e ${path_gwas}pre_imputation/QC/all_hce/logs/stderr_king_${i} \
-o ${path_gwas}pre_imputation/QC/all_hce/logs/stdout_king_${i} \
"/path/to/software/./king \
-b ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_noDuplicates.bed \
--related --degree 3 --cpus 8 --prefix ${path_gwas}pre_imputation/QC/pca_1000gp/king_relatedness_eur_noDuplicates_${i}"
done


#### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(viridis)
library(ggsci)

path<-"/path/to/ibdgwas/IIBDGC/"
array<-c("illumina370","illumina550","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")

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
#                    2nd   3rd    FS    PO    UN   4th
# affymetrix500        6     9     0     0     0     0
# affymetrix6          4     9     7    55     0     0
# gsa                259   320   407   812 52742     5
# humancoreexome      92   112   163   869 12679     6
# humanomni1           3    12    16     2     1     0
# humanomniexpress     1     1     1     2     0     0
# illumina370          1     6     3     1    18     0
# illumina550         10    39    65    13     0     0
# illuminaexome        1     2     4     6     0     0
# quad610              1    10     0     0     1     0


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
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/kinship_vs_propIBD_eur_nonDuplicates.pdf ~/tmp_plots/",sep=""))


pdf(paste(path,"pre_imputation/QC/pca_1000gp/kinship_histogram_eur_nonDuplicates.pdf",sep=""),width = 14, height = 7)
ggplot(all, aes(x=Kinship, color=array,fill=array)) +
  geom_density(alpha=0.4)
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/kinship_histogram_eur_nonDuplicates.pdf ~/tmp_plots/",sep=""))





rel_ids<-c(as.character(all$FID1),as.character(all$FID2))
length(rel_ids)
table(length(rel_ids)==(nrow(all)*2))
rel_ids<-rel_ids[!duplicated(rel_ids)]
length(rel_ids)
# [1] 45986

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
# [1] 45986     1

data_remove<-data_remove[!duplicated(data_remove$V1),,drop=F]
dim(data_remove)
# [1] 25298     1

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
#            15               72            20439             4582 
# humanomni1 humanomniexpress      illumina370      illumina550 
#         32                5               17              111 
# illuminaexome          quad610 
#            13               12 

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

# Uudate MHC based on https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38.p11
# MHC -- chr6 (NC_000006.12):28,510,120-33,480,577
# 28510120-500000 
# 33480577+500000
# 
# > 28510120-500000 
# [1] 28010120
# > 33480577+500000
# [1] 33980577


###########################################################################################################################################################################################

################################################################
# 2.0 EXCLUDE VARIANTS IN HIGH LD REGIONS, KEEP ONLY AUTOSOMAL #
################################################################

# source files, same as used in regenie stage1
array=(illumina370 illumina550 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 gsa illuminaexome)

for i in ${array[@]}
do
wc -l ${path_gwas}post_imputation/${i}/genotyped_data/${i}_all_studies_merged_eur_pruned.bim
done
# 191724 /path/to/ibdgwas/IIBDGC/post_imputation/illumina370/genotyped_data/illumina370_all_studies_merged_eur_pruned.bim
# 276246 /path/to/ibdgwas/IIBDGC/post_imputation/illumina550/genotyped_data/illumina550_all_studies_merged_eur_pruned.bim
# 265610 /path/to/ibdgwas/IIBDGC/post_imputation/affymetrix6/genotyped_data/affymetrix6_all_studies_merged_eur_pruned.bim
# 312099 /path/to/ibdgwas/IIBDGC/post_imputation/humanomniexpress/genotyped_data/humanomniexpress_all_studies_merged_eur_pruned.bim
# 192476 /path/to/ibdgwas/IIBDGC/post_imputation/affymetrix500/genotyped_data/affymetrix500_all_studies_merged_eur_pruned.bim
# 164254 /path/to/ibdgwas/IIBDGC/post_imputation/humancoreexome/genotyped_data/humancoreexome_all_studies_merged_eur_pruned.bim
# 352329 /path/to/ibdgwas/IIBDGC/post_imputation/humanomni1/genotyped_data/humanomni1_all_studies_merged_eur_pruned.bim
# 299132 /path/to/ibdgwas/IIBDGC/post_imputation/quad610/genotyped_data/quad610_all_studies_merged_eur_pruned.bim
# 329528 /path/to/ibdgwas/IIBDGC/post_imputation/gsa/genotyped_data/gsa_all_studies_merged_eur_pruned.bim
# 166704 /path/to/ibdgwas/IIBDGC/post_imputation/illuminaexome/genotyped_data/illuminaexome_all_studies_merged_eur_pruned.bim

for i in ${array[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}post_imputation/${i}/genotyped_data/${i}_all_studies_merged_eur_pruned \
--autosome \
--allow-no-sex \
--exclude range ${path_gwas}pre_imputation/QC/pca_1000gp/high-LD-regions-hg38-GRCh38.txt \
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
# 409 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina370_all_studies_merged_eur_list_associated_variants_toexclude
# 119 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina550_all_studies_merged_eur_list_associated_variants_toexclude
# 181 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix6_all_studies_merged_eur_list_associated_variants_toexclude
# 52 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomniexpress_all_studies_merged_eur_list_associated_variants_toexclude
# 86 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix500_all_studies_merged_eur_list_associated_variants_toexclude
# 361 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humancoreexome_all_studies_merged_eur_list_associated_variants_toexclude
# 83 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomni1_all_studies_merged_eur_list_associated_variants_toexclude
# 57 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/quad610_all_studies_merged_eur_list_associated_variants_toexclude
# 50089 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_list_associated_variants_toexclude
# 19 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illuminaexome_all_studies_merged_eur_list_associated_variants_toexclude

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
# 188410 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina370_all_studies_merged_eur_pruned_noHighLD_noassoclogP4.bim
# 265475 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina550_all_studies_merged_eur_pruned_noHighLD_noassoclogP4.bim
# 253600 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix6_all_studies_merged_eur_pruned_noHighLD_noassoclogP4.bim
# 298935 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomniexpress_all_studies_merged_eur_pruned_noHighLD_noassoclogP4.bim
# 185361 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix500_all_studies_merged_eur_pruned_noHighLD_noassoclogP4.bim
# 157207 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humancoreexome_all_studies_merged_eur_pruned_noHighLD_noassoclogP4.bim
# 336910 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomni1_all_studies_merged_eur_pruned_noHighLD_noassoclogP4.bim
# 286913 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/quad610_all_studies_merged_eur_pruned_noHighLD_noassoclogP4.bim
# 267397 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_pruned_noHighLD_noassoclogP4.bim
# 160037 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illuminaexome_all_studies_merged_eur_pruned_noHighLD_noassoclogP4.bim


###########################################################################################################################################################################################

#######################################################
# 4.0 ESTIMATE PC, EXPORT PCS, AND ALLELE FREQUENCIES #
#######################################################

# No duplicates or related

for i in ${array[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_noHighLD_noassoclogP4 \
--keep ${path_gwas}pre_imputation/QC/pca_1000gp/list_eur_ancestry_samples_noDuplicates \
--make-bed --out ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_noDuplicates
done

wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/*_all_studies_merged_eur_pruned_noDuplicates.fam
# 4637 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix500_all_studies_merged_eur_pruned_noDuplicates.fam
# 8238 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix6_all_studies_merged_eur_pruned_noDuplicates.fam
# 46860 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_pruned_noDuplicates.fam
# 21238 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humancoreexome_all_studies_merged_eur_pruned_noDuplicates.fam
# 2010 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomni1_all_studies_merged_eur_pruned_noDuplicates.fam
# 1236 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomniexpress_all_studies_merged_eur_pruned_noDuplicates.fam
# 3706 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina370_all_studies_merged_eur_pruned_noDuplicates.fam
# 7749 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina550_all_studies_merged_eur_pruned_noDuplicates.fam
# 635 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illuminaexome_all_studies_merged_eur_pruned_noDuplicates.fam
# 3400 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/quad610_all_studies_merged_eur_pruned_noDuplicates.fam
# 99709 total

for i in ${array[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_noDuplicates \
--remove ${path_gwas}pre_imputation/QC/pca_1000gp/list_eur_related_3rd_degree_samples \
--make-bed --out ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_noDuplicates_noRelated3rd
done

wc -l ${path_gwas}pre_imputation/QC/pca_1000gp/*_all_studies_merged_eur_pruned_noDuplicates_noRelated3rd.fam
# 4622 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix500_all_studies_merged_eur_pruned_noDuplicates_noRelated3rd.fam
# 8166 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/affymetrix6_all_studies_merged_eur_pruned_noDuplicates_noRelated3rd.fam
# 26421 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/gsa_all_studies_merged_eur_pruned_noDuplicates_noRelated3rd.fam
# 16656 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humancoreexome_all_studies_merged_eur_pruned_noDuplicates_noRelated3rd.fam
# 1978 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomni1_all_studies_merged_eur_pruned_noDuplicates_noRelated3rd.fam
# 1231 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/humanomniexpress_all_studies_merged_eur_pruned_noDuplicates_noRelated3rd.fam
# 3689 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina370_all_studies_merged_eur_pruned_noDuplicates_noRelated3rd.fam
# 7638 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illumina550_all_studies_merged_eur_pruned_noDuplicates_noRelated3rd.fam
# 622 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/illuminaexome_all_studies_merged_eur_pruned_noDuplicates_noRelated3rd.fam
# 3388 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pca_1000gp/quad610_all_studies_merged_eur_pruned_noDuplicates_noRelated3rd.fam
# 74411 total


# all samples
for i in ${array[@]}
do
/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_noHighLD_noassoclogP4 \
--keep ${path_gwas}pre_imputation/QC/pca_1000gp/list_eur_ancestry_samples_withDuplicates \
--make-bed --out ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_withDuplicates
done


# The following command exports PCs to project onto, along with the allele frequencies needed to calibrate the 'variance-standardize' operation:

array=(illumina370 illumina550 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 gsa illuminaexome)
path_gwas=/path/to/ibdgwas/IIBDGC/
  
MEM=20000

for i in ${array[@]}
do
bsub -J"bcftools" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal -n 8 \
-e ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stderr_${i}_pca_eur_nodup_pcs_for_analysis \
-o ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stdout_${i}_pca_eur_nodup_pcs_for_analysis \
"/path/to/software/./plink2  \
--bfile ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_noDuplicates_noRelated3rd \
--freq counts \
--pca '40' approx allele-wts \
--out ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_ref_pcs \
--threads 8 --memory $MEM"
done
# Job <501938..501947> is submitted to queue <normal>.




#########################################################
## 1.2 - PROJECT ONTO THOSE PCS ALL IIBDGC

array=(illumina370 illumina550 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 gsa illuminaexome)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=2000 # next time only this

for i in ${array[@]}
do
bsub -J"bcftools" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal -n 10 \
-e ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stderr_${i}_pca_eur_nodup_project_withDuplicates_pcs_for_analysis \
-o ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stdout_${i}_pca_eur_nodup_project_withDuplicates_pcs_for_analysis \
"/path/to/software/./plink2 \
--bfile ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_withDuplicates \
--autosome \
--read-freq ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_ref_pcs.acount \
--score ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_ref_pcs.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
--score-col-nums 6-45 \
--out ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_withDuplicates_ref_pcs_new_projection"
done
# Job <503949..503958> is submitted to queue <normal>.


###############################

### /software/R-4.3.1/bin/R


library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(viridis)
library(ggsci)

path<-"/path/to/ibdgwas/IIBDGC/"

# no german_illu
cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"chop_old_gwas","german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas"
           ,"swedish_uc_old_gwas","mccauley_gsa","ccfa_gsa","cedars_gsa")

length(cohorts)
# [1] 34

array<-c("illumina370","illumina550","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")


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
# [1] 117923      7

# add country of origin
fam$country[which(fam$cohort=="australia_omniexome")]<-"Australia"
fam$country[which(fam$cohort %in% c("gwas1","gwas2","all_hce"))]<-"UK"
fam$country[which(fam$cohort %in% c("spain_gsa","basque_gsa"))]<-"Spain"
fam$country[which(fam$cohort %in% c("kiel_austria_sibdcs_gsa","german_affy6_old_gwas"))]<-"Germany"                 
fam$country[which(fam$cohort %in% c("italy_gsa"))]<-"Italy"   
fam$country[which(fam$cohort %in% c("netherlands_gsa"))]<-"Netherlands"                     
fam$country[which(fam$cohort %in% c("slovenia_gsa"))]<-"Slovenia"                     
fam$country[which(fam$cohort %in% c("sweden_gsa","swedish_uc_old_gwas"))]<-"Sweden"                    
fam$country[which(fam$cohort %in% c("finland_illugwas"))]<-"Finland"
fam$country[which(fam$cohort %in% c("chop_old_gwas"))]<-"CHOP"
fam$country[which(fam$cohort %in% c("norway_affy6_old_gwas"))]<-"Norway"
fam$country[which(fam$cohort %in% c("belgium_inf1_old_gwas","belgium_inf2_old_gwas","belgium_louis_gsa","belgium_franchimont_gsa"
                                    ,"belgium_vermeire_gsa"))]<-"Belgium"
fam$country[which(fam$cohort %in% c("lithuania_gsa"))]<-"Lithuania"
fam$country[which(fam$cohort %in% c("prism_nfe_gwas","prism_nfe_gsa","cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","cedars_gsa",
                                    "niddk_broad_gsa","niddk_feinstein_gsa","niddk_old_gwas","pittsburgh_gsa","mccauley_gsa","ccfa_gsa"))]<-"USA"

table(fam$country,useNA="ifany")
fam$country<-as.factor(fam$country)
fam$country<-factor(fam$country, levels=c("Spain","Italy","Slovenia","Belgium","Germany","Netherlands","UK","Lithuania",
                                          "Sweden","Norway","Finland","Australia","CHOP","USA"))


# Add jewish ancestry
pheno<-read.csv("/path/to/ibdgwas/IIBDGC/pheno/gwas-mega-2-core-phenotypes.csv",head=T)
colnames(pheno)[14]<-"self_jewish"
pheno_cd<-fread(paste(path,"pheno/list_cedars_old_gwas_ids_TH.txt",sep=""),head=T,sep="\t")
pheno<-rbind(pheno[,c("sample_id","self_jewish")],pheno_cd[,c("sample_id","self_jewish")])
pheno$self_jewish<-as.character(pheno$self_jewish)
pheno$self_jewish[which(pheno$self_jewish %in% c("","Unknown"))]<-"Unknown"
pheno$self_jewish[which(pheno$self_jewish %in% c("No","Non-Jewish"))]<-"Non-Jewish"
pheno$self_jewish[which(pheno$self_jewish %in% c("Yes","Jewish"))]<-"Jewish"
table(pheno$self_jewish)
# Jewish Non-Jewish    Unknown 
# 5763      35488      39572 


fam$self_jewish<-"Unknown"
fam$self_jewish[which(fam$V1 %in% pheno$sample_id[which(pheno$self_jewish=="Non-Jewish")])]<-"Non-Jewish"
fam$self_jewish[which(fam$V1 %in% pheno$sample_id[which(pheno$self_jewish=="Jewish")])]<-"Jewish"
table(fam$self_jewish)
# Jewish Non-Jewish    Unknown 
# 5573      32381      79969 


tmp<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/iibdgc_EUR_withDuplicates_ref_pcs_new_projection_edited.sscore",sep=""),head=T,sep="\t")
fam$pca_jewish<-NA
fam$pca_jewish[which(fam$V1 %in% tmp$FID[which(tmp$pca_jewish=="Non-Jewish")])]<-"Non-Jewish"
fam$pca_jewish[which(fam$V1 %in% tmp$FID[which(tmp$pca_jewish=="Jewish")])]<-"Jewish"
rm(tmp)


# add pheno CD/UC/IBDu/Ctr:
pheno<-read.table(paste(path,"pheno/phenotype_data_for_all_cohorts_Dec2020_analysis.tsv",sep=""),head=T)

fam$pheno<-NA
fam$pheno[which(fam$V1 %in% pheno$FID[which(pheno$ibd==0)])]<-"Control"
fam$pheno[which(fam$V1 %in% pheno$FID[which(pheno$cd==0)])]<-"Control"
fam$pheno[which(fam$V1 %in% pheno$FID[which(pheno$uc==0)])]<-"Control"

fam$pheno[which(fam$V1 %in% pheno$FID[which(pheno$ibd==1)])]<-"IBDu"
fam$pheno[which(fam$V1 %in% pheno$FID[which(pheno$cd==1)])]<-"CD"
fam$pheno[which(fam$V1 %in% pheno$FID[which(pheno$uc==1)])]<-"UC"

fam$pheno<-factor(fam$pheno, levels=c("CD","IBDu","UC","Control"))

# array by array expore the PCS and estimate how many to include in the analysis

###############################################
# 1.- TOTAL VARIANCE EXPLAINED BY EACH PC:

for (i in 1:length(array)) {
  
  eigenval<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/",array[i],"_all_studies_merged_eur_pruned_ref_pcs.eigenval",sep=""),head=F)
  
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

pdf(paste(path,"pre_imputation/QC/pca_1000gp/cummulative_variance_explained_40PCs_for_analysis_per_array.pdf",sep=""),width = 10, height = 6)
print(p1)
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/cummulative_variance_explained_40PCs_for_analysis_per_array.pdf ~/tmp_plots/",sep=""))


###############################################
# 2.- EXPLORE EACH PC FOR:


data<-matrix(nrow=length(array),ncol=11)
data<-as.data.frame(data)
colnames(data)<-c("array","m0_vs_m1","m1_vs_m2","m2_vs_m3","m3_vs_m4","m4_vs_m5","m5_vs_m6","m6_vs_m7","m7_vs_m8","m8_vs_m9","m9_vs_m10")

for (ii in 1:length(array)) {
  
  print(array[ii])
  pca<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/",array[ii],"_all_studies_merged_eur_pruned_withDuplicates_ref_pcs_new_projection.sscore",sep=""),head=F)
  colnames(pca)<-c("FID","IID","PHENO1","ALLELE_CT","NAMED_ALLELE_DOSAGE_SUM",paste("PC",seq(1:40),sep=""))
  pca<-merge(pca,fam[,c("V1","cohort","country","self_jewish","pca_jewish","pheno")],by.x="FID",by.y="V1",all.x=T)
  
  # pca$country<-as.character(pca$country)
  # pca$pheno<-as.character(pca$pheno)
  # pca$self_jewish<-as.character(pca$self_jewish)
  
  # In Australia cohort only 1 self reported jewish, density needs at least 2 points
  if(ii==4) {
    pca$self_jewish<-"Unknown"
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
    
                    
    pdf(paste(path,"pre_imputation/QC/pca_1000gp/histogram_european_ancestry_samples_PC",j,"_per_array_",array[ii],".pdf",sep=""),width=14,height=25)
    print(annotate_figure(fig,top = text_grob(array[ii], face = "bold",size=14)))
    dev.off()
    
    system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/histogram_european_ancestry_samples_PC",j,"_per_array_",array[ii],".pdf ~/tmp_plots/",sep=""))
    
  }
  
  # create association report with phenotype, only nondup or non related samples:
  
  fam_analysis<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/",array[ii],"_all_studies_merged_eur_pruned_noDuplicates_noRelated3rd.fam",sep=""),head=F)
  
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

write.table(data,paste(path,"pre_imputation/QC/pca_1000gp/PC_selection_all_studies_perarray",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/PC_selection_all_studies_perarray ~/tmp_plots/",sep=""))


for (ii in 1:length(array)) {
  print(array[ii])
  print(all_eigen[which(all_eigen$array==array[ii] & all_eigen$var_exp>=3),])
}



p1<-ggplot(all_eigen, aes(x = PC, y=cumsum_var_expl, group = array, color = array)) + ylab("Cummulative Variance explained") +
  geom_point() + geom_line()+ theme(axis.text.x=element_text(angle=45, hjust=1))

pdf(paste(path,"pre_imputation/QC/pca_1000gp/cummulative_variance_explained_40PCs_for_analysis_per_array.pdf",sep=""),width = 10, height = 6)
print(p1)
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/cummulative_variance_explained_40PCs_for_analysis_per_array.pdf ~/tmp_plots/",sep=""))

p1<-ggplot(all_eigen, aes(x = PC, y=var_exp, group = array, color = array)) + ylab("Variance explained") +
  geom_point() + geom_line()+ theme(axis.text.x=element_text(angle=45, hjust=1)) 
# +scale_y_reverse()

pdf(paste(path,"pre_imputation/QC/pca_1000gp/variance_explained_40PCs_for_analysis_per_array.pdf",sep=""),width = 10, height = 6)
print(p1)
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/variance_explained_40PCs_for_analysis_per_array.pdf ~/tmp_plots/",sep=""))

tmp<-all_eigen[which( (all_eigen$array=="illumina370" & all_eigen$PC %in% paste("PC",seq(1:5),sep="")) |
                        (all_eigen$array=="illumina550" & all_eigen$PC %in% paste("PC",seq(1:2),sep="")) |
                        (all_eigen$array=="affymetrix6" & all_eigen$PC %in% paste("PC",seq(1),sep="")) |
                        (all_eigen$array=="humanomniexpress" & all_eigen$PC %in% paste("PC",seq(1),sep="")) |
                        (all_eigen$array=="affymetrix500" & all_eigen$PC %in% paste("PC",seq(1:3),sep="")) |
                        (all_eigen$array=="humancoreexome" & all_eigen$PC %in% paste("PC",seq(1:2),sep="")) |
                        (all_eigen$array=="humanomni1" & all_eigen$PC %in% paste("PC",seq(1:3),sep="")) |
                        (all_eigen$array=="quad610" & all_eigen$PC %in% paste("PC",seq(1),sep="")) |
                        (all_eigen$array=="gsa" & all_eigen$PC %in% paste("PC",seq(1:4),sep="")) |
                        (all_eigen$array=="illuminaexome" & all_eigen$PC %in% paste("PC",seq(1:4),sep=""))
                        ),]
p1<-ggplot(tmp, aes(x = PC, y=cumsum_var_expl, group = array, color = array)) + ylab("Cummulative Variance explained") +
  geom_point() + geom_line()+ theme(axis.text.x=element_text(angle=45, hjust=1))

pdf(paste(path,"pre_imputation/QC/pca_1000gp/cummulative_variance_explained_5PCs_for_analysis_per_array.pdf",sep=""),width = 5, height = 6)
print(p1)
dev.off()
system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/cummulative_variance_explained_5PCs_for_analysis_per_array.pdf ~/tmp_plots/",sep=""))


####################################################################################################################################################################
####################################################################################################################################################################
#### estimate contribution of variants to different PC:


for (ii in 1:length(array)) {
  
  print(array[ii])
  pcs<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/",array[ii],"_all_studies_merged_eur_pruned_ref_pcs.eigenvec.allele",sep=""),head=F)
  colnames(pcs)<-c("Chr","ID","Ref","Alt","Allele1",paste("PC",seq(1:40),sep=""))
  
  p1<-ggplot(pcs, aes(x=PC1)) + geom_histogram(color="black") + xlab("PC1")
  p2<-ggplot(pcs, aes(x=PC2)) + geom_histogram(color="black") + xlab("PC2")
  p3<-ggplot(pcs, aes(x=PC3)) + geom_histogram(color="black") + xlab("PC3")
  p4<-ggplot(pcs, aes(x=PC4)) + geom_histogram(color="black") + xlab("PC4")
  p5<-ggplot(pcs, aes(x=PC5)) + geom_histogram(color="black") + xlab("PC5")
  p6<-ggplot(pcs, aes(x=PC6)) + geom_histogram(color="black") + xlab("PC6")
  p7<-ggplot(pcs, aes(x=PC7)) + geom_histogram(color="black") + xlab("PC7")
  p8<-ggplot(pcs, aes(x=PC8)) + geom_histogram(color="black") + xlab("PC8")
  p9<-ggplot(pcs, aes(x=PC9)) + geom_histogram(color="black") + xlab("PC9")
  
  
  p<-ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,ncol=3,nrow=3)
  
  pdf(paste(path,"pre_imputation/QC/pca_1000gp/weights_per_allele_distribution_PC1_PC10_",array[ii],".pdf",sep=""),width = 15, height = 10)
  print(annotate_figure(p,
                        top = text_grob(array[ii], color = "black", face = "bold", size = 10)))
  dev.off()
  system(paste("cp ",path,"pre_imputation/QC/pca_1000gp/weights_per_allele_distribution_PC1_PC10_",array[ii],".pdf ~/tmp_plots/",sep=""))
  rm(p1,p2,p3,p4,p5,p6,p7,p8,p9)
  
}


# ###########################################################################################################################################################################
# ###########################################################################################################################################################################
# 
# # issue with heterogeneity when illumina550 is meta-analysed with other cohorts, PCs are showing jewish ancestry, what would happen if we exclude all jewish ancestry samples:
# 
# #### /software/R-4.3.1/bin/R
# 
# library(data.table)
# library(ggplot2)
# library(ggpubr)
# library(gtools)
# library(viridis)
# library(ggsci)
# library(lmtest)
# 
# path<-"/path/to/ibdgwas/IIBDGC/"
# 
# # no german_illu
# cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
#            ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
#            ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
#            ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
#            ,"prism_nfe_gwas","finland_illugwas"
#            ,"chop_old_gwas","german_affy6_old_gwas","norway_affy6_old_gwas"
#            ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
#            ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas"
#            ,"swedish_uc_old_gwas","mccauley_gsa","ccfa_gsa","cedars_gsa")
# 
# length(cohorts)
# # [1] 34
# 
# array<-c("illumina370","illumina550","affymetrix6","humanomniexpress","affymetrix500","humancoreexome","humanomni1","quad610","gsa","illuminaexome")
# 
# 
# for (i in 1:length(cohorts)){
#   a<-read.table(paste(path,"pre_imputation/QC/relatedness/",cohorts[i],"_subset_withDuplicates.fam",sep=""),head=F)
#   a$cohort<-cohorts[i]
#   if(i==1){
#     fam<-a
#   }else{
#     fam<-rbind(a,fam)
#   }
# }
# 
# dim(fam)
# # [1] 117923      7
# 
# # add country of origin
# fam$country[which(fam$cohort=="australia_omniexome")]<-"Australia"
# fam$country[which(fam$cohort %in% c("gwas1","gwas2","all_hce"))]<-"UK"
# fam$country[which(fam$cohort %in% c("spain_gsa","basque_gsa"))]<-"Spain"
# fam$country[which(fam$cohort %in% c("kiel_austria_sibdcs_gsa","german_affy6_old_gwas"))]<-"Germany"                 
# fam$country[which(fam$cohort %in% c("italy_gsa"))]<-"Italy"   
# fam$country[which(fam$cohort %in% c("netherlands_gsa"))]<-"Netherlands"                     
# fam$country[which(fam$cohort %in% c("slovenia_gsa"))]<-"Slovenia"                     
# fam$country[which(fam$cohort %in% c("sweden_gsa","swedish_uc_old_gwas"))]<-"Sweden"                    
# fam$country[which(fam$cohort %in% c("finland_illugwas"))]<-"Finland"
# fam$country[which(fam$cohort %in% c("chop_old_gwas"))]<-"CHOP"
# fam$country[which(fam$cohort %in% c("norway_affy6_old_gwas"))]<-"Norway"
# fam$country[which(fam$cohort %in% c("belgium_inf1_old_gwas","belgium_inf2_old_gwas","belgium_louis_gsa","belgium_franchimont_gsa"
#                                     ,"belgium_vermeire_gsa"))]<-"Belgium"
# fam$country[which(fam$cohort %in% c("lithuania_gsa"))]<-"Lithuania"
# fam$country[which(fam$cohort %in% c("prism_nfe_gwas","prism_nfe_gsa","cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","cedars_gsa",
#                                     "niddk_broad_gsa","niddk_feinstein_gsa","niddk_old_gwas","pittsburgh_gsa","mccauley_gsa","ccfa_gsa"))]<-"USA"
# 
# table(fam$country,useNA="ifany")
# fam$country<-as.factor(fam$country)
# fam$country<-factor(fam$country, levels=c("Spain","Italy","Slovenia","Belgium","Germany","Netherlands","UK","Lithuania",
#                                           "Sweden","Norway","Finland","Australia","CHOP","USA"))
# 
# 
# fam$pheno<-fam$V6
# 
# ii=2
# array[ii]
# 
# print(array[ii])
# pca<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/",array[ii],"_all_studies_merged_eur_pruned_withDuplicates_ref_pcs_new_projection.sscore",sep=""),head=F)
# colnames(pca)<-c("FID","IID","PHENO1","ALLELE_CT","NAMED_ALLELE_DOSAGE_SUM",paste("PC",seq(1:40),sep=""))
# pca<-merge(pca,fam[,c("V1","cohort","pheno")],by.x="FID",by.y="V1",all.x=T)
# 
# table(pca$pheno)
# #    1    2 
# # 6086 2364 
# 
# pca$pheno[which(pca$pheno==1)]<-0
# pca$pheno[which(pca$pheno==2)]<-1
# table(pca$pheno)
# # 0    1 
# # 6086 2364
# 
# data<-matrix(nrow=length(array),ncol=11)
# data<-as.data.frame(data)
# colnames(data)<-c("array","m0_vs_m1","m1_vs_m2","m2_vs_m3","m3_vs_m4","m4_vs_m5","m5_vs_m6","m6_vs_m7","m7_vs_m8","m8_vs_m9","m9_vs_m10")
# 
# 
# m0_ibd<-glm(pheno ~ 1, data = pca, family = "binomial")
# m1_ibd<-glm(pheno ~ PC1, data = pca, family = "binomial")
# m2_ibd<-glm(pheno ~ PC1 + PC2, data = pca, family = "binomial")
# m3_ibd<-glm(pheno ~ PC1 + PC2 + PC3, data = pca, family = "binomial")
# 
# j=1
# data$array[j]<-"illumina550"
# data$m0_vs_m1[j]<-lrtest(m0_ibd,m1_ibd)[["Pr(>Chisq)"]][2]
# data$m1_vs_m2[j]<-lrtest(m1_ibd,m2_ibd)[["Pr(>Chisq)"]][2]
# data$m2_vs_m3[j]<-lrtest(m2_ibd,m3_ibd)[["Pr(>Chisq)"]][2]
# 
# # see chop script to find definition:
# keep<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/chop_list_samples_homogeneous_central_EUR_to_keep",sep=""),head=T)
# 
# pca_2<-pca[which( (pca$FID %in% keep$FID) & !is.na(pca$pheno)),]
# table(pca_2$pheno)
# # 0    1 
# # 2687 1200    
# 
# m0_ibd<-glm(pheno ~ 1, data = pca_2, family = "binomial")
# m1_ibd<-glm(pheno ~ PC1, data = pca_2, family = "binomial")
# m2_ibd<-glm(pheno ~ PC1 + PC2, data = pca_2, family = "binomial")
# m3_ibd<-glm(pheno ~ PC1 + PC2 + PC3, data = pca_2, family = "binomial")
# 
# j=2
# data$array[j]<-"illumina550_non_jewish"
# data$m0_vs_m1[j]<-lrtest(m0_ibd,m1_ibd)[["Pr(>Chisq)"]][2]
# data$m1_vs_m2[j]<-lrtest(m1_ibd,m2_ibd)[["Pr(>Chisq)"]][2]
# data$m2_vs_m3[j]<-lrtest(m2_ibd,m3_ibd)[["Pr(>Chisq)"]][2]
# 
# 
# data
# # array     m0_vs_m1    m1_vs_m2  m2_vs_m3 m3_vs_m4 m4_vs_m5
# # 1             illumina550 4.121380e-03 0.002187803 0.1614308       NA       NA
# # 2  illumina550_non_jewish 1.601468e-09 0.215874319 0.7341893       NA       NA
# 
# 
# # still include these two PCs but exclude jewish ancestry samples:
# 
# # check out plot:
# 
# pca$homogeneous_ancestry<-"no"
# pca$homogeneous_ancestry[which(pca$FID %in% keep$FID)]<-"yes"
# 
# p1<-qplot(PC2,PC1, data = pca, 
#           colour = pheno)  + facet_grid( ~ homogeneous_ancestry)
# 
# pdf("~/tmp_plots/pca_PC1_PC2_included_as_covar_in_analysis_for_subset_homogeneous_chop_samples.pdf",width=21,height=7)
# p1
# dev.off()

# #########################################
# 
# # generate new PCs with this subset of chop:
# 
# path_gwas=/path/to/ibdgwas/IIBDGC/
# 
# 
# /path/to/software/./plink2  \
# --bfile ${path_gwas}pre_imputation/QC/pca_1000gp/illumina550_all_studies_merged_eur_pruned_noDuplicates_noRelated3rd \
# --keep ${path_gwas}pre_imputation/QC/pca_1000gp/chop_list_samples_homogeneous_central_EUR_to_keep \
# --freq \
# --make-bed --out ${path_gwas}pre_imputation/QC/pca_1000gp/illumina550_all_studies_merged_eur_pruned_noDuplicates_noRelated3rd_homogeneous_chop_samples
# # 3114 samples (1486 females, 1603 males, 25 ambiguous; 3114 founders) remaining
# 
# MEM=4000
# i=illumina550
# bsub -J"bcftools" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal -n 8 \
# -e ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stderr_${i}_pca_eur_nodup_pcs_for_analysis \
# -o ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stdout_${i}_pca_eur_nodup_pcs_for_analysis \
# "/path/to/software/./plink2  \
# --bfile ${path_gwas}pre_imputation/QC/pca_1000gp/illumina550_all_studies_merged_eur_pruned_noDuplicates_noRelated3rd_homogeneous_chop_samples \
# --freq counts \
# --pca allele-wts \
# --out ${path_gwas}pre_imputation/QC/pca_1000gp/illumina550_subset_all_studies_merged_eur_pruned_ref_pcs \
# --threads 8 --memory $MEM"
# 
# 
# 
# #########################################################
# ## 1.2 - PROJECT ONTO THOSE PCS ALL IIBDGC
# 
# array=(illumina550)
# 
# path_gwas=/path/to/ibdgwas/IIBDGC/
# MEM=2000 # next time only this
# i=illumina550
# bsub -J"bcftools" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal -n 10 \
# -e ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stderr_${i}_pca_eur_nodup_project_withDuplicates_pcs_for_analysis \
# -o ${path_gwas}pre_imputation/QC/pca_1000gp/logs/stdout_${i}_pca_eur_nodup_project_withDuplicates_pcs_for_analysis \
# "/path/to/software/./plink2 \
# --bfile ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_withDuplicates \
# --autosome \
# --read-freq ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_ref_pcs.acount \
# --score ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_ref_pcs.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
# --score-col-nums 6-15 \
# --out ${path_gwas}pre_imputation/QC/pca_1000gp/${i}_all_studies_merged_eur_pruned_withDuplicates_ref_pcs_new_projection"
# done
# 
# ###################################################################################################################################################
# 
# ########## /software/R-4.3.1/bin/R
# 
# library(data.table)
# library(ggplot2)
# library(ggpubr)
# library(gtools)
# library(viridis)
# 
# path<-"/path/to/ibdgwas/IIBDGC/"
# 
# # TOTAL VARIANCE EXPLAINED BY EACH PC - in CHOP:
# eigenval<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/illumina550_subset_all_studies_merged_eur_pruned_ref_pcs.eigenval",sep=""),head=F)
# 
# eigenval$var_exp<-NA
# for (i in 1:nrow(eigenval)){
#   eigenval$var_exp[i]<-(eigenval$V1[i] / sum(eigenval$V1))*100
# }
# 
# eigenval
# #         V1   var_exp
# # 1  2.39580 14.020119
# # 2  1.89048 11.063008
# # 3  1.72393 10.088365
# # 4  1.70080  9.953009
# # 5  1.57161  9.196994
# # 6  1.56583  9.163170
# # 7  1.56261  9.144327
# # 8  1.56074  9.133384
# # 9  1.55877  9.121855
# # 10 1.55773  9.115769
# 
# pca<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/illumina550_all_studies_merged_eur_pruned_withDuplicates_ref_pcs_new_projection.sscore",sep=""),head=F)
# colnames(pca)<-c("FID","IID","PHENO1","ALLELE_CT","NAMED_ALLELE_DOSAGE_SUM","PC1_AVG","PC2_AVG","PC3_AVG","PC4_AVG","PC5_AVG","PC6_AVG"
#                  ,"PC7_AVG","PC8_AVG","PC9_AVG","PC10_AVG")
# dim(pca)
# # [1] 8450   15
# 
# ibd<-read.table(paste(path,"post_imputation/illumina550/phenotype_data/illumina550_all_studies_merged_eur_phenotype_ibd",sep=""),head=T)
# 
# pca<-merge(pca,ibd[,c("FID","ibd")],by="FID",all.x=T,sort=F)
# table(pca$PHENO1)
# # 1    2
# # 6086 2364
# table(pca$ibd)
# # 0    1
# # 5984 1662
# 
# colnames(pca)<-gsub("_AVG","",colnames(pca))
# 
# data<-matrix(nrow=2,ncol=11)
# data<-as.data.frame(data)
# colnames(data)<-c("array","m0_vs_m1","m1_vs_m2","m2_vs_m3","m3_vs_m4","m4_vs_m5","m5_vs_m6","m6_vs_m7","m7_vs_m8","m8_vs_m9","m9_vs_m10")
# 
# m0_ibd<-glm(ibd ~ 1, data = pca, family = "binomial")
# m1_ibd<-glm(ibd ~ PC1, data = pca, family = "binomial")
# m2_ibd<-glm(ibd ~ PC1 + PC2, data = pca, family = "binomial")
# m3_ibd<-glm(ibd ~ PC1 + PC2 + PC3, data = pca, family = "binomial")
# m4_ibd<-glm(ibd ~ PC1 + PC2 + PC3 + PC4, data = pca, family = "binomial")
# m5_ibd<-glm(ibd ~ PC1 + PC2 + PC3 + PC4 + PC5, data = pca, family = "binomial")
# 
# j=1
# data$array[j]<-"illumina550"
# data$m0_vs_m1[j]<-lrtest(m0_ibd,m1_ibd)[["Pr(>Chisq)"]][2]
# data$m1_vs_m2[j]<-lrtest(m1_ibd,m2_ibd)[["Pr(>Chisq)"]][2]
# data$m2_vs_m3[j]<-lrtest(m2_ibd,m3_ibd)[["Pr(>Chisq)"]][2]
# data$m3_vs_m4[j]<-lrtest(m3_ibd,m4_ibd)[["Pr(>Chisq)"]][2]
# data$m4_vs_m5[j]<-lrtest(m4_ibd,m5_ibd)[["Pr(>Chisq)"]][2]
# 
# # keep only subset of samples:
# # see chop script to find definition:
# keep<-read.table(paste(path,"pre_imputation/QC/pca_1000gp/chop_list_samples_homogeneous_central_EUR_to_keep",sep=""),head=T)
# pca<-pca[which(pca$FID %in% keep$FID),]
# 
# table(pca$ibd)
# # 0    1
# # 3675  933
# 
# m0_ibd<-glm(ibd ~ 1, data = pca, family = "binomial")
# m1_ibd<-glm(ibd ~ PC1, data = pca, family = "binomial")
# m2_ibd<-glm(ibd ~ PC1 + PC2, data = pca, family = "binomial")
# m3_ibd<-glm(ibd ~ PC1 + PC2 + PC3, data = pca, family = "binomial")
# m4_ibd<-glm(ibd ~ PC1 + PC2 + PC3 + PC4, data = pca, family = "binomial")
# m5_ibd<-glm(ibd ~ PC1 + PC2 + PC3 + PC4 + PC5, data = pca, family = "binomial")
# 
# j=2
# data$array[j]<-"illumina550"
# data$m0_vs_m1[j]<-lrtest(m0_ibd,m1_ibd)[["Pr(>Chisq)"]][2]
# data$m1_vs_m2[j]<-lrtest(m1_ibd,m2_ibd)[["Pr(>Chisq)"]][2]
# data$m2_vs_m3[j]<-lrtest(m2_ibd,m3_ibd)[["Pr(>Chisq)"]][2]
# data$m3_vs_m4[j]<-lrtest(m3_ibd,m4_ibd)[["Pr(>Chisq)"]][2]
# data$m4_vs_m5[j]<-lrtest(m4_ibd,m5_ibd)[["Pr(>Chisq)"]][2]
# 
# data
# # array     m0_vs_m1     m1_vs_m2   m2_vs_m3   m3_vs_m4   m4_vs_m5
# # 1 illumina550 1.469111e-05 4.094027e-14 0.97891111 0.07460532 0.05517124
# # 2 illumina550 8.026732e-21 1.018882e-04 0.04543215 0.00403076 0.09515849
# 
# # keep the 4 first PCs



