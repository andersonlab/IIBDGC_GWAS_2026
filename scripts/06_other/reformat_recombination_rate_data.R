# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# reformat recombination rate files

##### get additional LD-related continuous annotations as in https://pubmed.ncbi.nlm.nih.gov/28892061/
# Recombination rate, nucleotide diversity, GC content, and CpG content were computed using windows of different sizes around each SNP: ±10 kb
# The genetic positions of surrounding windows were interpolated linearly from recombination maps using PLINK. We determined that the Oxford map provided the most 
# significant results (Supplementary Table 3), suggesting that the impact of recombination rate on trait heritability operates over a long time scale; 
# we thus used the Oxford map in all primary analyses.
  
# 10KB recombination rate:
path_gwas=/path/to/ibdgwas/IIBDGC/
cd ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/
mkdir oxford_recombination_rate/
# wget http://www.shapeit.fr/files/genetic_map_b37.tar.gz data not available anymore at this location, 
# neither at http://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html/files/genetic_map_b37.tar.gz
# get from SHAPEIT4
wget https://github.com/odelaneau/shapeit4/raw/master/maps/genetic_maps.b37.tar.gz
  
# update files to b38
# transform into bed files:

for chr in {1..22}
do
zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b37/chr${chr}.b37.gmap.gz | \
awk -v OFS="\t" '{$10=$2+1;print $2, $1, $1+1, $3}' | awk 'NR>1 {print}' \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b37/chr${chr}_gmap.bed
done

chr=X
zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b37/chr${chr}_par1.b37.gmap.gz | \
awk -v OFS="\t" '{$10=$2+1;print $2, $1, $1+1, $3}' | awk 'NR>1 {print}' \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b37/chr${chr}_par1_gmap.bed
zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b37/chr${chr}_par2.b37.gmap.gz | \
awk -v OFS="\t" '{$10=$2+1;print $2, $1, $1+1, $3}' | awk 'NR>1 {print}' \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b37/chr${chr}_par2_gmap.bed

cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b37/chr${chr}_par1_gmap.bed \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b37/chr${chr}_par2_gmap.bed > \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b37/chr${chr}_gmap.bed

for chr in {1..22}
do
sed 's/^/chr/' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b37/chr${chr}_gmap.bed \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b37/chr${chr}_gmap_ed.bed
done

# rm intermed files:
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b37/chr*_gmap.bed

for chr in {1..22}
do
${path_gwas}previous_qced_b38/liftover/liftOver \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b37/chr${chr}_gmap_ed.bed \
${path_gwas}previous_qced_b38/liftover/hg19ToHg38.over.chain.gz \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chr${chr}_gmap_b38.bed \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b37/chr${chr}_gmap_nolifted
done

for chr in {1..22}
do 
echo ${chr} && cat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b37/chr${chr}_gmap_nolifted | sed "/^#/d" | wc -l
done

# # edit to match SHAPEIT format:
# # zcat ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b37/chr22.b37.gmap.gz | head -10
# # # pos	chr	cM
# # # 16051347	22	0.000000
# # # 16052618	22	0.010291
# # # 16053624	22	0.01847

for chr in {1..22}
do
sed 's/^chr//' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chr${chr}_gmap_b38.bed \
| awk -v OFS="\t" '{print $2, $1, $4}' > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chr${chr}.b38.gmap
done

# sort by position:
for chr in {1..22}
do 
sort -n -k 1,1 chr${chr}.b38.gmap | sed  '1i pos\tchr\tcM' > chr${chr}_sorted.b38.gmap
done


# rm intermed files:
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b37/chr*_ed.bed

# rm intermed files:
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chr*_gmap_b38.bed
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chr*_gmap_ed.bed
rm ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chr{1..22}.b38.gmap


sed 's/X/23/g' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chrX_sorted.b38.gmap \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chr23_sorted.b38.gmap


# combine all files into one
# use plink to assign cm/mb data to variants:
# 1.- generate 2 bim files from vcf with -/+ 10kb from each snp:

# chromosome (1-22, X, Y or 0 if unplaced)
# rs# or snp identifier
# Genetic distance (morgans)
# Base-pair position (bp units)

# source
# ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis.bed

sed 's/:\+/\t/g' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bed | \
awk -v OFS="\t" '{print $6,$7}' > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp


awk -v OFS="\t" '{$10="0";print $1, $4, $10, $2-10000}' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bed \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/list_union_variants_cd_uc_ibd_metaanalysis_minus10kb.map

awk -v OFS="\t" '{$10="0";print $1, $4, $10, $2+10000}' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bed \
> ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/list_union_variants_cd_uc_ibd_metaanalysis_plus10kb.map


paste -d'\t' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/list_union_variants_cd_uc_ibd_metaanalysis_minus10kb.map \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp> \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/allchr_list_union_variants_cd_uc_ibd_metaanalysis_minus10kb.bim

paste -d'\t' ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/list_union_variants_cd_uc_ibd_metaanalysis_plus10kb.map \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/tmp> \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/allchr_list_union_variants_cd_uc_ibd_metaanalysis_plus10kb.bim

head ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/allchr_list_union_variants_cd_uc_ibd_metaanalysis_minus10kb.bim
head ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/allchr_list_union_variants_cd_uc_ibd_metaanalysis_plus10kb.bim

# split into chr, and add cM to both files using plink:

MEM=4000
MEM1=100000
for chr in {1..22} X
do
bsub -J"annota" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM:tmp=$MEM1]  span[hosts=1]" -G your_hpc_group -n 4 \
-e ${path_gwas}post_imputation/log/${chr}_recomb_bim_stderr \
-o ${path_gwas}post_imputation/log/${chr}_recomb_bim_stdout \
"/path/to/software/username/plink_linux_x86_64_20181202/./plink \
--bim ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/allchr_list_union_variants_cd_uc_ibd_metaanalysis_minus10kb.bim \
--chr $chr \
--allow-no-samples \
--memory $MEM1 \
--threads 4 \
--make-just-bim --out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_minus10kb"
done


MEM=800
MEM1=800
for chr in {1..22}
do
bsub -J"annota" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM:tmp=$MEM1]  span[hosts=1]" -G your_hpc_group -n 4 \
-e ${path_gwas}post_imputation/log/${chr}_recomb_bim_2_stderr \
-o ${path_gwas}post_imputation/log/${chr}_recomb_bim_2_stdout \
"/path/to/software/username/plink_linux_x86_64_20181202/./plink \
--bim ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_minus10kb.bim \
--cm-map ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chr23_sorted.b38.gmap 23 \
--allow-no-samples \
--make-just-bim --out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_minus10kb_withCM"
done

chr=X
bsub -J"annota" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM:tmp=$MEM1]  span[hosts=1]" -G your_hpc_group -n 4 \
-e ${path_gwas}post_imputation/log/${chr}_recomb_bim_2_stderr \
-o ${path_gwas}post_imputation/log/${chr}_recomb_bim_2_stdout \
"/path/to/software/username/plink_linux_x86_64_20181202/./plink \
--bim ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_minus10kb.bim \
--cm-map ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chr${chr}_sorted.b38.gmap $chr \
--allow-no-samples \
--output-chr M \
--make-just-bim --out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_minus10kb_withCM"


MEM=4000
MEM1=100000
for chr in {1..22}
do
bsub -J"annota" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM:tmp=$MEM1]  span[hosts=1]" -G your_hpc_group -n 4 \
-e ${path_gwas}post_imputation/log/${chr}_recomb_bim_3_stderr \
-o ${path_gwas}post_imputation/log/${chr}_recomb_bim_3_stdout \
"/path/to/software/username/plink_linux_x86_64_20181202/./plink \
--bim ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/allchr_list_union_variants_cd_uc_ibd_metaanalysis_plus10kb.bim \
--chr $chr \
--allow-no-samples \
--memory $MEM1 \
--threads 4 \
--make-just-bim --out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_plus10kb"
done


MEM=800
MEM1=800
for chr in {1..22}
do
bsub -J"annota" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM:tmp=$MEM1]  span[hosts=1]" -G your_hpc_group -n 4 \
-e ${path_gwas}post_imputation/log/${chr}_recomb_bim_4_stderr \
-o ${path_gwas}post_imputation/log/${chr}_recomb_bim_4_stdout \
"/path/to/software/username/plink_linux_x86_64_20181202/./plink \
--bim ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_plus10kb.bim \
--cm-map ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chr${chr}_sorted.b38.gmap ${chr} \
--allow-no-samples \
--make-just-bim --out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_plus10kb_withCM"
done

chr=X
bsub -J"annota" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM:tmp=$MEM1]  span[hosts=1]" -G your_hpc_group -n 4 \
-e ${path_gwas}post_imputation/log/${chr}_recomb_bim_4_stderr \
-o ${path_gwas}post_imputation/log/${chr}_recomb_bim_4_stdout \
"/path/to/software/username/plink_linux_x86_64_20181202/./plink \
--bim ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_plus10kb.bim \
--cm-map ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chr23_sorted.b38.gmap 23 \
--allow-no-samples \
--output-chr M \
--make-just-bim --out ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chr${chr}_list_union_variants_cd_uc_ibd_metaanalysis_plus10kb_withCM"


# # ESTIMATE RECOMBINATION RATE:

# # https://www.biostars.org/p/222697/
# # cM_position2 - cM_position1 / ( position2(mb)-position1(mb) )

  

# MEM=8000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q yesterday /software/bin/R-4.3.1

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"
  
for (chr in c(1:22,"X") {
  
  minus<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chr",chr,"_list_union_variants_cd_uc_ibd_metaanalysis_minus10kb_withCM.bim",sep=""),head=F)
  
  # anything below first row in map has negative cM -> set to 0
  minus$V3[which(minus$V3<0)]<-0
  
  minus<-minus[,c(2:3)]
  minus<-as.data.frame(minus)
  
  plus<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/chr",chr,"_list_union_variants_cd_uc_ibd_metaanalysis_plus10kb_withCM.bim",sep=""),head=F)

  # anything below first row in map has negative cM -> set to 0
  plus$V3[which(plus$V3<0)]<-0
  
  plus<-plus[,c(2:3)]
  plus<-as.data.frame(plus)
  
  print(table(minus$V2==plus$V2))
  
  # set to 0 any negative value
  
  tmp<-merge(minus,plus,by="V2",sort=F)
  tmp$recombRate_cMperMb_10kb<-(tmp$V3.y-tmp$V3.x)/0.02
  
  tmp<-tmp[,c(1,4)]
  
  if(chr==1) {
    rb<-tmp
  } else {
    rb<-rbind(rb,tmp)
  }
  
  rm(minus,plus,tmp)
  
}

colnames(rb)<-c("variant","recombRate_cMperMb_10kb")
write.table(rb,paste(path,"post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/allchr_recombination_rate_10kb_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed",sep=""),
            col.names=T,row.names=F,quote = F,sep="\t")

q("no")

# # sort by original file:
awk -F'\t' '{OFS=FS} NR==FNR {h[$1] = $0; next} {print h[$4]}'  \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/allchr_recombination_rate_10kb_b38_list_union_variants_cd_uc_ibd_metaanalysis_tmp.bed \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/list_union_variants_cd_uc_ibd_metaanalysis_sorted.bed \
| awk -F'\t' -v OFS="\t" 'BEGIN {print "variant","recombRate_cMperMb_10kb"}1' > ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/allchr_recombination_rate_10kb_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed

gzip ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/allchr_recombination_rate_10kb_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed

# mv data to final folder
mv ${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/oxford_recombination_rate/b38/allchr_recombination_rate_10kb_b38_list_union_variants_cd_uc_ibd_metaanalysis.bed.gz \
${path_gwas}post_imputation/2022/analysis/metaanalysis/annotation/matrix/individual_files/files_version_june2024/


