# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# DeCode data:

# OLD: path_gwas=/path/to/ibdgwas/IIBDGC/
# NEW: 
path_gwas=/path/to/ibdgwas/IIBDGC/

# /path/to/ibdgwas/IIBDGC/post_imputation/analysis/stage_2_summary_statistics/decode/



MEM=8000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1

library(data.table)
path<-"/path/to/ibdgwas/IIBDGC/"

pheno<-c("CrohnsDisease","InflammatoryBowelDisease","UlcerativeColitis")

for (i in 1:length(pheno)) {
  print(pheno[i])
  dat<-fread(paste(path,"post_imputation/analysis/stage_2_summary_statistics/decode/DECODE_",pheno[i],"_30062021.txt.gz",sep=""),head=T)
  print(nrow(dat))
  print(table(dat$CHROM))
  print(summary(dat$A1FREQ))
  print(summary(dat$INFO))
  x<-names(table(dat$ALLELE0))
  print(x[!grepl("[A-Z]",x)])
  x<-names(table(dat$ALLELE1))
  print(x[!grepl("[A-Z]",x)])
  rm(dat,x)
}

# [1] "CrohnsDisease"
# |--------------------------------------------------|
#   |==================================================|
#   [1] 24505450
# 
# chr1   chr10   chr11   chr12   chr13   chr14   chr15   chr16   chr17   chr18 
# 1945405 1208453 1242374 1153715  778904  725295  693991  783952  759217  686137 
# chr19    chr2   chr20   chr21   chr22    chr3    chr4    chr5    chr6    chr7 
# 564935 1963555  590729  303507  341320 1576143 1577225 1448452 1441346 1445923 
# chr8    chr9    chrX 
# 1232811 1173660  868401 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0010  0.0044  0.0333  0.1126  0.1886  0.5000 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8000  0.9930  0.9990  0.9782  1.0000  1.0000 
# [1] "!"
# [1] "!"
# [1] "InflammatoryBowelDisease"
# |--------------------------------------------------|
#   |==================================================|
#   [1] 24505450
# 
# chr1   chr10   chr11   chr12   chr13   chr14   chr15   chr16   chr17   chr18 
# 1945405 1208453 1242374 1153715  778904  725295  693991  783952  759217  686137 
# chr19    chr2   chr20   chr21   chr22    chr3    chr4    chr5    chr6    chr7 
# 564935 1963555  590729  303507  341320 1576143 1577225 1448452 1441346 1445923 
# chr8    chr9    chrX 
# 1232811 1173660  868401 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0010  0.0044  0.0333  0.1126  0.1886  0.5000 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8000  0.9930  0.9990  0.9782  1.0000  1.0000 
# [1] "!"
# [1] "!"
# [1] "UlcerativeColitis"
# |--------------------------------------------------|
# |==================================================|
#   [1] 24505450
# 
# chr1   chr10   chr11   chr12   chr13   chr14   chr15   chr16   chr17   chr18 
# 1945405 1208453 1242374 1153715  778904  725295  693991  783952  759217  686137 
# chr19    chr2   chr20   chr21   chr22    chr3    chr4    chr5    chr6    chr7 
# 564935 1963555  590729  303507  341320 1576143 1577225 1448452 1441346 1445923 
# chr8    chr9    chrX 
# 1232811 1173660  868401 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0010  0.0044  0.0333  0.1126  0.1886  0.5000 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8000  0.9930  0.9990  0.9782  1.0000  1.0000 


########################################################################################################################################################################
# how many of these variants are included in our current analysis:

# The ! symbol in allele annotation denotes other allele variations at this position combined when there are more than two alleles.
# It could be structural variants with complicated allele variation, or it could also be a triallelic SNP.
# 
# We only run association for variants with imputation info > 0.8 and as our association pipeline is set up to deal with very large
# number of association runs, it is not easy to change that threshold.

# The columns are wrongly labelled, A1 should be A0 and vice versa. A1 is the effect allele 
# This is true for all the markers, so you can simply exchange the columns for each other 
# A1 freq refers to the effect allele frequency, so when you have switched
# the columns it fits.


cd /path/to/ibdgwas/IIBDGC/post_imputation/analysis/stage_2_summary_statistics/decode/

# rename header to correct the issue (note:

# Ref allele in DeCode is Allele1, swap alleles and add header
zcat DECODE_CrohnsDisease_30062021.txt.gz | awk 'BEGIN {FS = " ";OFS = " "} {print $1,$2,$3,$5,$4,$6,$7,$8,$9,$10,$11,$12}' | tail -n +2 \
| cat <(zcat DECODE_CrohnsDisease_30062021.txt.gz | head -1) - | gzip > DECODE_cd_30062021_edited.txt.gz

zcat DECODE_InflammatoryBowelDisease_30062021.txt.gz | awk 'BEGIN {FS = " ";OFS = " "} {print $1,$2,$3,$5,$4,$6,$7,$8,$9,$10,$11,$12}' | tail -n +2 \
| cat <(zcat DECODE_InflammatoryBowelDisease_30062021.txt.gz | head -1) - | gzip > DECODE_ibd_30062021_edited.txt.gz

zcat DECODE_UlcerativeColitis_30062021.txt.gz | awk 'BEGIN {FS = " ";OFS = " "} {print $1,$2,$3,$5,$4,$6,$7,$8,$9,$10,$11,$12}' | tail -n +2 \
| cat <(zcat DECODE_UlcerativeColitis_30062021.txt.gz | head -1) - | gzip > DECODE_uc_30062021_edited.txt.gz


# split files in chr:

path_gwas=/path/to/ibdgwas/IIBDGC/
pheno=(ibd cd uc)
MEM=1000

for ph in ${pheno[@]}
do
for chr in X
do
bsub -J"decode_split" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/log/${ph}_${chr}_decode_split_stderr \
-o ${path_gwas}post_imputation/log/${ph}_${chr}_decode_split_stdout \
"zcat ${path_gwas}post_imputation/analysis/stage_2_summary_statistics/decode/edited_files/DECODE_${ph}_30062021_edited.txt.gz \
| grep 'chr${chr} ' | cat <(zcat ${path_gwas}post_imputation/analysis/stage_2_summary_statistics/decode/edited_files/DECODE_${ph}_30062021_edited.txt.gz | head -1) -  \
| gzip > ${path_gwas}post_imputation/analysis/stage_2_summary_statistics/decode/edited_files/chr${chr}_DECODE_${ph}_30062021_edited.txt.gz"
done
done



MEM=8000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1

library(data.table)
path<-"/path/to/ibdgwas/IIBDGC/"


# pheno<-c("cd","ibd","uc")
pheno<-c("ibd")
# list of variants used in our analysis:

# THIS SCRIPT ASSUMES ALL VARIANTS REFER TO + STRAND:

mat<-as.data.frame(matrix(ncol=6,nrow=22))
colnames(mat)<-c("chr","iibdgc","decode","shared","no_specific_allele_decode","no_specific_allele_decode_shared_position_and_one_allele")
mat$chr<-seq(1:22)

for (i in 1) {

  print(pheno[i])
  
  for (chr in 1:22) {
    
    print(chr)
    
    dat<-fread(paste(path,"post_imputation/analysis/stage_2_summary_statistics/decode/chr",chr,"_DECODE_",pheno[i],"_30062021_edited.txt.gz",sep=""),head=T)
    dat$X<-paste(dat$CHROM,dat$GENPOS_build38,dat$ALLELE0,dat$ALLELE1,sep=":")
    dat$Y<-paste(dat$CHROM,dat$GENPOS_build38,dat$ALLELE1,dat$ALLELE0,sep=":")
    
    tmp<-fread(paste(path,"post_imputation/analysis/metaanalysis/",pheno[i],"/",chr,"_",pheno[i],"_meta_noGC_PCs_firthse_no_illumina550_1_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""),head=T)
    # tmp<-tmp[which(tmp$rate_total_sample>=0.8),]
    
    tmp$chr<-paste("chr",chr,sep="")
    tmp$position<-gsub("chr[1-9]{1,2}:","",tmp$MarkerName)
    tmp$position<-gsub(":.*:.*$","",tmp$position)

    mat$iibdgc[which(mat$chr==chr)]<-nrow(tmp)
    mat$decode[which(mat$chr==chr)]<-nrow(dat)
    
    # N variants with same reference and alternative alleles:
    dat.1<-dat[which(dat$X %in% tmp$MarkerName),]

    # N variants with swapped ref and alt
    dat.2<-dat[which(dat$Y %in% tmp$MarkerName),]
    dat.2$BETA.Y1_ed<-(dat.2$BETA.Y1)*-1
    dat.2$A1FREQ_ed<-1-dat.2$A1FREQ
    dat.2<-dat.2[,c("CHROM","GENPOS_build38","ID","ALLELE1","ALLELE0","A1FREQ_ed","INFO","TEST","BETA.Y1_ed","SE.Y1","CHISQ.Y1","LOG10P.Y1")]
    colnames(dat.2)<-colnames(dat.1)[1:12]
    
    
    # N variants with one allele ! and the other available in IIBDGC data
    vec1<-c(paste(tmp$chr,tmp$position,tmp$A1,sep=":"),paste(tmp$chr,tmp$position,tmp$A2,sep=":"))
    
    dat.3<-dat[which( !(dat$Y %in% tmp$MarkerName) & !(dat$X %in% tmp$MarkerName)),]
    mat$no_specific_allele_decode<-nrow(dat.3)
    
    dat.3$xx<-gsub(":!","",dat.3$X)
    
    mat$no_specific_allele_decode_shared_position_and_one_allele[which(mat$chr==chr)]<-nrow(dat.3[which(dat.3$xx %in% vec1),])
    
    # some in tmp both directions,
    # dat.4<-dat[which( (dat$Y %in% tmp$MarkerName) & (dat$X %in% tmp$MarkerName)),]
    
    #### COMPARE EFFECT ALLELE FREQUENCY:
    
    dat<-rbind(dat.1[,1:12],dat.2)
    mat$shared[which(mat$chr==chr)]<-nrow(dat)
    
    dat$X<-paste(dat$CHROM,dat$GENPOS_build38,dat$ALLELE0,dat$ALLELE1,sep=":")
    dat<-merge(tmp,dat,by.x="MarkerName",by.y="X")
    
    if(chr==1){
      dat_final<-dat
      # dat.3_final<-dat.3
    } else {
      dat_final<-rbind(dat_final,dat)
      # dat.3_final<-rbind(dat.3_final,dat.3)
    }
    rm(tmp,dat,dat.1,dat.2,dat.3)
  }
}

write.table(mat,paste(path,"post_imputation/analysis/stage_2_summary_statistics/decode/number_variants_matching_iibdgc_meta_analysis.tsv",sep=""),
            # col.names=T,row.names=F,quote=F,sep="\t")

# write.table(mat,paste(path,"post_imputation/analysis/stage_2_summary_statistics/decode/number_variants_matching_iibdgc_meta_analysis_sample0.8.tsv",sep=""),
#             col.names=T,row.names=F,quote=F,sep="\t")

# plot effect allele frequencies to confirm match (no strand issues?)
dim(dat_final)
# [1] 11436286       30

sum(mat$shared)
# [1] 11436286

dim(dat_final[which(dat_final$rate_total_sample>=0.8),])
# [1] 9085628      30


sum(mat$shared)/sum(mat$decode)
# [1] 0.4838288

sum(mat$no_specific_allele_decode_shared_position)/sum(mat$decode)
# [1] 0.06262613

sum(mat$no_specific_allele_decode)/sum(mat$decode)
# [1] 0.316879

sum(mat$iibdgc)
# [1] 20895613

sum(mat$no_specific_allele_decode_shared_position)/sum(mat$iibdgc)
# [1] 0.07084248

dat_final$value<-((dat_final$avgA2FREQ-dat_final$A1FREQ)^2)/
  ((dat_final$avgA2FREQ+dat_final$A1FREQ)*
     (2-dat_final$avgA2FREQ-dat_final$A1FREQ))

dat_final$inconsistency_maf<-NA
dat_final$inconsistency_maf[which(dat_final$value>=0.125)]<-"1"
dat_final$inconsistency_maf[which(dat_final$value<0.125)]<-"0"

table(dat_final$inconsistency_maf)
# 0        1 
# 11435161     1125 

# plot those with iibdgc sample rate>0.8
library("ggpubr")

pdf("~/tmp_plots/effect_allele_decode_vs_iibdgc_metaanalysis.pdf")
ggscatter(dat_final[which(dat_final$rate_total_sample>=0.8),], x = "avgA2FREQ", y = "A1FREQ", fill="inconsistency_maf",
          shape = 21,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "IIBDGC", ylab = "DeCODE")
dev.off()


#### get summary table:

  


###########################################################################################################################################
# create files for MR-MEGA (see https://genomics.ut.ee/en/mr-mega) AND METAL
  
# Each GWA study file has mandatory column headers:
# 1) MARKERNAME – snp name
# 2) EA – effect allele
# 3) NEA – non effect allele
# 4) OR - odds ratio
# 5) OR_95L - lower confidence interval of OR
# 6) OR_95U - upper confidence interval of OR
# 7) EAF – effect allele frequency
# 8) N - sample size
# 9) CHROMOSOME  - chromosome of marker
# 10) POSITION - position of marker

MEM=8000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1

library(data.table)

pheno<-c("cd","ibd","uc")
path<-"/path/to/ibdgwas/IIBDGC/"

for (i in 1:length(pheno)) {
  
  # for (chr in 1:22) X {
    
  for (chr in "X") {
    
    print(chr)
    dat<-fread(paste(path,"post_imputation/analysis/stage_2_summary_statistics/decode/edited_files/chr",chr,"_DECODE_",pheno[i],"_30062021_edited.txt.gz",sep=""),head=T)
    
    # exclude any variant with no alleles
    dat<-dat[which(dat$ALLELE0!="!" & dat$ALLELE1!="!"),]
    dat$X<-paste(dat$CHROM,dat$GENPOS_build38,dat$ALLELE0,dat$ALLELE1,sep=":")
    dat$Y<-paste(dat$CHROM,dat$GENPOS_build38,dat$ALLELE1,dat$ALLELE0,sep=":")
    
    
    tmp<-fread(paste(path,"post_imputation/analysis/metaanalysis/",pheno[i],"/",chr,"_",pheno[i],"_meta_noGC_PCs_firthse_no_illumina550_1_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""),head=T)
    # tmp<-tmp[which(tmp$rate_total_sample>=0.8),]
    
    tmp$chr<-paste("chr",chr,sep="")
    tmp$position<-gsub("chr[1-9]{1,2}:","",tmp$MarkerName)
    tmp$position<-gsub(":.*:.*$","",tmp$position)
    
    # N variants with same reference and alternative alleles:
    dat.1<-dat[which(dat$X %in% tmp$MarkerName),]
    
    # N variants with swapped ref and alt
    dat.2<-dat[which( (dat$Y %in% tmp$MarkerName) & !(dat$X %in% tmp$MarkerName)),]
    dat.2$BETA.Y1_ed<-(dat.2$BETA.Y1)*-1
    dat.2$A1FREQ_ed<-1-dat.2$A1FREQ
    dat.2<-dat.2[,c("CHROM","GENPOS_build38","ID","ALLELE1","ALLELE0","A1FREQ_ed","INFO","TEST","BETA.Y1_ed","SE.Y1","CHISQ.Y1","LOG10P.Y1")]
    colnames(dat.2)<-colnames(dat.1)[1:12]
    
    # other
    dat.3<-dat[which(!(dat$Y %in% tmp$MarkerName) & !(dat$X %in% tmp$MarkerName)),]
    
    dat<-rbind(dat.1[,1:12],dat.2,dat.3[,1:12])
    rm(dat.1,dat.2,dat.3,tmp)
    
    dat<-dat[order(dat$GENPOS_build38,decreasing=F),]
    dat$ID_ed<-paste(dat$CHROM,dat$GENPOS_build38,dat$ALLELE0,dat$ALLELE1,sep=":")

    dat<-dat[,c("CHROM","GENPOS_build38","ID_ed","ALLELE0","ALLELE1","A1FREQ","INFO","TEST","BETA.Y1","SE.Y1","CHISQ.Y1","LOG10P.Y1")]
    colnames(dat)<-c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","TEST","BETA.Y1","SE.Y1","CHISQ.Y1","LOG10P.Y1")
    
    fwrite(dat,paste(path,"post_imputation/analysis/stage_2_summary_statistics/decode/",pheno[i],"/chr",chr,"_DECODE_",pheno[i],"_30062021_edited_noMult.regenie",sep=""),
           col.names=T,row.names=F,sep="\t",quote=F)
    
  }
  
}

###########################################################################################################################################
###########################################################################################################################################

# mv data to specific folders:

cd /path/to/ibdgwas/IIBDGC/post_imputation/analysis/stage_2_summary_statistics/decode/
pheno=(ibd cd uc)

for ph in ${pheno[@]}
do mkdir ${ph}
mv chr*_DECODE_${ph}_30062021_edited_noMult.regenie ${ph}/
mv chr*_DECODE_${ph}_30062021_MRMEGA_format.txt.gz ${ph}/
done

# MOVE INTERMEDIATE FILES AS WELL, REMOVE AFTER CHECKING ALL IS OK
mkdir edited_files/
mv *edited* edited_files/


  
###########################################################################################################################################
###########################################################################################################################################

# updated to match new format - 2022/2023 analysis:
  
MEM=8000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group /software/bin/R-4.3.1


pheno<-c("cd","ibd","uc")
path<-"/path/to/ibdgwas/IIBDGC/"

library(data.table)

for (i in 1:length(pheno)) {
  
  print(pheno[i])
  
  # for (chr in 1:22) {
  
  for (chr in "X") {
    
    print(chr)
    dat<-fread(paste(path,"post_imputation/analysis/stage_2_summary_statistics/decode/edited_files/chr",chr,"_DECODE_",pheno[i],"_30062021_edited.txt.gz",sep=""),head=T)
    
    # exclude any variant with no alleles
    dat<-dat[which(dat$ALLELE0!="!" & dat$ALLELE1!="!"),]
    dat$X<-paste(dat$CHROM,dat$GENPOS_build38,dat$ALLELE0,dat$ALLELE1,sep=":")
    dat$Y<-paste(dat$CHROM,dat$GENPOS_build38,dat$ALLELE1,dat$ALLELE0,sep=":")
    
    dat$XX<-paste(dat$CHROM,dat$GENPOS_build38,dat$ALLELE0,sep=":")
    dat$YY<-paste(dat$CHROM,dat$GENPOS_build38,dat$ALLELE1,sep=":")
    
    dat$Z<-paste(dat$CHROM,dat$GENPOS_build38,sep=":")
    
    tmp<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/chr",chr,"_",pheno[i],"_meta_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_with_rsid_dbsnp154.txt.gz",sep=""),head=T)
    # tmp<-tmp[which(tmp$rate_total_sample>=0.8),]
    
    tmp$chr<-paste("chr",chr,sep="")
    tmp$XX<-paste(tmp$chr,tmp$Position_b38,tmp$A1,sep=":")
    tmp$Z<-paste(tmp$chr,tmp$Position_b38,sep=":")
      
    
    # N variants with same reference and alternative alleles:
    dat.1<-dat[which(dat$X %in% tmp$MarkerName),]
    print(paste("N variants same ref and alt:",nrow(dat.1)))
    
    # N variants with swapped ref and alt
    dat.2<-dat[which( (dat$Y %in% tmp$MarkerName) & !(dat$X %in% tmp$MarkerName)),]
    dat.2$BETA.Y1_ed<-(dat.2$BETA.Y1)*-1
    dat.2$A1FREQ_ed<-1-dat.2$A1FREQ
    dat.2<-dat.2[,c("CHROM","GENPOS_build38","ID","ALLELE1","ALLELE0","A1FREQ_ed","INFO","TEST","BETA.Y1_ed","SE.Y1","CHISQ.Y1","LOG10P.Y1")]
    colnames(dat.2)<-colnames(dat.1)[1:12]
    print(paste("N variants swapped ref and alt:",nrow(dat.2)))
    
    
    # same reference (but not alt)
    dat.3<-dat[which( (dat$XX %in% tmp$XX) & !(dat$Y %in% tmp$MarkerName) & !(dat$X %in% tmp$MarkerName)),]
    print(paste("N variants same reference but swapped alt:",nrow(dat.3)))
    
    # swapped reference (but not alt)
    dat.4<-dat[which( (dat$YY %in% tmp$XX) & !(dat$XX %in% tmp$XX) & !(dat$Y %in% tmp$MarkerName) & !(dat$X %in% tmp$MarkerName)),]
    dat.4$BETA.Y1_ed<-(dat.4$BETA.Y1)*-1
    dat.4$A1FREQ_ed<-1-dat.4$A1FREQ
    dat.4<-dat.4[,c("CHROM","GENPOS_build38","ID","ALLELE1","ALLELE0","A1FREQ_ed","INFO","TEST","BETA.Y1_ed","SE.Y1","CHISQ.Y1","LOG10P.Y1")]
    colnames(dat.4)<-colnames(dat.1)[1:12]
    print(paste("N variants swapped reference different alt:",nrow(dat.4)))
    
    
    # no matching ref or alt, but same position - indels with +1nt in ref
    dat.5<-dat[which( (dat$Z %in% tmp$Z) & !(dat$YY %in% tmp$XX) & !(dat$XX %in% tmp$XX) & !(dat$Y %in% tmp$MarkerName) & !(dat$X %in% tmp$MarkerName)),]
    print(paste("N variants same position:",nrow(dat.5)))
    
    
    # rest, variants not in iibdgc meta:
    dat.6<-dat[which( !(dat$Z %in% tmp$Z) & !(dat$YY %in% tmp$XX) & !(dat$XX %in% tmp$XX) & !(dat$Y %in% tmp$MarkerName) & !(dat$X %in% tmp$MarkerName)),]
    print(paste("N variants not in iibdgc meta:",nrow(dat.6)))
    
    
    
    print(table(sum(nrow(dat.1),nrow(dat.2),nrow(dat.3),nrow(dat.4),nrow(dat.5),nrow(dat.6))==nrow(dat)))
    
    
    dat<-rbind(dat.1[,1:12],dat.2[,1:12],dat.3[,1:12],dat.4[,1:12],dat.5[,1:12],dat.6[,1:12])
    dat$ID_ed<-paste(dat$CHROM,dat$GENPOS_build38,dat$ALLELE0,dat$ALLELE1,sep=":")
    
    dat$A1FREQ_CASES<-NA
    dat$A1FREQ_CONTROLS<-NA
    dat$INFO<-1
    dat$EXTRA<-NA
    
    if(pheno[i]=="ibd"){
      dat$N<-sum(6069,357620)
      dat$N_CASES<-6069
      dat$N_CONTROLS<-357620 
    }

    if(pheno[i]=="cd"){
      dat$N<-sum(870,348243)
      dat$N_CASES<-870
      dat$N_CONTROLS<-348243 
    }
        
    if(pheno[i]=="uc"){
      dat$N<-sum(2548,369816)
      dat$N_CASES<-2548
      dat$N_CONTROLS<-369816 
    }
      
    dat<-dat[,c("CHROM","GENPOS_build38","ID_ed","ALLELE0","ALLELE1","A1FREQ","A1FREQ_CASES","A1FREQ_CONTROLS","INFO",
                "N","N_CASES","N_CONTROLS","TEST",
                "BETA.Y1","SE.Y1","CHISQ.Y1","LOG10P.Y1","EXTRA")]
    
    colnames(dat)<-c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","A1FREQ_CASES","A1FREQ_CONTROLS","INFO",
                     "N","N_CASES","N_CONTROLS","TEST",
                     "BETA","SE","CHISQ","LOG10P","EXTRA")
    dat<-dat[order(dat$GENPOS,decreasing=F),]
    

    fwrite(dat,paste(path,"post_imputation/analysis/stage_2_summary_statistics/decode/",pheno[i],"/chr",chr,"_DECODE_",pheno[i],"_30062021_edited_noMult.regenie",sep=""),
           col.names=T,row.names=F,sep="\t",quote=F,na="NA")
    
    rm(dat.1,dat.2,dat.3,dat.4,dat.5,dat.6,tmp,dat)
    
  }
  
}

#########################
#
# [1] "cd"
# [1] 1
# [1] "N variants same ref and alt: 778024"
# [1] "N variants swapped ref and alt: 121867"
# [1] "N variants same reference but swapped alt: 487"
# [1] "N variants swapped reference different alt: 85"
# [1] "N variants same position: 397"
# [1] "N variants not in iibdgc meta: 419764"
# 
# TRUE 
# 1 
# [1] 2
# [1] "N variants same ref and alt: 850527"
# [1] "N variants swapped ref and alt: 127298"
# [1] "N variants same reference but swapped alt: 564"
# [1] "N variants swapped reference different alt: 95"
# [1] "N variants same position: 363"
# [1] "N variants not in iibdgc meta: 363277"
# 
# TRUE 
# 1 
# [1] 3
# [1] "N variants same ref and alt: 713313"
# [1] "N variants swapped ref and alt: 110583"
# [1] "N variants same reference but swapped alt: 458"
# [1] "N variants swapped reference different alt: 76"
# [1] "N variants same position: 272"
# [1] "N variants not in iibdgc meta: 252757"
# 
# TRUE 
# 1 
# [1] 4
# |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 721511"
# [1] "N variants swapped ref and alt: 116731"
# [1] "N variants same reference but swapped alt: 446"
# [1] "N variants swapped reference different alt: 72"
# [1] "N variants same position: 358"
# [1] "N variants not in iibdgc meta: 250606"
# 
# TRUE 
# 1 
# [1] 5
# [1] "N variants same ref and alt: 654353"
# [1] "N variants swapped ref and alt: 97635"
# [1] "N variants same reference but swapped alt: 388"
# [1] "N variants swapped reference different alt: 60"
# [1] "N variants same position: 283"
# [1] "N variants not in iibdgc meta: 241181"
# 
# TRUE 
# 1 
# [1] 6
# [1] "N variants same ref and alt: 657120"
# [1] "N variants swapped ref and alt: 99632"
# [1] "N variants same reference but swapped alt: 409"
# [1] "N variants swapped reference different alt: 64"
# [1] "N variants same position: 285"
# [1] "N variants not in iibdgc meta: 225464"
# 
# TRUE 
# 1 
# [1] 7
# [1] "N variants same ref and alt: 590621"
# [1] "N variants swapped ref and alt: 88385"
# [1] "N variants same reference but swapped alt: 401"
# [1] "N variants swapped reference different alt: 60"
# [1] "N variants same position: 244"
# [1] "N variants not in iibdgc meta: 309777"
# 
# TRUE 
# 1 
# [1] 8
# |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 552784"
# [1] "N variants swapped ref and alt: 83237"
# [1] "N variants same reference but swapped alt: 322"
# [1] "N variants swapped reference different alt: 40"
# [1] "N variants same position: 240"
# [1] "N variants not in iibdgc meta: 216750"
# 
# TRUE 
# 1 
# [1] 9
# [1] "N variants same ref and alt: 437584"
# [1] "N variants swapped ref and alt: 64881"
# [1] "N variants same reference but swapped alt: 311"
# [1] "N variants swapped reference different alt: 35"
# [1] "N variants same position: 198"
# [1] "N variants not in iibdgc meta: 317109"
# 
# TRUE 
# 1 
# [1] 10
# [1] "N variants same ref and alt: 506730"
# [1] "N variants swapped ref and alt: 79077"
# [1] "N variants same reference but swapped alt: 275"
# [1] "N variants swapped reference different alt: 36"
# [1] "N variants same position: 198"
# [1] "N variants not in iibdgc meta: 242631"
# 
# TRUE 
# 1 
# [1] 11
# [1] "N variants same ref and alt: 489417"
# [1] "N variants swapped ref and alt: 80254"
# [1] "N variants same reference but swapped alt: 311"
# [1] "N variants swapped reference different alt: 47"
# [1] "N variants same position: 225"
# [1] "N variants not in iibdgc meta: 310254"
# 
# TRUE 
# 1 
# [1] 12
# [1] "N variants same ref and alt: 477756"
# [1] "N variants swapped ref and alt: 73548"
# [1] "N variants same reference but swapped alt: 314"
# [1] "N variants swapped reference different alt: 53"
# [1] "N variants same position: 214"
# [1] "N variants not in iibdgc meta: 229303"
# 
# TRUE 
# 1 
# [1] 13
# [1] "N variants same ref and alt: 362340"
# [1] "N variants swapped ref and alt: 60991"
# [1] "N variants same reference but swapped alt: 208"
# [1] "N variants swapped reference different alt: 32"
# [1] "N variants same position: 194"
# [1] "N variants not in iibdgc meta: 106171"
# 
# TRUE 
# 1 
# [1] 14
# [1] "N variants same ref and alt: 327524"
# [1] "N variants swapped ref and alt: 49282"
# [1] "N variants same reference but swapped alt: 200"
# [1] "N variants swapped reference different alt: 23"
# [1] "N variants same position: 160"
# [1] "N variants not in iibdgc meta: 111557"
# 
# TRUE 
# 1 
# [1] 15
# [1] "N variants same ref and alt: 284501"
# [1] "N variants swapped ref and alt: 43696"
# [1] "N variants same reference but swapped alt: 166"
# [1] "N variants swapped reference different alt: 25"
# [1] "N variants same position: 155"
# [1] "N variants not in iibdgc meta: 141160"
# 
# TRUE 
# 1 
# [1] 16
# [1] "N variants same ref and alt: 308210"
# [1] "N variants swapped ref and alt: 45226"
# [1] "N variants same reference but swapped alt: 201"
# [1] "N variants swapped reference different alt: 22"
# [1] "N variants same position: 101"
# [1] "N variants not in iibdgc meta: 169755"
# 
# TRUE 
# 1 
# [1] 17
# [1] "N variants same ref and alt: 270786"
# [1] "N variants swapped ref and alt: 39718"
# [1] "N variants same reference but swapped alt: 168"
# [1] "N variants swapped reference different alt: 19"
# [1] "N variants same position: 107"
# [1] "N variants not in iibdgc meta: 194082"
# 
# TRUE 
# 1 
# [1] 18
# [1] "N variants same ref and alt: 286650"
# [1] "N variants swapped ref and alt: 44109"
# [1] "N variants same reference but swapped alt: 189"
# [1] "N variants swapped reference different alt: 27"
# [1] "N variants same position: 148"
# [1] "N variants not in iibdgc meta: 145821"
# 
# TRUE 
# 1 
# [1] 19
# [1] "N variants same ref and alt: 223425"
# [1] "N variants swapped ref and alt: 31530"
# [1] "N variants same reference but swapped alt: 159"
# [1] "N variants swapped reference different alt: 25"
# [1] "N variants same position: 98"
# [1] "N variants not in iibdgc meta: 95050"
# 
# TRUE 
# 1 
# [1] 20
# [1] "N variants same ref and alt: 225069"
# [1] "N variants swapped ref and alt: 32563"
# [1] "N variants same reference but swapped alt: 123"
# [1] "N variants swapped reference different alt: 13"
# [1] "N variants same position: 96"
# [1] "N variants not in iibdgc meta: 151692"
# 
# TRUE 
# 1 
# [1] 21
# [1] "N variants same ref and alt: 132947"
# [1] "N variants swapped ref and alt: 22543"
# [1] "N variants same reference but swapped alt: 81"
# [1] "N variants swapped reference different alt: 17"
# [1] "N variants same position: 65"
# [1] "N variants not in iibdgc meta: 50721"
# 
# TRUE 
# 1 
# [1] 22
# [1] "N variants same ref and alt: 139022"
# [1] "N variants swapped ref and alt: 20073"
# [1] "N variants same reference but swapped alt: 102"
# [1] "N variants swapped reference different alt: 15"
# [1] "N variants same position: 69"
# [1] "N variants not in iibdgc meta: 67315"
# 
# TRUE 
# 1 
# [1] "ibd"
# [1] 1
# |--------------------------------------------------|
#   |==================================================|
#   |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 788429"
# [1] "N variants swapped ref and alt: 121909"
# [1] "N variants same reference but swapped alt: 494"
# [1] "N variants swapped reference different alt: 84"
# [1] "N variants same position: 391"
# [1] "N variants not in iibdgc meta: 409317"
# 
# TRUE 
# 1 
# [1] 2
# |--------------------------------------------------|
#   |==================================================|
#   |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 860731"
# [1] "N variants swapped ref and alt: 127314"
# [1] "N variants same reference but swapped alt: 572"
# [1] "N variants swapped reference different alt: 93"
# [1] "N variants same position: 374"
# [1] "N variants not in iibdgc meta: 353040"
# 
# TRUE 
# 1 
# [1] 3
# |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 722270"
# [1] "N variants swapped ref and alt: 110604"
# [1] "N variants same reference but swapped alt: 459"
# [1] "N variants swapped reference different alt: 78"
# [1] "N variants same position: 274"
# [1] "N variants not in iibdgc meta: 243774"
# 
# TRUE 
# 1 
# [1] 4
# |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 730386"
# [1] "N variants swapped ref and alt: 116739"
# [1] "N variants same reference but swapped alt: 455"
# [1] "N variants swapped reference different alt: 72"
# [1] "N variants same position: 357"
# [1] "N variants not in iibdgc meta: 241715"
# 
# TRUE 
# 1 
# [1] 5
# |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 662032"
# [1] "N variants swapped ref and alt: 97674"
# [1] "N variants same reference but swapped alt: 388"
# [1] "N variants swapped reference different alt: 62"
# [1] "N variants same position: 280"
# [1] "N variants not in iibdgc meta: 233464"
# 
# TRUE 
# 1 
# [1] 6
# |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 664838"
# [1] "N variants swapped ref and alt: 99640"
# [1] "N variants same reference but swapped alt: 410"
# [1] "N variants swapped reference different alt: 61"
# [1] "N variants same position: 287"
# [1] "N variants not in iibdgc meta: 217738"
# 
# TRUE 
# 1 
# [1] 7
# |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 597934"
# [1] "N variants swapped ref and alt: 88388"
# [1] "N variants same reference but swapped alt: 398"
# [1] "N variants swapped reference different alt: 56"
# [1] "N variants same position: 240"
# [1] "N variants not in iibdgc meta: 302472"
# 
# TRUE 
# 1 
# [1] 8
# [1] "N variants same ref and alt: 561395"
# [1] "N variants swapped ref and alt: 83249"
# [1] "N variants same reference but swapped alt: 349"
# [1] "N variants swapped reference different alt: 40"
# [1] "N variants same position: 244"
# [1] "N variants not in iibdgc meta: 208096"
# 
# TRUE 
# 1 
# [1] 9
# |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 444549"
# [1] "N variants swapped ref and alt: 64894"
# [1] "N variants same reference but swapped alt: 316"
# [1] "N variants swapped reference different alt: 36"
# [1] "N variants same position: 202"
# [1] "N variants not in iibdgc meta: 310121"
# 
# TRUE 
# 1 
# [1] 10
# |--------------------------------------------------|
#   |==================================================|
#   |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 513111"
# [1] "N variants swapped ref and alt: 79054"
# [1] "N variants same reference but swapped alt: 280"
# [1] "N variants swapped reference different alt: 32"
# [1] "N variants same position: 201"
# [1] "N variants not in iibdgc meta: 236269"
# 
# TRUE 
# 1 
# [1] 11
# |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 495213"
# [1] "N variants swapped ref and alt: 80266"
# [1] "N variants same reference but swapped alt: 317"
# [1] "N variants swapped reference different alt: 50"
# [1] "N variants same position: 229"
# [1] "N variants not in iibdgc meta: 304433"
# 
# TRUE 
# 1 
# [1] 12
# |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 484260"
# [1] "N variants swapped ref and alt: 73550"
# [1] "N variants same reference but swapped alt: 285"
# [1] "N variants swapped reference different alt: 54"
# [1] "N variants same position: 217"
# [1] "N variants not in iibdgc meta: 222822"
# 
# TRUE 
# 1 
# [1] 13
# |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 366802"
# [1] "N variants swapped ref and alt: 61060"
# [1] "N variants same reference but swapped alt: 204"
# [1] "N variants swapped reference different alt: 34"
# [1] "N variants same position: 195"
# [1] "N variants not in iibdgc meta: 101641"
# 
# TRUE 
# 1 
# [1] 14
# |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 331160"
# [1] "N variants swapped ref and alt: 49258"
# [1] "N variants same reference but swapped alt: 191"
# [1] "N variants swapped reference different alt: 24"
# [1] "N variants same position: 163"
# [1] "N variants not in iibdgc meta: 107950"
# 
# TRUE 
# 1 
# [1] 15
# |--------------------------------------------------|
#   |==================================================|
#   |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 287794"
# [1] "N variants swapped ref and alt: 43708"
# [1] "N variants same reference but swapped alt: 163"
# [1] "N variants swapped reference different alt: 25"
# [1] "N variants same position: 153"
# [1] "N variants not in iibdgc meta: 137860"
# 
# TRUE 
# 1 
# [1] 16
# |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 312170"
# [1] "N variants swapped ref and alt: 45226"
# [1] "N variants same reference but swapped alt: 184"
# [1] "N variants swapped reference different alt: 21"
# [1] "N variants same position: 100"
# [1] "N variants not in iibdgc meta: 165814"
# 
# TRUE 
# 1 
# [1] 17
# |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 274066"
# [1] "N variants swapped ref and alt: 39731"
# [1] "N variants same reference but swapped alt: 166"
# [1] "N variants swapped reference different alt: 19"
# [1] "N variants same position: 111"
# [1] "N variants not in iibdgc meta: 190787"
# 
# TRUE 
# 1 
# [1] 18
# |--------------------------------------------------|
#   |==================================================|
#   |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 289938"
# [1] "N variants swapped ref and alt: 44109"
# [1] "N variants same reference but swapped alt: 183"
# [1] "N variants swapped reference different alt: 25"
# [1] "N variants same position: 148"
# [1] "N variants not in iibdgc meta: 142541"
# 
# TRUE 
# 1 
# [1] 19
# |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 226002"
# [1] "N variants swapped ref and alt: 31535"
# [1] "N variants same reference but swapped alt: 164"
# [1] "N variants swapped reference different alt: 24"
# [1] "N variants same position: 100"
# [1] "N variants not in iibdgc meta: 92462"
# 
# TRUE 
# 1 
# [1] 20
# |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 227429"
# [1] "N variants swapped ref and alt: 32567"
# [1] "N variants same reference but swapped alt: 114"
# [1] "N variants swapped reference different alt: 14"
# [1] "N variants same position: 93"
# [1] "N variants not in iibdgc meta: 149339"
# 
# TRUE 
# 1 
# [1] 21
# [1] "N variants same ref and alt: 134473"
# [1] "N variants swapped ref and alt: 22539"
# [1] "N variants same reference but swapped alt: 80"
# [1] "N variants swapped reference different alt: 16"
# [1] "N variants same position: 64"
# [1] "N variants not in iibdgc meta: 49202"
# 
# TRUE 
# 1 
# [1] 22
# |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 140628"
# [1] "N variants swapped ref and alt: 20078"
# [1] "N variants same reference but swapped alt: 104"
# [1] "N variants swapped reference different alt: 16"
# [1] "N variants same position: 70"
# [1] "N variants not in iibdgc meta: 65700"
# 
# TRUE 
# 1 
# [1] "uc"
# [1] 1
# [1] "N variants same ref and alt: 792514"
# [1] "N variants swapped ref and alt: 121936"
# [1] "N variants same reference but swapped alt: 509"
# [1] "N variants swapped reference different alt: 87"
# [1] "N variants same position: 392"
# [1] "N variants not in iibdgc meta: 405186"
# 
# TRUE 
# 1 
# [1] 2
# [1] "N variants same ref and alt: 864931"
# [1] "N variants swapped ref and alt: 127270"
# [1] "N variants same reference but swapped alt: 625"
# [1] "N variants swapped reference different alt: 101"
# [1] "N variants same position: 377"
# [1] "N variants not in iibdgc meta: 348820"
# 
# TRUE 
# 1 
# [1] 3
# |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 725684"
# [1] "N variants swapped ref and alt: 110676"
# [1] "N variants same reference but swapped alt: 475"
# [1] "N variants swapped reference different alt: 80"
# [1] "N variants same position: 284"
# [1] "N variants not in iibdgc meta: 240260"
# 
# TRUE 
# 1 
# [1] 4
# [1] "N variants same ref and alt: 733963"
# [1] "N variants swapped ref and alt: 116776"
# [1] "N variants same reference but swapped alt: 468"
# [1] "N variants swapped reference different alt: 75"
# [1] "N variants same position: 363"
# [1] "N variants not in iibdgc meta: 238079"
# 
# TRUE 
# 1 
# [1] 5
# [1] "N variants same ref and alt: 665253"
# [1] "N variants swapped ref and alt: 97680"
# [1] "N variants same reference but swapped alt: 401"
# [1] "N variants swapped reference different alt: 59"
# [1] "N variants same position: 284"
# [1] "N variants not in iibdgc meta: 230223"
# 
# TRUE 
# 1 
# [1] 6
# [1] "N variants same ref and alt: 667572"
# [1] "N variants swapped ref and alt: 99639"
# [1] "N variants same reference but swapped alt: 437"
# [1] "N variants swapped reference different alt: 58"
# [1] "N variants same position: 288"
# [1] "N variants not in iibdgc meta: 214980"
# 
# TRUE 
# 1 
# [1] 7
# [1] "N variants same ref and alt: 600952"
# [1] "N variants swapped ref and alt: 88393"
# [1] "N variants same reference but swapped alt: 430"
# [1] "N variants swapped reference different alt: 58"
# [1] "N variants same position: 244"
# [1] "N variants not in iibdgc meta: 299411"
# 
# TRUE 
# 1 
# [1] 8
# [1] "N variants same ref and alt: 561834"
# [1] "N variants swapped ref and alt: 83258"
# [1] "N variants same reference but swapped alt: 337"
# [1] "N variants swapped reference different alt: 39"
# [1] "N variants same position: 242"
# [1] "N variants not in iibdgc meta: 207663"
# 
# TRUE 
# 1 
# [1] 9
# [1] "N variants same ref and alt: 444851"
# [1] "N variants swapped ref and alt: 64906"
# [1] "N variants same reference but swapped alt: 299"
# [1] "N variants swapped reference different alt: 35"
# [1] "N variants same position: 199"
# [1] "N variants not in iibdgc meta: 309828"
# 
# TRUE 
# 1 
# [1] 10
# [1] "N variants same ref and alt: 516013"
# [1] "N variants swapped ref and alt: 79068"
# [1] "N variants same reference but swapped alt: 304"
# [1] "N variants swapped reference different alt: 34"
# [1] "N variants same position: 205"
# [1] "N variants not in iibdgc meta: 233323"
# 
# TRUE 
# 1 
# [1] 11
# [1] "N variants same ref and alt: 497693"
# [1] "N variants swapped ref and alt: 80250"
# [1] "N variants same reference but swapped alt: 353"
# [1] "N variants swapped reference different alt: 50"
# [1] "N variants same position: 232"
# [1] "N variants not in iibdgc meta: 301930"
# 
# TRUE 
# 1 
# [1] 12
# [1] "N variants same ref and alt: 487496"
# [1] "N variants swapped ref and alt: 73568"
# [1] "N variants same reference but swapped alt: 299"
# [1] "N variants swapped reference different alt: 56"
# [1] "N variants same position: 215"
# [1] "N variants not in iibdgc meta: 219554"
# 
# TRUE 
# 1 
# [1] 13
# [1] "N variants same ref and alt: 368529"
# [1] "N variants swapped ref and alt: 61056"
# [1] "N variants same reference but swapped alt: 220"
# [1] "N variants swapped reference different alt: 36"
# [1] "N variants same position: 202"
# [1] "N variants not in iibdgc meta: 99893"
# 
# TRUE 
# 1 
# [1] 14
# [1] "N variants same ref and alt: 333173"
# [1] "N variants swapped ref and alt: 49334"
# [1] "N variants same reference but swapped alt: 208"
# [1] "N variants swapped reference different alt: 26"
# [1] "N variants same position: 163"
# [1] "N variants not in iibdgc meta: 105842"
# 
# TRUE 
# 1 
# [1] 15
# [1] "N variants same ref and alt: 289744"
# [1] "N variants swapped ref and alt: 43725"
# [1] "N variants same reference but swapped alt: 185"
# [1] "N variants swapped reference different alt: 25"
# [1] "N variants same position: 159"
# [1] "N variants not in iibdgc meta: 135865"
# 
# TRUE 
# 1 
# [1] 16
# [1] "N variants same ref and alt: 314278"
# [1] "N variants swapped ref and alt: 45234"
# [1] "N variants same reference but swapped alt: 207"
# [1] "N variants swapped reference different alt: 23"
# [1] "N variants same position: 103"
# [1] "N variants not in iibdgc meta: 163670"
# 
# TRUE 
# 1 
# [1] 17
# [1] "N variants same ref and alt: 275899"
# [1] "N variants swapped ref and alt: 39736"
# [1] "N variants same reference but swapped alt: 203"
# [1] "N variants swapped reference different alt: 21"
# [1] "N variants same position: 108"
# [1] "N variants not in iibdgc meta: 188913"
# 
# TRUE 
# 1 
# [1] 18
# |--------------------------------------------------|
#   |==================================================|
#   [1] "N variants same ref and alt: 291645"
# [1] "N variants swapped ref and alt: 44107"
# [1] "N variants same reference but swapped alt: 191"
# [1] "N variants swapped reference different alt: 30"
# [1] "N variants same position: 153"
# [1] "N variants not in iibdgc meta: 140818"
# 
# TRUE 
# 1 
# [1] 19
# [1] "N variants same ref and alt: 227270"
# [1] "N variants swapped ref and alt: 31541"
# [1] "N variants same reference but swapped alt: 172"
# [1] "N variants swapped reference different alt: 25"
# [1] "N variants same position: 104"
# [1] "N variants not in iibdgc meta: 91175"
# 
# TRUE 
# 1 
# [1] 20
# [1] "N variants same ref and alt: 228850"
# [1] "N variants swapped ref and alt: 32567"
# [1] "N variants same reference but swapped alt: 135"
# [1] "N variants swapped reference different alt: 18"
# [1] "N variants same position: 94"
# [1] "N variants not in iibdgc meta: 147892"
# 
# TRUE 
# 1 
# [1] 21
# [1] "N variants same ref and alt: 135436"
# [1] "N variants swapped ref and alt: 22544"
# [1] "N variants same reference but swapped alt: 90"
# [1] "N variants swapped reference different alt: 19"
# [1] "N variants same position: 65"
# [1] "N variants not in iibdgc meta: 48220"
# 
# TRUE 
# 1 
# [1] 22
# [1] "N variants same ref and alt: 141397"
# [1] "N variants swapped ref and alt: 20087"
# [1] "N variants same reference but swapped alt: 105"
# [1] "N variants swapped reference different alt: 17"
# [1] "N variants same position: 72"
# [1] "N variants not in iibdgc meta: 64918"
# 
# TRUE 
# 1 
# 
# [1] "cd"
# [1] "X"
# [1] "N variants same ref and alt: 4429"
# [1] "N variants swapped ref and alt: 878"
# [1] "N variants same reference but swapped alt: 1"
# [1] "N variants swapped reference different alt: 2"
# [1] "N variants same position: 2"
# [1] "N variants not in iibdgc meta: 581935"
# 
# TRUE 
# 1 
# 
# [1] "ibd"
# [1] "X"
# [1] "N variants same ref and alt: 4534"
# [1] "N variants swapped ref and alt: 883"
# [1] "N variants same reference but swapped alt: 1"
# [1] "N variants swapped reference different alt: 2"
# [1] "N variants same position: 2"
# [1] "N variants not in iibdgc meta: 581825"
# 
# TRUE 
# 1 
# [1] "uc"
# [1] "X"
# [1] "N variants same ref and alt: 4521"
# [1] "N variants swapped ref and alt: 884"
# [1] "N variants same reference but swapped alt: 1"
# [1] "N variants swapped reference different alt: 2"
# [1] "N variants same position: 2"
# [1] "N variants not in iibdgc meta: 581837"
# 
# TRUE 
# 1 




#########################
  
  
  