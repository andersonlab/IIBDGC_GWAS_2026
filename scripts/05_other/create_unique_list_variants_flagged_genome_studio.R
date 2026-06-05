# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
## create an unique list of variants flagged by Genome Studio:

library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"


br<-read.table("/path/to/ibdgwas/IIBDGC/from_kyle/Broad_flagged67262.txt",head=F)
dim(br)
# [1] 67262     1

fe<-read.table("/path/to/ibdgwas/IIBDGC/from_kyle/Feinstein_flagged54386.txt",head=F)
dim(fe)

all<-rbind(br,fe)
all<-all[!duplicated(all$V1),,drop=F]

dim(all)
# [1] 87859     1


## BROAD:

bim<-read.table(paste(path,"pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_hg19.bim",sep=""),head=F)
dim(br[which(br$V1 %in% bim$V2),,drop=F])
# 67262 2

dim(all[which(all$V1 %in% bim$V2),,drop=F])
# [1] 86874     1

all.1<-merge(all,bim,by.x="V1",by.y="V2",all.x=T)

## Feinstein:

bim<-read.table(paste(path,"pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_hg19.bim",sep=""),head=F)
dim(fe[which(fe$V1 %in% bim$V2),,drop=F])
# 54386    2

dim(all[which(all$V1 %in% bim$V2),,drop=F])
# [1] 87859     1

all.2<-merge(all,bim,by.x="V1",by.y="V2",all.x=T)


table(all.1$V5,all.1$V6)
#       -     A     C     D     G     I     T
# -  2509  7763  2249    62  7010   133   120
# A     0     0  8432     0 26167     0   152
# C     0  8340     0     0   169     0     0
# D     0     0     0     0     0   143     0
# G     0 23208   220     0     0     0     0
# I     0     0     0    59     0     0     0
# T     0   138     0     0     0     0     0

table(all.2$V5,all.2$V6)
#       -     A     C     D     G     I     T
# -  2187  8427  2243    70  7158   161   126
# A     0     0  8563     0 26442     0   172
# C     0  8326     0     0   205     0     0
# D     0     0     0     0     0   168     0
# G     0 23130   258     0     0     0     0
# I     0     0     0    83     0     0     0
# T     0   140     0     0     0     0     0

colnames(all.1)[2:ncol(all.1)]<-paste(colnames(all.1)[2:ncol(all.1)],"_br",sep="")
colnames(all.2)[2:ncol(all.2)]<-paste(colnames(all.2)[2:ncol(all.2)],"_fe",sep="")

all<-merge(all.1,all.2,by="V1")
dim(all)
# [1] 87859    11

table(all$V1.y_fe==all$V1.y_br)# same chr
# TRUE 
# 86874 

all$chr<-NA
all$pos<-NA
all$Ref<-NA
all$Alt<-NA


for (i in 1:nrow(all)) {
  
  all$chr[i]<-all$V1.y_fe[i]
  all$pos[i]<-all$V4_fe[i]
  
  tmp<-c(as.character(all$V5_br[i]),as.character(all$V5_fe[i]))
  tmp<-tmp[which(tmp!="-")]
  all$Ref[i]<-paste(tmp[!duplicated(tmp)],collapse="|")
  
  tmp<-c(as.character(all$V6_br[i]),as.character(all$V6_fe[i]))
  tmp<-tmp[which(tmp!="-")]
  all$Alt[i]<-paste(tmp[!duplicated(tmp)],collapse="|")
  
  if (nchar(all$Ref[i])>1 | nchar(all$Alt[i])>1) {
    tmp<-c(as.character(all$V5_br[i]),as.character(all$V5_fe[i]),as.character(all$V6_br[i]),as.character(all$V6_fe[i]))
    tmp<-tmp[which(tmp!="-")]
    tmp<-tmp[!duplicated(tmp)]
    if(length(tmp)>2) {
      print(i)
    }
    all$Ref[i]<-tmp[1]
    all$Alt[i]<-tmp[2]
  }
  
}

table(all$Ref,all$Alt)

all<-all[,c(1,12:15)]

write.table(all,"/path/to/ibdgwas/IIBDGC/from_kyle/list_variants_flagged_genome_studio_b37.tsv",col.names=T,row.names=F,sep="\t")
all<-read.table("/path/to/ibdgwas/IIBDGC/from_kyle/list_variants_flagged_genome_studio_b37.tsv",head=T)

table(all$chr)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 7387 5886 6920 6065 6338 5074 6129 4405 3490 3736 3401 3535 3489 2634 2447 2431 
# 16   17   18   19   20   21   22   25 
# 2561 2467 2277 2638 1682 1497 1320   50 

## chr25 is XY region, recode:
all$chr[which(all$chr==25)]<-"23"

dim(all[which(all$chr %in% 1:23),])
# [1] 80472     5
dim(all)
# [1] 87859     5

all<-all[which(all$chr %in% 1:23),]

all$chromStart<-all$pos
all$chromEnd<-all$pos+1

all$chromStart<-format(all$chromStart, scientific=F)
all$chromEnd<-format(all$chromEnd, scientific=F)

dim(all)
# [1] 80472      7

all$chro<-paste("chr",all$chr,sep="")
all$chro[which(all$chro=="chr23")]<-"chrX"

all$chr<-as.numeric(all$chr)
all<-all[order(all$chr,all$pos,decreasing=F),]

write.table(all[,c(8,6,7,1)],"/path/to/ibdgwas/IIBDGC/from_kyle/list_variants_flagged_genome_studio_b37.bed",col.names=F,row.names=F,quote=F,sep="\t")


#### lift positions

${path_gwas}previous_qced_b38/liftover/liftOver \
/path/to/ibdgwas/IIBDGC/from_kyle/list_variants_flagged_genome_studio_b37.bed \
${path_gwas}previous_qced_b38/liftover/hg19ToHg38.over.chain.gz \
/path/to/ibdgwas/IIBDGC/from_kyle/list_variants_flagged_genome_studio_b37_lifted_hg38 \
/path/to/ibdgwas/IIBDGC/from_kyle/list_variants_flagged_genome_studio_b37_no_lifted_hg38


wc -l /path/to/ibdgwas/IIBDGC/from_kyle/list_variants_flagged_genome_studio_b37_no_lifted_hg38
# 98 /path/to/ibdgwas/IIBDGC/from_kyle/list_variants_flagged_genome_studio_b37_no_lifted_hg38
# includes lines with comments, N variants = 98/2

#####

bed<-read.table("/path/to/ibdgwas/IIBDGC/from_kyle/list_variants_flagged_genome_studio_b37_lifted_hg38")
dim(bed)
# [1] 80423     4

dim(all)
# [1] 80472     8

98/2
# [1] 49

im(all)-dim(bed)
# [1] 49  4 ### OK

bed<-merge(bed,all,by.x="V4",by.y="V1",all.x=T)
table(bed$V1)
# chr1                   chr10                   chr11 
# 5746                    3398                    3535 
# chr12                   chr13                   chr14 
# 3489                    2634                    2442 
# chr14_GL000009v2_random                   chr15                   chr16 
# 56                    2427                    2561 
# chr17                   chr18                   chr19 
# 2466                    2277                    2630 
# chr19_KI270938v1_alt  chr1_KI270706v1_random     chr1_KI270766v1_alt 
# 4                      18                       2 
# chr2                   chr20                   chr21 
# 6912                    1682                    1499 
# chr22     chr2_KI270894v1_alt                    chr3 
# 1319                       2                    6064 
# chr4  chr4_GL000008v2_random                    chr5 
# 6337                      42                    5074 
# chr6                    chr7     chr7_KI270803v1_alt 
# 6128                    4377                      24 
# chr8     chr8_KI270821v1_alt                    chr9 
# 3486                       4                    3717 
# chrUn_KI270742v1                    chrX 
# 21                      50 

56+4+18+2+2+42+24+4+21
# [1] 173


bed<-bed[!grepl("_",bed$V1),]
dim(bed)
# [1] 80250     5

bed<-bed[,c(1,5,3,7,8)]

colnames(bed)<-c("var_ID","chr","pos","A0","A1")
write.table(bed,"/path/to/ibdgwas/IIBDGC/from_kyle/list_variants_flagged_genome_studio_b38.tsv",col.names=F,row.names=F,quote=F,sep="\t")





