# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
####################################
## DRAW A REGIONAL MANHATTAM PLOT:

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=4000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R


library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(ggbio)
library(biovizBase)
library(wesanderson)
# library(AnnotationHub)
library(karyoploteR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# hub <- AnnotationHub()
# query(hub, c("homo","orgdb"))
# file.copy(AnnotationHub::cache(hub["AH114084"]), "./org.Hs.eg.db.sqlite")
# makeAnnDbPkg(seed, "org.Hs.eg.db.sqlite")


pheno<-c("ibd","cd","uc")

path_gwas<-"/path/to/ibdgwas/IIBDGC/"
all<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_new_and_old_index_gwas_variants_conditional_on_known_index_gwas_fm_2023_with_iibdgc_tier1_tier2_eas_sas_summary_stats_plus_unsupervised_conditional.tsv",sep=""))
all<-as.data.frame(all)
dim(all)


new<-all[which(all$dbsnp154_eur_tier_1 %in% c("rs1887428")),]

# set Chr
chr<-new$chr

# set associated phenotype
# pheno<-c("ibd","cd","uc")
# cols<-c("#db7107","#004488","#BB5566")

colour_plot<-c("#004488")

# set intervals
# min_pos<-new$min
# max_pos<-new$max

min_pos<-new$pos-500000
max_pos<-new$pos+500000


# risk<-fread(paste(path,"post_imputation/analysis/metaanalysis/",pheno,"/",new$chr,"_",pheno,"_meta_noGC_PCs_firthse_no_illumina550_with_danish_ukbb_decode_finngen_1_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""),head=T)
i=2
risk_cd<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr,"_",pheno[i],"_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""),head=T)
risk_cd<-risk_cd[which(risk_cd$rate_total_sample>=0.4)]
risk_cd$position<-gsub("chr[0-9]{1,2}:","",risk_cd$MarkerName)
risk_cd$position<-as.numeric(gsub(":[A-Z]*:[A-Z]*","",risk_cd$position))
risk_cd$pvalue<-risk_cd$'P-value'
risk_cd<-risk_cd[which(risk_cd$position>=min_pos & risk_cd$position<=max_pos),]
risk_cd$chr<-paste("chr",chr,sep="")
risk_cd$variant<-risk_cd$MarkerName
risk_cd$pval<-risk_cd$'P-value'
risk_cd<-as.data.frame(risk_cd[,c("chr","position","variant","pval")])
risk_cd$strand<-"+"
head(risk_cd[order(risk_cd$pval,decreasing=F),])

i=3
risk_uc<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr,"_",pheno[i],"_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""),head=T)
risk_uc<-risk_uc[which(risk_uc$rate_total_sample>=0.39)]
risk_uc$position<-gsub("chr[0-9]{1,2}:","",risk_uc$MarkerName)
risk_uc$position<-as.numeric(gsub(":[A-Z]*:[A-Z]*","",risk_uc$position))
risk_uc$pvalue<-risk_uc$'P-value'
risk_uc<-risk_uc[which(risk_uc$position>=min_pos & risk_uc$position<=max_pos),]
risk_uc$chr<-paste("chr",chr,sep="")
risk_uc$variant<-risk_uc$MarkerName
risk_uc$pval<-risk_uc$'P-value'
risk_uc<-as.data.frame(risk_uc[,c("chr","position","variant","pval")])
risk_uc$strand<-"+"
head(risk_uc[order(risk_uc$pval,decreasing=F),])

i=1
risk_ibd<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/",chr,"_",pheno[i],"_meta_eur_tier_2_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""),head=T)
risk_ibd<-risk_ibd[which(risk_ibd$rate_total_sample>=0.4)]
risk_ibd$position<-gsub("chr[0-9]{1,2}:","",risk_ibd$MarkerName)
risk_ibd$position<-as.numeric(gsub(":[A-Z]*:[A-Z]*","",risk_ibd$position))
risk_ibd$pvalue<-risk_ibd$'P-value'
risk_ibd<-risk_ibd[which(risk_ibd$position>=min_pos & risk_ibd$position<=max_pos),]
risk_ibd$chr<-paste("chr",chr,sep="")
risk_ibd$variant<-risk_ibd$MarkerName
risk_ibd$pval<-risk_ibd$'P-value'
risk_ibd<-as.data.frame(risk_ibd[,c("chr","position","variant","pval")])
risk_ibd$strand<-"+"
head(risk_ibd[order(risk_ibd$pval,decreasing=F),])

# set max log pvalue
max_sig<-ceiling(max(-log10(risk_ibd$pval),-log10(risk_cd$pval),-log10(risk_uc$pval)) )

# convert to data ranges object
gr.risk_cd <- transformDfToGr(risk_cd, seqnames = "chr", start = "position", end = "position", strand="strand")
gr.risk_ibd <- transformDfToGr(risk_ibd, seqnames = "chr", start = "position", end = "position", strand="strand")
gr.risk_uc <- transformDfToGr(risk_uc, seqnames = "chr", start = "position", end = "position", strand="strand")

# set max log pvalue
max_sig<-ceiling(-log10(min(as.numeric(risk_cd$pval),as.numeric(risk_uc$pval),as.numeric(risk_ibd$pval))))

number_ticks<-max_sig+1
if(max_sig>10) {
  max_sig<-ceiling(max_sig/10)*10
  number_ticks<-1+(max_sig/10)
}


# cs<-read.table(paste("/path/to/ibdgwas/IIBDGC/post_imputation/analysis/finemapping/ibd/chr",chr,"_",min_pos,"_",max_pos,"_all_arrays.cred1",sep=""),head=T)
# # cs<-read.table(paste("/path/to/ibdgwas/IIBDGC/post_imputation/analysis/finemapping/ibd/chr",chr,"_",min_pos,"_",max_pos,"_all_arrays.cred2",sep=""),head=T)
# cred_set<-gr.risk[(elementMetadata(gr.risk)[, "variant"] %in% c(as.character(cs$cred1),as.character(cs$cred2)) )]

# cs_list<-c("chr16:50270774:G:T","chr16:50301163:C:A")

all<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_new_and_old_index_gwas_variants_conditional_on_known_index_gwas_fm_2023_with_iibdgc_tier1_tier2_eas_sas_summary_stats_plus_unsupervised_conditional.tsv",sep=""))
all<-as.data.frame(all)

cs_list<-all$MarkerName[which(all$updated_region==new$updated_region)]

# cs_list<-new$MarkerName
cred_set_cd<-gr.risk_cd[(elementMetadata(gr.risk_cd)[, "variant"] %in% cs_list)]
cred_set_uc<-gr.risk_uc[(elementMetadata(gr.risk_uc)[, "variant"] %in% cs_list)]
cred_set_ibd<-gr.risk_ibd[(elementMetadata(gr.risk_ibd)[, "variant"] %in% cs_list)]

col_cd<-"#004488"
col_ibd<-"#db7107"
col_uc<-"#BB5566"


pdf(paste("~/git/IIBDGC_GWAS/plots/chr",chr,"_",min_pos,"_",max_pos,"_",pheno,"_regional_manhattan_plot_with_danish_ukbb_finngen_decode_alt_highlight_credible_set.pdf",sep=""),height = 10, width = 10)

kp <- plotKaryotype(genome="hg38",main=new$closest_gene,plot.type=4,zoom=paste("chr",chr,":",min_pos,"-",max_pos,sep=""))

# PLOT 1 - CD
  
r0_pos<-0.77
r1_pos<-1

kpAddBaseNumbers(kp, add.units = TRUE, cex=1, tick.dist = 1e5)
kpPlotManhattan(kp, data=gr.risk_cd, suggestiveline = 0, genomewideline = 0,r0=r0_pos, r1=r1_pos, ymax=max_sig, points.col="grey",highlight=cred_set_cd,highlight.col=col_cd)
kpAddLabels(kp, labels = "CD -log10(Pval)", label.margin = 0.043,r0=r1_pos, r1=r1_pos,srt=90,cex=0.8)
kpAxis(kp, ymin = 0, ymax=max_sig, numticks = number_ticks, r0=r0_pos, r1=r1_pos,cex=0.8)

  
# PLOT 2 - UC

r0_pos<-0.52
r1_pos<-0.73

kpAddBaseNumbers(kp, add.units = TRUE, cex=1, tick.dist = 1e5)
kpPlotManhattan(kp, data=gr.risk_uc, suggestiveline = 0, genomewideline = 0,r0=r0_pos, r1=r1_pos, ymax=max_sig, points.col="grey",highlight=cred_set_uc,highlight.col=col_uc)
kpAddLabels(kp, labels = "UC -log10(Pval)", label.margin = 0.043,r0=r1_pos, r1=r1_pos,srt=90,cex=0.8)
kpAxis(kp, ymin = 0, ymax=max_sig, numticks = number_ticks, r0=r0_pos, r1=r1_pos,cex=0.8)

# PLOT 3 - IBD

r0_pos<-0.27
r1_pos<-0.48

kpAddBaseNumbers(kp, add.units = TRUE, cex=1, tick.dist = 1e5)
kpPlotManhattan(kp, data=gr.risk_ibd, suggestiveline = 0, genomewideline = 0,r0=r0_pos, r1=r1_pos, ymax=max_sig, points.col="grey",highlight=cred_set_ibd,highlight.col=col_ibd)
kpAddLabels(kp, labels = "IBD -log10(Pval)", label.margin = 0.043,r0=r1_pos, r1=r1_pos,srt=90,cex=0.8)
kpAxis(kp, ymin = 0, ymax=max_sig, numticks = number_ticks, r0=r0_pos, r1=r1_pos,cex=0.8)


# PLOT 4 - Genes
  
r0_pos<-0.0
r1_pos<-0.20

GenesData <- makeGenesDataFromTxDb(txdb = TxDb.Hsapiens.UCSC.hg38.knownGene, karyoplot = kp,plot.transcripts=TRUE)
GenesData <- addGeneNames(GenesData,keytype="ENTREZID", names="SYMBOL")
GenesData <- mergeTranscripts(GenesData)
kpPlotGenes(kp, data=GenesData, add.transcript.names = FALSE, cex=0.7, gene.name.position = "right", r0=r0_pos, r1=r1_pos)
  

#Close the device
dev.off()



