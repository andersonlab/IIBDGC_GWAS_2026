# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
####################################
## DRAW A REGIONAL MANHATTAM PLOT:

# # how to submit
# # singularity exec iibdgc_postprocess_10_singularity.sif

# MEM=8000
# release="eur_eas_sas_tier_2"
# region="1_780044_3066589"

# bsub -J"plot_imp1" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -m "modern_hardware" -G ibdgwas \
# -e ${path_gwas}scripts/logs/draw_manhattan_plots_stderr_${release} \
# -o ${path_gwas}scripts/logs/draw_manhattan_plots_stdout_${release} \
# "Rscript ~/git/IIBDGC_GWAS/scripts/other/draw_regional_manhattan_plot.R ${release} ${region} > \
# ${path_gwas}scripts/logs/draw_manhattan_plots_label_known_regions_${release}.Rout"

# MEM=6000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

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

rm(list=ls())

args <- commandArgs()
release<-args[6]
region<-args[7]

# options for release are: 
#   release<-"eur_tier_1"
#   release<-"eur_tier_2"
#   release<-"sas_tier_2"
#   release<-"eur_eas_sas_tier_2"

# options for regions:
# region<-"1_52976954_54020750"
# region<-"1_154561096_156553012"
# region<-"5_2840882_3846018"

pheno<-c("ibd","cd","uc")
path_gwas<-"/path/to/ibdgwas/IIBDGC/"

# set Chr
chr<-gsub("_.*","",region)

# set intervals
min_pos<-gsub("^[0-9]{1,2}_","",region)
max_pos<-as.numeric(gsub("^[0-9]*_","",min_pos))
min_pos<-as.numeric(gsub("_[0-9]*","",min_pos))

print(region)
print(release)

# extract summary statistics:
for (ph in pheno) {

  print(ph)

  if (release=="eur_tier_1") {
    risk<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",ph,"/",chr,"_",ph,"_meta_",release,"_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff_with_rsid_dbsnp154.txt.gz",sep=""),head=T)
  } else {
    risk<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",ph,"/",chr,"_",ph,"_meta_",release,"_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_withNeff.txt.gz",sep=""),head=T)
  }

  risk<-risk[which(risk$rate_Neff>=0.5)]

  # if (ph %in% c("cd","uc")) {
    risk<-risk[which(risk$HetPVal>1E-15)]
  # }

  risk$position<-gsub("chr[0-9]{1,2}:","",risk$MarkerName)
  risk$position<-as.numeric(gsub(":[A-Z]*:[A-Z]*","",risk$position))
  risk$pvalue<-risk$'P-value'
  risk<-risk[which(risk$position>=min_pos & risk$position<=max_pos),]
  risk$chr<-paste("chr",chr,sep="")
  risk$variant<-risk$MarkerName
  risk$pval<-as.numeric(risk$'P-value')
  risk<-as.data.frame(risk[,c("chr","position","variant","pval")])
  risk$strand<-"+"
  print(head(risk[order(risk$pval,decreasing=F),]))

  assign(paste0("risk_",ph),risk)
  
  rm(risk)

}

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

# define variants to highlight:
all<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_new_and_old_index_gwas_variants_conditional_on_known_index_gwas_fm_2023_with_iibdgc_tier1_tier2_eas_sas_summary_stats_plus_unsupervised_conditional_plus_new_fm_plus_pheno.tsv",sep=""),head=T)
all<-all[which(all$updated_region==region)]

cs_list<-all$MarkerName

# cs_list<-new$MarkerName
cred_set_cd<-gr.risk_cd[(elementMetadata(gr.risk_cd)[, "variant"] %in% cs_list)]
cred_set_uc<-gr.risk_uc[(elementMetadata(gr.risk_uc)[, "variant"] %in% cs_list)]
cred_set_ibd<-gr.risk_ibd[(elementMetadata(gr.risk_ibd)[, "variant"] %in% cs_list)]

col_cd<-"#004488"
col_ibd<-"#db7107"
col_uc<-"#BB5566"


pdf(paste("~/git/IIBDGC_GWAS/plots/regional_manhattan_plots/chr",chr,"_",min_pos,"_",max_pos,"_",release,"_regional_manhattan_plot_highlight_index.pdf",sep=""),height = 10, width = 10)

kp <- plotKaryotype(genome="hg38",main=all$updated_region[1],plot.type=4,zoom=paste("chr",chr,":",min_pos,"-",max_pos,sep=""))

# PLOT 1 - CD
  
r0_pos<-0.77
r1_pos<-1

kpAddBaseNumbers(kp, add.units = TRUE, cex=1, tick.dist = 1e5)
if (length(cred_set_cd$chr)>0) {
  kpPlotManhattan(kp, data=gr.risk_cd, suggestiveline = 0, genomewideline = 0,r0=r0_pos, r1=r1_pos, ymax=max_sig, points.col="grey",highlight=cred_set_cd,highlight.col=col_cd)
} else {
  kpPlotManhattan(kp, data=gr.risk_cd, suggestiveline = 0, genomewideline = 0,r0=r0_pos, r1=r1_pos, ymax=max_sig, points.col="grey")
}
kpAddLabels(kp, labels = "CD -log10(Pval)", label.margin = 0.043,r0=r1_pos, r1=r1_pos,srt=90,cex=0.8)
kpAxis(kp, ymin = 0, ymax=max_sig, numticks = number_ticks, r0=r0_pos, r1=r1_pos,cex=0.8)

  
# PLOT 2 - UC

r0_pos<-0.52
r1_pos<-0.73

kpAddBaseNumbers(kp, add.units = TRUE, cex=1, tick.dist = 1e5)
if(length(cred_set_uc$chr)>0) {
  kpPlotManhattan(kp, data=gr.risk_uc, suggestiveline = 0, genomewideline = 0,r0=r0_pos, r1=r1_pos, ymax=max_sig, points.col="grey",highlight=cred_set_uc,highlight.col=col_uc)
} else {
  kpPlotManhattan(kp, data=gr.risk_uc, suggestiveline = 0, genomewideline = 0,r0=r0_pos, r1=r1_pos, ymax=max_sig, points.col="grey")
}

kpAddLabels(kp, labels = "UC -log10(Pval)", label.margin = 0.043,r0=r1_pos, r1=r1_pos,srt=90,cex=0.8)
kpAxis(kp, ymin = 0, ymax=max_sig, numticks = number_ticks, r0=r0_pos, r1=r1_pos,cex=0.8)

# PLOT 3 - IBD

r0_pos<-0.27
r1_pos<-0.48

kpAddBaseNumbers(kp, add.units = TRUE, cex=1, tick.dist = 1e5)
if(length(cred_set_ibd$chr)>0) {
  kpPlotManhattan(kp, data=gr.risk_ibd, suggestiveline = 0, genomewideline = 0,r0=r0_pos, r1=r1_pos, ymax=max_sig, points.col="grey",highlight=cred_set_ibd,highlight.col=col_ibd)
} else {
  kpPlotManhattan(kp, data=gr.risk_ibd, suggestiveline = 0, genomewideline = 0,r0=r0_pos, r1=r1_pos, ymax=max_sig, points.col="grey")
}
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


q("no")
