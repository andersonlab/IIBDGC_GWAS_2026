# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
MEM=5000
bsub -Is -m "modern_hardware" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group -q yesterday /software/bin/R-4.3.1

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)
library(ggbio)
library(biovizBase)
library(wesanderson)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("karyoploteR")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

library(karyoploteR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

path<-"/path/to/ibdgwas/IIBDGC/"

all<-fread(paste(path,"summary_files/summary_stats_new_index_variants_accross_399_regions_post_forward_regression.tsv",sep=""),head=T)
new<-all[which(all$MarkerName=="chr7:117252792:T:C"),]



# set intervals
min_pos<-new$min
max_pos<-new$max

# set Gene ID
gene_name<-new$closest_gene

# set Chr
chr<-new$chr

# if (new$phenotype=="CD") {
#   pheno<-"cd"
#   colour_plot<-c("#004488")
# }
# if (new$phenotype=="UC") {
#   pheno<-"uc"
#   colour_plot<-c("#BB5566")
# }
# if (new$phenotype=="IBD_unsaturated") {
#   pheno<-"ibd"
#   colour_plot<-c("#db7107")
# }
# if (new$phenotype=="IBD_saturated") {
#   pheno<-"ibd"
#   colour_plot<-c("#db7107")
# }


pheno<-"ibd"
colour_plot<-c("#db7107")



  pheno<-"cd"
  colour_plot<-c("#004488")

risk<-fread(paste(path,"post_imputation/analysis/metaanalysis/",pheno,"/",new$chr,"_",pheno,"_meta_noGC_PCs_firthse_no_illumina550_with_danish_ukbb_decode_finngen_1_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq.txt.gz",sep=""),head=T)
risk$N<-rowSums(risk[,c("N_sample","danish_N_sample","ukbb_N_sample","finngen_N_sample","decode_N_sample")],na.rm=T)
risk$rate<-rowSums(risk[,c("N_sample","danish_N_sample","ukbb_N_sample","finngen_N_sample","decode_N_sample")],na.rm=T)/max(risk$N)
risk$rate_noukbb<-rowSums(risk[,c("N_sample","danish_N_sample","finngen_N_sample","decode_N_sample")],na.rm=T)/(max(risk$N)-361784)

risk$position<-gsub("chr[0-9]{1,2}:","",risk$MarkerName)
risk$position<-as.numeric(gsub(":[A-Z]*:[A-Z]*","",risk$position))


risk_2<-fread(paste(path,"post_imputation/analysis/metaanalysis/",pheno,"/",new$chr,"_meta_noGC_PCs_firthse_no_illumina550_with_danish_no_ukbb_with_decode_with_finngen_1.txt.gz",sep=""),head=T)
risk_2<-risk_2[,c("MarkerName","P-value")]
colnames(risk_2)[2]<-"P-value_noUKBB"

risk<-merge(risk,risk_2,by="MarkerName",all.x=T)

risk$pvalue<-risk$'P-value_noUKBB'

risk<-risk[which(risk$position>=new$min & risk$position<=new$max),]
risk$chr<-paste("chr",chr,sep="")
rm(risk_2)





# dim(risk[which(risk$rate_noukbb>=0.8),])
# # [1] 2410   39
# 
# dim(risk[which(risk$rate_noukbb>=0.7),])
# # [1] 2410   39
# dim(risk[which(risk$rate_noukbb>=0.6),])
# # [1] 2410   39

# dim(risk[which(risk$rate_noukbb>=0.5),])
# # [1] 3677   39

risk[which(risk$MarkerName=='chr7:117254861:A:T'),]

risk<-risk[which(risk$rate_noukbb>=0.5),]

risk$variant<-risk$MarkerName

risk<-as.data.frame(risk)
risk$pval<-risk$"P-value_noUKBB"

risk<-risk[,c("chr","position","variant","pval")]

risk$strand<-"+"

risk<-risk[order(risk$position,decreasing=F),]

# convert to data ranges object
gr.risk <- transformDfToGr(risk, seqnames = "chr", start = "position", end = "position", strand="strand")


# set max log pvalue
max_sig<-ceiling(max(-log10(risk$pval)) )


# update ranges if needed for aestetics

if(gene_name=="SHARPIN") {
  min_pos<-143700000
  max_pos<-144700000
}

# CARLMIL2, but closest to most signif is AGRP
if(gene_name=="AGRP") {
  min_pos<-67200000
  max_pos<-68200000
}

if(gene_name=="TSPAN32") {
  min_pos<-1352842
  max_pos<-2800421
}


pdf(paste("~/tmp_plots/chr",chr,"_",min_pos,"_",max_pos,"_",pheno,"_no_ukbb_regional_manhattan_plot.pdf",sep=""),height = 7, width = 10)

#Draw
kp <- plotKaryotype(genome="hg38",main=gene_name,plot.type=4,zoom=paste("chr",chr,":",min_pos,"-",max_pos,sep=""))

# plot 1
r0_pos<-0.45
r1_pos<-1

kpAddBaseNumbers(kp, add.units = TRUE, cex=1, tick.dist = 1e5)
kpPlotManhattan(kp, data=gr.risk, suggestiveline = 0, genomewideline = 0,r0=r0_pos, r1=r1_pos, ymax=max_sig, points.col=colour_plot)
kpAddLabels(kp, labels = paste(toupper(pheno)," risk -log10(Pval)",sep=""), label.margin = 0.043,r0=0.7, r1=r1_pos,srt=90,cex=0.8)
kpAxis(kp, ymin = 0, ymax=max_sig, numticks = max_sig+1, r0=r0_pos, r1=r1_pos,cex=0.8)

# plot 3
r0_pos<-0.0
r1_pos<-0.4

genes.data <- makeGenesDataFromTxDb(txdb = TxDb.Hsapiens.UCSC.hg38.knownGene, karyoplot = kp)
genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)
kpPlotGenes(kp, data=genes.data, add.transcript.names = FALSE, cex=0.7, gene.name.position = "right", r0=r0_pos, r1=r1_pos)

#Close the device
dev.off()






# # get gnomad data:
# 
# gnomad<-fread(paste(path,"resources/gnomad/gnomad_freq_edited_chr",new$chr,".tsv.gz",sep=""),head=F)
# colnames(gnomad)<-c("SNP","CHROM","POS","REF","ALT","AF","AF_nfe","AF_afr","AF_amr","AF_eas","AF_sas","AF_asj")
# class(gnomad$POS)
# # [1] "integer"
# 
# gnomad<-gnomad[which(gnomad$POS>=new$min & gnomad$POS<=new$max),]
# gnomad<-gnomad[which(gnomad$AF_nfe>=0.001 & gnomad$AF_nfe<=0.999),]
# dim(gnomad)
# # [1] 3768   12
# 
# 
# gnomad$MarkerName<-gsub("_",":",gnomad$SNP)
# gnomad$MarkerName<-paste("chr",gnomad$MarkerName,sep="")
# 
# dim(gnomad[which(gnomad$MarkerName %in% risk$MarkerName),])
# # [1] 3289   13
# 
# dim(gnomad[which(gnomad$MarkerName %in% risk$MarkerName[which(risk$rate_noukbb>=0.5)]),])
# # [1] 2575   13
# 
# 2575/3768
# # [1] 0.6833864
# 
# 




