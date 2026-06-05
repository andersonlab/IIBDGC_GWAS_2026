# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# test
# # singularity exec iibdgc_postprocess_10_singularity.sif
# MEM=2000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

# How to submit:

# cd /path/to/project
# files_ccre=($(ls | grep bed | sed 's/\.bed\.gz//g' ))

# MEM=1200
# for i in ${files_ccre[@]}
# do
# bsub -J"ccre" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 2 \
# -e ${path_gwas}post_imputation/log/${i}_ccre_stderr \
# -o ${path_gwas}post_imputation/log/${i}_ccre_stdout \
# "Rscript ~/git/IIBDGC_GWAS/scripts/other/split_ccre_annotation_per_type.R ${i} > \
# ${path_gwas}post_imputation/2022/log/split_ccre_annotation_per_type${i}.Rout"
# done



# MORE INFO ABOUT THE FORMAT
# https://omni-c.readthedocs.io/en/latest/epigenetics.html
# https://www.encodeproject.org/hic/

# .bedpe	loops	Pairs of loci that show significantly closer proximity with one another than with the loci lying between them	
# May include multiple files with this output depending on if the file was derived_from mapq1 or mapq30.
# The preferred files to use for analysis are those with mapq30 annotations. They are marked with preferred_default=true property. 

library(data.table)

# provide trait (in this case combination of database + trait) as input
args<-commandArgs(trailingOnly=TRUE)
bedpe_id<-args[1]

# for testing purposes:
# bedpe_id<-"ENCFF532RZR"
# bedpe_id<-"ENCFF992XHP"


path_files<-"/path/to/project"
path_gwas<-"/path/to/ibdgwas/IIBDGC/"

dat<-fread(paste0(path_files,bedpe_id,".bedpe.gz"))

# keep only loop files:

if (ncol(dat)==33) {

    colnames(dat)<-c("chr","x1","x2","chr","x1","x2","name","score","strand1","strand2","color","observed","expectedBL",
    "expectedDonut","expectedH","expectedV","fdrBL","fdrDonut","fdrH","fdrV","numCollapsed","centroid1","centroid2","radius","highRes_start_1",
    "highRes_end_1","highRes_start_2","highRes_end_2","localX","localY","localObserved","localPval","localPeakID")

    bed<-rbind(dat[,c(1:3)],dat[,c(4:6)])

    fwrite(bed,paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/annotation/encode_bedfiles/intact_hic/bed_b38/",bedpe_id,"_intact_hic.bed.gz"),col.names=F,row.names=F,quote=F,sep="\t")

}

q("no")
