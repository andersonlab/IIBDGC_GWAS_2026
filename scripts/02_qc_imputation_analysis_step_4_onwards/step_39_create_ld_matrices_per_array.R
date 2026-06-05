# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# Use ldstore to generate ld matrices for the different genotyping arrays:

# /path/to/software/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64


# see: http://www.christianbenner.com/#

# # Command-line arguments
# # 
# # --bcor-file		Option to specify a BCOR file		With --bcor-to-text
# # --bcor-to-text		Convert BCOR file to a text file		Subprogram
# # --bdose-version	  	Option to set the BDOSE file version (see below)	  	With --write-bdose
# # --compression	  	Option to specify the compression level (see below) of a BDOSE/BCOR file as 'ultra-low' (1 byte), 'low' (2 bytes), 'medium' (4 bytes) or 'high' (8 bytes)	  	Default is medium. With --write-bcor and --write-bdose
# # --dataset		Option to specify a delimiter-separated list of datasets as given in the master file (e.g. 1,2 or 1|2)		All datasets are processed by default
# # --in-files	  	Master file (see below)	  	With --bcor-to-text, --write-bcor, --write-bdose and --write-text
# # --ld-file		Option to specify a LD file		With --bcor-to-text
# # --memory		Option to limit the amount of memory in gigabyte that can be used during computation of SNP correlations		Default is 1Gb. With --read-bdose or --write-bdose when using either --write-bcor or --write-text
# # --n-threads		Option to set the number of parallel threads		Default is 1. With --write-bcor, --write-bdose or --write-text
# # --range		Option to specify a genomic range xx-yy to operate on when converting a BCOR file to a LD file where xx and yy are the start and end coordinates in base pairs		With --bcor-to-text
# # --read-bdose		Read dosage data from a BDOSE file		With --write-bcor or --write-text
# # --read-only-bgen		Read genotype probabilities from a BGEN file and store dosage data in memory		With --write-bcor or --write-text
# # --rsids		Option to sepcify a comma-separated list of SNP identifiers corresponding with the rsid column in a Z file (see below)		With --bcor-to-text or --write-text
# # --sample-miss		Option to set a missing data threshold between 0 and 1. If the missing data rate for a SNP is above the specified threshold, then the correlation of any SNP pair that includes this SNP is set to NA. If the missing data rate for a SNP is below the specified threshold, then missing data is mean-imputed		Default is 0.1. With --write-bcor or --write-text
# # --write-bcor		Write SNP correlations to a BCOR file		Subprogram
# # --write-bdose		Write dosage data to a BDOSE file		Subprogram and with --write-bcor or --write-text
# # --write-text		Write SNP correlations to a text file		Subprogram


# # # run test to understand inputs:
# # 
# # /path/to/software/bgen/gavinband-bgen-44fcabbc5c38/./build/apps/bgenix -g \
# # /path/to/software/ldstore_v2.0_x86_64/example/data.bgen \
# # -vcf  | head | cut -f1-10
# # 
# # # Welcome to bgenix
# # # (version: 1.1.6, revision )
# # # 
# # # (C) 2009-2017 University of Oxford
# # # 
# # # Building query                                              :  (55/?,0.0s,308801.4/s)
# # # ##fileformat=VCFv4.2
# # # ##FORMAT=<ID=GT,Type=String,Number=1,Description="Threshholded genotype call">
# # # ##FORMAT=<ID=GP,Type=Float,Number=G,Description="Genotype call probabilities">
# # # ##FORMAT=<ID=HP,Type=Float,Number=.,Description="Haplotype call probabilities">
# # # #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	(anonymous_sample_1)
# # # 01	1	rs1;---	A	G	.	.	.	GT:GP	0/1:0,1,0
# # # 01	2	rs2;---	A	G	.	.	.	GT:GP	0/1:0,1,0
# # # 01	3	rs3;---	A	G	.	.	.	GT:GP	1/1:0,0,1
# # # 01	4	rs4;---	A	G	.	.	.	GT:GP	1/1:0,0,1
# # # 01	5	rs5;---	A	G	.	.	.	GT:GP	0/0:1,0,0
# # 
# # 
# # # BGEN to BDOSE v1.1 file conversion
# # # Genotype data with 55 SNPs and 5363 individuals in BGEN format can be converted to dosage 
# # # data in BDOSE v1.1 format as follows.
# # cd
# # mkdir ~/test_ldose/
# # cd test_ldose/
# #   
# # cp /path/to/software/ldstore_v2.0_x86_64/example/data .
# # less data
# # # z;bgen;bgi;sample;bdose;bcor;ld;n_samples
# # # example/data.z;example/data.bgen;example/data.bgen.bgi;example/data.sample;example/data.bdose;example/data.bcor;example/data.ld;5363
# # 
# # mkdir example/
# # cd example/
# # cp /path/to/software/ldstore_v2.0_x86_64/example/data .
# # cp /path/to/software/ldstore_v2.0_x86_64/example/data.z .
# # cp /path/to/software/ldstore_v2.0_x86_64/example/data.bgen .
# # cp /path/to/software/ldstore_v2.0_x86_64/example/data.bgen.bgi .
# # cp /path/to/software/ldstore_v2.0_x86_64/example/data.sample .
# # 
# #   
# # /path/to/software/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 \
# # --in-files ~/test_ldose/example/data --write-bdose --bdose-version 1.1 --compression low
# # 
# # # Computation of SNP correlations
# # # SNP correlations for the same data as in the example above can be computed and written to a
# # # BCOR v1.1 file. 2) writing dosage data first to a BDOSE file
# # /path/to/software/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 \
# # --in-files example/data --write-bcor --write-bdose --bdose-version 1.1
# # 
# # # BCOR v1.1 to LD file conversion
# # /path/to/software/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 \
# # --in-files example/data --bcor-to-text


# ##################################################################################################################################################################
# ##################################################################################################################################################################

# # FIRST APPROACH, RELEASE LD MATRICES BY

# ###############################
# # 1. create bgen subset files #
# ###############################

# ##################################################################################################################################################################
# # 1.1 - create list of variants included in any analysis (ibd|cd|uc) - files cannot be subset when dosage is estimated, and bcor step requires lots of mem

# # singularity exec iibdgc_postprocess_10_singularity.sif

# ##################
# # MEM=20000
# # bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

# library(data.table)
# path_gwas<-"/path/to/ibdgwas/IIBDGC/"

# array<-c("illumina370","humanomniexpress","affymetrix500","humanomni1","quad610","illuminaexome","affymetrix6","humancoreexome","gsa")
# pheno<-c("ibd","cd","uc")

# regions<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_regions_425_cojo_supervised_plus_unsupervised.tsv.gz"))
# regions$min<-gsub("^[0-9]{1,2}_","",regions$updated_region)
# regions$max<-gsub("^[0-9]*_","",regions$min)
# regions$min<-gsub("_[0-9]*$","",regions$min)
# regions$chr<-gsub("_.*","",regions$updated_region)

# regions<-as.data.frame(regions)

# # create a list of SNPs per region:
# regions<-regions[!duplicated(regions),]

# table(regions$class)
# # new old 
# # 149 276

# # renames those regions:
# regions<-regions[,c("chr","updated_region","min","max")]
# regions<-as.data.frame(regions[!duplicated(regions),])

# regions<-regions[order(regions$min,decreasing="F"),]
# dim(regions)
# # [1] 425   4

# regions$min<-as.numeric(regions$min)
# regions$max<-as.numeric(regions$max)
# regions$chr<-as.numeric(regions$chr)

# # create a file with the chrs and regions to export later as env variables to submit the job:
# write.table(regions[,"updated_region",drop=F],paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/list_regions",sep=""),
#       col.names=F,row.names=F,quote=F,sep="\t")

# write.table(regions[,"chr",drop=F],paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/list_chromosomes",sep=""),
#       col.names=F,row.names=F,quote=F,sep="\t")



# for (i in 1:length(array)) {
  
#   print(array[i])
  
#   for (chr in c(1:22)){
  
#     print(chr)
    
#     for (j in 1:length(pheno)) {

#       print(pheno[j])

#       # if (chr==23){
#       #   chr<-"X"
#       # }
      
#       file<-paste(path_gwas,"post_imputation/2022/analysis/regenie/",array[i],"/",pheno[j],"/chr",chr,"_",array[i],"_eur_all_step2_",pheno[j],"_eur_sex_PCs_firthse_",pheno[j],".regenie.gz",sep="")
      
#       if (file.exists(file)) {

#         var_tokeep_tmp<-fread(file,sep=" ",head=T)
#         var_tokeep_tmp<-var_tokeep_tmp[which(var_tokeep_tmp$A1FREQ>=0.001 & var_tokeep_tmp$A1FREQ<=0.999 & var_tokeep_tmp$INFO>=0.4),c("ID","GENPOS")]
        
#         if(!exists("var_tokeep")) {
#           var_tokeep<-var_tokeep_tmp
#           print(nrow(var_tokeep))
#         } else {
#           var_tokeep_tmp<-var_tokeep_tmp[which(!var_tokeep_tmp$ID %in% var_tokeep$ID),c("ID","GENPOS")]
#           var_tokeep<-rbind(var_tokeep,var_tokeep_tmp)
#           print(nrow(var_tokeep))
#         }
        
#         rm(var_tokeep_tmp)
        
#       } else {
#         print("NO REGENIE FILE")
#       }

#     }

#     if (exists("var_tokeep")) {

#       var_tokeep<-as.data.frame(var_tokeep[order(var_tokeep$GENPOS,decreasing=F),])
#       regions_tmp<-as.data.frame(regions[which(regions$chr==chr),])

#       for (ii in 1:nrow(regions_tmp)) {

#         var_tokeep_tmp<-var_tokeep[which(var_tokeep$GENPOS>=regions_tmp$min[ii] & var_tokeep$GENPOS<=regions_tmp$max[ii]),]

#         write.table(var_tokeep_tmp[,"ID"],paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/",array[i],"/",array[i],"_chr",regions_tmp$updated_region[ii],"_list_variants_tokeep_in_analysis",sep=""),
#       col.names=F,row.names=F,quote=F,sep="\t")

#       }

#       rm(var_tokeep)

#     } else {
#       print(paste(chr,"file does not exist"))
#     }


#   }
# }

# q("no")


# ##################################################################
# # 1.2 - subset bgen files

# path_gwas=/path/to/ibdgwas/IIBDGC/
# array=(illumina370 humanomniexpress affymetrix500 humanomni1 quad610 illuminaexome affymetrix6 humancoreexome gsa)


# # export the env variables defined above
# mapfile -t region < <(cat /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/list_regions)
# mapfile -t chromosomes < <(cat /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/list_chromosomes)

# # export the env variables defined above - subset just the 6 new ones, plus the one that needs to be updated:
# mapfile -t region < <(cat /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/list_regions | grep -E "12_8422106_10257536|14_106197560_107197560|16_69127624_70127624|18_21986192_22986192|2_229716690_231482974|4_41896369_42896369|9_93754276_94754276")
# mapfile -t chromosomes < <(cat /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/list_chromosomes | grep -E "12|14|16|18|2|4|9")

# for i in ${array[@]}
# do
# for ii in {0..6} 
# do 
# ls -la ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${i}_chr${region[${ii}]}_list_variants_tokeep_in_analysis
# done
# done

# MEM=250

# for i in ${array[@]}
# do
# for ii in {0..6} 
# do 
# bsub -J"bgen_subset" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
# -o ${path_gwas}post_imputation/2022/log/${i}_chr${region[${ii}]}_bgen_subset_variants_stdout \
# -e ${path_gwas}post_imputation/2022/log/${i}_chr${region[${ii}]}_bgen_subset_variants_stderr \
# "/path/to/software/username/qctool_v2.2.0/./qctool \
# -g ${path_gwas}post_imputation/2022/${i}/imputed_data/${i}_chr${chromosomes[${ii}]}.dose.subset.bgen \
# -s ${path_gwas}post_imputation/2022/${i}/imputed_data/${i}_chr${chromosomes[${ii}]}.dose.subset.sample \
# -incl-rsids ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${i}_chr${region[${ii}]}_list_variants_tokeep_in_analysis \
# -os ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${i}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis.sample \
# -og ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${i}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis.bgen"
# done
# done


# for i in ${array[@]}; do ls -la ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${i}_chr*.dose.subset_list_variants_tokeep_in_analysis.bgen | wc -l; done
# for i in ${array[@]}; do for ii in {0..6}; do echo ${i} && tail -50 ${path_gwas}post_imputation/2022/log/${i}_chr${region[${ii}]}_bgen_subset_variants_stdout  | grep "Successfully"; done; done

# ############################
# # 1.3 - index files

# # Index bgen files

# MEM=100

# for i in ${array[@]}
# do
# for ii in {0..6}  
# do 
# bsub -J"bgi_subset" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
# -o ${path_gwas}post_imputation/2022/log/${i}_chr${region[${ii}]}_bgen_subset_variants_index_stdout \
# -e ${path_gwas}post_imputation/2022/log/${i}_chr${region[${ii}]}_bgen_subset_variants_index_stderr \
# "bgenix -index -g \
# ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${i}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis.bgen"
# done
# done

# for i in ${array[@]}; do ls -la ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${i}_chr*.dose.subset_list_variants_tokeep_in_analysis.bgen.bgi | wc -l; done
# for i in ${array[@]}; do for ii in {0..6}; do echo ${i} && tail -50 ${path_gwas}post_imputation/2022/log/${i}_chr${region[${ii}]}_bgen_subset_variants_index_stdout  | grep "Successfully"; done; done


# ########################################################
# # 1.4 - get the list of SNPs in the bgen file

# MEM=250

# # same variants between different phenotypes 
# for i in ${array[@]}
# do
# for ii in {0..6}   
# do 
# bsub -J"bgi" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -n 2 \
# -o ${path_gwas}post_imputation/2022/log/${i}_chr${region[${ii}]}_bgen_list_variants_stdout \
# -e ${path_gwas}post_imputation/2022/log/${i}_chr${region[${ii}]}_bgen_list_variants_stderr \
# "bgenix -g \
# ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${i}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis.bgen \
# -vcf | cut -f1-5 | gzip > \
# ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${i}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis.variants.gz"
# done
# done

# for i in ${array[@]}; do ls -la ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${i}_*.dose.subset_list_variants_tokeep_in_analysis.variants.gz | wc -l; done


# #################################################


# #################################
# # 2.- create z and master files #
# #################################

# # The master file is a semicolon-separated text file and contains no space. It contains the following mandatory 
# # column names and one dataset per line.
# # z column contains the names of Z files (input)
# # bgen column contains the names of BGEN files (input)
# # bgi column contains the names of BGI files (input)
# # bcor column contains the names of BCOR files when using --write-bcor (output)
# # ld column contains the names of LD files when using --bcor-to-text (output)
# # n_samples column contains the number of samples to include in any processing
# # bdose column contains the names of optional BDOSE files when using --write-bdose (optional output)
# # sample column contains the names of optional sample files when using --write-bdose and --bdose-version 1.1 (optional input)
# # incl column contains the names of optional sample inclusion files (optional input)

# # A master file with two datasets for writing SNP correlations to a BCOR or text file could look as follows.
# # z;bgen;bgi;bcor;ld;n_samples
# # dataset1.z;dataset1.bgen;dataset1.bgen.bgi;dataset1.bcor;;5363
# # dataset2.z;dataset2.bgen;dataset2.bgen.bgi;;dataset2.ld;5363

# # incl
# # The dataset.incl file is a text file to specify inclusion of samples in any processing. It contains one sample ID per line.


# # The dataset.z file is a space-delimited text file and contains meta information about the SNPs one SNP per line. 
# # It contains the mandatory column names in the following order.
# # A dataset.z file with three SNPs could look as follows.
# # rsid chromosome position allele1 allele2
# # rs1 10 1 T C
# # rs2 10 1 A G
# # rs3 10 1 G A

# # DOWNLOADED CENTROMERE POSITIONS TO SPLIT THE LD MATRICES INTO CHR ARMS
# https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=map&hgta_track=gap&hgta_table=gap&hgta_doSchema=describe+table+schema
# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/

  
  
# # array=(illumina370 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 gsa illuminaexome)
# # pheno=(ibd cd uc)
# # 
# # for i in ${array[@]}
# # do
# # for ph in ${pheno[@]}
# # do
# # mkdir -p ${i}/${ph}
# # done
# # done


# # MEM=4000
# # bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

# library(data.table)
# path_gwas<-"/path/to/ibdgwas/IIBDGC/"

# all<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_new_and_old_index_gwas_variants_conditional_on_known_index_gwas_fm_2023_with_iibdgc_tier1_tier2_eas_sas_summary_stats_plus_unsupervised_conditional.tsv",sep=""),head=T)

# # create a list of SNPs per region:
# regions<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_regions_425_cojo_supervised_plus_unsupervised.tsv.gz"))
# regions$min<-gsub("^[0-9]{1,2}_","",regions$updated_region)
# regions$max<-gsub("^[0-9]*_","",regions$min)
# regions$min<-gsub("_[0-9]*$","",regions$min)
# regions$chr<-gsub("_.*","",regions$updated_region)

# regions<-as.data.frame(regions)

# # create a list of SNPs per region:
# regions<-regions[!duplicated(regions),]

# table(regions$class)
# # new old 
# # 149 276

# # renames those regions:
# regions<-regions[,c("chr","updated_region","min","max")]
# regions<-as.data.frame(regions[!duplicated(regions),])

# regions<-regions[order(regions$min,decreasing="F"),]
# dim(regions)
# # [1] 425   4

# regions$min<-as.numeric(regions$min)
# regions$max<-as.numeric(regions$max)
# regions$chr<-as.numeric(regions$chr)

# # renames those regions:
# regions<-regions[,c("chr","updated_region","min","max")]
# regions<-as.data.frame(regions[!duplicated(regions),])

# regions<-regions[order(regions$min,decreasing="F"),]
# dim(regions)
# # [1] 425   4

# array<-c("illumina370","humanomniexpress","affymetrix500","humanomni1","quad610","illuminaexome","affymetrix6","humancoreexome","gsa")
# pheno<-c("ibd","cd","uc")

# for (i in 1:length(array)) {

#   print(array[i])
  
#   for (j in 1:length(pheno)) {
    
#     print(pheno[j])
    
#     # create .incl file:
    
#     ph<-fread(paste(path_gwas,"post_imputation/2022/",array[i],"/phenotype_data/",array[i],"_all_studies_merged_eur_all_phenotype_",pheno[j],sep=""),head=T)
#     ph<-as.data.frame(ph)
    
#     colnames(ph)[3]<-"ph"
#     ph<-ph[which(!is.na(ph$ph)),]
#     write.table(ph$IID,paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/",array[i],"/",pheno[j],"/individuals_dose.subset.incl",sep=""),
#                 col.names=F,row.names=F,quote=F,sep=" ")
  
      
#     for (chr in c(1:22)){
      
#       if (chr==23){
#         chr<-"X"
#       }
      
#       print(chr)
      
#       regions_tmp<-regions[which(regions$chr==chr),]
      
#       file<-paste(path_gwas,"post_imputation/2022/analysis/regenie/",array[i],"/",pheno[j],"/chr",chr,"_",array[i],"_eur_all_step2_",pheno[j],"_eur_sex_PCs_firthse_",pheno[j],".regenie.gz",sep="")

#       if (file.exists(file)) {

#         for (ii in 1:nrow(regions_tmp)) {

#           # list of variants in bgen file, that pass info and maf freq thresholds for ibd, cd and uc
#           z<-fread(paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/",array[i],"/",array[i],"_chr",regions_tmp$updated_region[ii],".dose.subset_list_variants_tokeep_in_analysis.variants.gz",sep=""),head=T)
          

#           if(nrow(z)>2) {
#             z<-z[,c("ID","#CHROM","POS","REF","ALT")]
#             colnames(z)<-c("rsid","chromosome","position","allele1","allele2")
            
#             fwrite(z,paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/",array[i],"/",array[i],"_chr",regions_tmp$updated_region[ii],".dose.subset.z",sep=""),
#                   col.names=T,row.names=F,quote=F,sep=" ")
#           }

#           # create .data file:
          
#           dat<-matrix(ncol=9,nrow=1)
#           dat<-as.data.frame(dat)
#           colnames(dat)<-c("z","bgen","bgi","bcor","ld","n_samples","bdose","sample","incl")
          
#           # variants in file:
#           dat$z<-paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/",array[i],"/",array[i],"_chr",regions_tmp$updated_region[ii],".dose.subset.z",sep="")

#           # individuals to keep to estimate the ld
#           dat$incl<-paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/",array[i],"/",pheno[j],"/individuals_dose.subset.incl",sep="")
          
#           dat$n_samples<-nrow(ph)
        
#           # input genotyping data

#           dat$bgen<-paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/",array[i],"/",array[i],"_chr",regions_tmp$updated_region[ii],".dose.subset_list_variants_tokeep_in_analysis.bgen",sep="")
#           dat$bgi<-paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/",array[i],"/",array[i],"_chr",regions_tmp$updated_region[ii],".dose.subset_list_variants_tokeep_in_analysis.bgen.bgi",sep="")
#           dat$sample<-paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/",array[i],"/",array[i],"_chr",regions_tmp$updated_region[ii],".dose.subset_list_variants_tokeep_in_analysis.sample",sep="")
          
#           # pre-defined outputs
          
#           dat$ld<-paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/",array[i],"/",pheno[j],"/",array[i],"_",pheno[j],"_chr",regions_tmp$updated_region[ii],".dose.subset.ld",sep="")
#           dat$bcor<-paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/",array[i],"/",pheno[j],"/",array[i],"_",pheno[j],"_chr",regions_tmp$updated_region[ii],".dose.subset.bcor",sep="")
#           dat$bdose<-paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/",array[i],"/",pheno[j],"/",array[i],"_",pheno[j],"_chr",regions_tmp$updated_region[ii],".dose.subset.bdose",sep="")
          
#           fwrite(dat,paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/",array[i],"/",pheno[j],"/",array[i],"_",pheno[j],"_chr",regions_tmp$updated_region[ii],".dose.subset.data",sep=""),
#                  col.names=T,row.names=F,quote=F,sep=";")
          
        
#         } 

        
#       } else {
#         print("NO REGENIE FILE")
#       }
#     }
    
#     rm(ph)
#   }
# }

# q("no")



# ####################################
# # 3.- Create bcor, bdose, ld files #
# ####################################

# ####################################
# # 3.1 - create bdose files


# path_gwas=/path/to/ibdgwas/IIBDGC/
# array=(illumina370 humanomniexpress affymetrix500 humanomni1 quad610 illuminaexome affymetrix6 humancoreexome gsa)
# pheno=(ibd cd uc)


# # export the env variables defined above
# mapfile -t region < <(cat /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/list_regions)
# mapfile -t chromosomes < <(cat /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/list_chromosomes)
# MEM=15000
# MEM1=15000

# # mapfile -t region < <(cat /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/list_regions | grep 6_29357657_37435948)
# # mapfile -t chromosomes < <(cat /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/list_chromosomes | grep 6_29357657_37435948)

# # MEM=95000
# # MEM1=95000

# # mapfile -t region < <(cat /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/list_regions | grep -E '1_91246790_93849934|19_1977822_4575286|19_48193018_50925886|4_37533101_40806352|9_3580755_6322802')
# # mapfile -t chromosomes < <(cat /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/list_chromosomes | grep -E '1_91246790_93849934|19_1977822_4575286|19_48193018_50925886|4_37533101_40806352|9_3580755_6322802')
# # MEM=25000
# # MEM1=25000

# array=(gsa)
# pheno=(uc)

# for i in ${array[@]}
# do
# for ph in ${pheno[@]}
# do
# for ii in {0..11}  
# do 
# bsub -J"bdose" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM:tmp=$MEM1] span[hosts=1]" -G ibdgwas -n 2 \
# -e ${path_gwas}post_imputation/log/${i}_chr${region[${ii}]}_${ph}_bdose_stderr \
# -o ${path_gwas}post_imputation/log/${i}_chr${region[${ii}]}_${ph}_bdose_stdout \
# "/path/to/software/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 \
# --in-files ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${ph}/${i}_${ph}_chr${region[${ii}]}.dose.subset.data \
# --read-only-bgen --write-text"
# done
# done
# done

# for i in ${array[@]}; do echo ${i} && for ph in ${pheno[@]}; do echo ${ph} && ls -la ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${ph}/${i}_${ph}_*.dose.subset.ld | wc -l; done; done

# pheno=(ibd cd uc)
# for i in ${array[@]}; do echo ${i} && for ph in ${pheno[@]}; do echo ${ph} && ls -lh ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${ph}/${i}_${ph}_*.dose.subset.ld | grep '[0-9]K'; done; done


# ## remove intermediate files files

# array=(illumina370 humanomniexpress affymetrix500 humanomni1 quad610 illuminaexome affymetrix6 humancoreexome gsa)

# for i in ${array[@]}
# do
# rm ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/*.dose.subset_list_variants_tokeep_in_analysis.bgen
# done

# for i in ${array[@]}
# do
# rm ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/*.dose.subset_list_variants_tokeep_in_analysis.bgen.bgi
# done

# for i in ${array[@]}
# do
# rm ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/*.dose.subset_list_variants_tokeep_in_analysis.sample
# done

# lfs quota -hg ibdgwas /path/to/project

# ## compress all existing files

# array=(illumina370 humanomniexpress affymetrix500 humanomni1 quad610 illuminaexome affymetrix6 humancoreexome gsa)
# array=(gsa)
# MEM=300
# MEM1=300

# for i in ${array[@]}
# do
# for ph in ${pheno[@]}
# do
# for ii in {0..419}  
# do 
# bsub -J"gzip" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM:tmp=$MEM1] span[hosts=1]" -G ibdgwas -n 2 \
# -e ${path_gwas}post_imputation/log/${i}_chr${region[${ii}]}_${ph}_bdose_ld_compress_stderr \
# -o ${path_gwas}post_imputation/log/${i}_chr${region[${ii}]}_${ph}_bdose_ld_compress_stdout \
# "gzip ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${ph}/${i}_${ph}_chr${region[${ii}]}.dose.subset.ld"
# done
# done
# done

# for i in ${array[@]}; do echo ${i} && for ph in ${pheno[@]}; do echo ${ph} && ls -lh ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${ph}/${i}_${ph}_*.dose.subset.ld.gz | wc -l; done; done
# for i in ${array[@]}; do echo ${i} && for ph in ${pheno[@]}; do echo ${ph} && ls -lh ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${ph}/${i}_${ph}_*.dose.subset.ld.gz | grep '[0-9]K'; done; done


# for i in ${array[@]}
# do
# rm ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/*subset_list_variants_tokeep_in_analysis.variants.gz
# done

# # for i in ${array[@]}
# # do
# # for ph in ${pheno[@]}
# # do
# # for ii in {0..419} 
# # do 
# # bsub -J"bdose" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM:tmp=$MEM1] span[hosts=1]" -G ibdgwas -n 2 \
# # -e ${path_gwas}post_imputation/log/${i}_chr${region[${ii}]}_${ph}_bcorr_stderr \
# # -o ${path_gwas}post_imputation/log/${i}_chr${region[${ii}]}_${ph}_bcorr_stdout \
# # "/path/to/software/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 \
# # --in-files ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${ph}/${i}_${ph}_chr${region[${ii}]}.dose.subset.data \
# # --read-only-bgen --write-bcor"
# # done
# # done
# # done

# # for i in ${array[@]}; do echo ${i} && for ph in ${pheno[@]}; do echo ${ph} && ls -la  ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${ph}/${i}_${ph}_*.dose.subset.bcor | wc -l; done; done

# #  UC
# # 16_78268798_81321307|16_81383307_83333852|17_38829907_43388282|19_16540599_20148843|4_120204800_123139785|5_171651212_174467164|6_19280778_22613885|7_25666136_28690783|7_3901126_7438156|7_97875685_101755911

# # ibd
# # 10_60240494_64101958|13_39138280_43328617|16_78268798_81321307|16_81383307_83333852|17_38829907_43388282|19_16540599_20148843|19_53703491_55371596|6_19280778_22613885|7_3901126_7438156|7_78257_1655687|7_97875685_101755911

# # cd
# # 13_39138280_43328617|16_78268798_81321307|16_81383307_83333852|17_38829907_43388282|17_66321594_68925063|19_16540599_20148843|4_37533101_40806352|5_171651212_174467164|7_25666136_28690783|7_3901126_7438156|7_97875685_101755911|9_132635004_134593510

# for i in ${array[@]}
# do
# rm ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/*subset_list_variants_tokeep_in_analysis.variants.gz
# done

# # remove ld files per array after sharing them wiht the IIBDGC
# for ph in ${pheno[@]}
# do
# for i in ${array[@]}
# do
# rm ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${ph}/${i}_${ph}_*.dose.subset.ld.gz
# done
# done






############################################################################################################################################################################################################
############################################################################################################################################################################################################

############################################################################################################################################################################################################
############################################################################################################################################################################################################


# SECOND APPROACH - COMBINE ALL GENOTYPING ARRAY DATA AND ESTIMATE EUR TIER 1 LD MATRICES:

# singularity exec iibdgc_postprocess_10_singularity.sif

path_gwas="/path/to/ibdgwas/IIBDGC/"

############################################################################
# 1.1 - Define the list of samples and variants to extract from the bgen files:

MEM=8000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
path_gwas<-"/path/to/ibdgwas/IIBDGC/"

array<-c("illumina370","humanomniexpress","affymetrix500","humanomni1","quad610","illuminaexome","affymetrix6","humancoreexome","gsa")
pheno<-c("ibd","cd","uc")


# create a list of SNPs per region:
regions<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_regions_420_cojo_supervised_plus_unsupervised.tsv.gz"))
regions$min<-gsub("^[0-9]{1,2}_","",regions$updated_region)
regions$max<-gsub("^[0-9]*_","",regions$min)
regions$min<-gsub("_[0-9]*$","",regions$min)
regions$chr<-gsub("_.*","",regions$updated_region)

regions<-as.data.frame(regions)

# create a list of SNPs per region:
regions<-regions[!duplicated(regions),]

table(regions$class)
# new old 
# 144 276

# renames those regions:
regions<-regions[,c("chr","updated_region","min","max")]
regions<-as.data.frame(regions[!duplicated(regions),])

regions<-regions[order(regions$min,decreasing="F"),]
dim(regions)
# [1] 420   4

regions$min<-as.numeric(regions$min)
regions$max<-as.numeric(regions$max)
regions$chr<-as.numeric(regions$chr)

# renames those regions:
regions<-regions[,c("chr","updated_region","min","max")]
regions<-as.data.frame(regions[!duplicated(regions),])

regions<-regions[order(regions$min,decreasing="F"),]
dim(regions)
# [1] 425   4

# create a file with the chrs and regions to export later as env variables to submit the job:
write.table(regions[,"updated_region",drop=F],paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/list_regions",sep=""),
      col.names=F,row.names=F,quote=F,sep="\t")

write.table(regions[,"chr",drop=F],paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/list_chromosomes",sep=""),
      col.names=F,row.names=F,quote=F,sep="\t")

# restrict this to the list of regions we need to update/create"
regions<-regions[which(regions$updated_region %in% c("12_8422106_10257536","14_106197560_107197560","16_69127624_70127624","18_21986192_22986192","2_229716690_231482974","4_41896369_42896369","9_93754276_94754276")),]

# define the list of samples per phenotype and genotyping array:


array<-c("illumina370","humanomniexpress","affymetrix500","humanomni1","quad610","illuminaexome","affymetrix6","humancoreexome","gsa")
pheno<-c("ibd","cd","uc")

for (i in 1:length(array)) {

  print(array[i])
  
 
  for (j in 1:length(pheno)) {
    
    print(pheno[j])
    
    if (array[i]=="humanomniexpress" & pheno[j]=="cd") {
      print(paste(array[i],"no",pheno[j],"analyisis"))
    } else {
          # create .incl file:
    
          ph<-fread(paste(path_gwas,"post_imputation/2022/",array[i],"/phenotype_data/",array[i],"_all_studies_merged_eur_all_phenotype_",pheno[j],sep=""),head=T)
          ph<-as.data.frame(ph)
          
          colnames(ph)[3]<-"ph"
          ph<-ph[which(!is.na(ph$ph)),]
          write.table(ph$IID,paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/",array[i],"/",pheno[j],"/individuals_dose.subset.incl",sep=""),
                col.names=F,row.names=F,quote=F,sep=" ")

    }

  
  }

}


# define the list of variants per phenotype - variants that contribute to tier_1 eur analyses:

for (i in 1:length(pheno)) {
  
  print(pheno[i])
   
  # for (chr in c(1:22)){
  for (chr in c(12,14,16,18,2,4,9)){
    
    tmp1<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/",pheno[i],"/chr",chr,"_",pheno[i],"_meta_noGC_PCs_firthse_formatted_A2_effect_and_alternative_with_per_sample_rate_and_avgA2freq_with_rsid_dbsnp154.txt.gz",sep=""),head=T)
    tmp1$chr<-chr

    tmp1<-as.data.frame(tmp1[order(tmp1$Position_b38,decreasing=F),])
    regions_tmp<-as.data.frame(regions[which(regions$chr==chr),])

    for (ii in 1:nrow(regions_tmp)) {

      var_tokeep<-tmp1[which(tmp1$Position_b38>=regions_tmp$min[ii] & tmp1$Position_b38<=regions_tmp$max[ii]),]

      write.table(var_tokeep[,"MarkerName"],paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/chr",regions_tmp$updated_region[ii],"_list_variants_tokeep_in_analysis_",pheno[i],sep=""),
      col.names=F,row.names=F,quote=F,sep="\t")
      rm(var_tokeep)

    }

    rm(regions_tmp)
  
  }
  rm(tmp1)
  
}

q("no")

##################



##################################################################
# 1.2 - subset bgen files

path_gwas=/path/to/ibdgwas/IIBDGC/
array=(illumina370 humanomniexpress affymetrix500 humanomni1 quad610 illuminaexome affymetrix6 humancoreexome gsa)

# export the env variables defined above
mapfile -t region < <(cat /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/list_regions)
mapfile -t chromosomes < <(cat /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/list_chromosomes)
pheno=(ibd cd uc)

# # export the env variables defined above - subset just the 6 new ones, plus the one that needs to be updated:
mapfile -t region < <(cat /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/list_regions | grep -E "12_8422106_10257536|14_106197560_107197560|16_69127624_70127624|18_21986192_22986192|2_229716690_231482974|4_41896369_42896369|9_93754276_94754276")
mapfile -t chromosomes < <(cat /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/list_regions | grep -E "12_8422106_10257536|14_106197560_107197560|16_69127624_70127624|18_21986192_22986192|2_229716690_231482974|4_41896369_42896369|9_93754276_94754276" | sed 's/\_.*//')

for ii in {0..6}
do 
echo ${region[${ii}]} && echo ${chromosomes[${ii}]}
done
done

# double check n samples per array and pheno
for i in ${array[@]}
do
echo ${i} && for ph in ${pheno[@]}
do
echo ${ph} && wc -l  ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${ph}/individuals_dose.subset.incl
done 
done


##########
# illumina370
# ibd
# 3263 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/illumina370/ibd/individuals_dose.subset.incl
# cd
# 1961 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/illumina370/cd/individuals_dose.subset.incl
# uc
# 2026 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/illumina370/uc/individuals_dose.subset.incl
# humanomniexpress
# ibd
# 1235 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/humanomniexpress/ibd/individuals_dose.subset.incl
# cd
# 589 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/humanomniexpress/cd/individuals_dose.subset.incl - TO REMOVE, NO CASES
# uc
# 1235 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/humanomniexpress/uc/individuals_dose.subset.incl
# affymetrix500
# ibd
# 4591 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/affymetrix500/ibd/individuals_dose.subset.incl
# cd
# 4542 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/affymetrix500/cd/individuals_dose.subset.incl
# uc
# 2934 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/affymetrix500/uc/individuals_dose.subset.incl
# humanomni1
# ibd
# 2034 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/humanomni1/ibd/individuals_dose.subset.incl
# cd
# 1659 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/humanomni1/cd/individuals_dose.subset.incl
# uc
# 1642 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/humanomni1/uc/individuals_dose.subset.incl
# quad610
# ibd
# 3364 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/quad610/ibd/individuals_dose.subset.incl
# cd
# 2583 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/quad610/cd/individuals_dose.subset.incl
# uc
# 2234 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/quad610/uc/individuals_dose.subset.incl
# illuminaexome
# ibd
# 1333 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/illuminaexome/ibd/individuals_dose.subset.incl
# cd
# 870 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/illuminaexome/cd/individuals_dose.subset.incl
# uc
# 745 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/illuminaexome/uc/individuals_dose.subset.incl
# affymetrix6
# ibd
# 8180 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/affymetrix6/ibd/individuals_dose.subset.incl
# cd
# 2818 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/affymetrix6/cd/individuals_dose.subset.incl
# uc
# 8121 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/affymetrix6/uc/individuals_dose.subset.incl
# humancoreexome
# ibd
# 20563 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/humancoreexome/ibd/individuals_dose.subset.incl
# cd
# 15802 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/humancoreexome/cd/individuals_dose.subset.incl
# uc
# 14866 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/humancoreexome/uc/individuals_dose.subset.incl
# gsa
# ibd
# 57998 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/gsa/ibd/individuals_dose.subset.incl
# cd
# 42599 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/gsa/cd/individuals_dose.subset.incl
# uc
# 30577 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/gsa/uc/individuals_dose.subset.incl

# file with only controls
rm /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/humanomniexpress/cd/individuals_dose.subset.incl

# double check the number of variants per region:
for ii in {0..6}
do 
echo ${region[${ii}]} && for i in ${array[@]}
do
echo ${i} && for ph in ${pheno[@]}
do
echo ${ph} && wc -l ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/chr${region[${ii}]}_list_variants_tokeep_in_analysis_${ph} 
done
done

##################

MEM=150
for i in ${array[@]}
do
for ii in {0..6} 
do 
for ph in ${pheno[@]}
do
bsub -J"bgen_subset" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-o ${path_gwas}post_imputation/2022/log/${i}_chr${region[${ii}]}_bgen_subset_variants_${ph}_stdout \
-e ${path_gwas}post_imputation/2022/log/${i}_chr${region[${ii}]}_bgen_subset_variants_${ph}_stderr \
"/path/to/software/username/qctool_v2.2.0/./qctool \
-g ${path_gwas}post_imputation/2022/${i}/imputed_data/${i}_chr${chromosomes[${ii}]}.dose.subset.bgen \
-s ${path_gwas}post_imputation/2022/${i}/imputed_data/${i}_chr${chromosomes[${ii}]}.dose.subset.sample \
-incl-rsids ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/chr${region[${ii}]}_list_variants_tokeep_in_analysis_${ph} \
-incl-samples ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${ph}/individuals_dose.subset.incl \
-og ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${i}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf && \
bgzip -f ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${i}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf && \
bcftools index ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${i}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz"
done
done
done

for i in ${array[@]}
do
echo ${i} && for ii in {0..6} 
do 
echo ${i} && for ph in ${pheno[@]}
do
echo ${ph} && tail -50 ${path_gwas}post_imputation/2022/log/${i}_chr${region[${ii}]}_bgen_subset_variants_${ph}_stdout  | grep "Successfully"
done
done
done



array=(illumina370 humanomniexpress affymetrix500 humanomni1 quad610 illuminaexome affymetrix6 humancoreexome gsa)

for i in ${array[@]}
do
echo ${i} && for ph in ${pheno[@]}
do
# echo ${ph} && ls -lh ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${i}_*.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz | grep ' 0 May'
echo ${ph} && ls -la ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${i}_*.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz | wc -l
done
done


##################################################################
# 1.3 - merge vcf files

array=(illumina370 humanomniexpress affymetrix500 humanomni1 quad610 illuminaexome affymetrix6 humancoreexome gsa)
MEM=150

pheno=(ibd uc)

for ii in {0..6} 
do 
for ph in ${pheno[@]}
do
bsub -J"reg_ibd_test" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -n 2 \
-o ${path_gwas}post_imputation/2022/log/chr${region[${ii}]}_merge_variants_${ph}_stdout \
-e ${path_gwas}post_imputation/2022/log/chr${region[${ii}]}_merge_variants_${ph}_stderr \
"bcftools merge \
${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${array[0]}/${array[0]}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz \
${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${array[1]}/${array[1]}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz \
${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${array[2]}/${array[2]}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz \
${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${array[3]}/${array[3]}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz \
${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${array[4]}/${array[4]}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz \
${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${array[5]}/${array[5]}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz \
${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${array[6]}/${array[6]}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz \
${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${array[7]}/${array[7]}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz \
${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${array[8]}/${array[8]}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz \
-Oz -o ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz \
-m none \
--threads 2"
done
done


pheno=(cd)

for ii in {0..6}
do 
for ph in ${pheno[@]}
do
bsub -J"merg_reg" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -n 2 \
-o ${path_gwas}post_imputation/2022/log/chr${region[${ii}]}_merge_variants_${ph}_stdout \
-e ${path_gwas}post_imputation/2022/log/chr${region[${ii}]}_merge_variants_${ph}_stderr \
"bcftools merge \
${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${array[0]}/${array[0]}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz \
${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${array[2]}/${array[2]}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz \
${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${array[3]}/${array[3]}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz \
${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${array[4]}/${array[4]}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz \
${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${array[5]}/${array[5]}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz \
${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${array[6]}/${array[6]}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz \
${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${array[7]}/${array[7]}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz \
${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${array[8]}/${array[8]}_chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz \
-Oz -o ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz \
-m none \
--threads 2"
done
done


for ii in {0..6} 
do 
echo ${region[${ii}]} && for ph in ${pheno[@]}
do
echo ${ph} && tail -50 ${path_gwas}post_imputation/2022/log/chr${region[${ii}]}_merge_variants_${ph}_stdout | grep "Successfully"
done
done



pheno=(ibd cd uc)
for ph in ${pheno[@]}
do
echo ${ph} && ls -la ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/chr*.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz | wc -l
done




# list number of samples

for ii in {0..6}
do
echo ${region[${ii}]} && for ph in ${pheno[@]}
do
echo ${ph} && bcftools query -l ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz | wc -l
done
done
# 102561 - ph = IBD OK
# 64380 - ph = UC OK
# 72834 - ph = CD OK

# list number of variants
for ii in {0..6}
do
echo ${region[${ii}]} && for ph in ${pheno[@]}
do
bcftools query -f '%CHROM %POS %REF %ALT\n' ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz | wc -l \
&& wc -l ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/chr${region[${ii}]}_list_variants_tokeep_in_analysis_${ph} | cut -d ' ' -f1
done
done


##################################################################
# 1.5 - convert into bgen

MEM=100
for ii in {0..6} 
do 
for ph in ${pheno[@]}
do
bsub -J"reg_ibd_test" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-o ${path_gwas}post_imputation/2022/log/chr${region[${ii}]}_vcf_to_bgen_${ph}_stdout \
-e ${path_gwas}post_imputation/2022/log/chr${region[${ii}]}_vcf_to_bgen_${ph}_stderr \
"/path/to/software/username/qctool_v2.2.0/./qctool \
-g ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_tmp.vcf.gz \
-filetype vcf -vcf-genotype-field GP \
-os ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}.dose.subset.sample \
-og ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}.dose.subset.bgen"
done
done


for ii in {0..6} 
do 
echo ${region[${ii}]} && for ph in ${pheno[@]}
do echo ${ph} && tail -50 ${path_gwas}post_imputation/2022/log/chr${region[${ii}]}_vcf_to_bgen_${ph}_stdout | grep -E "Successfully|Exit"
done
done

# remove intermed files:
rm ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/*_tmp.vcf.gz 

array=(illumina370 humanomniexpress affymetrix500 humanomni1 quad610 illuminaexome affymetrix6 humancoreexome gsa)
for i in ${array[@]}
do rm ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/${i}/${i}_*dose.subset_list_variants_tokeep_in_analysis_*_tmp.vcf.gz
done


############################
# 1.6 - index files

# Index bgen files

MEM=100
for ii in {0..6} 
do 
for ph in ${pheno[@]}
do
bsub -J"bgi_subset" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-o ${path_gwas}post_imputation/2022/log/chr${region[${ii}]}_bgen_index_${ph}_stdout \
-e ${path_gwas}post_imputation/2022/log/chr${region[${ii}]}_bgen_index_${ph}_stderr \
"bgenix -g ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}.dose.subset.bgen -index "
done
done

for ii in {0..6} 
do 
for ph in ${pheno[@]}
do
echo ${ph} && ls -la ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}.dose.subset.bgen.bgi
done
done

########################################################
# 1.4 - get the list of SNPs in the bgen file

for ph in ${pheno[@]}
do
mkdir ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/${ph}/
done

for ph in ${pheno[@]}
do
mv ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/*.dose.subset_list_variants_tokeep_in_analysis_${ph}.dose.subset.bgen* \
${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/${ph}/
done

for ph in ${pheno[@]}
do
mv ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/*.dose.subset_list_variants_tokeep_in_analysis_${ph}.dose.subset.sample \
${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/${ph}/
done

for ph in ${pheno[@]}
do
mv ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/*_list_variants_tokeep_in_analysis_${ph} \
${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/${ph}/
done

for ph in ${pheno[@]}
do
ls -la ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/${ph}/*.dose.subset_list_variants_tokeep_in_analysis_${ph}.dose.subset.bgen | wc -l
done

########################

MEM=100

# same variants between different phenotypes 
for ii in {0..6} 
do 
for ph in ${pheno[@]}
do
bsub -J"bgi" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-o ${path_gwas}post_imputation/2022/log/chr${region[${ii}]}_bgen_list_variants_${ph}_stdout \
-e ${path_gwas}post_imputation/2022/log/chr${region[${ii}]}_bgen_list_variants_${ph}_stderr \
"bgenix -g \
${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/${ph}/chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}.dose.subset.bgen \
-vcf | cut -f1-5 | gzip > \
${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/${ph}/chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}.variants.gz"
done
done

# CONTINUE HERE!

for ph in ${pheno[@]}
do
echo ${ph} && for ii in {0..6} 
do 
echo ${ii} && tail -50 ${path_gwas}post_imputation/2022/log/chr${region[${ii}]}_bgen_list_variants_${ph}_stdout | grep -E "Successfully|Exit"
done
done


## move each one to their pheno folder:

#################################
# 2.- create z and master files #
#################################

# The master file is a semicolon-separated text file and contains no space. It contains the following mandatory 
# column names and one dataset per line.
# z column contains the names of Z files (input)
# bgen column contains the names of BGEN files (input)
# bgi column contains the names of BGI files (input)
# bcor column contains the names of BCOR files when using --write-bcor (output)
# ld column contains the names of LD files when using --bcor-to-text (output)
# n_samples column contains the number of samples to include in any processing
# bdose column contains the names of optional BDOSE files when using --write-bdose (optional output)
# sample column contains the names of optional sample files when using --write-bdose and --bdose-version 1.1 (optional input)
# incl column contains the names of optional sample inclusion files (optional input)

# A master file with two datasets for writing SNP correlations to a BCOR or text file could look as follows.
# z;bgen;bgi;bcor;ld;n_samples
# dataset1.z;dataset1.bgen;dataset1.bgen.bgi;dataset1.bcor;;5363
# dataset2.z;dataset2.bgen;dataset2.bgen.bgi;;dataset2.ld;5363

# incl
# The dataset.incl file is a text file to specify inclusion of samples in any processing. It contains one sample ID per line.


# The dataset.z file is a space-delimited text file and contains meta information about the SNPs one SNP per line. 
# It contains the mandatory column names in the following order.
# A dataset.z file with three SNPs could look as follows.
# rsid chromosome position allele1 allele2
# rs1 10 1 T C
# rs2 10 1 A G
# rs3 10 1 G A

MEM=4000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

library(data.table)
path_gwas<-"/path/to/ibdgwas/IIBDGC/"

all<-fread(paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_2_new_and_old_index_gwas_variants_conditional_on_known_index_gwas_fm_2023_with_iibdgc_tier1_tier2_eas_sas_summary_stats_plus_unsupervised_conditional.tsv",sep=""),head=T)


# create a list of SNPs per region:
regions<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_regions_425_cojo_supervised_plus_unsupervised.tsv.gz"))
regions$min<-gsub("^[0-9]{1,2}_","",regions$updated_region)
regions$max<-gsub("^[0-9]*_","",regions$min)
regions$min<-gsub("_[0-9]*$","",regions$min)
regions$chr<-gsub("_.*","",regions$updated_region)

regions<-as.data.frame(regions)

# create a list of SNPs per region:
regions<-regions[!duplicated(regions),]

table(regions$class)
# new old 
# 149 276

# renames those regions:
regions<-regions[,c("chr","updated_region","min","max")]
regions<-as.data.frame(regions[!duplicated(regions),])

regions<-regions[order(regions$min,decreasing="F"),]
dim(regions)
# [1] 425   4

regions$min<-as.numeric(regions$min)
regions$max<-as.numeric(regions$max)
regions$chr<-as.numeric(regions$chr)

# renames those regions:
regions<-regions[,c("chr","updated_region","min","max")]
regions<-as.data.frame(regions[!duplicated(regions),])

regions<-regions[order(regions$min,decreasing="F"),]
dim(regions)
# [1] 425   4

pheno<-c("ibd","cd","uc")

for (j in 1:length(pheno)) {
    
  print(pheno[j])
    

  for (chr in c(1:22)){
      
    if (chr==23){
      chr<-"X"
    }
      
    print(chr)
      
    regions_tmp<-regions[which(regions$chr==chr),]
      
    for (ii in 1:nrow(regions_tmp)) {

      # list of variants in bgen file, that pass info and maf freq thresholds for ibd, cd and uc
      z<-fread(paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/",pheno[j],"/chr",regions_tmp$updated_region[ii],".dose.subset_list_variants_tokeep_in_analysis_",pheno[j],".variants.gz",sep=""),head=T)
          
      if(nrow(z)>2) {
        z<-z[,c("ID","#CHROM","POS","REF","ALT")]
        colnames(z)<-c("rsid","chromosome","position","allele1","allele2")
            
        fwrite(z,paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/",pheno[j],"/chr",regions_tmp$updated_region[ii],".dose.subset_list_variants_tokeep_in_analysis_",pheno[j],".z",sep=""),
              col.names=T,row.names=F,quote=F,sep=" ")
        rm(z)
      }

      # create .data file:
      dat<-matrix(ncol=7,nrow=1)
      dat<-as.data.frame(dat)
      colnames(dat)<-c("z","bgen","bgi","bcor","ld","n_samples","bdose")
          
      # variants in file - file created above:
      dat$z<-paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/",pheno[j],"/chr",regions_tmp$updated_region[ii],".dose.subset_list_variants_tokeep_in_analysis_",pheno[j],".z",sep="")

      # number of samples included:
      samples<-fread(paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/",pheno[j],"/chr",regions_tmp$updated_region[ii],".dose.subset_list_variants_tokeep_in_analysis_",pheno[j],".dose.subset.sample",sep=""),head=T)
      dat$n_samples<-nrow(samples)-1
      rm(samples)

      # input genotyping data
      dat$bgen<-paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/",pheno[j],"/chr",regions_tmp$updated_region[ii],".dose.subset_list_variants_tokeep_in_analysis_",pheno[j],".dose.subset.bgen",sep="")
      dat$bgi<-paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/",pheno[j],"/chr",regions_tmp$updated_region[ii],".dose.subset_list_variants_tokeep_in_analysis_",pheno[j],".dose.subset.bgen.bgi",sep="")
      # dat$sample<-paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/",pheno[j],"/chr",regions_tmp$updated_region[ii],".dose.subset_list_variants_tokeep_in_analysis_",pheno[j],".dose.subset.sample",sep="")
          
      # pre-defined outputs
          
      dat$ld<-paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/",pheno[j],"/chr",regions_tmp$updated_region[ii],".dose.subset_list_variants_tokeep_in_analysis_",pheno[j],"_all_arrays_combined.ld",sep="")
      dat$bcor<-paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/",pheno[j],"/chr",regions_tmp$updated_region[ii],".dose.subset_list_variants_tokeep_in_analysis_",pheno[j],"_all_arrays_combined.bcor",sep="")
      dat$bdose<-paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/",pheno[j],"/chr",regions_tmp$updated_region[ii],".dose.subset_list_variants_tokeep_in_analysis_",pheno[j],"_all_arrays_combined.bdose",sep="")
          
      fwrite(dat,paste(path_gwas,"post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/",pheno[j],"/chr",regions_tmp$updated_region[ii],".dose.subset_list_variants_tokeep_in_analysis_",pheno[j],"_all_arrays_combined.data",sep=""),
          col.names=T,row.names=F,quote=F,sep=";") 
        
    } 
    
  }
}

q("no")

# example of final data file:
# z;bgen;bgi;bcor;ld;n_samples;bdose;sample
# /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/ibd/chr1_780044_3066589.dose.subset_list_variants_tokeep_in_analysis_ibd.z;
# /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/ibd/chr1_780044_3066589.dose.subset_list_variants_tokeep_in_analysis_ibd.dose.subset.bgen
# /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/ibd/chr1_780044_3066589.dose.subset_list_variants_tokeep_in_analysis_ibd.dose.subset.bgen.bgi
# /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/ibd/chr1_780044_3066589.dose.subset_list_variants_tokeep_in_analysis_ibd_all_arrays_combined.bcor
# /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/ibd/chr1_780044_3066589.dose.subset_list_variants_tokeep_in_analysis_ibd_all_arrays_combined.ld
# 102562
# /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/ibd/chr1_780044_3066589.dose.subset_list_variants_tokeep_in_analysis_ibd_all_arrays_combined.bdose


###############################################################
# estimate LD

MEM=40000
MEM1=40000

# # export the env variables defined above
# mapfile -t region < <(cat /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/list_regions)
# mapfile -t chromosomes < <(cat /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/list_chromosomes)

# # export the env variables defined above - subset just the 6 new ones, plus the one that needs to be updated:
mapfile -t region < <(cat /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/list_regions | grep -E "12_8422106_10257536|14_106197560_107197560|16_69127624_70127624|18_21986192_22986192|2_229716690_231482974|4_41896369_42896369|9_93754276_94754276")
mapfile -t chromosomes < <(cat /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/finemapping/ld_matrices/list_regions | grep -E "12_8422106_10257536|14_106197560_107197560|16_69127624_70127624|18_21986192_22986192|2_229716690_231482974|4_41896369_42896369|9_93754276_94754276" | sed 's/\_.*//')

for ii in {0..6}
do 
echo ${region[${ii}]} && echo ${chromosomes[${ii}]}
done
done


pheno=(ibd cd uc)
for ph in ${pheno[@]}
do
for ii in {0..6} 
do 
bsub -J"bdose" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM:tmp=$MEM1] span[hosts=1]" -G ibdgwas -q week \
-e ${path_gwas}post_imputation/log/${ph}_chr${region[${ii}]}_${ph}_bdose_stderr \
-o ${path_gwas}post_imputation/log/${ph}_chr${region[${ii}]}_${ph}_bdose_stdout \
"/path/to/software/username/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 \
--in-files ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/${ph}/chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_all_arrays_combined.data \
--read-only-bgen --write-text"
done
done


for ph in ${pheno[@]}
do
for ii in {0..6} 
do echo ${ii} && tail -50 ${path_gwas}post_imputation/log/${ph}_chr${region[${ii}]}_${ph}_bdose_stdout | grep -E "Exit|Successfully"
done
done


for ph in ${pheno[@]}
do
ls -la ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/${ph}/chr*.dose.subset_list_variants_tokeep_in_analysis_${ph}_all_arrays_combined.ld | wc -l 
done



## compress all existing files


MEM=200
for ph in ${pheno[@]}
do
for ii in {0..6}  
do 
bsub -J"gzip" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/log/${i}_chr${region[${ii}]}_${ph}_bdose_ld_compress_stderr \
-o ${path_gwas}post_imputation/log/${i}_chr${region[${ii}]}_${ph}_bdose_ld_compress_stdout \
"gzip ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/${ph}/chr${region[${ii}]}.dose.subset_list_variants_tokeep_in_analysis_${ph}_all_arrays_combined.ld"
done
done

for ph in ${pheno[@]}
do
ls -la ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/${ph}/chr*.dose.subset_list_variants_tokeep_in_analysis_${ph}_all_arrays_combined.ld.gz | wc -l 
done
# 420
# 420
# 420

# list of variants
for ph in ${pheno[@]}
do
ls -la ${path_gwas}post_imputation/2022/analysis/finemapping/ld_matrices/all_arrays_combined/${ph}/chr*.dose.subset_list_variants_tokeep_in_analysis_${ph}.z | wc -l
done
# 420
# 420
# 420
