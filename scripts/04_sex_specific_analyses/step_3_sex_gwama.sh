#!/bin/bash
# Author: Talin Haritunians
# Institution: F. Widjaja Inflammatory Bowel Disease Institute
#

## sex-specific meta-analysis

## prepare regenie outputs for GWAMA input
## see regenie2GWAMA.R

## run GWAMA
# per {phenotype}


for chr in {1..22}; do
	
	GWAMA --filelist gwama_{pheno}_${chr}.in --sex --indel_alleles --output iibdgc_sex_{pheno}_${chr}

done


# gwama_{pheno}_${chr}.in = per chr per phenotype list of male and female GWAMA input files created with regenie2GWAMA.R