# Author: Talin Haritunians
# Institution: F. Widjaja Inflammatory Bowel Disease Institute
#

### regenie2GWAMA.R


library(data.table)
library(dplyr)
library(tidyverse)

chr <- commandArgs(trailingOnly = TRUE)

file_name <- chr[1]

	data <- fread(file_name, header=TRUE)

	data2 <- data %>%
  		dplyr::rename(MARKERNAME = ID, EA = ALLELE1, NEA = ALLELE0, N = N, EAF = A1FREQ)

	data2$range <- data2$SE * 1.96

	data2$L95 <- ifelse(data2$range < 0, 
                    data2$BETA + data2$range, 
                    data2$BETA - data2$range)

	data2$U95 <- ifelse(data2$range < 0, 
                    data2$BETA - data2$range, 
                    data2$BETA + data2$range)

	data2$OR = exp(data2$BETA)

	data2$OR_95L = exp(data2$L95)

	data2$OR_95U = exp(data2$U95)

	outgwama <- data2 %>%
    	select(MARKERNAME, EA, NEA, N, EAF, OR, OR_95L, OR_95U)
	
		
	output_name <- paste0(chr, ".{sex}_{pheno}_GWAMA.txt")
  		write.table(outgwama, output_name, quote = FALSE, sep = "\t", row.names = FALSE)




##
## run with gwama_prep.sh

#!/bin/bash

#for chr in {1..22}; do

#	Rscript regenie2GWAMA.R step2_{array}_chr${i}.{sex}.{pheno}_eur.regenie
	
#done
