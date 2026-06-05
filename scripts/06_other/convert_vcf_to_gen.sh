#!/bin/bash
# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# Author: https://orcid.org/0000-0002-7256-9752
#

# extract the variants form imputation

study=$1
imputation_panel=$2

path_gwas=/path/to/ibdgwas/IIBDGC/
pathway_script=${path_gwas}imputed/scripts/

wrapper_script=convert_vcf_to_gen_wrapper.sh

chr=23

bsub -n2 -J"vcfgen[1-${chr}]" -R"select[mem>500] rusage[mem=500] span[hosts=1]" -M500 -o ${path_gwas}imputed/${study}/logs/vcfgen_stdout_${study}${imputation_panel} -e ${path_gwas}imputed/${study}/logs/vcfgen_stderr_${study}${imputation_panel} bash ${pathway_script}${wrapper_script} ${study} ${imputation_panel} 

exit $?
