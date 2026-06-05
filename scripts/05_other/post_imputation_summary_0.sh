#!/bin/bash
# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# Author: https://orcid.org/0000-0002-7256-9752
#

# create summary plots and tables from imputation outputs

study=$1
chr=$2

MEM=2000

path_gwas=/path/to/ibdgwas/IIBDGC/
pathway_script=${path_gwas}scripts/

wrapper_script=${path_gwas}scripts/post_imputation_summary_0_wrapper.sh
r_script=${path_gwas}scripts/post_imputation_summary_0.R

bsub -J"postImp_${study}_${chr}" -q normal -G "ibdgwas" -m "modern_hardware" -M${MEM} -R"select[mem>${MEM}] rusage[mem=${MEM}] span[hosts=1]" -o ${pathway_script}logs/post_imputation_summary_0_stdout_${study}_${chr} -e ${pathway_script}logs/post_imputation_summary_stderr_0_${study}_${chr} bash ${wrapper_script} ${study} ${chr} ${r_script}  

exit $?