#!/bin/bash
# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# Author: https://orcid.org/0000-0002-7256-9752
#

# create fam files from final call report style files created by L. Philip Schumm <pschumm@uchicago.edu> 

study=$1
queue=$2
array=$3

path_gwas=/path/to/ibdgwas/IIBDGC/
pathway_script=${path_gwas}scripts/
  
wrapper_script=create_lgen_from_fcr_wrapper.sh

bsub -J"fcrlgen${study}" -q ${queue} -G "team152" -M500 -R"select[mem>500] rusage[mem=500] span[hosts=1]" -n2 -o ${pathway_script}logs/fcrlgen_stdout_${study} -e ${pathway_script}logs/fcrlgen_stderr_${study} bash ${pathway_script}${wrapper_script} ${study} ${path_gwas} ${array}

exit $?