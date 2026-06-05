#!/bin/bash
# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# Author: https://orcid.org/0000-0002-7256-9752
#

chr=$1

path=/path/to/ibdgwas/IIBDGC/
mem=8500

wrapper_script=~/git/IIBDGC_GWAS/scripts/other/estimate_CpG_wrapper.sh

r_script=~/git/IIBDGC_GWAS/scripts/other/estimate_CpG.R

bsub -R"select[mem>${mem}] rusage[mem=${mem}]" -M${mem} -q week -o ${path}scripts/logs/${chr}_CpG_stdout -e ${path}scripts/logs/${chr}_CpG_stderr bash ${wrapper_script} ${r_script} ${chr} 

exit $?