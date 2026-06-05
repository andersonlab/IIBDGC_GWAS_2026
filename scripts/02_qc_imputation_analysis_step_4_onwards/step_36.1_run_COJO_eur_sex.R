# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# double check we are not leaving out any region by using COJO to determine GWAS significant hits:

# give bfile as input

# --bfile test 
# Note: --bgen, --pfile, --bpfile, --mbfile, --mbgen, --mpfile and --mbpfile are currently only supported in the GRM and fastGWA analyses. 
# More functions will be available after rewriting some of the legacy codes. All the QC flags (e.g. --keep, --extract, --maf) in GCTA 1.92.4 are 
# currently applicable to these two formats.


# singularity exec iibdgc_postprocess_10_singularity.sif

# /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/metaanalysis/sex_specific_analyses/

# tar -xvzf cd_sexhet_cojo.tar.gz
# tar -xvzf uc_sexhet_cojo.tar.gz
# tar -xvzf ibd_sexhet_cojo.tar.gz
# mv IBD2clean_for_COJO/ ibd_sexhet_cojo
# for chr in {1..22} 
# do mv chr${chr}_sex_IBD2_gwama_metal_withHCE_linked_out_InfoMafFilt.txt chr${chr}_sex_IBD_gwama_metal_withHCE_linked_out_InfoMafFilt.txt
# done
# tar -xvzf aver_geno_rate_sex_male.tar.gz




path_gwas="/path/to/ibdgwas/IIBDGC/"
pheno=(ibd cd uc)
MEM=18000

for ph in ${pheno[@]}
do
bsub -J"ph" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}post_imputation/2022/log/reformat_${ph}_sex_specific_summary_stats_for_cojo_stderr \
-o ${path_gwas}post_imputation/2022/log/reformat_${ph}_sex_specific_summary_stats_for_cojo_stdout \
"Rscript ~/git/IIBDGC_GWAS/scripts/other/format_sex_specific_summary_stats_for_cojo.R ${ph} > \
${path_gwas}post_imputation/2022/log/format_sex_specific_summary_stats_for_cojo_${ph}.Rout"
done

for ph in ${pheno[@]}
do
echo ${ph} && tail -50 ${path_gwas}post_imputation/2022/log/reformat_${ph}_sex_specific_summary_stats_for_cojo_stdout | grep "Successfully"
done

sex_specific=(female male female_male)

for ph in ${pheno[@]}
do
echo ${ph} && for sex in ${sex_specific[@]}
do
echo ${sex} && wc -l ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/allchr_${ph}_meta_eur_tier_1_${sex}_sex_specific_gwama_cojo_format.ma
done
done

# ibd
# female
# 11179401 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/allchr_ibd_meta_eur_tier_1_female_sex_specific_gwama_cojo_format.ma
# male
# 11179432 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/allchr_ibd_meta_eur_tier_1_male_sex_specific_gwama_cojo_format.ma
# female_male
# 11179485 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/allchr_ibd_meta_eur_tier_1_female_male_sex_specific_gwama_cojo_format.ma
# cd
# female
# 11608952 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/allchr_cd_meta_eur_tier_1_female_sex_specific_gwama_cojo_format.ma
# male
# 11609760 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/allchr_cd_meta_eur_tier_1_male_sex_specific_gwama_cojo_format.ma
# female_male
# 11601477 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/allchr_cd_meta_eur_tier_1_female_male_sex_specific_gwama_cojo_format.ma
# uc
# female
# 11084830 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/allchr_uc_meta_eur_tier_1_female_sex_specific_gwama_cojo_format.ma
# male
# 11084592 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/allchr_uc_meta_eur_tier_1_male_sex_specific_gwama_cojo_format.ma
# female_male
# 11077678 /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/allchr_uc_meta_eur_tier_1_female_male_sex_specific_gwama_cojo_format.ma


#########################

path_gwas=/path/to/ibdgwas/IIBDGC/
pheno=(ibd cd uc)

MEM=45000

sex_specific=(female male female_male)

for ph in ${pheno[@]}
do
for chr in {1..22}
do
for sex in ${sex_specific[@]}
do
bsub -J"cojo_i" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-o ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_${sex}_no_conditioning_stdout \
-e ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_${sex}_no_conditioning_stderr \
"/path/to/software/username/gcta-1.94.1-linux-kernel-3-x86_64/./gcta64 \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_subset_included_in_${ph}_analysis \
--chr ${chr} \
--cojo-file ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/allchr_${ph}_meta_eur_tier_1_${sex}_sex_specific_gwama_cojo_format.ma \
--cojo-slct \
--cojo-p 5e-8 \
--out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/${chr}_${ph}_independent_index_variants_eur_tier_1_${sex}_sex_specific_gwama_ld_allarrays_no_conditioning"
done
done
done


for sex in ${sex_specific[@]}
do
echo ${sex} && for ph in ${pheno[@]}
do
echo ${ph} && for chr in {1..22} 
do
echo ${chr} && tail -30 ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_${sex}_no_conditioning_stdout | grep -E "Successfully|Exited"
done
done
done


for sex in ${sex_specific[@]}
do
echo ${sex} && for ph in ${pheno[@]}
do
echo ${ph} && wc -l ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/*_${ph}_independent_index_variants_eur_tier_1_${sex}_sex_specific_gwama_ld_allarrays_no_conditioning.jma.cojo
done
done



for sex in ${sex_specific[@]}
do
echo ${sex} && for ph in ${pheno[@]}
do
echo ${ph} && ls -la ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/*_${ph}_independent_index_variants_eur_tier_1_${sex}_sex_specific_gwama_ld_allarrays_no_conditioning.log | wc -l
done
done

for sex in ${sex_specific[@]}
do
echo ${sex} && for ph in ${pheno[@]}
do echo ${ph} && \
echo ${ph} && ls -la ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/*_${ph}_independent_index_variants_eur_tier_1_${sex}_sex_specific_gwama_ld_allarrays_no_conditioning.* | wc -l
done
done


for sex in ${sex_specific[@]}
do
echo ${sex} && for ph in ${pheno[@]}
do echo ${ph} && \
head -2 ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/${chr}_${ph}_independent_index_variants_eur_tier_1_${sex}_sex_specific_gwama_ld_allarrays_no_conditioning.cma.cojo
done
done


# female
# ibd
# Chr     SNP     bp      refA    freq    b       se      p       n       freq_geno       bC      bC_se   pC
# 1       chr1:1017466:G:A        1017466 A       0.003347        -0.0929621      0.140024        0.506754        45368.9 0.00220799   -0.0934353      0.140023        0.504591
# cd
# Chr     SNP     bp      refA    freq    b       se      p       n       freq_geno       bC      bC_se   pC
# 1       chr1:784474:CAG:C       784474  C       0.0075  0.0393303       0.151675        0.795398        17800.8 0.00211311  0.0393303        0.151671        0.795393
# uc
# Chr     SNP     bp      refA    freq    b       se      p       n       freq_geno       bC      bC_se   pC
# 1       chr1:981210:G:A 981210  A       0.176207        -0.00250814     0.0263898       0.924281        29572.1 0.16604 -0.00246254  0.0263894       0.925653
# male
# ibd
# Chr     SNP     bp      refA    freq    b       se      p       n       freq_geno       bC      bC_se   pC
# 1       chr1:1017466:G:A        1017466 A       0.003046        0.147847        0.154745        0.339362        41210.8 0.00220799   0.147847        0.154745        0.339361
# cd
# Chr     SNP     bp      refA    freq    b       se      p       n       freq_geno       bC      bC_se   pC
# 1       chr1:784474:CAG:C       784474  C       0.006905        0.495723        0.175338        0.00469495      14901.4 0.00211311   0.495723        0.175379        0.00470468
# uc
# Chr     SNP     bp      refA    freq    b       se      p       n       freq_geno       bC      bC_se   pC
# 1       chr1:981210:G:A 981210  A       0.175698        0.018331        0.0265449       0.48984 29196.5 0.16604 0.018454    0.0265447        0.486927
# female_male
# ibd
# Chr     SNP     bp      refA    freq    b       se      p       n       freq_geno       bC      bC_se   pC
# 1       chr1:1017466:G:A        1017466 A       0.003202        0.0154471       0.103828        0.881731        86640.4 0.00220799   0.0176469       0.103828        0.865039
# cd
# Chr     SNP     bp      refA    freq    b       se      p       n       freq_geno       bC      bC_se   pC
# 1       chr1:784474:CAG:C       784474  C       0.007224        0.234676        0.114712        0.0407768       32753.5 0.00211311   0.234676        0.114717        0.0407866
# uc
# Chr     SNP     bp      refA    freq    b       se      p       n       freq_geno       bC      bC_se   pC
# 1       chr1:981210:G:A 981210  A       0.175954        0.00785011      0.0187142       0.67487 58749.5 0.16604 0.00847372  0.0187141        0.650693


for sex in ${sex_specific[@]}
do
echo ${sex} && for ph in ${pheno[@]}
do echo ${ph} && \
grep 'chr1:1537160:T:G' ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/allchr_${ph}_meta_eur_tier_1_${sex}_sex_specific_gwama_cojo_format.ma
done
done

# female
# ibd
# chr1:1537160:T:G        G       T       0.338873        0.0363406080911177      0.0184207339232454      0.048534        52911
# cd
# chr1:1537160:T:G        G       T       0.338995        0.0478382827025523      0.0229617511375942      0.037235        36265
# uc
# chr1:1537160:T:G        G       T       0.335559        0.029373348152794       0.02329711190999        0.207355        32429
# male
# ibd
# chr1:1537160:T:G        G       T       0.338802        0.0313533076579957      0.0191782469188033      0.102081        49650
# cd
# chr1:1537160:T:G        G       T       0.339132        0.0354781415378224      0.0249689747318914      0.155334        31790
# uc
# chr1:1537160:T:G        G       T       0.336271        0.0375443087723367      0.023441370689747       0.109234        31951
# female_male
# ibd
# chr1:1537160:T:G        G       T       0.338839        0.0339472172996575      0.013285084616551       0.010629        102561
# cd
# chr1:1537160:T:G        G       T       0.339059        0.0421749987970228      0.0169015648813816      0.012604        68055
# uc
# chr1:1537160:T:G        G       T       0.335912        0.0334338089677812      0.0165243871118441      0.04306 64380


for sex in ${sex_specific[@]}
do
echo ${sex} && for ph in ${pheno[@]}
do echo ${ph} && \
head -1 ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/${chr}_${ph}_independent_index_variants_eur_tier_1_${sex}_sex_specific_gwama_ld_allarrays_no_conditioning.cma.cojo && \
grep 'chr1:1537160:T:G' ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/${chr}_${ph}_independent_index_variants_eur_tier_1_${sex}_sex_specific_gwama_ld_allarrays_no_conditioning.cma.cojo
done
done


# female
# ibd
# Chr     SNP     bp      refA    freq    b       se      p       n       freq_geno       bC      bC_se   pC
# 1       chr1:1537160:T:G        1537160 G       0.338873        0.0363406       0.0184207       0.0485172       39029.3 0.329723     0.0355209       0.0184214       0.0538251
# cd
# Chr     SNP     bp      refA    freq    b       se      p       n       freq_geno       bC      bC_se   pC
# 1       chr1:1537160:T:G        1537160 G       0.338995        0.0478383       0.0229618       0.0372157       25797   0.320821     0.0478383       0.0229632       0.037228
# uc
# Chr     SNP     bp      refA    freq    b       se      p       n       freq_geno       bC      bC_se   pC
# 1       chr1:1537160:T:G        1537160 G       0.335559        0.0293733       0.0232971       0.207376        24702.6 0.318319     0.0283122       0.0232974       0.224269
# male
# ibd
# Chr     SNP     bp      refA    freq    b       se      p       n       freq_geno       bC      bC_se   pC
# 1       chr1:1537160:T:G        1537160 G       0.338802        0.0313533       0.0191782       0.102083        36369   0.329723     0.0313533       0.0191787       0.102091
# cd
# Chr     SNP     bp      refA    freq    b       se      p       n       freq_geno       bC      bC_se   pC
# 1       chr1:1537160:T:G        1537160 G       0.339132        0.0354781       0.024969        0.155349        22492.1 0.320821     0.0354781       0.0249695       0.155358
# uc
# Chr     SNP     bp      refA    freq    b       se      p       n       freq_geno       bC      bC_se   pC
# 1       chr1:1537160:T:G        1537160 G       0.336271        0.0375443       0.0234414       0.109238        24292.1 0.318319     0.0368671       0.0234421       0.115792
# female_male
# ibd
# Chr     SNP     bp      refA    freq    b       se      p       n       freq_geno       bC      bC_se   pC
# 1       chr1:1537160:T:G        1537160 G       0.338839        0.0339472       0.0132851       0.01061 75390.3 0.329723    0.0332576        0.0132856       0.0123045
# cd
# Chr     SNP     bp      refA    freq    b       se      p       n       freq_geno       bC      bC_se   pC
# 1       chr1:1537160:T:G        1537160 G       0.339059        0.042175        0.0169016       0.012584        48284.2 0.320821     0.042175        0.0169025       0.0125888
# uc
# Chr     SNP     bp      refA    freq    b       se      p       n       freq_geno       bC      bC_se   pC
# 1       chr1:1537160:T:G        1537160 G       0.335912        0.0334338       0.0165244       0.0430421       48973.6 0.318319     0.0327857       0.0165249       0.0472542


ph=ibd
head -1 ${path_gwas}post_imputation/2022/analysis/metaanalysis/sex_specific_analyses/${ph}_sexhet_cojo/chr${chr}_sex_IBD_gwama_metal_withHCE_linked_out_InfoMafFilt.txt
grep 'chr1:1537160:T:G' ${path_gwas}post_imputation/2022/analysis/metaanalysis/sex_specific_analyses/${ph}_sexhet_cojo/chr${chr}_sex_IBD_gwama_metal_withHCE_linked_out_InfoMafFilt.txt
rs_number       reference_allele        other_allele    eaf     OR      OR_se   OR_95L  OR_95U  z       p-value _-log10_p-value      q_statistic     q_p-value       i2      n_studies       n_samples       effects male_eaf        male_OR male_OR_se  male_OR_95L      male_OR_95U     male_z  male_p-value    male_n_studies  male_n_samples  female_eaf      female_OR       female_OR_se female_OR_95L   female_OR_95U   female_z        female_p-value  female_n_studies        female_n_samples        gender_differentiated_p-value        gender_heterogeneity_p-value    Finfo   Fmaf    Minfo   Mmaf    Allele1_F_Metal Allele2_F_Metal      Effect_F_Metal  StdErr_F_Metal  P-value_F_Metal Direction_F_Metal       HetISq_F_Metal  HetChiSq_F_Metal        HetDf_F_Metal        HetPVal_F_Metal Allele1_M_Metal Allele2_M_Metal Effect_M_Metal  StdErr_M_Metal  P-value_M_Metal Direction_M_Metal    HetISq_M_Metal  HetChiSq_M_Metal        HetDf_M_Metal   HetPVal_M_Metal avg_rate_F      avg_N_F avg_rate_M      avg_N_M
chr1:1537160:T:G        G       T       0.338839        1.03453 0.013566        1.00794 1.061821        2.555288        0.010629     1.9735  11.089669       0.851865        0       18      102561  ++++++-+++++--++++      0.338802        1.03185 0.019422     0.993783        1.071374        1.634837        0.102081        9       49650   0.338873        1.037009        0.018762     1.000236        1.075133        1.97281 0.048534        9       52911   0.03754 0.85121 NA      NA      NA      NA  tg       0.0363  0.0184  0.04852 ++++++-++       0       6.257   8       0.6184  t       g       0.0314  0.0192  0.1021  +++--++++    0       4.797   8       0.779   1       5879    1       5516.66666666667

ph=cd
head -1 ${path_gwas}post_imputation/2022/analysis/metaanalysis/sex_specific_analyses/${ph}_sexhet_cojo/chr${chr}_sex_CD_gwama_metal_withHCE_linked_out_InfoMafFilt.txt
grep 'chr1:1537160:T:G' ${path_gwas}post_imputation/2022/analysis/metaanalysis/sex_specific_analyses/${ph}_sexhet_cojo/chr${chr}_sex_CD_gwama_metal_withHCE_linked_out_InfoMafFilt.txt
rs_number       reference_allele        other_allele    eaf     OR      OR_se   OR_95L  OR_95U  z       p-value _-log10_p-value      q_statistic     q_p-value       i2      n_studies       n_samples       effects male_eaf        male_OR male_OR_se  male_OR_95L      male_OR_95U     male_z  male_p-value    male_n_studies  male_n_samples  female_eaf      female_OR       female_OR_se female_OR_95L   female_OR_95U   female_z        female_p-value  female_n_studies        female_n_samples        gender_differentiated_p-value        gender_heterogeneity_p-value    Finfo   Fmaf    Minfo   Mmaf    Allele1_F_Metal Allele2_F_Metal      Effect_F_Metal  StdErr_F_Metal  P-value_F_Metal Direction_F_Metal       HetISq_F_Metal  HetChiSq_F_Metal        HetDf_F_Metal        HetPVal_F_Metal Allele1_M_Metal Allele2_M_Metal Effect_M_Metal  StdErr_M_Metal  P-value_M_Metal Direction_M_Metal    HetISq_M_Metal  HetChiSq_M_Metal        HetDf_M_Metal   HetPVal_M_Metal avg_rate_F      avg_N_F avg_rate_M      avg_N_M
chr1:1537160:T:G        G       T       0.339059        1.043077        0.017341        1.009089        1.078209        2.495331     0.012604        1.899486        10.306747       0.503024        0       12      68055   +++-+++-+-++    0.339132    1.036115 0.025248        0.986629        1.088083        1.420889        0.155334        6       31790   0.338995        1.049001     0.023553        1.002837        1.097289        2.08339 0.037235        6       36265   0.041597        0.715595    NA       NA      NA      NA      t       g       0.0478  0.023   0.03722 +++-++  15.9    5.944   5       0.3117  t       g   0.0355   0.025   0.1553  +-+-++  0       4.23    5       0.5167


ph=uc
head -1 ${path_gwas}post_imputation/2022/analysis/metaanalysis/sex_specific_analyses/${ph}_sexhet_cojo/chr${chr}_sex_UC_gwama_metal_withHCE_linked_out_InfoMafFilt.txt
grep 'chr1:1537160:T:G' ${path_gwas}post_imputation/2022/analysis/metaanalysis/sex_specific_analyses/${ph}_sexhet_cojo/chr${chr}_sex_UC_gwama_metal_withHCE_linked_out_InfoMafFilt.txt
rs_number       reference_allele        other_allele    eaf     OR      OR_se   OR_95L  OR_95U  z       p-value _-log10_p-value      q_statistic     q_p-value       i2      n_studies       n_samples       effects male_eaf        male_OR male_OR_se  male_OR_95L      male_OR_95U     male_z  male_p-value    male_n_studies  male_n_samples  female_eaf      female_OR       female_OR_se female_OR_95L   female_OR_95U   female_z        female_p-value  female_n_studies        female_n_samples        gender_differentiated_p-value        gender_heterogeneity_p-value    Finfo   Fmaf    Minfo   Mmaf    Allele1_F_Metal Allele2_F_Metal      Effect_F_Metal  StdErr_F_Metal  P-value_F_Metal Direction_F_Metal       HetISq_F_Metal  HetChiSq_F_Metal        HetDf_F_Metal        HetPVal_F_Metal Allele1_M_Metal Allele2_M_Metal Effect_M_Metal  StdErr_M_Metal  P-value_M_Metal Direction_M_Metal    HetISq_M_Metal  HetChiSq_M_Metal        HetDf_M_Metal   HetPVal_M_Metal avg_rate_F      avg_N_F avg_rate_M      avg_N_M
chr1:1537160:T:G        G       T       0.335912        1.033999        0.016812        1.001047        1.068036        2.023301     0.04306 1.365924        12.145113       0.791272        0       18      64380   -+--+++++-++--++++      0.336271    1.038258 0.023787        0.991634        1.087073        1.601626        0.109234        9       31951   0.335559        1.029809     0.023452        0.983843        1.077924        1.260815        0.207355        9       32429   0.125252        0.804738     NA      NA      NA      NA      t       g       0.0294  0.0233  0.2074  -+--+++++       0       5.889   8       0.6597       t       g       0.0375  0.0234  0.1092  -++--++++       0       6.195   8       0.6254  1       5879    1       5516.66666666667


#########################

path_gwas=/path/to/ibdgwas/IIBDGC/
pheno=(ibd cd uc)

MEM=45000

sex_specific=(female male female_male)

for ph in ${pheno[@]}
do
for chr in {1..22}
do
for sex in ${sex_specific[@]}
do
bsub -J"cojo_i" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-o ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_${sex}_no_conditioning_stdout \
-e ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_${sex}_no_conditioning_stderr \
"/path/to/software/username/gcta-1.94.1-linux-kernel-3-x86_64/./gcta64 \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_subset_included_in_${ph}_analysis \
--chr ${chr} \
--cojo-cond ${path_gwas}post_imputation/2022/analysis/conditional_analysis/final_independent_model/final_independent_model_final_join_model_no_subset_final.snplist \
--cojo-file ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/allchr_${ph}_meta_eur_tier_1_${sex}_sex_specific_gwama_cojo_format.ma \
--out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/${chr}_${ph}_independent_index_variants_eur_tier_1_${sex}_sex_specific_gwama_ld_allarrays__conditioning_on_final_independent_model_0.1_other_variants"
done
done
done

# CONTINUE HERE

for ph in ${pheno[@]}
do
echo ${ph} && for chr in {1..22}
do
echo ${chr} && for sex in ${sex_specific[@]}
do
echo  ${sex} && tail -50  ${path_gwas}post_imputation/2022/log/cojo_eur_${ph}_${chr}_${sex}_no_conditioning_stdout | grep "Successfully"
done
done
done

cd  /path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/

head -1 9_cd_independent_index_variants_eur_tier_1_female_male_sex_specific_gwama_ld_allarrays__conditioning_on_final_independent_model_0.1_other_variants.cma.cojo && grep chr9:114814020:G:A 9_*_independent_index_variants_eur_tier_1_*_sex_specific_gwama_ld_allarrays__conditioning_on_final_independent_model_0.1_other_variants.cma.cojo
# Chr     SNP     bp      refA    freq    b       se      p       n       freq_geno       bC      bC_se   pC
# 9_cd_independent_index_variants_eur_tier_1_female_male_sex_specific_gwama_ld_allarrays__conditioning_on_final_independent_model_0.1_other_variants.cma.cojo:9   chr9:114814020:G:A      114814020       A       0.136586        0.1332  0.0194803     8.04764e-12     69030.1 0.137105        0.0357181       0.0194868       0.0668116
# 9_cd_independent_index_variants_eur_tier_1_female_sex_specific_gwama_ld_allarrays__conditioning_on_final_independent_model_0.1_other_variants.cma.cojo:9        chr9:114814020:G:A      114814020       A       0.137216        0.119102     0.0265435        7.22127e-06     36524.5 0.137105        0.00595368      0.0265504       0.82257
# 9_cd_independent_index_variants_eur_tier_1_male_sex_specific_gwama_ld_allarrays__conditioning_on_final_independent_model_0.1_other_variants.cma.cojo:9  chr9:114814020:G:A      114814020       A       0.135869        0.149658        0.0286789     1.80465e-07     32520.7 0.137105        0.0703416       0.0286904       0.0142165
# 9_ibd_independent_index_variants_eur_tier_1_female_male_sex_specific_gwama_ld_allarrays__conditioning_on_final_independent_model_0.1_other_variants.cma.cojo:9  chr9:114814020:G:A      114814020       A       0.136689        0.107835     0.0156549        5.64761e-12     103034  0.139974        0.0314928       0.0156584       0.0443003
# 9_ibd_independent_index_variants_eur_tier_1_female_sex_specific_gwama_ld_allarrays__conditioning_on_final_independent_model_0.1_other_variants.cma.cojo:9       chr9:114814020:G:A      114814020       A       0.13694 0.0897219       0.0217657     3.75314e-05     52980   0.139974        0.00949809      0.021769        0.662609
# 9_ibd_independent_index_variants_eur_tier_1_male_sex_specific_gwama_ld_allarrays__conditioning_on_final_independent_model_0.1_other_variants.cma.cojo:9 chr9:114814020:G:A      114814020       A       0.136422        0.127247        0.0225329     1.63099e-08     50067.8 0.139974        0.0550793       0.0225399       0.0145399
# 9_uc_independent_index_variants_eur_tier_1_female_male_sex_specific_gwama_ld_allarrays__conditioning_on_final_independent_model_0.1_other_variants.cma.cojo:9   chr9:114814020:G:A      114814020       A       0.131478        0.0806727    0.0198874        4.98215e-05     66038.7 0.134025        0.0295972       0.0198898       0.136735
# 9_uc_independent_index_variants_eur_tier_1_female_sex_specific_gwama_ld_allarrays__conditioning_on_final_independent_model_0.1_other_variants.cma.cojo:9        chr9:114814020:G:A      114814020       A       0.131936        0.055377     0.0280169        0.0480914       33249.9 0.134025        0.0121074       0.0280181       0.665649
# 9_uc_independent_index_variants_eur_tier_1_male_sex_specific_gwama_ld_allarrays__conditioning_on_final_independent_model_0.1_other_variants.cma.cojo:9  chr9:114814020:G:A      114814020       A       0.131012        0.106363        0.0282345     0.000165147     32815.5 0.134025        0.047353        0.0282402       0.0935823


# combine in one tar file:

MEM=200
for sex in ${sex_specific[@]}
do
bsub -J"tar_files" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/log/${i}_${ph}_ld_tar_stderr \
-o ${path_gwas}post_imputation/log/${i}_${ph}_ld_tar_stdout \
"tar -cvzf ${path_gwas}to_upload_to_iibdgc_globus/${sex}_sex_specific_summary_stats_conditioned_619_indep_variants.tar.gz ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/*_independent_index_variants_eur_tier_1_${sex}_sex_specific_gwama_ld_allarrays__conditioning_on_final_independent_model_0.1_other_variants.cma.cojo"
done

for sex in ${sex_specific[@]}
do
tail -50 ${path_gwas}post_imputation/log/${i}_${ph}_ld_tar_stdout | grep "Successfully"
done

/path/to/ibdgwas/IIBDGC/to_upload_to_iibdgc_globus/female_male_sex_specific_summary_stats_conditioned_619_indep_variants.tar.gz



MEM=200

bsub -J"tar_files" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas \
-e ${path_gwas}post_imputation/log/${i}_${ph}_ld_tar_stderr \
-o ${path_gwas}post_imputation/log/${i}_${ph}_ld_tar_stdout \
"tar -cvzf ${path_gwas}to_upload_to_iibdgc_globus/allarrays_allchr_list_tier_1_sex_specific_conditional_variants_for_ld.tar.gz  ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr*_list_tier_1_sex_specific_conditional_variants_for_ld.ld"




########################################################################################
# Define independet sex-specific signals, and integrate those with our main analyses

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=15000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

path_gwas<-"/path/to/ibdgwas/IIBDGC/"

### is there any significant signal after conditioning on the 653?
pheno<-c("ibd","cd","uc")
sex_specific<-c("female","male","female_male")

for (sex in sex_specific) {

  for (ph in pheno) {

    for (chr in 1:22) {

      tmp<-read.table(paste0(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/",chr,"_",ph,"_independent_index_variants_eur_tier_1_",sex,"_sex_specific_gwama_ld_allarrays__conditioning_on_final_independent_model_0.1_other_variants.cma.cojo"),head=T)
      tmp<-tmp[which(tmp$p<5E-8 | tmp$pC<5E-8),]

      if(chr==1) {
        dat_tmp<-tmp
      } else {
        dat_tmp<-rbind(dat_tmp,tmp)
      }
      rm(tmp)
    }

    dat_tmp<-dat_tmp[,c("SNP","bC","bC_se","pC")]
    colnames(dat_tmp)[2:ncol(dat_tmp)]<-paste(colnames(dat_tmp)[2:ncol(dat_tmp)],ph,sep="_")
    
    print(dim(dat_tmp))

    if(ph=="ibd") {
      dat<-dat_tmp
    } else {
      dat<-merge(dat,dat_tmp,by="SNP",all=T)
    }
    rm(dat_tmp)
  } 

    dat$sex<-sex

    if (sex=="female") {
      dat_all<-dat
    } else {
      dat_all<-rbind(dat,dat_all)
    }

}

ids<-dat_all[which( (dat_all$sex %in% c("male","female")) & (dat_all$pC_ibd<5E-8 | dat_all$pC_uc<5E-8 | dat_all$pC_cd<5E-8) ),"SNP"]
ids<-ids[!duplicated(ids)]

ids[!grepl("chr6:",ids)]
#  [1] "chr10:62608500:T:C"   "chr10:62608958:T:C"   "chr10:62611228:C:G"  
#  [4] "chr10:62614600:G:A"   "chr10:62617743:G:A"   "chr16:50220349:CTT:C"
#  [7] "chr16:50278049:C:CT"  "chr10:62594503:C:T"   "chr10:62594693:G:A"  
# [10] "chr10:62597623:G:A"   "chr10:62605073:G:A"   "chr10:62609990:G:A"  
# [13] "chr10:72606376:T:G"   "chr10:72748205:G:A"   "chr10:72828909:T:G"  
# [16] "chr16:50377862:T:C"   "chr16:50712383:A:C"   "chr19:46636823:G:C"  
# [19] "chr19:46640496:C:T"   "chr19:46642401:T:C"   "chr2:240622406:G:C"  
# [22] "chr2:240627454:G:A"   "chr2:240630974:A:G"   "chr2:240631319:A:G"  
# [25] "chr2:240633278:C:T"   "chr2:240634984:G:A"   "chr2:240639120:C:G"  
# [28] "chr2:240639565:A:G"   "chr2:240639691:G:T"   "chr2:240639752:T:C"  
# [31] "chr2:240641033:C:T"   "chr2:240643822:G:A"   "chr2:240647279:C:T"  
# [34] "chr2:240649568:A:G"   "chr2:240653740:A:G"   "chr2:240657967:A:G"  

########################################################################################################################################################################################################################

# singularity exec iibdgc_postprocess_10_singularity.sif

MEM=8000
bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R

path_gwas<-"/path/to/ibdgwas/IIBDGC/"


library(data.table)
rm(list=ls())

path_gwas<-"/path/to/ibdgwas/IIBDGC/"



all<-fread(paste(path_gwas,"post_imputation/2022/analysis/final_tables/list_tier_2_unsupervised_conditional_frequency_with_old_data_with_new_class_signal_no_fm_signals_with_multiancestry_signals_final_join_model_with_exonic_variants_join_conditinal_model.tsv.gz",sep=""))
dim(all)
# [1] 653 163


## collect all independent signals:

# sex<-c("female","male","female_male")
sex<-c("female","male")
pheno<-c("ibd","cd","uc")

for (i in c(1:length(sex))) {

  print(sex[i])
  
  for (ph in pheno) {


    print(ph)

      for (chr in c(1:22)) {

        file_tmp<-paste0(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/",chr,"_",ph,"_independent_index_variants_eur_tier_1_",sex[i],"_sex_specific_gwama_ld_allarrays_no_conditioning.jma.cojo")
       
        if (file.exists(file_tmp)) {
          tmp<-fread(file_tmp)
        }

        if(chr==1) {
          dat<-tmp
        } else {
          dat<-rbind(dat,tmp)
        }
      }
    print(nrow(dat))

    dat<-dat[,c("SNP","b","se","p","bJ","bJ_se","pJ")]
    colnames(dat)[2:ncol(dat)]<-c("BETA","SE","Pvalue","BETA_cond","SE_cond","Pvalue_cond")
    colnames(dat)[2:ncol(dat)]<-paste0(colnames(dat)[2:ncol(dat)],"_",sex[i],"_",ph)

    if (ph=="ibd" && sex[i]=="female") {
      dat_final<-dat
    } else {
      dat_final<-merge(dat,dat_final,by="SNP",all=T)
    }

  }

}

dim(dat_final)
# [1] 343  37
dat_final<-dat_final[!duplicated(dat_final),]
dim(dat_final)
# [1] 318  37

dat_final[duplicated(dat_final$SNP),]
# 0


# find out LD between variants:

list_variants<-c(dat_final$SNP,all$MarkerName)
length(list_variants)
# [1] 971

list_variants<-list_variants[!duplicated(list_variants)]
length(list_variants)
# [1] 917

write.table(list_variants,paste(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_1_sex_specific_conditional_variants_for_ld.tsv",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

###############

# note this will creat hardcalls, and those will be used for ld1 estimation, all variants good imputation, it should1 not be biased...

MEM=2000
path_gwas="/path/to/ibdgwas/IIBDGC/"


for chr in {1..22} 
do
bsub -J"cojo_i" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-o ${path_gwas}post_imputation/2022/log/cojo_eur_sex_specific_${ph}_${chr}_ld11_stdout \
-e ${path_gwas}post_imputation/2022/log/cojo_eur_sex_specific_${ph}_${chr}_ld11_stderr \
"/path/to/software/username/./plink2 \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_subset_included_in_ibd_analysis \
--extract ${path_gwas}post_imputation/2022/analysis/metaanalysis/forward_regression/list_tier_1_sex_specific_conditional_variants_for_ld.tsv \
--make-bed --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_list_tier_1_sex_specific_conditional_variants_for_ld"
done

for chr in {1..22} 
do
echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/cojo_eur_sex_specific_${ph}_${chr}_ld11_stdout | grep "Successfully"
done

wc -l ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr*_list_tier_1_sex_specific_conditional_variants_for_ld.bim
# 917 total

for chr in {1..22} 
do
bsub -J"cojo_i" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-o ${path_gwas}post_imputation/2022/log/cojo_eur_sex_specific_${ph}_${chr}_ld12_stdout \
-e ${path_gwas}post_imputation/2022/log/cojo_eur_sex_specific_${ph}_${chr}_ld12_stderr \
"/path/to/software/username/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_list_tier_1_sex_specific_conditional_variants_for_ld \
--r2 --ld-window-kb '1000000000' --ld-window '1000000000' --ld-window-r2 '0' --out ${path_gwas}post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr${chr}_list_tier_1_sex_specific_conditional_variants_for_ld"
done

for chr in {1..22} 
do
echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/log/cojo_eur_sex_specific_${ph}_${chr}_ld12_stdout | grep "Successfully"
done


##############################################################################################################


for (chr in 1:22) {
  ld.tmp<-read.table(paste0(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr",chr,"_list_tier_1_sex_specific_conditional_variants_for_ld.ld"),head=T)
  ld.tmp<-ld.tmp[which(ld.tmp$R2>=0.1),]

  if(chr==1) {
    ld<-ld.tmp
  } else {
    ld<-rbind(ld,ld.tmp)
  }
  rm(ld.tmp)
}

dim(ld)
# [1] 431   7

ld1<-ld[which( (ld$SNP_A %in% dat_final$SNP) & (ld$SNP_B %in% dat_final$SNP)),]
dim(ld1)
# [1] 246   7


dim(dat_final)
# [1] 318  39

list_snps<-c(ld1$SNP_A,ld1$SNP_B)
list_snps<-list_snps[!duplicated(list_snps)]

dat_final$variants_in_ld_r2_sex_specific_0.1<-NA
dat_final$independent_sex_specific<-0

for (i in 1:length(list_snps)) {

  # extract all LD pairs
  tmp1<-ld1[which(ld1$SNP_A==list_snps[i] | ld1$SNP_B==list_snps[i]),]
  tmp1$variants_in_ld_r2_sex_specific_0.1<-paste(tmp1$SNP_A,tmp1$SNP_B,tmp1$R2,sep="|")

  dat_final$variants_in_ld_r2_sex_specific_0.1[which(dat_final$SNP %in% c(tmp1$SNP_A,tmp1$SNP_B))]<-paste(tmp1$variants_in_ld_r2_sex_specific_0.1,collapse=";")

  # retain only the most significant:
  tmp<-dat_final[which(dat_final$SNP %in% c(tmp1$SNP_A,tmp1$SNP_B)),]

    vec<-c("SNP",colnames(dat_final)[grep("Pvalue",colnames(dat_final))])
    tmp<-tmp[,..vec]
    minpval<-min(tmp[,2:ncol(tmp)],na.rm=t)

    vec<-c("SNP",names(tmp)[which(tmp == minpval, arr.ind = TRUE)[, "col"]])
    tmp<-tmp[,..vec]
    colnames(tmp)[2]<-"Pvalue"

    tmp<-tmp[which(tmp$Pvalue==minpval),]

  dat_final$independent_sex_specific[which(dat_final$SNP==tmp$SNP)]<-1

}

table(dat_final$independent_sex_specific)
  # 0   1 
# 230  88 


# re-label those variants with no friends in LD, as independent:
dim(dat_final[which(is.na(dat_final$variants_in_ld_r2_sex_specific_0.1)),])
# [1] 83 39

dat_final$independent_sex_specific[which(is.na(dat_final$variants_in_ld_r2_sex_specific_0.1))]<-1

table(dat_final$independent_sex_specific)
#   0   1 
# 147 171 

# fwrite(dat_final,
# paste0(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/list_tier_1_unsupervised_conditional_sex_specific_variants_female_male_only.tsv.gz"),
# col.names=T,row.names=F)

dat_final<-as.data.frame(dat_final)

ld1<-ld[which( (ld$SNP_A %in% dat_final$SNP & ld$SNP_B %in% all$MarkerName) | (ld$SNP_B %in% dat_final$SNP & ld$SNP_A %in% all$MarkerName)),]
dim(ld1)
# 234

dat_final$variants_in_ld_r2_sex_specific_vs_nonspecific_0.1<-NA
dat_final$consensus_index<-0

for (i in 1:nrow(ld1)) {

  tmp<-c(ld1$SNP_A[i],ld1$SNP_B[i])

  # tmp2<-dat_final$variants_in_ld_r2_sex_specific_vs_nonspecific_0.1[which(dat_final$SNP %in% tmp)]
  # if(length(tmp2)>1) {
  #   print(i)
  # }
  
  dat_final$variants_in_ld_r2_sex_specific_vs_nonspecific_0.1[which(dat_final$SNP %in% tmp)]<-rep(paste(c(as.character(ld1$SNP_A[i]),as.character(ld1$SNP_B[i]),ld1$R2[i]),collapse="|"),length(dat_final$variants_in_ld_r2_sex_specific_vs_nonspecific_0.1[which(dat_final$SNP %in% tmp)]))

  # retain the consensus Index SNP:
  
  tmp1<-paste(tmp[which(tmp %in% all$MarkerName)],collapse="|")
  # if(length(tmp[which(tmp %in% all$MarkerName)])>1) {
  #   print(i)
  # }
  dat_final$consensus_index[which(dat_final$SNP %in% tmp)]<-rep(tmp1,length(dat_final$consensus_index[which(dat_final$SNP %in% tmp)]))

}


# variants tagging an index SNP
dim(dat_final[which(dat_final$consensus_index!=0),])
# [1] 250  41


# add those sex_specific variants in the original
dim(dat_final[which(dat_final$independent_sex_specific==1 & dat_final$consensus_index==0),])
# [1] 49 41

dat_final$consensus_index[which(dat_final$SNP %in% all$MarkerName & dat_final$consensus_index==0)]<-dat_final$SNP[which(dat_final$SNP %in% all$MarkerName & dat_final$consensus_index==0)]

dim(dat_final[which(dat_final$independent_sex_specific==1 & dat_final$consensus_index==0),])
# [1] 24 45



# label HLA region to be excluded:
dat_final$HLA<-0
dat_final$Position<-gsub("chr[0-9]{1,2}:","",dat_final$SNP)
dat_final$Position<-as.numeric(gsub(":[A-Z]*:[A-Z]*$","",dat_final$Position))
dat_final$chr<-gsub(":[0-9]*:[A-Z]*:[A-Z]*$","",dat_final$SNP)
dat_final$chr<-as.numeric(gsub("chr","",dat_final$chr))

# remove variants in MHC: https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38.p11
# MHC -- chr6 (NC_000006.12):28,510,120-33,480,577
dat_final$HLA[which(dat_final$chr=="6" & dat_final$Position>=28510120-1500000 & dat_final$Position<=33480577+1500000)]<-1
#0


dim(dat_final[which(dat_final$independent_sex_specific==1 & dat_final$consensus_index==0 & dat_final$HLA==0),])
# [1]  12 45

## label those that fall within the consensus list of regions:
regions<-fread(paste0(path_gwas,"post_imputation/2022/analysis/metaanalysis/forward_regression/list_regions_420_cojo_supervised_plus_unsupervised.tsv.gz"))
regions$min<-gsub("^[0-9]{1,2}_","",regions$updated_region)
regions$max<-as.numeric(gsub("^[0-9]*_","",regions$min))
regions$min<-as.numeric(gsub("_[0-9]*$","",regions$min))
regions$chr<-gsub("_.*","",regions$updated_region)

regions<-as.data.frame(regions)

dat_final<-as.data.frame(dat_final)

dat_final$consensus_region<-NA
for (i in 1:nrow(dat_final)) {

  tmp<-regions[which(regions$chr==dat_final$chr[i] & regions$min<=dat_final$Position[i] & regions$max>=dat_final$Position[i]),]

  if(nrow(tmp)>0) {
    # print(i)
    dat_final$consensus_region[i]<-tmp$updated_region
  }
  
}

dat_final[which(is.na(dat_final$consensus_region) & dat_final$HLA==0),]
#    Position chr consensus_region
# 13 72606376  10             <NA>


dat_final[which(dat_final$independent_sex_specific==1 & dat_final$consensus_index==0 & dat_final$HLA==0 & !is.na(dat_final$consensus_region)),]
#                     SNP BETA_male_uc SE_male_uc Pvalue_male_uc
# 7    chr10:62617743:G:A           NA         NA             NA
# 30   chr11:64249945:A:G           NA         NA             NA
# 65   chr16:27399508:A:G           NA         NA             NA
# 74   chr16:50712383:A:C           NA         NA             NA
# 112  chr19:45898774:T:C           NA         NA             NA
# 114  chr19:46640496:C:T           NA         NA             NA
# 167 chr22:29426322:CA:C           NA         NA             NA
# 199   chr2:25097168:C:T           NA         NA             NA
# 238   chr5:40011635:C:T           NA         NA             NA
# 263   chr6:20809563:A:T           NA         NA             NA
# 310  chr9:114814020:G:A           NA         NA             NA




for (chr in 1:22) {
  ld.tmp<-read.table(paste0(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/allarrays_chr",chr,"_list_tier_1_sex_specific_conditional_variants_for_ld.ld"),head=T)
  # ld.tmp<-ld.tmp[which(ld.tmp$R2>=0.1),]

  if(chr==1) {
    ld<-ld.tmp
  } else {
    ld<-rbind(ld,ld.tmp)
  }
  rm(ld.tmp)
}


# chr9:114814020:G:A sex specific; chr9:114881082:C:T FM causal;  "chr9:114801407:T:G","chr9:114845889:A:G","chr9:114851353:GA:G" COJO indep 
vec<-c("chr9:114814020:G:A","chr9:114801407:T:G","chr9:114845889:A:G","chr9:114851353:GA:G","chr9:114881082:C:T")
ld[which(ld$SNP_A %in% vec & ld$SNP_B %in% vec),]
#       CHR_A      BP_A               SNP_A CHR_B      BP_B               SNP_B
# 17646     9 114801407  chr9:114801407:T:G     9 114814020  chr9:114814020:G:A
# 17648     9 114801407  chr9:114801407:T:G     9 114845889  chr9:114845889:A:G
# 17649     9 114801407  chr9:114801407:T:G     9 114851353 chr9:114851353:GA:G
# 17650     9 114801407  chr9:114801407:T:G     9 114881082  chr9:114881082:C:T
# 17681     9 114814020  chr9:114814020:G:A     9 114845889  chr9:114845889:A:G
# 17682     9 114814020  chr9:114814020:G:A     9 114851353 chr9:114851353:GA:G
# 17683     9 114814020  chr9:114814020:G:A     9 114881082  chr9:114881082:C:T
# 17711     9 114845889  chr9:114845889:A:G     9 114851353 chr9:114851353:GA:G
# 17712     9 114845889  chr9:114845889:A:G     9 114881082  chr9:114881082:C:T
# 17725     9 114851353 chr9:114851353:GA:G     9 114881082  chr9:114881082:C:T
#               R2
# 17646 0.05635010
# 17648 0.08910670
# 17649 0.00781178
# 17650 0.11597100
# 17681 0.06119440
# 17682 0.00380409
# 17683 0.05519570
# 17711 0.00884351
# 17712 0.12839100
# 17725 0.00724429


# chr10:62617743:G:A - SEX specific, chr10:62710915:C:T, FM causal; chr10:62710915:C:T COJO indep
vec<-c("chr10:62617743:G:A","chr10:62710915:C:T")
ld[which(ld$SNP_A %in% vec & ld$SNP_B %in% vec),]
#       CHR_A     BP_A              SNP_A CHR_B     BP_B              SNP_B
# 18484    10 62617743 chr10:62617743:G:A    10 62710915 chr10:62710915:C:T
#                R2
# 18484 0.000878164

# chr11:64249945:A:G sex specific; chr11:63796342:TC:T chr11:63929911:G:A COJO indep
vec<-c("chr11:64249945:A:G","chr11:63796342:TC:T","chr11:63929911:G:A")
ld[which(ld$SNP_A %in% vec & ld$SNP_B %in% vec),]
#       CHR_A     BP_A               SNP_A CHR_B     BP_B              SNP_B
# 19431    11 63796342 chr11:63796342:TC:T    11 63929911 chr11:63929911:G:A
# 19432    11 63796342 chr11:63796342:TC:T    11 64249945 chr11:64249945:A:G
# 19457    11 63929911  chr11:63929911:G:A    11 64249945 chr11:64249945:A:G
#                R2
# 19431 0.000393385
# 19432 0.010781400
# 19457 0.014317400

# chr16:27399508:A:G sex independent; chr16:27384341:C:CT COJO indep
vec<-c("chr16:27399508:A:G","chr16:27384341:C:CT")
ld[which(ld$SNP_A %in% vec & ld$SNP_B %in% vec),]
#       CHR_A     BP_A               SNP_A CHR_B     BP_B              SNP_B
# 21267    16 27384341 chr16:27384341:C:CT    16 27399508 chr16:27399508:A:G
#                R2
# 21267 0.000498677




# fwrite(dat_final,
# paste0(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/list_tier_1_unsupervised_conditional_sex_specific_variants.tsv.gz"),
# col.names=T,row.names=F)

fwrite(dat_final,
paste0(path_gwas,"post_imputation/2022/analysis/conditional_analysis/eur/sex_specific_analyses/list_tier_1_unsupervised_conditional_sex_specific_variants_female_male_only_202600310.tsv"),
col.names=T,row.names=F,sep="\t")


q("no")

################





