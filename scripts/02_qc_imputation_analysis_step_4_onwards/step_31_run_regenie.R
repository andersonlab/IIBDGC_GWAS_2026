# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
##########
# STEP 1 #
##########

path_gwas=/path/to/ibdgwas/IIBDGC/
  
MEM=4000 

array=(illumina370 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 gsa illuminaexome)
ancestry=(eur_nonjewish eur_all)

array=(illumina370 gsa)
ancestry=(eur_jewish)

pheno=(ibd cd uc)

for i in ${array[@]}
do
for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
bsub -J"reg_ibd" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -n 16 \
-e ${path_gwas}post_imputation/2022/analysis/regenie/log/${i}_${j}_${ph}_regenie_step1_sex_PCs_stderr \
-o ${path_gwas}post_imputation/2022/analysis/regenie/log/${i}_${j}_${ph}_regenie_step1_sex_PCs_stdout \
"/path/to/software/./regenie_v3.2.5.gz_x86_64_Linux \
--step 1 \
--bed ${path_gwas}post_imputation/2022/${i}/genotyped_data/${i}_all_studies_merged_eur_pruned \
--phenoFile ${path_gwas}post_imputation/2022/${i}/phenotype_data/${i}_all_studies_merged_${j}_phenotype_${ph} \
--covarFile ${path_gwas}post_imputation/2022/${i}/phenotype_data/${i}_all_studies_merged_${j}_covariate_sex_PCs \
--threads 16 \
--strict \
--bt \
--bsize 1000 \
--loocv \
--lowmem-prefix ${path_gwas}post_imputation/2022/analysis/regenie/${i}/${ph}/tmpdir/${i}_${j}_regenie_tmp_preds_${ph}_eur_sex_PCs \
--out ${path_gwas}post_imputation/2022/analysis/regenie/${i}/${ph}/${i}_${j}_step1_ibd_eur_sex_PCs"
done
done
done

ancestry=(eur_all eur_nonjewish eur_jewish)

for i in ${array[@]}
do echo ${i} && for j in ${ancestry[@]}
do echo ${j} && for ph in ${pheno[@]}
do echo ${ph} && tail -50 ${path_gwas}post_imputation/2022/analysis/regenie/log/${i}_${j}_${ph}_regenie_step1_sex_PCs_stdout | grep "Successfully"
done
done
done


############################

# illumina370 # NOT ENOUGH EUR JEWISH DATA FOR UC
# eur_all
# ibd
# Successfully completed.
# cd
# Successfully completed.
# uc
# Successfully completed.
# eur_nonjewish
# ibd
# Successfully completed.
# cd
# Successfully completed.
# uc
# Successfully completed.
# eur_jewish
# ibd
# Successfully completed.
# cd
# Successfully completed.
# uc

# affymetrix6
# eur_all
# ibd
# Successfully completed.
# cd
# Successfully completed.
# uc
# Successfully completed.
# eur_nonjewish
# ibd
# Successfully completed.
# cd
# Successfully completed.
# uc
# Successfully completed.
# eur_jewish
# ibd
# tail: cannot open '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/log/affymetrix6_eur_jewish_ibd_regenie_step1_sex_PCs_stdout' for reading: No such file or directory
# cd
# tail: cannot open '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/log/affymetrix6_eur_jewish_cd_regenie_step1_sex_PCs_stdout' for reading: No such file or directory
# uc
# tail: cannot open '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/log/affymetrix6_eur_jewish_uc_regenie_step1_sex_PCs_stdout' for reading: No such file or directory


# humanomniexpress # OK NO CD DATA
# eur_all
# ibd
# Successfully completed.
# cd
# uc
# Successfully completed.
# eur_nonjewish
# ibd
# Successfully completed.
# cd
# uc
# Successfully completed.
# eur_jewish
# ibd
# tail: cannot open '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/log/humanomniexpress_eur_jewish_ibd_regenie_step1_sex_PCs_stdout' for reading: No such file or directory
# cd
# tail: cannot open '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/log/humanomniexpress_eur_jewish_cd_regenie_step1_sex_PCs_stdout' for reading: No such file or directory
# uc
# tail: cannot open '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/log/humanomniexpress_eur_jewish_uc_regenie_step1_sex_PCs_stdout' for reading: No such file or directory


# quad610
# eur_all
# ibd
# Successfully completed.
# cd
# Successfully completed.
# uc
# Successfully completed.
# eur_nonjewish
# ibd
# Successfully completed.
# cd
# Successfully completed.
# uc
# Successfully completed.
# eur_jewish
# ibd
# tail: cannot open '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/log/quad610_eur_jewish_ibd_regenie_step1_sex_PCs_stdout' for reading: No such file or directory
# cd
# tail: cannot open '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/log/quad610_eur_jewish_cd_regenie_step1_sex_PCs_stdout' for reading: No such file or directory
# uc
# tail: cannot open '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/log/quad610_eur_jewish_uc_regenie_step1_sex_PCs_stdout' for reading: No such file or directory


# gsa
# eur_all
# ibd
# Successfully completed.
# cd
# Successfully completed.
# uc
# Successfully completed.
# eur_nonjewish
# ibd
# Successfully completed.
# cd
# Successfully completed.
# uc
# Successfully completed.
# eur_jewish
# ibd
# Successfully completed.
# cd
# Successfully completed.
# uc
# Successfully completed.

# illuminaexome
# eur_all
# ibd
# Successfully completed.
# cd
# Successfully completed.
# uc
# Successfully completed.
# eur_nonjewish
# ibd
# Successfully completed.
# cd
# Successfully completed.
# uc
# Successfully completed.
# eur_jewish
# ibd
# tail: cannot open '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/log/illuminaexome_eur_jewish_ibd_regenie_step1_sex_PCs_stdout' for reading: No such file or directory
# cd
# tail: cannot open '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/log/illuminaexome_eur_jewish_cd_regenie_step1_sex_PCs_stdout' for reading: No such file or directory
# uc
# tail: cannot open '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/log/illuminaexome_eur_jewish_uc_regenie_step1_sex_PCs_stdout' for reading: No such file or directory
############################


##########
# STEP 2 #
##########

## autosomal

MEM=8500
# array=(gsa) # queueu basement

# MEM=3000

# array=(illumina370 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 gsa illuminaexome gsa)
array=(gsa)
ancestry=(eur_nonjewish)

pheno=(cd ibd uc)

for i in ${array[@]}
do
for j in ${ancestry[@]}
do
for ph in ${pheno[@]}
do
for chr in 7
do 
bsub -J"reg_ibd" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -n 16 -q long \
-e ${path_gwas}post_imputation/2022/analysis/regenie/log/${i}_${j}_${ph}_chr${chr}_regenie_step2_sex_PCs_stderr \
-o ${path_gwas}post_imputation/2022/analysis/regenie/log/${i}_${j}_${ph}_chr${chr}_regenie_step2_sex_PCs_stdout \
"/path/to/software/./regenie_v3.2.5.gz_x86_64_Linux \
--step 2 \
--bgen ${path_gwas}post_imputation/2022/${i}/imputed_data/${i}_chr${chr}.dose.subset.bgen \
--sample ${path_gwas}post_imputation/2022/${i}/imputed_data/${i}_chr${chr}.dose.subset.sample \
--phenoFile ${path_gwas}post_imputation/2022/${i}/phenotype_data/${i}_all_studies_merged_${j}_phenotype_${ph} \
--covarFile ${path_gwas}post_imputation/2022/${i}/phenotype_data/${i}_all_studies_merged_${j}_covariate_sex_PCs \
--bt \
--ref-first \
--strict \
--firth \
--approx \
--firth-se \
--pThresh 0.1 \
--pred ${path_gwas}post_imputation/2022/analysis/regenie/${i}/${ph}/${i}_${j}_step1_${ph}_eur_sex_PCs_pred.list \
--threads 16 \
--bsize 400 \
--minMAC 1 \
--af-cc \
--out ${path_gwas}post_imputation/2022/analysis/regenie/${i}/${ph}/chr${chr}_${i}_${j}_step2_${ph}_eur_sex_PCs_firthse"
done
done
done
done





array=(illumina370 affymetrix6 humanomniexpress affymetrix500 humancoreexome humanomni1 quad610 gsa illuminaexome)
pheno=(ibd cd uc)

j=eur_all

for i in ${array[@]}
do echo ${i} && for ph in ${pheno[@]}
do echo ${ph} && for chr in {1..22} X
do echo ${chr} && tail -50 ${path_gwas}post_imputation/2022/analysis/regenie/log/${i}_${j}_${ph}_chr${chr}_regenie_step2_sex_PCs_stdout | grep "Successfully"
done
done
done


#############################

# illumina370
# ibd
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# cd
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# uc
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# affymetrix6
# ibd
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# Successfully completed.
# cd
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# Successfully completed.
# uc
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# Successfully completed.
# humanomniexpress
# ibd
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# Successfully completed.
# cd
# 1
# 2
# 3
# 4
# 5
# 6
# 7
# 8
# 9
# 10
# 11
# 12
# 13
# 14
# 15
# 16
# 17
# 18
# 19
# 20
# 21
# 22
# X
# uc
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# Successfully completed.
# affymetrix500
# ibd
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# tail: cannot open '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/log/affymetrix500_eur_all_ibd_chrX_regenie_step2_sex_PCs_stdout' for reading: No such file or directory
# cd
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# tail: cannot open '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/log/affymetrix500_eur_all_cd_chrX_regenie_step2_sex_PCs_stdout' for reading: No such file or directory
# uc
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# tail: cannot open '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/log/affymetrix500_eur_all_uc_chrX_regenie_step2_sex_PCs_stdout' for reading: No such file or directory
# humancoreexome
# ibd
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# tail: cannot open '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/log/humancoreexome_eur_all_ibd_chrX_regenie_step2_sex_PCs_stdout' for reading: No such file or directory
# cd
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# Successfully completed.
# uc
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# tail: cannot open '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/log/humancoreexome_eur_all_uc_chrX_regenie_step2_sex_PCs_stdout' for reading: No such file or directory
# humanomni1
# ibd
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# tail: cannot open '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/log/humanomni1_eur_all_ibd_chrX_regenie_step2_sex_PCs_stdout' for reading: No such file or directory
# cd
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# tail: cannot open '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/log/humanomni1_eur_all_cd_chrX_regenie_step2_sex_PCs_stdout' for reading: No such file or directory
# uc
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# tail: cannot open '/path/to/ibdgwas/IIBDGC/post_imputation/2022/analysis/regenie/log/humanomni1_eur_all_uc_chrX_regenie_step2_sex_PCs_stdout' for reading: No such file or directory
# quad610
# ibd
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# Successfully completed.
# cd
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# Successfully completed.
# uc
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# Successfully completed.
# gsa
# ibd
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# cd
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# uc
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# illuminaexome
# ibd
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# cd
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X
# uc
# 1
# Successfully completed.
# 2
# Successfully completed.
# 3
# Successfully completed.
# 4
# Successfully completed.
# 5
# Successfully completed.
# 6
# Successfully completed.
# 7
# Successfully completed.
# 8
# Successfully completed.
# 9
# Successfully completed.
# 10
# Successfully completed.
# 11
# Successfully completed.
# 12
# Successfully completed.
# 13
# Successfully completed.
# 14
# Successfully completed.
# 15
# Successfully completed.
# 16
# Successfully completed.
# 17
# Successfully completed.
# 18
# Successfully completed.
# 19
# Successfully completed.
# 20
# Successfully completed.
# 21
# Successfully completed.
# 22
# Successfully completed.
# X


#############################
