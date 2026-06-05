#!/bin/bash
# Author: Talin Haritunians
# Institution: F. Widjaja Inflammatory Bowel Disease Institute
#

## sex-specific analyses
## regenie step 1

## dataset filtered maf >5%; pruned in plink1.9 with --indep 50 5 10

for pheno in ibd cd uc; do

	regenie --step 1 --pgen ${array}_chrALL.dose_pruned_maf05 --phenoFile ${array}_eur_${pheno}.txt --covarFile ${array}_eur_covar.txt --keep ${array}_eur_females.txt --strict --bt --bsize 1000 --loocv --lowmem --out step1_${array}_chrALL.eur_female_${pheno}

done


for pheno in ibd cd uc; do

regenie --step 1 --pgen ${array}_chrALL.dose_pruned_maf05 --phenoFile ${array}_eur_${pheno}.txt --covarFile ${array}_eur_covar.txt --remove ${array}_eur_females.txt --strict --bt --bsize 1000 --loocv --lowmem --out step1_${array}_chrALL.eur_male_${pheno}

done




## regenie step 2

for pheno in ibd cd uc; do

	for i in {1..22}; do
	
		regenie --step 2 --pgen ${array}_chr${i}.dose.eur --phenoFile ${array}_eur_${pheno}.txt --covarFile ${array}_eur_covar.txt --keep ${array}_eur_females.txt --bt --firth --approx --firth-se --pThresh 0.1 --pred step1_${array}_chrALL.eur_female_${pheno}_pred.list --bsize 400 --minMAC 1 --write-samples --af-cc --out step2_${array}_chr${i}.female.${pheno}_eur

	done

done


for pheno in ibd cd uc; do

	for i in {1..22}; do
	
		regenie --step 2 --pgen ${array}_chr${i}.dose.eur --phenoFile ${array}_eur_${pheno}.txt --covarFile ${array}_eur_covar.txtt --remove ${array}_eur_females.txt --bt --firth --approx --firth-se --pThresh 0.1 --pred step1_${array}_chrALL.eur_male_${pheno}_pred.list --bsize 400 --minMAC 1 --write-samples --af-cc --out step2_${array}_chr${i}.male.${pheno}_eur

	done

done
