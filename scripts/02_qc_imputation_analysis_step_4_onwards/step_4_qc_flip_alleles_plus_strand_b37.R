# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
#########################
# 4.- ALIGN TO + STRAND #
#########################

################################################################################################
# STEP NO LONGER REQUIRED - NO CHR25 IN NEW BATCH OF SAMPLES - REQUIRED FOR:

##########################################################
# 4.1 GENERATE VCF (output-chr M - reformat chr23 to X)


path_gwas=/path/to/ibdgwas/IIBDGC/

studies=(australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

# recatch all_hce pre-check sex step
# recatch niddk_old_gwas  pre-check sex step

for i in ${studies[@]}
do
ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind.bed
done


##################################
# intermediate step to sort chr23:

MEM=1200

for i in ${studies[@]}
do
bsub -J"pl1" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_1_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_1_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind \
--allow-no-sex \
--merge-x no-fail \
--threads 2 \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_edited2"
done
# Job <1921299..1921360> is submitted to queue <normal>.


for i in ${studies[@]}
do
echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_1_${i} | grep "Successfully completed."
done

#####
# or i in ${studies[@]}
# > do echo ${i} && wc -l ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_edited2.bim
# > done
# australia_omniexome
# 807107 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/australia_omniexome/australia_omniexome_hg19_noind_edited2.bim
# gwas1
# 459426 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas1/gwas1_hg19_noind_edited2.bim
# gwas2
# 799264 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas2/gwas2_hg19_noind_edited2.bim
# pittsburgh_gsa
# 945375 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pittsburgh_gsa/pittsburgh_gsa_hg19_noind_edited2.bim
# spain_gsa
# 585734 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/spain_gsa/spain_gsa_hg19_noind_edited2.bim
# italy_gsa
# 586354 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/italy_gsa/italy_gsa_hg19_noind_edited2.bim
# kiel_austria_sibdcs_gsa
# 663705 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_hg19_noind_edited2.bim
# netherlands_gsa
# 684826 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/netherlands_gsa/netherlands_gsa_hg19_noind_edited2.bim
# slovenia_gsa
# 557690 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/slovenia_gsa/slovenia_gsa_hg19_noind_edited2.bim
# sweden_gsa
# 575192 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sweden_gsa/sweden_gsa_hg19_noind_edited2.bim
# niddk_broad_gsa
# 622147 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_hg19_noind_edited2.bim
# niddk_feinstein_gsa
# 645538 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_hg19_noind_edited2.bim
# basque_gsa
# 579360 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/basque_gsa/basque_gsa_hg19_noind_edited2.bim
# lithuania_gsa
# 580098 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lithuania_gsa/lithuania_gsa_hg19_noind_edited2.bim
# belgium_louis_gsa
# 605806 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_louis_gsa/belgium_louis_gsa_hg19_noind_edited2.bim
# belgium_franchimont_gsa
# 585266 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_franchimont_gsa/belgium_franchimont_gsa_hg19_noind_edited2.bim
# belgium_vermeire_gsa
# 674795 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_vermeire_gsa/belgium_vermeire_gsa_hg19_noind_edited2.bim
# prism_nfe_gsa
# 573760 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_hg19_noind_edited2.bim
# prism_nfe_gwas
# 243983 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gwas/prism_nfe_gwas_hg19_noind_edited2.bim
# finland_illugwas
# 239244 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/finland_illugwas/finland_illugwas_hg19_noind_edited2.bim
# german_affy6_old_gwas
# 893289 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/german_affy6_old_gwas/german_affy6_old_gwas_hg19_noind_edited2.bim
# norway_affy6_old_gwas
# 846288 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg19_noind_edited2.bim
# belgium_inf1_old_gwas
# 302700 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg19_noind_edited2.bim
# belgium_inf2_old_gwas
# 290730 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg19_noind_edited2.bim
# cedars_370k_old_gwas
# 339977 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg19_noind_edited2.bim
# cedars_610k_old_gwas
# 586380 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg19_noind_edited2.bim
# cedars_omni_old_gwas
# 723004 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg19_noind_edited2.bim
# swedish_uc_old_gwas
# 297004 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/swedish_uc_old_gwas_hg19_noind_edited2.bim
# mccauley_gsa
# 571173 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_gsa/mccauley_gsa_hg19_noind_edited2.bim
# ccfa_gsa
# 611354 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/ccfa_gsa/ccfa_gsa_hg19_noind_edited2.bim
# cedars_gsa
# 669422 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_gsa/cedars_gsa_hg19_noind_edited2.bim
# bernstein_gsa
# 565758 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/bernstein_gsa/bernstein_gsa_hg19_noind_edited2.bim
# farkkila_gsa
# 528311 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/farkkila_gsa/farkkila_gsa_hg19_noind_edited2.bim
# franchimont_gsa
# 671132 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franchimont_gsa/franchimont_gsa_hg19_noind_edited2.bim
# franke_gsa
# 571118 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franke_gsa/franke_gsa_hg19_noind_edited2.bim
# helmsley_prism_gsa
# 243983 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_prism_gsa/helmsley_prism_gsa_hg19_noind_edited2.bim
# helmsley_xavier_prism_gsa
# 243952 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_xavier_prism_gsa/helmsley_xavier_prism_gsa_hg19_noind_edited2.bim
# hyams_protect_gsa
# 580928 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/hyams_protect_gsa/hyams_protect_gsa_hg19_noind_edited2.bim
# lewis_sparc_gsa
# 619895 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lewis_sparc_gsa/lewis_sparc_gsa_hg19_noind_edited2.bim
# mccauley_new_gsa
# 591351 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_new_gsa/mccauley_new_gsa_hg19_noind_edited2.bim
# mcgovern_gsa
# 652046 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_hg19_noind_edited2.bim
# moayyedi_imagine_gsa
# 615521 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/moayyedi_imagine_gsa/moayyedi_imagine_gsa_hg19_noind_edited2.bim
# newberry_share_gsa
# 592708 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/newberry_share_gsa/newberry_share_gsa_hg19_noind_edited2.bim
# niddk_cho_gsa
# 608759 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_cho_gsa/niddk_cho_gsa_hg19_noind_edited2.bim
# niddk_duerr_gsa
# 596287 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_duerr_gsa/niddk_duerr_gsa_hg19_noind_edited2.bim
# niddk_rioux_gsa
# 589690 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_rioux_gsa/niddk_rioux_gsa_hg19_noind_edited2.bim
# niddk_silverberg_gsa
# 626861 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_silverberg_gsa/niddk_silverberg_gsa_hg19_noind_edited2.bim
# palotie_hus_gsa
# 565416 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/palotie_hus_gsa/palotie_hus_gsa_hg19_noind_edited2.bim
# pekow_share_gsa
# 581240 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pekow_share_gsa/pekow_share_gsa_hg19_noind_edited2.bim
# rioux_igenomed_gsa
# 551803 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_edited2.bim
# sands_msccr_gsa
# 604196 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sands_msccr_gsa/sands_msccr_gsa_hg19_noind_edited2.bim
# stampfer_gsa
# 672634 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/stampfer_gsa/stampfer_gsa_hg19_noind_edited2.bim
# vermeire_gsa
# 675377 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/vermeire_gsa/vermeire_gsa_hg19_noind_edited2.bim
# weersma_gsa
# 573894 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/weersma_gsa/weersma_gsa_hg19_noind_edited2.bim
# xavier_prism_gsa
# 584364 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_prism_gsa/xavier_prism_gsa_hg19_noind_edited2.bim
# xavier_share_gsa
# 584949 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_share_gsa/xavier_share_gsa_hg19_noind_edited2.bim
#####

##################################
# make vcf

MEM=1000
for i in ${studies[@]}
do
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_2_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_2_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_edited2 \
--allow-no-sex \
--threads 2 \
--output-chr M --recode vcf-iid --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19"
done
# Job <1922677..1922733> is submitted to queue <normal>.

for i in ${studies[@]}
do
echo ${i} && tail -200 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_2_${i} | grep "Successfully completed."
done

##################################
# compress file 

MEM=1000
for i in ${studies[@]}
do
bsub -J"bgz" -M"$MEM" -n 4 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_bgzip_1_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_bgzip_1_${i} \
"/path/to/software/htslib-1.11/./bgzip --threads 4 ${path_gwas}pre_imputation/QC/${i}/${i}_hg19.vcf"
done
# Job <1922752..1922808> is submitted to queue <normal>.

for i in ${studies[@]}
do
echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_bgzip_1_${i} | grep "Successfully completed."
done


for i in ${studies[@]}
do
ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_hg19.vcf.gz 
done

###########################################################
# 4.2 FLIP ALLELES TO + STRAND:

path_gwas=/path/to/ibdgwas/IIBDGC/
export BCFTOOLS_PLUGINS=/path/to/software/bcftools-1.16/plugins
MEM=500

# ALLL resubmitted - running now

studies=(australia_omniexome gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa)
studies=(niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas)
studies=(norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas xavier_prism_gsa xavier_share_gsa) 
studies=(cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa)
studies=(helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa)
studies=(niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa)


for i in ${studies[@]}
do
bsub -J"bcf" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_bcftools_1_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_bcftools_1_${i} \
"/path/to/software/bcftools-1.16/./bcftools +fixref ${path_gwas}pre_imputation/QC/${i}/${i}_hg19.vcf.gz \
--threads 1 -Oz -o ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_posstrandaligned.vcf.gz \
-- -f ${path_gwas}resources/1000GP/human_g1k_v37.fasta -m top"
done

for i in ${studies[@]}
do
echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_bcftools_1_${i} | grep "Successfully completed."
done

for i in ${studies[@]}
do
echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_posstrandaligned.vcf.gz
done

for i in ${studies[@]}
do
echo ${i} && tail -24 ${path_gwas}pre_imputation/QC/${i}/logs/stderr_bcftools_1_${i}
done

#################################
# australia_omniexome - OK
# ST	A>G	286896	35.5%
# ST	A>T	2012	0.2%
# ST	C>A	77261	9.6%
# ST	C>G	4822	0.6%
# ST	C>T	0	0.0%
# ST	G>A	360212	44.6%
# ST	G>C	5239	0.6%
# ST	G>T	0	0.0%
# ST	T>A	2636	0.3%
# ST	T>C	0	0.0%
# ST	T>G	0	0.0%
# # NS, Number of sites:
# NS	total        	807107
# NS	ref match    	299628	37.1%
# NS	ref mismatch 	507479	62.9%
# NS	flipped      	298646	37.0%
# NS	swapped      	104851	13.0%
# NS	flip+swap    	104715	13.0%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# gwas1 - OK
# ST	A>G	71627	15.6%
# ST	A>T	15147	3.3%
# ST	C>A	19087	4.2%
# ST	C>G	22552	4.9%
# ST	C>T	85316	18.6%
# ST	G>A	84978	18.5%
# ST	G>C	22502	4.9%
# ST	G>T	19011	4.1%
# ST	T>A	14880	3.2%
# ST	T>C	71301	15.5%
# ST	T>G	16562	3.6%
# # NS, Number of sites:
# NS	total        	459426
# NS	ref match    	320451	69.8%
# NS	ref mismatch 	138975	30.2%
# NS	flipped      	11516	2.5%
# NS	swapped      	127440	27.7%
# NS	flip+swap    	25982	5.7%
# NS	unresolved   	1	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# gwas2 - OK
# ST	A>G	124412	15.6%
# ST	A>T	26181	3.3%
# ST	C>A	37153	4.6%
# ST	C>G	36063	4.5%
# ST	C>T	156303	19.6%
# ST	G>A	145037	18.1%
# ST	G>C	36645	4.6%
# ST	G>T	30650	3.8%
# ST	T>A	24859	3.1%
# ST	T>C	124162	15.5%
# ST	T>G	24095	3.0%
# # NS, Number of sites:
# NS	total        	799264
# NS	ref match    	300526	37.6%
# NS	ref mismatch 	498737	62.4%
# NS	flipped      	271329	33.9%
# NS	swapped      	128672	16.1%
# NS	flip+swap    	128269	16.0%
# NS	unresolved   	1	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	1
# NS	non-ACGT     	1
# NS	non-SNP      	0
# NS	non-biallelic	0
# pittsburgh_gsa - OK
# ST	A>G	340623	36.0%
# ST	A>T	4294	0.5%
# ST	C>A	91657	9.7%
# ST	C>G	7255	0.8%
# ST	C>T	0	0.0%
# ST	G>A	405916	42.9%
# ST	G>C	8191	0.9%
# ST	G>T	0	0.0%
# ST	T>A	4922	0.5%
# ST	T>C	0	0.0%
# ST	T>G	0	0.0%
# # NS, Number of sites:
# NS	total        	945375
# NS	ref match    	349246	36.9%
# NS	ref mismatch 	596058	63.1%
# NS	flipped      	347834	36.8%
# NS	swapped      	125444	13.3%
# NS	flip+swap    	125516	13.3%
# NS	unresolved   	12	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	71
# NS	non-ACGT     	71
# NS	non-SNP      	0
# NS	non-biallelic	0
# spain_gsa - OK
# ST	A>G	220475	37.6%
# ST	A>T	581	0.1%
# ST	C>A	56920	9.7%
# ST	C>G	922	0.2%
# ST	C>T	0	0.0%
# ST	G>A	253231	43.2%
# ST	G>C	981	0.2%
# ST	G>T	0	0.0%
# ST	T>A	637	0.1%
# ST	T>C	0	0.0%
# ST	T>G	0	0.0%
# # NS, Number of sites:
# NS	total        	585734
# NS	ref match    	205824	35.1%
# NS	ref mismatch 	379910	64.9%
# NS	flipped      	206310	35.2%
# NS	swapped      	87213	14.9%
# NS	flip+swap    	86779	14.8%
# NS	unresolved   	4	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# italy_gsa - OK
# ST	A>G	198794	33.9%
# ST	A>T	945	0.2%
# ST	C>A	59417	10.1%
# ST	C>G	1637	0.3%
# ST	C>T	0	0.0%
# ST	G>A	273127	46.6%
# ST	G>C	1777	0.3%
# ST	G>T	0	0.0%
# ST	T>A	1238	0.2%
# ST	T>C	0	0.0%
# ST	T>G	0	0.0%
# # NS, Number of sites:
# NS	total        	586354
# NS	ref match    	249459	42.5%
# NS	ref mismatch 	336895	57.5%
# NS	flipped      	248750	42.4%
# NS	swapped      	44268	7.5%
# NS	flip+swap    	44319	7.6%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# kiel_austria_sibdcs_gsa
# ST	A>G	217399	32.8%
# ST	A>T	1220	0.2%
# ST	C>A	69063	10.4%
# ST	C>G	2201	0.3%
# ST	C>T	0	0.0%
# ST	G>A	314970	47.5%
# ST	G>C	2588	0.4%
# ST	G>T	0	0.0%
# ST	T>A	1661	0.3%
# ST	T>C	0	0.0%
# ST	T>G	0	0.0%
# # NS, Number of sites:
# NS	total        	663705
# NS	ref match    	286853	43.2%
# NS	ref mismatch 	376852	56.8%
# NS	flipped      	286469	43.2%
# NS	swapped      	45332	6.8%
# NS	flip+swap    	45511	6.9%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# netherlands_gsa
# ST	A>G	224706	32.8%
# ST	A>T	1442	0.2%
# ST	C>A	71450	10.4%
# ST	C>G	2768	0.4%
# ST	C>T	0	0.0%
# ST	G>A	322078	47.0%
# ST	G>C	2821	0.4%
# ST	G>T	0	0.0%
# ST	T>A	1908	0.3%
# ST	T>C	0	0.0%
# ST	T>G	0	0.0%
# # NS, Number of sites:
# NS	total        	684826
# NS	ref match    	298101	43.5%
# NS	ref mismatch 	386725	56.5%
# NS	flipped      	297620	43.5%
# NS	swapped      	44622	6.5%
# NS	flip+swap    	44883	6.6%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# slovenia_gsa - OK
# ST	A>G	191549	34.3%
# ST	A>T	687	0.1%
# ST	C>A	56458	10.1%
# ST	C>G	1086	0.2%
# ST	C>T	0	0.0%
# ST	G>A	258043	46.3%
# ST	G>C	1084	0.2%
# ST	G>T	0	0.0%
# ST	T>A	834	0.1%
# ST	T>C	0	0.0%
# ST	T>G	0	0.0%
# # NS, Number of sites:
# NS	total        	557690
# NS	ref match    	235723	42.3%
# NS	ref mismatch 	321967	57.7%
# NS	flipped      	235479	42.2%
# NS	swapped      	43342	7.8%
# NS	flip+swap    	43488	7.8%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# sweden_gsa
# ST	A>G	197840	34.4%
# ST	A>T	764	0.1%
# ST	C>A	57947	10.1%
# ST	C>G	1189	0.2%
# ST	C>T	0	0.0%
# ST	G>A	265857	46.2%
# ST	G>C	1229	0.2%
# ST	G>T	0	0.0%
# ST	T>A	902	0.2%
# ST	T>C	0	0.0%
# ST	T>G	0	0.0%
# # NS, Number of sites:
# NS	total        	575192
# NS	ref match    	244499	42.5%
# NS	ref mismatch 	330693	57.5%
# NS	flipped      	244364	42.5%
# NS	swapped      	43221	7.5%
# NS	flip+swap    	43440	7.6%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# niddk_broad_gsa
# ST	A>G	210759	33.9%
# ST	A>T	942	0.2%
# ST	C>A	63219	10.2%
# ST	C>G	1606	0.3%
# ST	C>T	0	0.0%
# ST	G>A	289763	46.6%
# ST	G>C	1674	0.3%
# ST	G>T	0	0.0%
# ST	T>A	1134	0.2%
# ST	T>C	0	0.0%
# ST	T>G	0	0.0%
# # NS, Number of sites:
# NS	total        	622147
# NS	ref match    	267818	43.0%
# NS	ref mismatch 	354329	57.0%
# NS	flipped      	267749	43.0%
# NS	swapped      	43390	7.0%
# NS	flip+swap    	43547	7.0%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# niddk_feinstein_gsa
# ST	A>G	215705	33.4%
# ST	A>T	1128	0.2%
# ST	C>A	65788	10.2%
# ST	C>G	2200	0.3%
# ST	C>T	0	0.0%
# ST	G>A	302502	46.9%
# ST	G>C	2343	0.4%
# ST	G>T	0	0.0%
# ST	T>A	1441	0.2%
# ST	T>C	0	0.0%
# ST	T>G	0	0.0%
# # NS, Number of sites:
# NS	total        	645538
# NS	ref match    	279650	43.3%
# NS	ref mismatch 	365888	56.7%
# NS	flipped      	279596	43.3%
# NS	swapped      	43205	6.7%
# NS	flip+swap    	43465	6.7%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# basque_gsa
# ST	A>G	197808	34.1%
# ST	A>T	937	0.2%
# ST	C>A	58305	10.1%
# ST	C>G	1609	0.3%
# ST	C>T	0	0.0%
# ST	G>A	268692	46.4%
# ST	G>C	1717	0.3%
# ST	G>T	0	0.0%
# ST	T>A	1172	0.2%
# ST	T>C	0	0.0%
# ST	T>G	0	0.0%
# # NS, Number of sites:
# NS	total        	579360
# NS	ref match    	245994	42.5%
# NS	ref mismatch 	333366	57.5%
# NS	flipped      	245546	42.4%
# NS	swapped      	44159	7.6%
# NS	flip+swap    	44085	7.6%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# lithuania_gsa
# ST	A>G	99502	17.2%
# ST	A>T	892	0.2%
# ST	C>A	29391	5.1%
# ST	C>G	1261	0.2%
# ST	C>T	134348	23.2%
# ST	G>A	133887	23.1%
# ST	G>C	1253	0.2%
# ST	G>T	29085	5.0%
# ST	T>A	866	0.1%
# ST	T>C	99698	17.2%
# ST	T>G	24950	4.3%
# # NS, Number of sites:
# NS	total        	580098
# NS	ref match    	490965	84.6%
# NS	ref mismatch 	89133	15.4%
# NS	flipped      	1040	0.2%
# NS	swapped      	88093	15.2%
# NS	flip+swap    	1144	0.2%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# belgium_louis_gsa
# ST	A>G	100309	16.6%
# ST	A>T	1042	0.2%
# ST	C>A	31555	5.2%
# ST	C>G	1607	0.3%
# ST	C>T	143570	23.7%
# ST	G>A	143194	23.6%
# ST	G>C	1585	0.3%
# ST	G>T	31101	5.1%
# ST	T>A	1018	0.2%
# ST	T>C	100570	16.6%
# ST	T>G	25119	4.1%
# # NS, Number of sites:
# NS	total        	605806
# NS	ref match    	514372	84.9%
# NS	ref mismatch 	91434	15.1%
# NS	flipped      	1229	0.2%
# NS	swapped      	90205	14.9%
# NS	flip+swap    	1443	0.2%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# belgium_franchimont_gsa
# ST	A>G	99519	17.0%
# ST	A>T	918	0.2%
# ST	C>A	29850	5.1%
# ST	C>G	1433	0.2%
# ST	C>T	136169	23.3%
# ST	G>A	135873	23.2%
# ST	G>C	1417	0.2%
# ST	G>T	29477	5.0%
# ST	T>A	956	0.2%
# ST	T>C	99701	17.0%
# ST	T>G	24998	4.3%
# # NS, Number of sites:
# NS	total        	585266
# NS	ref match    	498216	85.1%
# NS	ref mismatch 	87050	14.9%
# NS	flipped      	356	0.1%
# NS	swapped      	86694	14.8%
# NS	flip+swap    	2064	0.4%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# belgium_vermeire_gsa
# ST	A>G	109587	16.2%
# ST	A>T	1606	0.2%
# ST	C>A	35526	5.3%
# ST	C>G	2743	0.4%
# ST	C>T	160061	23.7%
# ST	G>A	159837	23.7%
# ST	G>C	2748	0.4%
# ST	G>T	35044	5.2%
# ST	T>A	1634	0.2%
# ST	T>C	109877	16.3%
# ST	T>G	27996	4.1%
# # NS, Number of sites:
# NS	total        	674795
# NS	ref match    	586861	87.0%
# NS	ref mismatch 	87934	13.0%
# NS	flipped      	378	0.1%
# NS	swapped      	87556	13.0%
# NS	flip+swap    	4047	0.6%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# prism_nfe_gsa
# ST	A>G	98474	17.2%
# ST	A>T	885	0.2%
# ST	C>A	29149	5.1%
# ST	C>G	1282	0.2%
# ST	C>T	132516	23.1%
# ST	G>A	132311	23.1%
# ST	G>C	1260	0.2%
# ST	G>T	28828	5.0%
# ST	T>A	836	0.1%
# ST	T>C	98705	17.2%
# ST	T>G	24744	4.3%
# # NS, Number of sites:
# NS	total        	573760
# NS	ref match    	486133	84.7%
# NS	ref mismatch 	87627	15.3%
# NS	flipped      	1060	0.2%
# NS	swapped      	86567	15.1%
# NS	flip+swap    	1147	0.2%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# prism_nfe_gwas
# ST	A>G	47719	19.6%
# ST	A>T	70	0.0%
# ST	C>A	11746	4.8%
# ST	C>G	111	0.0%
# ST	C>T	51368	21.1%
# ST	G>A	51072	20.9%
# ST	G>C	121	0.0%
# ST	G>T	11842	4.9%
# ST	T>A	80	0.0%
# ST	T>C	47361	19.4%
# ST	T>G	11168	4.6%
# # NS, Number of sites:
# NS	total        	243983
# NS	ref match    	167764	68.8%
# NS	ref mismatch 	76219	31.2%
# NS	flipped      	92	0.0%
# NS	swapped      	76127	31.2%
# NS	flip+swap    	94	0.0%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# finland_illugwas
# ST	A>G	46052	19.2%
# ST	A>T	70	0.0%
# ST	C>A	11734	4.9%
# ST	C>G	112	0.0%
# ST	C>T	51117	21.4%
# ST	G>A	50807	21.2%
# ST	G>C	119	0.0%
# ST	G>T	11867	5.0%
# ST	T>A	78	0.0%
# ST	T>C	45690	19.1%
# ST	T>G	10719	4.5%
# # NS, Number of sites:
# NS	total        	239244
# NS	ref match    	162387	67.9%
# NS	ref mismatch 	76857	32.1%
# NS	flipped      	91	0.0%
# NS	swapped      	76766	32.1%
# NS	flip+swap    	93	0.0%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# german_affy6_old_gwas
# ST	A>G	137268	15.4%
# ST	A>T	28511	3.2%
# ST	C>A	38193	4.3%
# ST	C>G	40636	4.5%
# ST	C>T	170767	19.1%
# ST	G>A	170283	19.1%
# ST	G>C	40579	4.5%
# ST	G>T	38083	4.3%
# ST	T>A	28106	3.1%
# ST	T>C	136927	15.3%
# ST	T>G	31978	3.6%
# # NS, Number of sites:
# NS	total        	893289
# NS	ref match    	646774	72.4%
# NS	ref mismatch 	246514	27.6%
# NS	flipped      	19245	2.2%
# NS	swapped      	227235	25.4%
# NS	flip+swap    	49759	5.6%
# NS	unresolved   	7	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	1
# NS	non-ACGT     	1
# NS	non-SNP      	0
# NS	non-biallelic	0
# norway_affy6_old_gwas
# ST	A>G	130056	15.4%
# ST	A>T	27045	3.2%
# ST	C>A	36285	4.3%
# ST	C>G	38670	4.6%
# ST	C>T	161207	19.0%
# ST	G>A	160960	19.0%
# ST	G>C	38443	4.5%
# ST	G>T	36193	4.3%
# ST	T>A	26718	3.2%
# ST	T>C	129946	15.4%
# ST	T>G	30440	3.6%
# # NS, Number of sites:
# NS	total        	846288
# NS	ref match    	606535	71.7%
# NS	ref mismatch 	239752	28.3%
# NS	flipped      	18685	2.2%
# NS	swapped      	221034	26.1%
# NS	flip+swap    	46810	5.5%
# NS	unresolved   	7	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	1
# NS	non-ACGT     	1
# NS	non-SNP      	0
# NS	non-biallelic	0
# belgium_inf1_old_gwas
# ST	A>G	58072	19.2%
# ST	A>T	0	0.0%
# ST	C>A	14360	4.7%
# ST	C>G	0	0.0%
# ST	C>T	66034	21.8%
# ST	G>A	65893	21.8%
# ST	G>C	0	0.0%
# ST	G>T	14189	4.7%
# ST	T>A	0	0.0%
# ST	T>C	57860	19.1%
# ST	T>G	13136	4.3%
# # NS, Number of sites:
# NS	total        	302700
# NS	ref match    	158926	52.5%
# NS	ref mismatch 	143774	47.5%
# NS	flipped      	42730	14.1%
# NS	swapped      	74315	24.6%
# NS	flip+swap    	26729	8.8%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# belgium_inf2_old_gwas
# ST	A>G	55207	19.0%
# ST	A>T	0	0.0%
# ST	C>A	14155	4.9%
# ST	C>G	0	0.0%
# ST	C>T	63482	21.8%
# ST	G>A	63475	21.8%
# ST	G>C	0	0.0%
# ST	G>T	13912	4.8%
# ST	T>A	0	0.0%
# ST	T>C	54963	18.9%
# ST	T>G	12755	4.4%
# # NS, Number of sites:
# NS	total        	290730
# NS	ref match    	152845	52.6%
# NS	ref mismatch 	137885	47.4%
# NS	flipped      	41263	14.2%
# NS	swapped      	70786	24.3%
# NS	flip+swap    	25836	8.9%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# cedars_370k_old_gwas
# ST	A>G	129007	37.9%
# ST	A>T	306	0.1%
# ST	C>A	32736	9.6%
# ST	C>G	436	0.1%
# ST	C>T	0	0.0%
# ST	G>A	146166	43.0%
# ST	G>C	445	0.1%
# ST	G>T	0	0.0%
# ST	T>A	336	0.1%
# ST	T>C	0	0.0%
# ST	T>G	0	0.0%
# # NS, Number of sites:
# NS	total        	339977
# NS	ref match    	114347	33.6%
# NS	ref mismatch 	225630	66.4%
# NS	flipped      	114396	33.6%
# NS	swapped      	55935	16.5%
# NS	flip+swap    	55462	16.3%
# NS	unresolved   	2	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# cedars_610k_old_gwas
# ST	A>G	220403	37.6%
# ST	A>T	595	0.1%
# ST	C>A	56898	9.7%
# ST	C>G	932	0.2%
# ST	C>T	0	0.0%
# ST	G>A	253767	43.3%
# ST	G>C	987	0.2%
# ST	G>T	0	0.0%
# ST	T>A	650	0.1%
# ST	T>C	0	0.0%
# ST	T>G	0	0.0%
# # NS, Number of sites:
# NS	total        	586380
# NS	ref match    	206587	35.2%
# NS	ref mismatch 	379793	64.8%
# NS	flipped      	207060	35.3%
# NS	swapped      	86784	14.8%
# NS	flip+swap    	86344	14.7%
# NS	unresolved   	4	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# cedars_omni_old_gwas
# ST	A>G	269224	37.2%
# ST	A>T	384	0.1%
# ST	C>A	71044	9.8%
# ST	C>G	668	0.1%
# ST	C>T	0	0.0%
# ST	G>A	316256	43.7%
# ST	G>C	872	0.1%
# ST	G>T	0	0.0%
# ST	T>A	556	0.1%
# ST	T>C	0	0.0%
# ST	T>G	0	0.0%
# # NS, Number of sites:
# NS	total        	723004
# NS	ref match    	247793	34.3%
# NS	ref mismatch 	475211	65.7%
# NS	flipped      	250247	34.6%
# NS	swapped      	111505	15.4%
# NS	flip+swap    	112310	15.5%
# NS	unresolved   	1350	0.2%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# swedish_uc_old_gwas
# ST	A>G	113460	38.2%
# ST	A>T	0	0.0%
# ST	C>A	28411	9.6%
# ST	C>G	0	0.0%
# ST	C>T	0	0.0%
# ST	G>A	128519	43.3%
# ST	G>C	0	0.0%
# ST	G>T	0	0.0%
# ST	T>A	0	0.0%
# ST	T>C	0	0.0%
# ST	T>G	0	0.0%
# # NS, Number of sites:
# NS	total        	297004
# NS	ref match    	98691	33.2%
# NS	ref mismatch 	198313	66.8%
# NS	flipped      	98866	33.3%
# NS	swapped      	49935	16.8%
# NS	flip+swap    	49512	16.7%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# mccauley_gsa
# ST	A>G	97703	17.1%
# ST	A>T	815	0.1%
# ST	C>A	28998	5.1%
# ST	C>G	1288	0.2%
# ST	C>T	132446	23.2%
# ST	G>A	132238	23.2%
# ST	G>C	1270	0.2%
# ST	G>T	28678	5.0%
# ST	T>A	873	0.2%
# ST	T>C	97832	17.1%
# ST	T>G	24496	4.3%
# # NS, Number of sites:
# NS	total        	571173
# NS	ref match    	484679	84.9%
# NS	ref mismatch 	86494	15.1%
# NS	flipped      	358	0.1%
# NS	swapped      	86136	15.1%
# NS	flip+swap    	1844	0.3%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# ccfa_gsa
# ST	A>G	102986	16.8%
# ST	A>T	1026	0.2%
# ST	C>A	31427	5.1%
# ST	C>G	1603	0.3%
# ST	C>T	142854	23.4%
# ST	G>A	142650	23.3%
# ST	G>C	1584	0.3%
# ST	G>T	30918	5.1%
# ST	T>A	1064	0.2%
# ST	T>C	103311	16.9%
# ST	T>G	25951	4.2%
# # NS, Number of sites:
# NS	total        	611354
# NS	ref match    	525069	85.9%
# NS	ref mismatch 	86285	14.1%
# NS	flipped      	358	0.1%
# NS	swapped      	85927	14.1%
# NS	flip+swap    	2331	0.4%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# cedars_gsa
# ST	A>G	108335	16.2%
# ST	A>T	1422	0.2%
# ST	C>A	35451	5.3%
# ST	C>G	2422	0.4%
# ST	C>T	159904	23.9%
# ST	G>A	159606	23.8%
# ST	G>C	2444	0.4%
# ST	G>T	35010	5.2%
# ST	T>A	1464	0.2%
# ST	T>C	108624	16.2%
# ST	T>G	27297	4.1%
# # NS, Number of sites:
# NS	total        	669422
# NS	ref match    	580282	86.7%
# NS	ref mismatch 	89140	13.3%
# NS	flipped      	381	0.1%
# NS	swapped      	88759	13.3%
# NS	flip+swap    	3537	0.5%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# bernstein_gsa
# ST	A>G	97240	17.2%
# ST	A>T	801	0.1%
# ST	C>A	28727	5.1%
# ST	C>G	1174	0.2%
# ST	C>T	130848	23.1%
# ST	G>A	130433	23.1%
# ST	G>C	1134	0.2%
# ST	G>T	28313	5.0%
# ST	T>A	815	0.1%
# ST	T>C	97505	17.2%
# ST	T>G	24358	4.3%
# # NS, Number of sites:
# NS	total        	565758
# NS	ref match    	479251	84.7%
# NS	ref mismatch 	86507	15.3%
# NS	flipped      	343	0.1%
# NS	swapped      	86164	15.2%
# NS	flip+swap    	1688	0.3%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# farkkila_gsa
# ST	A>G	92044	17.4%
# ST	A>T	679	0.1%
# ST	C>A	26694	5.1%
# ST	C>G	947	0.2%
# ST	C>T	121031	22.9%
# ST	G>A	120679	22.8%
# ST	G>C	910	0.2%
# ST	G>T	26311	5.0%
# ST	T>A	673	0.1%
# ST	T>C	92342	17.5%
# ST	T>G	23039	4.4%
# # NS, Number of sites:
# NS	total        	528311
# NS	ref match    	440617	83.4%
# NS	ref mismatch 	87694	16.6%
# NS	flipped      	790	0.1%
# NS	swapped      	86904	16.4%
# NS	flip+swap    	848	0.2%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# franchimont_gsa
# ST	A>G	107553	16.0%
# ST	A>T	1594	0.2%
# ST	C>A	35717	5.3%
# ST	C>G	2745	0.4%
# ST	C>T	160924	24.0%
# ST	G>A	160651	23.9%
# ST	G>C	2754	0.4%
# ST	G>T	35243	5.3%
# ST	T>A	1659	0.2%
# ST	T>C	107780	16.1%
# ST	T>G	27182	4.1%
# # NS, Number of sites:
# NS	total        	671132
# NS	ref match    	581524	86.6%
# NS	ref mismatch 	89608	13.4%
# NS	flipped      	401	0.1%
# NS	swapped      	89207	13.3%
# NS	flip+swap    	4010	0.6%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# franke_gsa
# ST	A>G	98085	17.2%
# ST	A>T	810	0.1%
# ST	C>A	28953	5.1%
# ST	C>G	1220	0.2%
# ST	C>T	132080	23.1%
# ST	G>A	131837	23.1%
# ST	G>C	1192	0.2%
# ST	G>T	28549	5.0%
# ST	T>A	845	0.1%
# ST	T>C	98354	17.2%
# ST	T>G	24606	4.3%
# # NS, Number of sites:
# NS	total        	571118
# NS	ref match    	484572	84.8%
# NS	ref mismatch 	86546	15.2%
# NS	flipped      	342	0.1%
# NS	swapped      	86204	15.1%
# NS	flip+swap    	1752	0.3%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# helmsley_prism_gsa
# ST	A>G	47335	19.4%
# ST	A>T	71	0.0%
# ST	C>A	11833	4.8%
# ST	C>G	111	0.0%
# ST	C>T	51767	21.2%
# ST	G>A	51456	21.1%
# ST	G>C	121	0.0%
# ST	G>T	11947	4.9%
# ST	T>A	79	0.0%
# ST	T>C	46962	19.2%
# ST	T>G	11063	4.5%
# # NS, Number of sites:
# NS	total        	243983
# NS	ref match    	167832	68.8%
# NS	ref mismatch 	76151	31.2%
# NS	flipped      	91	0.0%
# NS	swapped      	76060	31.2%
# NS	flip+swap    	95	0.0%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# helmsley_xavier_prism_gsa
# ST	A>G	46929	19.2%
# ST	A>T	70	0.0%
# ST	C>A	11951	4.9%
# ST	C>G	114	0.0%
# ST	C>T	52238	21.4%
# ST	G>A	51853	21.3%
# ST	G>C	118	0.0%
# ST	G>T	12041	4.9%
# ST	T>A	80	0.0%
# ST	T>C	46478	19.1%
# ST	T>G	10962	4.5%
# # NS, Number of sites:
# NS	total        	243952
# NS	ref match    	167809	68.8%
# NS	ref mismatch 	76143	31.2%
# NS	flipped      	92	0.0%
# NS	swapped      	76051	31.2%
# NS	flip+swap    	94	0.0%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# hyams_protect_gsa
# ST	A>G	99174	17.1%
# ST	A>T	935	0.2%
# ST	C>A	29602	5.1%
# ST	C>G	1293	0.2%
# ST	C>T	134841	23.2%
# ST	G>A	134464	23.1%
# ST	G>C	1279	0.2%
# ST	G>T	29208	5.0%
# ST	T>A	971	0.2%
# ST	T>C	99257	17.1%
# ST	T>G	24901	4.3%
# # NS, Number of sites:
# NS	total        	580928
# NS	ref match    	495253	85.3%
# NS	ref mismatch 	85675	14.7%
# NS	flipped      	357	0.1%
# NS	swapped      	85318	14.7%
# NS	flip+swap    	1923	0.3%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# lewis_sparc_gsa
# ST	A>G	103835	16.8%
# ST	A>T	1129	0.2%
# ST	C>A	31930	5.2%
# ST	C>G	1747	0.3%
# ST	C>T	145258	23.4%
# ST	G>A	145080	23.4%
# ST	G>C	1745	0.3%
# ST	G>T	31393	5.1%
# ST	T>A	1147	0.2%
# ST	T>C	104158	16.8%
# ST	T>G	26211	4.2%
# # NS, Number of sites:
# NS	total        	619895
# NS	ref match    	533530	86.1%
# NS	ref mismatch 	86365	13.9%
# NS	flipped      	360	0.1%
# NS	swapped      	86005	13.9%
# NS	flip+swap    	2578	0.4%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# mccauley_new_gsa
# ST	A>G	100026	16.9%
# ST	A>T	943	0.2%
# ST	C>A	30121	5.1%
# ST	C>G	1525	0.3%
# ST	C>T	137947	23.3%
# ST	G>A	137904	23.3%
# ST	G>C	1536	0.3%
# ST	G>T	29842	5.0%
# ST	T>A	988	0.2%
# ST	T>C	100197	16.9%
# ST	T>G	25116	4.2%
# # NS, Number of sites:
# NS	total        	591351
# NS	ref match    	505089	85.4%
# NS	ref mismatch 	86262	14.6%
# NS	flipped      	356	0.1%
# NS	swapped      	85906	14.5%
# NS	flip+swap    	2207	0.4%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# mcgovern_gsa
# ST	A>G	108778	16.7%
# ST	A>T	1230	0.2%
# ST	C>A	33731	5.2%
# ST	C>G	1994	0.3%
# ST	C>T	153024	23.5%
# ST	G>A	152711	23.4%
# ST	G>C	1978	0.3%
# ST	G>T	33264	5.1%
# ST	T>A	1237	0.2%
# ST	T>C	109054	16.7%
# ST	T>G	27429	4.2%
# # NS, Number of sites:
# NS	total        	652046
# NS	ref match    	563874	86.5%
# NS	ref mismatch 	88172	13.5%
# NS	flipped      	374	0.1%
# NS	swapped      	87798	13.5%
# NS	flip+swap    	2891	0.4%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# moayyedi_imagine_gsa
# ST	A>G	103918	16.9%
# ST	A>T	992	0.2%
# ST	C>A	31577	5.1%
# ST	C>G	1577	0.3%
# ST	C>T	143755	23.4%
# ST	G>A	143404	23.3%
# ST	G>C	1529	0.2%
# ST	G>T	31280	5.1%
# ST	T>A	1005	0.2%
# ST	T>C	104192	16.9%
# ST	T>G	26098	4.2%
# # NS, Number of sites:
# NS	total        	615521
# NS	ref match    	528936	85.9%
# NS	ref mismatch 	86585	14.1%
# NS	flipped      	350	0.1%
# NS	swapped      	86235	14.0%
# NS	flip+swap    	2272	0.4%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# newberry_share_gsa
# ST	A>G	99672	16.8%
# ST	A>T	924	0.2%
# ST	C>A	30494	5.1%
# ST	C>G	1425	0.2%
# ST	C>T	138993	23.5%
# ST	G>A	138626	23.4%
# ST	G>C	1418	0.2%
# ST	G>T	30187	5.1%
# ST	T>A	934	0.2%
# ST	T>C	99912	16.9%
# ST	T>G	25016	4.2%
# # NS, Number of sites:
# NS	total        	592708
# NS	ref match    	506608	85.5%
# NS	ref mismatch 	86100	14.5%
# NS	flipped      	346	0.1%
# NS	swapped      	85754	14.5%
# NS	flip+swap    	2060	0.3%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# niddk_cho_gsa
# ST	A>G	102857	16.9%
# ST	A>T	1027	0.2%
# ST	C>A	31124	5.1%
# ST	C>G	1587	0.3%
# ST	C>T	141966	23.3%
# ST	G>A	141667	23.3%
# ST	G>C	1553	0.3%
# ST	G>T	30879	5.1%
# ST	T>A	1031	0.2%
# ST	T>C	103141	16.9%
# ST	T>G	25920	4.3%
# # NS, Number of sites:
# NS	total        	608759
# NS	ref match    	522000	85.7%
# NS	ref mismatch 	86759	14.3%
# NS	flipped      	366	0.1%
# NS	swapped      	86393	14.2%
# NS	flip+swap    	2292	0.4%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# niddk_duerr_gsa
# ST	A>G	101205	17.0%
# ST	A>T	990	0.2%
# ST	C>A	30411	5.1%
# ST	C>G	1562	0.3%
# ST	C>T	138718	23.3%
# ST	G>A	138441	23.2%
# ST	G>C	1525	0.3%
# ST	G>T	30019	5.0%
# ST	T>A	1038	0.2%
# ST	T>C	101430	17.0%
# ST	T>G	25480	4.3%
# # NS, Number of sites:
# NS	total        	596287
# NS	ref match    	509346	85.4%
# NS	ref mismatch 	86941	14.6%
# NS	flipped      	355	0.1%
# NS	swapped      	86586	14.5%
# NS	flip+swap    	2273	0.4%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# niddk_rioux_gsa
# ST	A>G	100631	17.1%
# ST	A>T	949	0.2%
# ST	C>A	29976	5.1%
# ST	C>G	1394	0.2%
# ST	C>T	136593	23.2%
# ST	G>A	136524	23.2%
# ST	G>C	1406	0.2%
# ST	G>T	29638	5.0%
# ST	T>A	939	0.2%
# ST	T>C	100852	17.1%
# ST	T>G	25355	4.3%
# # NS, Number of sites:
# NS	total        	589690
# NS	ref match    	503027	85.3%
# NS	ref mismatch 	86663	14.7%
# NS	flipped      	350	0.1%
# NS	swapped      	86313	14.6%
# NS	flip+swap    	2053	0.3%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# niddk_silverberg_gsa
# ST	A>G	105242	16.8%
# ST	A>T	1090	0.2%
# ST	C>A	32213	5.1%
# ST	C>G	1696	0.3%
# ST	C>T	146838	23.4%
# ST	G>A	146539	23.4%
# ST	G>C	1699	0.3%
# ST	G>T	31877	5.1%
# ST	T>A	1098	0.2%
# ST	T>C	105419	16.8%
# ST	T>G	26531	4.2%
# # NS, Number of sites:
# NS	total        	626861
# NS	ref match    	539771	86.1%
# NS	ref mismatch 	87090	13.9%
# NS	flipped      	361	0.1%
# NS	swapped      	86729	13.8%
# NS	flip+swap    	2484	0.4%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# palotie_hus_gsa
# ST	A>G	97162	17.2%
# ST	A>T	818	0.1%
# ST	C>A	28663	5.1%
# ST	C>G	1177	0.2%
# ST	C>T	130775	23.1%
# ST	G>A	130411	23.1%
# ST	G>C	1149	0.2%
# ST	G>T	28327	5.0%
# ST	T>A	813	0.1%
# ST	T>C	97403	17.2%
# ST	T>G	24333	4.3%
# # NS, Number of sites:
# NS	total        	565416
# NS	ref match    	478295	84.6%
# NS	ref mismatch 	87121	15.4%
# NS	flipped      	347	0.1%
# NS	swapped      	86774	15.3%
# NS	flip+swap    	1673	0.3%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# pekow_share_gsa
# ST	A>G	99101	17.0%
# ST	A>T	886	0.2%
# ST	C>A	29658	5.1%
# ST	C>G	1328	0.2%
# ST	C>T	135030	23.2%
# ST	G>A	134682	23.2%
# ST	G>C	1330	0.2%
# ST	G>T	29251	5.0%
# ST	T>A	891	0.2%
# ST	T>C	99335	17.1%
# ST	T>G	24839	4.3%
# # NS, Number of sites:
# NS	total        	581240
# NS	ref match    	495532	85.3%
# NS	ref mismatch 	85708	14.7%
# NS	flipped      	345	0.1%
# NS	swapped      	85363	14.7%
# NS	flip+swap    	1946	0.3%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# rioux_igenomed_gsa
# ST	A>G	95234	17.3%
# ST	A>T	743	0.1%
# ST	C>A	27932	5.1%
# ST	C>G	1059	0.2%
# ST	C>T	127318	23.1%
# ST	G>A	126999	23.0%
# ST	G>C	1038	0.2%
# ST	G>T	27590	5.0%
# ST	T>A	771	0.1%
# ST	T>C	95452	17.3%
# ST	T>G	23795	4.3%
# # NS, Number of sites:
# NS	total        	551803
# NS	ref match    	465713	84.4%
# NS	ref mismatch 	86090	15.6%
# NS	flipped      	337	0.1%
# NS	swapped      	85753	15.5%
# NS	flip+swap    	1515	0.3%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# sands_msccr_gsa
# ST	A>G	102340	16.9%
# ST	A>T	940	0.2%
# ST	C>A	30898	5.1%
# ST	C>G	1463	0.2%
# ST	C>T	140851	23.3%
# ST	G>A	140612	23.3%
# ST	G>C	1498	0.2%
# ST	G>T	30554	5.1%
# ST	T>A	970	0.2%
# ST	T>C	102557	17.0%
# ST	T>G	25696	4.3%
# # NS, Number of sites:
# NS	total        	604196
# NS	ref match    	518318	85.8%
# NS	ref mismatch 	85878	14.2%
# NS	flipped      	364	0.1%
# NS	swapped      	85514	14.2%
# NS	flip+swap    	2130	0.4%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# stampfer_gsa
# ST	A>G	108187	16.1%
# ST	A>T	1598	0.2%
# ST	C>A	35720	5.3%
# ST	C>G	2754	0.4%
# ST	C>T	160722	23.9%
# ST	G>A	160520	23.9%
# ST	G>C	2772	0.4%
# ST	G>T	35219	5.2%
# ST	T>A	1649	0.2%
# ST	T>C	108535	16.1%
# ST	T>G	27432	4.1%
# # NS, Number of sites:
# NS	total        	672634
# NS	ref match    	582190	86.6%
# NS	ref mismatch 	90444	13.4%
# NS	flipped      	383	0.1%
# NS	swapped      	90061	13.4%
# NS	flip+swap    	4047	0.6%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# vermeire_gsa
# ST	A>G	109756	16.3%
# ST	A>T	1607	0.2%
# ST	C>A	35548	5.3%
# ST	C>G	2745	0.4%
# ST	C>T	160123	23.7%
# ST	G>A	159897	23.7%
# ST	G>C	2750	0.4%
# ST	G>T	35072	5.2%
# ST	T>A	1634	0.2%
# ST	T>C	110031	16.3%
# ST	T>G	28032	4.2%
# # NS, Number of sites:
# NS	total        	675377
# NS	ref match    	587367	87.0%
# NS	ref mismatch 	88010	13.0%
# NS	flipped      	381	0.1%
# NS	swapped      	87629	13.0%
# NS	flip+swap    	4046	0.6%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# weersma_gsa
# ST	A>G	98521	17.2%
# ST	A>T	825	0.1%
# ST	C>A	29248	5.1%
# ST	C>G	1222	0.2%
# ST	C>T	132731	23.1%
# ST	G>A	132377	23.1%
# ST	G>C	1229	0.2%
# ST	G>T	28807	5.0%
# ST	T>A	848	0.1%
# ST	T>C	98713	17.2%
# ST	T>G	24669	4.3%
# # NS, Number of sites:
# NS	total        	573894
# NS	ref match    	487583	85.0%
# NS	ref mismatch 	86311	15.0%
# NS	flipped      	334	0.1%
# NS	swapped      	85977	15.0%
# NS	flip+swap    	1793	0.3%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# xavier_prism_gsa
# ST	A>G	99907	17.1%
# ST	A>T	868	0.1%
# ST	C>A	29763	5.1%
# ST	C>G	1337	0.2%
# ST	C>T	135333	23.2%
# ST	G>A	135041	23.1%
# ST	G>C	1354	0.2%
# ST	G>T	29409	5.0%
# ST	T>A	912	0.2%
# ST	T>C	100113	17.1%
# ST	T>G	25163	4.3%
# # NS, Number of sites:
# NS	total        	584364
# NS	ref match    	497968	85.2%
# NS	ref mismatch 	86396	14.8%
# NS	flipped      	355	0.1%
# NS	swapped      	86041	14.7%
# NS	flip+swap    	1965	0.3%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0
# xavier_share_gsa
# ST	A>G	99702	17.0%
# ST	A>T	894	0.2%
# ST	C>A	29800	5.1%
# ST	C>G	1336	0.2%
# ST	C>T	135857	23.2%
# ST	G>A	135472	23.2%
# ST	G>C	1395	0.2%
# ST	G>T	29518	5.0%
# ST	T>A	915	0.2%
# ST	T>C	100004	17.1%
# ST	T>G	25029	4.3%
# # NS, Number of sites:
# NS	total        	584949
# NS	ref match    	498547	85.2%
# NS	ref mismatch 	86402	14.8%
# NS	flipped      	356	0.1%
# NS	swapped      	86046	14.7%
# NS	flip+swap    	1977	0.3%
# NS	unresolved   	0	0.0%
# NS	fixed pos    	0	0.0%
# NS	skipped      	0
# NS	non-ACGT     	0
# NS	non-SNP      	0
# NS	non-biallelic	0

#################################

# REMOVE MISSMATCH ALLELES, THEY  IMPUTATION SERVER TO CRASH - just for a subset of old studies

# studies=(australia_omniexome gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa)
# studies=(niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas)
# studies=(norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas xavier_prism_gsa xavier_share_gsa)
# studies=(cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa)
# studies=(helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa)
# studies=(niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa)

for i in ${studies[@]}
do
bsub -J"bcf1.2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_bcftools_1.2_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_bcftools_1.2_${i} \
"/path/to/software/bcftools-1.16/./bcftools norm \
--check-ref x -f ${path_gwas}resources/1000GP/human_g1k_v37.fasta \
${path_gwas}pre_imputation/QC/${i}/${i}_hg19_posstrandaligned.vcf.gz \
-Oz -o ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_posstrandaligned_2.vcf.gz"
done

# gwas2 - old error output
# [W::bcf_sr_add_reader] No BGZF EOF marker; file '/path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas2/gwas2_hg19_posstrandaligned.vcf.gz' may be truncated

for i in ${studies[@]}
do
echo ${i} && tail ${path_gwas}pre_imputation/QC/${i}/logs/stderr_bcftools_1.2_${i}
done

for i in ${studies[@]}
do
echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_posstrandaligned_2.vcf.gz
done

################
# australia_omniexome
# Lines   total/split/realigned/skipped:	807107/0/0/0
# gwas1
# Lines   total/split/realigned/skipped:	459426/0/0/1
# gwas2
# Lines   total/split/realigned/skipped:	799264/0/0/1
# pittsburgh_gsa
# Lines   total/split/realigned/skipped:	945375/0/0/83
# spain_gsa
# Lines   total/split/realigned/skipped:	585734/0/0/4
# italy_gsa
# Lines   total/split/realigned/skipped:	586354/0/0/0
# kiel_austria_sibdcs_gsa
# Lines   total/split/realigned/skipped:	663705/0/0/0
# netherlands_gsa
# Lines   total/split/realigned/skipped:	684826/0/0/0
# slovenia_gsa
# Lines   total/split/realigned/skipped:	557690/0/0/0
# sweden_gsa
# Lines   total/split/realigned/skipped:	575192/0/0/0
# niddk_broad_gsa
# Lines   total/split/realigned/skipped:	622147/0/0/0
# niddk_feinstein_gsa
# Lines   total/split/realigned/skipped:	645538/0/0/0
# basque_gsa
# Lines   total/split/realigned/skipped:	579360/0/0/0
# lithuania_gsa
# Lines   total/split/realigned/skipped:	580098/0/0/0
# belgium_louis_gsa
# Lines   total/split/realigned/skipped:	605806/0/0/0
# belgium_franchimont_gsa
# Lines   total/split/realigned/skipped:	585266/0/0/0
# belgium_vermeire_gsa
# Lines   total/split/realigned/skipped:	674795/0/0/0
# prism_nfe_gsa
# Lines   total/split/realigned/skipped:	573760/0/0/0
# prism_nfe_gwas
# Lines   total/split/realigned/skipped:	243983/0/0/0
# finland_illugwas
# Lines   total/split/realigned/skipped:	239244/0/0/0
# german_affy6_old_gwas
# Lines   total/split/realigned/skipped:	893289/0/0/8
# norway_affy6_old_gwas
# Lines   total/split/realigned/skipped:	846288/0/0/8
# belgium_inf1_old_gwas
# Lines   total/split/realigned/skipped:	302700/0/0/0
# belgium_inf2_old_gwas
# Lines   total/split/realigned/skipped:	290730/0/0/0
# cedars_370k_old_gwas
# Lines   total/split/realigned/skipped:	339977/0/0/2
# cedars_610k_old_gwas
# Lines   total/split/realigned/skipped:	586380/0/0/4
# cedars_omni_old_gwas
# Lines   total/split/realigned/skipped:	723004/0/0/1350
# swedish_uc_old_gwas
# Lines   total/split/realigned/skipped:	297004/0/0/0
# mccauley_gsa
# Lines   total/split/realigned/skipped:	571173/0/0/0
# ccfa_gsa
# Lines   total/split/realigned/skipped:	611354/0/0/0
# cedars_gsa
# Lines   total/split/realigned/skipped:	669422/0/0/0
# bernstein_gsa
# Lines   total/split/realigned/skipped:	565758/0/0/0
# farkkila_gsa
# Lines   total/split/realigned/skipped:	528311/0/0/0
# franchimont_gsa
# Lines   total/split/realigned/skipped:	671132/0/0/0
# franke_gsa
# Lines   total/split/realigned/skipped:	571118/0/0/0
# helmsley_prism_gsa
# Lines   total/split/realigned/skipped:	243983/0/0/0
# helmsley_xavier_prism_gsa
# Lines   total/split/realigned/skipped:	243952/0/0/0
# hyams_protect_gsa
# Lines   total/split/realigned/skipped:	580928/0/0/0
# lewis_sparc_gsa
# Lines   total/split/realigned/skipped:	619895/0/0/0
# mccauley_new_gsa
# Lines   total/split/realigned/skipped:	591351/0/0/0
# mcgovern_gsa
# Lines   total/split/realigned/skipped:	652046/0/0/0
# moayyedi_imagine_gsa
# Lines   total/split/realigned/skipped:	615521/0/0/0
# newberry_share_gsa
# Lines   total/split/realigned/skipped:	592708/0/0/0
# niddk_cho_gsa
# Lines   total/split/realigned/skipped:	608759/0/0/0
# niddk_duerr_gsa
# Lines   total/split/realigned/skipped:	596287/0/0/0
# niddk_rioux_gsa
# Lines   total/split/realigned/skipped:	589690/0/0/0
# niddk_silverberg_gsa
# Lines   total/split/realigned/skipped:	626861/0/0/0
# palotie_hus_gsa
# Lines   total/split/realigned/skipped:	565416/0/0/0
# pekow_share_gsa
# Lines   total/split/realigned/skipped:	581240/0/0/0
# rioux_igenomed_gsa
# Lines   total/split/realigned/skipped:	551803/0/0/0
# sands_msccr_gsa
# Lines   total/split/realigned/skipped:	604196/0/0/0
# stampfer_gsa
# Lines   total/split/realigned/skipped:	672634/0/0/0
# vermeire_gsa
# Lines   total/split/realigned/skipped:	675377/0/0/0
# weersma_gsa
# Lines   total/split/realigned/skipped:	573894/0/0/0
# xavier_prism_gsa
# Lines   total/split/realigned/skipped:	584364/0/0/0
# xavier_share_gsa
# Lines   total/split/realigned/skipped:	584949/0/0/0

################

###########################################################
# 4.3 DOUBLE-CHECK THE CONVERSION
gwas2 mcgovern_gsa

studies=(australia_omniexome pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa)
studies=(niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas)
studies=(norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas xavier_prism_gsa xavier_share_gsa)
studies=(cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa)
studies=(helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa  moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa)
studies=(niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa)


for i in ${studies[@]}
do
bsub -J"bcf2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_bcftools_2_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_bcftools_2_${i} \
"/path/to/software/bcftools-1.16/./bcftools +fixref \
${path_gwas}pre_imputation/QC/${i}/${i}_hg19_posstrandaligned_2.vcf.gz \
-- -f ${path_gwas}resources/1000GP/human_g1k_v37.fasta"
done


for i in ${studies[@]}
do
echo ${i} && tail -23 ${path_gwas}pre_imputation/QC/${i}/logs/stderr_bcftools_2_${i} | grep "ref match"
done
#######
# australia_omniexome
# NS	ref match    	807107	100.0%
# gwas1
# NS	ref match    	459425	100.0%
# gwas2
# NS	ref match    	799263	100.0%
# pittsburgh_gsa
# NS	ref match    	945292	100.0%
# spain_gsa
# NS	ref match    	585730	100.0%
# italy_gsa
# NS	ref match    	586354	100.0%
# kiel_austria_sibdcs_gsa
# NS	ref match    	663705	100.0%
# netherlands_gsa
# NS	ref match    	684826	100.0%
# slovenia_gsa
# NS	ref match    	557690	100.0%
# sweden_gsa
# NS	ref match    	575192	100.0%
# niddk_broad_gsa
# NS	ref match    	622147	100.0%
# niddk_feinstein_gsa
# NS	ref match    	645538	100.0%
# basque_gsa
# NS	ref match    	579360	100.0%
# lithuania_gsa
# NS	ref match    	580098	100.0%
# belgium_louis_gsa
# NS	ref match    	605806	100.0%
# belgium_franchimont_gsa
# NS	ref match    	585266	100.0%
# belgium_vermeire_gsa
# NS	ref match    	674795	100.0%
# prism_nfe_gsa
# NS	ref match    	573760	100.0%
# prism_nfe_gwas
# NS	ref match    	243983	100.0%
# finland_illugwas
# NS	ref match    	239244	100.0%
# german_affy6_old_gwas
# NS	ref match    	893281	100.0%
# norway_affy6_old_gwas
# NS	ref match    	846280	100.0%
# belgium_inf1_old_gwas
# NS	ref match    	302700	100.0%
# belgium_inf2_old_gwas
# NS	ref match    	290730	100.0%
# cedars_370k_old_gwas
# NS	ref match    	339975	100.0%
# cedars_610k_old_gwas
# NS	ref match    	586376	100.0%
# cedars_omni_old_gwas
# NS	ref match    	721654	100.0%
# swedish_uc_old_gwas
# NS	ref match    	297004	100.0%
# mccauley_gsa
# NS	ref match    	571173	100.0%
# ccfa_gsa
# NS	ref match    	611354	100.0%
# cedars_gsa
# NS	ref match    	669422	100.0%
# bernstein_gsa
# NS	ref match    	565758	100.0%
# farkkila_gsa
# NS	ref match    	528311	100.0%
# franchimont_gsa
# NS	ref match    	671132	100.0%
# franke_gsa
# NS	ref match    	571118	100.0%
# helmsley_prism_gsa
# NS	ref match    	243983	100.0%
# helmsley_xavier_prism_gsa
# NS	ref match    	243952	100.0%
# hyams_protect_gsa
# NS	ref match    	580928	100.0%
# lewis_sparc_gsa
# NS	ref match    	619895	100.0%
# mccauley_new_gsa
# NS	ref match    	591351	100.0%
# mcgovern_gsa
# NS	ref match    	652046	100.0%
# moayyedi_imagine_gsa
# NS	ref match    	615521	100.0%
# newberry_share_gsa
# NS	ref match    	592708	100.0%
# niddk_cho_gsa
# NS	ref match    	608759	100.0%
# niddk_duerr_gsa
# NS	ref match    	596287	100.0%
# niddk_rioux_gsa
# NS	ref match    	589690	100.0%
# niddk_silverberg_gsa
# NS	ref match    	626861	100.0%
# palotie_hus_gsa
# NS	ref match    	565416	100.0%
# pekow_share_gsa
# NS	ref match    	581240	100.0%
# rioux_igenomed_gsa
# NS	ref match    	551803	100.0%
# sands_msccr_gsa
# NS	ref match    	604196	100.0%
# stampfer_gsa
# NS	ref match    	672634	100.0%
# vermeire_gsa
# NS	ref match    	675377	100.0%
# weersma_gsa
# NS	ref match    	573894	100.0%
# xavier_prism_gsa
# NS	ref match    	584364	100.0%
# xavier_share_gsa
# NS	ref match    	584949	100.0%

#######

for i in ${studies[@]}
do ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_posstrandaligned_2.vcf.gz
done

# rm intermed files -
for i in ${studies[@]}
do ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_posstrandaligned.vcf.gz
done

for i in ${studies[@]}
do rm ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_posstrandaligned.vcf.gz
done

###########################################################
# 4.4 VCF to BED

MEM=500
for i in ${studies[@]}
do
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_3_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_3_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--vcf ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_posstrandaligned_2.vcf.gz \
--double-id \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr"
done

for i in ${studies[@]}
do
echo ${i} && tail -200 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_3_${i} | grep "Successfully completed."
done

for i in ${studies[@]}
do
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_4_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_4_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr \
--missing \
--out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr"
done

for i in ${studies[@]}
do
echo ${i} && tail -200 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_4_${i} | grep "Successfully completed."
done

for i in ${studies[@]}
do
echo ${i} && ls -la ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr.lmiss
done

