# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
##################################################################################################################################################################

##############################
# 23.- IMPUTATION SERVER QC  #
##############################


########################
# 23.1 Submit QC files #
########################

# create output folders

path_gwas=/path/to/ibdgwas/IIBDGC/
ancestry=(eur noneur)
  
i=
for j in ${ancestry[@]}
do echo ${i} && ls -la ${path_gwas}imputation_ready/${i}/${i}_${j}_chr*_2022.vcf.gz |  wc -l
done
for j in ${ancestry[@]}
do echo ${j} && python2 /path/to/user/scripts/IIBDGC/topmed_qc_${j}.py ${i}
done



##################

# i=swedish_uc_old_gwas
# # eur
# # swedish_uc_old_gwas
# 
# # QC  swedish_uc_old_gwas job_eur=job-20221028-165216-225
# # noneur
# # swedish_uc_old_gwas
# 
# # QC  swedish_uc_old_gwas job_noneur=job-20221028-165220-101
# job_eur=job-20221028-165216-225
# job_noneur=job-20221028-165220-101

# i=xavier_share_gsa
# # eur
# # xavier_share_gsa
# # QC  xavier_share_gsa job_eur=job-20221028-153042-630
# # noneur
# # xavier_share_gsa
# # QC  xavier_share_gsa job_noneur=job-20221028-153045-848
# job_eur=job-20221028-153042-630
# job_noneur=job-20221028-153045-848


# i=vermeire_gsa
# # eur
# # vermeire_gsa
# # QC  vermeire_gsa job_eur=job-20221028-152427-062
# # noneur
# # vermeire_gsa
# # QC  vermeire_gsa job_noneur=job-20221028-152436-985
# job_eur=job-20221028-152427-062
# job_noneur=job-20221028-152436-985

# i=xavier_prism_gsa
# # eur
# # xavier_prism_gsa
# # QC  xavier_prism_gsa job_eur=job-20221028-152912-241
# # noneur
# # xavier_prism_gsa
# # QC  xavier_prism_gsa job_noneur=job-20221028-152915-151
# job_eur=job-20221028-152912-241
# job_noneur=job-20221028-152915-151

# i=weersma_gsa
# # eur
# # weersma_gsa
# # QC  weersma_gsa job_eur=job-20221028-152552-423
# # noneur
# # weersma_gsa
# # QC  weersma_gsa job_noneur=job-20221028-152812-977
# job_eur=job-20221028-152552-423
# job_noneur=job-20221028-152812-977

# i=stampfer_gsa
# # eur
# # stampfer_gsa
# # QC  stampfer_gsa job_eur=job-20221028-152313-025
# # noneur
# # stampfer_gsa
# # QC  stampfer_gsa job_noneur=job-20221028-152317-332
# job_eur=job-20221028-152313-025
# job_noneur=job-20221028-152317-332

# i=sands_msccr_gsa
# # eur
# # sands_msccr_gsa
# # QC  sands_msccr_gsa job_eur=job-20221028-152122-714
# # noneur
# # sands_msccr_gsa
# # QC  sands_msccr_gsa job_noneur=job-20221028-152127-428
# job_eur=job-20221028-152122-714
# job_noneur=job-20221028-152127-428

# i=rioux_igenomed_gsa
# # eur
# # rioux_igenomed_gsa
# # QC  rioux_igenomed_gsa job_eur=job-20221028-151914-530
# # noneur
# # rioux_igenomed_gsa
# # QC  rioux_igenomed_gsa job_noneur=job-20221028-151916-543
# job_eur=job-20221028-151914-530
# job_noneur=job-20221028-151916-543

# i=pekow_share_gsa
# # pekow_share_gsa
# # QC  pekow_share_gsa job_eur=job-20221028-151536-142
# # noneur
# # pekow_share_gsa
# # QC  pekow_share_gsa job_noneur=job-20221028-151835-407
# job_eur=job-20221028-151536-142
# job_noneur=job-20221028-151835-407


# i=palotie_hus_gsa
# # eur
# # palotie_hus_gsa
# # QC  palotie_hus_gsa job_eur=job-20221028-151408-050
# # noneur
# # palotie_hus_gsa
# # QC  palotie_hus_gsa job_noneur=job-20221028-151411-924
# job_eur=job-20221028-151408-050
# job_noneur=job-20221028-151411-924


# i=niddk_silverberg_gsa
# # eur
# # niddk_silverberg_gsa
# # QC  niddk_silverberg_gsa job_eur=job-20221028-151249-933
# # noneur
# # niddk_silverberg_gsa
# # QC  niddk_silverberg_gsa job_noneur=job-20221028-151256-171
# job_eur=job-20221028-151249-933
# job_noneur=job-20221028-151256-171

# i=niddk_rioux_gsa
# # eur
# # niddk_rioux_gsa
# # QC  niddk_rioux_gsa job_eur=job-20221028-151132-249
# # noneur
# # niddk_rioux_gsa
# # QC  niddk_rioux_gsa job_noneur=job-20221028-151136-545
# job_eur=job-20221028-151132-249
# job_noneur=job-20221028-151136-545
# i=niddk_duerr_gsa
# # eur
# # niddk_duerr_gsa
# # QC  niddk_duerr_gsa job_eur=job-20221028-150619-477
# # noneur
# # niddk_duerr_gsa
# # QC  niddk_duerr_gsa job_noneur=job-20221028-150625-003
# job_eur=job-20221028-150619-477
# job_noneur=job-20221028-150625-003

# i=niddk_cho_gsa
# # eur
# # niddk_cho_gsa
# # QC  niddk_cho_gsa job_eur=job-20221028-150436-193
# # noneur
# # niddk_cho_gsa
# # QC  niddk_cho_gsa job_noneur=job-20221028-150442-146
# job_eur=job-20221028-150436-193
# job_noneur=job-20221028-150442-146

# i=newberry_share_gsa
# # eur
# # newberry_share_gsa
# # QC  newberry_share_gsa job_eur=job-20221028-150202-962
# # noneur
# # newberry_share_gsa
# # QC  newberry_share_gsa job_noneur=job-20221028-150206-296
# job_eur=job-20221028-150202-962
# job_noneur=job-20221028-150206-296


# i=moayyedi_imagine_gsa
# # eur
# # moayyedi_imagine_gsa
# # QC  moayyedi_imagine_gsa job_eur=job-20221028-145332-639
# # noneur
# # moayyedi_imagine_gsa
# # QC  moayyedi_imagine_gsa job_noneur=job-20221028-145336-297
# job_eur=job-20221028-145332-639
# job_noneur=job-20221028-145336-297

# i=mcgovern_gsa
# # eur
# # mcgovern_gsa
# # QC  mcgovern_gsa job_eur=job-20221028-145034-193
# # noneur
# # mcgovern_gsa
# # QC  mcgovern_gsa job_noneur=job-20221028-145045-450
# job_eur=job-20221028-145034-193
# job_noneur=job-20221028-145045-450
# i=mccauley_new_gsa
# # eur
# # mccauley_new_gsa
# # QC  mccauley_new_gsa job_eur=job-20221028-144428-575
# # noneur
# # mccauley_new_gsa
# # QC  mccauley_new_gsa job_noneur=job-20221028-144526-895
# job_eur=job-20221028-144428-575
# job_noneur=job-20221028-144526-895

# i=lewis_sparc_gsa
# # eur
# # lewis_sparc_gsa
# # QC  lewis_sparc_gsa job_eur=job-20221028-144321-016
# # noneur
# # lewis_sparc_gsa
# # QC  lewis_sparc_gsa job_noneur=job-20221028-144328-731
# job_eur=job-20221028-144321-016
# job_noneur=job-20221028-144328-731

# i=hyams_protect_gsa
# # eur
# # hyams_protect_gsa
# # QC  hyams_protect_gsa job_eur=job-20221028-144206-512
# # noneur
# # hyams_protect_gsa
# # QC  hyams_protect_gsa job_noneur=job-20221028-144209-106
# job_eur=job-20221028-144206-512
# job_noneur=job-20221028-144209-106

# i=helmsley_xavier_prism_gsa
# # eur
# # helmsley_xavier_prism_gsa
# # QC  helmsley_xavier_prism_gsa job_eur=job-20221028-144028-179
# # noneur
# # helmsley_xavier_prism_gsa
# # QC  helmsley_xavier_prism_gsa job_noneur=job-20221028-144031-888
# job_eur=job-20221028-144028-179
# job_noneur=job-20221028-144031-888

# i=helmsley_prism_gsa
# # eur
# # helmsley_prism_gsa
# # QC  helmsley_prism_gsa job_eur=job-20221028-143923-127
# # noneur
# # helmsley_prism_gsa
# # QC  helmsley_prism_gsa job_noneur=job-20221028-143925-706
# job_eur=job-20221028-143923-127
# job_noneur=job-20221028-143925-706

# i=franchimont_gsa
# # eur
# # franchimont_gsa
# # QC  franchimont_gsa job_eur=job-20221028-143404-938
# # noneur
# # franchimont_gsa
# # QC  franchimont_gsa job_noneur=job-20221028-143411-135
# job_eur=job-20221028-143404-938
# job_noneur=job-20221028-143411-135


# i=franke_gsa
# # eur
# # franke_gsa
# # QC  franke_gsa job_eur=job-20221028-143557-158
# # noneur
# # franke_gsa
# # QC  franke_gsa job_noneur=job-20221028-143600-465
# job_eur=job-20221028-143557-158
# job_noneur=job-20221028-143600-465

# i=farkkila_gsa
# # eur
# # farkkila_gsa
# 
# # QC  farkkila_gsa job_eur=job-20221028-143233-873
# job_eur=job-20221028-143233-873

# i=bernstein_gsa
# # eur
# # bernstein_gsa
# # QC  bernstein_gsa job_eur=job-20221028-143103-573
# # noneur
# # bernstein_gsa
# # QC  bernstein_gsa job_noneur=job-20221028-143105-854
# job_eur=job-20221028-143103-573
# job_noneur=job-20221028-143105-854


# i=cedars_gsa
# # eur
# # cedars_gsa
# # QC  cedars_gsa job_eur=job-20221028-142722-218
# # noneur
# # cedars_gsa
# # QC  cedars_gsa job_noneur=job-20221028-142732-420
# job_eur=job-20221028-142722-218
# job_noneur=job-20221028-142732-420

# i=ccfa_gsa
# # eur
# # ccfa_gsa
# # QC  ccfa_gsa job_eur=job-20221028-142536-627
# # noneur
# # ccfa_gsa
# # QC  ccfa_gsa job_noneur=job-20221028-142542-032
# job_eur=job-20221028-142536-627
# job_noneur=job-20221028-142542-032

# i=mccauley_gsa
# # eur
# # mccauley_gsa
# # QC  mccauley_gsa job_eur=job-20221028-142309-298
# # noneur
# # mccauley_gsa
# # QC  mccauley_gsa job_noneur=job-20221028-142312-565
# job_eur=job-20221028-142309-298
# job_noneur=job-20221028-142312-565

# i=cedars_omni_old_gwas
# # eur
# # cedars_omni_old_gwas
# # QC  cedars_omni_old_gwas job_eur=job-20221028-142100-583
# # noneur
# # cedars_omni_old_gwas
# # QC  cedars_omni_old_gwas job_noneur=job-20221028-142105-494
# job_eur=job-20221028-142100-583
# job_noneur=job-20221028-142105-494


# i=cedars_610k_old_gwas
# # eur
# # cedars_610k_old_gwas
# # QC  cedars_610k_old_gwas job_eur=job-20221028-141922-908
# # noneur
# # cedars_610k_old_gwas
# # QC  cedars_610k_old_gwas job_noneur=job-20221028-141930-599
# job_eur=job-20221028-141922-908
# job_noneur=job-20221028-141930-599

# i=cedars_370k_old_gwas
# # eur
# # cedars_370k_old_gwas
# # QC  cedars_370k_old_gwas job_eur=job-20221028-141728-153
# # noneur
# # cedars_370k_old_gwas
# # QC  cedars_370k_old_gwas job_noneur=job-20221028-141731-494
# job_eur=job-20221028-141728-153
# job_noneur=job-20221028-141731-494


# i=belgium_inf2_old_gwas
# # eur
# # belgium_inf2_old_gwas
# # QC  belgium_inf2_old_gwas job_eur=job-20221028-141555-705
# # noneur
# # belgium_inf2_old_gwas
# # QC  belgium_inf2_old_gwas job_noneur=job-20221028-141557-769
# job_eur=job-20221028-141555-705
# job_noneur=job-20221028-141557-769

# i=belgium_inf1_old_gwas
# # eur
# # belgium_inf1_old_gwas
# # QC  belgium_inf1_old_gwas job_eur=job-20221028-141328-933
# # noneur
# # belgium_inf1_old_gwas
# # QC  belgium_inf1_old_gwas job_noneur=job-20221028-141335-933
# job_eur=job-20221028-141328-933
# job_noneur=job-20221028-141335-933

# i=german_affy6_old_gwas
# # eur
# # german_affy6_old_gwas
# # QC  german_affy6_old_gwas job_eur=job-20221028-140606-591
# # noneur
# # german_affy6_old_gwas
# # QC  german_affy6_old_gwas job_noneur=job-20221028-140618-316
# job_eur=job-20221028-140606-591
# job_noneur=job-20221028-140618-316

# i=norway_affy6_old_gwas
# # eur
# # norway_affy6_old_gwas
# # QC  norway_affy6_old_gwas job_eur=job-20221028-141003-559
# # noneur
# # norway_affy6_old_gwas
# # QC  norway_affy6_old_gwas job_noneur=job-20221028-141007-042
# job_eur=job-20221028-141003-559
# job_noneur=job-20221028-141007-042


# i=finland_illugwas
# # eur
# # finland_illugwas
# # QC  finland_illugwas job_eur=job-20221028-140502-937
# # noneur
# # finland_illugwas
# # QC  finland_illugwas job_noneur=job-20221028-140505-075
# job_eur=job-20221028-140502-937
# job_noneur=job-20221028-140505-075

# i=prism_nfe_gwas
# # eur
# # prism_nfe_gwas
# # QC  prism_nfe_gwas job_eur=job-20221028-140320-463
# # noneur
# # prism_nfe_gwas
# # QC  prism_nfe_gwas job_noneur=job-20221028-140323-319
# job_eur=job-20221028-140320-463
# job_noneur=job-20221028-140323-319

# i=prism_nfe_gsa
# # eur
# # prism_nfe_gsa
# # QC  prism_nfe_gsa job_eur=job-20221028-140153-171
# # noneur
# # prism_nfe_gsa
# # QC  prism_nfe_gsa job_noneur=job-20221028-140155-951
# job_eur=job-20221028-140153-171
# job_noneur=job-20221028-140155-951

# i=belgium_vermeire_gsa
# # eur
# # belgium_vermeire_gsa
# # QC  belgium_vermeire_gsa job_eur=job-20221028-135934-257
# # noneur
# # belgium_vermeire_gsa
# # QC  belgium_vermeire_gsa job_noneur=job-20221028-135943-868
# job_eur=job-20221028-135934-257
# job_noneur=job-20221028-135943-868

# i=belgium_franchimont_gsa
# # eur
# # belgium_franchimont_gsa
# # QC  belgium_franchimont_gsa job_eur=job-20221028-134713-855
# # noneur
# # belgium_franchimont_gsa
# # QC  belgium_franchimont_gsa job_noneur=job-20221028-135903-226
# job_eur=job-20221028-134713-855
# job_noneur=job-20221028-135903-226

# i=belgium_louis_gsa
# # eur
# # belgium_louis_gsa
# 
# # QC  belgium_louis_gsa 
# # noneur
# # belgium_louis_gsa
# 
# # QC  belgium_louis_gsa 
# job_eur=job-20221028-134514-177
# job_noneur=job-20221028-134524-436

# i=niddk_feinstein_gsa
# # eur
# # niddk_feinstein_gsa
# # QC  niddk_feinstein_gsa job_eur=job-20221028-133402-226
# # noneur
# # niddk_feinstein_gsa
# # QC  niddk_feinstein_gsa job_noneur=job-20221028-133449-658
# job_eur=job-20221028-133402-226
# job_noneur=job-20221028-133449-658

# i=lithuania_gsa
# # eur
# # lithuania_gsa
# # QC  lithuania_gsa job_eur=job-20221028-134036-768
# # noneur
# # lithuania_gsa
# # QC  lithuania_gsa job_noneur=job-20221028-134320-607
# job_eur=job-20221028-134036-768
# job_noneur=job-20221028-134320-607

# i=basque_gsa
# # eur
# # basque_gsa
# 
# # QC  basque_gsa job_eur=job-20221028-133533-078
# # noneur
# # basque_gsa
# 
# # QC  basque_gsa job_noneur=job-20221028-134013-337
# job_eur=job-20221028-133533-078
# job_noneur=job-20221028-134013-337


# i=niddk_broad_gsa
# # eur
# # niddk_broad_gsa
# 
# # QC  niddk_broad_gsa job_eur=job-20221028-131918-632
# # noneur
# # niddk_broad_gsa
# 
# # QC  niddk_broad_gsa job_noneur=job-20221028-131929-940
# job_eur=job-20221028-131918-632
# job_noneur=job-20221028-131929-940

# i=sweden_gsa
# # eur
# # sweden_gsa
# 
# # QC  sweden_gsa job_eur=job-20221028-131846-401
# # noneur
# # sweden_gsa
# 
# # QC  sweden_gsa job_noneur=job-20221028-131850-883
# job_eur=job-20221028-131846-401
# job_noneur=job-20221028-131850-883

# 
# i=slovenia_gsa
# # eur
# # slovenia_gsa
# # QC  slovenia_gsa job_eur=job-20221028-130042-257
# job_eur=job-20221028-130042-257

# i=netherlands_gsa
# # eur
# # netherlands_gsa
# # QC  netherlands_gsa job_eur=job-20221028-125754-114
# # noneur
# # netherlands_gsa
# # QC  netherlands_gsa job_noneur=job-20221028-125806-822
# job_eur=job-20221028-125754-114
# job_noneur=job-20221028-125806-822

# i=kiel_austria_sibdcs_gsa
# # eur
# # kiel_austria_sibdcs_gsa
# # QC  kiel_austria_sibdcs_gsa job_eur=job-20221028-125312-215
# # noneur
# # kiel_austria_sibdcs_gsa
# # QC  kiel_austria_sibdcs_gsa job_noneur=job-20221028-125341-230
# job_eur=job-20221028-125312-215
# job_noneur=job-20221028-125341-230


# i=italy_gsa
# # eur
# # italy_gsa
# # QC  italy_gsa job_eur=job-20221028-124424-465
# # noneur
# # italy_gsa
# # QC  italy_gsa job_noneur=job-20221028-125227-049
# job_eur=job-20221028-124424-465
# job_noneur=job-20221028-125227-049

# i=spain_gsa
# # eur
# # spain_gsa
# # QC  spain_gsa job_eur=job-20221028-124320-928
# # noneur
# # spain_gsa
# # QC  spain_gsa job_noneur=job-20221028-124332-403
# job_eur=job-20221028-124320-928
# job_noneur=job-20221028-124332-403

# i=pittsburgh_gsa
# # eur
# # pittsburgh_gsa
# # QC  pittsburgh_gsa job_eur=job-20221028-123907-356
# # noneur
# # pittsburgh_gsa
# # QC  pittsburgh_gsa job_noneur=job-20221028-123920-150
# job_eur=job-20221028-123907-356
# job_noneur=job-20221028-123920-150


# i=gwas2
# eur
# gwas2
# QC  gwas2 job_eur=job-20221028-121008-738
# noneur
# gwas2
# QC  gwas2 job_eur=job-20221028-121954-782
# job_eur=job-20221028-121008-738
# job_noneur=job-20221028-121954-782

# i=gwas1
# eur
# gwas1
# QC  gwas1 job_eur=job-20221028-122132-247
# noneur
# gwas1
# QC  gwas1 job_eur=job-20221028-122143-438
# job_eur=job-20221028-122132-247
# job_noneur=job-20221028-122143-438

# i=australia_omniexome 
# eur
# australia_omniexome
# QC  australia_omniexome job_eur=job-20221028-122423-284
# noneur
# australia_omniexome
# QC  australia_omniexome job_noneur=job-20221028-122914-772
# job_eur=job-20221028-122423-284
# job_noneur=job-20221028-122914-772


# i=all_hce
# eur
# QC  all_hce eur job-20221027-214602-052
# noneur
# QC  all_hce noneur job-20221027-215014-014
# job_eur=job-20221027-214602-052
# job_noneur=job-20221027-215014-014


# 
# studies=(australia_omniexome gwas1 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa)
# for i in ${studies[@]}
# do
# for j in ${ancestry[@]}
# do echo ${j} && python2 /path/to/user/scripts/IIBDGC/topmed_qc_${j}.py ${i}
# done 
# done

##################


{
  echo ${i}
  
  j=eur
  TOKEN='YOUR_API_TOKEN';
  curl -H "X-Auth-Token: $TOKEN" https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/${job_eur} \
  > ${path_gwas}imputed/${i}/qc/2022/logs/status_${i}_${j}_qc.json
  
  j=noneur
  TOKEN='YOUR_API_TOKEN';
  curl -H "X-Auth-Token: $TOKEN" https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/${job_noneur} \
  > ${path_gwas}imputed/${i}/qc/2022/logs/status_${i}_${j}_qc.json
  
  /software/R-4.3.1/bin/Rscript \
  /path/to/user/scripts/IIBDGC/create_imputation_qc_download_commnads.R ${i} \
  > ${path_gwas}imputed/${i}/qc/2022/logs/create_imputation_qc_download_commnads_${i}.R_out
  
  ### edit files:
  j=eur
  sed -i '2 i\path=\/lustre\/scratch123\/hgi\/projects\/ibdgwas\/IIBDGC\/imputed\/${i}\/qc\/2022\/' ${path_gwas}imputed/${i}/qc/2022/logs/donwload_qc_${i}_${j}
  sed -i '2 i\j='${j} ${path_gwas}imputed/${i}/qc/2022/logs/donwload_qc_${i}_${j}
  sed -i '2 i\i='${i} ${path_gwas}imputed/${i}/qc/2022/logs/donwload_qc_${i}_${j}
  
  j=noneur
  sed -i -e '1,4d' ${path_gwas}imputed/${i}/qc/2022/logs/donwload_qc_${i}_${j}
  sed -i '2 i\j='${j} ${path_gwas}imputed/${i}/qc/2022/logs/donwload_qc_${i}_${j}
  
  cat  ${path_gwas}imputed/${i}/qc/2022/logs/donwload_qc_${i}_eur  ${path_gwas}imputed/${i}/qc/2022/logs/donwload_qc_${i}_noneur \
  > ${path_gwas}imputed/${i}/qc/2022/logs/donwload_qc_${i}
  rm ${path_gwas}imputed/${i}/qc/2022/logs/donwload_qc_${i}_*
    
    sed -i -e 's/-o .*.txt --create-dirs/-P ${path}\/${j}/g' ${path_gwas}imputed/${i}/qc/2022/logs/donwload_qc_${i}
  sed -i -e 's/curl -L/wget --tries=75 -c/g' ${path_gwas}imputed/${i}/qc/2022/logs/donwload_qc_${i}
  more ${path_gwas}imputed/${i}/qc/2022/logs/donwload_qc_${i}
  
  cd /path/to/ibdgwas/IIBDGC/imputed/${i}/qc/2022/logs/
    nohup bash /path/to/ibdgwas/IIBDGC/imputed/${i}/qc/2022/logs/donwload_qc_${i} &
    
    unset job_eur
  unset job_noneur
  
  # double check files:
  
  for j in ${ancestry[@]}
  do echo ${j} && wc -l /path/to/ibdgwas/IIBDGC/imputed/${i}/qc/2022/${j}/typed-only.txt
  done
  
}

for j in ${ancestry[@]}
do echo ${j} && wc -l /path/to/ibdgwas/IIBDGC/imputed/${i}/qc/2022/${j}/typed-only.txt
done
  

### EXCLUDE VARIANTS FLAGGED BY SERVER
 
## /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","swedish_uc_old_gwas"
           ,"mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

length(cohorts)
# [1] 58

ancestry<-c("eur","noneur")

for (j in 1:length(cohorts)) {
  
  print(cohorts[j])
  
  bim<-read.table(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim",sep=""),head=F)
  
  files_list<-list.files(paste(path,"imputed/",cohorts[j],"/qc/2022/eur/",sep=""))
  files_snp_exclude<-files_list[grep("snps-excluded",files_list)]
  
  for (i in 1:length(ancestry)) {
    file_eth<-paste(path,"imputed/",cohorts[j],"/qc/2022/",ancestry[i],"/",files_snp_exclude[1],sep="")
    if(file.exists(file_eth)) {
      tmp<-fread(file_eth,head=T,sep="\t",fill=TRUE)
      tmp<-as.data.frame(tmp)
      colnames(tmp)<-c("Position","FilterType","Info")
      if(i==1){
        snp_excl<-tmp
      } else {
        snp_excl<-rbind(snp_excl,tmp)
      }
    }
    rm(file_eth)
  }

  snp_excl<-snp_excl[!duplicated(snp_excl$Position),]
  
  snp_excl<-snp_excl[which(snp_excl$FilterType=="Allele mismatch"),]
  
  
  files_typed_only<-files_list[grep("typed-only",files_list)]

  for (i in 1:length(ancestry)) {
    file_eth<-paste(path,"imputed/",cohorts[j],"/qc/2022/",ancestry[i],"/",files_typed_only[1],sep="")
    if(file.exists(file_eth)) {
      tmp<-fread(file_eth,head=T,sep="\t",fill=TRUE)
      tmp<-as.data.frame(tmp)
      colnames(tmp)<-c("Position")
      if(i==1){
        snp_tyonly<-tmp
      } else {
        snp_tyonly<-rbind(snp_tyonly,tmp)
      }
    }
    rm(file_eth)
  }
  snp_tyonly<-snp_tyonly[!duplicated(snp_tyonly$Position),,drop=F]

  list_exclude<-rbind(snp_tyonly,snp_excl[,1,drop=F])
  
  list_exclude$SNP<-gsub(":","_",list_exclude$Position)
  list_exclude$SNP<-sub("_",":",list_exclude$SNP)
  list_exclude$SNP<-sub("chrX","chr23",list_exclude$SNP)
  list_exclude$SNP<-sub("chr","",list_exclude$SNP)
  
  eur<-read.table(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim",sep=""),head=F)
  
  print(paste("N variants exclude in EUR: ",nrow(list_exclude[which(list_exclude$SNP %in% eur$V2),]),sep=""))
  # [1] 101160     2
  
  file_noneur<-paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim",sep="")
  if(file.exists(file_noneur)) {
    noneur<-read.table(file_noneur,head=F)
    print(paste("N variants exclude in non-EUR: ",nrow(list_exclude[which(list_exclude$SNP %in% noneur$V2),]),sep=""))
    # [1] 123628     2
  }
  
  
  write.table(list_exclude[,"SNP",drop=F],paste(path,"imputed/",cohorts[j],"/qc/2022/list_variants_to_exclude",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  
  rm(list=ls()[!ls() %in% c("path","j","cohorts","ancestry")]) 

}

# double check why all_hce has >100K variants to be excluded

# frq<-read.table(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_",ancestry[i],"_ctr.frq.counts",sep=""),head=T)
# frq$maf<-frq$C1/(frq$C1+frq$C2)
# frq$maf[which(frq$maf>0.5)]<-1-frq$maf[which(frq$maf>0.5)]
# 
# eur<-merge(eur,frq[,c("SNP","maf")],by.x="V2",by.y="SNP")
# tmp<-eur[which(eur$V2 %in% list_exclude$SNP),]
# summary(tmp$maf)
# # Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# # 0.0000000 0.0000485 0.0000971 0.0006580 0.0004367 0.4938338 
# summary(frq$maf)
# # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# # 0.000000 0.000097 0.071477 0.149551 0.293495 0.500000 


############
# 
# [1] "australia_omniexome"
# [1] "N variants exclude in EUR: 43838"
# [1] "N variants exclude in non-EUR: 7640"
# [1] "gwas1"
# [1] "N variants exclude in EUR: 568"
# [1] "N variants exclude in non-EUR: 138"
# [1] "gwas2"
# [1] "N variants exclude in EUR: 683"
# [1] "N variants exclude in non-EUR: 303"
# [1] "all_hce"
# [1] "N variants exclude in EUR: 100966"
# [1] "N variants exclude in non-EUR: 123269"
# [1] "pittsburgh_gsa"
# [1] "N variants exclude in EUR: 6429"
# [1] "N variants exclude in non-EUR: 823"
# [1] "spain_gsa"
# [1] "N variants exclude in EUR: 794"
# [1] "N variants exclude in non-EUR: 328"
# [1] "italy_gsa"
# [1] "N variants exclude in EUR: 12776"
# [1] "N variants exclude in non-EUR: 551"
# [1] "kiel_austria_sibdcs_gsa"
# [1] "N variants exclude in EUR: 36715"
# [1] "N variants exclude in non-EUR: 25377"
# [1] "netherlands_gsa"
# [1] "N variants exclude in EUR: 20783"
# [1] "N variants exclude in non-EUR: 25901"
# [1] "slovenia_gsa"
# [1] "N variants exclude in EUR: 5452"
# [1] "sweden_gsa"
# [1] "N variants exclude in EUR: 13067"
# [1] "N variants exclude in non-EUR: 6965"
# [1] "niddk_broad_gsa"
# [1] "N variants exclude in EUR: 29755"
# [1] "N variants exclude in non-EUR: 35614"
# [1] "niddk_feinstein_gsa"
# [1] "N variants exclude in EUR: 39722"
# [1] "N variants exclude in non-EUR: 55658"
# [1] "basque_gsa"
# [1] "N variants exclude in EUR: 11787"
# [1] "N variants exclude in non-EUR: 662"
# [1] "prism_nfe_gsa"
# [1] "N variants exclude in EUR: 9658"
# [1] "N variants exclude in non-EUR: 8896"
# [1] "lithuania_gsa"
# [1] "N variants exclude in EUR: 21202"
# [1] "N variants exclude in non-EUR: 1353"
# [1] "belgium_louis_gsa"
# [1] "N variants exclude in EUR: 17843"
# [1] "N variants exclude in non-EUR: 1336"
# [1] "belgium_franchimont_gsa"
# [1] "N variants exclude in EUR: 16413"
# [1] "N variants exclude in non-EUR: 2542"
# [1] "belgium_vermeire_gsa"
# [1] "N variants exclude in EUR: 29483"
# [1] "N variants exclude in non-EUR: 68304"
# [1] "prism_nfe_gwas"
# [1] "N variants exclude in EUR: 143"
# [1] "N variants exclude in non-EUR: 144"
# [1] "finland_illugwas"
# [1] "N variants exclude in EUR: 141"
# [1] "N variants exclude in non-EUR: 117"
# [1] "german_affy6_old_gwas"
# [1] "N variants exclude in EUR: 1965"
# [1] "N variants exclude in non-EUR: 545"
# [1] "norway_affy6_old_gwas"
# [1] "N variants exclude in EUR: 1279"
# [1] "N variants exclude in non-EUR: 244"
# [1] "belgium_inf1_old_gwas"
# [1] "N variants exclude in EUR: 161"
# [1] "N variants exclude in non-EUR: 153"
# [1] "belgium_inf2_old_gwas"
# [1] "N variants exclude in EUR: 156"
# [1] "N variants exclude in non-EUR: 79"
# [1] "niddk_old_gwas"
# [1] "N variants exclude in EUR: 166"
# [1] "N variants exclude in non-EUR: 158"
# [1] "cedars_370k_old_gwas"
# [1] "N variants exclude in EUR: 328"
# [1] "N variants exclude in non-EUR: 319"
# [1] "cedars_610k_old_gwas"
# [1] "N variants exclude in EUR: 725"
# [1] "N variants exclude in non-EUR: 786"
# [1] "cedars_omni_old_gwas"
# [1] "N variants exclude in EUR: 3194"
# [1] "N variants exclude in non-EUR: 3499"
# [1] "swedish_uc_old_gwas"
# [1] "N variants exclude in EUR: 154"
# [1] "N variants exclude in non-EUR: 124"
# [1] "mccauley_gsa"
# [1] "N variants exclude in EUR: 14488"
# [1] "N variants exclude in non-EUR: 1212"
# [1] "ccfa_gsa"
# [1] "N variants exclude in EUR: 22848"
# [1] "N variants exclude in non-EUR: 28277"
# [1] "cedars_gsa"
# [1] "N variants exclude in EUR: 29434"
# [1] "N variants exclude in non-EUR: 46839"
# [1] "bernstein_gsa"
# [1] "N variants exclude in EUR: 8997"
# [1] "N variants exclude in non-EUR: 3458"
# [1] "farkkila_gsa"
# [1] "N variants exclude in EUR: 3845"
# [1] "franchimont_gsa"
# [1] "N variants exclude in EUR: 21922"
# [1] "N variants exclude in non-EUR: 12611"
# [1] "franke_gsa"
# [1] "N variants exclude in EUR: 12404"
# [1] "N variants exclude in non-EUR: 6306"
# [1] "helmsley_prism_gsa"
# [1] "N variants exclude in EUR: 144"
# [1] "N variants exclude in non-EUR: 145"
# [1] "helmsley_xavier_prism_gsa"
# [1] "N variants exclude in EUR: 146"
# [1] "N variants exclude in non-EUR: 146"
# [1] "hyams_protect_gsa"
# [1] "N variants exclude in EUR: 9058"
# [1] "N variants exclude in non-EUR: 13623"
# [1] "lewis_sparc_gsa"
# [1] "N variants exclude in EUR: 26764"
# [1] "N variants exclude in non-EUR: 32670"
# [1] "mccauley_new_gsa"
# [1] "N variants exclude in EUR: 20307"
# [1] "N variants exclude in non-EUR: 11273"
# [1] "mcgovern_gsa"
# [1] "N variants exclude in EUR: 39608"
# [1] "N variants exclude in non-EUR: 57246"
# [1] "moayyedi_imagine_gsa"
# [1] "N variants exclude in EUR: 14482"
# [1] "N variants exclude in non-EUR: 26993"
# [1] "newberry_share_gsa"
# [1] "N variants exclude in EUR: 10822"
# [1] "N variants exclude in non-EUR: 14446"
# [1] "niddk_cho_gsa"
# [1] "N variants exclude in EUR: 22837"
# [1] "N variants exclude in non-EUR: 24972"
# [1] "niddk_duerr_gsa"
# [1] "N variants exclude in EUR: 24095"
# [1] "N variants exclude in non-EUR: 10190"
# [1] "niddk_rioux_gsa"
# [1] "N variants exclude in EUR: 14912"
# [1] "N variants exclude in non-EUR: 11732"
# [1] "niddk_silverberg_gsa"
# [1] "N variants exclude in EUR: 27066"
# [1] "N variants exclude in non-EUR: 31814"
# [1] "palotie_hus_gsa"
# [1] "N variants exclude in EUR: 11346"
# [1] "N variants exclude in non-EUR: 2265"
# [1] "pekow_share_gsa"
# [1] "N variants exclude in EUR: 10510"
# [1] "N variants exclude in non-EUR: 11442"
# [1] "rioux_igenomed_gsa"
# [1] "N variants exclude in EUR: 4133"
# [1] "N variants exclude in non-EUR: 2635"
# [1] "sands_msccr_gsa"
# [1] "N variants exclude in EUR: 15564"
# [1] "N variants exclude in non-EUR: 29108"
# [1] "stampfer_gsa"
# [1] "N variants exclude in EUR: 22654"
# [1] "N variants exclude in non-EUR: 9973"
# [1] "vermeire_gsa"
# [1] "N variants exclude in EUR: 30356"
# [1] "N variants exclude in non-EUR: 14375"
# [1] "weersma_gsa"
# [1] "N variants exclude in EUR: 8644"
# [1] "N variants exclude in non-EUR: 7521"
# [1] "xavier_prism_gsa"
# [1] "N variants exclude in EUR: 12005"
# [1] "N variants exclude in non-EUR: 14229"
# [1] "xavier_share_gsa"
# [1] "N variants exclude in EUR: 13044"
# [1] "N variants exclude in non-EUR: 14059"


###########################################################################
# 23.2 Exclude variants from files and Convert ped/map files to VCF files #
###########################################################################


studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=500

ancestry=(eur noneur)

for chr in {1..23}
do
for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
bsub -J"pl2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_post_server_qc_${j}_${i}_chr${chr} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_post_server_qc_${j}_${i}_chr${chr} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_${j}_nomonom \
--allow-no-sex \
--keep-allele-order \
--exclude ${path_gwas}imputed/${i}/qc/2022/list_variants_to_exclude \
--chr ${chr} \
--output-chr chr26 --recode vcf-iid --out ${path_gwas}imputation_ready/${i}/${i}_${j}_chr${chr}_2022"
done
done
done


for i in ${studies[@]}
do 
echo ${i} && for j in ${ancestry[@]}
do 
ls -la ${path_gwas}imputation_ready/${i}/${i}_${j}_chr22_2022.vcf
done
done



###############################
# 23.3 Create a sorted vcf.gz #
###############################

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

ancestry=(eur noneur)
MEM=1200
for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
for chr in {1..22}
do
bsub -J"pl2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_2_post_server_qc_${j}_${i}_chr${chr} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_2_post_server_qc_${j}_${i}_chr${chr} \
"/path/to/software/bcftools-1.16/./bcftools sort ${path_gwas}imputation_ready/${i}/${i}_${j}_chr${chr}_2022.vcf \
-Oz -o ${path_gwas}imputation_ready/${i}/${i}_${j}_chr${chr}_2022.vcf.gz"
done
done
done

# rm intermediate files:

for i in ${studies[@]}
do rm ${path_gwas}imputation_ready/${i}/${i}_*_chr*.vcf
done
for i in ${studies[@]}
do rm ${path_gwas}imputation_ready/${i}/*.hh
done
for i in ${studies[@]}
do rm ${path_gwas}imputation_ready/${i}/*.log
done
for i in ${studies[@]}
do rm ${path_gwas}imputation_ready/${i}/*.nosex
done



##################################################################################################################################################################

# compare 2022 version (n samples) vs 2020:

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa )

ancestry=(eur noneur)
MEM=200

for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
bsub -J"pl2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_2.1_post_server_qc_${j}_${i}_chr${chr} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_2.1_post_server_qc_${j}_${i}_chr${chr} \
"/path/to/software/bcftools-1.16/./bcftools \
query -l ${path_gwas}imputation_ready/${i}/${i}_${j}_chr22_2022.vcf.gz > ${path_gwas}imputation_ready/test/${i}_${j}_2022_samples"
done
done

for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
bsub -J"pl2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_2.2_post_server_qc_${j}_${i}_chr${chr} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_2.2_post_server_qc_${j}_${i}_chr${chr} \
"/path/to/software/bcftools-1.16/./bcftools \
query -l ${path_gwas}imputation_ready/${i}/${i}_${j}_chr22.vcf.gz > ${path_gwas}imputation_ready/test/${i}_${j}_2020_samples"
done
done


####################

## /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","swedish_uc_old_gwas"
           ,"mccauley_gsa","ccfa_gsa","cedars_gsa")

length(cohorts)


ancestry<-c("eur","noneur")

for (j in 1:length(cohorts)) {
  
  print(cohorts[j])
  
  if (!cohorts[j] %in% c("slovenia_gsa","finland_illugwas")) {
    for (i in 1:length(ancestry) ) {
      s1<-read.table(paste(path,"imputation_ready/test/",cohorts[j],"_",ancestry[i],"_2020_samples",sep=""),head=F)
      s2<-read.table(paste(path,"imputation_ready/test/",cohorts[j],"_",ancestry[i],"_2022_samples",sep=""),head=F)
      
      if(nrow(s2[which(!s2$V1 %in% s1$V1),,drop=F])>0) {
        print(cohorts[j])
      }
      rm(s1,s2)
    }
  } else {
    i=1
    s1<-read.table(paste(path,"imputation_ready/test/",cohorts[j],"_",ancestry[i],"_2020_samples",sep=""),head=F)
    s2<-read.table(paste(path,"imputation_ready/test/",cohorts[j],"_",ancestry[i],"_2022_samples",sep=""),head=F)
    
    if(nrow(s2[which(!s2$V1 %in% s1$V1),,drop=F])>0) {
      print(cohorts[j])
    }
    rm(s1,s2)
  }
  
  
}


  
####################

##################################################################################################################################################################

###################################################
# 23.4 Submit files for first round of imputation #
###################################################


## submit first round of imputation to get rsqr and empirical rsq - for old studies keep old list - double check later that all is ok
 
studies=(bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

# some of the old studies (more samples now post-qc to be resubmitted)
studies=(gwas1 spain_gsa sweden_gsa)
studies=(all_hce niddk_old_gwas australia_omniexome gwas2 pittsburgh_gsa
         italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa
         niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa
         belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa
         prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas
         norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas
         cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas
         swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa)


path_gwas=/path/to/ibdgwas/IIBDGC/
ancestry=(eur noneur)

i=
for j in ${ancestry[@]}
do echo ${i} && ls -la ${path_gwas}imputation_ready/${i}/${i}_${j}_chr*_2022.vcf.gz |  wc -l
done
for j in ${ancestry[@]}
do echo ${j} && python2 /path/to/user/scripts/IIBDGC/topmed_imputation_${j}.py ${i}
done


###########

# i=spain_gsa - downloading 2
# # eur
# # spain_gsa
# # QC  spain_gsa job_eur=job-20221031-102053-143
# # noneur
# # spain_gsa
# # QC  spain_gsa job_noneur=job-20221031-103326-605
# job_eur=job-20221031-102053-143
# job_noneur=job-20221031-103326-605

# i=sweden_gsa - downloading 2
# # eur
# # sweden_gsa
# # Your job was successfully added to the job queue.
# # QC  sweden_gsa job_eur=job-20221031-114303-808
# # noneur
# # sweden_gsa
# # Your job was successfully added to the job queue.
# # QC  sweden_gsa job_noneur=job-20221031-114314-981
# job_eur=job-20221031-114303-808
# job_noneur=job-20221031-114314-981

# i=helmsley_prism_gsa - downloading 2
# # eur
# # helmsley_prism_gsa
# # QC  helmsley_prism_gsa job_eur=job-20221031-095649-051
# # noneur
# # helmsley_prism_gsa
# # QC  helmsley_prism_gsa job_noneur=job-20221031-095659-918
# job_eur=job-20221031-095649-051
# job_noneur=job-20221031-095659-918

# i=bernstein_gsa - downloading 2
# # eur
# # bernstein_gsa
# # QC  bernstein_gsa job_eur=job-20221031-070314-942
# # noneur
# # bernstein_gsa
# # QC  bernstein_gsa job_noneur=job-20221031-075433-968
# job_eur=job-20221031-070314-942
# job_noneur=job-20221031-075433-968

# i=newberry_share_gsa - downloading 2
# # eur
# # newberry_share_gsa
# # QC  newberry_share_gsa job_eur=job-20221031-044455-435
# # noneur
# # newberry_share_gsa
# # QC  newberry_share_gsa job_noneur=job-20221031-070247-080
# job_eur=job-20221031-044455-435
# job_noneur=job-20221031-070247-080

# i=farkkila_gsa - downloading 2
# # eur
# # farkkila_gsa
# # QC  farkkila_gsa job_eur=job-20221031-044001-863
# job_eur=job-20221031-044001-863

# i=franke_gsa - downloading 2
# # eur
# # franke_gsa
# # QC  franke_gsa job_eur=job-20221030-230034-501
# # noneur
# # franke_gsa
# # QC  franke_gsa job_noneur=job-20221031-043944-711
# job_eur=job-20221030-230034-501
# job_noneur=job-20221031-043944-711

# i=franchimont_gsa - downloading 2
# # eur
# # franchimont_gsa
# # QC  franchimont_gsa job_eur=job-20221030-220052-158
# # noneur
# # franchimont_gsa
# # QC  franchimont_gsa job_noneur=job-20221030-230005-586
# job_eur=job-20221030-220052-158
# job_noneur=job-20221030-230005-586

# i=helmsley_xavier_prism_gsa - downloading 2
# # eur
# # helmsley_xavier_prism_gsa
# # QC  helmsley_xavier_prism_gsa job_eur=job-20221030-210919-297
# # noneur
# # helmsley_xavier_prism_gsa
# # QC  helmsley_xavier_prism_gsa job_noneur=job-20221030-210939-348
# job_eur=job-20221030-210919-297
# job_noneur=job-20221030-210939-348

# i=lewis_sparc_gsa - downloading 2
# # eur
# # lewis_sparc_gsa
# # QC  lewis_sparc_gsa job_eur=job-20221030-193226-583
# # noneur
# # lewis_sparc_gsa
# # QC  lewis_sparc_gsa job_noneur=job-20221030-193233-823
# job_eur=job-20221030-193226-583
# job_noneur=job-20221030-193233-823


# i=hyams_protect_gsa - downloading 2
# # eur
# # hyams_protect_gsa
# # QC  hyams_protect_gsa job_eur=job-20221030-200731-062
# # noneur
# # hyams_protect_gsa
# # QC  hyams_protect_gsa job_noneur=job-20221030-202624-178
# job_eur=job-20221030-200731-062
# job_noneur=job-20221030-202624-178

# i=moayyedi_imagine_gsa - downloading 2
# # eur
# # moayyedi_imagine_gsa
# # QC  moayyedi_imagine_gsa job_eur=job-20221030-180107-779
# # noneur
# # moayyedi_imagine_gsa
# # QC  moayyedi_imagine_gsa job_noneur=job-20221030-184753-569
# job_eur=job-20221030-180107-779
# job_noneur=job-20221030-184753-569

# i=niddk_silverberg_gsa - downloading 2
# # eur
# # niddk_silverberg_gsa
# # QC  niddk_silverberg_gsa job_eur=job-20221030-140326-472
# # noneur
# # niddk_silverberg_gsa
# # QC  niddk_silverberg_gsa job_noneur=job-20221030-140333-044
# job_eur=job-20221030-140326-472
# job_noneur=job-20221030-140333-044


# i=niddk_duerr_gsa - downloading 2
# # eur
# # niddk_duerr_gsa
# # QC  niddk_duerr_gsa job_eur=job-20221030-164130-954
# # noneur
# # niddk_duerr_gsa
# # QC  niddk_duerr_gsa job_noneur=job-20221030-164143-293
# job_eur=job-20221030-164130-954
# job_noneur=job-20221030-164143-293

# i=niddk_rioux_gsa - downloading 2
# # eur
# # niddk_rioux_gsa
# # QC  niddk_rioux_gsa job_eur=job-20221030-140341-962
# # noneur
# # niddk_rioux_gsa
# # QC  niddk_rioux_gsa job_noneur=job-20221030-151948-570
# job_eur=job-20221030-140341-962
# job_noneur=job-20221030-151948-570


# i=niddk_rioux_gsa - downloading 2
# # eur
# # niddk_rioux_gsa
# # QC  niddk_rioux_gsa job_eur=job-20221030-140341-962
# # noneur
# # niddk_rioux_gsa
# # QC  niddk_rioux_gsa job_noneur=job-20221030-151948-570
# job_eur=job-20221030-140341-962
# job_noneur=job-20221030-151948-570

# i=palotie_hus_gsa - downloading 2
# # eur
# # palotie_hus_gsa
# # QC  palotie_hus_gsa job_eur=job-20221030-113111-956
# # noneur
# # palotie_hus_gsa
# # QC  palotie_hus_gsa job_noneur=job-20221030-115307-816
# job_eur=job-20221030-113111-956
# job_noneur=job-20221030-115307-816

# i=pekow_share_gsa - downloading 2
# # eur
# # pekow_share_gsa
# # QC  pekow_share_gsa job_eur=job-20221030-104328-299
# # noneur
# # pekow_share_gsa
# # QC  pekow_share_gsa job_noneur=job-20221030-104332-355
# job_eur=job-20221030-104328-299
# job_noneur=job-20221030-104332-355

# i=rioux_igenomed_gsa - downloading 2
# # eur
# # rioux_igenomed_gsa
# # QC  rioux_igenomed_gsa job_eur=job-20221030-080250-814
# # noneur
# # rioux_igenomed_gsa
# # QC  rioux_igenomed_gsa job_noneur=job-20221030-104301-984
# job_eur=job-20221030-080250-814
# job_noneur=job-20221030-104301-984

# i=sands_msccr_gsa - downloading 2
# # eur
# # sands_msccr_gsa
# # QC  sands_msccr_gsa job_eur=job-20221030-070343-813
# # noneur
# # sands_msccr_gsa
# # QC  sands_msccr_gsa job_noneur=job-20221030-070348-128
# job_eur=job-20221030-070343-813
# job_noneur=job-20221030-070348-128

# i=niddk_cho_gsa - downloading 2
# # eur
# # niddk_cho_gsa
# # QC  niddk_cho_gsa job_eur=job-20221030-032950-033
# # noneur
# # niddk_cho_gsa
# # QC  niddk_cho_gsa job_noneur=job-20221030-032955-561
# job_eur=job-20221030-032950-033
# job_noneur=job-20221030-032955-561


# i=mccauley_new_gsa - downloading 2
# # eur
# # mccauley_new_gsa
# # QC  mccauley_new_gsa job_eur=job-20221029-215856-971
# # noneur
# # mccauley_new_gsa
# # QC  mccauley_new_gsa job_noneur=job-20221030-032906-196
# job_eur=job-20221029-215856-971
# job_noneur=job-20221030-032906-196

# i=mcgovern_gsa - downloading 2
# # # eur
# # # mcgovern_gsa
# # # QC  mcgovern_gsa job_eur=job-20221029-151458-212
# # # noneur
# # # mcgovern_gsa
# # # QC  mcgovern_gsa job_noneur=job-20221029-215830-666
# job_eur=job-20221029-151458-212
# job_noneur=job-20221029-215830-666

i=vermeire_gsa
# # # eur
# # # vermeire_gsa
# # # QC  vermeire_gsa job_eur=job-20221029-092705-756
# # # noneur
# # # vermeire_gsa
# # # QC  vermeire_gsa job_noneur=job-20221029-092722-390
job_eur=job-20221029-092705-756
job_noneur=job-20221029-092722-390


# i=stampfer_gsa - downloading 2
# # # eur
# # # stampfer_gsa
# # # QC  stampfer_gsa job_eur=job-20221029-102249-793
# # # noneur
# # # stampfer_gsa
# # # QC  stampfer_gsa job_noneur=job-20221029-134725-792
# job_eur=job-20221029-102249-793
# job_noneur=job-20221029-134725-792

# i=weersma_gsa - downloading 2
# # eur
# # weersma_gsa
# # QC  weersma_gsa job_eur=job-20221029-080207-977
# # noneur
# # weersma_gsa
# # QC  weersma_gsa job_noneur=job-20221029-083859-382
# job_eur=job-20221029-080207-977
# job_noneur=job-20221029-083859-382

# i=xavier_prism_gsa - downloading 2
# # eur
# # xavier_prism_gsa
# # QC  xavier_prism_gsa job_eur=job-20221029-070125-496
# # noneur
# # xavier_prism_gsa
# # QC  xavier_prism_gsa job_noneur=job-20221029-070135-036
# job_eur=job-20221029-070125-496
# job_noneur=job-20221029-070135-036

# i=xavier_share_gsa - downloading 2
# # eur
# # xavier_share_gsa
# # QC  xavier_share_gsa job_eur=job-20221029-044915-186
# # noneur
# # xavier_share_gsa
# # QC  xavier_share_gsa job_noneur=job-20221029-044925-924
# job_eur=job-20221029-044915-186
# job_noneur=job-20221029-044925-924

# i=gwas1 - downloading 2
# # eur
# # gwas1
# # QC  gwas1 job_eur=job-20221028-171154-721
# # noneur
# # gwas1
# # QC  gwas1 job_noneur=job-20221028-171212-009
# job_eur=job-20221028-171154-721
# job_noneur=job-20221028-171212-009

###########

{
  echo ${i}
  
  j=eur
  TOKEN='YOUR_API_TOKEN';
  curl -H "X-Auth-Token: $TOKEN" https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/${job_eur} \
  > ${path_gwas}imputed/${i}/2022/logs/status_${i}_${j}_qc.json
  
  j=noneur
  TOKEN='YOUR_API_TOKEN';
  curl -H "X-Auth-Token: $TOKEN" https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/${job_noneur} \
  > ${path_gwas}imputed/${i}/2022/logs/status_${i}_${j}_qc.json
  
  /software/R-4.3.1/bin/Rscript \
  /path/to/user/scripts/IIBDGC/create_imputation_download_commnads.R ${i} \
  > ${path_gwas}imputed/${i}/2022/logs/create_imputation_download_commnads_${i}.R_out
  
  ### edit files:
  j=eur
  sed -i '2 i\path=\/lustre\/scratch123\/hgi\/projects\/ibdgwas\/IIBDGC\/imputed\/${i}\/2022\/' ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
  sed -i '2 i\j='${j} ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
  sed -i '2 i\i='${i} ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
  
  j=noneur
  sed -i -e '1,4d' ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
  sed -i '2 i\j='${j} ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
  
  cat  ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_eur  ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_noneur \
  > ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}
  rm ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_*
    
  sed -i -e 's/-o .*.zip --create-dirs/-P ${path}\/${j}/g' ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}
  sed -i -e 's/-o .*.md5 --create-dirs/-P ${path}\/${j}/g' ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}
  sed -i -e 's/curl -L/wget --tries=75 -c/g' ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}
  more ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}
  
  cd /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/logs/
  nohup bash /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/logs/donwload_imputed_results_${i} &
    
  unset job_eur
  unset job_noneur
  
  # double check files:
  
}

path_gwas=/path/to/ibdgwas/IIBDGC/
studies=(bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)
ancestry=(eur noneur)

for i in ${studies[@]}
do echo ${i} && for j in ${ancestry[@]}
do echo ${j} && ls -la /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/${j}/chr_*.zip | wc -l
done
done

for i in ${studies[@]}
do echo ${i} && tail -3 /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/logs/nohup.out
done


####################
# 23.5 Unzip files #
####################

studies=(gwas1 spain_gsa sweden_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa 
         helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa 
         newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa 
         rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)
ancestry=(eur noneur)
MEM=200

### CHR1 TO CHR22
for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
for chr in {1..22}
do
bsub -J"pl2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_3_post_server_qc_${j}_${i}_chr${chr} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_3_post_server_qc_${j}_${i}_chr${chr} \
"unzip -P PkSw6GlBg5VGuc ${path_gwas}imputed/${i}/2022/${j}/chr_${chr}.zip -d ${path_gwas}imputed/${i}/2022/${j}/"
done 
done
done

### CHR23
for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
bsub -J"pl2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_3_post_server_qc_${j}_${i}_chr23 \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_3_post_server_qc_${j}_${i}_chr23 \
"unzip -P PkSw6GlBg5VGuc ${path_gwas}imputed/${i}/2022/${j}/chr_X.zip -d ${path_gwas}imputed/${i}/2022/${j}/"
done
done


chr=1
for j in ${ancestry[@]}
do
for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_3_post_server_qc_${j}_${i}_chr${chr}| grep -E "completed"
done
done


###################################################
# 23.6 define variants to flip or exclude by empR #
###################################################

studies=(gwas1 spain_gsa sweden_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa 
         helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa 
         mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa 
         niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa 
         pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa 
         weersma_gsa xavier_prism_gsa xavier_share_gsa)

ancestry=(eur noneur)

for i in ${studies[@]}
do
for chr in {1..23}
do bash ${path_gwas}scripts/post_imputation_summary_0.sh ${i} ${chr}
done
done

for i in ${studies[@]}
do echo ${i} && for chr in {1..23}
do echo ${chr} && tail -50 ${path_gwas}/scripts/logs/post_imputation_summary_0_stdout_${i}_${chr} | grep "completed"
done
done

for i in ${studies[@]}
do 
echo ${i} && ls -la ${path_gwas}imputed/${i}/qc/2022/chr*_list_variants_negative_EmpR_toflip | wc -l
done

for i in ${studies[@]}
do 
echo ${i} && ls -la ${path_gwas}imputed/${i}/qc/2022/chr*_list_variants_low_EmpRsq_toexclude | wc -l
done


## variants to flip:
for i in ${studies[@]}
do cat ${path_gwas}imputed/${i}/qc/2022/chr*_list_variants_negative_EmpR_toflip > ${path_gwas}imputed/${i}/qc/2022/list_variants_negative_EmpR_toflip
done

for i in ${studies[@]}
do echo ${i} && wc -l ${path_gwas}imputed/${i}/qc/2022/chr*_list_variants_negative_EmpR_toflip && wc -l ${path_gwas}imputed/${i}/qc/2022/list_variants_negative_EmpR_toflip
done


for i in ${studies[@]}
do rm ${path_gwas}imputed/${i}/qc/2022/chr*_list_variants_negative_EmpR_toflip
done

# variants to exclude
for i in ${studies[@]}
do cat ${path_gwas}imputed/${i}/qc/2022/chr*_list_variants_low_EmpRsq_toexclude > ${path_gwas}imputed/${i}/qc/2022/list_variants_low_EmpRsq_toexclude
done

for i in ${studies[@]}
do echo ${i} && wc -l ${path_gwas}imputed/${i}/qc/2022/chr*_list_variants_low_EmpRsq_toexclude && wc -l ${path_gwas}imputed/${i}/qc/2022/list_variants_low_EmpRsq_toexclude
done

for i in ${studies[@]}
do rm ${path_gwas}imputed/${i}/qc/2022/chr*_list_variants_low_EmpRsq_toexclude
done

for i in ${studies[@]}
do cat ${path_gwas}imputed/${i}/qc/2022/list_variants_low_EmpRsq_toexclude ${path_gwas}imputed/${i}/qc/2022/list_variants_to_exclude  > \
${path_gwas}imputed/${i}/qc/2022/list_variants_to_exclude_2
done

for i in ${studies[@]}
do echo ${i} \
&& wc -l ${path_gwas}imputed/${i}/qc/2022/list_variants_low_EmpRsq_toexclude \
&& wc -l ${path_gwas}imputed/${i}/qc/2022/list_variants_to_exclude \
&& wc -l ${path_gwas}imputed/${i}/qc/2022/list_variants_to_exclude_2
done


studies=(all_hce niddk_old_gwas australia_omniexome gwas2 pittsburgh_gsa
         italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa
         niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa
         belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa
         prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas
         norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas
         cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas
         swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa)


for i in ${studies[@]}
do echo ${i} && wc -l ${path_gwas}imputed/${i}/qc/list_variants_low_EmpRsq_toexclude && \
wc -l ${path_gwas}imputed/${i}/qc/2022/list_variants_to_exclude && \
wc -l ${path_gwas}imputed/${i}/qc/2022/list_variants_to_exclude_2
done

for i in ${studies[@]}
do echo ${i} && wc -l ${path_gwas}imputed/${i}/qc/list_variants_negative_EmpR_toflip
done

for i in ${studies[@]}
do cat ${path_gwas}imputed/${i}/qc/list_variants_0.5_EmpRsq_toexclude ${path_gwas}imputed/${i}/qc/2022/list_variants_to_exclude  > \
${path_gwas}imputed/${i}/qc/2022/list_variants_to_exclude_2
done

##################
# get numbers by ancestry:

## /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","swedish_uc_old_gwas"
           ,"mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

length(cohorts)
# [1] 58

ancestry<-c("eur","noneur")

for (j in 1:length(cohorts)) {
  
  print(cohorts[j])
  
  bim<-read.table(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged.bim",sep=""),head=F)
  
  files_snp_flip<-paste(path,"imputed/",cohorts[j],"/qc/2022/list_variants_negative_EmpR_toflip",sep="")
  if(file.exists(files_snp_flip)){
    snp_flip<-try(read.table(files_snp_flip),silent=TRUE)
  }else{
    files_snp_flip<-paste(path,"imputed/",cohorts[j],"/qc/list_variants_negative_EmpR_toflip",sep="")
    snp_flip<-try(read.table(files_snp_flip),silent=TRUE)
  }
  if(is.null(nrow(snp_flip))) {
    snp_flip<-as.data.frame(matrix(ncol=1,nrow=0))
  }
  
  files_snp_exclude<-paste(path,"imputed/",cohorts[j],"/qc/2022/list_variants_low_EmpRsq_toexclude",sep="")
  if(file.exists(files_snp_exclude)){
    snp_exclude<-read.table(files_snp_exclude)
  }else{
    files_snp_exclude<-paste(path,"imputed/",cohorts[j],"/qc/list_variants_0.5_EmpRsq_toexclude",sep="")
    snp_exclude<-read.table(files_snp_exclude)
  }
  
  
  eur<-read.table(paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim",sep=""),head=F)
  print(paste("N variants exclude in EUR: ",nrow(eur[which(eur$V2 %in% snp_exclude$V1),]),sep=""))
  print(paste("N variants flip in EUR: ",nrow(eur[which(eur$V2 %in% snp_flip$V1),]),sep=""))
  
  
  file_noneur<-paste(path,"pre_imputation/QC/",cohorts[j],"/",cohorts[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim",sep="")
  if(file.exists(file_noneur)) {
    noneur<-read.table(file_noneur,head=F)
    print(paste("N variants exclude in nonEUR: ",nrow(noneur[which(noneur$V2 %in% snp_exclude$V1),]),sep=""))
    print(paste("N variants flip in nonEUR: ",nrow(noneur[which(noneur$V2 %in% snp_flip$V1),]),sep=""))
    
  }
  
  rm(list=ls()[!ls() %in% c("path","j","cohorts","ancestry")]) 
  
}

##################
# [1] "australia_omniexome"
# [1] "N variants exclude in EUR: 1453"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 809"
# [1] "N variants flip in nonEUR: 0"
# [1] "gwas1"
# [1] "N variants exclude in EUR: 9308"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 2233"
# [1] "N variants flip in nonEUR: 0"
# [1] "gwas2"
# [1] "N variants exclude in EUR: 4437"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 775"
# [1] "N variants flip in nonEUR: 0"
# [1] "all_hce"
# [1] "N variants exclude in EUR: 5198"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 5190"
# [1] "N variants flip in nonEUR: 0"
# [1] "pittsburgh_gsa"
# [1] "N variants exclude in EUR: 4483"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 1331"
# [1] "N variants flip in nonEUR: 0"
# [1] "spain_gsa"
# [1] "N variants exclude in EUR: 876"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 506"
# [1] "N variants flip in nonEUR: 0"
# [1] "italy_gsa"
# [1] "N variants exclude in EUR: 10169"
# [1] "N variants flip in EUR: 3"
# [1] "N variants exclude in nonEUR: 3090"
# [1] "N variants flip in nonEUR: 3"
# [1] "kiel_austria_sibdcs_gsa"
# [1] "N variants exclude in EUR: 3292"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 2894"
# [1] "N variants flip in nonEUR: 0"
# [1] "netherlands_gsa"
# [1] "N variants exclude in EUR: 2154"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 1861"
# [1] "N variants flip in nonEUR: 0"
# [1] "slovenia_gsa"
# [1] "N variants exclude in EUR: 4298"
# [1] "N variants flip in EUR: 0"
# [1] "sweden_gsa"
# [1] "N variants exclude in EUR: 1223"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 809"
# [1] "N variants flip in nonEUR: 0"
# [1] "niddk_broad_gsa"
# [1] "N variants exclude in EUR: 1850"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 1699"
# [1] "N variants flip in nonEUR: 0"
# [1] "niddk_feinstein_gsa"
# [1] "N variants exclude in EUR: 1528"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 1482"
# [1] "N variants flip in nonEUR: 0"
# [1] "basque_gsa"
# [1] "N variants exclude in EUR: 2040"
# [1] "N variants flip in EUR: 1"
# [1] "N variants exclude in nonEUR: 722"
# [1] "N variants flip in nonEUR: 1"
# [1] "prism_nfe_gsa"
# [1] "N variants exclude in EUR: 1465"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 878"
# [1] "N variants flip in nonEUR: 0"
# [1] "lithuania_gsa"
# [1] "N variants exclude in EUR: 3468"
# [1] "N variants flip in EUR: 1"
# [1] "N variants exclude in nonEUR: 752"
# [1] "N variants flip in nonEUR: 1"
# [1] "belgium_louis_gsa"
# [1] "N variants exclude in EUR: 3874"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 1639"
# [1] "N variants flip in nonEUR: 0"
# [1] "belgium_franchimont_gsa"
# [1] "N variants exclude in EUR: 3241"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 1926"
# [1] "N variants flip in nonEUR: 0"
# [1] "belgium_vermeire_gsa"
# [1] "N variants exclude in EUR: 4014"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 3532"
# [1] "N variants flip in nonEUR: 0"
# [1] "prism_nfe_gwas"
# [1] "N variants exclude in EUR: 3415"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 3415"
# [1] "N variants flip in nonEUR: 0"
# [1] "finland_illugwas"
# [1] "N variants exclude in EUR: 358"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 185"
# [1] "N variants flip in nonEUR: 0"
# [1] "german_affy6_old_gwas"
# [1] "N variants exclude in EUR: 29701"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 13880"
# [1] "N variants flip in nonEUR: 0"
# [1] "norway_affy6_old_gwas"
# [1] "N variants exclude in EUR: 9911"
# [1] "N variants flip in EUR: 2"
# [1] "N variants exclude in nonEUR: 2052"
# [1] "N variants flip in nonEUR: 2"
# [1] "belgium_inf1_old_gwas"
# [1] "N variants exclude in EUR: 667"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 621"
# [1] "N variants flip in nonEUR: 0"
# [1] "belgium_inf2_old_gwas"
# [1] "N variants exclude in EUR: 702"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 366"
# [1] "N variants flip in nonEUR: 0"
# [1] "niddk_old_gwas"
# [1] "N variants exclude in EUR: 38"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 33"
# [1] "N variants flip in nonEUR: 0"
# [1] "cedars_370k_old_gwas"
# [1] "N variants exclude in EUR: 108"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 102"
# [1] "N variants flip in nonEUR: 0"
# [1] "cedars_610k_old_gwas"
# [1] "N variants exclude in EUR: 320"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 304"
# [1] "N variants flip in nonEUR: 0"
# [1] "cedars_omni_old_gwas"
# [1] "N variants exclude in EUR: 1121"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 1111"
# [1] "N variants flip in nonEUR: 0"
# [1] "swedish_uc_old_gwas"
# [1] "N variants exclude in EUR: 101"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 79"
# [1] "N variants flip in nonEUR: 0"
# [1] "mccauley_gsa"
# [1] "N variants exclude in EUR: 1193"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 429"
# [1] "N variants flip in nonEUR: 0"
# [1] "ccfa_gsa"
# [1] "N variants exclude in EUR: 1847"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 1506"
# [1] "N variants flip in nonEUR: 0"
# [1] "cedars_gsa"
# [1] "N variants exclude in EUR: 2336"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 1726"
# [1] "N variants flip in nonEUR: 0"
# [1] "bernstein_gsa"
# [1] "N variants exclude in EUR: 1996"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 1153"
# [1] "N variants flip in nonEUR: 0"
# [1] "farkkila_gsa"
# [1] "N variants exclude in EUR: 3023"
# [1] "N variants flip in EUR: 4"
# [1] "franchimont_gsa"
# [1] "N variants exclude in EUR: 3750"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 3048"
# [1] "N variants flip in nonEUR: 0"
# [1] "franke_gsa"
# [1] "N variants exclude in EUR: 1996"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 796"
# [1] "N variants flip in nonEUR: 0"
# [1] "helmsley_prism_gsa"
# [1] "N variants exclude in EUR: 825"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 823"
# [1] "N variants flip in nonEUR: 0"
# [1] "helmsley_xavier_prism_gsa"
# [1] "N variants exclude in EUR: 522"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 522"
# [1] "N variants flip in nonEUR: 0"
# [1] "hyams_protect_gsa"
# [1] "N variants exclude in EUR: 1866"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 1461"
# [1] "N variants flip in nonEUR: 0"
# [1] "lewis_sparc_gsa"
# [1] "N variants exclude in EUR: 2051"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 1712"
# [1] "N variants flip in nonEUR: 0"
# [1] "mccauley_new_gsa"
# [1] "N variants exclude in EUR: 1242"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 984"
# [1] "N variants flip in nonEUR: 0"
# [1] "mcgovern_gsa"
# [1] "N variants exclude in EUR: 2771"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 2215"
# [1] "N variants flip in nonEUR: 0"
# [1] "moayyedi_imagine_gsa"
# [1] "N variants exclude in EUR: 1440"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 1175"
# [1] "N variants flip in nonEUR: 0"
# [1] "newberry_share_gsa"
# [1] "N variants exclude in EUR: 1006"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 762"
# [1] "N variants flip in nonEUR: 0"
# [1] "niddk_cho_gsa"
# [1] "N variants exclude in EUR: 1711"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 1268"
# [1] "N variants flip in nonEUR: 0"
# [1] "niddk_duerr_gsa"
# [1] "N variants exclude in EUR: 2297"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 1480"
# [1] "N variants flip in nonEUR: 0"
# [1] "niddk_rioux_gsa"
# [1] "N variants exclude in EUR: 2111"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 1248"
# [1] "N variants flip in nonEUR: 0"
# [1] "niddk_silverberg_gsa"
# [1] "N variants exclude in EUR: 2311"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 1776"
# [1] "N variants flip in nonEUR: 0"
# [1] "palotie_hus_gsa"
# [1] "N variants exclude in EUR: 2044"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 765"
# [1] "N variants flip in nonEUR: 0"
# [1] "pekow_share_gsa"
# [1] "N variants exclude in EUR: 1161"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 896"
# [1] "N variants flip in nonEUR: 0"
# [1] "rioux_igenomed_gsa"
# [1] "N variants exclude in EUR: 2766"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 933"
# [1] "N variants flip in nonEUR: 0"
# [1] "sands_msccr_gsa"
# [1] "N variants exclude in EUR: 1135"
# [1] "N variants flip in EUR: 1"
# [1] "N variants exclude in nonEUR: 988"
# [1] "N variants flip in nonEUR: 1"
# [1] "stampfer_gsa"
# [1] "N variants exclude in EUR: 2714"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 1360"
# [1] "N variants flip in nonEUR: 0"
# [1] "vermeire_gsa"
# [1] "N variants exclude in EUR: 4093"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 2835"
# [1] "N variants flip in nonEUR: 0"
# [1] "weersma_gsa"
# [1] "N variants exclude in EUR: 1568"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 853"
# [1] "N variants flip in nonEUR: 0"
# [1] "xavier_prism_gsa"
# [1] "N variants exclude in EUR: 1320"
# [1] "N variants flip in EUR: 1"
# [1] "N variants exclude in nonEUR: 900"
# [1] "N variants flip in nonEUR: 1"
# [1] "xavier_share_gsa"
# [1] "N variants exclude in EUR: 1512"
# [1] "N variants flip in EUR: 0"
# [1] "N variants exclude in nonEUR: 1001"
# [1] "N variants flip in nonEUR: 0"
##################


###########################################
# 23.7 Convert ped/map files to VCF files #
###########################################


studies=(gwas1 spain_gsa sweden_gsa bernstein_gsa farkkila_gsa franchimont_gsa 
         franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa 
         lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa 
         newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa 
         niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa 
         sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)
path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=1000

ancestry=(eur noneur)

for chr in {1..23}
do
for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_4_post_server_qc_${j}_${i}_chr${chr} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_4_post_server_qc_${j}_${i}_chr${chr} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_${j}_nomonom \
--allow-no-sex \
--keep-allele-order \
--flip ${path_gwas}imputed/${i}/qc/2022/list_variants_negative_EmpR_toflip \
--exclude ${path_gwas}imputed/${i}/qc/2022/list_variants_to_exclude_2 \
--chr ${chr} \
--output-chr chr26 --recode vcf-iid --out ${path_gwas}imputation_ready/${i}/${i}_${j}_chr${chr}_2022"
done
done
done

# to the same with old studies pointing to old flip file
studies=(all_hce niddk_old_gwas australia_omniexome gwas2 pittsburgh_gsa
         italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa
         niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa
         belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa
         prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas
         norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas
         cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas
         swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa)

for chr in {1..23}
do
for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_4_post_server_qc_${j}_${i}_chr${chr} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_4_post_server_qc_${j}_${i}_chr${chr} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_${j}_nomonom \
--allow-no-sex \
--keep-allele-order \
--flip ${path_gwas}imputed/${i}/qc/list_variants_negative_EmpR_toflip \
--exclude ${path_gwas}imputed/${i}/qc/2022/list_variants_to_exclude_2 \
--chr ${chr} \
--output-chr chr26 --recode vcf-iid --out ${path_gwas}imputation_ready/${i}/${i}_${j}_chr${chr}_2022"
done
done
done

# studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa 
#          netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa 
#          belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas 
#          finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas 
#          belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas 
#          swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

studies=(all_hce niddk_old_gwas australia_omniexome gwas2 pittsburgh_gsa
         italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa
         niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa
         belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa
         prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas
         norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas
         cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas
         swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa)



for i in ${studies[@]}
do
echo ${i} && for j in ${ancestry[@]}
do
echo ${j} && tail -40 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_4_post_server_qc_${j}_${i}_chr${chr} | grep -E "completed"
done
done

for i in ${studies[@]}
do ls -la ${path_gwas}imputation_ready/${i}/${i}_${j}_chr${chr}_2022.vcf
done


########################################################
# 23.8 Create a sorted vcf.gz - final round imputation #
########################################################


# studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

studies=(all_hce niddk_old_gwas australia_omniexome gwas2 pittsburgh_gsa
         italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa
         niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa
         belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa
         prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas
         norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas
         cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas
         swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa gwas1)
ancestry=(eur noneur)

for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
for chr in {1..23}
do
bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_5_post_server_qc_${j}_${i}_chr${chr} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_5_post_server_qc_${j}_${i}_chr${chr} \
"/path/to/software/bcftools-1.16/./bcftools sort ${path_gwas}imputation_ready/${i}/${i}_${j}_chr${chr}_2022.vcf \
-Oz -o ${path_gwas}imputation_ready/${i}/${i}_${j}_chr${chr}_2022.vcf.gz"
done
done
done

# rm intermediate files:
for i in ${studies[@]}
do rm ${path_gwas}imputation_ready/${i}/${i}_*_chr*.vcf
done
for i in ${studies[@]}
do rm ${path_gwas}imputation_ready/${i}/*.hh
done
for i in ${studies[@]}
do rm ${path_gwas}imputation_ready/${i}/*.log
done
for i in ${studies[@]}
do rm ${path_gwas}imputation_ready/${i}/*.nosex
done

for i in ${studies[@]}
do for j in ${ancestry[@]}
do echo ${i} && ls -la ${path_gwas}imputation_ready/${i}/${i}_${j}_chr*_2022.vcf.gz |  wc -l
done
done

for i in ${studies[@]}
do for j in ${ancestry[@]}
do echo ${i} && ls -la ${path_gwas}imputation_ready/${i}/${i}_${j}_chr1_2022.vcf.gz
done
done


## get total number of variants per final file:
for i in ${studies[@]}
do
rm ${path_gwas}imputation_ready/test/count_n_variants_${i}_*
done

for i in ${studies[@]}
do 
for j in ${ancestry[@]}
do
for chr in {1..23}
do
/path/to/software/bcftools-1.16/./bcftools stats ${path_gwas}imputation_ready/${i}/${i}_${j}_chr${chr}_2022.vcf.gz | grep "number of records" | \
cut -f4 | awk 'NR>1' >> ${path_gwas}imputation_ready/test/count_n_variants_${i}_${j}
done
done
done

###### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","swedish_uc_old_gwas"
           ,"mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

length(cohorts)
# [1] 58

ancestry<-c("eur","noneur")

for (i in 1:length(cohorts)) {
  
  print(cohorts[i])
  
  for (j in 1:length(ancestry)) {
    print(ancestry[j])
    a<-tryCatch(read.table(paste(path,"imputation_ready/test/count_n_variants_",cohorts[i],"_",ancestry[j],sep=""),head=F), error=function(e) NULL)
    if(!is.null(a)) {
      print(sum(a$V1))
    } else {
      print("No variants")
    }

  }
  print("           ")
  rm(a)
}

#############################
# [1] "australia_omniexome"
# [1] "eur"
# [1] 669551
# [1] "noneur"
# [1] 653177
# [1] "           "
# [1] "gwas1"
# [1] "eur"
# [1] 427055
# [1] "noneur"
# [1] 359278
# [1] "           "
# [1] "gwas2"
# [1] "eur"
# [1] 752046
# [1] "noneur"
# [1] 693961
# [1] "           "
# [1] "all_hce"
# [1] "eur"
# [1] 253045
# [1] "noneur"
# [1] 257243
# [1] "           "
# [1] "pittsburgh_gsa"
# [1] "eur"
# [1] 859591
# [1] "noneur"
# [1] 741763
# [1] "           "
# [1] "spain_gsa"
# [1] "eur"
# [1] 565580
# [1] "noneur"
# [1] 523752
# [1] "           "
# [1] "italy_gsa"
# [1] "eur"
# [1] 478301
# [1] "noneur"
# [1] 286577
# [1] "           "
# [1] "kiel_austria_sibdcs_gsa"
# [1] "eur"
# [1] 454949
# [1] "noneur"
# [1] 452229
# [1] "           "
# [1] "netherlands_gsa"
# [1] "eur"
# [1] 476837
# [1] "noneur"
# [1] 475696
# [1] "           "
# [1] "slovenia_gsa"
# [1] "eur"
# [1] 463105
# [1] "noneur"
# [1] "No variants"
# [1] "           "
# [1] "sweden_gsa"
# [1] "eur"
# [1] 474118
# [1] "noneur"
# [1] 410797
# [1] "           "
# [1] "niddk_broad_gsa"
# [1] "eur"
# [1] 480273
# [1] "noneur"
# [1] 479360
# [1] "           "
# [1] "niddk_feinstein_gsa"
# [1] "eur"
# [1] 474867
# [1] "noneur"
# [1] 477813
# [1] "           "
# [1] "basque_gsa"
# [1] "eur"
# [1] 488406
# [1] "noneur"
# [1] 346900
# [1] "           "
# [1] "prism_nfe_gsa"
# [1] "eur"
# [1] 465896
# [1] "noneur"
# [1] 403050
# [1] "           "
# [1] "lithuania_gsa"
# [1] "eur"
# [1] 472001
# [1] "noneur"
# [1] 332575
# [1] "           "
# [1] "belgium_louis_gsa"
# [1] "eur"
# [1] 474878
# [1] "noneur"
# [1] 342079
# [1] "           "
# [1] "belgium_franchimont_gsa"
# [1] "eur"
# [1] 473238
# [1] "noneur"
# [1] 413262
# [1] "           "
# [1] "belgium_vermeire_gsa"
# [1] "eur"
# [1] 475896
# [1] "noneur"
# [1] 475524
# [1] "           "
# [1] "prism_nfe_gwas"
# [1] "eur"
# [1] 234933
# [1] "noneur"
# [1] 236178
# [1] "           "
# [1] "finland_illugwas"
# [1] "eur"
# [1] 230208
# [1] "noneur"
# [1] 205910
# [1] "           "
# [1] "german_affy6_old_gwas"
# [1] "eur"
# [1] 729844
# [1] "noneur"
# [1] 668972
# [1] "           "
# [1] "norway_affy6_old_gwas"
# [1] "eur"
# [1] 696669
# [1] "noneur"
# [1] 447207
# [1] "           "
# [1] "belgium_inf1_old_gwas"
# [1] "eur"
# [1] 297400
# [1] "noneur"
# [1] 288013
# [1] "           "
# [1] "belgium_inf2_old_gwas"
# [1] "eur"
# [1] 285920
# [1] "noneur"
# [1] 164935
# [1] "           "
# [1] "niddk_old_gwas"
# [1] "eur"
# [1] 296499
# [1] "noneur"
# [1] 294374
# [1] "           "
# [1] "cedars_370k_old_gwas"
# [1] "eur"
# [1] 332036
# [1] "noneur"
# [1] 332758
# [1] "           "
# [1] "cedars_610k_old_gwas"
# [1] "eur"
# [1] 568736
# [1] "noneur"
# [1] 572323
# [1] "           "
# [1] "cedars_omni_old_gwas"
# [1] "eur"
# [1] 585711
# [1] "noneur"
# [1] 597174
# [1] "           "
# [1] "swedish_uc_old_gwas"
# [1] "eur"
# [1] 292950
# [1] "noneur"
# [1] 239937
# [1] "           "
# [1] "mccauley_gsa"
# [1] "eur"
# [1] 477342
# [1] "noneur"
# [1] 376449
# [1] "           "
# [1] "ccfa_gsa"
# [1] "eur"
# [1] 479094
# [1] "noneur"
# [1] 475801
# [1] "           "
# [1] "cedars_gsa"
# [1] "eur"
# [1] 478990
# [1] "noneur"
# [1] 478315
# [1] "           "
# [1] "bernstein_gsa"
# [1] "eur"
# [1] 408001
# [1] "noneur"
# [1] 336502
# [1] "           "
# [1] "farkkila_gsa"
# [1] "eur"
# [1] 428819
# [1] "noneur"
# [1] "No variants"
# [1] "           "
# [1] "franchimont_gsa"
# [1] "eur"
# [1] 475399
# [1] "noneur"
# [1] 473029
# [1] "           "
# [1] "franke_gsa"
# [1] "eur"
# [1] 474496
# [1] "noneur"
# [1] 330722
# [1] "           "
# [1] "helmsley_prism_gsa"
# [1] "eur"
# [1] 237607
# [1] "noneur"
# [1] 238957
# [1] "           "
# [1] "helmsley_xavier_prism_gsa"
# [1] "eur"
# [1] 238652
# [1] "noneur"
# [1] 239300
# [1] "           "
# [1] "hyams_protect_gsa"
# [1] "eur"
# [1] 474800
# [1] "noneur"
# [1] 453498
# [1] "           "
# [1] "lewis_sparc_gsa"
# [1] "eur"
# [1] 479720
# [1] "noneur"
# [1] 480424
# [1] "           "
# [1] "mccauley_new_gsa"
# [1] "eur"
# [1] 481384
# [1] "noneur"
# [1] 480928
# [1] "           "
# [1] "mcgovern_gsa"
# [1] "eur"
# [1] 480322
# [1] "noneur"
# [1] 483357
# [1] "           "
# [1] "moayyedi_imagine_gsa"
# [1] "eur"
# [1] 477498
# [1] "noneur"
# [1] 457525
# [1] "           "
# [1] "newberry_share_gsa"
# [1] "eur"
# [1] 476835
# [1] "noneur"
# [1] 451803
# [1] "           "
# [1] "niddk_cho_gsa"
# [1] "eur"
# [1] 478691
# [1] "noneur"
# [1] 464548
# [1] "           "
# [1] "niddk_duerr_gsa"
# [1] "eur"
# [1] 477968
# [1] "noneur"
# [1] 439240
# [1] "           "
# [1] "niddk_rioux_gsa"
# [1] "eur"
# [1] 476716
# [1] "noneur"
# [1] 411207
# [1] "           "
# [1] "niddk_silverberg_gsa"
# [1] "eur"
# [1] 479082
# [1] "noneur"
# [1] 470980
# [1] "           "
# [1] "palotie_hus_gsa"
# [1] "eur"
# [1] 471792
# [1] "noneur"
# [1] 364674
# [1] "           "
# [1] "pekow_share_gsa"
# [1] "eur"
# [1] 475925
# [1] "noneur"
# [1] 442332
# [1] "           "
# [1] "rioux_igenomed_gsa"
# [1] "eur"
# [1] 468358
# [1] "noneur"
# [1] 329256
# [1] "           "
# [1] "sands_msccr_gsa"
# [1] "eur"
# [1] 478163
# [1] "noneur"
# [1] 480530
# [1] "           "
# [1] "stampfer_gsa"
# [1] "eur"
# [1] 474703
# [1] "noneur"
# [1] 381192
# [1] "           "
# [1] "vermeire_gsa"
# [1] "eur"
# [1] 475923
# [1] "noneur"
# [1] 457700
# [1] "           "
# [1] "weersma_gsa"
# [1] "eur"
# [1] 473561
# [1] "noneur"
# [1] 377280
# [1] "           "
# [1] "xavier_prism_gsa"
# [1] "eur"
# [1] 475370
# [1] "noneur"
# [1] 438553
# [1] "           "
# [1] "xavier_share_gsa"
# [1] "eur"
# [1] 475786
# [1] "noneur"
# [1] 428075
# [1] "           "
#############################


  
##################################################################################################################################################################

####################################################
# 23.9 submit second and FINAL round of imputation #
####################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
ancestry=(eur)

studies=()


i=niddk_broad_gsa
for j in ${ancestry[@]}
do echo ${i} && ls -la ${path_gwas}imputation_ready/${i}/${i}_${j}_chr*_2022.vcf.gz |  wc -l
done
for j in ${ancestry[@]}
do echo ${j} && python2 /path/to/user/scripts/IIBDGC/topmed_imputation_${j}.py ${i}
done



#############
# downloading

# i=swedish_uc_old_gwas
# # eur
# # swedish_uc_old_gwas
# # QC  swedish_uc_old_gwas job_eur=job-20221123-125735-493
# # noneur
# # swedish_uc_old_gwas
# # QC  swedish_uc_old_gwas job_noneur=job-20221123-125739-425 - downloaded
# job_eur=job-20221123-125735-493
# job_noneur=job-20221123-125739-425

# i=german_affy6_old_gwas
# # eur
# # german_affy6_old_gwas
# # QC  german_affy6_old_gwas job_eur=job-20221123-074250-648 - downloaded
# # noneur
# # german_affy6_old_gwas
# # QC  german_affy6_old_gwas job_noneur=job-20221123-090813-725 - downloaded
# job_eur=job-20221123-074250-648
# job_noneur=job-20221123-090813-725


# i=mccauley_gsa
# # eur
# # mccauley_gsa
# # QC  mccauley_gsa job_eur=job-20221123-113202-478 - downloaded
# # noneur
# # mccauley_gsa
# # QC  mccauley_gsa job_noneur=job-20221123-115341-593 - downloaded
# job_eur=job-20221123-113202-478
# job_noneur=job-20221123-115341-593


# i=norway_affy6_old_gwas
# # eur
# # norway_affy6_old_gwas
# # QC  norway_affy6_old_gwas job_eur=job-20221123-102046-753 - downloaded
# # noneur
# # norway_affy6_old_gwas
# # QC  norway_affy6_old_gwas job_noneur=job-20221123-102110-372 - downloaded
# job_eur=job-20221123-102046-753
# job_noneur=job-20221123-102110-372


# i=prism_nfe_gwas
# # eur
# # prism_nfe_gwas
# # QC  prism_nfe_gwas job_eur=job-20221123-070624-730 - downloaded
# # noneur
# # prism_nfe_gwas
# # QC  prism_nfe_gwas job_noneur=job-20221123-074212-532 - downloaded
# job_eur=job-20221123-070624-730
# job_noneur=job-20221123-074212-532

# i=niddk_broad_gsa
# # eur
# # niddk_broad_gsa
# # QC  niddk_broad_gsa job_eur=job-20221122-200603-312 - downloaded
# # noneur
# # niddk_broad_gsa
# # QC  niddk_broad_gsa job_noneur=job-20221122-200633-187 - downloaded
# job_eur=job-20221122-200603-312
# job_noneur=job-20221122-200633-187


# i=prism_nfe_gsa
# # eur
# # prism_nfe_gsa
# # QC  prism_nfe_gsa job_eur=job-20221123-040106-395 - downloaded
# # noneur
# # prism_nfe_gsa
# # QC  prism_nfe_gsa job_noneur=job-20221123-070556-654 - downloaded
# job_eur=job-20221123-040106-395
# job_noneur=job-20221123-070556-654

# i=belgium_louis_gsa
# # eur
# # belgium_louis_gsa
# # QC  belgium_louis_gsa job_eur=job-20221123-022600-353 - downloaded
# # noneur
# # belgium_louis_gsa
# # QC  belgium_louis_gsa job_noneur=job-20221123-022606-563 - downloaded
# job_eur=job-20221123-022600-353
# job_noneur=job-20221123-022606-563


# i=netherlands_gsa
# # eur
# # netherlands_gsa
# # QC  netherlands_gsa job_eur=job-20221122-115945-718 - downloaded
# # noneur
# # netherlands_gsa
# # QC  netherlands_gsa job_noneur=job-20221122-160612-194 - downloaded
# job_eur=job-20221122-115945-718
# job_noneur=job-20221122-160612-194

# i=basque_gsa
# # eur
# # basque_gsa
# # QC  basque_gsa job_eur=job-20221122-214524-622 - downloaded
# # noneur
# # basque_gsa
# # QC  basque_gsa job_noneur=job-20221122-214657-285 - downloaded
# job_eur=job-20221122-214524-622
# job_noneur=job-20221122-214657-285

# i=cedars_omni_old_gwas
# # eur
# # cedars_omni_old_gwas
# # QC  cedars_omni_old_gwas job_eur=job-20221122-170556-790 - downloaded
# # noneur
# # cedars_omni_old_gwas
# # QC  cedars_omni_old_gwas job_noneur=job-20221122-180810-808 - downloaded
# job_eur=job-20221122-170556-790
# job_noneur=job-20221122-180810-808

# i=cedars_610k_old_gwas
# # eur
# # cedars_610k_old_gwas
# # QC  cedars_610k_old_gwas job_eur=job-20221122-084839-520 - downloaded
# # noneur
# # QC  cedars_610k_old_gwas job_noneur=job-20221122-143412-620 - downloaded
# # job_eur=job-20221122-084839-520
# # job_noneur=job-20221122-143412-620

# i=gwas2
# # eur
# # gwas2
# # QC  gwas2 job_eur=job-20221121-225313-998 - downloaded
# # noneur
# # gwas2
# # QC gwas 2 job_noneur=job-20221121-215430-097 - downloaded
# job_eur=job-20221121-225313-998
# job_noneur=job-20221121-215430-097


# i=finland_illugwas
# # eur
# # finland_illugwas
# # QC  finland_illugwas job_eur=job-20221119-191237-859 - downloaded
# # noneur
# # finland_illugwas
# # QC  finland_illugwas job_noneur=job-20221122-132503-256 - downloaded
# # job_eur=job-20221119-191237-859
# job_noneur=job-20221122-132503-256


# i=cedars_370k_old_gwas
# # eur
# # cedars_370k_old_gwas
# # QC  cedars_370k_old_gwas job_eur=job-20221122-084226-820 - downloaded
# # noneur
# # belgium_franchimont_gsa
# # QC  belgium_franchimont_gsa job_noneur=job-20221122-070654-150 - downloaded
# job_eur=job-20221122-084226-820
# job_noneur=job-20221122-070654-150

# i=vermeire_gsa
# # eur
# # vermeire_gsa
# # QC  vermeire_gsa job_eur=job-20221121-190312-351 - downloaded
# # noneur
# # vermeire_gsa
# # QC  vermeire_gsa job_noneur=job-20221121-195641-005 - downloaded
# job_eur=job-20221121-190312-351
# job_noneur=job-20221121-195641-005

# i=belgium_franchimont_gsa
# # eur
# # belgium_franchimont_gsa
# # QC  belgium_franchimont_gsa job_eur=job-20221122-020344-367 - downloaded
# # noneur
# # belgium_franchimont_gsa
# # QC  belgium_franchimont_gsa job_noneur=job-20221122-065644-700 - downloaded
# job_eur=job-20221122-020344-367
# job_noneur=job-20221122-065644-700

# i=lithuania_gsa
# # eur
# # lithuania_gsa
# # QC  lithuania_gsa job_eur=job-20221121-205839-354 - downloaded
# # noneur
# # lithuania_gsa
# # QC  lithuania_gsa job_noneur=job-20221121-205905-324 - downloaded
# job_eur=job-20221121-205839-354
# job_noneur=job-20221121-205905-324

# i=gwas1
# # eur
# # gwas1
# # QC  gwas1 job_eur=job-20221120-201723-720 - downloaded
# # noneur
# # gwas1
# # QC  gwas1 job_noneur=job-20221121-190242-381 - downloaded
# job_eur=job-20221120-201723-720
# job_noneur=job-20221121-190242-381

# i=pittsburgh_gsa
# # eur
# # pittsburgh_gsa
# # QC  pittsburgh_gsa job_eur=job-20221120-122403-793 - downloaded
# # noneur
# # pittsburgh_gsa
# # QC  pittsburgh_gsa job_noneur=job-20221120-122416-947 - downloaded
# job_eur=job-20221120-122403-793
# job_noneur=job-20221120-122416-947

# i=italy_gsa
# # eur
# # italy_gsa
# # QC  italy_gsa job_eur=job-20221120-201640-524 - downloaded
# # noneur
# # italy_gsa
# # QC  pittsburgh_gsa job_noneur=job-20221120-122943-909 - downloaded
# job_eur=job-20221120-201640-524
# job_noneur=job-20221120-122943-909

# i=belgium_vermeire_gsa
# # eur
# # belgium_vermeire_gsa
# # QC  belgium_vermeire_gsa job_eur=job-20221119-230306-840 - downloaded
# # noneur
# # belgium_vermeire_gsa
# # QC  belgium_vermeire_gsa job_noneur=job-20221119-230318-275 - downloaded
# job_eur=job-20221119-230306-840
# job_noneur=job-20221119-230318-275

# i=belgium_inf1_old_gwas
# # eur
# # belgium_inf1_old_gwas
# # QC  belgium_inf1_old_gwas job_eur=job-20221119-203119-528 - downloaded
# # noneur
# # belgium_inf1_old_gwas
# # QC  belgium_inf1_old_gwas job_noneur=job-20221119-203259-294 - downloaded
# job_eur=job-20221119-203119-528
# job_noneur=job-20221119-203259-294

# i=cedars_gsa
# # eur
# # cedars_gsa
# # QC  cedars_gsa job_eur=job-20221119-162511-748 -downloaded
# # noneur
# # cedars_gsa
# # QC  cedars_gsa job_noneur=job-20221119-175420-201 - downloaded
# job_eur=job-20221119-162511-748
# job_noneur=job-20221119-175420-201


# i=all_hce
# # eur
# # all_hce
# # QC  all_hce job_eur=job-20221117-190946-778 - DOWNLOADED
# # noneur
# # all_hce
# # QC  all_hce job_noneur=job-20221117-195209-577 - DOWNLOADED
# job_noneur=job-20221117-195209-577
# job_eur=job-20221117-190946-778


# i=belgium_inf2_old_gwas
# # eur
# # belgium_inf2_old_gwas
# # QC  belgium_inf2_old_gwas job_eur=job-20221119-214809-174 - downloaded
# # noneur
# # belgium_inf2_old_gwas
# # QC  belgium_inf2_old_gwas job_noneur=job-20221119-214821-161 - downloaded
# job_eur=job-20221119-214809-174
# job_noneur=job-20221119-214821-161

# i=slovenia_gsa
# # eur
# # slovenia_gsa
# # QC  slovenia_gsa job_eur=job-20221119-191206-168 - downloaded
# job_eur=job-20221119-191206-168


# i=ccfa_gsa
# # eur
# # ccfa_gsa
# # QC  ccfa_gsa job_eur=job-20221119-152544-508 - downloaded
# # noneur
# # ccfa_gsa
# # QC  ccfa_gsa job_noneur=job-20221119-161733-605 - downloaded
# job_eur=job-20221119-152544-508
# job_noneur=job-20221119-161733-605

# i=australia_omniexome
# # eur
# # australia_omniexome
# # QC  australia_omniexome job_eur=job-20221119-110203-777 - DOWNLOADED
# # noneur
# # australia_omniexome
# # QC  australia_omniexome job_noneur=job-20221119-145507-438 - DOWNLOADED
# job_eur=job-20221119-110203-777
# job_noneur=job-20221119-145507-438

# i=niddk_feinstein_gsa
# # eur
# # niddk_feinstein_gsa
# # QC  niddk_feinstein_gsa job_eur=job-20221118-215831-792  - DOWNLOADED
# # noneur
# # niddk_feinstein_gsa
# # QC  niddk_feinstein_gsa job_noneur=job-20221119-072834-286 - DOWNLOADED
# # job_eur=job-20221118-215831-792
# # job_noneur=job-20221119-072834-286


# i=kiel_austria_sibdcs_gsa
# # eur
# # kiel_austria_sibdcs_gsa
# # QC  kiel_austria_sibdcs_gsa job_eur=job-20221117-195248-602 - DOWNLOADED
# # noneur
# # kiel_austria_sibdcs_gsa
# # QC  kiel_austria_sibdcs_gsa job_noneur=job-20221118-032544-186 - DOWNLOADED
# job_eur=job-20221117-195248-602
# job_noneur=job-20221118-032544-186

# i=niddk_old_gwas
# # eur
# # niddk_old_gwas
# # QC  niddk_old_gwas job_eur=job-20221118-152004-114 - DOWNLOADED
# # noneur
# # niddk_old_gwas
# # QC  niddk_old_gwas job_noneur=job-20221119-072747-453  - DOWNLOADED
# job_eur=job-20221118-152004-114
# job_noneur=job-20221119-072747-453

# i=weersma_gsa - donwloaded
# # eur
# # weersma_gsa
# # QC  weersma_gsa job_eur=job-20221114-113120-312
# # noneur
# # weersma_gsa
# # QC  weersma_gsa job_noneur=job-20221114-135016-862
# job_eur=job-20221114-113120-312
# job_noneur=job-20221114-135016-862


# i=swedish_uc_old_gwas
# # eur
# # swedish_uc_old_gwas
# # QC  swedish_uc_old_gwas job_eur=job-20221112-210359-349
# # noneur
# # swedish_uc_old_gwas
# # QC  swedish_uc_old_gwas job_noneur=job-20221112-210403-261
# job_eur=job-20221112-210359-349
# job_noneur=job-20221112-210403-261


# i=ccfa_gsa - downloaded
# # eur
# # ccfa_gsa
# # QC  ccfa_gsa job_eur=job-20221112-155239-418
# # noneur
# # ccfa_gsa
# # QC  ccfa_gsa job_noneur=job-20221112-155244-825
# job_eur=job-20221112-155239-418
# job_noneur=job-20221112-155244-825

# i=cedars_gsa
# # eur
# # cedars_gsa
# # QC  cedars_gsa job_eur=job-20221112-083227-378
# # noneur
# # cedars_gsa
# # QC  cedars_gsa job_noneur=job-20221112-155213-777
# job_eur=job-20221112-083227-378
# job_noneur=job-20221112-155213-777

# i=cedars_610k_old_gwas - downloaded
# # eur
# # cedars_610k_old_gwas
# # QC  cedars_610k_old_gwas job_eur=job-20221112-082854-637
# # noneur
# # cedars_610k_old_gwas
# # QC  cedars_610k_old_gwas job_noneur=job-20221112-082859-105
# job_eur=job-20221112-082854-637
# job_noneur=job-20221112-082859-105

# i=cedars_omni_old_gwas - downloaded
# # eur
# # cedars_omni_old_gwas
# # QC  cedars_omni_old_gwas job_eur=job-20221111-204521-713
# # noneur
# # cedars_omni_old_gwas
# # QC  cedars_omni_old_gwas job_noneur=job-20221111-204538-104
# job_eur=job-20221111-204521-713
# job_noneur=job-20221111-204538-104

# i=cedars_370k_old_gwas - downloaded
# # eur
# # cedars_370k_old_gwas
# # QC  cedars_370k_old_gwas job_eur=job-20221111-172118-270
# # noneur
# # cedars_370k_old_gwas
# # QC  cedars_370k_old_gwas job_noneur=job-20221111-204424-168
# job_eur=job-20221111-172118-270
# job_noneur=job-20221111-204424-168

# i=belgium_inf2_old_gwas - downloaded
# # eur
# # belgium_inf2_old_gwas
# # QC  belgium_inf2_old_gwas job_eur=job-20221111-172040-089
# # noneur
# # belgium_inf2_old_gwas
# # QC  belgium_inf2_old_gwas job_noneur=job-20221111-172043-641
# job_eur=job-20221111-172040-089
# job_noneur=job-20221111-172043-641


# i=mccauley_gsa - downloaded
# # eur
# # mccauley_gsa
# # QC  mccauley_gsa job_eur=job-20221111-101140-993
# # noneur
# # mccauley_gsa
# # QC  mccauley_gsa job_noneur=job-20221111-101145-892
# job_eur=job-20221111-101140-993
# job_noneur=job-20221111-101145-892

# i=belgium_inf1_old_gwas - downloaded
# # eur
# # belgium_inf1_old_gwas
# # QC  belgium_inf1_old_gwas job_eur=job-20221111-052614-183
# # noneur
# # belgium_inf1_old_gwas
# # QC  belgium_inf1_old_gwas job_noneur=job-20221111-101222-115
# job_noneur=job-20221111-101222-115
# job_eur=job-20221111-052614-183


# i=norway_affy6_old_gwas - downloaded
# # eur
# # norway_affy6_old_gwas
# # QC  norway_affy6_old_gwas job_eur=job-20221110-215501-310
# # noneur
# # norway_affy6_old_gwas
# # QC  norway_affy6_old_gwas job_noneur=job-20221111-052544-435
# job_eur=job-20221110-215501-310
# job_noneur=job-20221111-052544-435

# i=belgium_franchimont_gsa - downloaded
# # eur
# # belgium_franchimont_gsa
# # QC  belgium_franchimont_gsa job_eur=job-20221110-175048-544
# # noneur
# # belgium_franchimont_gsa
# # QC  belgium_franchimont_gsa job_noneur=job-20221111-052558-105
# job_eur=job-20221110-175048-544
# job_noneur=job-20221111-052558-105

# i=belgium_louis_gsa - downloaded
# # eur
# # belgium_louis_gsa
# # QC  belgium_louis_gsa job_eur=job-20221110-174921-994
# # noneur
# # belgium_louis_gsa
# # QC  belgium_louis_gsa job_noneur=job-20221110-155957-434
# job_eur=job-20221110-174921-994
# job_noneur=job-20221110-155957-434

# i=lithuania_gsa - downloaded
# # eur
# # lithuania_gsa
# # QC  lithuania_gsa job_eur=job-20221110-155805-916
# # noneur
# # lithuania_gsa
# # QC  lithuania_gsa job_noneur=job-20221110-155812-341
# job_eur=job-20221110-155805-916
# job_noneur=job-20221110-155812-341

# i=sweden_gsa - downloaded
# # eur
# # sweden_gsa
# # QC  sweden_gsa job_eur=job-20221110-075233-749
# # noneur
# # sweden_gsa
# # QC  sweden_gsa job_noneur=job-20221110-075238-219
# job_eur=job-20221110-075233-749
# job_noneur=job-20221110-075238-219

# i=italy_gsa - downloaded
# # eur
# # italy_gsa
# # QC  italy_gsa job_eur=job-20221109-225258-710
# # noneur
# # italy_gsa
# # QC italy_gsa job_noneur=job-20221110-075153-750
# job_eur=job-20221109-225258-710
# job_noneur=job-20221110-075153-750


# i=basque_gsa - downloaded
# # eur
# # basque_gsa
# # QC  basque_gsa job_eur=job-20221109-171117-379
# # noneur
# # basque_gsa
# # QC  basque_gsa job_noneur=job-20221109-225226-121
# job_eur=job-20221109-171117-379
# job_noneur=job-20221109-225226-121

# i=spain_gsa - downloaded
# # eur
# # spain_gsa
# # QC  spain_gsa job_eur=job-20221109-171015-160
# # noneur
# # spain_gsa
# # QC  spain_gsa job_noneur=job-20221109-171048-430
# job_eur=job-20221109-171015-160
# job_noneur=job-20221109-171048-430

# i=pittsburgh_gsa - donwloaded
# # eur
# # pittsburgh_gsa
# # QC  pittsburgh_gsa job_eur=job-20221109-062023-412
# # noneur
# # pittsburgh_gsa
# # QC  pittsburgh_gsa job_noneur=job-20221109-070829-116
# job_eur=job-20221109-062023-412
# job_noneur=job-20221109-070829-116

# i=gwas1 - downloaded
# # eur
# # gwas1
# # QC  gwas1 job_eur=job-20221109-002032-881
# # noneur
# # gwas1
# # QC  gwas1 job_noneur=job-20221109-062002-482
# job_eur=job-20221109-002032-881
# job_noneur=job-20221109-062002-482

# i=german_affy6_old_gwas
# # eur
# # german_affy6_old_gwas
# # QC  german_affy6_old_gwas job_eur=job-20221108-231351-762
# # noneur
# # german_affy6_old_gwas
# # QC  german_affy6_old_gwas job_noneur=job-20221108-231412-492
# job_eur=job-20221108-231351-762
# job_noneur=job-20221108-231412-492

# i=belgium_vermeire_gsa - downloaded
# # eur
# # belgium_vermeire_gsa
# # QC  belgium_vermeire_gsa job_eur=job-20221108-215757-683
# # noneur
# # belgium_vermeire_gsa
# # QC  belgium_vermeire_gsa job_noneur=job-20221108-215806-821
# job_eur=job-20221108-215757-683
# job_noneur=job-20221108-215806-821

# i=franchimont_gsa - downloaded
# # eur
# # franchimont_gsa
# # QC  franchimont_gsa job_eur=job-20221108-172154-599
# # noneur
# # franchimont_gsa
# # QC  franchimont_gsa job_noneur=job-20221108-172200-755
# job_eur=job-20221108-172154-599
# job_noneur=job-20221108-172200-755

# i=farkkila_gsa - downloaded
# # eur
# # farkkila_gsa
# # QC  farkkila_gsa job_eur=job-20221108-202557-269
# job_eur=job-20221108-202557-269

# i=franke_gsa - downloaded
# # eur
# # franke_gsa
# # QC  franke_gsa job_eur=job-20221108-172429-139
# # noneur
# # franke_gsa
# # QC  franke_gsa job_noneur=job-20221108-202210-077
# job_eur=job-20221108-172429-139
# job_noneur=job-20221108-202210-077

# i=niddk_rioux_gsa - downloaded
# # eur
# # niddk_rioux_gsa
# # QC  niddk_rioux_gsa job_eur=job-20221108-141156-826
# # noneur
# # niddk_rioux_gsa
# # QC  niddk_rioux_gsa job_noneur=job-20221108-152601-578
# job_eur=job-20221108-141156-826
# job_noneur=job-20221108-152601-578

# i=niddk_duerr_gsa - downloaded
# # eur
# # niddk_duerr_gsa
# # QC  niddk_duerr_gsa job_eur=job-20221108-121415-038
# # noneur
# # niddk_duerr_gsa
# # QC  niddk_duerr_gsa job_noneur=job-20221108-130034-912
# job_eur=job-20221108-121415-038
# job_noneur=job-20221108-130034-912

# i=niddk_cho_gsa - downloaded
# # eur
# # niddk_cho_gsa
# # QC  niddk_cho_gsa job_eur=job-20221108-121349-371
# # noneur
# # niddk_cho_gsa
# # QC  niddk_cho_gsa job_noneur=job-20221108-121354-511
# job_eur=job-20221108-121349-371
# job_noneur=job-20221108-121354-511

# i=moayyedi_imagine_gsa - downloaded
# # eur
# # moayyedi_imagine_gsa
# # QC  moayyedi_imagine_gsa job_eur=job-20221108-102036-436
# # noneur
# # moayyedi_imagine_gsa
# # QC  moayyedi_imagine_gsa job_noneur=job-20221108-105039-008
# job_eur=job-20221108-102036-436
# job_noneur=job-20221108-105039-008

# i=mccauley_new_gsa - downloaded
# # eur
# # mccauley_new_gsa
# # QC  mccauley_new_gsa job_eur=job-20221108-080942-825
# # noneur
# # mccauley_new_gsa
# # QC  mccauley_new_gsa job_noneur=job-20221108-100122-472
# job_eur=job-20221108-080942-825
# job_noneur=job-20221108-100122-472


# i=vermeire_gsa
# # eur
# # vermeire_gsa
# # QC  vermeire_gsa job_eur=job-20221108-022855-809
# # noneur
# # vermeire_gsa
# # QC  vermeire_gsa job_noneur=job-20221108-072138-725
# job_eur=job-20221108-022855-809
# job_noneur=job-20221108-072138-725


# i=netherlands_gsa - downloaded
# # eur
# # netherlands_gsa
# # QC  netherlands_gsa job_eur=job-20221107-220102-481
# # noneur
# # netherlands_gsa
# # QC  netherlands_gsa job_noneur=job-20221108-063521-637
# job_eur=job-20221107-220102-481
# job_noneur=job-20221108-063521-637

# i=niddk_silverberg_gsa - downloaded
# # eur
# # niddk_silverberg_gsa
# # QC  niddk_silverberg_gsa job_eur=job-20221107-195918-491
# # noneur
# # niddk_silverberg_gsa
# # QC  niddk_silverberg_gsa job_noneur=job-20221107-195923-797
# job_eur=job-20221107-195918-491
# job_noneur=job-20221107-195923-797

# i=newberry_share_gsa - downloaded
# # eur
# # newberry_share_gsa
# # QC  newberry_share_gsa job_eur=job-20221107-200223-729
# # noneur
# # newberry_share_gsa
# # QC  newberry_share_gsa job_noneur=job-20221107-212401-449
# job_eur=job-20221107-200223-729
# job_noneur=job-20221107-212401-449

# i=palotie_hus_gsa - downloaded
# # eur
# # palotie_hus_gsa
# # QC  palotie_hus_gsa job_eur=job-20221107-183613-675
# # noneur
# # palotie_hus_gsa
# # QC  palotie_hus_gsa job_noneur=job-20221107-183627-087
# job_eur=job-20221107-183613-675
# job_noneur=job-20221107-183627-087

# i=stampfer_gsa - downloaded
# # eur
# # stampfer_gsa
# # QC  stampfer_gsa job_eur=job-20221107-174329-567
# # noneur
# # stampfer_gsa
# # QC  stampfer_gsa job_noneur=job-20221107-180641-452
# job_eur=job-20221107-174329-567
# job_noneur=job-20221107-180641-452


# i=australia_omniexome - downloaded
# # eur
# # australia_omniexome
# # QC  australia_omniexome job_eur=job-20221107-162140-089
# # noneur
# # australia_omniexome
# # QC  australia_omniexome job_noneur=job-20221107-162552-719
# job_eur=job-20221107-162140-089
# job_noneur=job-20221107-162552-719

# i=prism_nfe_gsa - downloaded
# # eur
# # prism_nfe_gsa
# # QC  prism_nfe_gsa job_eur=job-20221107-162631-438
# # noneur
# # prism_nfe_gsa
# # QC  prism_nfe_gsa job_noneur=job-20221107-172106-173
# job_eur=job-20221107-162631-438
# job_noneur=job-20221107-172106-173

# i=slovenia_gsa - downloaded
# # eur
# # slovenia_gsa
# # QC  slovenia_gsa job_eur=job-20221107-143558-671
# job_eur=job-20221107-143558-671

# i=pekow_share_gsa - downloaded
# # eur
# # pekow_share_gsa
# # QC  pekow_share_gsa job_eur=job-20221107-140832-433
# # noneur
# # pekow_share_gsa
# # QC  pekow_share_gsa job_noneur=job-20221107-140835-390
# job_eur=job-20221107-140832-433
# job_noneur=job-20221107-140835-390


# i=mcgovern_gsa - downloaded
# # eur
# # mcgovern_gsa
# # QC  mcgovern_gsa job_eur=job-20221106-213914-184
# # noneur
# # mcgovern_gsa
# # QC  mcgovern_gsa job_noneur=job-20221107-080521-517
# job_eur=job-20221106-213914-184
# job_noneur=job-20221107-080521-517

# i=finland_illugwas - downloaded
# # eur
# # finland_illugwas
# # QC  finland_illugwas job_eur=job-20221107-103806-116
# job_eur=job-20221107-103806-116

# i=niddk_broad_gsa - downloaded
# # eur
# # niddk_broad_gsa
# # QC  niddk_broad_gsa job_eur=job-20221106-224840-628
# # noneur
# # niddk_broad_gsa
# # QC  niddk_broad_gsa job_noneur=job-20221107-080508-105
# job_eur=job-20221106-224840-628
# job_noneur=job-20221107-080508-105


# i=prism_nfe_gwas - downloaded
# # eur
# # prism_nfe_gwas
# # QC  prism_nfe_gwas job_eur=job-20221106-202749-956
# # noneur
# # prism_nfe_gwas
# # QC  prism_nfe_gwas job_noneur=job-20221106-202753-255
# job_eur=job-20221106-202749-956
# job_noneur=job-20221106-202753-255

# i=rioux_igenomed_gsa - downloaded
# # eur
# # rioux_igenomed_gsa
# # QC  rioux_igenomed_gsa job_eur=job-20221106-190947-897
# # noneur
# # rioux_igenomed_gsa
# # QC  rioux_igenomed_gsa job_noneur=job-20221106-190958-404
# job_eur=job-20221106-190947-897
# job_noneur=job-20221106-190958-404

# i=xavier_prism_gsa - downloaded
# # eur
# # xavier_prism_gsa
# # QC  xavier_prism_gsa job_eur=job-20221106-174415-867
# # noneur
# # xavier_prism_gsa
# # QC  xavier_prism_gsa job_noneur=job-20221106-175920-717
# job_eur=job-20221106-174415-867
# job_noneur=job-20221106-175920-717

# i=xavier_share_gsa - downloaded
# # eur
# # xavier_share_gsa
# # QC  xavier_share_gsa job_eur=job-20221106-160327-938
# # noneur
# # xavier_share_gsa
# # QC  xavier_share_gsa job_noneur=job-20221106-172925-698
# job_eur=job-20221106-160327-938
# job_noneur=job-20221106-172925-698

# i=lewis_sparc_gsa - downloaded
# # eur
# # lewis_sparc_gsa
# # QC  lewis_sparc_gsa job_eur=job-20221106-141026-232
# # noneur
# # lewis_sparc_gsa
# # QC  lewis_sparc_gsa job_noneur=job-20221106-142109-584
# job_eur=job-20221106-141026-232
# job_noneur=job-20221106-142109-584

# i=hyams_protect_gsa - downloaded
# # eur
# # hyams_protect_gsa
# # QC  hyams_protect_gsa job_eur=job-20221106-131830-519
# # noneur
# # hyams_protect_gsa
# # QC  hyams_protect_gsa job_noneur=job-20221106-133925-133
# job_eur=job-20221106-131830-519
# job_noneur=job-20221106-133925-133

# i=helmsley_xavier_prism_gsa - downloaded
# # eur
# # helmsley_xavier_prism_gsa
# # QC  helmsley_xavier_prism_gsa job_eur=job-20221106-105308-444
# # noneur
# # helmsley_xavier_prism_gsa
# # QC  helmsley_xavier_prism_gsa job_noneur=job-20221106-125510-399
# job_eur=job-20221106-105308-444
# job_noneur=job-20221106-125510-399

# i=niddk_feinstein_gsa - downloaded
# # eur
# # niddk_feinstein_gsa
# # QC  niddk_feinstein_gsa job_eur=job-20221105-224236-528
# # noneur
# # niddk_feinstein_gsa
# # QC  niddk_feinstein_gsa job_noneur=job-20221106-084549-050
# job_eur=job-20221105-224236-528
# job_noneur=job-20221106-084549-050

# i=niddk_old_gwas - completed
# # eur
# # niddk_old_gwas
# # QC  niddk_old_gwas job_eur=job-20221105-184744-458
# # noneur
# # niddk_old_gwas
# # QC  niddk_old_gwas job_noneur=job-20221106-081558-847
# job_eur=job-20221105-184744-458
# job_noneur=job-20221106-081558-847

# i=helmsley_prism_gsa - Completed
# # eur
# # helmsley_prism_gsa
# # QC  helmsley_prism_gsa job_eur=job-20221105-104614-542
# # noneur
# # helmsley_prism_gsa
# # QC  helmsley_prism_gsa job_noneur=job-20221105-170248-057
# job_eur=job-20221105-104614-542
# job_noneur=job-20221105-170248-057

# i=sands_msccr_gsa - Completed
# # eur
# # sands_msccr_gsa
# # QC  sands_msccr_gsa job_eur=job-20221105-093623-304
# # noneur
# # sands_msccr_gsa
# # QC  sands_msccr_gsa job_noneur=job-20221105-100727-785
# job_eur=job-20221105-093623-304
# job_noneur=job-20221105-100727-785

# i=bernstein_gsa - Completed
# # # eur
# # # bernstein_gsa
# # # QC  bernstein_gsa job_eur=job-20221105-092512-514
# # # noneur
# # # bernstein_gsa
# # # QC  bernstein_gsa job_noneur=job-20221105-092514-974
# job_eur=job-20221105-092512-514
# job_noneur=job-20221105-092514-974

# i=rioux_igenomed_gsa - Completed
# # eur
# # rioux_igenomed_gsa
# # QC  rioux_igenomed_gsa job_eur=job-20221105-084243-338
# # noneur
# # rioux_igenomed_gsa
# # QC  rioux_igenomed_gsa job_noneur=job-20221105-084245-200
# job_eur=job-20221105-084243-338
# job_noneur=job-20221105-084245-200

# i=farkkila_gsa - Completed
# # eur
# # farkkila_gsa
# # QC  farkkila_gsa job_eur=job-20221105-084154-330
# job_eur=job-20221105-084154-330

############# 

for j in ${ancestry[@]}
do rm /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/${j}/*.vcf.gz
done

for j in ${ancestry[@]}
do rm /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/${j}/chr*.info.gz
done

for j in ${ancestry[@]}
do rm /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/${j}/chr_*.zip
done

for j in ${ancestry[@]}
do rm /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/${j}/results.md5
done

rm /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/logs/nohup.out
rm /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/logs/nohup_ch23.out
rm /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/logs/donwload_imputed_results_chr23_*
rm /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/logs/status_*
  
{
  echo ${i}
  
  j=eur
  TOKEN='YOUR_API_TOKEN';
  curl -H "X-Auth-Token: $TOKEN" https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/${job_eur} \
  > ${path_gwas}imputed/${i}/2022/logs/status_${i}_${j}_qc.json
  
  j=noneur
  TOKEN='YOUR_API_TOKEN';
  curl -H "X-Auth-Token: $TOKEN" https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/${job_noneur} \
  > ${path_gwas}imputed/${i}/2022/logs/status_${i}_${j}_qc.json
  
  /software/R-4.3.1/bin/Rscript \
  /path/to/user/scripts/IIBDGC/create_imputation_download_commnads.R ${i} \
  > ${path_gwas}imputed/${i}/2022/logs/create_imputation_download_commnads_${i}.R_out
  
  ### edit files:
  j=eur
  sed -i '2 i\path=\/lustre\/scratch123\/hgi\/projects\/ibdgwas\/IIBDGC\/imputed\/${i}\/2022\/' ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
  sed -i '2 i\j='${j} ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
  sed -i '2 i\i='${i} ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
  
  j=noneur
  sed -i -e '1,4d' ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
  sed -i '2 i\j='${j} ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
  
  cat  ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_eur  ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_noneur \
  > ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}
  rm ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_*
    
  # # if noneur only:
  # j=noneur
  # TOKEN='YOUR_API_TOKEN';
  # curl -H "X-Auth-Token: $TOKEN" https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/${job_noneur} \
  # > ${path_gwas}imputed/${i}/2022/logs/status_${i}_${j}_qc.json
  # 
  # /software/R-4.3.1/bin/Rscript \
  # /path/to/user/scripts/IIBDGC/create_imputation_download_commnads.R ${i} \
  # > ${path_gwas}imputed/${i}/2022/logs/create_imputation_download_commnads_${i}.R_out
  # 
  # j=noneur
  # sed -i '2 i\path=\/lustre\/scratch123\/hgi\/projects\/ibdgwas\/IIBDGC\/imputed\/${i}\/2022\/' ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
  # sed -i '2 i\j='${j} ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
  # sed -i '2 i\i='${i} ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
  # 
  # cat ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_noneur \
  # > ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}
  # rm ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_*

  sed -i -e 's/-o .*.zip --create-dirs/-P ${path}\/${j}/g' ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}
  sed -i -e 's/-o .*.md5 --create-dirs/-P ${path}\/${j}/g' ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}
  sed -i -e 's/curl -L/wget --tries=75 -c/g' ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}
  # more ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}
  
  cd /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/logs/
  nohup bash /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/logs/donwload_imputed_results_${i} > nohup_new.out &
    
  unset job_eur
  unset job_noneur

  
}


##############
##############


studies=(australia_omniexome niddk_feinstein_gsa)
studies=(belgium_franchimont_gsa lithuania_gsa gwas1 pittsburgh_gsa italy_gsa)
studies=(belgium_vermeire_gsa belgium_inf1_old_gwas cedars_gsa all_hce belgium_inf2_old_gwas slovenia_gsa finland_illugwas ccfa_gsa weersma_gsa)
studies=(kiel_austria_sibdcs_gsa vermeire_gsa niddk_old_gwas cedars_370k_old_gwas)
studies=(cedars_610k_old_gwas gwas2 finland_illugwas cedars_omni_old_gwas basque_gsa netherlands_gsa prism_nfe_gsa belgium_louis_gsa)
studies=(prism_nfe_gwas niddk_broad_gsa norway_affy6_old_gwas)

studies=(swedish_uc_old_gwas german_affy6_old_gwas mccauley_gsa)


for i in ${studies[@]}
do echo ${i} && for j in ${ancestry[@]}
do echo ${j} && ls -la ${path_gwas}imputed/${i}/2022/${j}/chr_9.zip
done
done

##############
##############

####################################################
# 23.10  Unzip files:

path_gwas=/path/to/ibdgwas/IIBDGC/
  
ancestry=(eur noneur)
MEM=200

### CHR1 TO CHR22
for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
for chr in {1..22}
do
bsub -J"pl2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_6_post_server_qc_${j}_${i}_chr${chr} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_6_post_server_qc_${j}_${i}_chr${chr} \
"unzip -P PkSw6GlBg5VGuc ${path_gwas}imputed/${i}/2022/${j}/chr_${chr}.zip -d ${path_gwas}imputed/${i}/2022/${j}/"
done 
done
done

### CHR23
for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
bsub -J"pl2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_6_post_server_qc_${j}_${i}_chr23 \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_6_post_server_qc_${j}_${i}_chr23 \
"unzip -P PkSw6GlBg5VGuc ${path_gwas}imputed/${i}/2022/${j}/chr_X.zip -d ${path_gwas}imputed/${i}/2022/${j}/"
done
done


chr=1
for i in ${studies[@]}
do echo ${i} && for j in ${ancestry[@]}
do echo ${j} && ls -la ${path_gwas}imputed/${i}/2022/${j}/chr${chr}.info.gz
done
done

###### remove the .zip file for the studies unzipped - double check first unzipped files exist

for i in ${studies[@]}
do echo ${i} && for j in ${ancestry[@]}
do echo ${j} && ls -la ${path_gwas}imputed/${i}/2022/${j}/chr*.vcf.gz | wc -l && \
ls -la ${path_gwas}imputed/${i}/2022/${j}/chr*.info.gz | wc -l && \
ls -la ${path_gwas}imputed/${i}/2022/${j}/chr9.info.gz && \
ls -la ${path_gwas}imputed/${i}/2022/${j}/chrX.dose.vcf.gz && echo "        "
done
done

for i in ${studies[@]}
do echo ${i} && for j in ${ancestry[@]}
do echo ${j} && ls -la ${path_gwas}imputed/${i}/2022/${j}/chr_X.zip && ls -la ${path_gwas}imputed/${i}/2022/${j}/chr_9.zip
done
done

for i in ${studies[@]}
do echo ${i} && for j in ${ancestry[@]}
do rm ${path_gwas}imputed/${i}/2022/${j}/chr_*.zip
done
done



####################################################
# 23.11 create imputation summary files and plots:

for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
mkdir -p ${path_gwas}imputed/${i}/2022/${j}/plots/
done
done

path_gwas="/path/to/ibdgwas/IIBDGC/"

# CREATE final rainfall plots:
for i in ${studies[@]}
do for chr in {1..23}
do bash ${path_gwas}scripts/post_imputation_summary_1.sh ${i} ${chr}
done
done


## double check results
for i in ${studies[@]}
do echo ${i} && for chr in {1..23}
do echo ${chr} && tail -50 ${path_gwas}scripts/logs/post_imputation_summary_stdout_${i}_${chr} | grep -E "completed"
done
done

# CREATE final summary tables of imputation accuracy:
for i in ${studies[@]}
do bash ${path_gwas}scripts/post_imputation_summary_3.sh ${i}
done

### TO CHECK RUNNING FOR 


for i in ${studies[@]}
do echo ${i} && tail -50 /path/to/ibdgwas/IIBDGC/scripts/logs/post_imputation_summary_stdout_${i}_allchr | grep -E "completed"
done

# ### visually inspect plots:
# for i in ${studies[@]}
# done


### output number of samples in post-imputation files:

MEM=500
for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
bsub -J"pl2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_7_post_server_imputation_${j}_${i}_chr23 \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_7_post_server_imputation_${j}_${i}_chr23 \
"/path/to/software/bcftools-1.16/./bcftools \
query -l ${path_gwas}imputed/${i}/2022/${j}/chrX.dose.vcf.gz > ${path_gwas}imputed/test/${i}_${j}_2020_samples_imputed_files_chr23"
done
done

MEM=500
for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
bsub -J"pl2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_7_post_server_imputation_${j}_${i}_chr9 \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_7_post_server_imputation_${j}_${i}_chr9 \
"/path/to/software/bcftools-1.16/./bcftools \
query -l ${path_gwas}imputed/${i}/2022/${j}/chr9.dose.vcf.gz > ${path_gwas}imputed/test/${i}_${j}_2020_samples_imputed_files_chr9"
done
done


for i in ${studies[@]}
do echo ${i} && for j in ${ancestry[@]}
do tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_7_post_server_imputation_${j}_${i}_chr23 | grep -E "completed" && \
tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_7_post_server_imputation_${j}_${i}_chr9 | grep -E "completed" 
done
done

###### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","swedish_uc_old_gwas"
           ,"mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

cohorts<-c("bernstein_gsa","farkkila_gsa","franchimont_gsa",
"franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
"mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
"niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
"rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
"xavier_share_gsa")


length(cohorts)
# [1] 58

ancestry<-c("eur","noneur")

for (i in 1:length(cohorts)) {
  
  print(cohorts[i])
  
  for (j in 1:length(ancestry)) {
   
    print(ancestry[j])
    
    file_tmp1<-paste(path,"imputed/test/",cohorts[i],"_",ancestry[j],"_2020_samples_imputed_files_chr9",sep="")
    file_tmp2<-paste(path,"imputed/test/",cohorts[i],"_",ancestry[j],"_2020_samples_imputed_files_chr23",sep="")
    
    if(file.exists(file_tmp1)) {
      tmp<-fread(file_tmp1,head=F)
      print(nrow(tmp))
      tmp<-fread(file_tmp2,head=F)
      print(nrow(tmp))
    } else{
      print("File not ready")
    }
    rm(tmp)
  }
}

################
# 
# [1] "australia_omniexome"
# [1] "eur"
# [1] 1245
# [1] 1245
# [1] "noneur"
# [1] 54
# [1] 54

# [1] "gwas1"
# [1] "eur"
# [1] 4652
# [1] 4652
# [1] "noneur"
# [1] 7
# [1] 7
# [1] "gwas2"
# [1] "eur"
# [1] 7758
# [1] 7758
# [1] "noneur"
# [1] 15
# [1] 15
# [1] "all_hce"
# [1] "eur"
# [1] 22588
# [1] 22588
# [1] "noneur"
# [1] 1175
# [1] 1175

# [1] "pittsburgh_gsa"
# [1] "eur"
# [1] 2701
# [1] 2701
# [1] "noneur"
# [1] 8
# [1] 8
# [1] "spain_gsa"
# [1] "eur"
# [1] 3396
# [1] 3396
# [1] "noneur"
# [1] 12
# [1] 12
# [1] "italy_gsa"
# [1] "eur"
# [1] 938
# [1] 938
# [1] "noneur"
# [1] 5
# [1] 5
# [1] "kiel_austria_sibdcs_gsa"
# [1] "eur"
# [1] 14086
# [1] 14086
# [1] "noneur"
# [1] 243
# [1] 243


# [1] "netherlands_gsa"
# [1] "eur"
# [1] 4315
# [1] 4315
# [1] "noneur"
# [1] 222
# [1] 222
# [1] "slovenia_gsa"
# [1] "eur"
# [1] 263
# [1] 263
# [1] "noneur"
# [1] 0
# [1] 0
# [1] "sweden_gsa"
# [1] "eur"
# [1] 1364
# [1] 1364
# [1] "noneur"
# [1] 32
# [1] 32

# [1] "niddk_broad_gsa"
# [1] "eur"
# [1] 5130
# [1] 5130
# [1] "noneur"
# [1] 296
# [1] 296

# [1] "niddk_feinstein_gsa"
# [1] "eur"
# [1] 6965
# [1] 6965
# [1] "noneur"
# [1] 1159
# [1] 1159
# [1] "basque_gsa"
# [1] "eur"
# [1] 1487
# [1] 1487
# [1] "noneur"
# [1] 11
# [1] 11
# [1] "prism_nfe_gsa"
# [1] "eur"
# [1] 426
# [1] 426
# [1] "noneur"
# [1] 36
# [1] 36
# [1] "lithuania_gsa"
# [1] "eur"
# [1] 2213
# [1] 2213
# [1] "noneur"
# [1] 12
# [1] 12
# [1] "belgium_louis_gsa"
# [1] "eur"
# [1] 1499
# [1] 1499
# [1] "noneur"
# [1] 12
# [1] 12
# [1] "belgium_franchimont_gsa"
# [1] "eur"
# [1] 1459
# [1] 1459
# [1] "noneur"
# [1] 32
# [1] 32
# [1] "belgium_vermeire_gsa"
# [1] "eur"
# [1] 3894
# [1] 3894
# [1] "noneur"
# [1] 80
# [1] 80

# [1] "prism_nfe_gwas"
# [1] "eur"
# [1] 717
# [1] 717
# [1] "noneur"
# [1] 93
# [1] 93

# [1] "finland_illugwas"
# [1] "eur"
# [1] 440
# [1] 440
# [1] "noneur"
# [1] 4
# [1] 4

# [1] "german_affy6_old_gwas"
# [1] "eur"
# [1] 2787
# [1] 2787
# [1] "noneur"
# [1] 10
# [1] 10


# [1] "norway_affy6_old_gwas"
# [1] "eur"
# [1] 542
# [1] 542
# [1] "noneur"
# [1] 2
# [1] 2

# [1] "belgium_inf1_old_gwas"
# [1] "eur"
# [1] 1396
# [1] 1396
# [1] "noneur"
# [1] 10
# [1] 10
# [1] "belgium_inf2_old_gwas"
# [1] "eur"
# [1] 270
# [1] 270
# [1] "noneur"
# [1] 1
# [1] 1

# [1] "niddk_old_gwas"
# [1] "eur"
# [1] 2712
# [1] 2712
# [1] "noneur"
# [1] 20
# [1] 20
# [1] "cedars_370k_old_gwas"
# [1] "eur"
# [1] 509
# [1] 509
# [1] "noneur"
# [1] 88
# [1] 88
# [1] "cedars_610k_old_gwas"
# [1] "eur"
# [1] 825
# [1] 825
# [1] "noneur"
# [1] 56
# [1] 56
# [1] "cedars_omni_old_gwas"
# [1] "eur"
# [1] 1104
# [1] 1104
# [1] "noneur"
# [1] 103
# [1] 103

# [1] "swedish_uc_old_gwas"
# [1] "eur"
# [1] 1255
# [1] 0
# [1] "noneur"
# [1] 3
# [1] 0
# [1] "mccauley_gsa"
# [1] "eur"
# [1] 758
# [1] 758
# [1] "noneur"
# [1] 18
# [1] 18

# [1] "ccfa_gsa"
# [1] "eur"
# [1] 1932
# [1] 1932
# [1] "noneur"
# [1] 235
# [1] 235
# [1] "cedars_gsa"
# [1] "eur"
# [1] 2699
# [1] 2699
# [1] "noneur"
# [1] 360
# [1] 360
# [1] "bernstein_gsa"
# [1] "eur"
# [1] 482
# [1] 482
# [1] "noneur"
# [1] 27
# [1] 27
# [1] "farkkila_gsa"
# [1] "eur"
# [1] 68
# [1] 68
# [1] "noneur"
# [1] 0
# [1] 0
# [1] "franchimont_gsa"
# [1] "eur"
# [1] 2360
# [1] 2360
# [1] "noneur"
# [1] 409
# [1] 409
# [1] "franke_gsa"
# [1] "eur"
# [1] 853
# [1] 853
# [1] "noneur"
# [1] 11
# [1] 11
# [1] "helmsley_prism_gsa"
# [1] "eur"
# [1] 669
# [1] 669
# [1] "noneur"
# [1] 84
# [1] 84
# [1] "helmsley_xavier_prism_gsa"
# [1] "eur"
# [1] 1171
# [1] 1171
# [1] "noneur"
# [1] 110
# [1] 110
# [1] "hyams_protect_gsa"
# [1] "eur"
# [1] 323
# [1] 323
# [1] "noneur"
# [1] 93
# [1] 93
# [1] "lewis_sparc_gsa"
# [1] "eur"
# [1] 2529
# [1] 2529
# [1] "noneur"
# [1] 317
# [1] 317
# [1] "mccauley_new_gsa"
# [1] "eur"
# [1] 1323
# [1] 1323
# [1] "noneur"
# [1] 287
# [1] 287
# [1] "mcgovern_gsa"
# [1] "eur"
# [1] 5149
# [1] 5149
# [1] "noneur"
# [1] 829
# [1] 829
# [1] "moayyedi_imagine_gsa"
# [1] "eur"
# [1] 1011
# [1] 1011
# [1] "noneur"
# [1] 109
# [1] 109
# [1] "newberry_share_gsa"
# [1] "eur"
# [1] 747
# [1] 747
# [1] "noneur"
# [1] 111
# [1] 111
# [1] "niddk_cho_gsa"
# [1] "eur"
# [1] 1611
# [1] 1611
# [1] "noneur"
# [1] 143
# [1] 143
# [1] "niddk_duerr_gsa"
# [1] "eur"
# [1] 1868
# [1] 1868
# [1] "noneur"
# [1] 69
# [1] 69
# [1] "niddk_rioux_gsa"
# [1] "eur"
# [1] 871
# [1] 871
# [1] "noneur"
# [1] 42
# [1] 42
# [1] "niddk_silverberg_gsa"
# [1] "eur"
# [1] 2112
# [1] 2112
# [1] "noneur"
# [1] 235
# [1] 235
# [1] "palotie_hus_gsa"
# [1] "eur"
# [1] 853
# [1] 853
# [1] "noneur"
# [1] 16
# [1] 16
# [1] "pekow_share_gsa"
# [1] "eur"
# [1] 538
# [1] 538
# [1] "noneur"
# [1] 94
# [1] 94
# [1] "rioux_igenomed_gsa"
# [1] "eur"
# [1] 171
# [1] 171
# [1] "noneur"
# [1] 10
# [1] 10
# [1] "sands_msccr_gsa"
# [1] "eur"
# [1] 1120
# [1] 1120
# [1] "noneur"
# [1] 299
# [1] 299
# [1] "stampfer_gsa"
# [1] "eur"
# [1] 1439
# [1] 1439
# [1] "noneur"
# [1] 25
# [1] 25

# [1] "vermeire_gsa"
# [1] "eur"
# [1] 4540
# [1] 4540
# [1] "noneur"
# [1] 119
# [1] 119
# [1] "weersma_gsa"
# [1] "eur"
# [1] 682
# [1] 682
# [1] "noneur"
# [1] 23
# [1] 23

# [1] "xavier_prism_gsa"
# [1] "eur"
# [1] 619
# [1] 619
# [1] "noneur"
# [1] 67
# [1] 67
# [1] "xavier_share_gsa"
# [1] "eur"
# [1] 637
# [1] 637
# [1] "noneur"
# [1] 55
# [1] 55



################

####################################################
# 23.12 Final double check to EmpRsq

# still some genotyped variants with EmpRsq<0.5, double check distribution of these per study:

# final rainfall plots:
for i in ${studies[@]}
do bash ${path_gwas}scripts/post_imputation_summary_4.sh ${i}
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}scripts/logs/post_imputation_summary_4_stdout_${i}_allchr  | grep -E "completed"
done

###### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("australia_omniexome","gwas1"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","swedish_uc_old_gwas"
           ,"mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","xavier_prism_gsa",
           "xavier_share_gsa","gwas2","all_hce","kiel_austria_sibdcs_gsa","vermeire_gsa","weersma_gsa")

new_cohorts<-c("bernstein_gsa","farkkila_gsa","franchimont_gsa",
               "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
               "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
               "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
               "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","xavier_prism_gsa",
               "xavier_share_gsa")


length(cohorts)

# [1] 58

ancestry<-"eur"
j=1

dat<-as.data.frame(matrix(ncol=4,nrow=length(cohorts)))
colnames(dat)<-c("cohort","n_genotyped_variants","n_variants_empRsq_lower_0.5","n_variants_empRsq_lower_0.4")
dat$cohort<-cohorts

for (i in 1:length(cohorts)) {
  
  print(cohorts[i])
  
  file_tmp<-paste(path,"imputed/",cohorts[i],"/2022/eur/plots/",cohorts[i],"_eur_low_EmpRsq_after_first_roudn_exclussions.tsv",sep="")
  
  if(file.exists(file_tmp)) {
    
    tmp<-read.table(file_tmp,head=T)
    tmp$study<-cohorts[i]
    
    if(!exists("all_tmp")) {
      all_tmp<-tmp
    }else{
      all_tmp<-rbind(all_tmp,tmp)
    }
    tmp_nvar<-read.table(paste(path,"imputation_ready/test/count_n_variants_",cohorts[i],"_",ancestry[j],sep=""),head=F)
    
    dat$n_genotyped_variants[which(dat$cohort==cohorts[i])]<-sum(tmp_nvar$V1)
    dat$n_variants_empRsq_lower_0.5[which(dat$cohort==cohorts[i])]<-nrow(tmp)
    dat$n_variants_empRsq_lower_0.4[which(dat$cohort==cohorts[i])]<-nrow(tmp[which(tmp$EmpRsq<0.4),])
    rm(tmp,tmp_nvar)
  
  }
}

dat$percentage_empRsq_lower_0.5<-dat$n_variants_empRsq_lower_0.5/dat$n_genotyped_variants
dat<-dat[order(dat$percentage_empRsq_lower_0.5,decreasing=T),]

dat[which( (dat$n_variants_empRsq_lower_0.5/dat$n_genotyped_variants)>=0.001),]
#                     cohort n_genotyped_variants n_variants_empRsq_lower_0.5
# 5                italy_gsa               475735                        2566
# 55                 all_hce               252293                         752
# 19   german_affy6_old_gwas               729844                        1973
# 7             slovenia_gsa               462215                         890
# 2                    gwas1               426348                         707
# 20   norway_affy6_old_gwas               695725                         944
# 14       belgium_louis_gsa               474265                         613
# 15 belgium_franchimont_gsa               472743                         495
# 33         franchimont_gsa               474922                         477
#    n_variants_empRsq_lower_0.4 percentage_empRsq_lower_0.5
# 5                          409                 0.005393759
# 55                         314                 0.002980661
# 19                        1403                 0.002703317
# 7                          197                 0.001925511
# 2                          354                 0.001658270
# 20                         838                 0.001356858
# 14                         200                 0.001292526
# 15                         119                 0.001047081
# 33                         133                 0.001004375




cohorts_2<-c("italy_gsa","all_hce","german_affy6_old_gwas","slovenia_gsa","gwas1","norway_affy6_old_gwas","belgium_louis_gsa","belgium_franchimont_gsa","franchimont_gsa")

pdf(paste(path,"imputed/plots/low_EmpRsq_after_first_roudn_exclussions_2022.tsv.pdf",sep=""),width = 60, height = 7,onefile=FALSE)
ggplot(all_tmp[which(all_tmp$MAF>0.0001 & !(all_tmp$study %in% cohorts_2)),], aes(x=study, y=EmpRsq)) + geom_violin(trim=T) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25)
dev.off()

p1<-ggplot(all_tmp, aes(x=EmpRsq,color=study)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_density() + ylim(0,100)
p2<-ggplot(all_tmp[which(all_tmp$study %in% cohorts_2),], aes(x=EmpRsq,color=study)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_density() + ylim(0,100)

pdf(paste(path,"imputed/plots/low_EmpRsq_after_first_roudn_exclussions_2022_density.pdf",sep=""),width = 20, height = 14,onefile=FALSE)
ggarrange(p1,p2,nrow=2,legend="right")
dev.off()

write.table(dat,paste(path,"imputed/plots/low_EmpRsq_after_first_roudn_exclussions_2022.tsv",sep=""),col.names=T,row.names=F,quote=F,sep="\t")



all_tmp<-all_tmp[which(all_tmp$study %in% cohorts_2),]

for (j in 1:length(cohorts_2)) {
  
  print(cohorts_2[j])
  
  tmp1<-read.table(paste(path,"imputed/",cohorts_2[j],"/qc/2022/list_variants_to_exclude_2",sep=""))
  tmp2<-all_tmp[which(all_tmp$study==cohorts_2[j]),"SNP",drop=F]


  tmp2$SNP<-gsub(":","_",tmp2$SNP)
  tmp2$SNP<-sub("_",":",tmp2$SNP)
  tmp2$SNP<-sub("chrX","chr23",tmp2$SNP)
  tmp2$SNP<-sub("chr","",tmp2$SNP)
  
  colnames(tmp2)<-"V1"
  
  eur<-read.table(paste(path,"pre_imputation/QC/",cohorts_2[j],"/",cohorts_2[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_eur_nomonom.bim",sep=""),head=F)
  print(paste("N variants exclude in EUR: ",nrow(eur[which(eur$V2 %in% tmp2$V1),]),sep=""))
  
  file_noneur<-paste(path,"pre_imputation/QC/",cohorts_2[j],"/",cohorts_2[j],"_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_noneur_nomonom.bim",sep="")
  if(file.exists(file_noneur)) {
    noneur<-read.table(file_noneur,head=F)
    print(paste("N variants exclude in nonEUR: ",nrow(noneur[which(noneur$V2 %in% tmp2$V1),]),sep=""))
    
  }
  
  tmp<-rbind(tmp1,tmp2)
  
  write.table(tmp,paste(path,"imputed/",cohorts_2[j],"/qc/2022/list_variants_to_exclude_3",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
  
}

q("no")

####################
 
# [1] "italy_gsa"
# [1] "N variants exclude in EUR: 2566"
# [1] "N variants exclude in nonEUR: 735"
# [1] "all_hce"
# [1] "N variants exclude in EUR: 752"
# [1] "N variants exclude in nonEUR: 748"
# [1] "german_affy6_old_gwas"
# [1] "N variants exclude in EUR: 1973"
# [1] "N variants exclude in nonEUR: 938"
# [1] "slovenia_gsa"
# [1] "N variants exclude in EUR: 890"
# [1] "gwas1"
# [1] "N variants exclude in EUR: 707"
# [1] "N variants exclude in nonEUR: 212"
# [1] "norway_affy6_old_gwas"
# [1] "N variants exclude in EUR: 944"
# [1] "N variants exclude in nonEUR: 175"
# [1] "belgium_louis_gsa"
# [1] "N variants exclude in EUR: 613"
# [1] "N variants exclude in nonEUR: 234"
# [1] "belgium_franchimont_gsa"
# [1] "N variants exclude in EUR: 495"
# [1] "N variants exclude in nonEUR: 304"
# [1] "franchimont_gsa"
# [1] "N variants exclude in EUR: 477"
# [1] "N variants exclude in nonEUR: 414"

####################


##########################################################################
# 23.13 Convert ped/map files to VCF files - subset of studies to re-run #
##########################################################################


studies=(franchimont_gsa)
path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=1000

ancestry=(eur noneur)

for chr in {1..23}
do
for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_4_post_server_qc_${j}_${i}_chr${chr} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_4_post_server_qc_${j}_${i}_chr${chr} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_${j}_nomonom \
--allow-no-sex \
--keep-allele-order \
--flip ${path_gwas}imputed/${i}/qc/2022/list_variants_negative_EmpR_toflip \
--exclude ${path_gwas}imputed/${i}/qc/2022/list_variants_to_exclude_3 \
--chr ${chr} \
--output-chr chr26 --recode vcf-iid --out ${path_gwas}imputation_ready/${i}/${i}_${j}_chr${chr}_2022"
done
done
done

# to the same with old studies pointing to old flip file
studies=(german_affy6_old_gwas italy_gsa all_hce slovenia_gsa gwas1 norway_affy6_old_gwas belgium_louis_gsa belgium_franchimont_gsa)

for chr in {1..23}
do
for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_4_post_server_qc_${j}_${i}_chr${chr} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_4_post_server_qc_${j}_${i}_chr${chr} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_${j}_nomonom \
--allow-no-sex \
--keep-allele-order \
--flip ${path_gwas}imputed/${i}/qc/list_variants_negative_EmpR_toflip \
--exclude ${path_gwas}imputed/${i}/qc/2022/list_variants_to_exclude_3 \
--chr ${chr} \
--output-chr chr26 --recode vcf-iid --out ${path_gwas}imputation_ready/${i}/${i}_${j}_chr${chr}_2022"
done
done
done

studies=(german_affy6_old_gwas italy_gsa all_hce slovenia_gsa gwas1 norway_affy6_old_gwas belgium_louis_gsa belgium_franchimont_gsa franchimont_gsa)

for i in ${studies[@]}
do
echo ${i} && for j in ${ancestry[@]}
do
echo ${j} && tail -40 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_4_post_server_qc_${j}_${i}_chr${chr} | grep -E "completed"
done
done

for i in ${studies[@]}
do
for j in ${ancestry[@]}
do ls -la ${path_gwas}imputation_ready/${i}/${i}_${j}_chr${chr}_2022.vcf
done
done

########################################################
# 23.14 Create a sorted vcf.gz - final round imputation #
########################################################

studies=(german_affy6_old_gwas italy_gsa all_hce slovenia_gsa gwas1 norway_affy6_old_gwas belgium_louis_gsa belgium_franchimont_gsa franchimont_gsa)

ancestry=(eur noneur)

for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
for chr in {1..23}
do
bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_5_post_server_qc_${j}_${i}_chr${chr} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_5_post_server_qc_${j}_${i}_chr${chr} \
"/path/to/software/bcftools-1.16/./bcftools sort ${path_gwas}imputation_ready/${i}/${i}_${j}_chr${chr}_2022.vcf \
-Oz -o ${path_gwas}imputation_ready/${i}/${i}_${j}_chr${chr}_2022.vcf.gz"
done
done
done


for i in ${studies[@]}
do for j in ${ancestry[@]}
do echo ${i} && ls -la ${path_gwas}imputation_ready/${i}/${i}_${j}_chr*_2022.vcf.gz |  wc -l
done
done


for i in ${studies[@]}
do for j in ${ancestry[@]}
do echo ${i} && ls -la ${path_gwas}imputation_ready/${i}/${i}_${j}_chr1_2022.vcf.gz
done
done

# rm intermediate files:
for i in ${studies[@]}
do rm ${path_gwas}imputation_ready/${i}/${i}_*_chr*.vcf
done
for i in ${studies[@]}
do rm ${path_gwas}imputation_ready/${i}/*.hh
done
for i in ${studies[@]}
do rm ${path_gwas}imputation_ready/${i}/*.log
done
for i in ${studies[@]}
do rm ${path_gwas}imputation_ready/${i}/*.nosex
done


## get total number of variants per final file for all studies
studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)
for i in ${studies[@]}
do
rm ${path_gwas}imputation_ready/test/count_n_variants_${i}_*
done

for i in ${studies[@]}
do 
for j in ${ancestry[@]}
do
for chr in {1..23}
do
/path/to/software/bcftools-1.16/./bcftools stats ${path_gwas}imputation_ready/${i}/${i}_${j}_chr${chr}_2022.vcf.gz | grep "number of records" | \
cut -f4 | awk 'NR>1' >> ${path_gwas}imputation_ready/test/count_n_variants_${i}_${j}
done
done
done

###### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("german_affy6_old_gwas","italy_gsa","all_hce","slovenia_gsa","gwas1","norway_affy6_old_gwas","belgium_louis_gsa","belgium_franchimont_gsa","franchimont_gsa")

length(cohorts)
# [1] 9

ancestry<-c("eur","noneur")

for (i in 1:length(cohorts)) {
  
  print(cohorts[i])
  
  for (j in 1:length(ancestry)) {
    print(ancestry[j])
    a<-tryCatch(read.table(paste(path,"imputation_ready/test/count_n_variants_",cohorts[i],"_",ancestry[j],sep=""),head=F), error=function(e) NULL)
    if(!is.null(a)) {
      print(sum(a$V1))
    } else {
      print("No variants")
    }
    
  }
  print("           ")
  rm(a)
}

#############################

# [1] "german_affy6_old_gwas"
# [1] "eur"
# [1] 727871
# [1] "noneur"
# [1] 668034
# [1] "           "
# [1] "italy_gsa"
# [1] "eur"
# [1] 475735
# [1] "noneur"
# [1] 285842
# [1] "           "
# [1] "all_hce"
# [1] "eur"
# [1] 252293
# [1] "noneur"
# [1] 256495
# [1] "           "
# [1] "slovenia_gsa"
# [1] "eur"
# [1] 462215
# [1] "noneur"
# [1] "No variants"
# [1] "           "
# [1] "gwas1"
# [1] "eur"
# [1] 426348
# [1] "noneur"
# [1] 359066
# [1] "           "
# [1] "norway_affy6_old_gwas"
# [1] "eur"
# [1] 695725
# [1] "noneur"
# [1] 447032
# [1] "           "
# [1] "belgium_louis_gsa"
# [1] "eur"
# [1] 474265
# [1] "noneur"
# [1] 341845
# [1] "           "
# [1] "belgium_franchimont_gsa"
# [1] "eur"
# [1] 472743
# [1] "noneur"
# [1] 412958
# [1] "           "
# [1] "franchimont_gsa"
# [1] "eur"
# [1] 474922
# [1] "noneur"
# [1] 472615
# [1] "           "

#############################

##################################################################################################################################################################

####################################################################################################
# 23.15 submit third round of imputation  for the subset of studies with large % variants empR>0.5 #
####################################################################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
studies=(gwas1)
ancestry=(eur)

i=gwas1
for j in ${ancestry[@]}
do echo ${i} && ls -la ${path_gwas}imputation_ready/${i}/${i}_${j}_chr*_2022.vcf.gz |  wc -l
done
for j in ${ancestry[@]}
do echo ${j} && python2 /path/to/user/scripts/IIBDGC/topmed_imputation_${j}.py ${i}
done

i=gwas1
job_eur=job-20221219-083616-516


#############
# downloading

# i=belgium_franchimont_gsa
# # eur
# # belgium_franchimont_gsa
# # QC  belgium_franchimont_gsa job_eur=job-20221128-093233-984
# # noneur
# # belgium_franchimont_gsa
# # QC  belgium_franchimont_gsa job_noneur=job-20221128-093257-802
# job_eur=job-20221128-093233-984
# job_noneur=job-20221128-093257-802


# i=gwas1
# # # eur
# # # gwas1
# # # QC  gwas1 job_eur=job-20221124-230856-642 - downloaded
# # # noneur
# # # gwas1
# # # QC  gwas1 job_noneur=job-20221124-213347-378 - downloaded
# job_eur=job-20221124-230856-642
# job_noneur=job-20221124-213347-378

# i=franchimont_gsa
# # eur
# # franchimont_gsa
# # QC  franchimont_gsa job_eur=job-20221124-104054-779 - downloaded
# # noneur
# # franchimont_gsa
# # QC  franchimont_gsa job_noneur=job-20221124-081041-223 - downloaded
# job_noneur=job-20221124-081041-223
# job_eur=job-20221124-104054-779

# i=belgium_louis_gsa
# # eur
# # belgium_louis_gsa
# # QC  belgium_louis_gsa job_eur=job-20221125-073204-448 - downloaded
# # noneur
# # belgium_louis_gsa
# # QC  belgium_louis_gsa job_noneur=job-20221125-063124-013 - downloaded
# job_eur=job-20221125-073204-448
# job_noneur=job-20221125-063124-013

# i=slovenia_gsa
# # eur
# # slovenia_gsa
# # QC  slovenia_gsa job_eur=job-20221124-112409-027 - downloaded
# job_eur=job-20221124-112409-027

# i=german_affy6_old_gwas
# # eur
# # german_affy6_old_gwas
# # QC  german_affy6_old_gwas job_eur=job-20221123-193927-526 - downloaded
# # german_affy6_old_gwas
# # QC  german_affy6_old_gwas job_noneur=job-20221123-214101-568 - downloaded
# job_eur=job-20221123-193927-526
# job_noneur=job-20221123-214101-568


# i=norway_affy6_old_gwas
# # eur
# # norway_affy6_old_gwas
# # QC  norway_affy6_old_gwas job_eur=job-20221124-035203-109 - downloaded
# # noneur
# # norway_affy6_old_gwas
# # QC  norway_affy6_old_gwas job_noneur=job-20221124-035107-102 - downloaded
# job_eur=job-20221124-035203-109
# job_noneur=job-20221124-035107-102


# i=italy_gsa
# # eur
# # italy_gsa
# # QC  italy_gsa job_eur=job-20221123-160534-810 - downloaded
# # noneur
# # italy_gsa
# # QC  italy_gsa job_noneur=job-20221123-160538-030 - downloaded
# job_eur=job-20221123-160534-810
# job_noneur=job-20221123-160538-030


############# 

for j in ${ancestry[@]}
do rm /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/${j}/*.vcf.gz
done

for j in ${ancestry[@]}
do rm /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/${j}/chr*.info.gz
done

for j in ${ancestry[@]}
do rm /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/${j}/chr_*.zip
done

for j in ${ancestry[@]}
do rm /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/${j}/results.md5
done

rm /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/logs/nohup.out
rm /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/logs/nohup_ch23.out
rm /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/logs/donwload_imputed_results_chr23_*
rm /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/logs/status_*

{
  echo ${i}
  
  j=eur
  TOKEN='YOUR_API_TOKEN';
  curl -H "X-Auth-Token: $TOKEN" https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/${job_eur} \
  > ${path_gwas}imputed/${i}/2022/logs/status_${i}_${j}_qc.json
  
  j=noneur
  TOKEN='YOUR_API_TOKEN';
  curl -H "X-Auth-Token: $TOKEN" https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/${job_noneur} \
  > ${path_gwas}imputed/${i}/2022/logs/status_${i}_${j}_qc.json
  
  /software/R-4.3.1/bin/Rscript \
  /path/to/user/scripts/IIBDGC/create_imputation_download_commnads.R ${i} \
  > ${path_gwas}imputed/${i}/2022/logs/create_imputation_download_commnads_${i}.R_out
  
  ### edit files:
  j=eur
  sed -i '2 i\path=\/lustre\/scratch123\/hgi\/projects\/ibdgwas\/IIBDGC\/imputed\/${i}\/2022\/' ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
  sed -i '2 i\j='${j} ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
  sed -i '2 i\i='${i} ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
  
  j=noneur
  sed -i -e '1,4d' ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
  sed -i '2 i\j='${j} ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
  
  cat  ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_eur  ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_noneur \
  > ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}
  rm ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_*

  # # if noneur only:
  # j=noneur
  # TOKEN='YOUR_API_TOKEN';
  # curl -H "X-Auth-Token: $TOKEN" https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/${job_noneur} \
  # > ${path_gwas}imputed/${i}/2022/logs/status_${i}_${j}_qc.json
  # 
  # /software/R-4.3.1/bin/Rscript \
  # /path/to/user/scripts/IIBDGC/create_imputation_download_commnads.R ${i} \
  # > ${path_gwas}imputed/${i}/2022/logs/create_imputation_download_commnads_${i}.R_out
  # 
  # j=noneur
  # sed -i '2 i\path=\/lustre\/scratch123\/hgi\/projects\/ibdgwas\/IIBDGC\/imputed\/${i}\/2022\/' ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
  # sed -i '2 i\j='${j} ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
  # sed -i '2 i\i='${i} ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
  # 
  # cat ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_noneur \
  # > ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}
  # rm ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_*
  
  sed -i -e 's/-o .*.zip --create-dirs/-P ${path}\/${j}/g' ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}
  sed -i -e 's/-o .*.md5 --create-dirs/-P ${path}\/${j}/g' ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}
  sed -i -e 's/curl -L/wget --tries=75 -c/g' ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}
  # more ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}
  
  cd /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/logs/
  nohup bash /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/logs/donwload_imputed_results_${i} > nohup_new.out &
    
  unset job_eur
  unset job_noneur
  
  
}


####################################################
# 23.16  Unzip files:

studies=(all_hce)

path_gwas=/path/to/ibdgwas/IIBDGC/
  
ancestry=(eur noneur)
MEM=200
MEM=800

### CHR1 TO CHR22
for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
for chr in {1..22}
do
bsub -J"pl2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_6_post_server_qc_${j}_${i}_chr${chr} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_6_post_server_qc_${j}_${i}_chr${chr} \
"unzip -P PkSw6GlBg5VGuc ${path_gwas}imputed/${i}/2022/${j}/chr_${chr}.zip -d ${path_gwas}imputed/${i}/2022/${j}/"
done 
done
done

### CHR23
for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
bsub -J"pl2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_6_post_server_qc_${j}_${i}_chr23 \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_6_post_server_qc_${j}_${i}_chr23 \
"unzip -P PkSw6GlBg5VGuc ${path_gwas}imputed/${i}/2022/${j}/chr_X.zip -d ${path_gwas}imputed/${i}/2022/${j}/"
done
done

chr=1
for i in ${studies[@]}
do echo ${i} && for j in ${ancestry[@]}
do echo ${j} && ls -la ${path_gwas}imputed/${i}/2022/${j}/chr${chr}.info.gz
done
done

###### remove the .zip file for the studies unzipped - double check first unzipped files exist

for i in ${studies[@]}
do echo ${i} && for j in ${ancestry[@]}
do echo ${j} && ls -la ${path_gwas}imputed/${i}/2022/${j}/chr*.vcf.gz | wc -l && \
ls -la ${path_gwas}imputed/${i}/2022/${j}/chr*.info.gz | wc -l && \
ls -la ${path_gwas}imputed/${i}/2022/${j}/chr9.info.gz && \
ls -la ${path_gwas}imputed/${i}/2022/${j}/chrX.dose.vcf.gz && echo "        "
done
done

for i in ${studies[@]}
do echo ${i} && for j in ${ancestry[@]}
do echo ${j} && ls -la ${path_gwas}imputed/${i}/2022/${j}/chr_X.zip && ls -la ${path_gwas}imputed/${i}/2022/${j}/chr_9.zip
done
done

for i in ${studies[@]}
do echo ${i} && for j in ${ancestry[@]}
do rm ${path_gwas}imputed/${i}/2022/${j}/chr_*.zip
done
done



####################################################
# 23.11 create imputation summary files and plots:

for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
mkdir -p ${path_gwas}imputed/${i}/2022/${j}/plots/
done
done

path_gwas="/path/to/ibdgwas/IIBDGC/"

# CREATE final rainfall plots:
for i in ${studies[@]}
do for chr in {1..2}
do bash ${path_gwas}scripts/post_imputation_summary_1.sh ${i} ${chr}
done
done

## double check results
for i in ${studies[@]}
do echo ${i} && for chr in {1..23}
do echo ${chr} && tail -50 ${path_gwas}scripts/logs/post_imputation_summary_stdout_${i}_${chr} | grep -E "completed"
done
done

# CREATE final summary tables of imputation accuracy:
for i in ${studies[@]}
do bash ${path_gwas}scripts/post_imputation_summary_3.sh ${i}
done


### TO CHECK RUNNING FOR 

for i in ${studies[@]}
do echo ${i} && tail -50 /path/to/ibdgwas/IIBDGC/scripts/logs/post_imputation_summary_stdout_${i}_allchr | grep -E "completed"
done

# ### visually inspect plots:
# for i in ${studies[@]}
# done


### output number of samples in post-imputation files:

MEM=500
for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
bsub -J"pl2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_7_post_server_imputation_${j}_${i}_chr23 \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_7_post_server_imputation_${j}_${i}_chr23 \
"/path/to/software/bcftools-1.16/./bcftools \
query -l ${path_gwas}imputed/${i}/2022/${j}/chrX.dose.vcf.gz > ${path_gwas}imputed/test/${i}_${j}_2020_samples_imputed_files_chr23"
done
done

MEM=500
for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
bsub -J"pl2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_7_post_server_imputation_${j}_${i}_chr9 \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_7_post_server_imputation_${j}_${i}_chr9 \
"/path/to/software/bcftools-1.16/./bcftools \
query -l ${path_gwas}imputed/${i}/2022/${j}/chr9.dose.vcf.gz > ${path_gwas}imputed/test/${i}_${j}_2020_samples_imputed_files_chr9"
done
done


for i in ${studies[@]}
do echo ${i} && for j in ${ancestry[@]}
do tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_7_post_server_imputation_${j}_${i}_chr23 | grep -E "completed" && \
tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_7_post_server_imputation_${j}_${i}_chr9 | grep -E "completed" 
done
done

###### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("australia_omniexome","gwas1","gwas2","all_hce"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa","kiel_austria_sibdcs_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","swedish_uc_old_gwas"
           ,"mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

length(cohorts)
# [1] 58

ancestry<-c("eur","noneur")

dat<-as.data.frame(matrix(nrow=length(cohorts),ncol=5))
colnames(dat)<-c("cohort","eur_nsamples","eur_nvariants","noneur_nsamples","noneur_nvariants")
dat$cohort<-cohorts

for (i in 1:length(cohorts)) {
  
  print(cohorts[i])
  
  for (j in 1:length(ancestry)) {
    
    print(ancestry[j])
    
    file_tmp1<-paste(path,"imputed/test/",cohorts[i],"_",ancestry[j],"_2020_samples_imputed_files_chr9",sep="")
    # file_tmp2<-paste(path,"imputed/test/",cohorts[i],"_",ancestry[j],"_2020_samples_imputed_files_chr23",sep="")
    
    if(file.exists(file_tmp1)) {
      tmp<-fread(file_tmp1,head=F)
      print(nrow(tmp))
      # tmp<-fread(file_tmp2,head=F)
      # print(nrow(tmp))
      
      dat[which(dat$cohort==cohorts[i]),paste(ancestry[j],"nsamples",sep="_")]<-nrow(tmp)
      
    } else{
      print("File not ready")
    }
    rm(tmp)
    
    print(ancestry[j])
    a<-tryCatch(read.table(paste(path,"imputation_ready/test/count_n_variants_",cohorts[i],"_",ancestry[j],sep=""),head=F), error=function(e) NULL)
    if(!is.null(a)) {
      print(sum(a$V1))
      dat[which(dat$cohort==cohorts[i]),paste(ancestry[j],"nvariants",sep="_")]<-sum(a$V1)
      
    } else {
      print("No variants")
    }
    rm(a)

  }
}

write.table(dat,paste(path,"imputed/plots/number_samples_variants_per_ancestry_preqc_2022.tsv",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
# NOTE: file shared with Phil - updated to just show number of cases from all_hce

################
# [1] "australia_omniexome"
# [1] "eur"
# [1] 1245
# [1] 1245
# [1] "noneur"
# [1] 54
# [1] 54
# [1] "gwas1"
# [1] "eur"
# [1] 4652
# [1] 4652
# [1] "noneur"
# [1] 7
# [1] 7
# [1] "gwas2"
# [1] "eur"
# [1] 7758
# [1] 7758
# [1] "noneur"
# [1] 15
# [1] 15
# [1] "all_hce"
# [1] "eur"
# [1] 22588
# [1] 22588
# [1] "noneur"
# [1] 1175
# [1] 1175
# [1] "pittsburgh_gsa"
# [1] "eur"
# [1] 2701
# [1] 2701
# [1] "noneur"
# [1] 8
# [1] 8
# [1] "spain_gsa"
# [1] "eur"
# [1] 3396
# [1] 3396
# [1] "noneur"
# [1] 12
# [1] 12
# [1] "italy_gsa"
# [1] "eur"
# [1] 938
# [1] 938
# [1] "noneur"
# [1] 5
# [1] 5
# [1] "kiel_austria_sibdcs_gsa"
# [1] "eur"
# [1] 14086
# [1] 14086
# [1] "noneur"
# [1] 243
# [1] 243
# [1] "netherlands_gsa"
# [1] "eur"
# [1] 4315
# [1] 4315
# [1] "noneur"
# [1] 222
# [1] 222
# [1] "slovenia_gsa"
# [1] "eur"
# [1] 263
# [1] 263
# [1] "noneur"
# [1] 0
# [1] 0
# [1] "sweden_gsa"
# [1] "eur"
# [1] 1364
# [1] 1364
# [1] "noneur"
# [1] 32
# [1] 32
# [1] "niddk_broad_gsa"
# [1] "eur"
# [1] 5130
# [1] 5130
# [1] "noneur"
# [1] 296
# [1] 296
# [1] "niddk_feinstein_gsa"
# [1] "eur"
# [1] 6965
# [1] 6965
# [1] "noneur"
# [1] 1159
# [1] 1159
# [1] "basque_gsa"
# [1] "eur"
# [1] 1487
# [1] 1487
# [1] "noneur"
# [1] 11
# [1] 11
# [1] "prism_nfe_gsa"
# [1] "eur"
# [1] 426
# [1] 426
# [1] "noneur"
# [1] 36
# [1] 36
# [1] "lithuania_gsa"
# [1] "eur"
# [1] 2213
# [1] 2213
# [1] "noneur"
# [1] 12
# [1] 12
# [1] "belgium_louis_gsa"
# [1] "eur"
# [1] 1499
# [1] 1499
# [1] "noneur"
# [1] 12
# [1] 12
# [1] "belgium_franchimont_gsa"
# [1] "eur"
# [1] 1459
# [1] 1459
# [1] "noneur"
# [1] 32
# [1] 32
# [1] "belgium_vermeire_gsa"
# [1] "eur"
# [1] 3894
# [1] 3894
# [1] "noneur"
# [1] 80
# [1] 80
# [1] "prism_nfe_gwas"
# [1] "eur"
# [1] 717
# [1] 717
# [1] "noneur"
# [1] 93
# [1] 93
# [1] "finland_illugwas"
# [1] "eur"
# [1] 440
# [1] 440
# [1] "noneur"
# [1] 4
# [1] 4
# [1] "german_affy6_old_gwas"
# [1] "eur"
# [1] 2787
# [1] 2787
# [1] "noneur"
# [1] 10
# [1] 10
# [1] "norway_affy6_old_gwas"
# [1] "eur"
# [1] 542
# [1] 542
# [1] "noneur"
# [1] 2
# [1] 2
# [1] "belgium_inf1_old_gwas"
# [1] "eur"
# [1] 1396
# [1] 1396
# [1] "noneur"
# [1] 10
# [1] 10
# [1] "belgium_inf2_old_gwas"
# [1] "eur"
# [1] 270
# [1] 270
# [1] "noneur"
# [1] 1
# [1] 1
# [1] "niddk_old_gwas"
# [1] "eur"
# [1] 2712
# [1] 2712
# [1] "noneur"
# [1] 20
# [1] 20
# [1] "cedars_370k_old_gwas"
# [1] "eur"
# [1] 509
# [1] 509
# [1] "noneur"
# [1] 88
# [1] 88
# [1] "cedars_610k_old_gwas"
# [1] "eur"
# [1] 825
# [1] 825
# [1] "noneur"
# [1] 56
# [1] 56
# [1] "cedars_omni_old_gwas"
# [1] "eur"
# [1] 1104
# [1] 1104
# [1] "noneur"
# [1] 103
# [1] 103
# [1] "swedish_uc_old_gwas"
# [1] "eur"
# [1] 1255
# [1] 0
# [1] "noneur"
# [1] 3
# [1] 0
# [1] "mccauley_gsa"
# [1] "eur"
# [1] 758
# [1] 758
# [1] "noneur"
# [1] 18
# [1] 18
# [1] "ccfa_gsa"
# [1] "eur"
# [1] 1932
# [1] 1932
# [1] "noneur"
# [1] 235
# [1] 235
# [1] "cedars_gsa"
# [1] "eur"
# [1] 2699
# [1] 2699
# [1] "noneur"
# [1] 360
# [1] 360
# [1] "bernstein_gsa"
# [1] "eur"
# [1] 482
# [1] 482
# [1] "noneur"
# [1] 27
# [1] 27
# [1] "farkkila_gsa"
# [1] "eur"
# [1] 68
# [1] 68
# [1] "noneur"
# [1] 0
# [1] 0
# [1] "franchimont_gsa"
# [1] "eur"
# [1] 2360
# [1] 2360
# [1] "noneur"
# [1] 409
# [1] 409
# [1] "franke_gsa"
# [1] "eur"
# [1] 853
# [1] 853
# [1] "noneur"
# [1] 11
# [1] 11
# [1] "helmsley_prism_gsa"
# [1] "eur"
# [1] 669
# [1] 669
# [1] "noneur"
# [1] 84
# [1] 84
# [1] "helmsley_xavier_prism_gsa"
# [1] "eur"
# [1] 1171
# [1] 1171
# [1] "noneur"
# [1] 110
# [1] 110
# [1] "hyams_protect_gsa"
# [1] "eur"
# [1] 323
# [1] 323
# [1] "noneur"
# [1] 93
# [1] 93
# [1] "lewis_sparc_gsa"
# [1] "eur"
# [1] 2529
# [1] 2529
# [1] "noneur"
# [1] 317
# [1] 317
# [1] "mccauley_new_gsa"
# [1] "eur"
# [1] 1323
# [1] 1323
# [1] "noneur"
# [1] 287
# [1] 287
# [1] "mcgovern_gsa"
# [1] "eur"
# [1] 5149
# [1] 5149
# [1] "noneur"
# [1] 829
# [1] 829
# [1] "moayyedi_imagine_gsa"
# [1] "eur"
# [1] 1011
# [1] 1011
# [1] "noneur"
# [1] 109
# [1] 109
# [1] "newberry_share_gsa"
# [1] "eur"
# [1] 747
# [1] 747
# [1] "noneur"
# [1] 111
# [1] 111
# [1] "niddk_cho_gsa"
# [1] "eur"
# [1] 1611
# [1] 1611
# [1] "noneur"
# [1] 143
# [1] 143
# [1] "niddk_duerr_gsa"
# [1] "eur"
# [1] 1868
# [1] 1868
# [1] "noneur"
# [1] 69
# [1] 69
# [1] "niddk_rioux_gsa"
# [1] "eur"
# [1] 871
# [1] 871
# [1] "noneur"
# [1] 42
# [1] 42
# [1] "niddk_silverberg_gsa"
# [1] "eur"
# [1] 2112
# [1] 2112
# [1] "noneur"
# [1] 235
# [1] 235
# [1] "palotie_hus_gsa"
# [1] "eur"
# [1] 853
# [1] 853
# [1] "noneur"
# [1] 16
# [1] 16
# [1] "pekow_share_gsa"
# [1] "eur"
# [1] 538
# [1] 538
# [1] "noneur"
# [1] 94
# [1] 94
# [1] "rioux_igenomed_gsa"
# [1] "eur"
# [1] 171
# [1] 171
# [1] "noneur"
# [1] 10
# [1] 10
# [1] "sands_msccr_gsa"
# [1] "eur"
# [1] 1120
# [1] 1120
# [1] "noneur"
# [1] 299
# [1] 299
# [1] "stampfer_gsa"
# [1] "eur"
# [1] 1439
# [1] 1439
# [1] "noneur"
# [1] 25
# [1] 25
# [1] "vermeire_gsa"
# [1] "eur"
# [1] 4540
# [1] 4540
# [1] "noneur"
# [1] 119
# [1] 119
# [1] "weersma_gsa"
# [1] "eur"
# [1] 682
# [1] 682
# [1] "noneur"
# [1] 23
# [1] 23
# [1] "xavier_prism_gsa"
# [1] "eur"
# [1] 619
# [1] 619
# [1] "noneur"
# [1] 67
# [1] 67
# [1] "xavier_share_gsa"
# [1] "eur"
# [1] 637
# [1] 637
# [1] "noneur"
# [1] 55
# [1] 55
################

####################################################
# 23.12 Final double check to EmpRsq

# still some genotyped variants with EmpRsq<0.5, double check distribution of these per study:

# final rainfall plots:
for i in ${studies[@]}
do bash ${path_gwas}scripts/post_imputation_summary_4.sh ${i}
done

for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}scripts/logs/post_imputation_summary_4_stdout_${i}_allchr  | grep -E "completed"
done

###### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"

cohorts<-c("australia_omniexome","gwas1"
           ,"pittsburgh_gsa","spain_gsa","italy_gsa"
           ,"netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa","niddk_feinstein_gsa","basque_gsa","prism_nfe_gsa"
           ,"lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa","belgium_vermeire_gsa"
           ,"prism_nfe_gwas","finland_illugwas"
           ,"german_affy6_old_gwas","norway_affy6_old_gwas"
           ,"belgium_inf1_old_gwas","belgium_inf2_old_gwas","niddk_old_gwas"
           ,"cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas","swedish_uc_old_gwas"
           ,"mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","xavier_prism_gsa",
           "xavier_share_gsa","gwas2","all_hce","kiel_austria_sibdcs_gsa","vermeire_gsa","weersma_gsa")

new_cohorts<-c("bernstein_gsa","farkkila_gsa","franchimont_gsa",
               "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
               "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
               "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
               "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","xavier_prism_gsa",
               "xavier_share_gsa")


length(cohorts)

# [1] 58

ancestry<-"eur"
j=1

dat<-as.data.frame(matrix(ncol=4,nrow=length(cohorts)))
colnames(dat)<-c("cohort","n_genotyped_variants","n_variants_empRsq_lower_0.5","n_variants_empRsq_lower_0.4")
dat$cohort<-cohorts

for (i in 1:length(cohorts)) {
  
  print(cohorts[i])
  
  file_tmp<-paste(path,"imputed/",cohorts[i],"/2022/eur/plots/",cohorts[i],"_eur_low_EmpRsq_after_first_roudn_exclussions.tsv",sep="")
  
  if(file.exists(file_tmp)) {
    
    tmp<-read.table(file_tmp,head=T)
    tmp$study<-cohorts[i]
    
    if(!exists("all_tmp")) {
      all_tmp<-tmp
    }else{
      all_tmp<-rbind(all_tmp,tmp)
    }
    tmp_nvar<-read.table(paste(path,"imputation_ready/test/count_n_variants_",cohorts[i],"_",ancestry[j],sep=""),head=F)
    
    dat$n_genotyped_variants[which(dat$cohort==cohorts[i])]<-sum(tmp_nvar$V1)
    dat$n_variants_empRsq_lower_0.5[which(dat$cohort==cohorts[i])]<-nrow(tmp)
    dat$n_variants_empRsq_lower_0.4[which(dat$cohort==cohorts[i])]<-nrow(tmp[which(tmp$EmpRsq<0.4),])
    rm(tmp,tmp_nvar)
    
  }
}

dat$percentage_empRsq_lower_0.5<-dat$n_variants_empRsq_lower_0.5/dat$n_genotyped_variants
dat<-dat[order(dat$percentage_empRsq_lower_0.5,decreasing=T),]

dat[which( (dat$n_variants_empRsq_lower_0.5/dat$n_genotyped_variants)>=0.001),]
# cohort n_genotyped_variants n_variants_empRsq_lower_0.5
# 5 italy_gsa               475735                        1152
# n_variants_empRsq_lower_0.4 percentage_empRsq_lower_0.5
# 5                         126                 0.002421516

write.table(dat,paste(path,"imputed/plots/low_EmpRsq_after_second_final_roudn_exclussions_2022.tsv",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

q("no")



#####################

# remove old imputed files - there is a copy archived in irods:

studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

for i in ${studies[@]}
do rm -r ${path_gwas}imputed/${i}/noneur/
done

################################################################################################################################################
################################################################################################################################################

# one sample with low call rate at PAR2 within niddk_broad_gsa is causing that chunk to be excluded from imputation - re-do (see 'pseudoautosomal_region_coverage')

i=niddk_broad_gsa
emacs  ${path_gwas}pre_imputation/QC/${i}/${i}_list_samples_exclude_from_chr23
# sample_id sample_id



###########################################
# 23.13 Convert ped/map files to VCF files #
###########################################


studies=(niddk_broad_gsa)

path_gwas=/path/to/ibdgwas/IIBDGC/
MEM=1000

ancestry=(eur)


chr=23

for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_8_post_server_qc_${j}_${i}_chr${chr} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_8_post_server_qc_${j}_${i}_chr${chr} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_postqc_lifted_hg38_nomonom_RefAlt_posstrandaligned_updated_TOPMed_flip2_freqCheck_noflagged_${j}_nomonom \
--allow-no-sex \
--keep-allele-order \
--flip ${path_gwas}imputed/${i}/qc/list_variants_negative_EmpR_toflip \
--exclude ${path_gwas}imputed/${i}/qc/2022/list_variants_to_exclude_2 \
--remove ${path_gwas}pre_imputation/QC/${i}/${i}_list_samples_exclude_from_chr23 \
--chr ${chr} \
--output-chr chr26 --recode vcf-iid --out ${path_gwas}imputation_ready/${i}/${i}_${j}_chr${chr}_2022"
done
done



for i in ${studies[@]}
do
echo ${i} && for j in ${ancestry[@]}
do
echo ${j} && tail -40 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_8_post_server_qc_${j}_${i}_chr${chr} | grep -E "completed"
done
done

for i in ${studies[@]}
do ls -la ${path_gwas}imputation_ready/${i}/${i}_${j}_chr${chr}_2022.vcf
done


########################################################
# 23.14 Create a sorted vcf.gz - final round imputation #
########################################################


for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
for chr in 23
do
bsub -J"gnomad" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_9_post_server_qc_${j}_${i}_chr${chr} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_9_post_server_qc_${j}_${i}_chr${chr} \
"/path/to/software/bcftools-1.16/./bcftools sort ${path_gwas}imputation_ready/${i}/${i}_${j}_chr${chr}_2022.vcf \
-Oz -o ${path_gwas}imputation_ready/${i}/${i}_${j}_chr${chr}_2022.vcf.gz"
done
done
done

# rm intermediate files:
for i in ${studies[@]}
do rm ${path_gwas}imputation_ready/${i}/${i}_*_chr*.vcf
done
for i in ${studies[@]}
do rm ${path_gwas}imputation_ready/${i}/*.hh
done
for i in ${studies[@]}
do rm ${path_gwas}imputation_ready/${i}/*.log
done
for i in ${studies[@]}
do rm ${path_gwas}imputation_ready/${i}/*.nosex
done

for i in ${studies[@]}
do for j in ${ancestry[@]}
do echo ${i} && ls -la ${path_gwas}imputation_ready/${i}/${i}_${j}_chr*_2022.vcf.gz |  wc -l
done
done

for i in ${studies[@]}
do for j in ${ancestry[@]}
do echo ${i} && ls -la ${path_gwas}imputation_ready/${i}/${i}_${j}_chr1_2022.vcf.gz
done
done

##################################################################################################################################################################
##################################################################################################################################################################

####################################################
# 23.15 submit second and FINAL round of imputation #
####################################################

path_gwas=/path/to/ibdgwas/IIBDGC/
ancestry=(eur)

studies=()

# new server - updated to v1.7.1
i=niddk_broad_gsa
for j in ${ancestry[@]}
do echo ${i} && ls -la ${path_gwas}imputation_ready/${i}/${i}_${j}_chr*_2022.vcf.gz |  wc -l
done
for j in ${ancestry[@]}
do echo ${j} && python2 /path/to/user/scripts/IIBDGC/topmed_imputation_${j}.py ${i}
done


# eur
# niddk_broad_gsa
# Your job was successfully added to the job queue.
# QC  niddk_broad_gsa job_eur=job-20230228-094347-334
i=niddk_broad_gsa
j=eur
job_eur=job-20230228-094347-334

for j in ${ancestry[@]}
do rm /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/${j}/chrX.info.gz
done

for j in ${ancestry[@]}
do rm /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/${j}/chr_X*.zip
done


rm /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/logs/nohup_*
rm /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/logs/donwload_imputed_results_*
rm /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/logs/status_*
  
{
    echo ${i}
    
    j=eur
    TOKEN='YOUR_API_TOKEN';
    curl -H "X-Auth-Token: $TOKEN" https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/${job_eur} \
    > ${path_gwas}imputed/${i}/2022/logs/status_${i}_${j}_qc.json
    
    j=noneur
    TOKEN='YOUR_API_TOKEN';
    curl -H "X-Auth-Token: $TOKEN" https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/${job_noneur} \
    > ${path_gwas}imputed/${i}/2022/logs/status_${i}_${j}_qc.json
    
    /software/R-4.3.1/bin/Rscript \
    /path/to/user/scripts/IIBDGC/create_imputation_download_commnads.R ${i} \
    > ${path_gwas}imputed/${i}/2022/logs/create_imputation_download_commnads_${i}.R_out
    
    ### edit files:
    j=eur
    sed -i '2 i\path=\/lustre\/scratch123\/hgi\/projects\/ibdgwas\/IIBDGC\/imputed\/${i}\/2022\/' ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
    sed -i '2 i\j='${j} ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
    sed -i '2 i\i='${i} ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
    
    j=noneur
    sed -i -e '1,4d' ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
    sed -i '2 i\j='${j} ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
    
    cat  ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_eur  ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_noneur \
    > ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}
    rm ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_*
      
      # # if noneur only:
      # j=noneur
      # TOKEN='YOUR_API_TOKEN';
      # curl -H "X-Auth-Token: $TOKEN" https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/${job_noneur} \
      # > ${path_gwas}imputed/${i}/2022/logs/status_${i}_${j}_qc.json
      # 
      # /software/R-4.3.1/bin/Rscript \
      # /path/to/user/scripts/IIBDGC/create_imputation_download_commnads.R ${i} \
      # > ${path_gwas}imputed/${i}/2022/logs/create_imputation_download_commnads_${i}.R_out
      # 
    # j=noneur
    # sed -i '2 i\path=\/lustre\/scratch123\/hgi\/projects\/ibdgwas\/IIBDGC\/imputed\/${i}\/2022\/' ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
    # sed -i '2 i\j='${j} ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
    # sed -i '2 i\i='${i} ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_${j}
    # 
    # cat ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_noneur \
    # > ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}
    # rm ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}_*
    
    sed -i -e 's/-o .*.zip --create-dirs/-P ${path}\/${j}/g' ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}
    sed -i -e 's/-o .*.md5 --create-dirs/-P ${path}\/${j}/g' ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}
    sed -i -e 's/curl -L/wget --tries=75 -c/g' ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}
    # more ${path_gwas}imputed/${i}/2022/logs/donwload_imputed_results_${i}
    
    cd /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/logs/
      nohup bash /path/to/ibdgwas/IIBDGC/imputed/${i}/2022/logs/donwload_imputed_results_${i} > nohup_new.out &
      
      unset job_eur
    unset job_noneur
    
    
  }

##############
##############


studies=(niddk_broad_gsa)

for i in ${studies[@]}
do echo ${i} && for j in ${ancestry[@]}
do echo ${j} && ls -la ${path_gwas}imputed/${i}/2022/${j}/chr_X.zip
done
done

##############
##############

####################################################
# 23.10  Unzip files:

path_gwas=/path/to/ibdgwas/IIBDGC/
  
studies=(niddk_broad_gsa)

ancestry=(eur)
MEM=200


### CHR23
for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
bsub -J"pl2" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G ibdgwas -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_imputation_ready_10_post_server_qc_${j}_${i}_chr23 \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_imputation_ready_10_post_server_qc_${j}_${i}_chr23 \
"unzip -P PkSw6GlBg5VGuc ${path_gwas}imputed/${i}/2022/${j}/chr_X.zip -d ${path_gwas}imputed/${i}/2022/${j}/"
done
done


###### remove the .zip file for the studies unzipped - double check first unzipped files exist

for i in ${studies[@]}
do echo ${i} && for j in ${ancestry[@]}
do echo ${j} && ls -la ${path_gwas}imputed/${i}/2022/${j}/chr*.vcf.gz | wc -l && \
ls -la ${path_gwas}imputed/${i}/2022/${j}/chrX.info.gz | wc -l && \
ls -la ${path_gwas}imputed/${i}/2022/${j}/chrX.dose.vcf.gz && echo "        "
done
done

for i in ${studies[@]}
do echo ${i} && for j in ${ancestry[@]}
do echo ${j} && ls -la ${path_gwas}imputed/${i}/2022/${j}/chr_X.zip
done
done

for i in ${studies[@]}
do echo ${i} && for j in ${ancestry[@]}
do rm ${path_gwas}imputed/${i}/2022/${j}/chr_*.zip
done
done



####################################################
# 23.11 create imputation summary files and plots:

for i in ${studies[@]}
do
for j in ${ancestry[@]}
do
mkdir -p ${path_gwas}imputed/${i}/2022/${j}/plots/
  done
done

path_gwas="/path/to/ibdgwas/IIBDGC/"

# CREATE final rainfall plots:
for i in ${studies[@]}
do for chr in 23
do bash ${path_gwas}scripts/post_imputation_summary_1.sh ${i} ${chr}
done
done


## double check results
for i in ${studies[@]}
do echo ${i} && for chr in 23
do echo ${chr} && tail -50 ${path_gwas}scripts/logs/post_imputation_summary_stdout_${i}_${chr} | grep -E "completed"
done
done

# CREATE final summary tables of imputation accuracy:
for i in ${studies[@]}
do bash ${path_gwas}scripts/post_imputation_summary_3.sh ${i}
done

### TO CHECK RUNNING FOR 


for i in ${studies[@]}
do echo ${i} && tail -50 /path/to/ibdgwas/IIBDGC/scripts/logs/post_imputation_summary_stdout_${i}_allchr | grep -E "completed"
done

# ### visually inspect plots:
# for i in ${studies[@]}
# done

