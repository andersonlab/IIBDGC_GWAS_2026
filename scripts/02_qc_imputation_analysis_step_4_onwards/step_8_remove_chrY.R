# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
###################
# 8.- REMOVE ChrY #
###################

### R



##############################
### /software/R-4.3.1/bin/R

library(data.table)
library(ggplot2)
library(ggpubr)
library(gtools)

path<-"/path/to/ibdgwas/IIBDGC/"
studies<-c("all_hce","niddk_old_gwas","australia_omniexome","gwas1","gwas2","pittsburgh_gsa","spain_gsa","italy_gsa",
           "kiel_austria_sibdcs_gsa","netherlands_gsa","slovenia_gsa","sweden_gsa","niddk_broad_gsa",
           "niddk_feinstein_gsa","basque_gsa","lithuania_gsa","belgium_louis_gsa","belgium_franchimont_gsa",
           "belgium_vermeire_gsa","prism_nfe_gsa","prism_nfe_gwas","finland_illugwas","german_affy6_old_gwas",
           "norway_affy6_old_gwas","belgium_inf1_old_gwas","belgium_inf2_old_gwas",
           "cedars_370k_old_gwas","cedars_610k_old_gwas","cedars_omni_old_gwas",
           "swedish_uc_old_gwas",
           "mccauley_gsa","ccfa_gsa","cedars_gsa","bernstein_gsa","farkkila_gsa","franchimont_gsa",
           "franke_gsa","helmsley_prism_gsa","helmsley_xavier_prism_gsa","hyams_protect_gsa","lewis_sparc_gsa",
           "mccauley_new_gsa","mcgovern_gsa","moayyedi_imagine_gsa","newberry_share_gsa","niddk_cho_gsa",
           "niddk_duerr_gsa","niddk_rioux_gsa","niddk_silverberg_gsa","palotie_hus_gsa","pekow_share_gsa",
           "rioux_igenomed_gsa","sands_msccr_gsa","stampfer_gsa","vermeire_gsa","weersma_gsa","xavier_prism_gsa",
           "xavier_share_gsa")

for (j in 1:length(studies)) {
  
  print(studies[j])
  bim<-read.table(paste(path,"pre_imputation/QC/",studies[j],"/",studies[j],"_hg19_noind_posstr_nodup_flip_sexcheck_missinghh.bim",sep=""),head=F)
  
  print(table(bim$V1))
  list_toremove<-bim[which(bim$V1==24),"V2",drop=F]
  
  write.table(list_toremove,paste(path,"pre_imputation/QC/",studies[j],"/list_chry_variants_toremove",sep=""),col.names=F,row.names=F,sep="\t",quote=F)
  
  rm(list=ls()[!ls() %in% c("j","studies","path")])
}

q("no")

############
# [1] "all_hce"
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 45564 39110 32968 27279 27371 31804 25706 22831 21606 22942 28232 25154 14366 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 16196 16123 18450 19993 11488 20284 12781  6598  8743 12477  1504   148 
# [1] "niddk_old_gwas"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 23045 25082 21377 18907 19055 20577 16514 18099 15683 15422 14491 14853 11418 
# 14    15    16    17    18    19    20    21    22    23    25 
# 9748  8811  8898  8234 10389  5845  7767  5439  5417  8968     2 
# [1] "australia_omniexome"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 65585 61855 50833 43260 44955 52436 40910 39302 35706 41373 40872 38678 28511 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 25482 23999 25687 23829 22246 19526 19997 11079 12056 18586   934   465 
# [1] "gwas1"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 36370 37719 30904 29279 29458 28834 23595 25050 20990 26109 23908 22771 17558 
# 14    15    16    17    18    19    20    21    22    23 
# 14297 13032 13856 10226 13541  5775 11311  6503  5616  9746 
# [1] "gwas2"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 61139 63444 52679 48260 49484 49509 41386 42181 36541 41830 38596 37071 29569 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 24488 22450 23851 17917 22761 10455 19813 11004 10128 36382   206   410 
# [1] "pittsburgh_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 77688 73255 59825 55835 54124 71132 49577 49261 43620 49681 46645 45205 33195 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 28706 27878 29759 26310 26407 20823 25491 13375 13450 22362   999   374 
# [1] "spain_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 44608 47320 39393 35250 35682 39984 31768 32400 27486 30432 28372 28220 21887 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 19058 17348 17381 15123 17237 10123 14517  8455  8748 14735    93    41 
# [1] "italy_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 44589 46739 38604 36143 33399 40656 31380 29079 24095 28119 28015 26965 20089 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 18168 17148 18633 16775 16269 13037 13654  7861  8244 24945  1698   644 
# [1] "kiel_austria_sibdcs_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 51329 53385 43960 40468 37669 45277 35263 32323 27354 31455 32224 30616 22052 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 20534 19665 21293 19895 18047 15917 15351  8905  9542 26756  2007   697 
# [1] "netherlands_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 54420 55672 46392 42043 39607 47118 37208 33811 28783 32783 34163 32269 22819 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 21699 20863 22461 21446 18747 17148 16120  9370 10006 17412  1346   526 
# [1] "slovenia_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 42972 44882 37838 35450 32783 39939 30363 28275 23611 27301 26971 26085 19482 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 17779 16693 18055 16155 16081 12486 13352  7777  7946 13858   764   421 
# [1] "sweden_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 44487 46337 39032 36476 33652 40986 31245 29183 24400 28161 27929 26901 19946 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 18313 17241 18667 16731 16576 12954 13788  8012  8222 14315   806   464 
# [1] "niddk_broad_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 48741 50232 42306 39225 36371 43838 33821 31369 26317 30226 30581 29220 21260 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 19726 18662 20231 18544 17570 14652 14804  8579  8987 15302   685   500 
# [1] "niddk_feinstein_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 50960 52115 43745 40240 37553 45109 34975 32364 27250 31183 31864 30436 21778 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 20480 19438 21189 19623 18036 15825 15330  8879  9447 15451  1002   496 
# [1] "basque_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 44005 46047 38135 35741 33016 40247 30884 28662 23800 27799 27696 26617 19770 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 17982 16973 18465 16560 16124 12932 13521  7799  8125 24702  1713   649 
# [1] "lithuania_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 44914 46790 39435 36774 34019 41159 31520 29331 24654 28366 28177 27199 20135 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 18455 17413 18830 16893 16634 13128 13828  8088  8317 14423   721   442 
# [1] "belgium_louis_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 47232 48681 41116 38038 35317 42570 32853 30383 25704 29408 29585 28367 20651 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 19200 18282 19755 18144 17024 14312 14456  8419  8827 15664   850   471 
# [1] "belgium_franchimont_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 45355 47151 39771 36965 34183 41437 31836 29514 24835 28613 28528 27422 20223 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 18619 17524 18958 17203 16651 13454 13972  8168  8407 14798   735   476 
# [1] "belgium_vermeire_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 53689 54864 45666 41424 38991 46368 36585 33355 28297 32399 33716 31880 22567 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 21403 20591 22182 21134 18599 16845 15932  9173  9897 16934  1154   515 
# [1] "prism_nfe_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 44302 46161 38958 36452 33557 40800 31122 28993 24350 28097 27902 26918 19876 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 18275 17220 18630 16675 16493 12926 13710  7985  8227 14522   718   471 
# [1] "prism_nfe_gwas"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 18859 19813 16749 15681 14850 15223 13350 13042 10743 12503 11734 11627  9433 
# 14    15    16    17    18    19    20    21    22    23 
# 7841  7090  7450  6428  7261  4938  6047  3579  3356  6371 
# [1] "finland_illugwas"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 18416 19404 16406 15403 14610 14999 13180 12826 10594 12252 11498 11412  9263 
# 14    15    16    17    18    19    20    21    22    23 
# 7677  6969  7294  6273  7094  4856  5933  3525  3289  6056 
# [1] "german_affy6_old_gwas"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 69755 72246 59323 54669 55187 55129 46029 47422 40683 47149 43527 41684 33525 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 27451 25474 27047 20207 25971 11692 22338 12360 11323 34415   236   356 
# [1] "norway_affy6_old_gwas"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 65961 68490 56415 51885 52762 52815 44077 45193 38852 44875 41377 39723 31908 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 26198 24126 25514 19048 24548 11129 21144 11790 10746 31902   183   345 
# [1] "belgium_inf1_old_gwas"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 22412 24256 20632 18070 18404 19914 15865 17501 15217 14960 14013 14414 10961 
# 14    15    16    17    18    19    20    21    22    23 
# 9430  8549  8663  8022 10043  5664  7575  5264  5303  7568 
# [1] "belgium_inf2_old_gwas"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 21525 23278 19734 17188 17554 18988 15126 16771 14634 14400 13450 13777 10459 
# 14    15    16    17    18    19    20    21    22    23 
# 9096  8239  8354  7750  9615  5405  7353  5029  5076  7929 
# [1] "cedars_370k_old_gwas"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 25183 26979 23064 20607 20489 24255 18218 19055 16459 16733 15732 15905 12448 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 10471  9502  9398  8868 10789  6451  8035  5601  5689  9924    65    25 
# [1] "cedars_610k_old_gwas"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 44653 47367 39427 35301 35707 40001 31840 32422 27504 30472 28405 28262 21919 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 19092 17388 17387 15131 17245 10161 14528  8467  8763 14735    77    39 
# [1] "cedars_omni_old_gwas"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 58586 57570 47105 40423 42030 47719 37842 37012 32766 39009 36438 35349 27828 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 23344 21582 22503 20018 21686 14849 18290 10253 10470 17670   870   434 
# [1] "mccauley_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 44155 45905 38756 36141 33374 40627 30977 28835 24291 28023 27775 26752 19850 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 18222 17091 18556 16664 16409 12941 13644  7983  8203 14417   694   470 
# [1] "ccfa_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 47843 49313 41554 38467 35653 42959 33227 30757 25929 29740 29950 28747 20921 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 19473 18367 19852 18180 17343 14394 14552  8447  8819 15164   751   488 
# [1] "cedars_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 53098 54302 45392 41363 38806 46201 36306 33254 28189 32223 33222 31580 22396 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 21189 20331 22013 20739 18474 16610 15845  9167  9845 16799  1009   504 
# [1] "bernstein_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 43675 45462 38396 35963 33136 40338 30724 28654 24019 27715 27473 26556 19668 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 18058 16962 18342 16374 16340 12674 13554  7866  8066 14244   630   451 
# [1] "farkkila_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 40543 42436 35894 33742 31115 37510 28765 26837 22519 25988 25498 24769 18579 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 16897 15810 17067 15178 15380 11594 12657  7369  7594 13238   541   402 
# [1] "franchimont_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 53049 54304 45665 41430 39002 46349 36554 33330 28306 32117 33260 31498 22287 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 21219 20361 22057 20909 18388 16737 15894  9182  9884 17028  1179   513 
# [1] "franke_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 44044 46031 38791 36257 33478 40650 31050 28980 24250 28002 27740 26754 19860 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 18208 17100 18475 16566 16450 12880 13678  7944  8169 14212   684   445 
# [1] "helmsley_prism_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 18859 19813 16749 15681 14850 15223 13351 13042 10743 12503 11734 11627  9433 
# 14    15    16    17    18    19    20    21    22    23 
# 7841  7090  7450  6428  7261  4938  6047  3579  3356  6371 
# [1] "helmsley_xavier_prism_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 18860 19809 16749 15680 14848 15219 13345 13041 10743 12503 11733 11625  9432 
# 14    15    16    17    18    19    20    21    22    23 
# 7841  7089  7449  6428  7261  4936  6045  3579  3357  6371 
# [1] "hyams_protect_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 45113 46709 39436 36769 33944 41238 31475 29323 24629 28442 28293 27243 20054 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 18530 17370 18853 17045 16630 13228 13910  8062  8312 14714   689   482 
# [1] "lewis_sparc_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 48643 50087 42090 38899 36109 43427 33708 31091 26257 30104 30445 29184 21156 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 19702 18659 20148 18565 17520 14697 14731  8540  8957 15381   828   491 
# [1] "mccauley_new_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 46009 47575 40048 37220 34382 41763 32035 29701 25067 28891 29026 27748 20417 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 18817 17790 19238 17554 16878 13758 14093  8207  8535 14945   726   482 
# [1] "mcgovern_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 51347 52856 44276 40723 38050 45408 35495 32727 27545 31538 32139 30696 22062 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 20669 19667 21301 19761 18195 15811 15442  8958  9484 15958   912   498 
# [1] "moayyedi_imagine_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 48116 49890 41893 38815 36031 43319 33493 31145 26037 29972 30151 28928 21178 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 19502 18469 19962 18244 17464 14265 14625  8525  8833 14969   751   489 
# [1] "newberry_share_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 46154 47665 40127 37431 34646 41847 32120 29809 25186 28877 28961 27797 20391 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 18822 17778 19263 17507 16851 13845 14132  8226  8561 15057   727   485 
# [1] "niddk_cho_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 47595 49101 41342 38379 35585 42831 33099 30649 25722 29621 29898 28563 20866 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 19322 18261 19779 18133 17251 14215 14517  8434  8775 15106   782   484 
# [1] "niddk_duerr_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 46380 48023 40464 37554 34800 42030 32429 29998 25215 29011 29196 28057 20481 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 18938 17917 19425 17686 16959 13901 14168  8314  8640 14975   780   486 
# [1] "niddk_rioux_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 45872 47559 39990 37277 34457 41637 32076 29811 24972 28794 28773 27702 20347 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 18745 17709 19086 17306 16850 13553 14107  8182  8453 14769   747   480 
# [1] "niddk_silverberg_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 49186 50767 42569 39378 36640 43898 34105 31575 26530 30444 30828 29490 21389 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 19914 18863 20449 18695 17689 14794 14886  8641  9032 15313   822   490 
# [1] "palotie_hus_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 43694 45459 38378 35925 33155 40240 30691 28534 24049 27684 27422 26531 19680 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 18048 16926 18332 16388 16247 12697 13532  7882  8118 14252   686   452 
# [1] "pekow_share_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 45136 46765 39414 36741 33978 41197 31527 29323 24665 28358 28349 27329 20076 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 18482 17475 18893 17026 16606 13317 13857  8081  8362 14678   685   480 
# [1] "rioux_igenomed_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 42515 44287 37447 35047 32394 39462 29923 27909 23456 27115 26680 25938 19259 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 17586 16546 17884 15912 15920 12298 13248  7706  7847 13960   608   446 
# [1] "sands_msccr_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 47141 48824 40944 38072 35298 42615 32844 30472 25617 29455 29567 28394 20701 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 19228 18138 19629 17902 17165 14095 14426  8366  8729 14937   706   484 
# [1] "stampfer_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 53442 54540 45567 41401 38917 46294 36514 33241 28286 32272 33426 31714 22421 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 21333 20447 22076 21031 18485 16852 15867  9217  9881 17138  1111   512 
# [1] "vermeire_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 53722 54906 45707 41484 39029 46406 36640 33387 28337 32421 33736 31895 22583 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 21409 20606 22198 21150 18610 16855 15942  9188  9905 16946  1157   515 
# [1] "weersma_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 44394 46203 38962 36389 33661 40772 31217 29006 24343 28100 27928 26879 19896 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 18281 17192 18578 16707 16496 12999 13759  7987  8184 14441   639   462 
# [1] "xavier_prism_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 45263 47067 39668 37060 34204 41432 31748 29507 24786 28570 28470 27427 20193 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 18597 17539 18975 17071 16731 13297 13939  8113  8372 14668   739   481 
# [1] "xavier_share_gsa"
# 
# 1     2     3     4     5     6     7     8     9    10    11    12    13 
# 45401 47182 39722 37001 34172 41453 31751 29516 24779 28619 28505 27503 20194 
# 14    15    16    17    18    19    20    21    22    23    24    25 
# 18600 17550 18968 17112 16745 13381 13947  8104  8370 14720   729   480 

############

path_gwas=/path/to/ibdgwas/IIBDGC/
studies=(all_hce niddk_old_gwas australia_omniexome gwas1 gwas2 pittsburgh_gsa spain_gsa italy_gsa kiel_austria_sibdcs_gsa netherlands_gsa slovenia_gsa sweden_gsa niddk_broad_gsa niddk_feinstein_gsa basque_gsa lithuania_gsa belgium_louis_gsa belgium_franchimont_gsa belgium_vermeire_gsa prism_nfe_gsa prism_nfe_gwas finland_illugwas german_affy6_old_gwas norway_affy6_old_gwas belgium_inf1_old_gwas belgium_inf2_old_gwas cedars_370k_old_gwas cedars_610k_old_gwas cedars_omni_old_gwas swedish_uc_old_gwas mccauley_gsa ccfa_gsa cedars_gsa bernstein_gsa farkkila_gsa franchimont_gsa franke_gsa helmsley_prism_gsa helmsley_xavier_prism_gsa hyams_protect_gsa lewis_sparc_gsa mccauley_new_gsa mcgovern_gsa moayyedi_imagine_gsa newberry_share_gsa niddk_cho_gsa niddk_duerr_gsa niddk_rioux_gsa niddk_silverberg_gsa palotie_hus_gsa pekow_share_gsa rioux_igenomed_gsa sands_msccr_gsa stampfer_gsa vermeire_gsa weersma_gsa xavier_prism_gsa xavier_share_gsa)

MEM=500

for i in ${studies[@]}
do 
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_21_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_21_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh \
--allow-no-sex \
--exclude ${path_gwas}pre_imputation/QC/${i}/list_chry_variants_toremove \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry"
done

i=swedish_uc_old_gwas
bsub -J"pl2" -M"$MEM" -n 2 -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G your_hpc_group -q normal \
-e ${path_gwas}pre_imputation/QC/${i}/logs/stderr_plink_21_${i} \
-o ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_21_${i} \
"/path/to/software/plink_linux_x86_64_20181202/./plink \
--bfile ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh \
--allow-no-sex \
--make-bed --out ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry"


for i in ${studies[@]}
do echo ${i} && tail -50 ${path_gwas}pre_imputation/QC/${i}/logs/stdout_plink_21_${i} | grep -E "completed"
done

for i in ${studies[@]}
do echo ${i} \
&& wc -l ${path_gwas}pre_imputation/QC/${i}/list_chry_variants_toremove \
&& wc -l ${path_gwas}pre_imputation/QC/${i}/${i}_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
done

###########
# all_hce
# 1504 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/all_hce/list_chry_variants_toremove
# 508214 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/all_hce/all_hce_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# niddk_old_gwas
# 0 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_old_gwas/list_chry_variants_toremove
# 314041 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_old_gwas/niddk_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# australia_omniexome
# 934 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/australia_omniexome/list_chry_variants_toremove
# 787228 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/australia_omniexome/australia_omniexome_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# gwas1
# 0 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas1/list_chry_variants_toremove
# 456448 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas1/gwas1_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# gwas2
# 206 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas2/list_chry_variants_toremove
# 791348 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/gwas2/gwas2_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# pittsburgh_gsa
# 999 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pittsburgh_gsa/list_chry_variants_toremove
# 943978 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pittsburgh_gsa/pittsburgh_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# spain_gsa
# 93 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/spain_gsa/list_chry_variants_toremove
# 585568 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/spain_gsa/spain_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# italy_gsa
# 1698 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/italy_gsa/list_chry_variants_toremove
# 583250 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/italy_gsa/italy_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# kiel_austria_sibdcs_gsa
# 2007 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/kiel_austria_sibdcs_gsa/list_chry_variants_toremove
# 659977 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/kiel_austria_sibdcs_gsa/kiel_austria_sibdcs_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# netherlands_gsa
# 1346 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/netherlands_gsa/list_chry_variants_toremove
# 682886 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/netherlands_gsa/netherlands_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# slovenia_gsa
# 764 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/slovenia_gsa/list_chry_variants_toremove
# 556555 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/slovenia_gsa/slovenia_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# sweden_gsa
# 806 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sweden_gsa/list_chry_variants_toremove
# 574018 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sweden_gsa/sweden_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# niddk_broad_gsa
# 685 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_broad_gsa/list_chry_variants_toremove
# 621064 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_broad_gsa/niddk_broad_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# niddk_feinstein_gsa
# 1002 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_feinstein_gsa/list_chry_variants_toremove
# 643766 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_feinstein_gsa/niddk_feinstein_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# basque_gsa
# 1713 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/basque_gsa/list_chry_variants_toremove
# 576251 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/basque_gsa/basque_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# lithuania_gsa
# 721 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lithuania_gsa/list_chry_variants_toremove
# 578924 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lithuania_gsa/lithuania_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# belgium_louis_gsa
# 850 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_louis_gsa/list_chry_variants_toremove
# 604459 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_louis_gsa/belgium_louis_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# belgium_franchimont_gsa
# 735 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_franchimont_gsa/list_chry_variants_toremove
# 584063 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_franchimont_gsa/belgium_franchimont_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# belgium_vermeire_gsa
# 1154 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_vermeire_gsa/list_chry_variants_toremove
# 673006 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_vermeire_gsa/belgium_vermeire_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# prism_nfe_gsa
# 718 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gsa/list_chry_variants_toremove
# 572622 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gsa/prism_nfe_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# prism_nfe_gwas
# 0 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gwas/list_chry_variants_toremove
# 243968 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/prism_nfe_gwas/prism_nfe_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# finland_illugwas
# 0 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/finland_illugwas/list_chry_variants_toremove
# 239229 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/finland_illugwas/finland_illugwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# german_affy6_old_gwas
# 236 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/german_affy6_old_gwas/list_chry_variants_toremove
# 884962 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/german_affy6_old_gwas/german_affy6_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# norway_affy6_old_gwas
# 183 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/norway_affy6_old_gwas/list_chry_variants_toremove
# 840823 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/norway_affy6_old_gwas/norway_affy6_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# belgium_inf1_old_gwas
# 0 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf1_old_gwas/list_chry_variants_toremove
# 302700 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf1_old_gwas/belgium_inf1_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# belgium_inf2_old_gwas
# 0 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf2_old_gwas/list_chry_variants_toremove
# 290730 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/belgium_inf2_old_gwas/belgium_inf2_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# cedars_370k_old_gwas
# 65 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_370k_old_gwas/list_chry_variants_toremove
# 339880 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_370k_old_gwas/cedars_370k_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# cedars_610k_old_gwas
# 77 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_610k_old_gwas/list_chry_variants_toremove
# 586216 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_610k_old_gwas/cedars_610k_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# cedars_omni_old_gwas
# 870 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_omni_old_gwas/list_chry_variants_toremove
# 720776 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_omni_old_gwas/cedars_omni_old_gwas_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# swedish_uc_old_gwas
# wc: /path/to/ibdgwas/IIBDGC/pre_imputation/QC/swedish_uc_old_gwas/list_chry_variants_toremove: No such file or directory
# mccauley_gsa
# 694 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_gsa/list_chry_variants_toremove
# 570061 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_gsa/mccauley_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# ccfa_gsa
# 751 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/ccfa_gsa/list_chry_variants_toremove
# 610139 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/ccfa_gsa/ccfa_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# cedars_gsa
# 1009 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_gsa/list_chry_variants_toremove
# 667848 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/cedars_gsa/cedars_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# bernstein_gsa
# 630 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/bernstein_gsa/list_chry_variants_toremove
# 564710 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/bernstein_gsa/bernstein_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# farkkila_gsa
# 541 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/farkkila_gsa/list_chry_variants_toremove
# 527381 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/farkkila_gsa/farkkila_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# franchimont_gsa
# 1179 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franchimont_gsa/list_chry_variants_toremove
# 669323 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franchimont_gsa/franchimont_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# franke_gsa
# 684 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franke_gsa/list_chry_variants_toremove
# 570014 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/franke_gsa/franke_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# helmsley_prism_gsa
# 0 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_prism_gsa/list_chry_variants_toremove
# 243969 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_prism_gsa/helmsley_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# helmsley_xavier_prism_gsa
# 0 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_xavier_prism_gsa/list_chry_variants_toremove
# 243943 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/helmsley_xavier_prism_gsa/helmsley_xavier_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# hyams_protect_gsa
# 689 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/hyams_protect_gsa/list_chry_variants_toremove
# 579804 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/hyams_protect_gsa/hyams_protect_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# lewis_sparc_gsa
# 828 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lewis_sparc_gsa/list_chry_variants_toremove
# 618591 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/lewis_sparc_gsa/lewis_sparc_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# mccauley_new_gsa
# 726 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_new_gsa/list_chry_variants_toremove
# 590179 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mccauley_new_gsa/mccauley_new_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# mcgovern_gsa
# 912 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mcgovern_gsa/list_chry_variants_toremove
# 650606 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/mcgovern_gsa/mcgovern_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# moayyedi_imagine_gsa
# 751 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/moayyedi_imagine_gsa/list_chry_variants_toremove
# 614315 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/moayyedi_imagine_gsa/moayyedi_imagine_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# newberry_share_gsa
# 727 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/newberry_share_gsa/list_chry_variants_toremove
# 591538 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/newberry_share_gsa/newberry_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# niddk_cho_gsa
# 782 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_cho_gsa/list_chry_variants_toremove
# 607528 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_cho_gsa/niddk_cho_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# niddk_duerr_gsa
# 780 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_duerr_gsa/list_chry_variants_toremove
# 595047 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_duerr_gsa/niddk_duerr_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# niddk_rioux_gsa
# 747 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_rioux_gsa/list_chry_variants_toremove
# 588507 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_rioux_gsa/niddk_rioux_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# niddk_silverberg_gsa
# 822 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_silverberg_gsa/list_chry_variants_toremove
# 625565 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/niddk_silverberg_gsa/niddk_silverberg_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# palotie_hus_gsa
# 686 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/palotie_hus_gsa/list_chry_variants_toremove
# 564316 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/palotie_hus_gsa/palotie_hus_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# pekow_share_gsa
# 685 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pekow_share_gsa/list_chry_variants_toremove
# 580115 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/pekow_share_gsa/pekow_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# rioux_igenomed_gsa
# 608 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/rioux_igenomed_gsa/list_chry_variants_toremove
# 550785 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/rioux_igenomed_gsa/rioux_igenomed_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# sands_msccr_gsa
# 706 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sands_msccr_gsa/list_chry_variants_toremove
# 603043 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/sands_msccr_gsa/sands_msccr_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# stampfer_gsa
# 1111 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/stampfer_gsa/list_chry_variants_toremove
# 670874 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/stampfer_gsa/stampfer_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# vermeire_gsa
# 1157 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/vermeire_gsa/list_chry_variants_toremove
# 673577 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/vermeire_gsa/vermeire_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# weersma_gsa
# 639 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/weersma_gsa/list_chry_variants_toremove
# 572836 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/weersma_gsa/weersma_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# xavier_prism_gsa
# 739 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_prism_gsa/list_chry_variants_toremove
# 583178 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_prism_gsa/xavier_prism_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim
# xavier_share_gsa
# 729 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_share_gsa/list_chry_variants_toremove
# 583775 /path/to/ibdgwas/IIBDGC/pre_imputation/QC/xavier_share_gsa/xavier_share_gsa_hg19_noind_posstr_nodup_flip_sexcheck_missinghh_nochry.bim

###########