# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
## /software/R-4.3.1/bin/R

library("rjson")

rm(list=ls())

path<-"/path/to/ibdgwas/IIBDGC/"

args<-commandArgs()

cohorts<-args[6]
print(cohorts)

ancestry<-c("eur","noneur")

for (i in 1:length(cohorts)) {
  
  print(cohorts[i])

  for (j in 1:length(ancestry)) {

    file_json<-paste(path,"imputed/",cohorts[i],"/2022/logs/status_",cohorts[i],"_",ancestry[j],"_qc.json",sep="")

    if(file.exists(file_json)) {

      print(ancestry[j])
      dat<-fromJSON(file=file_json)

      commnad_download<-paste("curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/",dat$outputParams[[3]]$id,"/",
                              dat$outputParams[[3]]$hash," > ",path,"imputed/",cohorts[i],"/2022/logs/donwload_imputed_results_",cohorts[i],"_",ancestry[j],sep="")
      system(commnad_download)
    }
  }
}
q("no")
d