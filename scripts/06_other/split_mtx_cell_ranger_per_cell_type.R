# Author: Laura Fachal
# Institution: Wellcome Sanger Institute
# ORCID: https://orcid.org/0000-0002-7256-9752
#
# MEM=8000
# bsub -Is -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM]" -G your_hpc_group R \

library(data.table)

args<-commandArgs(trailingOnly=TRUE)
tissue<-args[1]
# tissue<-"non_multiome_duodenum_epithelial"

path<-"/path/to/ibdgwas/IIBDGC/"

mtx<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/annotation/hickey_2023/",tissue,"_hickney_peak_matrix_final_matrix.csv",sep=""),head=F)
colnames(mtx)<-c("cell_barcode","chr","start","end","cut_sites")

# 10x cellular barcode
# Peak chromosome
# Peak start position
# Peak end position
# # of cut sites within peak

dict<-fread(paste(path,"post_imputation/2022/analysis/metaanalysis/annotation/hickey_2023/scATAC_cell_types_",tissue,".tsv",sep=""),head=T)
# print(table(dict$CellType))
dict$CellType<-gsub(" ","_",dict$CellType)
dict$CellType<-gsub("/","_",dict$CellType)
dict$CellType<-gsub("\\+","",dict$CellType)
dict$CellType<-gsub("-","_",dict$CellType)

# combine with dictionary, and split by cell type:
mtx<-merge(mtx,dict,by.x="cell_barcode",by.y="Cell")

ids<-names(table(mtx$CellType))
print(ids)
class(mtx$start)
class(mtx$end)

# save one file per cell type
for (i in 1:length(ids)) {

    print(ids[i])
    tmp<-mtx[which(mtx$CellType==ids[i]),]
    tmp$chr_ed<-as.numeric(gsub("chr","",tmp$chr))
    tmp$chr_ed[which(is.na(tmp$chr_ed))]<-23

    tmp<-tmp[order(tmp$chr,tmp$start,decreasing=F),c("chr","start","end")]
    tmp<-tmp[!duplicated(tmp),]
    print(dim(tmp))

    fwrite(tmp[,c("chr","start","end")],paste(path,"post_imputation/2022/analysis/metaanalysis/annotation/hickey_2023/",tissue,"_",ids[i],".bed",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

}

q("no")



