# Author: Rui Zhang
# Institution: Analytic and Translational Genetics Unit, Department of Medicine,
#              Massachusetts General Hospital, Boston, MA, USA
#
library(susieR)
library(data.table)

args <- commandArgs(T)
Type <- as.character(args[1])
chr <- as.character(args[2])
start <- as.character(args[3])
end <- as.character(args[4])

proj_path <- "/stanley/huang_lab/home/rzhang/11_IBD_FM/"

dtSumS <- fread(paste0(proj_path,"00_data_Aug25/",Type,"_sumstats/",chr,"_",start,"_",end,".sumstats.txt"))
dtLD <- as.matrix(fread(paste0(proj_path,"00_data_Aug25/",Type,"_ld/",chr,"_",start,"_",end,".dose.subset_list_variants_tokeep_in_analysis_",Type,"_all_arrays_combined.ld.gz")))
dtLD[is.na(dtLD)] <- 0

if (Type=="cd"){
	N = 65266
}else if (Type=="uc"){
	N = 58246
}else if (Type=="ibd"){
	N = 91113
}
LL = 10
finish = 0

result <- susie_rss(z=dtSumS$BETA/dtSumS$SE, R=dtLD, n=N, L=LL, estimate_residual_variance=TRUE)

if (result$converged == TRUE) {
	if (length(result$sets$coverage) >= 5) {
		LL = 20
		result <- susie_rss(z=dtSumS$BETA/dtSumS$SE, R=dtLD, n=N, L=LL, estimate_residual_variance=TRUE)

		if (result$converged == FALSE) {
			while (finish == 0) {
				LL = LL-1
				result <- susie_rss(z=dtSumS$BETA/dtSumS$SE, R=dtLD, n=N, L=LL, estimate_residual_variance=TRUE)

				if (result$converged == TRUE) {
					finish = 1
				}	
			}
		}
	}
} else {
	while (finish == 0) {
		LL = LL-1
		result <- susie_rss(z=dtSumS$BETA/dtSumS$SE, R=dtLD, n=N, L=LL, estimate_residual_variance=TRUE)

		if (result$converged == TRUE) {
			finish = 1
		}
	}
}

var1 = names(result)[1:9]
for (var in var1) {
	index = which(var==var1)
	write.table(result[[var]], file=paste0(proj_path,"11_FM/02_results/",Type,"/0",index,"_",var,"/",chr,"_",start,"_",end,"_",var,".tsv"),row.names=FALSE,col.names=FALSE,sep="\t")	
}

var2 = names(result)[10:16]
if (var2[3] == "elbo") {
	var2[3:5] = c("converged", "elbo", "niter")
}
for (var in var2) {
	index = which(var==var2)+9
	write.table(result[[var]], file=paste0(proj_path,"11_FM/02_results/",Type,"/",index,"_",var,"/",chr,"_",start,"_",end,"_",var,".tsv"),row.names=FALSE,col.names=FALSE,sep="\t")
}

capture.output(print(result$sets), file=paste0(proj_path,"11_FM/02_results/",Type,"/17_sets/",chr,"_",start,"_",end,"_sets.txt"))

nVari = nrow(dtSumS)

CreIndex <- rep(0, nVari)
CreSetL <- names(result$sets$cs)
for (CreSet in CreSetL) {
	for (indexL in result$sets$cs[CreSet]) {
		for (index in indexL) {
			CreIndex[index] = 1
		}
	}
}

CS_N_purity <- rep(0, nVari)
CS_alpha <- rep(0, nVari)
purity_min <- rep(0, nVari)
purity_mean <- rep(0, nVari)
purity_median <- rep(0, nVari)
CS_coverage <- rep(0, nVari)
label <- 0
for (CreSet in CreSetL) {
	label <- label+1
	nL <- result$sets$cs_index[label]
	for (indexL in result$sets$cs[CreSet]) {
		for (index in indexL) {
			CS_N_purity[index] = label
			CS_alpha[index] = result$alpha[nL, index]
			purity_min[index] = result$sets$purity[label,1]
			purity_mean[index] = result$sets$purity[label,2]
			purity_median[index] = result$sets$purity[label,3]
			CS_coverage[index] = result$sets$coverage[label]
		}
	}
}

CS_N_cov <- rep(0, nVari)
if (length(result$sets$cs) != 0) {
	label <- 0
	sorted_cov <- names(result$sets$cs)[order(result$sets$coverage, decreasing=TRUE)]
	for (CreSet in sorted_cov) {
		label <- label+1
		for (indexL in result$sets$cs[CreSet]) {
			for (index in indexL) {
				CS_N_cov[index] = label
			}
		}
	}
}

PIP <- cbind(dtSumS[,c(1:8)], result$pip, CS_alpha, CreIndex, CS_N_purity, purity_min, purity_mean, purity_median, CS_N_cov, CS_coverage, dtSumS[,c(9:15,17:21,23)])
PIP[dbsnp154=="", dbsnp154:=MarkerName]
names(PIP)[8] <- "P_value"
names(PIP)[9] <- "pip"
names(PIP)[10] <- "CS_alpha"
names(PIP)[11] <- "CreIndex"
names(PIP)[12] <- "CS_N_purity"
names(PIP)[13] <- "purity_min"
names(PIP)[14] <- "purity_mean"
names(PIP)[15] <- "purity_median"
names(PIP)[16] <- "CS_N_cov"
names(PIP)[17] <- "CS_coverage"

fwrite(PIP, file=paste0(proj_path,"11_FM/02_results/",Type,"/18_pip/",chr,"_",start,"_",end,"_pip.tsv"), sep="\t", quote=FALSE)

