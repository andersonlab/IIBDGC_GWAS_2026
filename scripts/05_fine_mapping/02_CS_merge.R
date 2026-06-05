# Author: Rui Zhang
# Institution: Analytic and Translational Genetics Unit, Department of Medicine,
#              Massachusetts General Hospital, Boston, MA, USA
#
library(data.table)
library(readxl)
library(dplyr)
library(tidyr)
library(Matrix)
library(writexl)

data_CreSet <- data.frame(
	chr = character(),
	region_start = integer(),
	region_end = integer(),
	region_len = integer(),
	region_label = character(),
	trait = character(),
	trait_CS = character(),
	CS_index = integer(),
	CS_start = integer(),
	CS_end = integer(),
	CS_len = integer(),
	nSNP_CS = integer(),
	purity_min = numeric(),
	variant_lead = character(),
	position_lead = integer(),
	MarkerName_lead = character(),
	alpha_lead = numeric(),
	BETA_lead = numeric(),
	SE_lead = numeric(),
	p_value_lead = numeric(),
	Func_lead = character(),
	Gene_lead = character(),
	GeneDetail_lead = character(),
	Exonic_Func_lead = character(),
	AAChange_lead = character(),
	Canonical_lead = character(),
	AAChange_lead_VEP = character(),
	N_lead = integer(),
	A2freq_lead = numeric(),
	MAF_lead = numeric(),
	N_cases_lead = integer(),
	A2freq_cases_lead = numeric(),
	N_controls_lead = numeric(),
	A2freq_controls_lead = numeric(),
	rate_total_sample_lead = numeric(),
	N_effective_lead = integer(),
	rate_total_effective_sample_lead = numeric(),
	MAF_mean = numeric(),
	rate_total_sample_mean = numeric(),
	rate_total_effective_sample_mean = numeric(),
	all_CS_variant = character())

data_Vari <- data.frame(
	chr = character(),
	region_mb = character(),
	region_label = character(),
	trait = character(),
	trait_CS = character(),
	trait_variant = character(),
	variant = character(),
	position = integer(),
	MarkerName = character(),
	CS_index = integer(),
	nSNP_CS = integer(),
	purity_min = numeric(),
	alpha_lead = numeric(),
	BETA_lead = numeric(),
	SE_lead = numeric(),
	p_value_lead = numeric(),
	MAF_lead = numeric(),
	rate_total_sample_lead = numeric(),
	rate_total_effective_sample_lead = numeric(),
	MAF_mean = numeric(),
	rate_total_sample_mean = numeric(),
	rate_total_effective_sample_mean = numeric(),
	lead = logical(),
	alpha = numeric(),
	BETA = numeric(),
	SE = numeric(),
	p_value = numeric(),
	novel = logical(),
	pre_Pro = numeric(),
	pre_p_multi = numeric(),
	pre_tier2 = character(),
	pre_trait = character(),
	pre_trait_reassigned = character(),
	Func = character(),
	Gene = character(),
	GeneDetail = character(),
	Exonic_Func = character(),
	AAChange = character(),
	Canonical = character(),
	AAChange_VEP = character(),
	N = integer(),
	A2freq = numeric(),
	MAF = numeric(),
	N_cases = integer(),
	A2freq_cases = numeric(),
	N_controls = integer(),
	A2freq_controls = numeric(),
	rate_total_sample = numeric(),
	N_effective = integer(),
	rate_total_effective_sample = numeric())

proj_path <- "/stanley/huang_lab/home/rzhang/11_IBD_FM/"

jaccard <- function(input_m) {
	#Non zero values
	A <- tcrossprod(t(input_m))
	#Indices for non-zero values, pair of CS need to be compared
	im <- which(as.matrix(A)>0, arr.ind=TRUE)
	#Weighted jaccard similarity
	wj <- apply(im, 1, function(x) {
		x2 <- matrixStats::rowRanges(as.matrix(input_m[, x, drop=F]), cols=c(1,2))
		s <- matrixStats::colSums2(x2)
		s[1]/s[2]
	})

	J <- sparseMatrix(
		i = im[,1],
		j = im[,2],
		x = wj,
		dims = dim(A),
		dimnames = dimnames(A))

	return(J)
}

find_overlaps <- function(CS_HC) {
	if (length(unique(CS_HC$id)) == 1) {
		CS_HC <- dplyr::mutate(CS_HC, csm_id = 1) %>%
				select(id, trait, cs_id, pip, variant, csm_id) 
		return(CS_HC)
	}
		
	#Convert to sparse matrix, N_unique_variant * N_CS_total
	CS_HC_mat <- CS_HC %>%
		tidytext::cast_sparse(variant, id, pip) %>%
		as(., "dgCMatrix")

	#Compute weighted Jaccard similarity matrix
	CS_HC_jsm <- jaccard(CS_HC_mat)
	#Get dissimilarity matrix
	CS_HC_dist <- as.dist(1-CS_HC_jsm)

	#Hierarchical clustering using Complete Linkage
	CS_HC_hc1 <- hclust(CS_HC_dist, method="complete")

	#Give cutoff, assign IDs to 95% CSs
	CS_HC_hc1 <- cutree(CS_HC_hc1, h=0.9) %>%
		data.frame() %>%
		tibble::rownames_to_column() %>%
		tibble::as_tibble() %>%
		dplyr::rename("id"="rowname", "csm_id"=".")

	CS_HC <- CS_HC %>%
		merge(., CS_HC_hc1, by="id", sort=FALSE) %>%
		tibble::as_tibble() %>%
		tidyr::separate(id, into=c("trait", "cs_id"), sep=";", remove=F) %>%
		select(id, trait, cs_id, pip, variant, csm_id)

	return(CS_HC)
}

CS_merge <- function(CS_HC, high_LD_list) {
	CSm <- data.frame(new_trait = character(),
			CS_trait = character(),
			variant_trait = character(),
			pip = numeric(),
			variant = character(),
			csm_id = integer())

	# further merge if the LD between two lead variants is larger than 0.8	
	lead_vari_df = data.frame(vari=character(), index=integer())
	further_merge_id = c()
	for (CS_N in unique(CS_HC$csm_id)) {
		temp <- CS_HC[CS_HC$csm_id==CS_N,]
		N_class = length(unique(temp$id))
		temp <- temp %>% 
			group_by(variant) %>%
			mutate(mean_pip = sum(pip)/N_class)	

		lead_vari <- temp[which.max(temp$mean_pip),]$variant
		temp_data <- data.frame(vari = lead_vari, index=CS_N)
		lead_vari_df <- rbind(lead_vari_df, temp_data)
	}
	n_lead = nrow(lead_vari_df)
	for (lead_index1 in 1:(n_lead-1)) {
		for (lead_index2 in (lead_index1+1):n_lead) {
			vari_comb <- paste(lead_vari_df[lead_index1,]$vari, lead_vari_df[lead_index2,]$vari, sep="_")
			if (vari_comb %in% high_LD_list) {
				id_comb <- paste(lead_index1, lead_index2, sep="_")
				further_merge_id = c(further_merge_id, id_comb)
			}
		}
	}

	CS_HC <- CS_HC %>% mutate(id_change = 0)

	for (id_comb in further_merge_id) {
		id_comb <- strsplit(id_comb, "_")[[1]]
		id1 <- as.integer(id_comb[1])
		id2 <- as.integer(id_comb[2])

		new_id_change = CS_HC[CS_HC$csm_id==id1,]$id_change[1]
		
		CS_HC <- CS_HC %>%
			mutate(id_change = if_else(csm_id==id2, new_id_change, id_change)) %>%
			mutate(csm_id = if_else(csm_id==id2, id1, csm_id))

		if (id2 < n_lead) {
			for (id_behind in (id2+1):n_lead) {
				CS_HC <- CS_HC %>%
					mutate(id_change = if_else(csm_id==id_behind, id_change+1, id_change))
			}
		}
	}		

	CS_HC <- CS_HC %>%
		mutate(csm_id = csm_id-id_change) %>%
		select(-id_change)

	for (CS_N in unique(CS_HC$csm_id)) {
		temp <- CS_HC[CS_HC$csm_id==CS_N,]
		phe <- unique(temp$trait)
		N_class <- length(unique(temp$id))

		if (N_class == 1) {
			new_trait <- rep(phe[1], nrow(temp))
			CSm_part <- cbind(new_trait, new_trait, new_trait, temp[,4:6])
		} else {
			CS_trait <- paste(phe, collapse="_")

			CSm_part <- data.frame()
			for (vari in unique(temp$variant)) {
				temp_vari <- temp[temp$variant == vari,]
				variant_trait <- paste(unique(temp_vari$trait), collapse="_")
				mean_PIP <- sum(temp_vari$pip)/N_class

				final_trait = "ibd"				
				if (CS_trait == "ibd_cd") {
					final_trait = "cd"	
				} else if (CS_trait == "ibd_uc") {
					final_trait = "uc"
				}

				CSm_part_row <- data.frame(final_trait, CS_trait, variant_trait, mean_PIP, vari, CS_N)
				CSm_part <- rbind(CSm_part, CSm_part_row)
			}
		}
		names(CSm_part) <- names(CSm)
		CSm_part <- CSm_part %>% arrange(desc(pip))
		CSm <- rbind(CSm, CSm_part)
	}

	return(CSm)
}

#2017 nature paper
nature <- fread(paste0(proj_path,"02_preResult/17_Huang_results/all_hg38.tsv"))
nature$chr <- as.character(nature$chr)
nature$position_b38 <- as.character(nature$position_b38)
nature$MarkerName <- apply(nature[,c(2,16,17,18)], 1, function(x) paste0("chr",paste(x,collapse=":")))
#MarkerName, pro, p_multi, tier2, trait...8, trait_reassigned
nature <- nature[,c(36,26,14,6,8,9)]
nature$tier2[nature$tier2 == "No"] <- "three_method"

fsvari <- nature[nature$MarkerName=="chr16:50729870:I:D",]
fsvari$MarkerName = "chr16:50729867:G:GC"
nature <- rbind(nature, fsvari)

#regions <- read.table(paste0(proj_path,"00_data_Aug25/01_region/sig_class/IL23R_region.txt"), header=F)
regions <- read.table(paste0(proj_path,"00_data_Aug25/01_region/sig_class/regions_significance.txt"), header=F)

regions <- regions %>%
	mutate(region_name = paste0(substr(V1,4,nchar(V1)),"_",V2,"_",V3))

for (i in 1:nrow(regions)) {
	phe <- c("ibd","cd","uc")

	Chr = regions[i,1]
	start = as.numeric(regions[i,2])
	end = as.numeric(regions[i,3])
	region_name = regions[i,4]
	
	print(region_name)

	input_all <- data.frame()

	for (Type in phe) {
		pips <- fread(paste0(proj_path,"11_FM/02_results/",Type,"/18_pip/",Chr,"_",start,"_",end,"_pip.tsv"))
		
		dtSumS <- fread(paste0(proj_path,"00_data_Aug25/",Type,"_sumstats/",Chr,"_",start,"_",end,".sumstats.txt"))
		dtSumS <- dtSumS %>%
			select(MarkerName, MAF, Neff, rate_Neff)

		#add MAF
		pips <- pips %>%
			left_join(dtSumS %>% select(MarkerName, MAF), by="MarkerName")	

		#add annotation info
		annotation <- read.table(paste0(proj_path,"01_anno_AV/",Type,"/final_result/",Type,"_",Chr,"_",start,"_",end,".anno.txt"), header =TRUE)
		pips <- cbind(pips, annotation[, c(2:8)])

		#add effective sample size
		pips <- pips %>%
			left_join(dtSumS %>% select(MarkerName,Neff,rate_Neff), by="MarkerName") %>%
			rename("N_effective"="Neff","rate_total_effective_sample"="rate_Neff")

		#filter CS based on MAF & rate_total_effective_sample of lead variant
		#filter CS based on marginal p-value of lead variant
		N_CS <- max(pips$CS_N_purity)
		if (N_CS != 0) {
			list1 <- unique(pips$CS_N_purity)
			for (n_CS in 1:N_CS) {
				pips_CS <- pips[pips$CS_N_purity == n_CS,]
				lead_vari <- pips_CS[pips_CS$CS_alpha == max(pips_CS$CS_alpha),]
				if ((lead_vari[1,]$MAF<0.01) & (lead_vari[1,]$rate_total_effective_sample<0.95)) {
					pips <- pips[pips$CS_N_purity != n_CS,]
				}
				if (lead_vari[1,]$P_value >= 5E-8) {
					pips <- pips[pips$CS_N_purity != n_CS,]
				}
			}
			list2 <- unique(pips$CS_N_purity)
			diff <- setdiff(list1, list2)
			pips$n_diff <- rep(0, nrow(pips))
			for (CS in list2) {
				n = sum(diff < CS)
				pips[pips$CS_N_purity == CS,]$n_diff <- n
			}
			pips$CS_N_purity <- pips$CS_N_purity-pips$n_diff
			pips <- pips[,1:(ncol(pips)-1)]
		}

		Trait <- rep(Type, nrow(pips))
		pips <- cbind(Trait, pips)

		input_all <- rbind(input_all, pips)
	}

	input_all_CS <- input_all[input_all$CS_N_purity != 0, ]

	region_no_signal = "FALSE"
	if (nrow(input_all_CS) == 0) {
		region_no_signal = "TRUE"
		allCS_MarkerName = c()
	}

	if (nrow(input_all_CS) != 0) {	
		#merge LD information of lead variant
		high_LD_list <- c()	
		for (Type in phe) {
			LD_path <- paste0(proj_path,"11_FM/04_CS_merge/01_lead_LD/02_LD_result/",Type,"/",Type,"_",Chr,"_",start,"_",end,".ld")
			if (!(file.exists(LD_path))) {
				next
			}	
			LD <- read.table(LD_path, header = TRUE)
			var_list <- apply(LD, 1, function(row) paste(row[3], row[6], sep="_"))
			high_LD_list <- c(high_LD_list, var_list)
		}
		high_LD_list <- unique(high_LD_list)

		#CS alignment across three traits
		input <- input_all_CS[,c(1,13,11,2)]
		names(input) <- c("trait", "cs_id", "pip", "variant")

		CSm <- dplyr::mutate(input, id=paste(trait, cs_id, sep=";")) %>% 
			find_overlaps() %>% 
			CS_merge(., high_LD_list)
 
		#extract other info from pips & nature_result based on CSm
		#generate data_list_of_CreSet & data_list_of_variant
		N_cs = as.integer(max(CSm$csm_id))
		allCS_MarkerName = CSm$variant

		region_CreSet <- data.frame()
		region_Vari <- data.frame()

		for (Index in 1:N_cs) {
			CSm_part <- CSm[CSm$csm_id == Index, ]
			CSTrait <- CSm_part$CS_trait[1]
			newTrait <- CSm_part$new_trait[1]

			pips_part <- data.frame()
			for (vari in c(CSm_part$variant)) {
				CSm_part_row <- CSm_part[CSm_part$variant==vari, ]
				TraitList <- strsplit(CSm_part_row$variant_trait, "_")[[1]]
				if (length(TraitList) == 2) {
					VariTrait = TraitList[2]
				} else {
					VariTrait = TraitList[1]
				}
	
				pips_row <- input_all_CS[(input_all_CS$Trait==VariTrait) & (input_all_CS$MarkerName==vari),]
				pips_row$CS_alpha <- CSm_part_row$pip
				pips_row$variant_trait <- CSm_part_row$variant_trait
				pips_part <- rbind(pips_part, pips_row)
			}

			cs_start = min(pips_part$Position_b38)
			cs_end = max(pips_part$Position_b38)
			n_cs = nrow(pips_part)

			#list of credible set
			new_row_CreSet_part <- data.frame(
				chr = Chr,
				region_start = start,
				region_end = end,
				region_len = end-start+1,
				region_label = region_name,
				trait = toupper(newTrait),
				trait_CS <- toupper(CSTrait),
				CS_index = Index,
				CS_start = cs_start,
				CS_end = cs_end,
				CS_len = cs_end-cs_start+1,
				nSNP_CS = n_cs)

			MAF_mean = mean(pips_part$MAF)
			rate_mean = mean(pips_part$rate_total_sample)
			rate_effective_mean = mean(pips_part$rate_total_effective_sample)

			pips_part$purity_min <- min(pips_part$purity_min)

			new_row_CreSet <- cbind(new_row_CreSet_part,
				pips_part[1,c(14,3,4,2,11,7:9,33:39,25:26,32,27:31,40,41)],
				MAF_mean, rate_mean, rate_effective_mean,
				all_CS_variant = paste(pips_part$dbsnp154, collapse=", "))

			names(new_row_CreSet) <- names(data_CreSet)
			region_CreSet <- rbind(region_CreSet, new_row_CreSet)

			#list of variants
			chr_list = rep(Chr, n_cs)
			re_mb = rep(paste(round(start/1000000,2),round(end/1000000,2),sep="-"), n_cs)
			region_name_list = rep(region_name, n_cs)			

			trait_list = rep(toupper(newTrait), n_cs)
			trait_CS_list = rep(toupper(CSTrait), n_cs)
			CS_index = rep(Index, n_cs)
			lead_list = c("TRUE", rep("FALSE", n_cs-1))
			novel_list = !(pips_part$MarkerName %in% nature$MarkerName)

			nSNP_CS_list = rep(n_cs, n_cs)
			alpha_lead_list = rep(pips_part$CS_alpha[1], n_cs)
			BETA_lead_list = rep(pips_part$BETA[1], n_cs)
			SE_lead_list = rep(pips_part$SE[1], n_cs)
			p_value_lead_list = rep(pips_part$P_value[1], n_cs)
			MAF_lead_list = rep(pips_part$MAF[1], n_cs)	
			rate_lead_list = rep(pips_part$rate_total_sample[1], n_cs)
			rate_effective_lead_list = rep(pips_part$rate_total_effective_sample[1], n_cs)
			MAF_mean_list = rep(MAF_mean, n_cs)
			rate_mean_list = rep(rate_mean, n_cs)
			rate_effective_mean_list = rep(rate_effective_mean, n_cs)
		
			pips_part$variant_trait <- toupper(pips_part$variant_trait)	


			#extract information from 2017 nature paper
			pre_data <- data.frame(
				pre_Pro = rep(0, n_cs),
				pre_p_multi = rep(0, n_cs),
				pre_tier2 = rep(".", n_cs),
				pre_trait = rep(".", n_cs),
				pre_trait_reassigned = rep(".", n_cs))

			overlap_index = which(pips_part$MarkerName %in% nature$MarkerName)
			for (overlap in overlap_index) {
				pre_index = which(nature$MarkerName == pips_part$MarkerName[overlap])
				pre_data[overlap,1:5] = nature[pre_index,2:6]
			}

			CS_Vari <- cbind(chr_list, re_mb, region_name_list, trait_list, trait_CS_list, pips_part[,c(42,3,4,2)], CS_index, nSNP_CS_list, pips_part[,14], alpha_lead_list, BETA_lead_list, SE_lead_list, p_value_lead_list, MAF_lead_list, rate_lead_list, rate_effective_lead_list, MAF_mean_list, rate_mean_list, rate_effective_mean_list, lead_list, pips_part[,c(11,7:9)], novel_list, pre_data, pips_part[,c(33:39,25:26,32,27:31,40,41)])
			names(CS_Vari) <- names(data_Vari)
			region_Vari <- rbind(region_Vari, CS_Vari)
		}

		names(region_CreSet) <- names(data_CreSet)
		names(region_Vari) <- names(data_Vari)
 
		region_CreSet <- region_CreSet %>%
			group_by(CS_index) %>%
			arrange(desc(purity_min), desc(alpha_lead)) %>%
			ungroup() %>%
			mutate(CS_index_new = dense_rank(factor(CS_index, levels=unique(CS_index)))) %>%
			mutate(CS_index = CS_index_new) %>%
			select(-CS_index_new)

		region_Vari <- region_Vari %>%
			group_by(CS_index) %>%
			mutate(first_alpha = first(alpha)) %>%
			arrange(desc(purity_min), desc(first_alpha)) %>%
			ungroup() %>%
			mutate(CS_index_new = dense_rank(factor(CS_index, levels=unique(CS_index)))) %>%
			mutate(CS_index = CS_index_new) %>%
			select(-CS_index_new, -first_alpha)

		data_CreSet <- rbind(data_CreSet, region_CreSet)
		data_Vari <- rbind(data_Vari, region_Vari)
	}
}

write_xlsx(data_CreSet, "merge_results/list_of_credible_sets.xlsx")
write_xlsx(data_Vari, "merge_results/list_of_variants.xlsx")

