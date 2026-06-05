# IIBDGC GWAS Pipeline

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

Scripts for quality control, imputation, association analysis, meta-analysis, and functional annotation of the International IBD Genetics Consortium (IIBDGC) genome-wide association study.

## Overview

This repository contains the full analytical pipeline for a large-scale, multi-cohort GWAS of inflammatory bowel disease (IBD) and its subtypes — Crohn's disease (CD) and ulcerative colitis (UC) — across multiple ancestry groups. The pipeline covers everything from per-cohort genotype QC through to functional annotation, and pathway analysis of GWAS signals.

### Study design
- **Phenotypes:** IBD, CD, UC
- **Ancestry groups:** European (all), European non-Jewish, European Jewish, South Asian, East Asian, multi-ancestry
- **Tier 1 cohorts:** 50+ cohorts across multiple genotyping arrays (Illumina 370K, Affymetrix 6.0, HumanOmniExpress, GSA, and others)
- **Tier 2 cohorts:** FinnGen, deCODE (Iceland), Danish cohort, East Asian IIBDGC

---

## Repository structure

```
scripts/
├── 01_qc_per_cohort_step_1_to_step_3/       # Cohort-specific QC (Steps 1–3)
├── 02_qc_imputation_analysis_step_4_onwards/ # Shared QC, analysis & meta-analysis (Steps 4–64)
├── 03_tier_2_datasets/                       # Reformatting external summary statistics
├── 04_sex_specific_analyses/                 # Sex-specific GWAS and meta-analysis
└── 05_other/                                 # Supplementary analyses and visualizations
```

---

## Pipeline

### Phase 1 — Cohort-specific QC (Steps 1–3)

Scripts in `01_qc_per_cohort_step_1_to_step_3/` are named `[COHORT]_qc_imputation_README.R` and perform:

- Sample QC: sex checks, call rates, heterozygosity
- Variant QC: allele frequency filtering, platform-specific exclusions
- Preparation of files for imputation

### Phase 2 — Shared QC, analysis, and meta-analysis (Steps 4–64)

Scripts in `02_qc_imputation_analysis_step_4_onwards/` are applied uniformly across all Tier 1 cohorts:

| Steps | Description |
|-------|-------------|
| 4–10  | Data harmonisation: strand flips, variant renaming, allele frequency comparison against 1000GP, sex inference, Y-chromosome removal |
| 11–15 | Sample and variant QC: call rate comparisons, duplicate removal, PCA-based ancestry assignment, HWE filtering, heterozygosity |
| 16–23 | Build 38 conversion, TopMed variant checks, AT/CG allele handling, ancestry splitting, imputation input preparation |
| 24–29 | Post-imputation QC, relatedness inference, per-array PCA |
| 30–31 | Phenotype file creation; association testing with REGENIE (mixed-model, per ancestry and per array) |
| 32–35 | Meta-analysis: EUR with METAL; SAS with METAL; multi-ancestry with METAL and MR-MEGA; integration of Tier 2 datasets |
| 36–43 | Fine-mapping: COJO conditional analysis, independent signal definition, LD matrix generation, closest gene annotation |
| 44–45 | WES results integration; GSMR direction-of-effect analysis |
| 46–50 | Enrichment analyses: PolyFun, heritability (LDSC), genetic correlations (IBD/CD/UC), latent causal model, local genetic correlations |
| 51–54 | eQTL and pQTL colocalization; effector gene consensus; direction of effect validation |
| 55–59 | Functional annotation: receptor-ligand pairs, protein complexes, monogenic disease genes, OTAR drug targets, pathway analysis |
| 60–64 | Final signal table compilation; comparative analyses; regional Manhattan plots |

### Tier 2 dataset integration

Scripts in `03_tier_2_datasets/` reformat external summary statistics (FinnGen, deCODE, Danish, East Asian) into the format required by METAL for meta-analysis.

### Sex-specific analyses

Scripts in `04_sex_specific_analyses/` run sex-stratified GWAS with REGENIE and perform sex-specific meta-analysis using GWAMA. These analyses were developed and implemented by Talin Haritunians (F. Widjaja Inflammatory Bowel Disease Institute).

### Supplementary analyses

Scripts in `05_other/` include:
- Colocalization weight matrix generation and result organisation
- Comparison of results against prior IIBDGC and external meta-analyses (GWAS3, Liu, Huang, de Lange)
- Visualisation: Manhattan plots, locus zoom plots, colocalization heatmaps, effect size comparisons
- eQTL data reformatting for multiple sources (eQTL Catalogue, GTEx, IBDverse, macromap, hu_2021)
- Burden heritability enrichment analyses

---

## Tools and software

| Tool | Purpose |
|------|---------|
| PLINK 1.9 / 2.0 | Genotype QC and data management |
| REGENIE v3.2.5 | Mixed-model GWAS association testing |
| METAL | Fixed/random-effects meta-analysis |
| MR-MEGA | Multi-ancestry meta-analysis |
| GCTA-COJO | Conditional and joint analysis for fine-mapping |
| LDSC | SNP heritability and genetic correlation |
| PolyFun | Functionally-informed enrichment and fine-mapping priors |
| coloc | Bayesian colocalization with eQTL/pQTL data |
| GSMR | Direction-of-effect testing |
| BCFtools / SAMtools | VCF handling and variant annotation |

**R packages** (key): `data.table`, `ggplot2`, `coloc`, `qusage` (pathway analysis), and others listed within individual scripts.

---

## Functional annotation data sources

**eQTL:**
- eQTL Catalogue (multiple tissues and cell types)
- GTEx
- hu_2021 (Human Cell Atlas single-cell eQTL)
- macromap (macrophage single-cell eQTL)
- IBDverse (IBD-specific single-cell eQTL)

**pQTL:**
- deCODE
- SPARC

**Gene and pathway annotation:**
- CORUM (protein complexes)
- OTAR (drug target evidence)
- Baderlab gene sets (pathway analysis)
- Ligand-receptor pair databases

---

## Reference panels and key QC thresholds

- **Reference panels:** 1000 Genomes Project, UK Biobank, TopMed
- **Sample call rate:** >80% (samples below this threshold removed)
- **Variant call rate:** >95–98%
- **Imputation INFO score:** >0.4 (minimum); higher thresholds applied for specific analyses
- **MAF:** >0.001
- **HWE p-value:** standard thresholds applied in controls

---

## Software environment

Two Singularity containers are available on Zenodo:

| Container | Used for | Zenodo |
|-----------|----------|--------|
| `iibdgc_postprocess_10_singularity.sif` | Main pipeline (QC, association, meta-analysis, annotation) | https://zenodo.org/records/20557380 |
| `polyfun_2_singularity.sif` | Enrichment and heritability analyses (step 46–48) | https://zenodo.org/records/20557380 |

To run:
```bash
singularity shell iibdgc_postprocess_10_singularity.sif
```

---

## Citation

A preprint is available on medRxiv:

> Fachal L, et al. (2026). Resolving inflammatory bowel disease risk variants to genes and cell types. medRxiv. https://www.medrxiv.org/content/10.64898/2026.05.13.26352926v2
