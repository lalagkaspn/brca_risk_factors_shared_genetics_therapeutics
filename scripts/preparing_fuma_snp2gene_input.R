
### Creating input files for FUMA SNP2GENE to get the protein-coding genes located within the LOGODetect positively correlated regions (for each brca-predisposing disease pair) ###

library(dplyr)
library(data.table)

## -- load GWAS summary statistics used as LOGODetect input -- ##
bc = fread("preprocessed_data/icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics_LOGODetect.txt")
depression = fread("preprocessed_data/PGC_UKB_depression_genome-wide_LOGODetect.txt")
hdl = fread("preprocessed_data/HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results_LOGODetect.txt")
ldl = fread("preprocessed_data/LDL_GLGC_EUR_2021/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results_LOGODetect.txt")
prostate = fread("preprocessed_data/PROSTATE_CANCER_PRACTICAL_EUR_2018_LOGODetect.txt")
schizophrenia = fread("preprocessed_data/PGC3_SCZ_wave3.european.autosome.public.v3.vcf_LOGODetect.txt")
t2dm_diamante = fread("preprocessed_data/DIAMANTE-EUR.sumstat_LOGODetect.txt")

## -- LOGODetect regions -- ##
bc_dep_regions = fread("logodetect_regions/brca_depression_logodetect_regions.txt")
bc_hdl_regions = fread("logodetect_regions/brca_hdl_logodetect_regions.txt")
bc_ldl_regions = fread("logodetect_regions/brca_ldl_logodetect_regions.txt")
bc_prostate_regions = fread("logodetect_regions/brca_prostate_logodetect_regions.txt")
bc_schizophrenia_regions = fread("logodetect_regions/brca_schizophrenia_regions.txt")
bc_t2dm_diamante_regions = fread("logodetect_regions/brca_t2dm_logodetect_regions.txt")

## -- Filter for positively correlated regions -- ##
bc_dep_regions = bc_dep_regions %>% filter(stat > 0) %>% distinct()
bc_hdl_regions = bc_hdl_regions %>% filter(stat > 0) %>% distinct()
bc_ldl_regions = bc_ldl_regions %>% filter(stat > 0) %>% distinct()
bc_prostate_regions = bc_prostate_regions %>% filter(stat > 0) %>% distinct()
bc_schizophrenia_regions = bc_schizophrenia_regions %>% filter(stat > 0) %>% distinct()
bc_t2dm_diamante_regions = bc_t2dm_diamante_regions %>% filter(stat > 0) %>% distinct()

## -- GWAS pre-processing for FUMA -- ##

## Keep common SNPs in the LOGODetect positively correlated regions for each pair of GWAS
fuma_gwas_input = function(gwas1, gwas2, logodetect_pos_regions) {
  
  # common SNPs
  common_snps = intersect(gwas1$SNP, gwas2$SNP)
  gwas1_common = gwas1 %>% filter(SNP %in% common_snps) %>% dplyr::select(CHR, SNP, BP, P) %>% distinct()
  gwas2_common = gwas2 %>% filter(SNP %in% common_snps) %>% dplyr::select(CHR, SNP, BP, P) %>% distinct()
  
  # create empty data frame to populate
  gwas_common_snps_pos_reg = data.frame(CHR = NA, SNP = NA, BP = NA, P = NA)
  # positively correlated regions
  for (i in 1:nrow(logodetect_pos_regions)) {
    # subset gwas (doesn't matter if its gwas1_common or gwas2_common as both contain the same SNPs)
    # also, p-value doesn't matter in SNP2GENE FUMA as we will pre-define the SNPs of interest (and not choose them based on p-value)
    gwas_subset = gwas1_common %>% filter(CHR == logodetect_pos_regions[i, chr] &
                                            BP >= logodetect_pos_regions[i, begin_pos] &
                                            BP <= logodetect_pos_regions[i, stop_pos])
    gwas_common_snps_pos_reg = rbind(gwas_common_snps_pos_reg, gwas_subset)
  }
  
  # remove first row
  gwas_common_snps_pos_reg = gwas_common_snps_pos_reg[-1, ]
  rownames(gwas_common_snps_pos_reg) = NULL
  
  # sort
  gwas_common_snps_pos_reg = gwas_common_snps_pos_reg %>% arrange(CHR, BP)
  
  # return
  return(gwas_common_snps_pos_reg)
}
fuma_snp_input = function(gwas_common_snps_pos_reg) {
  gwas_common_snps_pos_reg = gwas_common_snps_pos_reg %>% dplyr::select(rsID = SNP,
                                                                        chromosome = CHR, 
                                                                        position = BP)
  return(gwas_common_snps_pos_reg)
}

# BC-Depression
bc_dep_fuma_gwas_input = fuma_gwas_input(bc, depression, bc_dep_regions)
bc_dep_fuma_snp_input = fuma_snp_input(bc_dep_fuma_gwas_input)
# BC-HDL
bc_hdl_fuma_gwas_input = fuma_gwas_input(bc, hdl, bc_hdl_regions)
bc_hdl_fuma_snp_input = fuma_snp_input(bc_hdl_fuma_gwas_input)
# BC-LDL
bc_ldl_fuma_gwas_input = fuma_gwas_input(bc, ldl, bc_ldl_regions)
bc_ldl_fuma_snp_input = fuma_snp_input(bc_ldl_fuma_gwas_input)
# BC-Prostate cancer
bc_prostate_fuma_gwas_input = fuma_gwas_input(bc, prostate, bc_prostate_regions)
bc_prostate_fuma_snp_input = fuma_snp_input(bc_prostate_fuma_gwas_input)
# BC-Schizophrenia
bc_schizophrenia_fuma_gwas_input = fuma_gwas_input(bc, schizophrenia, bc_schizophrenia_regions)
bc_schizophrenia_fuma_snp_input = fuma_snp_input(bc_schizophrenia_fuma_gwas_input)
# BC-T2DM_DIAMANTE
bc_t2dm_diamante_fuma_gwas_input = fuma_gwas_input(bc, t2dm_diamante, bc_t2dm_diamante_regions)
bc_t2dm_diamante_fuma_snp_input = fuma_snp_input(bc_t2dm_diamante_fuma_gwas_input)

## Save
dir.create("fuma_snp2gene")
dir.create("fuma_snp2gene/input")
fwrite(bc_dep_fuma_gwas_input, "fuma_snp2gene_input/bc_dep_fuma_gwas_input.txt", sep = "\t", row.names = FALSE)
fwrite(bc_dep_fuma_snp_input, "fuma_snp2gene_input/bc_dep_fuma_snp_input.txt", sep = "\t", row.names = FALSE)

fwrite(bc_hdl_fuma_gwas_input, "fuma_snp2gene_input/bc_hdl_fuma_gwas_input.txt", sep = "\t", row.names = FALSE)
fwrite(bc_hdl_fuma_snp_input, "fuma_snp2gene_input/bc_hdl_fuma_snp_input.txt", sep = "\t", row.names = FALSE)

fwrite(bc_ldl_fuma_gwas_input, "fuma_snp2gene_input/bc_ldl_fuma_gwas_input.txt", sep = "\t", row.names = FALSE)
fwrite(bc_ldl_fuma_snp_input, "fuma_snp2gene_input/bc_ldl_fuma_snp_input.txt", sep = "\t", row.names = FALSE)

fwrite(bc_prostate_fuma_gwas_input, "fuma_snp2gene_input/bc_prostate_fuma_gwas_input.txt", sep = "\t", row.names = FALSE)
fwrite(bc_prostate_fuma_snp_input, "fuma_snp2gene_input/bc_prostate_fuma_snp_input.txt", sep = "\t", row.names = FALSE)

fwrite(bc_schizophrenia_fuma_gwas_input, "fuma_snp2gene_input/bc_schizophrenia_fuma_gwas_input.txt", sep = "\t", row.names = FALSE)
fwrite(bc_schizophrenia_fuma_snp_input, "fuma_snp2gene_input/bc_schizophrenia_fuma_snp_input.txt", sep = "\t", row.names = FALSE)

fwrite(bc_t2dm_diamante_fuma_gwas_input, "fuma_snp2gene_input/bc_t2dm_diamante_fuma_gwas_input.txt", sep = "\t", row.names = FALSE)
fwrite(bc_t2dm_diamante_fuma_snp_input, "fuma_snp2gene_input/bc_t2dm_diamante_fuma_snp_input.txt", sep = "\t", row.names = FALSE)

## Use the above files as input to FUMA SNP2GENE (https://fuma.ctglab.nl/) to get the genes positionally located within the positively correlated regions of each BRCA-predisposing disease pair
## The mapped genes, for each disease pair, can be found in the the "fuma_snp2gene/mapped_genes" directory --> these are the shared genes for each disease pair
## FUMA SNP2GENE was run using default settings - example of parameters used can be found below (params.config FUMA SNP2GENE file):
# [version]
# FUMA = v1.5.3
# MAGMA = v1.08
# GWAScatalog = e0_r2022-11-29
# ANNOVAR = 2017-07-17
# 
# [inputfiles]
# gwasfile = bc_dep_fuma_gwas_input.txt
# chrcol = CHR
# poscol = BP
# rsIDcol = SNP
# pcol = P
# eacol = NA
# neacol = NA
# orcol = NA
# becol = NA
# secol = NA
# leadSNPsfile = bc_dep_fuma_snp_input.txt
# addleadSNPs = 0
# regionsfile = NA
# 
# [params]
# N = 200000
# Ncol = NA
# exMHC = 1
# MHCopt = annot
# extMHC = NA
# ensembl = v92
# genetype = protein_coding
# leadP = 5e-8
# gwasP = 0.05
# r2 = 0.6
# r2_2 = 0.1
# refpanel = 1KG/Phase3
# pop = EUR
# MAF = 0
# refSNPs = 0
# mergeDist = 250
# 
# [magma]
# magma = 0
# magma_window = NA
# magma_exp = NA
# 
# [posMap]
# posMap = 1
# posMapWindowSize = 10
# posMapAnnot = NA
# posMapCADDth = 0
# posMapRDBth = NA
# posMapChr15 = NA
# posMapChr15Max = NA
# posMapChr15Meth = NA
# posMapAnnoDs = NA
# posMapAnnoMeth = NA
# 
# [eqtlMap]
# eqtlMap = 0
# eqtlMaptss = NA
# eqtlMapSig = 1
# eqtlMapP = 1
# eqtlMapCADDth = 0
# eqtlMapRDBth = NA
# eqtlMapChr15 = NA
# eqtlMapChr15Max = NA
# eqtlMapChr15Meth = NA
# eqtlMapAnnoDs = NA
# eqtlMapAnnoMeth = NA
# 
# [ciMap]
# ciMap = 0
# ciMapBuiltin = NA
# ciMapFileN = 0
# ciMapFiles = NA
# ciMapFDR = NA
# ciMapPromWindow = NA
# ciMapRoadmap = NA
# ciMapEnhFilt = 0
# ciMapPromFilt = 0
# ciMapCADDth = 0
# ciMapRDBth = NA
# ciMapChr15 = NA
# ciMapChr15Max = NA
# ciMapChr15Meth = NA
# ciMapAnnoDs = NA
# ciMapAnnoMeth = NA
