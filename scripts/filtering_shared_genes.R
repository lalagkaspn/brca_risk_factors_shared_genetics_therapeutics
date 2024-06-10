
### Filtering shared genes using MAGMA and S-MultiXcan results ###

## S-MultiXcan ##
## - For each disease, we used its GWAS summary statistics to run S-PrediXcan (default settings)
##    - Instructions about how to run S-PrediXcan can be found here: https://github.com/hakyimlab/MetaXcan/wiki/Tutorial:-GTEx-v8-MASH-models-integration-with-a-Coronary-Artery-Disease-GWAS
## - Then, using the S-PrediXcan results for each tissue, we ran S-MultiXcan (default settings)
##    - Instructions about how to run S-MultiXcan can be found here: https://github.com/hakyimlab/MetaXcan/tree/master
## - We provide the S-MultiXcan results for each disease in the directory "smultixcan_results"
## - If you want to reproduce them, you can follow the instructions in the links above to install and run S-PrediXcan/S-MultiXcan
##    - The input GWAS can be generated with the script "/smultixcan_results/gwas_preprocessing_for_spredixcan_input.R"

## MAGMA ##
## - For each disease, we used its GWAS summary statistics to run MAGMA
## - We provide the MAGMA results for each disease in the directory "magma_results".
## - If you want to reproduce them, you can follow the instructions in this website (https://cncr.nl/research/magma/) to download, install and run MAGMA

library(dplyr)
library(data.table)
library(ggplot2)

## -- load shared genes -- ##
# filtered for those in STRING network
string_high_conf = fread("preprocessed_data/string_high_conf_processed.txt") # this file can be created by running the "string_pre_processing.R" script
genes_in_string_high_conf = unique(c(string_high_conf$from, string_high_conf$to))

bc_ldl = fread("fuma_snp2gene/mapped_genes/brca_ldl_fuma_snp2gene_genes_logodetect_qvalues.txt") %>% filter(entrezID %in% genes_in_string_high_conf)
bc_hdl = fread("fuma_snp2gene/mapped_genes/brca_hdl_fuma_snp2gene_genes_logodetect_qvalues.txt") %>% filter(entrezID %in% genes_in_string_high_conf)
bc_depression = fread("fuma_snp2gene/mapped_genes/brca_depression_fuma_snp2gene_genes_logodetect_qvalues.txt") %>% filter(entrezID %in% genes_in_string_high_conf)
bc_prostate = fread("fuma_snp2gene/mapped_genes/brca_prostate_fuma_snp2gene_genes_logodetect_qvalues.txt") %>% filter(entrezID %in% genes_in_string_high_conf)
bc_scz = fread("fuma_snp2gene/mapped_genes/brca_schizophrenia_fuma_snp2gene_genes_logodetect_qvalues.txt") %>% filter(entrezID %in% genes_in_string_high_conf)
bc_t2dm_diamante = fread("fuma_snp2gene/mapped_genes/brca_t2dm_fuma_snp2gene_genes_logodetect_qvalues.txt") %>% filter(entrezID %in% genes_in_string_high_conf)

## -- load S-MultiXcan results -- ##
bc_smultixcan = fread("smultixcan_results/BRCA_BCAC_EUR_2020_imputed_smultixcan.txt")
depression_smultixcan = fread("smultixcan_results/DEPRESSION_PGC_UKBB_EUR_2019_imputed_smultixcan.txt")
hdl_smultixcan = fread("smultixcan_results/HDL_GLGC_EUR_2021_imputed_smultixcan.txt")
ldl_smultixcan = fread("smultixcan_results/LDL_GLGC_EUR_2021_imputed_smultixcan.txt")
prostate_smultixcan = fread("smultixcan_results/PROSTATE_CANCER_PRACTICAL_EUR_2018_imputed_smultixcan.txt")
scz_smultixcan = fread("smultixcan_results/SCHIZOPHRENIA_PGC_EUR_2022_imputed_smultixcan.txt")
t2dm_diamante_smultixcan = fread("smultixcan_results/T2DM_DIAMANTE_2022_imputed_smultixcan.txt")

## -- load MAGMA results -- ##
bc_magma = fread("magma_results/bc.genes.out")
depression_magma = fread("magma_results/depression.genes.out")
hdl_magma = fread("magma_results/hdl.genes.out")
ldl_magma = fread("magma_results/ldl.genes.out")
prostate_magma = fread("magma_results/prostate.genes.out")
scz_magma = fread("magma_results/schizophrenia.genes.out")
t2dm_diamante_magma = fread("magma_results/t2dm_diamante.genes.out")
# # adjust MAGMA p-values for multiple hypothesis testing
# bc_magma$P_adj_BH = p.adjust(bc_magma$P, method = "BH", n = nrow(bc_magma))
# depression_magma$P_adj_BH = p.adjust(depression_magma$P, method = "BH", n = nrow(depression_magma))
# hdl_magma$P_adj_BH = p.adjust(hdl_magma$P, method = "BH", n = nrow(hdl_magma))
# ldl_magma$P_adj_BH = p.adjust(ldl_magma$P, method = "BH", n = nrow(ldl_magma))
# prostate_magma$P_adj_BH = p.adjust(prostate_magma$P, method = "BH", n = nrow(prostate_magma))
# scz_magma$P_adj_BH = p.adjust(scz_magma$P, method = "BH", n = nrow(scz_magma))
# t2dm_diamante_magma$P_adj_BH = p.adjust(t2dm_diamante_magma$P, method = "BH", n = nrow(t2dm_diamante_magma))

## -- filter for MAGMA significant genes for the BC-related disease and same S-MultiXcan direction -- ##

## BC-LDL
# filter for ldl-MAGMA significant genes (after adjusting for the number of original LOGODetect genes)
ldl_magma_shared = ldl_magma %>% 
  filter(GENE %in% bc_ldl$entrezID) %>%
  mutate(P.adj = p.adjust(P, method = "BH", n = length(P))) %>%
  filter(P.adj < 0.05) %>%
  dplyr::select(GENE) %>%
  distinct()
# filter bc-smultixcan genes for the LOGODetect ones and keep the direction of z
bc_smultixcan_shared = bc_smultixcan %>%
  filter(gene_name %in% bc_ldl$gene_name) %>%
  dplyr::select(gene_name, z_mean) %>%
  distinct() %>%
  mutate(z_sign = if_else(z_mean > 0, 1, 0)) %>%
  dplyr::select(-z_mean) %>%
  distinct()
# filter ldl-smultixcan genes for the LOGODetect ones and keep the direction of z
ldl_smultixcan_shared = ldl_smultixcan %>% 
  filter(gene_name %in% bc_ldl$gene_name) %>%
  dplyr::select(gene_name, z_mean) %>%
  distinct() %>%
  mutate(z_sign = if_else(z_mean > 0, 1, 0)) %>%
  dplyr::select(-z_mean) %>%
  distinct()
# combine bc and ldl smultixcan genes and keep the ones that have the same direction (either both positive or both negative)
bc_ldl_smultixcan_shared_combined = left_join(bc_smultixcan_shared, ldl_smultixcan_shared, by = "gene_name") %>%
  filter(z_sign.x == 1 & z_sign.y == 1 | z_sign.x == 0 & z_sign.y == 0) ; rm(bc_smultixcan_shared, ldl_smultixcan_shared)
bc_ldl_magma_smultixcan_dir_shared = bc_ldl %>%
  filter(entrezID %in% ldl_magma_shared$GENE) %>%
  left_join(bc_ldl_smultixcan_shared_combined, by = "gene_name") %>%
  na.omit() %>%
  dplyr::select(-z_sign.x, -z_sign.y) %>%
  distinct() 
rm(ldl_magma_shared, bc_ldl_smultixcan_shared_combined)
# 60 / 73 genes are ldl MAGMA significant and have the same direction (positive or negative) as BC in S-MultiXcan
fwrite(bc_ldl_magma_smultixcan_dir_shared, "magma_smultixcan_filtered_shared_genes/bc_ldl_filtered_shared_genes.txt", sep = "\t", row.names = FALSE)

## BC-HDL
# filter for hdl-MAGMA significant genes (after adjusting for the number of original LOGODetect genes)
hdl_magma_shared = hdl_magma %>% 
  filter(GENE %in% bc_hdl$entrezID) %>%
  mutate(P.adj = p.adjust(P, method = "BH", n = length(P))) %>%
  filter(P.adj < 0.05) %>%
  dplyr::select(GENE) %>%
  distinct()
# filter bc-smultixcan genes for the LOGODetect ones and keep the direction of z
bc_smultixcan_shared = bc_smultixcan %>%
  filter(gene_name %in% bc_hdl$gene_name) %>%
  dplyr::select(gene_name, z_mean) %>%
  distinct() %>%
  mutate(z_sign = if_else(z_mean > 0, 1, 0)) %>%
  dplyr::select(-z_mean) %>%
  distinct()
# filter hdl-smultixcan genes for the LOGODetect ones and keep the direction of z
hdl_smultixcan_shared = hdl_smultixcan %>% 
  filter(gene_name %in% bc_hdl$gene_name) %>%
  dplyr::select(gene_name, z_mean) %>%
  distinct() %>%
  mutate(z_sign = if_else(z_mean > 0, 1, 0)) %>%
  dplyr::select(-z_mean) %>%
  distinct()
# combine bc and hdl smultixcan genes and keep the ones that have the same direction (either both positive or both negative)
bc_hdl_smultixcan_shared_combined = left_join(bc_smultixcan_shared, hdl_smultixcan_shared, by = "gene_name") %>%
  filter(z_sign.x == 1 & z_sign.y == 1 | z_sign.x == 0 & z_sign.y == 0) ; rm(bc_smultixcan_shared, hdl_smultixcan_shared)
bc_hdl_magma_smultixcan_dir_shared = bc_hdl %>%
  filter(entrezID %in% hdl_magma_shared$GENE) %>%
  left_join(bc_hdl_smultixcan_shared_combined, by = "gene_name") %>%
  na.omit() %>%
  dplyr::select(-z_sign.x, -z_sign.y) %>%
  distinct() 
rm(hdl_magma_shared, bc_hdl_smultixcan_shared_combined)
# 3 / 7 genes are hdl MAGMA significant and have the same direction (positive or negative) as BC in S-MultiXcan
fwrite(bc_hdl_magma_smultixcan_dir_shared, "magma_smultixcan_filtered_shared_genes/bc_hdl_filtered_shared_genes.txt", sep = "\t", row.names = FALSE)

## BC-DEPRESSION
# filter for depression-MAGMA significant genes (after adjusting for the number of original LOGODetect genes)
depression_magma_shared = depression_magma %>% 
  filter(GENE %in% bc_depression$entrezID) %>%
  mutate(P.adj = p.adjust(P, method = "BH", n = length(P))) %>%
  filter(P.adj < 0.05) %>%
  dplyr::select(GENE) %>%
  distinct()
# filter bc-smultixcan genes for the LOGODetect ones and keep the direction of z
bc_smultixcan_shared = bc_smultixcan %>%
  filter(gene_name %in% bc_depression$gene_name) %>%
  dplyr::select(gene_name, z_mean) %>%
  distinct() %>%
  mutate(z_sign = if_else(z_mean > 0, 1, 0)) %>%
  dplyr::select(-z_mean) %>%
  distinct()
# filter depression-smultixcan genes for the LOGODetect ones and keep the direction of z
depression_smultixcan_shared = depression_smultixcan %>% 
  filter(gene_name %in% bc_depression$gene_name) %>%
  dplyr::select(gene_name, z_mean) %>%
  distinct() %>%
  mutate(z_sign = if_else(z_mean > 0, 1, 0)) %>%
  dplyr::select(-z_mean) %>%
  distinct()
# combine bc and depression smultixcan genes and keep the ones that have the same direction (either both positive or both negative)
bc_depression_smultixcan_shared_combined = left_join(bc_smultixcan_shared, depression_smultixcan_shared, by = "gene_name") %>%
  filter(z_sign.x == 1 & z_sign.y == 1 | z_sign.x == 0 & z_sign.y == 0) ; rm(bc_smultixcan_shared, depression_smultixcan_shared)
bc_depression_magma_smultixcan_dir_shared = bc_depression %>%
  filter(entrezID %in% depression_magma_shared$GENE) %>%
  left_join(bc_depression_smultixcan_shared_combined, by = "gene_name") %>%
  na.omit() %>%
  dplyr::select(-z_sign.x, -z_sign.y) %>%
  distinct() 
rm(depression_magma_shared, bc_depression_smultixcan_shared_combined)
# 3 / 10 genes are depression MAGMA significant and have the same direction (positive or negative) as BC in S-MultiXcan
fwrite(bc_depression_magma_smultixcan_dir_shared, "magma_smultixcan_filtered_shared_genes//bc_depression_filtered_shared_genes.txt", sep = "\t", row.names = FALSE)

## BC-PROSTATE
# filter for prostate-MAGMA significant genes (after adjusting for the number of original LOGODetect genes)
prostate_magma_shared = prostate_magma %>% 
  filter(GENE %in% bc_prostate$entrezID) %>%
  mutate(P.adj = p.adjust(P, method = "BH", n = length(P))) %>%
  filter(P.adj < 0.05) %>%
  dplyr::select(GENE) %>%
  distinct()
# filter bc-smultixcan genes for the LOGODetect ones and keep the direction of z
bc_smultixcan_shared = bc_smultixcan %>%
  filter(gene_name %in% bc_prostate$gene_name) %>%
  dplyr::select(gene_name, z_mean) %>%
  distinct() %>%
  mutate(z_sign = if_else(z_mean > 0, 1, 0)) %>%
  dplyr::select(-z_mean) %>%
  distinct()
# filter prostate-smultixcan genes for the LOGODetect ones and keep the direction of z
prostate_smultixcan_shared = prostate_smultixcan %>% 
  filter(gene_name %in% bc_prostate$gene_name) %>%
  dplyr::select(gene_name, z_mean) %>%
  distinct() %>%
  mutate(z_sign = if_else(z_mean > 0, 1, 0)) %>%
  dplyr::select(-z_mean) %>%
  distinct()
# combine bc and prostate smultixcan genes and keep the ones that have the same direction (either both positive or both negative)
bc_prostate_smultixcan_shared_combined = left_join(bc_smultixcan_shared, prostate_smultixcan_shared, by = "gene_name") %>%
  filter(z_sign.x == 1 & z_sign.y == 1 | z_sign.x == 0 & z_sign.y == 0) ; rm(bc_smultixcan_shared, prostate_smultixcan_shared)
bc_prostate_magma_smultixcan_dir_shared = bc_prostate %>%
  filter(entrezID %in% prostate_magma_shared$GENE) %>%
  left_join(bc_prostate_smultixcan_shared_combined, by = "gene_name") %>%
  na.omit() %>%
  dplyr::select(-z_sign.x, -z_sign.y) %>%
  distinct() 
rm(prostate_magma_shared, bc_prostate_smultixcan_shared_combined)
# 38 / 46 genes are prostate MAGMA significant and have the same direction (positive or negative) as BC in S-MultiXcan
fwrite(bc_prostate_magma_smultixcan_dir_shared, "magma_smultixcan_filtered_shared_genes/bc_prostate_filtered_shared_genes.txt", sep = "\t", row.names = FALSE)

## BC-SCZ
# filter for scz-MAGMA significant genes (after adjusting for the number of original LOGODetect genes)
scz_magma_shared = scz_magma %>% 
  filter(GENE %in% bc_scz$entrezID) %>%
  mutate(P.adj = p.adjust(P, method = "BH", n = length(P))) %>%
  filter(P.adj < 0.05) %>%
  dplyr::select(GENE) %>%
  distinct()
# filter bc-smultixcan genes for the LOGODetect ones and keep the direction of z
bc_smultixcan_shared = bc_smultixcan %>%
  filter(gene_name %in% bc_scz$gene_name) %>%
  dplyr::select(gene_name, z_mean) %>%
  distinct() %>%
  mutate(z_sign = if_else(z_mean > 0, 1, 0)) %>%
  dplyr::select(-z_mean) %>%
  distinct()
# filter scz-smultixcan genes for the LOGODetect ones and keep the direction of z
scz_smultixcan_shared = scz_smultixcan %>% 
  filter(gene_name %in% bc_scz$gene_name) %>%
  dplyr::select(gene_name, z_mean) %>%
  distinct() %>%
  mutate(z_sign = if_else(z_mean > 0, 1, 0)) %>%
  dplyr::select(-z_mean) %>%
  distinct()
# combine bc and scz smultixcan genes and keep the ones that have the same direction (either both positive or both negative)
bc_scz_smultixcan_shared_combined = left_join(bc_smultixcan_shared, scz_smultixcan_shared, by = "gene_name") %>%
  filter(z_sign.x == 1 & z_sign.y == 1 | z_sign.x == 0 & z_sign.y == 0) ; rm(bc_smultixcan_shared, scz_smultixcan_shared)
bc_scz_magma_smultixcan_dir_shared = bc_scz %>%
  filter(entrezID %in% scz_magma_shared$GENE) %>%
  left_join(bc_scz_smultixcan_shared_combined, by = "gene_name") %>%
  na.omit() %>%
  dplyr::select(-z_sign.x, -z_sign.y) %>%
  distinct() 
rm(scz_magma_shared, bc_scz_smultixcan_shared_combined)
# 13 / 21 genes are scz MAGMA significant and have the same direction (positive or negative) as BC in S-MultiXcan
fwrite(bc_scz_magma_smultixcan_dir_shared, "magma_smultixcan_filtered_shared_genes/bc_schizophrenia_filtered_shared_genes.txt", sep = "\t", row.names = FALSE)

## BC-T2DM_DIAMANTE
# filter for t2dm_diamante-MAGMA significant genes (after adjusting for the number of original LOGODetect genes)
t2dm_diamante_magma_shared = t2dm_diamante_magma %>% 
  filter(GENE %in% bc_t2dm_diamante$entrezID) %>%
  mutate(P.adj = p.adjust(P, method = "BH", n = length(P))) %>%
  filter(P.adj < 0.05) %>%
  dplyr::select(GENE) %>%
  distinct()
# filter bc-smultixcan genes for the LOGODetect ones and keep the direction of z
bc_smultixcan_shared = bc_smultixcan %>%
  filter(gene_name %in% bc_t2dm_diamante$gene_name) %>%
  dplyr::select(gene_name, z_mean) %>%
  distinct() %>%
  mutate(z_sign = if_else(z_mean > 0, 1, 0)) %>%
  dplyr::select(-z_mean) %>%
  distinct()
# filter t2dm_diamante-smultixcan genes for the LOGODetect ones and keep the direction of z
t2dm_diamante_smultixcan_shared = t2dm_diamante_smultixcan %>% 
  filter(gene_name %in% bc_t2dm_diamante$gene_name) %>%
  dplyr::select(gene_name, z_mean) %>%
  distinct() %>%
  mutate(z_sign = if_else(z_mean > 0, 1, 0)) %>%
  dplyr::select(-z_mean) %>%
  distinct()
# combine bc and t2dm_diamante smultixcan genes and keep the ones that have the same direction (either both positive or both negative)
bc_t2dm_diamante_smultixcan_shared_combined = left_join(bc_smultixcan_shared, t2dm_diamante_smultixcan_shared, by = "gene_name") %>%
  filter(z_sign.x == 1 & z_sign.y == 1 | z_sign.x == 0 & z_sign.y == 0) ; rm(bc_smultixcan_shared, t2dm_diamante_smultixcan_shared)
bc_t2dm_diamante_magma_smultixcan_dir_shared = bc_t2dm_diamante %>%
  filter(entrezID %in% t2dm_diamante_magma_shared$GENE) %>%
  left_join(bc_t2dm_diamante_smultixcan_shared_combined, by = "gene_name") %>%
  na.omit() %>%
  dplyr::select(-z_sign.x, -z_sign.y) %>%
  distinct() 
rm(t2dm_diamante_magma_shared, bc_t2dm_diamante_smultixcan_shared_combined)
# 52 / 68 genes are t2dm_diamante MAGMA significant and have the same direction (positive or negative) as BC in S-MultiXcan
fwrite(bc_t2dm_diamante_magma_smultixcan_dir_shared, "magma_smultixcan_filtered_shared_genes/bc_t2dm_diamante_filtered_shared_genes.txt", sep = "\t", row.names = FALSE)
