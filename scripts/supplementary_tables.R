
### Create supplementary tables of the manuscript

library(dplyr)
library(data.table)
library(openxlsx)

# we will only provide genes that are in the STRING network
string_high_conf = fread("preprocessed_data/string_high_conf_processed.txt")
genes_in_string = unique(c(string_high_conf$from, string_high_conf$to))

## Supplementary Table 1-2 - Shared loci between BRCA and predisposing diseases ##
# logodetect loci
bc_depression_loci = fread("logodetect_regions/brca_depression_logodetect_regions.txt") %>% dplyr::rename(begin_pos_GRCh37 = begin_pos, stop_pos_GRCh37 = stop_pos)
bc_hdl_loci = fread("logodetect_regions/brca_hdl_logodetect_regions.txt") %>% dplyr::rename(begin_pos_GRCh37 = begin_pos, stop_pos_GRCh37 = stop_pos)
bc_ldl_loci = fread("logodetect_regions/brca_ldl_logodetect_regions.txt") %>% dplyr::rename(begin_pos_GRCh37 = begin_pos, stop_pos_GRCh37 = stop_pos)
bc_prostate_loci = fread("logodetect_regions/brca_prostate_logodetect_regions.txt") %>% dplyr::rename(begin_pos_GRCh37 = begin_pos, stop_pos_GRCh37 = stop_pos)
bc_scz_loci = fread("logodetect_regions/brca_schizophrenia_regions.txt") %>% dplyr::rename(begin_pos_GRCh37 = begin_pos, stop_pos_GRCh37 = stop_pos)
bc_t2dm_loci = fread("logodetect_regions/brca_t2dm_logodetect_regions.txt") %>% dplyr::rename(begin_pos_GRCh37 = begin_pos, stop_pos_GRCh37 = stop_pos)

# genes in loci
bc_depression_genes = fread("fuma_snp2gene/mapped_genes/brca_depression_fuma_snp2gene_genes.txt") %>% dplyr::select(symbol, entrezID, chr, start, end) %>% distinct() %>% na.omit() %>% mutate(disease = "Depression") %>% filter(entrezID %in% genes_in_string)
bc_hdl_genes = fread("fuma_snp2gene/mapped_genes/brca_hdl_fuma_snp2gene_genes.txt") %>% dplyr::select(symbol, entrezID, chr, start, end) %>% distinct() %>% na.omit() %>% mutate(disease = "HDL") %>% filter(entrezID %in% genes_in_string)
bc_ldl_genes = fread("fuma_snp2gene/mapped_genes/brca_ldl_fuma_snp2gene_genes.txt") %>% dplyr::select(symbol, entrezID, chr, start, end) %>% distinct() %>% na.omit() %>% mutate(disease = "LDL") %>% filter(entrezID %in% genes_in_string)
bc_prostate_genes = fread("fuma_snp2gene/mapped_genes/brca_prostate_fuma_snp2gene_genes.txt") %>% dplyr::select(symbol, entrezID, chr, start, end) %>% distinct() %>% na.omit() %>% mutate(disease = "Prostate") %>% filter(entrezID %in% genes_in_string)
bc_scz_genes = fread("fuma_snp2gene/mapped_genes/brca_schizophrenia_fuma_snp2gene_genes.txt") %>% dplyr::select(symbol, entrezID, chr, start, end) %>% distinct() %>% na.omit() %>% mutate(disease = "Scz") %>% filter(entrezID %in% genes_in_string)
bc_t2dm_genes = fread("fuma_snp2gene/mapped_genes/brca_t2dm_fuma_snp2gene_genes.txt") %>% dplyr::select(symbol, entrezID, chr, start, end) %>% distinct() %>% na.omit() %>% mutate(disease = "T2DM") %>% filter(entrezID %in% genes_in_string)

# positively correlated loci
bc_depression_loci_pos = bc_depression_loci %>% filter(stat > 0)
bc_hdl_loci_pos = bc_hdl_loci %>% filter(stat > 0)
bc_ldl_loci_pos = bc_ldl_loci %>% filter(stat > 0)
bc_prostate_loci_pos = bc_prostate_loci %>% filter(stat > 0)
bc_scz_loci_pos = bc_scz_loci %>% filter(stat > 0)
bc_t2dm_loci_pos = bc_t2dm_loci %>% filter(stat > 0)

# negatively correlated loci
bc_depression_loci_neg = bc_depression_loci %>% filter(stat < 0)
bc_hdl_loci_neg = bc_hdl_loci %>% filter(stat < 0)
bc_ldl_loci_neg = bc_ldl_loci %>% filter(stat < 0)
bc_prostate_loci_neg = bc_prostate_loci %>% filter(stat < 0)
bc_scz_loci_neg = bc_scz_loci %>% filter(stat < 0)
bc_t2dm_loci_neg = bc_t2dm_loci %>% filter(stat < 0)

# ultimate file for positively correlated loci
pos_loci_with_genes = function(pos_loci, genes, disease) {
  pos_loci = pos_loci %>% mutate("Disease 1" = rep("Breast cancer", nrow(pos_loci)), "Disease 2" = rep(disease, nrow(pos_loci)), .before = chr)
  # add genes
  for (i in 1:nrow(pos_loci)) {
    chr_temp = pos_loci[i, chr]
    begin_temp = pos_loci[i, begin_pos_GRCh37] - 10000 # default from FUMA SNP2GENE
    stop_temp = pos_loci[i, stop_pos_GRCh37] + 10000 # default from FUMA SNP2GENE
    
    genes_temp = genes %>% filter(chr == chr_temp, 
                                  start <= begin_temp & end >= begin_temp | # part of gene is within the start of the locus and rest of it is outsise of the locus
                                    start >= begin_temp & end <= stop_temp | # all gene is within the locus
                                    start <= stop_temp & end >= stop_temp) # part of the gene is close to the end of the locus and rest of it outside of locus
    
    if (nrow(genes_temp) > 0) {
      genes_temp$symbol_entrez = paste0(genes_temp$symbol, " (", genes_temp$entrezID, ")")
      pos_loci[i, "Gene in locus (entrezID)"] = paste0(genes_temp$symbol_entrez, collapse = " | ")
    }
  }
  colnames(pos_loci) = c("Disease 1", "Disease 2", "Chromosome", "Begin_pos_GRCh37", "Stop_pos_GRCh37", "LOGODetect_stat", "LOGODetect_pvalue", "LOGODetect_qvalue", "Genes in locus (entrezID)")
  return(pos_loci)
}

bc_depression_loci_pos = pos_loci_with_genes(pos_loci = bc_depression_loci_pos, genes = bc_depression_genes, disease = "Depression")
bc_hdl_loci_pos = pos_loci_with_genes(pos_loci = bc_hdl_loci_pos, genes = bc_hdl_genes, disease = "high HDL")
bc_ldl_loci_pos = pos_loci_with_genes(pos_loci = bc_ldl_loci_pos, genes = bc_ldl_genes, disease = "high LDL")
bc_prostate_loci_pos = pos_loci_with_genes(pos_loci = bc_prostate_loci_pos, genes = bc_prostate_genes, disease = "Prostate cancer")
bc_scz_loci_pos = pos_loci_with_genes(pos_loci = bc_scz_loci_pos, genes = bc_scz_genes, disease = "Schizophrenia")
bc_t2dm_loci_pos = pos_loci_with_genes(pos_loci = bc_t2dm_loci_pos, genes = bc_t2dm_genes, disease = "Type 2 diabetes")

pos_loci_ultimate = rbind(bc_depression_loci_pos, bc_hdl_loci_pos, bc_ldl_loci_pos, bc_prostate_loci_pos, bc_scz_loci_pos, bc_t2dm_loci_pos)

# ultimate file for negatively correlated loci
neg_loci = function(neg_loci, disease) {
  neg_loci = neg_loci %>% mutate("Disease 1" = rep("Breast cancer", nrow(neg_loci)), "Disease 2" = rep(disease, nrow(neg_loci)), .before = chr)
  colnames(neg_loci) = c("Disease 1", "Disease 2", "Chromosome", "Begin_pos_GRCh37", "Stop_pos_GRCh37", "LOGODetect_stat", "LOGODetect_pvalue", "LOGODetect_qvalue")
  return(neg_loci)
}
bc_depression_loci_neg = neg_loci(neg_loci = bc_depression_loci_neg, disease = "Depression")
bc_hdl_loci_neg = neg_loci(neg_loci = bc_hdl_loci_neg, disease = "high HDL")
bc_ldl_loci_neg = neg_loci(neg_loci = bc_ldl_loci_neg, disease = "high LDL")
bc_prostate_loci_neg = neg_loci(neg_loci = bc_prostate_loci_neg, disease = "Prostate cancer")
bc_scz_loci_neg = neg_loci(neg_loci = bc_scz_loci_neg, disease = "Schizophrenia")
bc_t2dm_loci_neg = neg_loci(neg_loci = bc_t2dm_loci_neg, disease = "Type 2 diabetes")

neg_loci_ultimate = rbind(bc_depression_loci_neg, bc_hdl_loci_neg, bc_ldl_loci_neg, bc_prostate_loci_neg, bc_scz_loci_neg, bc_t2dm_loci_neg)

# save
write.xlsx(pos_loci_ultimate, file = "supplementary_tables/S1_positively_correlated_loci.xlsx")
write.xlsx(neg_loci_ultimate, file = "supplementary_tables/S2_negatively_correlated_loci.xlsx")


## Supplementary Table 3 - genes in positively correlated loci + MAGMA/S-PrediXcan + canonical pathways they are connected to
# MAGMA results
depression_magma = fread("magma_results/depression.genes.out") %>% filter(GENE %in% bc_depression_genes$entrezID) %>% dplyr::select(GENE, P) %>% filter(GENE %in% genes_in_string)
hdl_magma = fread("magma_results/hdl.genes.out") %>% filter(GENE %in% bc_hdl_genes$entrezID) %>% dplyr::select(GENE, P) %>% filter(GENE %in% genes_in_string)
ldl_magma = fread("magma_results/ldl.genes.out") %>% filter(GENE %in% bc_ldl_genes$entrezID) %>% dplyr::select(GENE, P) %>% filter(GENE %in% genes_in_string)
prostate_magma = fread("magma_results/prostate.genes.out") %>% filter(GENE %in% bc_prostate_genes$entrezID) %>% dplyr::select(GENE, P) %>% filter(GENE %in% genes_in_string)
scz_magma = fread("magma_results/schizophrenia.genes.out") %>% filter(GENE %in% bc_scz_genes$entrezID) %>% dplyr::select(GENE, P) %>% filter(GENE %in% genes_in_string)
t2dm_magma = fread("magma_results/t2dm_diamante.genes.out") %>% filter(GENE %in% bc_t2dm_genes$entrezID) %>% dplyr::select(GENE, P) %>% filter(GENE %in% genes_in_string)
# S-MultiXcan results
bc_smultixcan = fread("smultixcan_results/BRCA_BCAC_EUR_2020_imputed_smultixcan.txt") %>% dplyr::select(gene_name, pvalue, z_mean) %>% distinct()
depression_smultixcan = fread("smultixcan_results/DEPRESSION_PGC_UKBB_EUR_2019_imputed_smultixcan.txt") %>% dplyr::select(gene_name, pvalue, z_mean) %>% distinct() %>% filter(gene_name %in% bc_depression_genes$symbol) %>% left_join(bc_depression_genes[, c(1,2)], by = c("gene_name" = "symbol"))
hdl_smultixcan = fread("smultixcan_results/HDL_GLGC_EUR_2021_imputed_smultixcan.txt") %>% dplyr::select(gene_name, pvalue, z_mean) %>% distinct() %>% filter(gene_name %in% bc_hdl_genes$symbol) %>% left_join(bc_hdl_genes[, c(1,2)], by = c("gene_name" = "symbol"))
ldl_smultixcan = fread("smultixcan_results/LDL_GLGC_EUR_2021_imputed_smultixcan.txt") %>% dplyr::select(gene_name, pvalue, z_mean) %>% distinct() %>% filter(gene_name %in% bc_ldl_genes$symbol) %>% left_join(bc_ldl_genes[, c(1,2)], by = c("gene_name" = "symbol"))
prostate_smultixcan = fread("smultixcan_results/PROSTATE_CANCER_PRACTICAL_EUR_2018_imputed_smultixcan.txt") %>% dplyr::select(gene_name, pvalue, z_mean) %>% distinct() %>% filter(gene_name %in% bc_prostate_genes$symbol) %>% left_join(bc_prostate_genes[, c(1,2)], by = c("gene_name" = "symbol"))
scz_smultixcan = fread("smultixcan_results/SCHIZOPHRENIA_PGC_EUR_2022_imputed_smultixcan.txt") %>% dplyr::select(gene_name, pvalue, z_mean) %>% distinct() %>% filter(gene_name %in% bc_scz_genes$symbol) %>% left_join(bc_scz_genes[, c(1,2)], by = c("gene_name" = "symbol"))
t2dm_smultixcan = fread("smultixcan_results/T2DM_DIAMANTE_2022_imputed_smultixcan.txt") %>% dplyr::select(gene_name, pvalue, z_mean) %>% distinct() %>% filter(gene_name %in% bc_t2dm_genes$symbol) %>% left_join(bc_t2dm_genes[, c(1,2)], by = c("gene_name" = "symbol"))

magma_smultixcan_data = function(disease1_smultixcan, disease2_smultixcan, disease2_magma, disease2) {
  
  # annotate genes with entrezIDs
  genes_with_entrez = disease2_smultixcan %>% dplyr::select(gene_name, entrezID)
  disease1_smultixcan = left_join(disease1_smultixcan, genes_with_entrez, by = c("gene_name")) %>% na.omit()
  disease2_magma = left_join(disease2_magma, genes_with_entrez, by = c("GENE" = "entrezID")) %>% na.omit()
  
  # adjust MAGMA p-values
  disease2_magma$P_adj_BH = p.adjust(disease2_magma$P, method = "BH")
  
  disease2_magma = disease2_magma %>%
    dplyr::select("Shared gene symbol" = gene_name, "Shared gene entrezID" = GENE, MAGMA_pvalue_adj_disease2 = "P_adj_BH") %>%
    mutate("Disease 1" = rep("Breast cancer", nrow(disease2_magma)),
           "Disease 2" = rep(disease2, nrow(disease2_magma)), .before = "Shared gene symbol")
  # add S-MultiXcan information
  disease1_smultixcan = disease1_smultixcan %>% dplyr::select(entrezID, z_mean) %>% distinct()
  disease2_smultixcan = disease2_smultixcan %>% dplyr::select(entrezID, z_mean) %>% distinct()
  disease2_magma = left_join(disease2_magma, disease1_smultixcan, by = c("Shared gene entrezID" = "entrezID"))
  disease2_magma = left_join(disease2_magma, disease2_smultixcan, by = c("Shared gene entrezID" = "entrezID"))
  disease2_magma = disease2_magma %>% dplyr::rename("S-MultiXcan disease 1 (z-score)" = "z_mean.x",
                                                    "S-MultiXcan disease 2 (z-score)" = "z_mean.y")
  # add column for passing or not filters
  disease2_magma$`Pass MAGMA filter` = ifelse(disease2_magma$MAGMA_pvalue_adj_disease2 < 0.05, 1, 0)
  disease2_magma$`Pass S-MultiXcan filter` = ifelse(disease2_magma$`S-MultiXcan disease 1 (z-score)` * disease2_magma$`S-MultiXcan disease 2 (z-score)` > 0, 1, 0)
  disease2_magma$`Pass both filters` = ifelse(disease2_magma$`Pass MAGMA filter` * disease2_magma$`Pass S-MultiXcan filter` == 1, 1, 0)
  
  return(disease2_magma)
}

bc_depression_genes_filtered = magma_smultixcan_data(disease1_smultixcan = bc_smultixcan, disease2_smultixcan = depression_smultixcan, disease2_magma = depression_magma, disease2 = "Depression")
bc_hdl_genes_filtered = magma_smultixcan_data(disease1_smultixcan = bc_smultixcan, disease2_smultixcan = hdl_smultixcan, disease2_magma = hdl_magma, disease2 = "High HDL")
bc_ldl_genes_filtered = magma_smultixcan_data(disease1_smultixcan = bc_smultixcan, disease2_smultixcan = ldl_smultixcan, disease2_magma = ldl_magma, disease2 = "High LDL")
bc_prostate_genes_filtered = magma_smultixcan_data(disease1_smultixcan = bc_smultixcan, disease2_smultixcan = prostate_smultixcan, disease2_magma = prostate_magma, disease2 = "Prostate cancer")
bc_scz_genes_filtered = magma_smultixcan_data(disease1_smultixcan = bc_smultixcan, disease2_smultixcan = scz_smultixcan, disease2_magma = scz_magma, disease2 = "Schizophrenia")
bc_t2dm_genes_filtered = magma_smultixcan_data(disease1_smultixcan = bc_smultixcan, disease2_smultixcan = t2dm_smultixcan, disease2_magma = t2dm_magma, disease2 = "Type 2 diabetes")

genes_pos_loci_magma_smultixcan = rbind(bc_depression_genes_filtered, bc_hdl_genes_filtered, bc_ldl_genes_filtered, bc_prostate_genes_filtered, bc_scz_genes_filtered, bc_t2dm_genes_filtered)
genes_pos_loci_magma_smultixcan$`Pass MAGMA filter` = ifelse(genes_pos_loci_magma_smultixcan$`Pass MAGMA filter` == 1, "Yes", "No")
genes_pos_loci_magma_smultixcan$`Pass S-MultiXcan filter` = ifelse(genes_pos_loci_magma_smultixcan$`Pass S-MultiXcan filter` == 1, "Yes", "No")
genes_pos_loci_magma_smultixcan$`Pass both filters` = ifelse(genes_pos_loci_magma_smultixcan$`Pass both filters` == 1, "Yes", "No")

write.xlsx(genes_pos_loci_magma_smultixcan, file = "supplementary_tables/S3_magma_smultixcan_filters.xlsx")


## Supplementary Table 4 - Filtered shared genes - canonical pathways ##
shared_genes_to_canonical_pathways = function(permuted_values_path, shared_genes, disease) {
  permuted_values = readRDS(permuted_values_path)
  sig_cp = permuted_values
  for (i in 1:length(sig_cp)) {
    sig_cp[[i]] = sig_cp[[i]] %>%
      mutate(ppr_pvalue_adj = p.adjust(perm_pvalue, "bonferroni", length(perm_pvalue))) %>%
      mutate(ppr_pvalue_adj = min(ppr_pvalue_adj)) %>%
      dplyr::select(ppr_pvalue_adj) %>%
      distinct()
    cat(i, "\n")
  } ; rm(i)
  sig_cp = 
    lapply(names(sig_cp), function(pathway) {
      sig_cp[[pathway]] %>% 
        mutate(pathway = pathway)
    }) %>%
    bind_rows() %>%
    select(pathway, shared_gene, ppr_pvalue_adj) %>%
    group_by(shared_gene) %>%
    mutate(ppr_pvalue_adj = p.adjust(ppr_pvalue_adj, method = "bonferroni")) %>%
    ungroup() %>%
    filter(ppr_pvalue_adj < 0.05)
  sig_cp = sig_cp %>%
    mutate("Disease 1" = rep("Breast cancer", nrow(sig_cp)),
           "Disease 2" = rep(disease , nrow(sig_cp)), .before = pathway)  %>%
    dplyr::select("Disease 1", "Disease 2", "Shared gene entrezID" = shared_gene, "Connected canonical pathway" = pathway)
  # add gene symbol
  sig_cp = left_join(sig_cp, shared_genes[, c(1,2)], by = c("Shared gene entrezID" = "entrezID"))
  sig_cp = sig_cp %>% dplyr::select("Disease 1", "Disease 2", "Shared gene symbol" = symbol, "Shared gene entrezID", "Connected canonical pathway")
  return(sig_cp)
}

bc_depression_sig_cp = shared_genes_to_canonical_pathways(permuted_values_path = "shared_genes_to_canonical_pathways/permuted_BC_DEPRESSION.RDS", 
                                                          shared_genes = bc_depression_genes, disease = "Depression")
bc_hdl_sig_cp = shared_genes_to_canonical_pathways(permuted_values_path = "shared_genes_to_canonical_pathways/permuted_BC_HDL.RDS", 
                                                   shared_genes = bc_hdl_genes, disease = "high HDL")
bc_ldl_sig_cp = shared_genes_to_canonical_pathways(permuted_values_path = "shared_genes_to_canonical_pathways/permuted_BC_LDL.RDS", 
                                                   shared_genes = bc_ldl_genes, disease = "high LDL")
bc_prostate_sig_cp = shared_genes_to_canonical_pathways(permuted_values_path = "shared_genes_to_canonical_pathways/permuted_BC_PROSTATE.RDS", 
                                                        shared_genes = bc_prostate_genes, disease = "Prostate cancer")
bc_scz_sig_cp = shared_genes_to_canonical_pathways(permuted_values_path = "shared_genes_to_canonical_pathways/permuted_BC_SCZ.RDS", 
                                                   shared_genes = bc_scz_genes, disease = "Schizophrenia")
bc_t2dm_sig_cp = shared_genes_to_canonical_pathways(permuted_values_path = "shared_genes_to_canonical_pathways/permuted_BC_T2DM_DIAMANTE.RDS", 
                                                    shared_genes = bc_t2dm_genes, disease = "Type 2 diabetes")
# sort based on shared gene
bc_depression_sig_cp = bc_depression_sig_cp %>% arrange(`Shared gene symbol`)
bc_hdl_sig_cp = bc_hdl_sig_cp %>% arrange(`Shared gene symbol`)
bc_ldl_sig_cp = bc_ldl_sig_cp %>% arrange(`Shared gene symbol`)
bc_prostate_sig_cp = bc_prostate_sig_cp %>% arrange(`Shared gene symbol`)
bc_scz_sig_cp = bc_scz_sig_cp %>% arrange(`Shared gene symbol`)
bc_t2dm_sig_cp = bc_t2dm_sig_cp %>% arrange(`Shared gene symbol`)

filt_shared_genes_connected_cps = rbind(bc_depression_sig_cp, bc_hdl_sig_cp, bc_ldl_sig_cp, bc_prostate_sig_cp, bc_scz_sig_cp, bc_t2dm_sig_cp)

write.xlsx(filt_shared_genes_connected_cps, file = "supplementary_tables/S4_prioritized_shared_genes_pathways.xlsx")


## Supplementary Table 5 - Drugs - significant canonical pathways ##
depression_drugs = fread("drugs_to_shared_canonical_pathways/BC_DEPRESSION_drug_targets_perm_pvalues_bonferroni.txt")
hdl_drugs = fread("drugs_to_shared_canonical_pathways/BC_HDL_drug_targets_perm_pvalues_bonferroni.txt")
ldl_drugs = fread("drugs_to_shared_canonical_pathways/BC_LDL_drug_targets_perm_pvalues_bonferroni.txt")
prostate_drugs = fread("drugs_to_shared_canonical_pathways/BC_PROSTATE_drug_targets_perm_pvalues_bonferroni.txt")
scz_drugs = fread("drugs_to_shared_canonical_pathways/BC_SCZ_drug_targets_perm_pvalues_bonferroni.txt")
t2dm_drugs = fread("drugs_to_shared_canonical_pathways/BC_T2DM_drug_targets_perm_pvalues_bonferroni.txt")

## After you get the appropriate license from DrugBank, you can download the necessary file from here:
## - Link between DrugBank IDs and drug names: https://go.drugbank.com/releases/latest#external-links --> "External Drug Links" --> Drug Group "All"
## NOTE: download and place the files in the "raw_data" directory
db_drug = fread("raw_data/drug_links_all.csv") %>% dplyr::select(drug_name = Name, "DrugBank ID") %>% distinct()
bc_inv_ind_drugs = fread("preprocessed_data/bc_indicated_investigated_drugs_of_predisposing_diseases.txt")

drugs_sig_cps = function(drugs_cps_pvalues, disease2, drug_dbids = db_drug) {
  drugs_cps_pvalues$drug = stringr::str_to_title(drugs_cps_pvalues$drug)
  drugs_cps_pvalues = drugs_cps_pvalues %>% arrange(drug)
  temp = data.frame(Drug = unique(drugs_cps_pvalues$drug),
                    Indication = disease2)
  
  for (i in 1:nrow(temp)) {
    drug_temp = temp[i, "Drug"]
    drug_cp_temp = drugs_cps_pvalues %>% filter(drug == drug_temp)
    drug_cp_temp = drug_cp_temp %>% filter(perm_pvalue_adj < 0.05)
    if (nrow(drug_cp_temp) == 0) {
      temp[i, "Connected to shared canonical pathway"] = ""
      temp[i, "Recommended for breast cancer repurposing"] = "No"
    } else {
      temp[i, "Connected to shared canonical pathway"] = paste0(drug_cp_temp$canonical_pathway, collapse = " | ")
      temp[i, "Recommended for breast cancer repurposing"] = "Yes"
    }
    if (length(intersect(bc_inv_ind_drugs$drug_name, drug_temp)) == 0) {
      temp[i, "Current drug development phase for breast cancer"] = "Not investigated/approved"
    } else {
      x = bc_inv_ind_drugs %>% filter(drug_name == drug_temp)
      if (nrow(x) > 0) {
        temp[i, "Current drug development phase for breast cancer"] = x$max_phase
      }
    }
  }
  temp = left_join(temp, drug_dbids, by = c("Drug" = "drug_name"))
  temp = temp %>% dplyr::select("Drug Name" = Drug, `DrugBank ID`, everything())
  
  return(temp)
}

depression_drugs = drugs_sig_cps(drugs_cps_pvalues = depression_drugs, disease2 = "Depression", drug_dbids = db_drug)
hdl_drugs = drugs_sig_cps(drugs_cps_pvalues = hdl_drugs, disease2 = "high HDL", drug_dbids = db_drug)
ldl_drugs = drugs_sig_cps(drugs_cps_pvalues = ldl_drugs, disease2 = "high LDL", drug_dbids = db_drug)
prostate_drugs = drugs_sig_cps(drugs_cps_pvalues = prostate_drugs, disease2 = "Prostate cancer", drug_dbids = db_drug)
scz_drugs = drugs_sig_cps(drugs_cps_pvalues = scz_drugs, disease2 = "Schizophrenia", drug_dbids = db_drug)
t2dm_drugs = drugs_sig_cps(drugs_cps_pvalues = t2dm_drugs, disease2 = "Type 2 diabetes", drug_dbids = db_drug)

drugs_cps_recommended = rbind(depression_drugs, hdl_drugs, ldl_drugs, prostate_drugs, scz_drugs, t2dm_drugs)
drugs_cps_recommended$`Current drug development phase for breast cancer` = ifelse(drugs_cps_recommended$`Current drug development phase for breast cancer` == "4", "Approved", drugs_cps_recommended$`Current drug development phase for breast cancer`)
drugs_cps_recommended$`Current drug development phase for breast cancer` = ifelse(drugs_cps_recommended$`Current drug development phase for breast cancer` == "1", "Phase I", drugs_cps_recommended$`Current drug development phase for breast cancer`)
drugs_cps_recommended$`Current drug development phase for breast cancer` = ifelse(drugs_cps_recommended$`Current drug development phase for breast cancer` == "2", "Phase II", drugs_cps_recommended$`Current drug development phase for breast cancer`)
drugs_cps_recommended$`Current drug development phase for breast cancer` = ifelse(drugs_cps_recommended$`Current drug development phase for breast cancer` == "3", "Phase II", drugs_cps_recommended$`Current drug development phase for breast cancer`)
# fix some drugbank IDs
# drugs_cps_recommended[which(is.na(drugs_cps_recommended$`DrugBank ID`)), ]
drugs_cps_recommended[35, "DrugBank ID"] = "DB11936"
drugs_cps_recommended[41, "DrugBank ID"] = "DB13873"
drugs_cps_recommended[63, "DrugBank ID"] = "DB16778"
drugs_cps_recommended[103, "DrugBank ID"] = "DB00030"

drugs_cps_recommended = drugs_cps_recommended %>% dplyr::select(`Drug Name`, `DrugBank ID`, Indication, `Recommended for breast cancer repurposing`, `Current drug development phase for breast cancer`, `Connected to shared canonical pathway`)
write.xlsx(drugs_cps_recommended, file = "supplementary_tables/supplementary_tables/S5_drugs_to_shared_genetics.xlsx")
