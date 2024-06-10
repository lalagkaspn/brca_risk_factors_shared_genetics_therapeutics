
### Connecting drug targets to shared canonical pathways, for each disease pair ###

library(dplyr)
library(data.table)
library(GSA)
library(igraph)

## all canonical pathways with their genes
msigdb_cp = readRDS("preprocessed_data/c2.cp.v2023.2.Hs.entrez_STRING_filtered.RDS")

## STRING network
string_high_conf = fread("preprocessed_data/string_high_conf_processed.txt")
string_genes = unique(c(string_high_conf$from, string_high_conf$to))
g = graph_from_data_frame(string_high_conf, directed = FALSE)

### OBSERVED PPR SCORES ###

## BC-DEPRESSION
# significant canonical pathways
bc_depression_sig_cp = fread("shared_genes_to_canonical_pathways/BC_DEPRESSION_CPs_pvalues_bonferroni.txt") %>% dplyr::select(-sig) %>% filter(pvalue_adj < 0.05) %>% distinct()
# drug targets
depression_drug_targets = fread("preprocessed_data/depression_drug_targets.txt") %>%
  filter(entrezID %in% string_genes)
unique_depression_drug_targets = unique(depression_drug_targets$entrezID)
# observed average PPR score of genes in significant canonical pathways when using each drug target as seed gene
sig_cp_depression_drug_targets_observed = vector("list", length = length(unique_depression_drug_targets))
names(sig_cp_depression_drug_targets_observed) = unique_depression_drug_targets
for (i in 1:length(sig_cp_depression_drug_targets_observed)) {
  # drug target
  drug_target = names(sig_cp_depression_drug_targets_observed)[i]
  # dataframe to populate later
  ppr_avg_scores = data.frame(canonical_pathway = bc_depression_sig_cp$canonical_pathway,
                              average_ppr_score = 0)
  for (z in 1:nrow(bc_depression_sig_cp)) {
    # canonical pathway
    cp = bc_depression_sig_cp[z, canonical_pathway]
    # genes in canonical pathway
    cp_genes = msigdb_cp[[cp]]
    # PPR personalization vector
    p = c(rep(0, length(V(g))))
    p[which(V(g)$name == drug_target)] = 1
    # PPR
    pr_result = page_rank(g, algo = "prpack", damping = 0.85, personalized = p, weights = NA)
    pr_result = pr_result$vector
    pr_result = as.data.frame(pr_result)
    pr_result$entrezID = rownames(pr_result)
    rownames(pr_result) = NULL
    pr_result = pr_result %>%
      dplyr::select(entrezID, propagation_score = pr_result) %>%
      filter(entrezID %in% cp_genes) %>%
      mutate(avg_propagation_score = mean(propagation_score)) %>%
      dplyr::select(avg_propagation_score) %>%
      distinct()
    # populate dataframe
    ppr_avg_scores[z, "average_ppr_score"] = pr_result$avg_propagation_score
  }
  sig_cp_depression_drug_targets_observed[[i]] = ppr_avg_scores
  
  # track progress
  cat(i, "-", length(sig_cp_depression_drug_targets_observed), "\n")
}
rm(i, drug_target, z, cp, cp_genes, p, pr_result, ppr_avg_scores)
saveRDS(sig_cp_depression_drug_targets_observed, "drugs_to_shared_canonical_pathways/BC_DEPRESSION_drug_targets_sigCP_avg_ppr_scores.RDS")

## BC-HDL
# significant canonical pathways
bc_hdl_sig_cp = fread("shared_genes_to_canonical_pathways/BC_HDL_CPs_pvalues_bonferroni.txt") %>% dplyr::select(-sig) %>% filter(pvalue_adj < 0.05) %>% distinct()
# drug targets
hdl_drug_targets = fread("preprocessed_data/hdl_drug_targets.txt") %>%
  filter(entrezID %in% string_genes)
unique_hdl_drug_targets = unique(hdl_drug_targets$entrezID)
# observed average PPR score of genes in significant canonical pathways when using each drug target as seed gene
sig_cp_hdl_drug_targets_observed = vector("list", length = length(unique_hdl_drug_targets))
names(sig_cp_hdl_drug_targets_observed) = unique_hdl_drug_targets
for (i in 1:length(sig_cp_hdl_drug_targets_observed)) {
  # drug target
  drug_target = names(sig_cp_hdl_drug_targets_observed)[i]
  # dataframe to populate later
  ppr_avg_scores = data.frame(canonical_pathway = bc_hdl_sig_cp$canonical_pathway,
                              average_ppr_score = 0)
  for (z in 1:nrow(bc_hdl_sig_cp)) {
    # canonical pathway
    cp = bc_hdl_sig_cp[z, canonical_pathway]
    # genes in canonical pathway
    cp_genes = msigdb_cp[[cp]]
    # PPR personalization vector
    p = c(rep(0, length(V(g))))
    p[which(V(g)$name == drug_target)] = 1
    # PPR
    pr_result = page_rank(g, algo = "prpack", damping = 0.85, personalized = p, weights = NA)
    pr_result = pr_result$vector
    pr_result = as.data.frame(pr_result)
    pr_result$entrezID = rownames(pr_result)
    rownames(pr_result) = NULL
    pr_result = pr_result %>%
      dplyr::select(entrezID, propagation_score = pr_result) %>%
      filter(entrezID %in% cp_genes) %>%
      mutate(avg_propagation_score = mean(propagation_score)) %>%
      dplyr::select(avg_propagation_score) %>%
      distinct()
    # populate dataframe
    ppr_avg_scores[z, "average_ppr_score"] = pr_result$avg_propagation_score
  }
  sig_cp_hdl_drug_targets_observed[[i]] = ppr_avg_scores
  
  # track progress
  cat(i, "-", length(sig_cp_hdl_drug_targets_observed), "\n")
}
rm(i, drug_target, z, cp, cp_genes, p, pr_result, ppr_avg_scores)
saveRDS(sig_cp_hdl_drug_targets_observed, "drugs_to_shared_canonical_pathways/BC_HDL_drug_targets_sigCP_avg_ppr_scores.RDS")

## BC-LDL
# significant canonical pathways
bc_ldl_sig_cp = fread("shared_genes_to_canonical_pathways/BC_LDL_CPs_pvalues_bonferroni.txt") %>% dplyr::select(-sig) %>% filter(pvalue_adj < 0.05) %>% distinct()
# drug targets
ldl_drug_targets = fread("preprocessed_data/ldl_drug_targets.txt") %>%
  filter(entrezID %in% string_genes)
unique_ldl_drug_targets = unique(ldl_drug_targets$entrezID)
# observed average PPR score of genes in significant canonical pathways when using each drug target as seed gene
sig_cp_ldl_drug_targets_observed = vector("list", length = length(unique_ldl_drug_targets))
names(sig_cp_ldl_drug_targets_observed) = unique_ldl_drug_targets
for (i in 1:length(sig_cp_ldl_drug_targets_observed)) {
  # drug target
  drug_target = names(sig_cp_ldl_drug_targets_observed)[i]
  # dataframe to populate later
  ppr_avg_scores = data.frame(canonical_pathway = bc_ldl_sig_cp$canonical_pathway,
                              average_ppr_score = 0)
  for (z in 1:nrow(bc_ldl_sig_cp)) {
    # canonical pathway
    cp = bc_ldl_sig_cp[z, canonical_pathway]
    # genes in canonical pathway
    cp_genes = msigdb_cp[[cp]]
    # PPR personalization vector
    p = c(rep(0, length(V(g))))
    p[which(V(g)$name == drug_target)] = 1
    # PPR
    pr_result = page_rank(g, algo = "prpack", damping = 0.85, personalized = p, weights = NA)
    pr_result = pr_result$vector
    pr_result = as.data.frame(pr_result)
    pr_result$entrezID = rownames(pr_result)
    rownames(pr_result) = NULL
    pr_result = pr_result %>%
      dplyr::select(entrezID, propagation_score = pr_result) %>%
      filter(entrezID %in% cp_genes) %>%
      mutate(avg_propagation_score = mean(propagation_score)) %>%
      dplyr::select(avg_propagation_score) %>%
      distinct()
    # populate dataframe
    ppr_avg_scores[z, "average_ppr_score"] = pr_result$avg_propagation_score
  }
  sig_cp_ldl_drug_targets_observed[[i]] = ppr_avg_scores
  
  # track progress
  cat(i, "-", length(sig_cp_ldl_drug_targets_observed), "\n")
}
rm(i, drug_target, z, cp, cp_genes, p, pr_result, ppr_avg_scores)
saveRDS(sig_cp_ldl_drug_targets_observed, "drugs_to_shared_canonical_pathways/BC_LDL_drug_targets_sigCP_avg_ppr_scores.RDS")

## BC-PROSTATE
# significant canonical pathways
bc_prostate_sig_cp = fread("shared_genes_to_canonical_pathways/BC_PROSTATE_CPs_pvalues_bonferroni.txt") %>% dplyr::select(-sig) %>% filter(pvalue_adj < 0.05) %>% distinct()
# drug targets
prostate_drug_targets = fread("preprocessed_data/prostate_drug_targets.txt") %>%
  filter(entrezID %in% string_genes)
unique_prostate_drug_targets = unique(prostate_drug_targets$entrezID)
# observed average PPR score of genes in significant canonical pathways when using each drug target as seed gene
sig_cp_prostate_drug_targets_observed = vector("list", length = length(unique_prostate_drug_targets))
names(sig_cp_prostate_drug_targets_observed) = unique_prostate_drug_targets
for (i in 1:length(sig_cp_prostate_drug_targets_observed)) {
  # drug target
  drug_target = names(sig_cp_prostate_drug_targets_observed)[i]
  # dataframe to populate later
  ppr_avg_scores = data.frame(canonical_pathway = bc_prostate_sig_cp$canonical_pathway,
                              average_ppr_score = 0)
  for (z in 1:nrow(bc_prostate_sig_cp)) {
    # canonical pathway
    cp = bc_prostate_sig_cp[z, canonical_pathway]
    # genes in canonical pathway
    cp_genes = msigdb_cp[[cp]]
    # PPR personalization vector
    p = c(rep(0, length(V(g))))
    p[which(V(g)$name == drug_target)] = 1
    # PPR
    pr_result = page_rank(g, algo = "prpack", damping = 0.85, personalized = p, weights = NA)
    pr_result = pr_result$vector
    pr_result = as.data.frame(pr_result)
    pr_result$entrezID = rownames(pr_result)
    rownames(pr_result) = NULL
    pr_result = pr_result %>%
      dplyr::select(entrezID, propagation_score = pr_result) %>%
      filter(entrezID %in% cp_genes) %>%
      mutate(avg_propagation_score = mean(propagation_score)) %>%
      dplyr::select(avg_propagation_score) %>%
      distinct()
    # populate dataframe
    ppr_avg_scores[z, "average_ppr_score"] = pr_result$avg_propagation_score
  }
  sig_cp_prostate_drug_targets_observed[[i]] = ppr_avg_scores
  
  # track progress
  cat(i, "-", length(sig_cp_prostate_drug_targets_observed), "\n")
}
rm(i, drug_target, z, cp, cp_genes, p, pr_result, ppr_avg_scores)
saveRDS(sig_cp_prostate_drug_targets_observed, "drugs_to_shared_canonical_pathways/BC_PROSTATE_drug_targets_sigCP_avg_ppr_scores.RDS")

## BC-SCZ
# significant canonical pathways
bc_scz_sig_cp = fread("shared_genes_to_canonical_pathways/BC_SCZ_CPs_pvalues_bonferroni.txt") %>% dplyr::select(-sig) %>% filter(pvalue_adj < 0.05) %>% distinct()
# drug targets
scz_drug_targets = fread("preprocessed_data/schizophrenia_drug_targets.txt") %>%
  filter(entrezID %in% string_genes)
unique_scz_drug_targets = unique(scz_drug_targets$entrezID)
# observed average PPR score of genes in significant canonical pathways when using each drug target as seed gene
sig_cp_scz_drug_targets_observed = vector("list", length = length(unique_scz_drug_targets))
names(sig_cp_scz_drug_targets_observed) = unique_scz_drug_targets
for (i in 1:length(sig_cp_scz_drug_targets_observed)) {
  # drug target
  drug_target = names(sig_cp_scz_drug_targets_observed)[i]
  # dataframe to populate later
  ppr_avg_scores = data.frame(canonical_pathway = bc_scz_sig_cp$canonical_pathway,
                              average_ppr_score = 0)
  for (z in 1:nrow(bc_scz_sig_cp)) {
    # canonical pathway
    cp = bc_scz_sig_cp[z, canonical_pathway]
    # genes in canonical pathway
    cp_genes = msigdb_cp[[cp]]
    # PPR personalization vector
    p = c(rep(0, length(V(g))))
    p[which(V(g)$name == drug_target)] = 1
    # PPR
    pr_result = page_rank(g, algo = "prpack", damping = 0.85, personalized = p, weights = NA)
    pr_result = pr_result$vector
    pr_result = as.data.frame(pr_result)
    pr_result$entrezID = rownames(pr_result)
    rownames(pr_result) = NULL
    pr_result = pr_result %>%
      dplyr::select(entrezID, propagation_score = pr_result) %>%
      filter(entrezID %in% cp_genes) %>%
      mutate(avg_propagation_score = mean(propagation_score)) %>%
      dplyr::select(avg_propagation_score) %>%
      distinct()
    # populate dataframe
    ppr_avg_scores[z, "average_ppr_score"] = pr_result$avg_propagation_score
  }
  sig_cp_scz_drug_targets_observed[[i]] = ppr_avg_scores
  
  # track progress
  cat(i, "-", length(sig_cp_scz_drug_targets_observed), "\n")
}
rm(i, drug_target, z, cp, cp_genes, p, pr_result, ppr_avg_scores)
saveRDS(sig_cp_scz_drug_targets_observed, "drugs_to_shared_canonical_pathways/BC_SCZ_drug_targets_sigCP_avg_ppr_scores.RDS")

## BC-T2DM_DIAMANTE
# significant canonical pathways
bc_t2dm_sig_cp = fread("shared_genes_to_canonical_pathways/BC_T2DM_DIAMANTE_CPs_pvalues_bonferroni.txt") %>% dplyr::select(-sig) %>% filter(pvalue_adj < 0.05) %>% distinct()
# drug targets
t2dm_drug_targets = fread("preprocessed_data/t2dm_drug_targets.txt") %>%
  filter(entrezID %in% string_genes)
unique_t2dm_drug_targets = unique(t2dm_drug_targets$entrezID)
# observed average PPR score of genes in significant canonical pathways when using each drug target as seed gene
sig_cp_t2dm_drug_targets_observed = vector("list", length = length(unique_t2dm_drug_targets))
names(sig_cp_t2dm_drug_targets_observed) = unique_t2dm_drug_targets
for (i in 1:length(sig_cp_t2dm_drug_targets_observed)) {
  # drug target
  drug_target = names(sig_cp_t2dm_drug_targets_observed)[i]
  # dataframe to populate later
  ppr_avg_scores = data.frame(canonical_pathway = bc_t2dm_sig_cp$canonical_pathway,
                              average_ppr_score = 0)
  for (z in 1:nrow(bc_t2dm_sig_cp)) {
    # canonical pathway
    cp = bc_t2dm_sig_cp[z, canonical_pathway]
    # genes in canonical pathway
    cp_genes = msigdb_cp[[cp]]
    # PPR personalization vector
    p = c(rep(0, length(V(g))))
    p[which(V(g)$name == drug_target)] = 1
    # PPR
    pr_result = page_rank(g, algo = "prpack", damping = 0.85, personalized = p, weights = NA)
    pr_result = pr_result$vector
    pr_result = as.data.frame(pr_result)
    pr_result$entrezID = rownames(pr_result)
    rownames(pr_result) = NULL
    pr_result = pr_result %>%
      dplyr::select(entrezID, propagation_score = pr_result) %>%
      filter(entrezID %in% cp_genes) %>%
      mutate(avg_propagation_score = mean(propagation_score)) %>%
      dplyr::select(avg_propagation_score) %>%
      distinct()
    # populate dataframe
    ppr_avg_scores[z, "average_ppr_score"] = pr_result$avg_propagation_score
  }
  sig_cp_t2dm_drug_targets_observed[[i]] = ppr_avg_scores
  
  # track progress
  cat(i, "-", length(sig_cp_t2dm_drug_targets_observed), "\n")
}
rm(i, drug_target, z, cp, cp_genes, p, pr_result, ppr_avg_scores)
saveRDS(sig_cp_t2dm_drug_targets_observed, "drugs_to_shared_canonical_pathways/BC_T2DM_DIAMANTE_drug_targets_sigCP_avg_ppr_scores.RDS")

rm(list = ls())
gc()

### PERMUTATIONS ###

### Connect drug targets to the canonical pathways that are significantly connected to the shared genes for each
### BC and BC-related disease pair
## permutations

## STRING network
string_high_conf = fread("preprocessed_data/string_high_conf_processed.txt")
string_genes = unique(c(string_high_conf$from, string_high_conf$to))
g = graph_from_data_frame(string_high_conf, directed = FALSE)
node_degrees = data.frame(degree = degree(g))
node_degrees$entrezID = rownames(node_degrees) ; rownames(node_degrees) = NULL
# degree bins for permutations
bins = list(first = seq(1,4,1), second = seq(5,13,1), third = seq(14,36, 1), fourth = seq(37,729, 1))

## average PPR scores of canonical pathways to all genes in STRING
# This file was created by running Personalized PageRank using a gene as seed gene and genes in one canonical pathway as target --> to speed up permutation analysis
all_genes_cp_avg_ppr_scores = fread("/project/pi_rachel_melamed_uml_edu/Panos/GWAS/msigdb_canonical_pathways/all_genes_cp_ppr_scores.txt")

## BC-DEPRESSION
bc_depression_observed = readRDS("drugs_to_shared_canonical_pathways/BC_DEPRESSION_drug_targets_sigCP_avg_ppr_scores.RDS")
bc_depression_sig_cps = fread("shared_genes_to_canonical_pathways/BC_DEPRESSION_CPs_pvalues_bonferroni.txt") %>%
  filter(pvalue_adj < 0.05) %>% 
  dplyr::select(canonical_pathway) %>%
  distinct()
bc_depression_permuted = matrix(ncol = nrow(bc_depression_sig_cps),
                                nrow = length(bc_depression_observed))
colnames(bc_depression_permuted) = bc_depression_sig_cps$canonical_pathway
rownames(bc_depression_permuted) = names(bc_depression_observed)
bc_depression_permuted = as.data.frame(bc_depression_permuted)

for (i in 1:nrow(bc_depression_permuted)) {
  # drug gene-target
  gene_temp = rownames(bc_depression_permuted)[i]
  # degree of the drug gene-target
  degree_gene_temp = node_degrees %>% filter(entrezID == gene_temp) ; degree_gene_temp = degree_gene_temp$degree
  # degree bin that the drug gene-target belongs to
  degree_bin = ifelse(degree_gene_temp %in% bins$first, "first", 
                      ifelse(degree_gene_temp %in% bins$second, "second",
                             ifelse(degree_gene_temp %in% bins$third, "third", "fourth")))
  # all genes in this degree bin
  all_genes_in_degree_bin = node_degrees %>% filter(degree %in% bins[[degree_bin]])
  # choose 1,00 random genes from the above list
  random_genes_in_degree_bin = sample(all_genes_in_degree_bin$entrezID, size = 1000, replace = FALSE)
  # average PPR score of the canonical pathway when using the above genes as seed genes (one at a time)
  random_genes_cp_avg_ppr_scores = all_genes_cp_avg_ppr_scores[which(all_genes_cp_avg_ppr_scores$random_gene %in% random_genes_in_degree_bin), ]
  for (z in 1:ncol(bc_depression_permuted)) {
    random_genes_cp_avg_ppr_scores_temp = random_genes_cp_avg_ppr_scores %>% filter(canonical_pathway == colnames(bc_depression_permuted)[[z]])
    observed_cp_avg_ppr_scores = bc_depression_observed[[gene_temp]][which(bc_depression_observed[[gene_temp]]$canonical_pathway == colnames(bc_depression_permuted)[[z]]), "average_ppr_score"]
    # permuted p-value
    bc_depression_permuted[i, z] = sum(observed_cp_avg_ppr_scores <= random_genes_cp_avg_ppr_scores_temp$average_ppr_score) / 1000
  }
  # track progress
  cat(i, "\n")
} ; rm(i, z, gene_temp, degree_gene_temp, degree_bin, all_genes_in_degree_bin, random_genes_in_degree_bin, random_genes_cp_avg_ppr_scores, random_genes_cp_avg_ppr_scores_temp, observed_cp_avg_ppr_scores)
bc_depression_permuted$drug_target = rownames(bc_depression_permuted) ; rownames(bc_depression_permuted) = NULL
bc_depression_permuted = bc_depression_permuted %>% dplyr::select(drug_target, everything())
fwrite(bc_depression_permuted, "drugs_to_shared_canonical_pathways/BC_DEPRESSION_drug_targets_perm_pvalues_unadjusted.txt", sep = "\t", row.names = FALSE)
rm(bc_depression_observed, bc_depression_permuted, bc_depression_sig_cps)

## BC-HDL
bc_hdl_observed = readRDS("drugs_to_shared_canonical_pathways/BC_HDL_drug_targets_sigCP_avg_ppr_scores.RDS")
bc_hdl_sig_cps = fread("shared_genes_to_canonical_pathways/BC_HDL_CPs_pvalues_bonferroni.txt") %>%
  filter(pvalue_adj < 0.05) %>% 
  dplyr::select(canonical_pathway) %>%
  distinct()
bc_hdl_permuted = matrix(ncol = nrow(bc_hdl_sig_cps),
                         nrow = length(bc_hdl_observed))
colnames(bc_hdl_permuted) = bc_hdl_sig_cps$canonical_pathway
rownames(bc_hdl_permuted) = names(bc_hdl_observed)
bc_hdl_permuted = as.data.frame(bc_hdl_permuted)
for (i in 1:nrow(bc_hdl_permuted)) {
  # drug gene-target
  gene_temp = rownames(bc_hdl_permuted)[i]
  # degree of the drug gene-target
  degree_gene_temp = node_degrees %>% filter(entrezID == gene_temp) ; degree_gene_temp = degree_gene_temp$degree
  # degree bin that the drug gene-target belongs to
  degree_bin = ifelse(degree_gene_temp %in% bins$first, "first", 
                      ifelse(degree_gene_temp %in% bins$second, "second",
                             ifelse(degree_gene_temp %in% bins$third, "third", "fourth")))
  # all genes in this degree bin
  all_genes_in_degree_bin = node_degrees %>% filter(degree %in% bins[[degree_bin]])
  # choose 1,00 random genes from the above list
  random_genes_in_degree_bin = sample(all_genes_in_degree_bin$entrezID, size = 1000, replace = FALSE)
  # average PPR score of the canonical pathway when using the above genes as seed genes (one at a time)
  random_genes_cp_avg_ppr_scores = all_genes_cp_avg_ppr_scores[which(all_genes_cp_avg_ppr_scores$random_gene %in% random_genes_in_degree_bin), ]
  for (z in 1:ncol(bc_hdl_permuted)) {
    random_genes_cp_avg_ppr_scores_temp = random_genes_cp_avg_ppr_scores %>% filter(canonical_pathway == colnames(bc_hdl_permuted)[[z]])
    observed_cp_avg_ppr_scores = bc_hdl_observed[[gene_temp]][which(bc_hdl_observed[[gene_temp]]$canonical_pathway == colnames(bc_hdl_permuted)[[z]]), "average_ppr_score"]
    # permuted p-value
    bc_hdl_permuted[i, z] = sum(observed_cp_avg_ppr_scores <= random_genes_cp_avg_ppr_scores_temp$average_ppr_score) / 1000
  }
  # track progress
  cat(i, "\n")
} ; rm(i, z, gene_temp, degree_gene_temp, degree_bin, all_genes_in_degree_bin, random_genes_in_degree_bin, random_genes_cp_avg_ppr_scores, random_genes_cp_avg_ppr_scores_temp, observed_cp_avg_ppr_scores)
bc_hdl_permuted$drug_target = rownames(bc_hdl_permuted) ; rownames(bc_hdl_permuted) = NULL
bc_hdl_permuted = bc_hdl_permuted %>% dplyr::select(drug_target, everything())
fwrite(bc_hdl_permuted, "drugs_to_shared_canonical_pathways/BC_HDL_drug_targets_perm_pvalues_unadjusted.txt", sep = "\t", row.names = FALSE)
rm(bc_hdl_observed, bc_hdl_permuted, bc_hdl_sig_cps)

## BC-LDL
bc_ldl_observed = readRDS("drugs_to_shared_canonical_pathways/BC_LDL_drug_targets_sigCP_avg_ppr_scores.RDS")
bc_ldl_sig_cps = fread("shared_genes_to_canonical_pathways/BC_LDL_CPs_pvalues_bonferroni.txt") %>%
  filter(pvalue_adj < 0.05) %>%
  dplyr::select(canonical_pathway) %>%
  distinct()
bc_ldl_permuted = matrix(ncol = nrow(bc_ldl_sig_cps),
                         nrow = length(bc_ldl_observed))
colnames(bc_ldl_permuted) = bc_ldl_sig_cps$canonical_pathway
rownames(bc_ldl_permuted) = names(bc_ldl_observed)
bc_ldl_permuted = as.data.frame(bc_ldl_permuted)
for (i in 1:nrow(bc_ldl_permuted)) {
  # drug gene-target
  gene_temp = rownames(bc_ldl_permuted)[i]
  # degree of the drug gene-target
  degree_gene_temp = node_degrees %>% filter(entrezID == gene_temp) ; degree_gene_temp = degree_gene_temp$degree
  # degree bin that the drug gene-target belongs to
  degree_bin = ifelse(degree_gene_temp %in% bins$first, "first", 
                      ifelse(degree_gene_temp %in% bins$second, "second",
                             ifelse(degree_gene_temp %in% bins$third, "third", "fourth")))
  # all genes in this degree bin
  all_genes_in_degree_bin = node_degrees %>% filter(degree %in% bins[[degree_bin]])
  # choose 1,00 random genes from the above list
  random_genes_in_degree_bin = sample(all_genes_in_degree_bin$entrezID, size = 1000, replace = FALSE)
  # average PPR score of the canonical pathway when using the above genes as seed genes (one at a time)
  random_genes_cp_avg_ppr_scores = all_genes_cp_avg_ppr_scores[which(all_genes_cp_avg_ppr_scores$random_gene %in% random_genes_in_degree_bin), ]
  for (z in 1:ncol(bc_ldl_permuted)) {
    random_genes_cp_avg_ppr_scores_temp = random_genes_cp_avg_ppr_scores %>% filter(canonical_pathway == colnames(bc_ldl_permuted)[[z]])
    observed_cp_avg_ppr_scores = bc_ldl_observed[[gene_temp]][which(bc_ldl_observed[[gene_temp]]$canonical_pathway == colnames(bc_ldl_permuted)[[z]]), "average_ppr_score"]
    # permuted p-value
    bc_ldl_permuted[i, z] = sum(observed_cp_avg_ppr_scores <= random_genes_cp_avg_ppr_scores_temp$average_ppr_score) / 1000
  }
  # track progress
  cat(i, "\n")
} ; rm(i, z, gene_temp, degree_gene_temp, degree_bin, all_genes_in_degree_bin, random_genes_in_degree_bin, random_genes_cp_avg_ppr_scores, random_genes_cp_avg_ppr_scores_temp, observed_cp_avg_ppr_scores)
bc_ldl_permuted$drug_target = rownames(bc_ldl_permuted) ; rownames(bc_ldl_permuted) = NULL
bc_ldl_permuted = bc_ldl_permuted %>% dplyr::select(drug_target, everything())
fwrite(bc_ldl_permuted, "drugs_to_shared_canonical_pathways/BC_LDL_drug_targets_perm_pvalues_unadjusted.txt", sep = "\t", row.names = FALSE)
rm(bc_ldl_observed, bc_ldl_permuted, bc_ldl_sig_cps)

## BC-PROSTATE
bc_prostate_observed = readRDS("drugs_to_shared_canonical_pathways//BC_PROSTATE_drug_targets_sigCP_avg_ppr_scores.RDS")
bc_prostate_sig_cps = fread("shared_genes_to_canonical_pathways/BC_PROSTATE_CPs_pvalues_bonferroni.txt") %>%
  filter(pvalue_adj < 0.05) %>%
  dplyr::select(canonical_pathway) %>%
  distinct()
bc_prostate_permuted = matrix(ncol = nrow(bc_prostate_sig_cps),
                              nrow = length(bc_prostate_observed))
colnames(bc_prostate_permuted) = bc_prostate_sig_cps$canonical_pathway
rownames(bc_prostate_permuted) = names(bc_prostate_observed)
bc_prostate_permuted = as.data.frame(bc_prostate_permuted)
for (i in 1:nrow(bc_prostate_permuted)) {
  # drug gene-target
  gene_temp = rownames(bc_prostate_permuted)[i]
  # degree of the drug gene-target
  degree_gene_temp = node_degrees %>% filter(entrezID == gene_temp) ; degree_gene_temp = degree_gene_temp$degree
  # degree bin that the drug gene-target belongs to
  degree_bin = ifelse(degree_gene_temp %in% bins$first, "first", 
                      ifelse(degree_gene_temp %in% bins$second, "second",
                             ifelse(degree_gene_temp %in% bins$third, "third", "fourth")))
  # all genes in this degree bin
  all_genes_in_degree_bin = node_degrees %>% filter(degree %in% bins[[degree_bin]])
  # choose 1,00 random genes from the above list
  random_genes_in_degree_bin = sample(all_genes_in_degree_bin$entrezID, size = 1000, replace = FALSE)
  # average PPR score of the canonical pathway when using the above genes as seed genes (one at a time)
  random_genes_cp_avg_ppr_scores = all_genes_cp_avg_ppr_scores[which(all_genes_cp_avg_ppr_scores$random_gene %in% random_genes_in_degree_bin), ]
  for (z in 1:ncol(bc_prostate_permuted)) {
    random_genes_cp_avg_ppr_scores_temp = random_genes_cp_avg_ppr_scores %>% filter(canonical_pathway == colnames(bc_prostate_permuted)[[z]])
    observed_cp_avg_ppr_scores = bc_prostate_observed[[gene_temp]][which(bc_prostate_observed[[gene_temp]]$canonical_pathway == colnames(bc_prostate_permuted)[[z]]), "average_ppr_score"]
    # permuted p-value
    bc_prostate_permuted[i, z] = sum(observed_cp_avg_ppr_scores <= random_genes_cp_avg_ppr_scores_temp$average_ppr_score) / 1000
  }
  # track progress
  cat(i, "\n")
} ; rm(i, z, gene_temp, degree_gene_temp, degree_bin, all_genes_in_degree_bin, random_genes_in_degree_bin, random_genes_cp_avg_ppr_scores, random_genes_cp_avg_ppr_scores_temp, observed_cp_avg_ppr_scores)
bc_prostate_permuted$drug_target = rownames(bc_prostate_permuted) ; rownames(bc_prostate_permuted) = NULL
bc_prostate_permuted = bc_prostate_permuted %>% dplyr::select(drug_target, everything())
fwrite(bc_prostate_permuted, "drugs_to_shared_canonical_pathways/BC_PROSTATE_drug_targets_perm_pvalues_unadjusted.txt", sep = "\t", row.names = FALSE)
rm(bc_prostate_observed, bc_prostate_permuted, bc_prostate_sig_cps)

## BC-SCZ
bc_scz_observed = readRDS("drugs_to_shared_canonical_pathways/BC_SCZ_drug_targets_sigCP_avg_ppr_scores.RDS")
bc_scz_sig_cps = fread("shared_genes_to_canonical_pathways/BC_SCZ_CPs_pvalues_bonferroni.txt") %>%
  filter(pvalue_adj < 0.05) %>% 
  dplyr::select(canonical_pathway) %>%
  distinct()
bc_scz_permuted = matrix(ncol = nrow(bc_scz_sig_cps),
                         nrow = length(bc_scz_observed))
colnames(bc_scz_permuted) = bc_scz_sig_cps$canonical_pathway
rownames(bc_scz_permuted) = names(bc_scz_observed)
bc_scz_permuted = as.data.frame(bc_scz_permuted)
for (i in 1:nrow(bc_scz_permuted)) {
  # drug gene-target
  gene_temp = rownames(bc_scz_permuted)[i]
  # degree of the drug gene-target
  degree_gene_temp = node_degrees %>% filter(entrezID == gene_temp) ; degree_gene_temp = degree_gene_temp$degree
  # degree bin that the drug gene-target belongs to
  degree_bin = ifelse(degree_gene_temp %in% bins$first, "first", 
                      ifelse(degree_gene_temp %in% bins$second, "second",
                             ifelse(degree_gene_temp %in% bins$third, "third", "fourth")))
  # all genes in this degree bin
  all_genes_in_degree_bin = node_degrees %>% filter(degree %in% bins[[degree_bin]])
  # choose 1,00 random genes from the above list
  random_genes_in_degree_bin = sample(all_genes_in_degree_bin$entrezID, size = 1000, replace = FALSE)
  # average PPR score of the canonical pathway when using the above genes as seed genes (one at a time)
  random_genes_cp_avg_ppr_scores = all_genes_cp_avg_ppr_scores[which(all_genes_cp_avg_ppr_scores$random_gene %in% random_genes_in_degree_bin), ]
  for (z in 1:ncol(bc_scz_permuted)) {
    random_genes_cp_avg_ppr_scores_temp = random_genes_cp_avg_ppr_scores %>% filter(canonical_pathway == colnames(bc_scz_permuted)[[z]])
    observed_cp_avg_ppr_scores = bc_scz_observed[[gene_temp]][which(bc_scz_observed[[gene_temp]]$canonical_pathway == colnames(bc_scz_permuted)[[z]]), "average_ppr_score"]
    # permuted p-value
    bc_scz_permuted[i, z] = sum(observed_cp_avg_ppr_scores <= random_genes_cp_avg_ppr_scores_temp$average_ppr_score) / 1000
  }
  # track progress
  cat(i, "\n")
} ; rm(i, z, gene_temp, degree_gene_temp, degree_bin, all_genes_in_degree_bin, random_genes_in_degree_bin, random_genes_cp_avg_ppr_scores, random_genes_cp_avg_ppr_scores_temp, observed_cp_avg_ppr_scores)
bc_scz_permuted$drug_target = rownames(bc_scz_permuted) ; rownames(bc_scz_permuted) = NULL
bc_scz_permuted = bc_scz_permuted %>% dplyr::select(drug_target, everything())
fwrite(bc_scz_permuted, "drugs_to_shared_canonical_pathways/BC_SCZ_drug_targets_perm_pvalues_unadjusted.txt", sep = "\t", row.names = FALSE)
rm(bc_scz_observed, bc_scz_permuted, bc_scz_sig_cps)

## BC-T2DM_DIAMANTE
bc_t2dm_observed = readRDS("drugs_to_shared_canonical_pathways/BC_T2DM_DIAMANTE_drug_targets_sigCP_avg_ppr_scores.RDS")
bc_t2dm_sig_cps = fread("shared_genes_to_canonical_pathways/BC_T2DM_DIAMANTE_CPs_pvalues_bonferroni.txt") %>%
  filter(pvalue_adj < 0.05) %>% 
  dplyr::select(canonical_pathway) %>%
  distinct()
bc_t2dm_permuted = matrix(ncol = nrow(bc_t2dm_sig_cps),
                          nrow = length(bc_t2dm_observed))
colnames(bc_t2dm_permuted) = bc_t2dm_sig_cps$canonical_pathway
rownames(bc_t2dm_permuted) = names(bc_t2dm_observed)
bc_t2dm_permuted = as.data.frame(bc_t2dm_permuted)
for (i in 1:nrow(bc_t2dm_permuted)) {
  # drug gene-target
  gene_temp = rownames(bc_t2dm_permuted)[i]
  # degree of the drug gene-target
  degree_gene_temp = node_degrees %>% filter(entrezID == gene_temp) ; degree_gene_temp = degree_gene_temp$degree
  # degree bin that the drug gene-target belongs to
  degree_bin = ifelse(degree_gene_temp %in% bins$first, "first", 
                      ifelse(degree_gene_temp %in% bins$second, "second",
                             ifelse(degree_gene_temp %in% bins$third, "third", "fourth")))
  # all genes in this degree bin
  all_genes_in_degree_bin = node_degrees %>% filter(degree %in% bins[[degree_bin]])
  # choose 1,00 random genes from the above list
  random_genes_in_degree_bin = sample(all_genes_in_degree_bin$entrezID, size = 1000, replace = FALSE)
  # average PPR score of the canonical pathway when using the above genes as seed genes (one at a time)
  random_genes_cp_avg_ppr_scores = all_genes_cp_avg_ppr_scores[which(all_genes_cp_avg_ppr_scores$random_gene %in% random_genes_in_degree_bin), ]
  for (z in 1:ncol(bc_t2dm_permuted)) {
    random_genes_cp_avg_ppr_scores_temp = random_genes_cp_avg_ppr_scores %>% filter(canonical_pathway == colnames(bc_t2dm_permuted)[[z]])
    observed_cp_avg_ppr_scores = bc_t2dm_observed[[gene_temp]][which(bc_t2dm_observed[[gene_temp]]$canonical_pathway == colnames(bc_t2dm_permuted)[[z]]), "average_ppr_score"]
    # permuted p-value
    bc_t2dm_permuted[i, z] = sum(observed_cp_avg_ppr_scores <= random_genes_cp_avg_ppr_scores_temp$average_ppr_score) / 1000
  }
  # track progress
  cat(i, "\n")
} ; rm(i, z, gene_temp, degree_gene_temp, degree_bin, all_genes_in_degree_bin, random_genes_in_degree_bin, random_genes_cp_avg_ppr_scores, random_genes_cp_avg_ppr_scores_temp, observed_cp_avg_ppr_scores)
bc_t2dm_permuted$drug_target = rownames(bc_t2dm_permuted) ; rownames(bc_t2dm_permuted) = NULL
bc_t2dm_permuted = bc_t2dm_permuted %>% dplyr::select(drug_target, everything())
fwrite(bc_t2dm_permuted, "drugs_to_shared_canonical_pathways/BC_T2DM_DIAMANTE_drug_targets_perm_pvalues_unadjusted.txt", sep = "\t", row.names = FALSE)
rm(bc_t2dm_observed, bc_t2dm_permuted, bc_t2dm_sig_cps)

rm(list = ls())
gc()

##### CALCULATE PVALUE FOR EACH DRUG-SHARED CANONICAL PATHWAY PAIR #####
## For each drug, we adjust for the number of gene-targets (bonferonni)

## STRING network
string_high_conf = fread("preprocessed_data/string_high_conf_processed.txt")
string_genes = unique(c(string_high_conf$from, string_high_conf$to))

## BC-DEPRESSION
cp_depression_drug_targets = fread("drugs_to_shared_canonical_pathways/BC_DEPRESSION_drug_targets_perm_pvalues_unadjusted.txt")
cp_depression_drug_targets = tidyr::gather(cp_depression_drug_targets, key = "canonical_pathway", value = "perm_pvalue", -"drug_target")
depression_drug_targets = fread("preprocessed_data/depression_drug_targets.txt") %>%
  filter(entrezID %in% string_genes) %>%
  dplyr::select(drug, entrezID) %>% 
  distinct()
depression_drug_targets = left_join(depression_drug_targets, cp_depression_drug_targets, by = c("entrezID" = "drug_target"))
depression_drug_targets = depression_drug_targets %>%
  dplyr::select(-entrezID) %>%
  group_by(drug, canonical_pathway) %>%
  mutate(perm_pvalue_adj = p.adjust(perm_pvalue, method = "bonferroni")) %>%
  mutate(perm_pvalue_adj = min(perm_pvalue_adj)) %>%
  ungroup() %>%
  dplyr::select(-perm_pvalue) %>%
  distinct()
rm(cp_depression_drug_targets)

## BC-HDL
cp_hdl_drug_targets = fread("drugs_to_shared_canonical_pathways/BC_HDL_drug_targets_perm_pvalues_unadjusted.txt")
cp_hdl_drug_targets = tidyr::gather(cp_hdl_drug_targets, key = "canonical_pathway", value = "perm_pvalue", -"drug_target")
hdl_drug_targets = fread("preprocessed_data/hdl_drug_targets.txt") %>%
  filter(entrezID %in% string_genes) %>%
  dplyr::select(drug, entrezID) %>% 
  distinct()
hdl_drug_targets = left_join(hdl_drug_targets, cp_hdl_drug_targets, by = c("entrezID" = "drug_target"))
hdl_drug_targets = hdl_drug_targets %>%
  dplyr::select(-entrezID) %>%
  group_by(drug, canonical_pathway) %>%
  mutate(perm_pvalue_adj = p.adjust(perm_pvalue, method = "bonferroni")) %>%
  mutate(perm_pvalue_adj = min(perm_pvalue_adj)) %>%
  ungroup() %>%
  dplyr::select(-perm_pvalue) %>%
  distinct()
rm(cp_hdl_drug_targets)

## BC-LDL
cp_ldl_drug_targets = fread("drugs_to_shared_canonical_pathways/BC_LDL_drug_targets_perm_pvalues_unadjusted.txt")
cp_ldl_drug_targets = tidyr::gather(cp_ldl_drug_targets, key = "canonical_pathway", value = "perm_pvalue", -"drug_target")
ldl_drug_targets = fread("preprocessed_data/ldl_drug_targets.txt") %>%
  filter(entrezID %in% string_genes) %>%
  dplyr::select(drug, entrezID) %>% 
  distinct()
ldl_drug_targets = left_join(ldl_drug_targets, cp_ldl_drug_targets, by = c("entrezID" = "drug_target"))
ldl_drug_targets = ldl_drug_targets %>%
  dplyr::select(-entrezID) %>%
  group_by(drug, canonical_pathway) %>%
  mutate(perm_pvalue_adj = p.adjust(perm_pvalue, method = "bonferroni")) %>%
  mutate(perm_pvalue_adj = min(perm_pvalue_adj)) %>%
  ungroup() %>%
  dplyr::select(-perm_pvalue) %>%
  distinct()
rm(cp_ldl_drug_targets)

## BC-PROSTATE
cp_prostate_drug_targets = fread("drugs_to_shared_canonical_pathways/BC_PROSTATE_drug_targets_perm_pvalues_unadjusted.txt")
cp_prostate_drug_targets = tidyr::gather(cp_prostate_drug_targets, key = "canonical_pathway", value = "perm_pvalue", -"drug_target")
prostate_drug_targets = fread("preprocessed_data/prostate_drug_targets.txt") %>%
  filter(entrezID %in% string_genes) %>%
  dplyr::select(drug, entrezID) %>% 
  distinct()
prostate_drug_targets = left_join(prostate_drug_targets, cp_prostate_drug_targets, by = c("entrezID" = "drug_target"))
prostate_drug_targets = prostate_drug_targets %>%
  dplyr::select(-entrezID) %>%
  group_by(drug, canonical_pathway) %>%
  mutate(perm_pvalue_adj = p.adjust(perm_pvalue, method = "bonferroni")) %>%
  mutate(perm_pvalue_adj = min(perm_pvalue_adj)) %>%
  ungroup() %>%
  dplyr::select(-perm_pvalue) %>%
  distinct()
rm(cp_prostate_drug_targets)

## BC-SCZ
cp_scz_drug_targets = fread("drugs_to_shared_canonical_pathways/BC_SCZ_drug_targets_perm_pvalues_unadjusted.txt")
cp_scz_drug_targets = tidyr::gather(cp_scz_drug_targets, key = "canonical_pathway", value = "perm_pvalue", -"drug_target")
scz_drug_targets = fread("preprocessed_data/schizophrenia_drug_targets.txt") %>%
  filter(entrezID %in% string_genes) %>%
  dplyr::select(drug, entrezID) %>% 
  distinct()
scz_drug_targets = left_join(scz_drug_targets, cp_scz_drug_targets, by = c("entrezID" = "drug_target"))
scz_drug_targets = scz_drug_targets %>%
  dplyr::select(-entrezID) %>%
  group_by(drug, canonical_pathway) %>%
  mutate(perm_pvalue_adj = p.adjust(perm_pvalue, method = "bonferroni")) %>%
  mutate(perm_pvalue_adj = min(perm_pvalue_adj)) %>%
  ungroup() %>%
  dplyr::select(-perm_pvalue) %>%
  distinct()
rm(cp_scz_drug_targets)

## BC-T2DM_DIAMANTE
cp_t2dm_drug_targets = fread("drugs_to_shared_canonical_pathways/BC_T2DM_DIAMANTE_drug_targets_perm_pvalues_unadjusted.txt")
cp_t2dm_drug_targets = tidyr::gather(cp_t2dm_drug_targets, key = "canonical_pathway", value = "perm_pvalue", -"drug_target")
t2dm_drug_targets = fread("preprocessed_data/t2dm_drug_targets.txt") %>%
  filter(entrezID %in% string_genes) %>%
  dplyr::select(drug, entrezID) %>% 
  distinct()
t2dm_drug_targets = left_join(t2dm_drug_targets, cp_t2dm_drug_targets, by = c("entrezID" = "drug_target"))
t2dm_drug_targets = t2dm_drug_targets %>%
  dplyr::select(-entrezID) %>%
  group_by(drug, canonical_pathway) %>%
  mutate(perm_pvalue_adj = p.adjust(perm_pvalue, method = "bonferroni")) %>%
  mutate(perm_pvalue_adj = min(perm_pvalue_adj)) %>%
  ungroup() %>%
  dplyr::select(-perm_pvalue) %>%
  distinct()
rm(cp_t2dm_drug_targets)

# save
fwrite(depression_drug_targets, "drugs_to_shared_canonical_pathways/BC_DEPRESSION_drug_targets_perm_pvalues_bonferroni.txt", sep = "\t", row.names = FALSE)
fwrite(hdl_drug_targets, "drugs_to_shared_canonical_pathways/BC_HDL_drug_targets_perm_pvalues_bonferroni.txt", sep = "\t", row.names = FALSE)
fwrite(ldl_drug_targets, "drugs_to_shared_canonical_pathways/BC_LDL_drug_targets_perm_pvalues_bonferroni.txt", sep = "\t", row.names = FALSE)
fwrite(prostate_drug_targets, "drugs_to_shared_canonical_pathways/BC_PROSTATE_drug_targets_perm_pvalues_bonferroni.txt", sep = "\t", row.names = FALSE)
fwrite(scz_drug_targets, "drugs_to_shared_canonical_pathways/BC_SCZ_drug_targets_perm_pvalues_bonferroni.txt", sep = "\t", row.names = FALSE)
fwrite(t2dm_drug_targets, "drugs_to_shared_canonical_pathways/BC_T2DM_drug_targets_perm_pvalues_bonferroni.txt", sep = "\t", row.names = FALSE)
