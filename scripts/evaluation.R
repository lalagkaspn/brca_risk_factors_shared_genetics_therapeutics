
### Evaluating if the number of prioritized candidate drugs for repurposing are enriched for drugs currently investigated or indicated for breast cancer ###

library(dplyr)
library(data.table)

## drug pvalues
depression_drug_targets = fread("drugs_to_shared_canonical_pathways/BC_DEPRESSION_drug_targets_perm_pvalues_bonferroni.txt")
hdl_drug_targets = fread("drugs_to_shared_canonical_pathways//BC_HDL_drug_targets_perm_pvalues_bonferroni.txt")
ldl_drug_targets = fread("drugs_to_shared_canonical_pathways//BC_LDL_drug_targets_perm_pvalues_bonferroni.txt")
prostate_drug_targets = fread("drugs_to_shared_canonical_pathways//BC_PROSTATE_drug_targets_perm_pvalues_bonferroni.txt")
scz_drug_targets = fread("drugs_to_shared_canonical_pathways//BC_SCZ_drug_targets_perm_pvalues_bonferroni.txt")
t2dm_drug_targets = fread("drugs_to_shared_canonical_pathways//BC_T2DM_drug_targets_perm_pvalues_bonferroni.txt")

## Find candidate drugs for repurposing by looking if a drug is significantly connected to at least one shared canonical pathway
bc_depression_recommended_drugs = data.frame(drug = unique(depression_drug_targets$drug), recommended = 2)
for (i in 1:nrow(bc_depression_recommended_drugs)) {
  drug_sig_cps = depression_drug_targets %>% filter(drug == bc_depression_recommended_drugs[i, "drug"]) %>% filter(perm_pvalue_adj < 0.05)
  if (nrow(drug_sig_cps) == 0) {
    bc_depression_recommended_drugs[i, "recommended"] = 0
  } else {
    bc_depression_recommended_drugs[i, "recommended"] = 1
  }
}
bc_hdl_recommended_drugs = data.frame(drug = unique(hdl_drug_targets$drug), recommended = 2)
for (i in 1:nrow(bc_hdl_recommended_drugs)) {
  drug_sig_cps = hdl_drug_targets %>% filter(drug == bc_hdl_recommended_drugs[i, "drug"]) %>% filter(perm_pvalue_adj < 0.05)
  if (nrow(drug_sig_cps) == 0) {
    bc_hdl_recommended_drugs[i, "recommended"] = 0
  } else {
    bc_hdl_recommended_drugs[i, "recommended"] = 1
  }
}
bc_ldl_recommended_drugs = data.frame(drug = unique(ldl_drug_targets$drug), recommended = 2)
for (i in 1:nrow(bc_ldl_recommended_drugs)) {
  drug_sig_cps = ldl_drug_targets %>% filter(drug == bc_ldl_recommended_drugs[i, "drug"]) %>% filter(perm_pvalue_adj < 0.05)
  if (nrow(drug_sig_cps) == 0) {
    bc_ldl_recommended_drugs[i, "recommended"] = 0
  } else {
    bc_ldl_recommended_drugs[i, "recommended"] = 1
  }
}
bc_prostate_recommended_drugs = data.frame(drug = unique(prostate_drug_targets$drug), recommended = 2)
for (i in 1:nrow(bc_prostate_recommended_drugs)) {
  drug_sig_cps = prostate_drug_targets %>% filter(drug == bc_prostate_recommended_drugs[i, "drug"]) %>% filter(perm_pvalue_adj < 0.05)
  if (nrow(drug_sig_cps) == 0) {
    bc_prostate_recommended_drugs[i, "recommended"] = 0
  } else {
    bc_prostate_recommended_drugs[i, "recommended"] = 1
  }
}
bc_scz_recommended_drugs = data.frame(drug = unique(scz_drug_targets$drug), recommended = 2)
for (i in 1:nrow(bc_scz_recommended_drugs)) {
  drug_sig_cps = scz_drug_targets %>% filter(drug == bc_scz_recommended_drugs[i, "drug"]) %>% filter(perm_pvalue_adj < 0.05)
  if (nrow(drug_sig_cps) == 0) {
    bc_scz_recommended_drugs[i, "recommended"] = 0
  } else {
    bc_scz_recommended_drugs[i, "recommended"] = 1
  }
}
bc_t2dm_recommended_drugs = data.frame(drug = unique(t2dm_drug_targets$drug), recommended = 2)
for (i in 1:nrow(bc_t2dm_recommended_drugs)) {
  drug_sig_cps = t2dm_drug_targets %>% filter(drug == bc_t2dm_recommended_drugs[i, "drug"]) %>% filter(perm_pvalue_adj < 0.05)
  if (nrow(drug_sig_cps) == 0) {
    bc_t2dm_recommended_drugs[i, "recommended"] = 0
  } else {
    bc_t2dm_recommended_drugs[i, "recommended"] = 1
  }
}

### EVALUATION ###

### All predisposing diseases ###
# drugs tested
drugs_PAall = rbind(bc_depression_recommended_drugs, bc_hdl_recommended_drugs, bc_ldl_recommended_drugs, bc_prostate_recommended_drugs, bc_scz_recommended_drugs, bc_t2dm_recommended_drugs)

## -- evaluation using drug names -- ##
bc_drugs = fread("preprocessed_data/bc_indicated_investigated_drugs_of_predisposing_diseases.txt")

## Keep minimum p-value for each drug
drugs_PAall = drugs_PAall %>%
  group_by(drug) %>% 
  mutate(recommended = max(recommended)) %>%
  ungroup() %>%
  distinct()

# PAall
contingency_table = matrix(ncol = 2, nrow = 2)
colnames(contingency_table) = c("BC_inv_ind", "BC_NOT_inv_ind")
rownames(contingency_table) = c("Sig", "NOT sig")
temp = drugs_PAall
temp = na.omit(temp) ; rownames(temp) = NULL
temp_sig = temp %>% filter(recommended == 1)
temp_notsig = temp %>% filter(recommended == 0)
contingency_table[1, 1] = length(intersect(temp_sig$drug, tolower(bc_drugs$drug_name)))
contingency_table[1, 2] = length(setdiff(temp_sig$drug, tolower(bc_drugs$drug_name)))
contingency_table[2, 1] = length(intersect(temp_notsig$drug, tolower(bc_drugs$drug_name)))
contingency_table[2, 2] = nrow(temp) - contingency_table[1, 1] - contingency_table[1, 2] - contingency_table[2, 1]
fisher.test(contingency_table, alternative = "greater")
contingency_table

### All diseases but LDL ###
# drugs tested
drugs_PAall = rbind(bc_depression_recommended_drugs, 
                    bc_hdl_recommended_drugs, 
                    # bc_ldl_recommended_drugs,
                    bc_prostate_recommended_drugs,
                    bc_scz_recommended_drugs,
                    bc_t2dm_recommended_drugs)

## Keep minimum p-value for each drug
drugs_PAall = drugs_PAall %>%
  group_by(drug) %>% 
  mutate(recommended = max(recommended)) %>%
  ungroup() %>%
  distinct()

# PAall
contingency_table = matrix(ncol = 2, nrow = 2)
colnames(contingency_table) = c("BC_inv_ind", "BC_NOT_inv_ind")
rownames(contingency_table) = c("Sig", "NOT sig")
temp = drugs_PAall
temp = na.omit(temp) ; rownames(temp) = NULL
temp_sig = temp %>% filter(recommended == 1)
temp_notsig = temp %>% filter(recommended == 0)
contingency_table[1, 1] = length(intersect(temp_sig$drug, tolower(bc_drugs$drug_name)))
contingency_table[1, 2] = length(setdiff(temp_sig$drug, tolower(bc_drugs$drug_name)))
contingency_table[2, 1] = length(intersect(temp_notsig$drug, tolower(bc_drugs$drug_name)))
contingency_table[2, 2] = nrow(temp) - contingency_table[1, 1] - contingency_table[1, 2] - contingency_table[2, 1]
fisher.test(contingency_table, alternative = "greater")
contingency_table
