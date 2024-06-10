
### Add LOGODetect region q-value to shared genes for each BRCA-predisposing disease pair ###

library(dplyr)

### BC - Depression
# LOGODetect positively correlated regions
bc_depression_regions = data.table::fread("logodetect_regions/brca_depression_logodetect_regions.txt") %>%
  filter(stat > 0)
# mapped genes
bc_depression_genes_pos_mapping = data.table::fread("fuma_snp2gene/mapped_genes/brca_depression_fuma_snp2gene_genes.txt") %>%
  dplyr::select(gene_name = symbol, entrezID, chr, start, end) %>%
  na.omit()
# add q-values to genes
bc_depression_genes_pos_mapping$qvalue = 0
for (i in 1:nrow(bc_depression_genes_pos_mapping)) {
  chromosome = bc_depression_genes_pos_mapping[i, chr]
  start = bc_depression_genes_pos_mapping[i, start]
  stop = bc_depression_genes_pos_mapping[i, end]
  value = bc_depression_regions %>% filter(chr == chromosome & ((begin_pos %in% (start - 10000):(stop + 10000) | stop_pos %in% (start - 10000):(stop + 10000)) | (begin_pos <= start & stop_pos >= stop))) %>% dplyr::select(qval)
  value = value$qval
  bc_depression_genes_pos_mapping[i, "qvalue"] = value
}
# filter for needed columns
bc_depression_genes_pos_mapping = bc_depression_genes_pos_mapping %>% dplyr::select(gene_name, entrezID, qvalue) %>% distinct()
# save
data.table::fwrite(bc_depression_genes_pos_mapping, "fuma_snp2gene/mapped_genes/brca_depression_fuma_snp2gene_genes_logodetect_qvalues.txt", sep = "\t", row.names = FALSE)

### BC - HDL
# LOGODetect regions positively correlated regions
bc_hdl_regions = data.table::fread("logodetect_regions/brca_hdl_logodetect_regions.txt") %>%
  filter(stat > 0)
# mapped genes
bc_hdl_genes_pos_mapping = data.table::fread("fuma_snp2gene/mapped_genes/brca_hdl_fuma_snp2gene_genes.txt") %>%
  dplyr::select(gene_name = symbol, entrezID, chr, start, end) %>%
  na.omit()
# add q-values to genes
bc_hdl_genes_pos_mapping$qvalue = 0
for (i in 1:nrow(bc_hdl_genes_pos_mapping)) {
  chromosome = bc_hdl_genes_pos_mapping[i, chr]
  start = bc_hdl_genes_pos_mapping[i, start]
  stop = bc_hdl_genes_pos_mapping[i, end]
  value = bc_hdl_regions %>% filter(chr == chromosome & ((begin_pos %in% (start - 10000):(stop + 10000) | stop_pos %in% (start - 10000):(stop + 10000)) | (begin_pos <= start & stop_pos >= stop))) %>% dplyr::select(qval)
  value = value$qval
  bc_hdl_genes_pos_mapping[i, "qvalue"] = value
}
# filter for needed columns
bc_hdl_genes_pos_mapping = bc_hdl_genes_pos_mapping %>% dplyr::select(gene_name, entrezID, qvalue) %>% distinct()
# save
data.table::fwrite(bc_hdl_genes_pos_mapping, "fuma_snp2gene/mapped_genes/brca_hdl_fuma_snp2gene_genes_logodetect_qvalues.txt", sep = "\t", row.names = FALSE)

### BC - LDL
# LOGODetect positively correlated regions
bc_ldl_regions = data.table::fread("logodetect_regions/brca_ldl_logodetect_regions.txt") %>%
  filter(stat > 0)
# mapped genes
bc_ldl_genes_pos_mapping = data.table::fread("fuma_snp2gene/mapped_genes/brca_ldl_fuma_snp2gene_genes.txt") %>%
  dplyr::select(gene_name = symbol, entrezID, chr, start, end) %>%
  na.omit()
# add q-values to genes
bc_ldl_genes_pos_mapping$qvalue = 0
for (i in 1:nrow(bc_ldl_genes_pos_mapping)) {
  chromosome = bc_ldl_genes_pos_mapping[i, chr]
  start = bc_ldl_genes_pos_mapping[i, start]
  stop = bc_ldl_genes_pos_mapping[i, end]
  value = bc_ldl_regions %>% filter(chr == chromosome & ((begin_pos %in% (start - 10000):(stop + 10000) | stop_pos %in% (start - 10000):(stop + 10000)) | (begin_pos <= start & stop_pos >= stop))) %>% dplyr::select(qval)
  value = value$qval
  bc_ldl_genes_pos_mapping[i, "qvalue"] = value
}
# filter for needed columns
bc_ldl_genes_pos_mapping = bc_ldl_genes_pos_mapping %>% dplyr::select(gene_name, entrezID, qvalue) %>% distinct()
# save
data.table::fwrite(bc_ldl_genes_pos_mapping, "fuma_snp2gene/mapped_genes/brca_ldl_fuma_snp2gene_genes_logodetect_qvalues.txt", sep = "\t", row.names = FALSE)

### BC - Prostate
# LOGODetect positively correlated regions
bc_prostate_regions = data.table::fread("logodetect_regions/brca_prostate_logodetect_regions.txt") %>%
  filter(stat > 0)
# mapped genes
bc_prostate_genes_pos_mapping = data.table::fread("fuma_snp2gene/mapped_genes/brca_prostate_fuma_snp2gene_genes.txt") %>%
  dplyr::select(gene_name = symbol, entrezID, chr, start, end) %>%
  na.omit()
# add q-values to genes
bc_prostate_genes_pos_mapping$qvalue = 0
for (i in 1:nrow(bc_prostate_genes_pos_mapping)) {
  chromosome = bc_prostate_genes_pos_mapping[i, chr]
  start = bc_prostate_genes_pos_mapping[i, start]
  stop = bc_prostate_genes_pos_mapping[i, end]
  value = bc_prostate_regions %>% filter(chr == chromosome & ((begin_pos %in% (start - 10000):(stop + 10000) | stop_pos %in% (start - 10000):(stop + 10000)) | (begin_pos <= start & stop_pos >= stop))) %>% dplyr::select(qval)
  value = value$qval
  bc_prostate_genes_pos_mapping[i, "qvalue"] = value
}
# filter for needed columns
bc_prostate_genes_pos_mapping = bc_prostate_genes_pos_mapping %>% dplyr::select(gene_name, entrezID, qvalue) %>% distinct()
# save
data.table::fwrite(bc_prostate_genes_pos_mapping, "fuma_snp2gene/mapped_genes/brca_prostate_fuma_snp2gene_genes_logodetect_qvalues.txt", sep = "\t", row.names = FALSE)

### BC - Schizophrenia
# LOGODetect positively correlated regions
bc_scz_regions = data.table::fread("logodetect_regions/brca_schizophrenia_regions.txt") %>%
  filter(stat > 0)
# mapped genes
bc_scz_genes_pos_mapping = data.table::fread("fuma_snp2gene/mapped_genes/brca_schizophrenia_fuma_snp2gene_genes.txt") %>%
  dplyr::select(gene_name = symbol, entrezID, chr, start, end) %>%
  na.omit()
# add q-values to genes
bc_scz_genes_pos_mapping$qvalue = 0
for (i in 1:nrow(bc_scz_genes_pos_mapping)) {
  chromosome = bc_scz_genes_pos_mapping[i, chr]
  start = bc_scz_genes_pos_mapping[i, start]
  stop = bc_scz_genes_pos_mapping[i, end]
  value = bc_scz_regions %>% filter(chr == chromosome & ((begin_pos %in% (start - 10000):(stop + 10000) | stop_pos %in% (start - 10000):(stop + 10000)) | (begin_pos <= start & stop_pos >= stop))) %>% dplyr::select(qval)
  value = value$qval
  bc_scz_genes_pos_mapping[i, "qvalue"] = value
}
# filter for needed columns
bc_scz_genes_pos_mapping = bc_scz_genes_pos_mapping %>% dplyr::select(gene_name, entrezID, qvalue) %>% distinct()
# save
data.table::fwrite(bc_scz_genes_pos_mapping, "fuma_snp2gene/mapped_genes/brca_schizophrenia_fuma_snp2gene_genes_logodetect_qvalues.txt", sep = "\t", row.names = FALSE)

### BC - T2DM_DIAMANTE
# LOGODetect positively correlated regions
bc_t2dm_diamante_regions = data.table::fread("logodetect_regions/brca_t2dm_logodetect_regions.txt") %>%
  filter(stat > 0)
# mapped genes
bc_t2dm_diamante_genes_pos_mapping = data.table::fread("fuma_snp2gene/mapped_genes/brca_t2dm_fuma_snp2gene_genes.txt") %>%
  dplyr::select(gene_name = symbol, entrezID, chr, start, end) %>%
  na.omit()
# add q-values to genes
bc_t2dm_diamante_genes_pos_mapping$qvalue = 0
for (i in 1:nrow(bc_t2dm_diamante_genes_pos_mapping)) {
  chromosome = bc_t2dm_diamante_genes_pos_mapping[i, chr]
  start = bc_t2dm_diamante_genes_pos_mapping[i, start]
  stop = bc_t2dm_diamante_genes_pos_mapping[i, end]
  value = bc_t2dm_diamante_regions %>% filter(chr == chromosome & ((begin_pos %in% (start - 10000):(stop + 10000) | stop_pos %in% (start - 10000):(stop + 10000)) | (begin_pos <= start & stop_pos >= stop))) %>% dplyr::select(qval)
  value = value$qval
  bc_t2dm_diamante_genes_pos_mapping[i, "qvalue"] = value
}
# filter for needed columns
bc_t2dm_diamante_genes_pos_mapping = bc_t2dm_diamante_genes_pos_mapping %>% dplyr::select(gene_name, entrezID, qvalue) %>% distinct()
# save
data.table::fwrite(bc_t2dm_diamante_genes_pos_mapping, "fuma_snp2gene/mapped_genes/brca_t2dm_fuma_snp2gene_genes_logodetect_qvalues.txt", sep = "\t", row.names = FALSE)
