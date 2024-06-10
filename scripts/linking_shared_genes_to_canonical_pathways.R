
### Find canonical pathways closely connected to filtered shared genes for each disease pair ###
### NOTE: this script requires >64GB RAM to run

library(GSA)
library(data.table)
library(dplyr)
library(parallel)
library(igraph)
library(doParallel)
library(foreach)
# HPC test
print("Packages successfully loaded!")

## load MSigDB canonical pathways with genes in STRING ##
# msigdb = GSA.read.gmt("https://zenodo.org/records/11540678/files/c2.cp.v2023.2.Hs.entrez.gmt?download=1")
# # name the lists
# names(msigdb$genesets) = msigdb$geneset.names
# msigdb = msigdb$genesets # 3,795 pathways
# ## filter for genes found in the STRING network
# string_high_conf = fread("preprocessed_data/string_high_conf_processed.txt")
# genes_in_string = unique(c(string_high_conf$from, string_high_conf$to))
# # Set up parallel processing
# num_cores = 4
# cl = makeCluster(num_cores)
# clusterExport(cl, c("msigdb", "genes_in_string"))
# # Parallel intersection using parLapply
# msigdb <- parLapply(cl, msigdb, function(x) {
#   intersect(x, genes_in_string)
# })
# # Stop the cluster
# stopCluster(cl)
# saveRDS(msigdb, "preprocessed_data/c2.cp.v2023.2.Hs.entrez_STRING_filtered.RDS")
msigdb = readRDS("preprocessed_data/c2.cp.v2023.2.Hs.entrez_STRING_filtered.RDS")


##### OBSERVED PERSONALIZED PAGE RANK SCORES #####

### BC-Depression ###
## -- load shared genes -- ##
## MAGMA significant for depression and S-MultiXcan same direction for BC and depression
bc_depression_shared_genes = fread("magma_smultixcan_filtered_shared_genes/bc_depression_filtered_shared_genes.txt")
## find average PPR score of the genes in a canonical pathway using the shared genes as seed genes ##
string_high_conf = fread("preprocessed_data/string_high_conf_processed.txt")
g = graph_from_data_frame(string_high_conf, directed = FALSE)

## observed
avg_ppr_obs = vector(mode = "list", length = length(msigdb))
names(avg_ppr_obs) = names(msigdb)

# Set up parallel processing
cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) # detectCores() doesn't work with the way the Slurm limits the number of cores available.
print(paste0("Number of cores used for parallelization: ", cores))
registerDoParallel(cores)

# Parallelized loop
avg_ppr_obs = foreach(cp = 1:length(avg_ppr_obs), .packages = c("igraph", "dplyr")) %dopar% {
  genes_in_cp = msigdb[[cp]]
  obs_avg_ppr_score_per_shared_gene = data.frame(shared_gene = unique(bc_depression_shared_genes$entrezID),
                                                 average_ppr_score_cp = 0,
                                                 cp = names(avg_ppr_obs)[cp])
  for (i in 1:nrow(obs_avg_ppr_score_per_shared_gene)) {
    shared_gene_temp = as.character(obs_avg_ppr_score_per_shared_gene[i, "shared_gene"])
    p = c(rep(0, length(V(g))))
    p[which(V(g)$name == shared_gene_temp)] = 1
    pr_result = page_rank(g, algo = "prpack", damping = 0.85, personalized = p, weights = NA)
    pr_result = pr_result$vector
    pr_result = as.data.frame(pr_result)
    pr_result$entrezID = rownames(pr_result)
    rownames(pr_result) = NULL
    pr_result = pr_result %>%
      dplyr::select(entrezID, propagation_score = pr_result) %>%
      filter(entrezID %in% genes_in_cp) %>%
      mutate(avg_propagation_score = mean(propagation_score)) %>%
      dplyr::select(avg_propagation_score) %>%
      distinct()
    obs_avg_ppr_score_per_shared_gene[i, "average_ppr_score_cp"] = pr_result
  }
  return(obs_avg_ppr_score_per_shared_gene)
}
names(avg_ppr_obs) = names(msigdb)
saveRDS(avg_ppr_obs, "shared_genes_to_canonical_pathways/observed_BC_DEPRESSION.RDS")

### BC-HDL ###
## -- load shared genes -- ##
bc_hdl_shared_genes = fread("magma_smultixcan_filtered_shared_genes/bc_hdl_filtered_shared_genes.txt")
## observed
avg_ppr_obs = vector(mode = "list", length = length(msigdb))
names(avg_ppr_obs) = names(msigdb)
# Set up parallel processing
cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) # detectCores() doesn't work with the way the Slurm limits the number of cores available.
print(paste0("Number of cores used for parallelization: ", cores))
registerDoParallel(cores)
# Parallelized loop
avg_ppr_obs = foreach(cp = 1:length(avg_ppr_obs), .packages = c("igraph", "dplyr")) %dopar% {
  genes_in_cp = msigdb[[cp]]
  obs_avg_ppr_score_per_shared_gene = data.frame(shared_gene = unique(bc_hdl_shared_genes$entrezID),
                                                 average_ppr_score_cp = 0,
                                                 cp = names(avg_ppr_obs)[cp])
  for (i in 1:nrow(obs_avg_ppr_score_per_shared_gene)) {
    shared_gene_temp = as.character(obs_avg_ppr_score_per_shared_gene[i, "shared_gene"])
    p = c(rep(0, length(V(g))))
    p[which(V(g)$name == shared_gene_temp)] = 1
    pr_result = page_rank(g, algo = "prpack", damping = 0.85, personalized = p, weights = NA)
    pr_result = pr_result$vector
    pr_result = as.data.frame(pr_result)
    pr_result$entrezID = rownames(pr_result)
    rownames(pr_result) = NULL
    pr_result = pr_result %>%
      dplyr::select(entrezID, propagation_score = pr_result) %>%
      filter(entrezID %in% genes_in_cp) %>%
      mutate(avg_propagation_score = mean(propagation_score)) %>%
      dplyr::select(avg_propagation_score) %>%
      distinct()
    obs_avg_ppr_score_per_shared_gene[i, "average_ppr_score_cp"] = pr_result
  }
  return(obs_avg_ppr_score_per_shared_gene)
}
names(avg_ppr_obs) = names(msigdb)
saveRDS(avg_ppr_obs, "shared_genes_to_canonical_pathways/observed_BC_HDL.RDS")

### BC-LDL ###
## -- load shared genes -- ##
bc_ldl_shared_genes = fread("magma_smultixcan_filtered_shared_genes/bc_ldl_filtered_shared_genes.txt")
## observed
avg_ppr_obs = vector(mode = "list", length = length(msigdb))
names(avg_ppr_obs) = names(msigdb)
# Set up parallel processing
cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) # detectCores() doesn't work with the way the Slurm limits the number of cores available.
print(paste0("Number of cores used for parallelization: ", cores))
registerDoParallel(cores)
# Parallelized loop
avg_ppr_obs = foreach(cp = 1:length(avg_ppr_obs), .packages = c("igraph", "dplyr")) %dopar% {
  genes_in_cp = msigdb[[cp]]
  obs_avg_ppr_score_per_shared_gene = data.frame(shared_gene = unique(bc_ldl_shared_genes$entrezID),
                                                 average_ppr_score_cp = 0,
                                                 cp = names(avg_ppr_obs)[cp])
  for (i in 1:nrow(obs_avg_ppr_score_per_shared_gene)) {
    shared_gene_temp = as.character(obs_avg_ppr_score_per_shared_gene[i, "shared_gene"])
    p = c(rep(0, length(V(g))))
    p[which(V(g)$name == shared_gene_temp)] = 1
    pr_result = page_rank(g, algo = "prpack", damping = 0.85, personalized = p, weights = NA)
    pr_result = pr_result$vector
    pr_result = as.data.frame(pr_result)
    pr_result$entrezID = rownames(pr_result)
    rownames(pr_result) = NULL
    pr_result = pr_result %>%
      dplyr::select(entrezID, propagation_score = pr_result) %>%
      filter(entrezID %in% genes_in_cp) %>%
      mutate(avg_propagation_score = mean(propagation_score)) %>%
      dplyr::select(avg_propagation_score) %>%
      distinct()
    obs_avg_ppr_score_per_shared_gene[i, "average_ppr_score_cp"] = pr_result
  }
  return(obs_avg_ppr_score_per_shared_gene)
}
names(avg_ppr_obs) = names(msigdb)
saveRDS(avg_ppr_obs, "shared_genes_to_canonical_pathways/observed_BC_LDL.RDS")

### BC-Prostate cancer ###
## -- load shared genes -- ##
bc_prostate_shared_genes = fread("magma_smultixcan_filtered_shared_genes/bc_prostate_filtered_shared_genes.txt")
## observed
avg_ppr_obs = vector(mode = "list", length = length(msigdb))
names(avg_ppr_obs) = names(msigdb)
# Set up parallel processing
cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) # detectCores() doesn't work with the way the Slurm limits the number of cores available.
print(paste0("Number of cores used for parallelization: ", cores))
registerDoParallel(cores)
# Parallelized loop
avg_ppr_obs = foreach(cp = 1:length(avg_ppr_obs), .packages = c("igraph", "dplyr")) %dopar% {
  genes_in_cp = msigdb[[cp]]
  obs_avg_ppr_score_per_shared_gene = data.frame(shared_gene = unique(bc_prostate_shared_genes$entrezID),
                                                 average_ppr_score_cp = 0,
                                                 cp = names(avg_ppr_obs)[cp])
  for (i in 1:nrow(obs_avg_ppr_score_per_shared_gene)) {
    shared_gene_temp = as.character(obs_avg_ppr_score_per_shared_gene[i, "shared_gene"])
    p = c(rep(0, length(V(g))))
    p[which(V(g)$name == shared_gene_temp)] = 1
    pr_result = page_rank(g, algo = "prpack", damping = 0.85, personalized = p, weights = NA)
    pr_result = pr_result$vector
    pr_result = as.data.frame(pr_result)
    pr_result$entrezID = rownames(pr_result)
    rownames(pr_result) = NULL
    pr_result = pr_result %>%
      dplyr::select(entrezID, propagation_score = pr_result) %>%
      filter(entrezID %in% genes_in_cp) %>%
      mutate(avg_propagation_score = mean(propagation_score)) %>%
      dplyr::select(avg_propagation_score) %>%
      distinct()
    obs_avg_ppr_score_per_shared_gene[i, "average_ppr_score_cp"] = pr_result
  }
  return(obs_avg_ppr_score_per_shared_gene)
}
names(avg_ppr_obs) = names(msigdb)
saveRDS(avg_ppr_obs, "shared_genes_to_canonical_pathways/observed_BC_PROSTATE.RDS")

### BC-Schizophrenia ###
## -- load shared genes -- ##
bc_scz_shared_genes = fread("magma_smultixcan_filtered_shared_genes/bc_schizophrenia_filtered_shared_genes.txt")
## observed
avg_ppr_obs = vector(mode = "list", length = length(msigdb))
names(avg_ppr_obs) = names(msigdb)
# Set up parallel processing
cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) # detectCores() doesn't work with the way the Slurm limits the number of cores available.
print(paste0("Number of cores used for parallelization: ", cores))
registerDoParallel(cores)
# Parallelized loop
avg_ppr_obs = foreach(cp = 1:length(avg_ppr_obs), .packages = c("igraph", "dplyr")) %dopar% {
  genes_in_cp = msigdb[[cp]]
  obs_avg_ppr_score_per_shared_gene = data.frame(shared_gene = unique(bc_scz_shared_genes$entrezID),
                                                 average_ppr_score_cp = 0,
                                                 cp = names(avg_ppr_obs)[cp])
  for (i in 1:nrow(obs_avg_ppr_score_per_shared_gene)) {
    shared_gene_temp = as.character(obs_avg_ppr_score_per_shared_gene[i, "shared_gene"])
    p = c(rep(0, length(V(g))))
    p[which(V(g)$name == shared_gene_temp)] = 1
    pr_result = page_rank(g, algo = "prpack", damping = 0.85, personalized = p, weights = NA)
    pr_result = pr_result$vector
    pr_result = as.data.frame(pr_result)
    pr_result$entrezID = rownames(pr_result)
    rownames(pr_result) = NULL
    pr_result = pr_result %>%
      dplyr::select(entrezID, propagation_score = pr_result) %>%
      filter(entrezID %in% genes_in_cp) %>%
      mutate(avg_propagation_score = mean(propagation_score)) %>%
      dplyr::select(avg_propagation_score) %>%
      distinct()
    obs_avg_ppr_score_per_shared_gene[i, "average_ppr_score_cp"] = pr_result
  }
  return(obs_avg_ppr_score_per_shared_gene)
}
names(avg_ppr_obs) = names(msigdb)
saveRDS(avg_ppr_obs, "shared_genes_to_canonical_pathways/observed_BC_SCZ.RDS")

### BC-T2DM ###
## -- load shared genes -- ##
bc_t2dm_shared_genes = fread("magma_smultixcan_filtered_shared_genes/bc_t2dm_diamante_filtered_shared_genes.txt")
## observed
avg_ppr_obs = vector(mode = "list", length = length(msigdb))
names(avg_ppr_obs) = names(msigdb)
# Set up parallel processing
cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) # detectCores() doesn't work with the way the Slurm limits the number of cores available.
print(paste0("Number of cores used for parallelization: ", cores))
registerDoParallel(cores)
# Parallelized loop
avg_ppr_obs = foreach(cp = 1:length(avg_ppr_obs), .packages = c("igraph", "dplyr")) %dopar% {
  genes_in_cp = msigdb[[cp]]
  obs_avg_ppr_score_per_shared_gene = data.frame(shared_gene = unique(bc_t2dm_shared_genes$entrezID),
                                                 average_ppr_score_cp = 0,
                                                 cp = names(avg_ppr_obs)[cp])
  for (i in 1:nrow(obs_avg_ppr_score_per_shared_gene)) {
    shared_gene_temp = as.character(obs_avg_ppr_score_per_shared_gene[i, "shared_gene"])
    p = c(rep(0, length(V(g))))
    p[which(V(g)$name == shared_gene_temp)] = 1
    pr_result = page_rank(g, algo = "prpack", damping = 0.85, personalized = p, weights = NA)
    pr_result = pr_result$vector
    pr_result = as.data.frame(pr_result)
    pr_result$entrezID = rownames(pr_result)
    rownames(pr_result) = NULL
    pr_result = pr_result %>%
      dplyr::select(entrezID, propagation_score = pr_result) %>%
      filter(entrezID %in% genes_in_cp) %>%
      mutate(avg_propagation_score = mean(propagation_score)) %>%
      dplyr::select(avg_propagation_score) %>%
      distinct()
    obs_avg_ppr_score_per_shared_gene[i, "average_ppr_score_cp"] = pr_result
  }
  return(obs_avg_ppr_score_per_shared_gene)
}
names(avg_ppr_obs) = names(msigdb)
saveRDS(avg_ppr_obs, "shared_genes_to_canonical_pathways/observed_BC_T2DM_DIAMANTE.RDS")



##### PERMUTATIONS - REPLACING THE SHARED GENES WITH RANDOM GENES, DEGREE MATCHED, TO THE PERSONALIZED PAGE RANK #####

node_deg = degree(g)
bins = list(first = seq(1,4,1), second = seq(5,13,1), third = seq(14,36, 1), fourth = seq(37,729, 1)) # quantiles

## BC-Depression ##
avg_ppr_perm = vector(mode = "list", length = 1000)
names(avg_ppr_perm) = paste0("permutation_", 1:1000)

# Set up parallel processing
cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) # detectCores() doesn't work with the way the Slurm limits the number of cores available.
print(paste0("Number of cores used for parallelization: ", cores))
registerDoParallel(cores)
# Parallelized loop
for (perm in 1:length(avg_ppr_perm)) {
  
  ## -------------------------------------- shuffle shared genes based on degree bins --------------------------------------------- ##
  perm_bc_depression_shared_genes = data.frame(original_shared_entrezID = NA, perm_shared_entrezID = NA)
  for (i in 1:length(bins)) {
    # find entrezIDs in this degree bin
    nodes_within_bin = as.data.frame(names(node_deg[which(node_deg %in% bins[[i]])]))
    colnames(nodes_within_bin) = "entrezID"
    # annotate the ones that are shared genes
    nodes_within_bin$original_entrezID = ifelse(nodes_within_bin$entrezID %in% bc_depression_shared_genes$entrezID, 1, 0)
    # shuffle the shared genes
    random_seed_nodes = sample(nodes_within_bin$entrezID, size = sum(nodes_within_bin$original_entrezID))
    # annotate the permuted shared genes
    nodes_within_bin$perm_shared_genes = ifelse(nodes_within_bin$entrezID %in% random_seed_nodes, 1, 0)
    new_shared_genes = nodes_within_bin %>% filter(original_entrezID == 1) %>% dplyr::select(original_shared_entrezID = entrezID)
    new_shared_genes = cbind(new_shared_genes,
                             nodes_within_bin %>% filter(perm_shared_genes == 1) %>% dplyr::select(perm_shared_entrezID = entrezID))
    perm_bc_depression_shared_genes = rbind(perm_bc_depression_shared_genes, new_shared_genes)
  }
  perm_bc_depression_shared_genes = perm_bc_depression_shared_genes[-1, ] ; rownames(perm_bc_depression_shared_genes) = NULL
  
  ## --------------------- parallelize the PPR runs for each canonical pathway per gene for each permutation ---------------------- ##
  cp_perm_genes = vector(mode = "list", length = length(msigdb))
  names(cp_perm_genes) = names(msigdb)
  
  cp_perm_genes = foreach(cp = 1:length(msigdb), .packages = c("igraph", "dplyr")) %dopar% {
    genes_in_cp = msigdb[[cp]]
    perm_avg_ppr_score_per_shared_gene = data.frame(shared_gene = unique(perm_bc_depression_shared_genes$perm_shared_entrezID),
                                                    average_ppr_score_cp = 0,
                                                    cp = names(msigdb)[cp])
    for (i in 1:nrow(perm_avg_ppr_score_per_shared_gene)) {
      shared_gene_temp = as.character(perm_avg_ppr_score_per_shared_gene[i, "shared_gene"])
      p = c(rep(0, length(V(g))))
      p[which(V(g)$name == shared_gene_temp)] = 1
      pr_result = page_rank(g, algo = "prpack", damping = 0.85, personalized = p, weights = NA)
      pr_result = pr_result$vector
      pr_result = as.data.frame(pr_result)
      pr_result$entrezID = rownames(pr_result)
      rownames(pr_result) = NULL
      pr_result = pr_result %>%
        dplyr::select(entrezID, propagation_score = pr_result) %>%
        filter(entrezID %in% genes_in_cp) %>%
        mutate(avg_propagation_score = mean(propagation_score)) %>%
        dplyr::select(avg_propagation_score) %>%
        distinct()
      perm_avg_ppr_score_per_shared_gene[i, "average_ppr_score_cp"] = pr_result
    }
    perm_avg_ppr_score_per_shared_gene = left_join(perm_avg_ppr_score_per_shared_gene, perm_bc_depression_shared_genes, by = c("shared_gene" = "perm_shared_entrezID"))
    perm_avg_ppr_score_per_shared_gene = perm_avg_ppr_score_per_shared_gene %>%
      dplyr::select(original_shared_entrezID, perm_shared_gene = shared_gene, average_ppr_score_cp_perm_gene = average_ppr_score_cp, cp)
    return(perm_avg_ppr_score_per_shared_gene)
  }
  names(cp_perm_genes) = names(msigdb)
  # populate list with the results from one permutation
  avg_ppr_perm[[perm]] = cp_perm_genes
  # track progress
  cat(perm, "- 1000", "\n")
}
saveRDS(avg_ppr_perm, "shared_genes_to_canonical_pathways/permuted_BC_DEPRESSION.RDS")

## BC-HDL ##
avg_ppr_perm = vector(mode = "list", length = 1000)
names(avg_ppr_perm) = paste0("permutation_", 1:1000)
# Set up parallel processing
cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) # detectCores() doesn't work with the way the Slurm limits the number of cores available.
print(paste0("Number of cores used for parallelization: ", cores))
registerDoParallel(cores)
# Parallelized loop
for (perm in 1:length(avg_ppr_perm)) {
  ## -------------------------------------- shuffle shared genes based on degree bins --------------------------------------------- ##
  perm_bc_hdl_shared_genes = data.frame(original_shared_entrezID = NA, perm_shared_entrezID = NA)
  for (i in 1:length(bins)) {
    # find entrezIDs in this degree bin
    nodes_within_bin = as.data.frame(names(node_deg[which(node_deg %in% bins[[i]])]))
    colnames(nodes_within_bin) = "entrezID"
    # annotate the ones that are shared genes
    nodes_within_bin$original_entrezID = ifelse(nodes_within_bin$entrezID %in% bc_hdl_shared_genes$entrezID, 1, 0)
    # shuffle the shared genes
    random_seed_nodes = sample(nodes_within_bin$entrezID, size = sum(nodes_within_bin$original_entrezID))
    # annotate the permuted shared genes
    nodes_within_bin$perm_shared_genes = ifelse(nodes_within_bin$entrezID %in% random_seed_nodes, 1, 0)
    new_shared_genes = nodes_within_bin %>% filter(original_entrezID == 1) %>% dplyr::select(original_shared_entrezID = entrezID)
    new_shared_genes = cbind(new_shared_genes,
                             nodes_within_bin %>% filter(perm_shared_genes == 1) %>% dplyr::select(perm_shared_entrezID = entrezID))
    perm_bc_hdl_shared_genes = rbind(perm_bc_hdl_shared_genes, new_shared_genes)
  }
  perm_bc_hdl_shared_genes = perm_bc_hdl_shared_genes[-1, ] ; rownames(perm_bc_hdl_shared_genes) = NULL
  
  ## --------------------- parallelize the PPR runs for each canonical pathway per gene for each permutation ---------------------- ##
  cp_perm_genes = vector(mode = "list", length = length(msigdb))
  names(cp_perm_genes) = names(msigdb)
  
  cp_perm_genes = foreach(cp = 1:length(msigdb), .packages = c("igraph", "dplyr")) %dopar% {
    genes_in_cp = msigdb[[cp]]
    perm_avg_ppr_score_per_shared_gene = data.frame(shared_gene = unique(perm_bc_hdl_shared_genes$perm_shared_entrezID),
                                                    average_ppr_score_cp = 0,
                                                    cp = names(msigdb)[cp])
    for (i in 1:nrow(perm_avg_ppr_score_per_shared_gene)) {
      shared_gene_temp = as.character(perm_avg_ppr_score_per_shared_gene[i, "shared_gene"])
      p = c(rep(0, length(V(g))))
      p[which(V(g)$name == shared_gene_temp)] = 1
      pr_result = page_rank(g, algo = "prpack", damping = 0.85, personalized = p, weights = NA)
      pr_result = pr_result$vector
      pr_result = as.data.frame(pr_result)
      pr_result$entrezID = rownames(pr_result)
      rownames(pr_result) = NULL
      pr_result = pr_result %>%
        dplyr::select(entrezID, propagation_score = pr_result) %>%
        filter(entrezID %in% genes_in_cp) %>%
        mutate(avg_propagation_score = mean(propagation_score)) %>%
        dplyr::select(avg_propagation_score) %>%
        distinct()
      perm_avg_ppr_score_per_shared_gene[i, "average_ppr_score_cp"] = pr_result
    }
    perm_avg_ppr_score_per_shared_gene = left_join(perm_avg_ppr_score_per_shared_gene, perm_bc_hdl_shared_genes, by = c("shared_gene" = "perm_shared_entrezID"))
    perm_avg_ppr_score_per_shared_gene = perm_avg_ppr_score_per_shared_gene %>%
      dplyr::select(original_shared_entrezID, perm_shared_gene = shared_gene, average_ppr_score_cp_perm_gene = average_ppr_score_cp, cp)
    return(perm_avg_ppr_score_per_shared_gene)
  }
  names(cp_perm_genes) = names(msigdb)
  # populate list with the results from one permutation
  avg_ppr_perm[[perm]] = cp_perm_genes
  # track progress
  cat(perm, "- 1000", "\n")
}
saveRDS(avg_ppr_perm, "shared_genes_to_canonical_pathways/permuted_BC_HDL.RDS")

## BC-LDL ##
avg_ppr_perm = vector(mode = "list", length = 1000)
names(avg_ppr_perm) = paste0("permutation_", 1:1000)
# Set up parallel processing
cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) # detectCores() doesn't work with the way the Slurm limits the number of cores available.
print(paste0("Number of cores used for parallelization: ", cores))
registerDoParallel(cores)
# Parallelized loop
for (perm in 1:length(avg_ppr_perm)) {
  ## -------------------------------------- shuffle shared genes based on degree bins --------------------------------------------- ##
  perm_bc_ldl_shared_genes = data.frame(original_shared_entrezID = NA, perm_shared_entrezID = NA)
  for (i in 1:length(bins)) {
    # find entrezIDs in this degree bin
    nodes_within_bin = as.data.frame(names(node_deg[which(node_deg %in% bins[[i]])]))
    colnames(nodes_within_bin) = "entrezID"
    # annotate the ones that are shared genes
    nodes_within_bin$original_entrezID = ifelse(nodes_within_bin$entrezID %in% bc_ldl_shared_genes$entrezID, 1, 0)
    # shuffle the shared genes
    random_seed_nodes = sample(nodes_within_bin$entrezID, size = sum(nodes_within_bin$original_entrezID))
    # annotate the permuted shared genes
    nodes_within_bin$perm_shared_genes = ifelse(nodes_within_bin$entrezID %in% random_seed_nodes, 1, 0)
    new_shared_genes = nodes_within_bin %>% filter(original_entrezID == 1) %>% dplyr::select(original_shared_entrezID = entrezID)
    new_shared_genes = cbind(new_shared_genes,
                             nodes_within_bin %>% filter(perm_shared_genes == 1) %>% dplyr::select(perm_shared_entrezID = entrezID))
    perm_bc_ldl_shared_genes = rbind(perm_bc_ldl_shared_genes, new_shared_genes)
  }
  perm_bc_ldl_shared_genes = perm_bc_ldl_shared_genes[-1, ] ; rownames(perm_bc_ldl_shared_genes) = NULL
  
  ## --------------------- parallelize the PPR runs for each canonical pathway per gene for each permutation ---------------------- ##
  cp_perm_genes = vector(mode = "list", length = length(msigdb))
  names(cp_perm_genes) = names(msigdb)
  
  cp_perm_genes = foreach(cp = 1:length(msigdb), .packages = c("igraph", "dplyr")) %dopar% {
    genes_in_cp = msigdb[[cp]]
    perm_avg_ppr_score_per_shared_gene = data.frame(shared_gene = unique(perm_bc_ldl_shared_genes$perm_shared_entrezID),
                                                    average_ppr_score_cp = 0,
                                                    cp = names(msigdb)[cp])
    for (i in 1:nrow(perm_avg_ppr_score_per_shared_gene)) {
      shared_gene_temp = as.character(perm_avg_ppr_score_per_shared_gene[i, "shared_gene"])
      p = c(rep(0, length(V(g))))
      p[which(V(g)$name == shared_gene_temp)] = 1
      pr_result = page_rank(g, algo = "prpack", damping = 0.85, personalized = p, weights = NA)
      pr_result = pr_result$vector
      pr_result = as.data.frame(pr_result)
      pr_result$entrezID = rownames(pr_result)
      rownames(pr_result) = NULL
      pr_result = pr_result %>%
        dplyr::select(entrezID, propagation_score = pr_result) %>%
        filter(entrezID %in% genes_in_cp) %>%
        mutate(avg_propagation_score = mean(propagation_score)) %>%
        dplyr::select(avg_propagation_score) %>%
        distinct()
      perm_avg_ppr_score_per_shared_gene[i, "average_ppr_score_cp"] = pr_result
    }
    perm_avg_ppr_score_per_shared_gene = left_join(perm_avg_ppr_score_per_shared_gene, perm_bc_ldl_shared_genes, by = c("shared_gene" = "perm_shared_entrezID"))
    perm_avg_ppr_score_per_shared_gene = perm_avg_ppr_score_per_shared_gene %>%
      dplyr::select(original_shared_entrezID, perm_shared_gene = shared_gene, average_ppr_score_cp_perm_gene = average_ppr_score_cp, cp)
    return(perm_avg_ppr_score_per_shared_gene)
  }
  names(cp_perm_genes) = names(msigdb)
  
  # populate list with the results from one permutation
  avg_ppr_perm[[perm]] = cp_perm_genes
  
  # track progress
  cat(perm, "- 1000", "\n")
}
saveRDS(avg_ppr_perm, "shared_genes_to_canonical_pathways/permuted_BC_LDL.RDS")

## BC-Prostate cancer ##
avg_ppr_perm = vector(mode = "list", length = 1000)
names(avg_ppr_perm) = paste0("permutation_", 1:1000)
# Set up parallel processing
cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) # detectCores() doesn't work with the way the Slurm limits the number of cores available.
print(paste0("Number of cores used for parallelization: ", cores))
registerDoParallel(cores)
# Parallelized loop
for (perm in 1:length(avg_ppr_perm)) {
  ## -------------------------------------- shuffle shared genes based on degree bins --------------------------------------------- ##
  perm_bc_prostate_shared_genes = data.frame(original_shared_entrezID = NA, perm_shared_entrezID = NA)
  for (i in 1:length(bins)) {
    # find entrezIDs in this degree bin
    nodes_within_bin = as.data.frame(names(node_deg[which(node_deg %in% bins[[i]])]))
    colnames(nodes_within_bin) = "entrezID"
    # annotate the ones that are shared genes
    nodes_within_bin$original_entrezID = ifelse(nodes_within_bin$entrezID %in% bc_prostate_shared_genes$entrezID, 1, 0)
    # shuffle the shared genes
    random_seed_nodes = sample(nodes_within_bin$entrezID, size = sum(nodes_within_bin$original_entrezID))
    # annotate the permuted shared genes
    nodes_within_bin$perm_shared_genes = ifelse(nodes_within_bin$entrezID %in% random_seed_nodes, 1, 0)
    new_shared_genes = nodes_within_bin %>% filter(original_entrezID == 1) %>% dplyr::select(original_shared_entrezID = entrezID)
    new_shared_genes = cbind(new_shared_genes,
                             nodes_within_bin %>% filter(perm_shared_genes == 1) %>% dplyr::select(perm_shared_entrezID = entrezID))
    perm_bc_prostate_shared_genes = rbind(perm_bc_prostate_shared_genes, new_shared_genes)
  }
  perm_bc_prostate_shared_genes = perm_bc_prostate_shared_genes[-1, ] ; rownames(perm_bc_prostate_shared_genes) = NULL
  ## --------------------- parallelize the PPR runs for each canonical pathway per gene for each permutation ---------------------- ##
  cp_perm_genes = vector(mode = "list", length = length(msigdb))
  names(cp_perm_genes) = names(msigdb)
  
  cp_perm_genes = foreach(cp = 1:length(msigdb), .packages = c("igraph", "dplyr")) %dopar% {
    genes_in_cp = msigdb[[cp]]
    perm_avg_ppr_score_per_shared_gene = data.frame(shared_gene = unique(perm_bc_prostate_shared_genes$perm_shared_entrezID),
                                                    average_ppr_score_cp = 0,
                                                    cp = names(msigdb)[cp])
    for (i in 1:nrow(perm_avg_ppr_score_per_shared_gene)) {
      shared_gene_temp = as.character(perm_avg_ppr_score_per_shared_gene[i, "shared_gene"])
      p = c(rep(0, length(V(g))))
      p[which(V(g)$name == shared_gene_temp)] = 1
      pr_result = page_rank(g, algo = "prpack", damping = 0.85, personalized = p, weights = NA)
      pr_result = pr_result$vector
      pr_result = as.data.frame(pr_result)
      pr_result$entrezID = rownames(pr_result)
      rownames(pr_result) = NULL
      pr_result = pr_result %>%
        dplyr::select(entrezID, propagation_score = pr_result) %>%
        filter(entrezID %in% genes_in_cp) %>%
        mutate(avg_propagation_score = mean(propagation_score)) %>%
        dplyr::select(avg_propagation_score) %>%
        distinct()
      perm_avg_ppr_score_per_shared_gene[i, "average_ppr_score_cp"] = pr_result
    }
    perm_avg_ppr_score_per_shared_gene = left_join(perm_avg_ppr_score_per_shared_gene, perm_bc_prostate_shared_genes, by = c("shared_gene" = "perm_shared_entrezID"))
    perm_avg_ppr_score_per_shared_gene = perm_avg_ppr_score_per_shared_gene %>%
      dplyr::select(original_shared_entrezID, perm_shared_gene = shared_gene, average_ppr_score_cp_perm_gene = average_ppr_score_cp, cp)
    return(perm_avg_ppr_score_per_shared_gene)
  }
  names(cp_perm_genes) = names(msigdb)
  # populate list with the results from one permutation
  avg_ppr_perm[[perm]] = cp_perm_genes
  # track progress
  cat(perm, "- 1000", "\n")
}
saveRDS(avg_ppr_perm, "shared_genes_to_canonical_pathways/permuted_BC_PROSTATE.RDS")

## BC-Schizophrenia ##
avg_ppr_perm = vector(mode = "list", length = 1000)
names(avg_ppr_perm) = paste0("permutation_", 1:1000)
# Set up parallel processing
cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) # detectCores() doesn't work with the way the Slurm limits the number of cores available.
print(paste0("Number of cores used for parallelization: ", cores))
registerDoParallel(cores)
# Parallelized loop
for (perm in 1:length(avg_ppr_perm)) {
  ## -------------------------------------- shuffle shared genes based on degree bins --------------------------------------------- ##
  perm_bc_scz_shared_genes = data.frame(original_shared_entrezID = NA, perm_shared_entrezID = NA)
  for (i in 1:length(bins)) {
    # find entrezIDs in this degree bin
    nodes_within_bin = as.data.frame(names(node_deg[which(node_deg %in% bins[[i]])]))
    colnames(nodes_within_bin) = "entrezID"
    # annotate the ones that are shared genes
    nodes_within_bin$original_entrezID = ifelse(nodes_within_bin$entrezID %in% bc_scz_shared_genes$entrezID, 1, 0)
    # shuffle the shared genes
    random_seed_nodes = sample(nodes_within_bin$entrezID, size = sum(nodes_within_bin$original_entrezID))
    # annotate the permuted shared genes
    nodes_within_bin$perm_shared_genes = ifelse(nodes_within_bin$entrezID %in% random_seed_nodes, 1, 0)
    new_shared_genes = nodes_within_bin %>% filter(original_entrezID == 1) %>% dplyr::select(original_shared_entrezID = entrezID)
    new_shared_genes = cbind(new_shared_genes,
                             nodes_within_bin %>% filter(perm_shared_genes == 1) %>% dplyr::select(perm_shared_entrezID = entrezID))
    perm_bc_scz_shared_genes = rbind(perm_bc_scz_shared_genes, new_shared_genes)
  }
  perm_bc_scz_shared_genes = perm_bc_scz_shared_genes[-1, ] ; rownames(perm_bc_scz_shared_genes) = NULL
  ## --------------------- parallelize the PPR runs for each canonical pathway per gene for each permutation ---------------------- ##
  cp_perm_genes = vector(mode = "list", length = length(msigdb))
  names(cp_perm_genes) = names(msigdb)
  
  cp_perm_genes = foreach(cp = 1:length(msigdb), .packages = c("igraph", "dplyr")) %dopar% {
    genes_in_cp = msigdb[[cp]]
    perm_avg_ppr_score_per_shared_gene = data.frame(shared_gene = unique(perm_bc_scz_shared_genes$perm_shared_entrezID),
                                                    average_ppr_score_cp = 0,
                                                    cp = names(msigdb)[cp])
    for (i in 1:nrow(perm_avg_ppr_score_per_shared_gene)) {
      shared_gene_temp = as.character(perm_avg_ppr_score_per_shared_gene[i, "shared_gene"])
      p = c(rep(0, length(V(g))))
      p[which(V(g)$name == shared_gene_temp)] = 1
      pr_result = page_rank(g, algo = "prpack", damping = 0.85, personalized = p, weights = NA)
      pr_result = pr_result$vector
      pr_result = as.data.frame(pr_result)
      pr_result$entrezID = rownames(pr_result)
      rownames(pr_result) = NULL
      pr_result = pr_result %>%
        dplyr::select(entrezID, propagation_score = pr_result) %>%
        filter(entrezID %in% genes_in_cp) %>%
        mutate(avg_propagation_score = mean(propagation_score)) %>%
        dplyr::select(avg_propagation_score) %>%
        distinct()
      perm_avg_ppr_score_per_shared_gene[i, "average_ppr_score_cp"] = pr_result
    }
    perm_avg_ppr_score_per_shared_gene = left_join(perm_avg_ppr_score_per_shared_gene, perm_bc_scz_shared_genes, by = c("shared_gene" = "perm_shared_entrezID"))
    perm_avg_ppr_score_per_shared_gene = perm_avg_ppr_score_per_shared_gene %>%
      dplyr::select(original_shared_entrezID, perm_shared_gene = shared_gene, average_ppr_score_cp_perm_gene = average_ppr_score_cp, cp)
    return(perm_avg_ppr_score_per_shared_gene)
  }
  names(cp_perm_genes) = names(msigdb)
  # populate list with the results from one permutation
  avg_ppr_perm[[perm]] = cp_perm_genes
  # track progress
  cat(perm, "- 1000", "\n")
}
saveRDS(avg_ppr_perm, "shared_genes_to_canonical_pathways/permuted_BC_SCZ.RDS")

## BC-T2DM ##
avg_ppr_perm = vector(mode = "list", length = 1000)
names(avg_ppr_perm) = paste0("permutation_", 1:1000)
# Set up parallel processing
cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) # detectCores() doesn't work with the way the Slurm limits the number of cores available.
print(paste0("Number of cores used for parallelization: ", cores))
registerDoParallel(cores)
# Parallelized loop
for (perm in 1:length(avg_ppr_perm)) {
  ## -------------------------------------- shuffle shared genes based on degree bins --------------------------------------------- ##
  perm_bc_t2dm_shared_genes = data.frame(original_shared_entrezID = NA, perm_shared_entrezID = NA)
  for (i in 1:length(bins)) {
    # find entrezIDs in this degree bin
    nodes_within_bin = as.data.frame(names(node_deg[which(node_deg %in% bins[[i]])]))
    colnames(nodes_within_bin) = "entrezID"
    # annotate the ones that are shared genes
    nodes_within_bin$original_entrezID = ifelse(nodes_within_bin$entrezID %in% bc_t2dm_shared_genes$entrezID, 1, 0)
    # shuffle the shared genes
    random_seed_nodes = sample(nodes_within_bin$entrezID, size = sum(nodes_within_bin$original_entrezID))
    # annotate the permuted shared genes
    nodes_within_bin$perm_shared_genes = ifelse(nodes_within_bin$entrezID %in% random_seed_nodes, 1, 0)
    new_shared_genes = nodes_within_bin %>% filter(original_entrezID == 1) %>% dplyr::select(original_shared_entrezID = entrezID)
    new_shared_genes = cbind(new_shared_genes,
                             nodes_within_bin %>% filter(perm_shared_genes == 1) %>% dplyr::select(perm_shared_entrezID = entrezID))
    perm_bc_t2dm_shared_genes = rbind(perm_bc_t2dm_shared_genes, new_shared_genes)
  }
  perm_bc_t2dm_shared_genes = perm_bc_t2dm_shared_genes[-1, ] ; rownames(perm_bc_t2dm_shared_genes) = NULL
  ## --------------------- parallelize the PPR runs for each canonical pathway per gene for each permutation ---------------------- ##
  cp_perm_genes = vector(mode = "list", length = length(msigdb))
  names(cp_perm_genes) = names(msigdb)
  
  cp_perm_genes = foreach(cp = 1:length(msigdb), .packages = c("igraph", "dplyr")) %dopar% {
    genes_in_cp = msigdb[[cp]]
    perm_avg_ppr_score_per_shared_gene = data.frame(shared_gene = unique(perm_bc_t2dm_shared_genes$perm_shared_entrezID),
                                                    average_ppr_score_cp = 0,
                                                    cp = names(msigdb)[cp])
    for (i in 1:nrow(perm_avg_ppr_score_per_shared_gene)) {
      shared_gene_temp = as.character(perm_avg_ppr_score_per_shared_gene[i, "shared_gene"])
      p = c(rep(0, length(V(g))))
      p[which(V(g)$name == shared_gene_temp)] = 1
      pr_result = page_rank(g, algo = "prpack", damping = 0.85, personalized = p, weights = NA)
      pr_result = pr_result$vector
      pr_result = as.data.frame(pr_result)
      pr_result$entrezID = rownames(pr_result)
      rownames(pr_result) = NULL
      pr_result = pr_result %>%
        dplyr::select(entrezID, propagation_score = pr_result) %>%
        filter(entrezID %in% genes_in_cp) %>%
        mutate(avg_propagation_score = mean(propagation_score)) %>%
        dplyr::select(avg_propagation_score) %>%
        distinct()
      perm_avg_ppr_score_per_shared_gene[i, "average_ppr_score_cp"] = pr_result
    }
    perm_avg_ppr_score_per_shared_gene = left_join(perm_avg_ppr_score_per_shared_gene, perm_bc_t2dm_shared_genes, by = c("shared_gene" = "perm_shared_entrezID"))
    perm_avg_ppr_score_per_shared_gene = perm_avg_ppr_score_per_shared_gene %>%
      dplyr::select(original_shared_entrezID, perm_shared_gene = shared_gene, average_ppr_score_cp_perm_gene = average_ppr_score_cp, cp)
    return(perm_avg_ppr_score_per_shared_gene)
  }
  names(cp_perm_genes) = names(msigdb)
  # populate list with the results from one permutation
  avg_ppr_perm[[perm]] = cp_perm_genes
  # track progress
  cat(perm, "- 1000", "\n")
}
saveRDS(avg_ppr_perm, "shared_genes_to_canonical_pathways/permuted_BC_T2DM_DIAMANTE.RDS")

rm(list = ls())
gc()

##### CALCULATING P-VALUES FOR EACH CANONICAL PATHWAY, FOR EACH DISEASE PAIR #####
## Adjusting for the number of shared genes and the number of canonical pathways tested (bonferonni)

## BC-DEPRESSION
# observed
bc_depression_observed = readRDS("shared_genes_to_canonical_pathways/observed_BC_DEPRESSION.RDS")
# permuted
bc_depression_permuted = readRDS("shared_genes_to_canonical_pathways/permuted_BC_DEPRESSION.RDS")
# create list to populate later with permuted pvalues
bc_depression_cp_pvalues = vector("list", length(bc_depression_observed))
names(bc_depression_cp_pvalues) = names(bc_depression_observed)
for (i in 1:length(bc_depression_cp_pvalues)) {
  # For a canonical pathway...
  cp_temp = names(bc_depression_cp_pvalues)[i]
  ## observed values for that canonical pathway
  obs_value = bc_depression_observed[[cp_temp]]
  # Data frame with permuted p-values --> populate later
  perm_pvalues = data.frame(shared_gene = obs_value$shared_gene, perm_pvalue = 99)
  # for each shared gene...
  for (z in 1:nrow(perm_pvalues)) {
    # find the shared gene
    shared_gene_temp = perm_pvalues[z, "shared_gene"]
    # observed PPR score for that gene
    obs_value_temp = obs_value %>% filter(shared_gene == shared_gene_temp)
    # permuted PPR scores for that gene
    perm_value = sapply(bc_depression_permuted, function(perm) {
      x = perm[[cp_temp]] %>% filter(original_shared_entrezID == shared_gene_temp)
      x$average_ppr_score_cp_perm_gene
    })
    perm_pvalues[z , "perm_pvalue"] = sum(obs_value_temp$average_ppr_score_cp <= perm_value) / 1000
  }
  bc_depression_cp_pvalues[[i]] = perm_pvalues
} ; rm(i, cp_temp, obs_value, perm_pvalues, z, shared_gene_temp, obs_value_temp)
# canonical pathways significantly connected to the shared genes
bc_depression_sig_cp = bc_depression_cp_pvalues
for (i in 1:length(bc_depression_sig_cp)) {
  bc_depression_sig_cp[[i]] = bc_depression_sig_cp[[i]] %>%
    # adjust for the number of shared genes tested for each canonical pathway
    mutate(ppr_pvalue_adj = p.adjust(perm_pvalue, "bonferroni", length(perm_pvalue))) %>%
    mutate(ppr_pvalue_adj = min(ppr_pvalue_adj)) %>%
    dplyr::select(ppr_pvalue_adj) %>%
    distinct()
  cat(i, "\n")
} ; rm(i)
bc_depression_sig_cp = as.data.frame(unlist(bc_depression_sig_cp))
bc_depression_sig_cp$canonical_pathway = rownames(bc_depression_sig_cp) ; rownames(bc_depression_sig_cp) = NULL
bc_depression_sig_cp$canonical_pathway = gsub(".ppr_pvalue_adj", "", bc_depression_sig_cp$canonical_pathway)
bc_depression_sig_cp = bc_depression_sig_cp %>% dplyr::select(canonical_pathway, p_adjusted = "unlist(bc_depression_sig_cp)") %>% distinct()
bc_depression_sig_cp$sig = ifelse(bc_depression_sig_cp$p_adjusted < 0.05, 1, 0) 
# adjust for the number of canonical pathways tested
bc_depression_sig_cp$p_adjusted_final = p.adjust(bc_depression_sig_cp$p_adjusted, method = "bonferroni")
bc_depression_sig_cp = bc_depression_sig_cp %>%
  dplyr::select(canonical_pathway, pvalue_adj = p_adjusted_final) %>%
  distinct()
bc_depression_sig_cp$sig = ifelse(bc_depression_sig_cp$pvalue_adj < 0.05, 1, 0)
fwrite(bc_depression_sig_cp, "shared_genes_to_canonical_pathways/BC_DEPRESSION_CPs_pvalues_bonferroni.txt", sep = "\t", row.names = FALSE)

## BC-HDL
# observed
bc_hdl_observed = readRDS("shared_genes_to_canonical_pathways/observed_BC_HDL.RDS")
# permuted
bc_hdl_permuted = readRDS("shared_genes_to_canonical_pathways/permuted_BC_HDL.RDS")
# create list to populate later with permuted pvalues
bc_hdl_cp_pvalues = vector("list", length(bc_hdl_observed))
names(bc_hdl_cp_pvalues) = names(bc_hdl_observed)
for (i in 1:length(bc_hdl_cp_pvalues)) {
  # For a canonical pathway...
  cp_temp = names(bc_hdl_cp_pvalues)[i]
  ## observed values for that canonical pathway
  obs_value = bc_hdl_observed[[cp_temp]]
  # Data frame with permuted p-values --> populate later
  perm_pvalues = data.frame(shared_gene = obs_value$shared_gene, perm_pvalue = 99)
  # for each shared gene...
  for (z in 1:nrow(perm_pvalues)) {
    # find the shared gene
    shared_gene_temp = perm_pvalues[z, "shared_gene"]
    # observed PPR score for that gene
    obs_value_temp = obs_value %>% filter(shared_gene == shared_gene_temp)
    # permuted PPR scores for that gene
    perm_value = sapply(bc_hdl_permuted, function(perm) {
      x = perm[[cp_temp]] %>% filter(original_shared_entrezID == shared_gene_temp)
      x$average_ppr_score_cp_perm_gene
    })
    perm_pvalues[z , "perm_pvalue"] = sum(obs_value_temp$average_ppr_score_cp <= perm_value) / 1000
  }
  bc_hdl_cp_pvalues[[i]] = perm_pvalues
} ; rm(i, cp_temp, obs_value, perm_pvalues, z, shared_gene_temp, obs_value_temp)

# canonical pathways significantly connected to the shared genes
bc_hdl_sig_cp = bc_hdl_cp_pvalues
for (i in 1:length(bc_hdl_sig_cp)) {
  bc_hdl_sig_cp[[i]] = bc_hdl_sig_cp[[i]] %>%
    # adjust for the number of shared genes tested for each canonical pathway
    mutate(ppr_pvalue_adj = p.adjust(perm_pvalue, "bonferroni", length(perm_pvalue))) %>%
    mutate(ppr_pvalue_adj = min(ppr_pvalue_adj)) %>%
    dplyr::select(ppr_pvalue_adj) %>%
    distinct()
  cat(i, "\n")
} ; rm(i)
bc_hdl_sig_cp = as.data.frame(unlist(bc_hdl_sig_cp))
bc_hdl_sig_cp$canonical_pathway = rownames(bc_hdl_sig_cp) ; rownames(bc_hdl_sig_cp) = NULL
bc_hdl_sig_cp$canonical_pathway = gsub(".ppr_pvalue_adj", "", bc_hdl_sig_cp$canonical_pathway)
bc_hdl_sig_cp = bc_hdl_sig_cp %>% dplyr::select(canonical_pathway, p_adjusted = "unlist(bc_hdl_sig_cp)") %>% distinct()
bc_hdl_sig_cp$sig = ifelse(bc_hdl_sig_cp$p_adjusted < 0.05, 1, 0) 
# adjust for the number of canonical pathways tested
bc_hdl_sig_cp$p_adjusted_final = p.adjust(bc_hdl_sig_cp$p_adjusted, method = "bonferroni")
bc_hdl_sig_cp = bc_hdl_sig_cp %>%
  dplyr::select(canonical_pathway, pvalue_adj = p_adjusted_final) %>%
  distinct()
bc_hdl_sig_cp$sig = ifelse(bc_hdl_sig_cp$pvalue_adj < 0.05, 1, 0)
fwrite(bc_hdl_sig_cp, "shared_genes_to_canonical_pathways/BC_HDL_CPs_pvalues_bonferroni.txt", sep = "\t", row.names = FALSE)

## BC-LDL
# observed
bc_ldl_observed = readRDS("shared_genes_to_canonical_pathways/observed_BC_LDL.RDS")
# permuted
bc_ldl_permuted = readRDS("shared_genes_to_canonical_pathways/permuted_BC_LDL.RDS")
# create list to populate later with permuted pvalues
bc_ldl_cp_pvalues = vector("list", length(bc_ldl_observed))
names(bc_ldl_cp_pvalues) = names(bc_ldl_observed)
for (i in 1:length(bc_ldl_cp_pvalues)) {
  # For a canonical pathway...
  cp_temp = names(bc_ldl_cp_pvalues)[i]
  ## observed values for that canonical pathway
  obs_value = bc_ldl_observed[[cp_temp]]
  # Data frame with permuted p-values --> populate later
  perm_pvalues = data.frame(shared_gene = obs_value$shared_gene, perm_pvalue = 99)
  # for each shared gene...
  for (z in 1:nrow(perm_pvalues)) {
    # find the shared gene
    shared_gene_temp = perm_pvalues[z, "shared_gene"]
    # observed PPR score for that gene
    obs_value_temp = obs_value %>% filter(shared_gene == shared_gene_temp)
    # permuted PPR scores for that gene
    perm_value = sapply(bc_ldl_permuted, function(perm) {
      x = perm[[cp_temp]] %>% filter(original_shared_entrezID == shared_gene_temp)
      x$average_ppr_score_cp_perm_gene
    })
    perm_pvalues[z , "perm_pvalue"] = sum(obs_value_temp$average_ppr_score_cp <= perm_value) / 1000
  }
  bc_ldl_cp_pvalues[[i]] = perm_pvalues
} ; rm(i, cp_temp, obs_value, perm_pvalues, z, shared_gene_temp, obs_value_temp)

# canonical pathways significantly connected to the shared genes
bc_ldl_sig_cp = bc_ldl_cp_pvalues
for (i in 1:length(bc_ldl_sig_cp)) {
  bc_ldl_sig_cp[[i]] = bc_ldl_sig_cp[[i]] %>%
    # adjust for the number of shared genes tested for each canonical pathway
    mutate(ppr_pvalue_adj = p.adjust(perm_pvalue, "bonferroni", length(perm_pvalue))) %>%
    mutate(ppr_pvalue_adj = min(ppr_pvalue_adj)) %>%
    dplyr::select(ppr_pvalue_adj) %>%
    distinct()
  cat(i, "\n")
} ; rm(i)
bc_ldl_sig_cp = as.data.frame(unlist(bc_ldl_sig_cp))
bc_ldl_sig_cp$canonical_pathway = rownames(bc_ldl_sig_cp) ; rownames(bc_ldl_sig_cp) = NULL
bc_ldl_sig_cp$canonical_pathway = gsub(".ppr_pvalue_adj", "", bc_ldl_sig_cp$canonical_pathway)
bc_ldl_sig_cp = bc_ldl_sig_cp %>% dplyr::select(canonical_pathway, p_adjusted = "unlist(bc_ldl_sig_cp)") %>% distinct()
bc_ldl_sig_cp$sig = ifelse(bc_ldl_sig_cp$p_adjusted < 0.05, 1, 0) 
# adjust for the number of canonical pathways tested
bc_ldl_sig_cp$p_adjusted_final = p.adjust(bc_ldl_sig_cp$p_adjusted, method = "bonferroni")
bc_ldl_sig_cp = bc_ldl_sig_cp %>%
  dplyr::select(canonical_pathway, pvalue_adj = p_adjusted_final) %>%
  distinct()
bc_ldl_sig_cp$sig = ifelse(bc_ldl_sig_cp$pvalue_adj < 0.05, 1, 0)
fwrite(bc_ldl_sig_cp, "shared_genes_to_canonical_pathways/BC_LDL_CPs_pvalues_bonferroni.txt", sep = "\t", row.names = FALSE)

## BC-PROSTATE
# observed
bc_prostate_observed = readRDS("shared_genes_to_canonical_pathways/observed_BC_PROSTATE.RDS")
# permuted
bc_prostate_permuted = readRDS("shared_genes_to_canonical_pathways/permuted_BC_PROSTATE.RDS")
# create list to populate later with permuted pvalues
bc_prostate_cp_pvalues = vector("list", length(bc_prostate_observed))
names(bc_prostate_cp_pvalues) = names(bc_prostate_observed)
for (i in 1:length(bc_prostate_cp_pvalues)) {
  # For a canonical pathway...
  cp_temp = names(bc_prostate_cp_pvalues)[i]
  ## observed values for that canonical pathway
  obs_value = bc_prostate_observed[[cp_temp]]
  # Data frame with permuted p-values --> populate later
  perm_pvalues = data.frame(shared_gene = obs_value$shared_gene, perm_pvalue = 99)
  # for each shared gene...
  for (z in 1:nrow(perm_pvalues)) {
    # find the shared gene
    shared_gene_temp = perm_pvalues[z, "shared_gene"]
    # observed PPR score for that gene
    obs_value_temp = obs_value %>% filter(shared_gene == shared_gene_temp)
    # permuted PPR scores for that gene
    perm_value = sapply(bc_prostate_permuted, function(perm) {
      x = perm[[cp_temp]] %>% filter(original_shared_entrezID == shared_gene_temp)
      x$average_ppr_score_cp_perm_gene
    })
    perm_pvalues[z , "perm_pvalue"] = sum(obs_value_temp$average_ppr_score_cp <= perm_value) / 1000
  }
  bc_prostate_cp_pvalues[[i]] = perm_pvalues
} ; rm(i, cp_temp, obs_value, perm_pvalues, z, shared_gene_temp, obs_value_temp)

# canonical pathways significantly connected to the shared genes
bc_prostate_sig_cp = bc_prostate_cp_pvalues
for (i in 1:length(bc_prostate_sig_cp)) {
  bc_prostate_sig_cp[[i]] = bc_prostate_sig_cp[[i]] %>%
    # adjust for the number of shared genes tested for each canonical pathway
    mutate(ppr_pvalue_adj = p.adjust(perm_pvalue, "bonferroni", length(perm_pvalue))) %>%
    mutate(ppr_pvalue_adj = min(ppr_pvalue_adj)) %>%
    dplyr::select(ppr_pvalue_adj) %>%
    distinct()
  cat(i, "\n")
} ; rm(i)
bc_prostate_sig_cp = as.data.frame(unlist(bc_prostate_sig_cp))
bc_prostate_sig_cp$canonical_pathway = rownames(bc_prostate_sig_cp) ; rownames(bc_prostate_sig_cp) = NULL
bc_prostate_sig_cp$canonical_pathway = gsub(".ppr_pvalue_adj", "", bc_prostate_sig_cp$canonical_pathway)
bc_prostate_sig_cp = bc_prostate_sig_cp %>% dplyr::select(canonical_pathway, p_adjusted = "unlist(bc_prostate_sig_cp)") %>% distinct()
bc_prostate_sig_cp$sig = ifelse(bc_prostate_sig_cp$p_adjusted < 0.05, 1, 0) 
# adjust for the number of canonical pathways tested
bc_prostate_sig_cp$p_adjusted_final = p.adjust(bc_prostate_sig_cp$p_adjusted, method = "bonferroni")
bc_prostate_sig_cp = bc_prostate_sig_cp %>%
  dplyr::select(canonical_pathway, pvalue_adj = p_adjusted_final) %>%
  distinct()
bc_prostate_sig_cp$sig = ifelse(bc_prostate_sig_cp$pvalue_adj < 0.05, 1, 0)
fwrite(bc_prostate_sig_cp, "shared_genes_to_canonical_pathways/BC_PROSTATE_CPs_pvalues_bonferroni.txt", sep = "\t", row.names = FALSE)

## BC-SCZ
# observed
bc_scz_observed = readRDS("shared_genes_to_canonical_pathways/observed_BC_SCZ.RDS")
# permuted
bc_scz_permuted = readRDS("shared_genes_to_canonical_pathways/permuted_BC_SCZ.RDS")
# create list to populate later with permuted pvalues
bc_scz_cp_pvalues = vector("list", length(bc_scz_observed))
names(bc_scz_cp_pvalues) = names(bc_scz_observed)
for (i in 1:length(bc_scz_cp_pvalues)) {
  # For a canonical pathway...
  cp_temp = names(bc_scz_cp_pvalues)[i]
  ## observed values for that canonical pathway
  obs_value = bc_scz_observed[[cp_temp]]
  # Data frame with permuted p-values --> populate later
  perm_pvalues = data.frame(shared_gene = obs_value$shared_gene, perm_pvalue = 99)
  # for each shared gene...
  for (z in 1:nrow(perm_pvalues)) {
    # find the shared gene
    shared_gene_temp = perm_pvalues[z, "shared_gene"]
    # observed PPR score for that gene
    obs_value_temp = obs_value %>% filter(shared_gene == shared_gene_temp)
    # permuted PPR scores for that gene
    perm_value = sapply(bc_scz_permuted, function(perm) {
      x = perm[[cp_temp]] %>% filter(original_shared_entrezID == shared_gene_temp)
      x$average_ppr_score_cp_perm_gene
    })
    perm_pvalues[z , "perm_pvalue"] = sum(obs_value_temp$average_ppr_score_cp <= perm_value) / 1000
  }
  bc_scz_cp_pvalues[[i]] = perm_pvalues
} ; rm(i, cp_temp, obs_value, perm_pvalues, z, shared_gene_temp, obs_value_temp)

# canonical pathways significantly connected to the shared genes
bc_scz_sig_cp = bc_scz_cp_pvalues
for (i in 1:length(bc_scz_sig_cp)) {
  bc_scz_sig_cp[[i]] = bc_scz_sig_cp[[i]] %>%
    # adjust for the number of shared genes tested for each canonical pathway
    mutate(ppr_pvalue_adj = p.adjust(perm_pvalue, "bonferroni", length(perm_pvalue))) %>%
    mutate(ppr_pvalue_adj = min(ppr_pvalue_adj)) %>%
    dplyr::select(ppr_pvalue_adj) %>%
    distinct()
  cat(i, "\n")
} ; rm(i)
bc_scz_sig_cp = as.data.frame(unlist(bc_scz_sig_cp))
bc_scz_sig_cp$canonical_pathway = rownames(bc_scz_sig_cp) ; rownames(bc_scz_sig_cp) = NULL
bc_scz_sig_cp$canonical_pathway = gsub(".ppr_pvalue_adj", "", bc_scz_sig_cp$canonical_pathway)
bc_scz_sig_cp = bc_scz_sig_cp %>% dplyr::select(canonical_pathway, p_adjusted = "unlist(bc_scz_sig_cp)") %>% distinct()
bc_scz_sig_cp$sig = ifelse(bc_scz_sig_cp$p_adjusted < 0.05, 1, 0) 
# adjust for the number of canonical pathways tested
bc_scz_sig_cp$p_adjusted_final = p.adjust(bc_scz_sig_cp$p_adjusted, method = "bonferroni")
bc_scz_sig_cp = bc_scz_sig_cp %>%
  dplyr::select(canonical_pathway, pvalue_adj = p_adjusted_final) %>%
  distinct()
bc_scz_sig_cp$sig = ifelse(bc_scz_sig_cp$pvalue_adj < 0.05, 1, 0)
fwrite(bc_scz_sig_cp, "shared_genes_to_canonical_pathways/BC_SCZ_CPs_pvalues_bonferroni.txt", sep = "\t", row.names = FALSE)

## BC-T2DM_DIAMANTE
# observed
bc_t2dm_diamante_observed = readRDS("shared_genes_to_canonical_pathways/observed_BC_T2DM_DIAMANTE.RDS")
# permuted
bc_t2dm_diamante_permuted = readRDS("shared_genes_to_canonical_pathways/permuted_BC_T2DM_DIAMANTE.RDS")
# create list to populate later with permuted pvalues
bc_t2dm_diamante_cp_pvalues = vector("list", length(bc_t2dm_diamante_observed))
names(bc_t2dm_diamante_cp_pvalues) = names(bc_t2dm_diamante_observed)
for (i in 1:length(bc_t2dm_diamante_cp_pvalues)) {
  # For a canonical pathway...
  cp_temp = names(bc_t2dm_diamante_cp_pvalues)[i]
  ## observed values for that canonical pathway
  obs_value = bc_t2dm_diamante_observed[[cp_temp]]
  # Data frame with permuted p-values --> populate later
  perm_pvalues = data.frame(shared_gene = obs_value$shared_gene, perm_pvalue = 99)
  # for each shared gene...
  for (z in 1:nrow(perm_pvalues)) {
    # find the shared gene
    shared_gene_temp = perm_pvalues[z, "shared_gene"]
    # observed PPR score for that gene
    obs_value_temp = obs_value %>% filter(shared_gene == shared_gene_temp)
    # permuted PPR scores for that gene
    perm_value = sapply(bc_t2dm_diamante_permuted, function(perm) {
      x = perm[[cp_temp]] %>% filter(original_shared_entrezID == shared_gene_temp)
      x$average_ppr_score_cp_perm_gene
    })
    perm_pvalues[z , "perm_pvalue"] = sum(obs_value_temp$average_ppr_score_cp <= perm_value) / 1000
  }
  bc_t2dm_diamante_cp_pvalues[[i]] = perm_pvalues
} ; rm(i, cp_temp, obs_value, perm_pvalues, z, shared_gene_temp, obs_value_temp)

# canonical pathways significantly connected to the shared genes
bc_t2dm_diamante_sig_cp = bc_t2dm_diamante_cp_pvalues
for (i in 1:length(bc_t2dm_diamante_sig_cp)) {
  bc_t2dm_diamante_sig_cp[[i]] = bc_t2dm_diamante_sig_cp[[i]] %>%
    # adjust for the number of shared genes tested for each canonical pathway
    mutate(ppr_pvalue_adj = p.adjust(perm_pvalue, "bonferroni", length(perm_pvalue))) %>%
    mutate(ppr_pvalue_adj = min(ppr_pvalue_adj)) %>%
    dplyr::select(ppr_pvalue_adj) %>%
    distinct()
  cat(i, "\n")
} ; rm(i)
bc_t2dm_diamante_sig_cp = as.data.frame(unlist(bc_t2dm_diamante_sig_cp))
bc_t2dm_diamante_sig_cp$canonical_pathway = rownames(bc_t2dm_diamante_sig_cp) ; rownames(bc_t2dm_diamante_sig_cp) = NULL
bc_t2dm_diamante_sig_cp$canonical_pathway = gsub(".ppr_pvalue_adj", "", bc_t2dm_diamante_sig_cp$canonical_pathway)
bc_t2dm_diamante_sig_cp = bc_t2dm_diamante_sig_cp %>% dplyr::select(canonical_pathway, p_adjusted = "unlist(bc_t2dm_diamante_sig_cp)") %>% distinct()
bc_t2dm_diamante_sig_cp$sig = ifelse(bc_t2dm_diamante_sig_cp$p_adjusted < 0.05, 1, 0) 
# adjust for the number of canonical pathways tested
bc_t2dm_diamante_sig_cp$p_adjusted_final = p.adjust(bc_t2dm_diamante_sig_cp$p_adjusted, method = "bonferroni")
bc_t2dm_diamante_sig_cp = bc_t2dm_diamante_sig_cp %>%
  dplyr::select(canonical_pathway, pvalue_adj = p_adjusted_final) %>%
  distinct()
bc_t2dm_diamante_sig_cp$sig = ifelse(bc_t2dm_diamante_sig_cp$pvalue_adj < 0.05, 1, 0)
fwrite(bc_t2dm_diamante_sig_cp, "shared_genes_to_canonical_pathways/BC_T2DM_DIAMANTE_CPs_pvalues_bonferroni.txt", sep = "\t", row.names = FALSE)
