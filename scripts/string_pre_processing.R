
### Pre-processing of STRING protein-protein interaction network (v11.5) ###

library(dplyr)
library(igraph)
library(knitr)

### STRING network was downloaded from here: https://string-db.org/cgi/download?sessionId=bIAz6gR8tk72
### The version we used is v11.5 (the latest at the moment of the study)
## load STRING data
string_net = data.table::fread("https://zenodo.org/records/11540071/files/9606.protein.links.detailed.v11.5.txt?download=1")
string_net_aliasses = data.table::fread("https://zenodo.org/records/11540071/files/9606.protein.aliases.v11.5.txt?download=1") %>%
  filter(source %in% c("Ensembl_HGNC_Entrez_Gene_ID", "Ensembl_HGNC_Entrez_Gene_ID(supplied_by_NCBI)")) %>% 
  dplyr::select(protein_id = `#string_protein_id`, entrezID = alias) %>%
  distinct()

## filter for high confidence score
string_net_high_conf = string_net %>% filter(combined_score >= 700) %>% dplyr::select(protein1, protein2)

## convert string protein IDs to entrezID
string_net_high_conf = left_join(string_net_high_conf, string_net_aliasses, by = c("protein1" = "protein_id"))
string_net_high_conf = left_join(string_net_high_conf, string_net_aliasses, by = c("protein2" = "protein_id"))
string_net_high_conf = string_net_high_conf %>% 
  na.omit() %>%
  dplyr::select(gene1 = entrezID.x, gene2 = entrezID.y) %>%
  distinct()

## create network
g = graph_from_data_frame(string_net_high_conf, directed = FALSE)
g = simplify(g, remove.multiple = TRUE, remove.loops = TRUE) # remove duplicate nodes and loops
string_net_high_conf = get.data.frame(g)
data.table::fwrite(string_net_high_conf, "preprocessed_data/string_high_conf_processed.txt", sep = "\t", row.names = FALSE)
