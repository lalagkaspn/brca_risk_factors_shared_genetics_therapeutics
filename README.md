# Shared genetics between breast cancer and predisposing diseases identifies novel breast cancer treatment candidates

## Panagiotis N. Lalagkas & Rachel D. Melamed

In this work, we propose a novel use of genetics by leveraging pleiotropy to facilitate drug repurposing. As a test case, we focus on breast cancer, the most common type of cancer in women worldwide. Breast cancer has known predisposing health conditions, such as high cholesterol and type 2 diabetes. These predisposing health conditions increase the risk for breast cancer because they share genetic factors with it. Therefore, we hypothesize that by dissecting the shared genetic etiology of breast cancer and a predisposing health condition, we can identify drugs that currently treat the predisposing health condition but also impact shared pathways with cancer. To this end, we analyze publicly available Genome-Wide Association Studies (GWAS) summary statistics data for breast cancer and each predisposing disease to identify genes shared between them. Next, we use a network biology method to link the identified shared genes to canonical pathways. We do the same for all drugs treating the predisposing health condition by linking their targets to these pathways. By finding drugs that target shared pathways, we both prioritize candidate drugs for repurposing for breast cancer and provide biologica linsights that support their effect in disease treatment. Finally, we validdate our list of candidate drugs for repurposing using breast cancer clinical trials data and currently approved indicated drugs.  

There are ten scripts in this repository:

- `gwas_preprocessing.R` is the script to prepare the GWAS summary statistics data for all diseases for local genetic correlation analysis (LOGODetect input).
- `string_pre_processing.R` is the script to preprocess the STRING protein-protein interaction network.
- `preparing_fuma_snp2gene_input.R` is the script to prepare the input files for identifying the genes located within the positivelly correlated loci (FUMA SNP2GENE input).
- `adding_logodetect_qvalues_to_shared_genes.R` is the script to match each gene with the LOGODetect q-value of the locus within which the gene is located.
- `filtering_shared_genes.R` is the script to filter the identifed shared genes using the MAGMA and S-MultiXcan results.
- `linking_shared_genes_to_canonical_pathways.R` is the script to link the filtered shared genes to canonical pathways.
- `linking_drug_targets_to_shared_canonical_pathways.R` is the script to link drugs approved for a predisposing disease to its shared canonical pathways with breast cancer and prioritize candidate drugs for repurposing for breast cancer.
- `evaluation.R` is the script to evaluate the prioritized candidate drugs for repurposing using drugs approved or investigated for breast cancer at the time of the study.
- `figures.R` is the script to reproduce all the figures in the manuscript.
- `supplementary_tables.R` is the script to reproduce all the supplementary tables in the manuscript.

All scripts are written in **R**.

If you have any questions, please reach out to [panagiotis.lalagkas@gmail.com](mailto:panagiotis.lalagkas@gmail.com)
