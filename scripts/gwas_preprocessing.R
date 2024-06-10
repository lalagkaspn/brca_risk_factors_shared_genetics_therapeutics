
### Pre-processing of GWAS summary statistics data to use as input to LOGODetect ###

## LOGODetect GitHub repository: https://github.com/ghm17/LOGODetect

## LOGODetect requires the input data to have the following columns:
# - CHR (chromosome)
# - SNP (rsID)
# - BP (hg19)
# - A1 (effect allele)
# - A2 (alternate allele)
# - BETA or OR (effect size)
# - P (p-value)
# NOTE: I get the rsIDs from the 1K reference panel, as provided by the LOGODetect developers to prevent confusion with different versions of dbSNP

library(dplyr)

### 1K Reference Genome panel (from LOGODetect GitHub)
ref_panel = data.table::fread("https://zenodo.org/records/11537759/files/1000G_EUR_QC.bim?download=1") 
ref_panel = ref_panel %>% 
  dplyr::select(CHR = V1, BP = V4, A1 = V5, A2 = V6, SNP = V2) %>%
  mutate(CHR_BP = paste0(CHR, "_", BP)) %>% 
  dplyr::select(CHR_BP, A1_ref = A1, A2_ref = A2, SNP) %>% 
  distinct()

## To run this script, you need GWAS summary statistics data for breast cancer and its predisposing health conditions. Below, you find the links from where to download the GWAS data
## Breast cancer
## - Website: https://bcac.ccge.medschl.cam.ac.uk/bcacdata/oncoarray/oncoarray-and-combined-summary-result/gwas-summary-associations-breast-cancer-risk-2020/
## - File name: icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt
## Depression: 
## - Website: https://pgc.unc.edu/for-researchers/download-results/
## - Look for "mdd2019edinburgh" --> Genome-wide summary statistics from a meta-analysis of PGC and UK Biobank (347.4Mb)
## HDL:
## - Website: https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/
## - File name: HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz
## LDL:
## - Website: https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/
## - File name: LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.txt
## Prostate cancer
## - Website: http://practical.icr.ac.uk/blog/?page_id=8164
## - File name: PROSTATE_CANCER_PRACTICAL_EUR_2018.txt
## Schizophrenia
## - Website: https://pgc.unc.edu/for-researchers/download-results/
## - Look for "scz2022" --> PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz
## Type 2 diabetes:
## - Website: https://t2d.hugeamp.org/dinspector.html?dataset=Mahajan2022_T2D_EU
## Please, download the above GWAS summary statistics data and place them into the "raw_data" directory.

dir.create("raw_data/")
dir.create("preprocessed_data/")

#### Breast Cancer ####
## load data
bc = data.table::fread("data/icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt") %>% 
  dplyr::select(CHR = chr.iCOGs,
                BP = Position.iCOGs,
                A1 = Effect.Meta,
                A2 = Baseline.Meta,
                BETA = Beta.meta,
                P = p.meta) %>% 
  distinct()

## keep bi-allelic SNPs only
bc = bc %>% filter(A1 %in% c("A", "T", "C", "G"),
                   A2 %in% c("A", "T", "C", "G")) %>% distinct()
## add snp rsIDs
bc$CHR_BP = paste0(bc$CHR, "_", bc$BP)
bc = left_join(bc, ref_panel, by = "CHR_BP")

## remove no matches
bc = na.omit(bc)
rownames(bc) = NULL

## check A1/A2 concordance
same = bc$A1 == bc$A1_ref & bc$A2 == bc$A2_ref | bc$A1 == bc$A2_ref & bc$A2 == bc$A1_ref
sum(same)
not_same = !same
sum(not_same)

## keep LOGODetect required columns
bc = bc %>% dplyr::select(CHR, SNP, BP, A1, A2, BETA, P) %>% distinct()

## save
data.table::fwrite(bc, "preprocessed_data/icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics_LOGODetect.txt", sep = "\t", row.names = FALSE)

#### LDL ####
## load data
ldl = data.table::fread("raw_data/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.txt")

## keep bi-allelic SNPs only
ldl = ldl %>% filter(REF %in% c("A", "T", "C", "G"),
                     ALT %in% c("A", "T", "C", "G")) %>% distinct()

## keep needed columns
ldl = ldl %>% 
  dplyr::select(CHR = CHROM, BP = POS_b37, A1 = ALT, A2 = REF, BETA = EFFECT_SIZE, P = pvalue, N) %>% 
  distinct()

## add SNP rsIDs
ldl$CHR_BP = paste0(ldl$CHR, "_", ldl$BP)
ldl = left_join(ldl, ref_panel, by = "CHR_BP")

## remove no matches
ldl = na.omit(ldl)
rownames(ldl) = NULL

## check A1/A2 concordance
same = ldl$A1 == ldl$A1_ref & ldl$A2 == ldl$A2_ref | ldl$A1 == ldl$A2_ref & ldl$A2 == ldl$A1_ref
sum(same)
not_same = !same
sum(not_same)
not_same = which(not_same)
ldl = ldl[-not_same, ]
rownames(ldl) = NULL

## sample size
summary(ldl$N) # mean = 1199124

## keep LOGODetect required columns
ldl = ldl %>% dplyr::select(CHR, SNP, BP, A1, A2, BETA, P) %>% distinct()

## convert p-value to numeric
ldl$P = as.numeric(ldl$P)

## save
data.table::fwrite(ldl, "preprocessed_data/LDL_GLGC_EUR_2021/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results_LOGODetect.txt", sep = "\t", row.names = FALSE)

#### Depression ####
## load data
depression = data.table::fread("raw_data/PGC_UKB_depression_genome-wide.txt")
## sample size: 500,199

## capitalize A1 and A2
depression$A1 = toupper(depression$A1)
depression$A2 = toupper(depression$A2)

## keep bi-allelic SNPs only
depression = depression %>% filter(A1 %in% c("A", "T", "C", "G"),
                                   A2 %in% c("A", "T", "C", "G")) %>% 
  distinct()

## they do not provide CHR or BP --> I add this information using the LOGODetect 1K reference panel
depression = left_join(depression, ref_panel, by = c("MarkerName" = "SNP"))
depression = na.omit(depression)
rownames(depression) = NULL

## keep rows with A1/A2 concordance
same = depression$A1 == depression$A1_ref & depression$A2 == depression$A2_ref | depression$A1 == depression$A2_ref & depression$A2 == depression$A1_ref
sum(same)
not_same = !same
sum(not_same)
not_same = which(not_same)
depression = depression[-not_same, ]
rownames(depression) = NULL

## filter columns
depression = depression %>% dplyr::select(CHR_BP, SNP = MarkerName, A1, A2, BETA = LogOR, P) %>% distinct()

## split CHR_BP column
temp = stringr::str_split_fixed(depression$CHR_BP, "_", n = Inf)
depression = cbind(depression, temp) ; rm(temp)
depression$V1 = as.numeric(depression$V1)
depression$V2 = as.numeric(depression$V2)
depression = depression %>% dplyr::select(CHR = V1, SNP, BP = V2, A1, A2, BETA, P) %>% distinct() %>% arrange(CHR, BP)

## save
data.table::fwrite(depression, "preprocessed_data/PGC_UKB_depression_genome-wide_LOGODetect.txt", sep = "\t", row.names = FALSE)

#### HDL ####
## load data
hdl = data.table::fread("raw_data/HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz")

## keep bi-allelic SNPs only
hdl = hdl %>% filter(REF %in% c("A", "T", "C", "G"),
                     ALT %in% c("A", "T", "C", "G")) %>% 
  distinct()

## filter columns
hdl = hdl %>% dplyr::select(CHR = CHROM, BP = POS_b37, A1 = ALT, A2 = REF, BETA = EFFECT_SIZE, P = pvalue, N) %>% distinct()

## add SNP rsIDs
hdl$CHR_BP = paste0(hdl$CHR, "_", hdl$BP)
hdl = left_join(hdl, ref_panel, by = "CHR_BP")

## remove no matches
hdl = na.omit(hdl)
rownames(hdl) = NULL

## check A1/A2 concordance
same = hdl$A1 == hdl$A1_ref & hdl$A2 == hdl$A2_ref | hdl$A1 == hdl$A2_ref & hdl$A2 == hdl$A1_ref
sum(same)
not_same = !same
sum(not_same)
not_same = which(not_same)
hdl = hdl[-not_same, ]
rownames(hdl) = NULL

## sample size
summary(hdl$N) # mean = 1209154

## keep LOGOGOetect required columns
hdl = hdl %>% dplyr::select(CHR, SNP, BP, A1, A2, BETA, P)

## convert p-value to numeric
hdl$P = as.numeric(hdl$P)

## save
data.table::fwrite(hdl, "preprocessed_data/HDL_GLGC_EUR_2021/HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results_LOGODetect.txt", sep = "\t", row.names = FALSE)

#### Schizophrenia ####
## load data
schizophrenia = data.table::fread("raw_data/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz")

## keep bi-allelic SNPs only
schizophrenia = schizophrenia %>% filter(A1 %in% c("A", "T", "C", "G"),
                                         A2 %in% c("A", "T", "C", "G")) %>% distinct()

## filter for imputation score > 0.3
schizophrenia = schizophrenia %>% filter(IMPINFO >= 0.3)

## total sample size
schizophrenia$N = schizophrenia$NCAS + schizophrenia$NCON

## filter columns
schizophrenia = schizophrenia %>% dplyr::select(CHR = CHROM, BP = POS, A1, A2, BETA, P = PVAL, N) %>% distinct()

## add SNP rsIDs
schizophrenia$CHR_BP = paste0(schizophrenia$CHR, "_", schizophrenia$BP)
schizophrenia = left_join(schizophrenia, ref_panel, by = "CHR_BP")

## remove no matches
schizophrenia = na.omit(schizophrenia)
rownames(schizophrenia) = NULL

## check A1/A2 concordance
same = schizophrenia$A1 == schizophrenia$A1_ref & schizophrenia$A2 == schizophrenia$A2_ref | schizophrenia$A1 == schizophrenia$A2_ref & schizophrenia$A2 == schizophrenia$A1_ref
sum(same)
not_same = !same
sum(not_same)
not_same = which(not_same)
schizophrenia = schizophrenia[-not_same, ]
rownames(schizophrenia) = NULL

## sample size
summary(schizophrenia$N) # mean = 130263

## keep LOGOGOetect required columns
schizophrenia = schizophrenia %>% 
  dplyr::select(CHR, SNP, BP, A1, A2, BETA, P) %>%
  arrange(CHR, BP) %>% distinct()

## save
data.table::fwrite(schizophrenia, "preprocessed_data/PGC3_SCZ_wave3.european.autosome.public.v3.vcf_LOGODetect.txt", sep = "\t", row.names = FALSE)

#### Type 2 Diabetes - DIAMANTE ####
## load data
t2dm_diamante = data.table::fread("raw_data/DIAMANTE-EUR.sumstat.txt.gz")

## keep bi-allelic SNPs only
t2dm_diamante$effect_allele = toupper(t2dm_diamante$effect_allele)
t2dm_diamante$other_allele = toupper(t2dm_diamante$other_allele)
t2dm_diamante = t2dm_diamante %>% 
  filter(effect_allele %in% c("A", "T", "C", "G"),
         other_allele %in% c("A", "T", "C", "G")) %>%
  distinct()

## keep LOGOGOetect required columns
t2dm_diamante = t2dm_diamante %>%
  dplyr::select(CHR = "chromosome(b37)", SNP = rsID, BP = "position(b37)", A1 = effect_allele, A2 = other_allele, BETA = "Fixed-effects_beta", P = "Fixed-effects_p-value") %>%
  distinct() %>%
  arrange(CHR, BP)

## sample size
# 80,154 T2D cases and 853,816
# n = 933970

str(t2dm_diamante)
t2dm_diamante$P = as.numeric(t2dm_diamante$P)

## save
data.table::fwrite(t2dm_diamante, "preprocessed_data/DIAMANTE-EUR.sumstat_LOGODetect.txt", sep = "\t", row.names = FALSE)

#### Prostate cancer ####
## load data
prostate = data.table::fread("raw_data/PROSTATE_CANCER_PRACTICAL_EUR_2018.txt")

## capitalize
prostate$Allele1 = toupper(prostate$Allele1)
prostate$Allele2 = toupper(prostate$Allele2)

## keep bi-allelic SNPs only
prostate = prostate %>% filter(Allele1 %in% c("A", "T", "C", "G"),
                               Allele2 %in% c("A", "T", "C", "G"))

## filter columns
prostate = prostate %>% dplyr::select(CHR = Chr, BP = position, A1 = Allele1, A2 = Allele2, BETA = Effect, P = Pvalue) %>% distinct()

## add SNP rsIDs
prostate$CHR_BP = paste0(prostate$CHR, "_", prostate$BP)
prostate = left_join(prostate, ref_panel, by = "CHR_BP")

## remove no matches
prostate = na.omit(prostate)
rownames(prostate) = NULL

## check A1/A2 concordance
same = prostate$A1 == prostate$A1_ref & prostate$A2 == prostate$A2_ref | prostate$A1 == prostate$A2_ref & prostate$A2 == prostate$A1_ref
sum(same)
not_same = !same
sum(not_same) # 0
not_same = which(not_same)
prostate = prostate[-not_same, ]
rownames(prostate) = NULL

## keep LOGOGOetect required columns
prostate = prostate %>% 
  dplyr::select(CHR, SNP, BP, A1, A2, BETA, P) %>%
  arrange(CHR, BP)

## save
data.table::fwrite(prostate, "preprocessed_data/PROSTATE_CANCER_PRACTICAL_EUR_2018_LOGODetect.txt", sep = "\t", row.names = FALSE)

### Then, use the saved data as input to LOGODetect (default settings) to find positively and negatively correlated loci between breast cancer and each predisposing health condition
## Breast cancer - Depression
## Breast cancer - HDL
## Breast cancer - LDL
## Breast cancer - Prostate cancer
## Breast cancer - Schizophrenia
## Breast cancer - Type 2 diabetes

## The output of LOGODetect for each pair of diseases is saved in "/logodetect_regions" directory (but you can run it by yourself if you install LOGODetect locally)
## LOGODetect was run with default settings
