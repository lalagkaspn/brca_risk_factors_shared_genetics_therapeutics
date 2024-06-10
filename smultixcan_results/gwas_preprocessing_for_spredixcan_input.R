
### Prepare GWAS for predixcan input (harmonization step) ####

library(dplyr)

## -- Breast cancer -- ##
bc = data.table::fread("data/icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt")
# keep predixcan needed columns
bc = bc %>%
  dplyr::select(noneffect_allele = Baseline.Meta,
                effect_allele = Effect.Meta,
                beta = Beta.meta,
                pvalue = p.meta,
                chromosome = chr.iCOGs,
                position = Position.iCOGs) %>%
  distinct()
# save
data.table::fwrite(bc, "BC_unharmonized.txt", sep = "\t", row.names = FALSE)

## -- LDL -- ##
ldl = data.table::fread("data/LDL_GLGC_EUR_2021/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.txt")
# keep predixcan needed columns
ldl = ldl %>%
  dplyr::select(variant_id = rsID,
                noneffect_allele = REF,
                effect_allele = ALT,
                beta = EFFECT_SIZE,
                pvalue = pvalue,
                chromosome = CHROM,
                position = POS_b37,
                effect_allele_frequency = POOLED_ALT_AF) %>%
  distinct()
str(ldl)
ldl$pvalue = as.numeric(ldl$pvalue)
# save
data.table::fwrite(ldl, "LDL_unharmonized.txt", sep = "\t", row.names = FALSE)

## -- HDL -- ##
hdl = data.table::fread("data/HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz")
# keep predixcan needed columns
hdl = hdl %>%
  dplyr::select(variant_id = rsID,
                noneffect_allele = REF,
                effect_allele = ALT,
                beta = EFFECT_SIZE,
                pvalue = pvalue,
                chromosome = CHROM,
                position = POS_b37,
                effect_allele_frequency = POOLED_ALT_AF) %>%
  distinct()

str(hdl)
hdl$pvalue = as.numeric(hdl$pvalue)
# save
data.table::fwrite(hdl, "HDL_unharmonized.txt", sep = "\t", row.names = FALSE)

## -- Depression -- ##
depression = data.table::fread("data/PGC_UKB_depression_genome-wide.txt")
# capitalize A1/A2
depression$A1 = toupper(depression$A1)
depression$A2 = toupper(depression$A2)
# add chromosome and position information using the 1K Reference Genome panel (this is already in b38 --> no need for liftover!)
ref_panel = data.table::fread("https://zenodo.org/records/11540342/files/variant_metadata.txt?download=1")
depression = left_join(depression, ref_panel[, c(1,2,7)], by = c("MarkerName" = "rsid"))
# remove unmatched
depression = na.omit(depression) ; rownames(depression) = NULL
# keep predixcan needed columns
depression = depression %>%
  dplyr::select(variant_id = MarkerName,
                noneffect_allele = A2,
                effect_allele = A1,
                beta = LogOR,
                pvalue = P,
                chromosome = chromosome,
                position = position) %>%
  distinct()
str(depression)
# save
data.table::fwrite(depression, "depression_unharmonized.txt", sep = "\t", row.names = FALSE)

## -- Prostate cancer -- ##

# load data
prostate = data.table::fread("data/PROSTATE_CANCER_PRACTICAL_EUR_2018.txt")

# capitalize A1/A2
prostate$Allele1 = toupper(prostate$Allele1)
prostate$Allele2 = toupper(prostate$Allele2)

# keep predixcan needed columns
prostate = prostate %>%
  dplyr::select(noneffect_allele = Allele2,
                effect_allele = Allele1,
                beta = Effect,
                pvalue = Pvalue,
                chromosome = Chr,
                position = position,
                effect_allele_frequency = Freq1) %>%
  distinct()

str(prostate)

# save
data.table::fwrite(prostate, "prostate_unharmonized.txt", sep = "\t", row.names = FALSE)
rm(prostate)


## -- T2DM DIAMANTE -- ##

# load data
t2dm = data.table::fread("data/DIAMANTE-EUR.sumstat.txt.gz")

# capitalize A1/A2
t2dm$effect_allele = toupper(t2dm$effect_allele)
t2dm$other_allele = toupper(t2dm$other_allele)

# keep predixcan needed columns
t2dm = t2dm %>%
  dplyr::select(variant_id = rsID,
                noneffect_allele = other_allele,
                effect_allele = effect_allele,
                beta = "Fixed-effects_beta",
                pvalue = "Fixed-effects_p-value",
                chromosome = "chromosome(b37)",
                position = "position(b37)",
                effect_allele_frequency = effect_allele_frequency) %>%
  distinct()

str(t2dm)
t2dm$pvalue = as.numeric(t2dm$pvalue)

# save
data.table::fwrite(t2dm, "t2dm_diamante_unharmonized.txt", sep = "\t", row.names = FALSE)
rm(t2dm)

## -- Schizophrenia -- ##

# load data
scz = data.table::fread("data/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz")

# Average A1 frequency in controls and cases
scz$FREQ = rowMeans(scz[, c(6,7)])

# keep predixcan needed columns
scz = scz %>%
  dplyr::select(variant_id = ID,
                noneffect_allele = A2,
                effect_allele = A1,
                beta = BETA,
                pvalue = PVAL,
                chromosome = CHROM,
                position = POS,
                effect_allele_frequency = FREQ) %>%
  distinct()

str(scz)

# save
data.table::fwrite(scz, "scz_unharmonized.txt", sep = "\t", row.names = FALSE)
rm(scz)

