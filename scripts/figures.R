
### Create figures in the manuscript

library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(gwaRs)
# BiocManager::install("EnsDb.Hsapiens.v75")
# install.packages("locuszoomr")
library(locuszoomr) # Vignette: https://cran.r-project.org/web/packages/locuszoomr/vignettes/locuszoomr.html
# remotes::install_github("mrcieu/ieugwasr")
library(ieugwasr)
# devtools::install_github("explodecomputer/plinkbinr")
library(plinkbinr)
library(AnnotationHub)
library(rtracklayer)
library(igraph)
library(ComplexHeatmap)
library(circlize)

dir.create("figures")

# ------------------------------------------------------------------------------------------------------------ #

### Figure 1

# Figure 1B-1 - mirrored manhattan plot
# https://rdrr.io/github/LindoNkambule/gwaRs/src/R/mirrored_man_plot.R

## data
bc_gwas = data.table::fread("preprocessed_data/icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics_LOGODetect.txt")
bc_gwas = bc_gwas %>% dplyr::select(CHR, SNP, BP, P) %>% mutate(Trait = "Breast cancer")
bc_gwas$P = ifelse(bc_gwas$P == 0, 1e-319, bc_gwas$P)
hdl_gwas = data.table::fread("preprocessed_data/HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results_LOGODetect.txt")
hdl_gwas = hdl_gwas %>% dplyr::select(CHR, SNP, BP, P) %>% mutate(Trait = "HDL")
hdl_gwas$P = ifelse(hdl_gwas$P == 0, 1e-319, hdl_gwas$P)
bc_hdl_logodetect = data.table::fread("logodetect_regions/brca_hdl_logodetect_regions.txt")

## define some variables for mirrored manhattan plot
trait1 = "Breast cancer"
trait2 = "HDL"
trait1_chromCols = c("gray66", "grey36")
trait2_chromCols = c("steelblue1", "steelblue4")

## create input data
input_data = rbind(bc_gwas, hdl_gwas)

## calculate variables to be used for the manhattan plot
input_data = input_data %>%
  dplyr::group_by(CHR) %>%
  dplyr::summarise(chr_len=max(BP), .groups = 'drop') %>%
  dplyr::mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  dplyr::left_join(input_data, ., by=c("CHR"="CHR")) %>%
  dplyr::arrange(CHR, BP) %>%
  dplyr::mutate(BPcum=BP+tot)

# sort the chromosomes in ascending order
input_data = input_data[order(as.numeric(as.character(input_data$CHR))), ]

## filter data
# trait 1
t1 = input_data %>% dplyr::filter(Trait == trait1)
t1 = t1 %>% mutate(log = -log10(P))
t1df = with(t1, data.frame(CHR = unique(t1$CHR), col = rep(trait1_chromCols, length(unique(t1$CHR)))))
t1colors = t1df[!duplicated(t1df$CHR), ] # trait1 chromosome colors
t1input_data = merge(t1,t1colors, by  = "CHR")
# trait 2
t2 = input_data %>% dplyr::filter(Trait == trait2)
t2 = t2 %>% mutate(log = -(-log10(P)))
t2df = with(t2, data.frame(CHR = unique(t2$CHR), col = rep(trait2_chromCols, length(unique(t2$CHR)))))
t2colors = t2df[!duplicated(t2df$CHR), ]
t2input_data = merge(t2,t2colors, by  = "CHR")

## input data for the manhattan plot
maninput_data = rbind(t1input_data, t2input_data)

# total numebr of chromosomes
nCHR = length(unique(maninput_data$CHR))
# define x-axis position (middle of each chromosome range)
axis.set = maninput_data %>% 
  dplyr::group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 , .groups = 'drop')
# define x-axis labels
axis_label = axis.set[order(as.numeric(as.character(axis.set$CHR))), ]
# define y-axis limits
ylim = abs(floor(log10(min(maninput_data$P)))) + 1
# define position of trait annotations on the manhattan plot (top and bottom left corners)
annotations = data.frame(x = c(1,1),
                         y = c(ylim, -ylim),
                         label = c(trait1, trait2))

# common SNPs in genome-wide significant level - to be used for annotation on the manhattan plot
trait1_gwas = maninput_data %>% dplyr::filter(Trait == trait1 & P <= 5e-08)
trait2_gwas = maninput_data %>% dplyr::filter(Trait == trait2 & P <= 5e-08)
common_genomewide_sig = intersect(trait1_gwas$SNP, trait2_gwas$SNP)

# SNPs in LOGODetect positively correlated regions
# GRCh37
trait1_gwas_logodetectsnps = maninput_data %>% 
  dplyr::filter(Trait == "Breast cancer") %>% 
  dplyr::filter(CHR == 6 & BP >= 43756863 & BP <= 43765902 |
                  CHR == 7 & BP >= 72823777 & BP <= 73063515 |
                  CHR == 7 & BP >= 130423972 & BP <= 130467272 |
                  CHR == 9 & BP >= 107625803 & BP <= 107678916 |
                  CHR == 18 & BP >= 57852948 & BP <= 57965193)
trait2_gwas_logodetectsnps = maninput_data %>% 
  dplyr::filter(Trait == "HDL") %>% 
  dplyr::filter(CHR == 6 & BP >= 43756863 & BP <= 43765902 |
                  CHR == 7 & BP >= 72823777 & BP <= 73063515 |
                  CHR == 7 & BP >= 130423972 & BP <= 130467272 |
                  CHR == 9 & BP >= 107625803 & BP <= 107678916 |
                  CHR == 18 & BP >= 57852948 & BP <= 57965193)

maninput_data_bc = maninput_data %>% filter(Trait == "Breast cancer")
ggplot(maninput_data_bc, aes(x=BPcum, y=log)) +
  
  # Show all points
  geom_point(aes(color=as.factor(CHR)), size=0.5) +
  scale_color_manual(values = rep(c("gray66", "grey36"), 22 )) +
  
  # custom axes
  scale_x_continuous(labels = axis_label$CHR, breaks = axis_label$center) +
  scale_y_continuous(breaks = seq(0, 350, 50), labels = abs(seq(0, 350, 50)), limits = c(0, 350)) +
  labs(x = "Chromosome", y = bquote(-log[10](p))) +
  # custom theme
  theme_classic() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(size = 16, family = "Arial", color = "black", angle = 0, vjust = 0.5),
        axis.text.y = element_text(size = 16, family = "Arial", color = "black"),
        axis.title.y = element_text(size = 18, family = "Arial", face = "bold", color = "black", 
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size = 18, family = "Arial", face = "bold", color = "black",
                                    margin = margin(t = 10, r = 0, b = 0, l = 0)))
ggsave("figure_1B_1_bc.png", device = "png", path = "figures/", width = 17, height = 5)

maninput_data_hdl = maninput_data %>% filter(Trait == "HDL")
maninput_data_hdl$log = -1 * maninput_data_hdl$log
ggplot(maninput_data_hdl, aes(x=BPcum, y=log)) +
  
  # Show all points
  geom_point(aes(color=as.factor(CHR)), size=0.5) +
  scale_color_manual(values = rep(c("steelblue1", "steelblue4"), 22 )) +
  
  # custom axes
  scale_x_continuous(labels = axis_label$CHR, breaks = axis_label$center) +
  scale_y_continuous(breaks = seq(0, 350, 50), labels = abs(seq(0, 350, 50)), limits = c(0, 350)) +
  labs(x = "Chromosome", y = bquote(-log[10](p))) +
  # custom theme
  theme_classic() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(size = 16, family = "Arial", color = "black", angle = 0, vjust = 0.5),
        axis.text.y = element_text(size = 16, family = "Arial", color = "black"),
        axis.title.y = element_text(size = 18, family = "Arial", face = "bold", color = "black", 
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size = 18, family = "Arial", face = "bold", color = "black",
                                    margin = margin(t = 10, r = 0, b = 0, l = 0)))
ggsave("figure_1B_1_hdl.png", device = "png", path = "figures/", width = 17, height = 5)

figure_1B_2 = ggplot2::ggplot(maninput_data, aes(x = BPcum, y = log, color = factor(col), label = SNP)) +
  # data points
  geom_point(size = 0.5) +
  # genome-wide line
  geom_hline(yintercept = 5e-08, color = "red", linetype = "dashed") +
  geom_hline(yintercept = -(5e-08), color = "red", linetype = "dashed") +
  # suggestive line
  geom_hline(yintercept = 0.05, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = -(0.05), color = "blue", linetype = "dashed") +
  # trait1-specific annotations
  # geom_point(data = maninput_data %>% dplyr::filter(Trait == trait1 & SNP %in% common_genomewide_sig), color = "red", size = 1) +
  geom_point(data = trait1_gwas_logodetectsnps, color = "red", size = 0.8) +
  
  # ggrepel::geom_label_repel(data = maninput_data %>% dplyr::filter(SNP %in% common_genomewide_sig), label.size = 0.1, size = 3, color = "black") +
  # trait2-specific annotations
  # geom_point(data = maninput_data %>% dplyr::filter(Trait == trait2 & SNP %in% common_genomewide_sig), color = "red", size = 1) +
  geom_point(data = trait2_gwas_logodetectsnps, color = "red", size = 0.8) +
  
  # ggrepel::geom_label_repel(data = maninput_data %>% dplyr::filter(Trait == trait2 & P <= 1e-300), label.size = 0.1, size = 3, color = "black") +
  # axes and titles
  scale_x_continuous(labels = axis_label$CHR, breaks = axis_label$center) +
  scale_y_continuous(breaks = seq(-350, 350, 50), labels = abs(seq(-350, 350, 50)), limits = c(-350, 350)) +
  scale_colour_identity() +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = "Chromosome", y = bquote(-log[10](p))) +
  theme_classic() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(size = 16, family = "Arial", color = "black", angle = 0, vjust = 0.5),
        axis.text.y = element_text(size = 16, family = "Arial", color = "black"),
        axis.title.y = element_text(size = 18, family = "Arial", face = "bold", color = "black", 
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size = 18, family = "Arial", face = "bold", color = "black",
                                    margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ggtitle("") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_hline(yintercept = 0, lwd = 0.5, colour = "black")
# geom_label(data = annotations, aes(x = x, y = y, label = label),
#            color = "black", hjust = 0,
#            fill = c(trait1_chromCols[1], trait2_chromCols[1]),
#            size = 6, fontface = "bold")
figure_1B_2
ggsave("figure_1B_2.png", device = "png", path = "figures/", width = 17, height = 8)

# clean environment
rm(annotations, axis_label, axis.set, input_data, maninput_data, t1, t1colors, t1df, t1input_data, t2, t2colors, t2df, t2input_data, trait1_gwas, trait1_gwas_logodetectsnps, trait2_gwas, trait2_gwas_logodetectsnps, common_genomewide_sig, nCHR, trait1, trait1_chromCols, trait2, trait2_chromCols, ylim,
   bc_gwas, hdl_gwas, bc_hdl_logodetect)

# ### Figure 1B-2 - locuszoom plot
# ## data
# bc_gwas = data.table::fread("data/gwas_summary_statistics/BRCA_BCAC_EUR_2020/icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics_LOGODetect.txt")
# bc_gwas = bc_gwas %>% dplyr::select(CHR, SNP, BP, BETA, P)
# bc_gwas$P = ifelse(bc_gwas$P == 0, 1e-319, bc_gwas$P)
# hdl_gwas = data.table::fread("data/gwas_summary_statistics/HDL_GLGC_EUR_2021/HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results_LOGODetect.txt")
# hdl_gwas = hdl_gwas %>% dplyr::select(CHR, SNP, BP, BETA, P)
# hdl_gwas$P = ifelse(hdl_gwas$P == 0, 1e-319, hdl_gwas$P)
# bc_hdl_logodetect = data.table::fread("LOGODetect/results/BC_HDL/LOGODetect_regions.txt")
# 
# ## filter for LOGODetect regions
# bc_gwas_logodetect = bc_gwas %>% dplyr::filter(CHR == 6 & BP >= 43756863 & BP <= 43765902 |
#                                                  CHR == 7 & BP >= 72823777 & BP <= 73063515 |
#                                                  CHR == 7 & BP >= 130423972 & BP <= 130467272 |
#                                                  CHR == 9 & BP >= 107625803 & BP <= 107678916 |
#                                                  CHR == 18 & BP >= 57852948 & BP <= 57965193)
# hdl_gwas_logodetect = hdl_gwas %>% dplyr::filter(CHR == 6 & BP >= 43756863 & BP <= 43765902 |
#                                                    CHR == 7 & BP >= 72823777 & BP <= 73063515 |
#                                                    CHR == 7 & BP >= 130423972 & BP <= 130467272 |
#                                                    CHR == 9 & BP >= 107625803 & BP <= 107678916 |
#                                                    CHR == 18 & BP >= 57852948 & BP <= 57965193)
# 
# ## Liftover to GRCh38
# # Chromosome column must be in the format "chr#"
# bc_gwas_logodetect$CHR = paste0("chr", bc_gwas_logodetect$CHR)
# hdl_gwas_logodetect$CHR = paste0("chr", hdl_gwas_logodetect$CHR)
# # import chain object for liftover
# # can be found here: 
# # IMPORTANT: You have to import the unzipped file
# # hg17 -> hg38 (hg17ToHg38.over.chain)
# # hg18 -> hg19 (hg18ToHg19.over.chain)
# # hg18 -> hg38 (hg18ToHg38.over.chain)
# # hg19 -> hg38 (hg19ToHg38.over.chain)
# hg19tohg38_chain = import.chain("/project/pi_rachel_melamed_uml_edu/rachel/hold/predixcan/data/liftover/hg19ToHg38.over.chain")
# # Create GRanges from SNP positions
# bc_granges_data = makeGRangesFromDataFrame(bc_gwas_logodetect,
#                                            keep.extra.columns = TRUE,
#                                            seqnames.field = "CHR",
#                                            start.field = "BP",
#                                            end.field = "BP",
#                                            ignore.strand = TRUE)
# hdl_granges_data = makeGRangesFromDataFrame(hdl_gwas_logodetect,
#                                             keep.extra.columns = TRUE,
#                                             seqnames.field = "CHR",
#                                             start.field = "BP",
#                                             end.field = "BP",
#                                             ignore.strand = TRUE)
# # Liftover from GRCh37 to GRCh38 - necessary as locuszoom plot works with GRCh38
# bc_gwas_hg38 = liftOver(x = bc_granges_data, chain = hg19tohg38_chain)
# hdl_gwas_hg38 = liftOver(x = hdl_granges_data, chain = hg19tohg38_chain)
# # Transform GRanges object to data.frame
# bc_gwas_hg38 = as.data.frame(bc_gwas_hg38)
# bc_gwas_hg38 = bc_gwas_hg38 %>% dplyr::select(CHR = seqnames, BP = start, SNP:P)
# hdl_gwas_hg38 = as.data.frame(hdl_gwas_hg38)
# hdl_gwas_hg38 = hdl_gwas_hg38 %>% dplyr::select(CHR = seqnames, BP = start, SNP:P)
# 
# bc_gwas_hg38$CHR = gsub("chr", "", bc_gwas_hg38$CHR)
# hdl_gwas_hg38$CHR = gsub("chr", "", hdl_gwas_hg38$CHR)
# 
# # query ensembl for the most recent database
# ah = AnnotationHub()
# ensmbl_query = query(ah, c("EnsDb", "Homo sapiens")) ; ensmbl_query
# ensDb_v109 = ah[["AH109606"]] # GRCh38 - cannot use more recent due to version of R used
# 
# # estimate LD matrix
# plink = get_plink_exe()
# unique(bc_gwas_hg38$CHR)
# # bc_ld_matrix_chr6 = ld_matrix_local(variants = bc_gwas_hg38$SNP,
# #                                     with_alleles = FALSE,
# #                                     bfile = paste0("/project/pi_rachel_melamed_uml_edu/ref_data/1K_Genomes_phase3_EUR_Reference_Panel/1K_Genomes_phase3_EUR_Reference_panel_MAF_greater_0.01/1kg_eur_1pct_chr", 6),
# #                                     plink_bin = plink)
# bc_ld_matrix_chr7 = ld_matrix_local(variants = bc_gwas_hg38$SNP,
#                                     with_alleles = FALSE,
#                                     bfile = paste0("/project/pi_rachel_melamed_uml_edu/ref_data/1K_Genomes_phase3_EUR_Reference_Panel/1K_Genomes_phase3_EUR_Reference_panel_MAF_greater_0.01/1kg_eur_1pct_chr", 7),
#                                     plink_bin = plink)
# # bc_ld_matrix_chr9 = ld_matrix_local(variants = bc_gwas_hg38$SNP,
# #                                     with_alleles = FALSE,
# #                                     bfile = paste0("/project/pi_rachel_melamed_uml_edu/ref_data/1K_Genomes_phase3_EUR_Reference_Panel/1K_Genomes_phase3_EUR_Reference_panel_MAF_greater_0.01/1kg_eur_1pct_chr", 9),
# #                                     plink_bin = plink)
# # bc_ld_matrix_chr18 = ld_matrix_local(variants = bc_gwas_hg38$SNP,
# #                                      with_alleles = FALSE,
# #                                      bfile = paste0("/project/pi_rachel_melamed_uml_edu/ref_data/1K_Genomes_phase3_EUR_Reference_Panel/1K_Genomes_phase3_EUR_Reference_panel_MAF_greater_0.01/1kg_eur_1pct_chr", 18),
# #                                      plink_bin = plink)
# 
# # hdl_ld_matrix_chr6 = ld_matrix_local(variants = hdl_gwas_hg38$SNP,
# #                                      with_alleles = FALSE,
# #                                      bfile = paste0("/project/pi_rachel_melamed_uml_edu/ref_data/1K_Genomes_phase3_EUR_Reference_Panel/1K_Genomes_phase3_EUR_Reference_panel_MAF_greater_0.01/1kg_eur_1pct_chr", 6),
# #                                      plink_bin = plink)
# hdl_ld_matrix_chr7 = ld_matrix_local(variants = hdl_gwas_hg38$SNP,
#                                      with_alleles = FALSE,
#                                      bfile = paste0("/project/pi_rachel_melamed_uml_edu/ref_data/1K_Genomes_phase3_EUR_Reference_Panel/1K_Genomes_phase3_EUR_Reference_panel_MAF_greater_0.01/1kg_eur_1pct_chr", 7),
#                                      plink_bin = plink)
# # hdl_ld_matrix_chr9 = ld_matrix_local(variants = hdl_gwas_hg38$SNP,
# #                                      with_alleles = FALSE,
# #                                      bfile = paste0("/project/pi_rachel_melamed_uml_edu/ref_data/1K_Genomes_phase3_EUR_Reference_Panel/1K_Genomes_phase3_EUR_Reference_panel_MAF_greater_0.01/1kg_eur_1pct_chr", 9),
# #                                      plink_bin = plink)
# # hdl_ld_matrix_chr18 = ld_matrix_local(variants = hdl_gwas_hg38$SNP,
# #                                       with_alleles = FALSE,
# #                                       bfile = paste0("/project/pi_rachel_melamed_uml_edu/ref_data/1K_Genomes_phase3_EUR_Reference_Panel/1K_Genomes_phase3_EUR_Reference_panel_MAF_greater_0.01/1kg_eur_1pct_chr", 18),
# #                                       plink_bin = plink)
# # remove SNPs from GWAS that have no LD r2
# # bc_snps_w_r2 = c(colnames(bc_ld_matrix_chr6), colnames(bc_ld_matrix_chr7), colnames(bc_ld_matrix_chr9), colnames(bc_ld_matrix_chr18))
# bc_snps_w_r2 = colnames(bc_ld_matrix_chr7)
# bc_gwas_hg38 = bc_gwas_hg38 %>% dplyr::filter(SNP %in% bc_snps_w_r2)
# hdl_snps_w_r2 = colnames(hdl_ld_matrix_chr7)
# hdl_gwas_hg38 = hdl_gwas_hg38 %>% dplyr::filter(SNP %in% hdl_snps_w_r2)
# 
# # add LD r2
# bc_gwas_hg38$r2 = bc_ld_matrix_chr7[, bc_gwas_hg38[which(bc_gwas_hg38$P == min(bc_gwas_hg38$P)), "SNP"]]
# hdl_gwas_hg38$r2 = hdl_ld_matrix_chr7[, hdl_gwas_hg38[which(hdl_gwas_hg38$P == min(hdl_gwas_hg38$P)), "SNP"]]
# 
# bc_gwas_hg38$CHR = gsub("chr", "", bc_gwas_hg38$CHR)
# hdl_gwas_hg38$CHR = gsub("chr", "", hdl_gwas_hg38$CHR)
# # create locus zoom plot
# bc_gwas_hg38$P_neglog10 = -log10(bc_gwas_hg38$P)
# loc_bc = locus(data = bc_gwas_hg38,
#                # gene = 'BCL7B',
#                # fix_window = 1e6,
#                seqname = 7, xrange = c(73409447,73649185), # must be GRCh38! # this is from the 1st LOGODetect region in Chr7
#                ens_db = ensDb_v109,
#                chrom = "CHR", pos = "BP", labs = "SNP", p = "P",
#                # yvar = "BETA", # change the variable on the y-axis - if this is specified, then you cannot specify the P argument (column with pvalues)
#                LD = "r2") # this is the column that contains the LD r2 in reference to the selected SNP (column from LD matrix)
# summary(loc_bc)
# # customize the plot
# locus_plot(loc_bc,
#            # labels = c("index"), # mention specific SNPs - index is the SNP with the lowest p-value in the locus
#            # label_x = c(3, -2), # position of the SNP name indexed above
#            # filter_gene_name = c("FZD9", "BAZ1B"), # show only specific genes
#            filter_gene_biotype = "protein_coding", # show only protein coding genes
#            # xlab = "Chromosome",
#            # ylab = "Effect size\n",
#            cex = 1, # size of points
#            cex.axis = 1, # axis.text font size    
#            cex.lab = 1.2, # axis.title font size
#            cex.text = 1,  # gene names font size
#            heights = 2, # ratio top-bottom (top=data points; bottom = gene_names)
#            showExons = TRUE, # show exons or only gene (just a box without exon-intron information)
#            xticks = "bottom", # chromosome position should be right below the data points ("top") or below the gene names ("bottom")
#            border = FALSE, # bounding box plotted to separate top and bottom parts of the figure 
#            gene_col = "blue", # color for gene lines
#            exon_col = "blue",  # fill color for exons
#            exon_border = TRUE, # border for exon or not
#            text_pos = "top", # name of genes to be "top" or "left" from their boxes
#            # recomb_col = "recomb" # column of recombination rate
#            legend_pos = NA,
#            align = TRUE, # to align the plots if not alignedf
#            yzero = TRUE # force y-axis to include 0,
# )
# # manually save - figure_1B_2a.png
# 
# loc_hdl = locus(data = hdl_gwas_hg38,
#                 # gene = 'hdlL7B',
#                 # fix_window = 1e6,
#                 seqname = 7, xrange = c(73409447,73649185), # must be GRCh38! # this is from the 1st LOGODetect region in Chr7
#                 ens_db = ensDb_v109,
#                 chrom = "CHR", pos = "BP", labs = "SNP", p = "P",
#                 # yvar = "BETA", # change the variable on the y-axis - if this is specified, then you cannot specify the P argument (column with pvalues)
#                 LD = "r2") # this is the column that contains the LD r2 in reference to the selected SNP (column from LD matrix)
# summary(loc_hdl)
# # customize the plot
# locus_plot(loc_hdl,
#            # labels = c("index"), # mention specific SNPs - index is the SNP with the lowest p-value in the locus
#            # label_x = c(3, -2), # position of the SNP name indexed above
#            # filter_gene_name = c("FZD9", "BAZ1B"), # show only specific genes
#            filter_gene_biotype = "protein_coding", # show only protein coding genes
#            # xlab = "Chromosome",
#            # ylab = "Effect size\n",
#            cex = 1, # size of points
#            cex.axis = 1, # axis.text font size    
#            cex.lab = 1.2, # axis.title font size
#            cex.text = 1,  # gene names font size
#            heights = 2, # ratio top-bottom (top=data points; bottom = gene_names)
#            showExons = TRUE, # show exons or only gene (just a box without exon-intron information)
#            xticks = "bottom", # chromosome position should be right below the data points ("top") or below the gene names ("bottom")
#            border = FALSE, # bounding box plotted to separate top and bottom parts of the figure 
#            gene_col = "blue", # color for gene lines
#            exon_col = "blue",  # fill color for exons
#            exon_border = TRUE, # border for exon or not
#            text_pos = "top", # name of genes to be "top" or "left" from their boxes
#            # recomb_col = "recomb" # column of recombination rate
#            legend_pos = "topleft",
#            align = TRUE, # to align the plots if not alignedf
#            yzero = TRUE # force y-axis to include 0,
# )
# # manually save - figure_1B_2b.png
# 
# rm(ah, bc_granges_data, bc_gwas, bc_gwas_hg38, bc_gwas_logodetect, bc_hdl_logodetect, bc_ld_matrix_chr18, bc_ld_matrix_chr6, bc_ld_matrix_chr7, bc_ld_matrix_chr9, ensDb_v109, ensmbl_query, hdl_granges_data, hdl_gwas, hdl_gwas_hg38, hdl_gwas_logodetect, hdl_ld_matrix_chr18, hdl_ld_matrix_chr6, hdl_ld_matrix_chr7, hdl_ld_matrix_chr9, hg19tohg38_chain, x,
#    bc_snps_w_r2, hdl_snps_w_r2, plink)

### Figure 1B-3 - protein-coding genes per positively correlated LOGODetect loci
# BC-HDL positively correlated loci shared loci
bc_hdl_logodetect = data.table::fread("logodetect_regions/brca_hdl_logodetect_regions.txt")
# BRCA-HDL shared genes in STRING high confidence network
string_high_conf = fread("preprocessed_data/string_high_conf_processed.txt")
genes_in_string_high_conf = unique(c(string_high_conf$from, string_high_conf$to))
bc_hdl_genes_in_pos_reg = fread("fuma_snp2gene/mapped_genes/brca_hdl_fuma_snp2gene_genes.txt") %>% filter(entrezID %in% genes_in_string_high_conf) %>% dplyr::select(entrezID, chr, start, end) %>% distinct()
fig_1b_3_data = data.frame(locus = c("chr18:57852948-57965193","chr9:107625803-107678916","chr7:130423972-130467272","chr7:72823777-73063515","chr6:43756863-43765902"),
                           mapped_genes = c(0,1,0,5,1))
fig_1b_3_data$locus = factor(fig_1b_3_data$locus, levels = fig_1b_3_data$locus, labels = fig_1b_3_data$locus)
figure_1B_3 = ggplot(fig_1b_3_data, aes(x = mapped_genes, y = locus)) +
  geom_col(fill = "darkorange3", width = 0.6) +
  xlab("Number of mapped genes") +
  ylab("Genomic locus") +
  theme_classic() +
  theme(axis.text = element_text(size = 16, family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", color = "black", face = "bold"))
figure_1B_3
ggsave("figures/figure_1B_3.png", device = "png", width = 7, height = 3, dpi = 600)

rm(bc_hdl_logodetect, string_high_conf, genes_in_string_high_conf, bc_hdl_genes_in_pos_reg, fig_1b_3_data)

### Figure 1B-4 - applying MAGMA and S-MultiXcan filters
## load BRCA-HDL LOGODetect genes that are in the STRING network
string_high_conf = fread("preprocessed_data/string_high_conf_processed.txt")
genes_in_string_high_conf = unique(c(string_high_conf$from, string_high_conf$to))
bc_hdl = fread("fuma_snp2gene/mapped_genes/brca_hdl_fuma_snp2gene_genes_logodetect_qvalues.txt") %>%
  dplyr::filter(entrezID %in% genes_in_string_high_conf)
## load S-MultiXcan results
bc_smultixcan = fread("smultixcan_results/BRCA_BCAC_EUR_2020_imputed_smultixcan.txt")
hdl_smultixcan = fread("smultixcan_results/HDL_GLGC_EUR_2021_imputed_smultixcan.txt")

## -- load MAGMA results -- ##
hdl_magma = fread("magma_results/hdl.genes.out")

## -- filter for MAGMA significant genes for the BC-related disease and same S-MultiXcan direction -- ##
## BC-HDL
# filter for hdl-MAGMA significant genes (after adjusting for the number of original LOGODetect genes)
hdl_magma_shared = hdl_magma %>% 
  dplyr::filter(GENE %in% bc_hdl$entrezID) %>%
  mutate(P.adj = p.adjust(P, method = "BH", n = length(P))) %>%
  mutate(hdl_magma_sig = if_else(P.adj < 0.05, 1, 0)) %>%
  # filter(P.adj < 0.05) %>%
  dplyr::select(GENE, P.adj, hdl_magma_sig) %>%
  distinct()
# filter bc-smultixcan genes for the LOGODetect ones and keep the direction of z
bc_smultixcan_shared = bc_smultixcan %>%
  dplyr::filter(gene_name %in% bc_hdl$gene_name) %>%
  dplyr::select(gene_name, z_mean) %>%
  distinct() %>%
  # mutate(z_sign = if_else(z_mean > 0, 1, 0)) %>%
  # dplyr::select(-z_mean) %>%
  distinct()
# filter hdl-smultixcan genes for the LOGODetect ones and keep the direction of z
hdl_smultixcan_shared = hdl_smultixcan %>% 
  dplyr::filter(gene_name %in% bc_hdl$gene_name) %>%
  dplyr::select(gene_name, z_mean) %>%
  distinct() %>%
  # mutate(z_sign = if_else(z_mean > 0, 1, 0)) %>%
  # dplyr::select(-z_mean) %>%
  distinct()
# combine bc and hdl smultixcan genes and keep the ones that have the same direction (either both positive or both negative)
bc_hdl_smultixcan_shared_combined = left_join(bc_smultixcan_shared, hdl_smultixcan_shared, by = "gene_name")
bc_hdl_magma_smultixcan_dir_shared = bc_hdl %>%
  left_join(hdl_magma_shared, by = c("entrezID" = "GENE")) %>%
  left_join(bc_hdl_smultixcan_shared_combined, by = "gene_name") %>%
  dplyr::rename(z_mean_hdl = z_mean.y, z_mean_bc = z_mean.x) %>%
  na.omit() %>%
  distinct() 
figure_1B_4 = ggplot(bc_hdl_magma_smultixcan_dir_shared, aes(x = z_mean_hdl, y = z_mean_bc, color = factor(hdl_magma_sig, levels = c(0, 1), labels = c("≥0.05", "<0.05")))) +
  # green rectangle for the selected genes after applying MAGMA + S-MultiXcan filters
  geom_rect(aes(xmin = 0, xmax = 14.6, ymin = 0, ymax = 3), color = "green", fill = "green", alpha = 0.02) +
  geom_rect(aes(xmin = -3, xmax = 0, ymin = -1, ymax = 0), color = "green", fill = "green", alpha = 0.02) +
  geom_point(size = 6) +
  labs(color = "Predisposing disease MAGMA P.adjusted") +
  xlab("Predisposing disease\n S-MultiXcan z-mean") +
  ylab("Breast cancer\n S-MultiXcan z-mean") +
  scale_x_continuous(breaks = seq(-3, 14, 1), limits = c(-3, 14.6)) +
  scale_y_continuous(breaks = seq(-1, 3, 1), limits = c(-1, 3)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_text_repel(data=bc_hdl_magma_smultixcan_dir_shared,
                  aes(label=gene_name),
                  fontface = "italic",
                  family = "Arial",
                  color = "black",
                  direction = "y",
                  size = 8, nudge_x = 0.5, nudge_y = -0.1) +
  geom_label(label = "S-MultiXcan same direction", x=3.8, y=2.7, size = 8,
             label.padding = unit(0.4, "lines"), # Rectangle size around label
             color = "black",
             fill = "green", alpha = 0.05
  ) +
  theme_classic() +
  theme(axis.text = element_text(size = 20, color = "black", family = "Arial"),
        axis.title = element_text(size = 24, color = "black", family = "Arial"),
        legend.position = "top",
        legend.text = element_text(size = 24, color = "black", family = "Arial"),
        legend.title = element_text(size = 24, color = "black", family = "Arial"),
        legend.key.size = unit(10, 'mm'))
figure_1B_4
ggsave("/project/pi_rachel_melamed_uml_edu/Panos/GWAS/figures/figure_1B_4.png", device = "png", width = 12, height = 7)

rm(bc_hdl, bc_hdl_magma_smultixcan_dir_shared, bc_hdl_smultixcan_shared_combined, bc_smultixcan, bc_smultixcan_shared, hdl_magma, hdl_magma_shared, hdl_smultixcan, hdl_smultixcan_shared, string_high_conf, genes_in_string_high_conf)

### Figure 1B-5 - PPI network, shared genes, canonical pathways and drug targets
## In this part of the figure, I will plot a small network that shows how we connect shared genes to drug targets through canonical pathways
## Chosen example: BRCA-HDL MILXP shared gene; REACTOME_AMPK_INHIBITS_CHREBP_TRANSCRIPTIONAL_ACTIVATION_ACTIVITY pathway; HMGCR drug target

## HDL drugs and connected CPs
hdl_drugs_cps = fread("drugs_to_shared_canonical_pathways/BC_HDL_drug_targets_perm_pvalues_bonferroni.txt")
## drugs and to which pathways they are connected
hdl_drugs_cps %>% dplyr::filter(perm_pvalue_adj < 0.05)
## BRCA-HDL filtered shared genes
bc_hdl_logodetect_genes_filtered = fread("magma_smultixcan_filtered_shared_genes/bc_hdl_filtered_shared_genes.txt")
## Genes in canonical pathway
msigdb_pathways = readRDS("preprocessed_data/c2.cp.v2023.2.Hs.entrez_STRING_filtered.RDS")
msigdb_pathways[["REACTOME_AMPK_INHIBITS_CHREBP_TRANSCRIPTIONAL_ACTIVATION_ACTIVITY"]]
## BRCA-HDL connected canonical pathways to filtered shared genes
bc_hdl_cps = fread("shared_genes_to_canonical_pathways/BC_HDL_CPs_pvalues_bonferroni.txt")
## STRING network
string_net = fread("preprocessed_data/string_high_conf_processed.txt")

# plot network that highlights shared gene and genes in the pathway
figure_1B_5_net = string_net %>% dplyr::filter(from %in% c(51085, 3156, msigdb_pathways[["REACTOME_AMPK_INHIBITS_CHREBP_TRANSCRIPTIONAL_ACTIVATION_ACTIVITY"]], 10608, 60491, 10628),
                                               to %in% c(51085, 3156, msigdb_pathways[["REACTOME_AMPK_INHIBITS_CHREBP_TRANSCRIPTIONAL_ACTIVATION_ACTIVITY"]], 10608, 60491, 10628))
figure_1B_5_net_igraph = graph_from_data_frame(figure_1B_5_net, directed = FALSE)
V(figure_1B_5_net_igraph)$label = c("PRKAB2","PRKAG2","HMGCR","MLXIPL","STK11","MXD4","ADIPOR1","ADIPOR2","PRKAA2","ADIPOQ","NIF3L1","TXNIP")
png("figures/figure_1B_5.png", res = 300, width = 1800, height = 1800)
set.seed(1122334455)
plot(figure_1B_5_net_igraph, layout = layout_with_fr,
     # customize vertices
     vertex.shape="circle", # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
     vertex.frame.color = "black",
     vertex.color = c("steelblue","steelblue","green4","orange","steelblue", "gray","steelblue","steelblue","steelblue","steelblue","gray","gray"),
     vertex.size = 20,
     # customize vertex labels
     # vertex.label=LETTERS[1:10], # Character vector used to label the nodes
     vertex.label.color = "black",
     vertex.label.family = "Arial",                  # Font family of the label (e.g.“Times”, “Helvetica”)
     vertex.label.font = 4,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
     vertex.label.dist = c(2.6,-2.6,2.6,2.6,-2.6,3,2.6,2.6,-2.6,2.6,-2.6,-2.6), 
     vertex.label.degree = c(300,250,300,300,250,180,300,300,250,250,250,250) , # The position of the label in relation to the vertex (use pi)
     vertex.label.cex = 0.9, # Font size (multiplication factor, device-dependent)
     # customize edges
     edge.color = c("gray60","darkred","darkred","gray60","gray60","gray60","gray60","gray60","darkred","darkred","gray60","gray60","gray60","darkred",
                    "gray60","gray60","gray60","darkred","gray60","gray60","gray60","gray60","gray60","gray60","gray60","gray60","gray60","gray60"), # Edge color
     edge.width = 1.8, # Edge width, defaults to 1
     # edge.arrow.size = 0.2, # Arrow size, defaults to 1
     # edge.arrow.width = 0.2, # Arrow width, defaults to 1
     edge.lty = c("solid","dashed","dashed","solid","solid","solid","solid","solid","dashed","dashed","solid","solid","solid","dashed",
                  "solid","solid","solid","dashed","solid","solid","solid","solid","solid","solid","solid","solid","solid","solid"), # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
     edge.curved = 0) # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
# customize title
# main = expression(bold("Connecting shared gene"~bolditalic(MILXP)~"to pathway")))
dev.off()

rm(bc_hdl_cps, bc_hdl_logodetect_genes_filtered, hdl_drugs_cps, msigdb_pathways, string_net)

# ------------------------------------------------------------------------------------------------------------ #

### Figure 2 ###

# 2A: Shared loci between BRCA and clinically associated diseases - quantitatively
bc_depression_shared_loci = fread("logodetect_regions/brca_depression_logodetect_regions.txt") %>% mutate(direction = if_else(stat > 0, "positive", "negative"), disease = "Depression")
bc_hdl_shared_loci = fread("logodetect_regions/brca_hdl_logodetect_regions.txt") %>% mutate(direction = if_else(stat > 0, "positive", "negative"), disease = "high HDL")
bc_ldl_shared_loci = fread("logodetect_regions/brca_ldl_logodetect_regions.txt") %>% mutate(direction = if_else(stat > 0, "positive", "negative"), disease = "high LDL")
bc_prostate_shared_loci = fread("logodetect_regions/brca_prostate_logodetect_regions.txt") %>% mutate(direction = if_else(stat > 0, "positive", "negative"), disease = "Prostate cancer")
bc_scz_shared_loci = fread("logodetect_regions/brca_schizophrenia_regions.txt") %>% mutate(direction = if_else(stat > 0, "positive", "negative"), disease = "Schizophrenia")
bc_t2dm_shared_loci = fread("logodetect_regions/brca_t2dm_logodetect_regions.txt") %>% mutate(direction = if_else(stat > 0, "positive", "negative"), disease = "Type 2 diabetes")

shared_loci_ultimate = rbind(bc_depression_shared_loci, bc_hdl_shared_loci, bc_ldl_shared_loci, bc_prostate_shared_loci, bc_scz_shared_loci, bc_t2dm_shared_loci)
shared_loci_ultimate = shared_loci_ultimate %>%
  group_by(disease) %>%
  mutate(Negative = sum(direction == "negative"),
         Positive = sum(direction == "positive")) %>%
  ungroup() %>%
  dplyr::select(disease, Negative, Positive) %>% 
  distinct()
shared_loci_ultimate = reshape2::melt(shared_loci_ultimate, "disease", colnames(shared_loci_ultimate)[2:ncol(shared_loci_ultimate)])

shared_loci_ultimate_positive = shared_loci_ultimate %>% dplyr::filter(variable == "Positive") %>% arrange(value)
shared_loci_ultimate_negative = shared_loci_ultimate %>% dplyr::filter(variable == "Negative") %>% arrange(value)

shared_loci_ultimate_positive$disease = factor(shared_loci_ultimate_positive$disease, levels = c("high HDL", "Schizophrenia", "Depression", "high LDL", "Prostate cancer", "Type 2 diabetes"), labels = c("high HDL", "Schizophrenia", "Depression", "high LDL", "Prostate cancer", "Type 2 diabetes"))
shared_loci_ultimate_negative$disease = factor(shared_loci_ultimate_negative$disease, levels = c("high HDL", "Schizophrenia", "Depression", "high LDL", "Prostate cancer", "Type 2 diabetes"), labels = c("high HDL", "Schizophrenia", "Depression", "high LDL", "Prostate cancer", "Type 2 diabetes"))

p1 = ggplot(shared_loci_ultimate_negative, aes(x = value, y = disease)) +
  geom_col(fill = "#F8766D") + 
  labs(title = bquote(atop(bold("Negatively"), "correlated loci with breast cancer"))) +
  xlab("Number of loci") +
  ylab("") +
  scale_x_reverse(breaks = seq(0, 40, 5), limits = c(40, 0)) +
  scale_y_discrete(position = "right") +
  theme_bw() +
  theme(plot.title = element_text(size = 60, color = "black", family = "Arial", margin = margin(t = 40, b = 40)),
        panel.grid.major.x = element_line(color = "lightgrey"),
        panel.grid.major.y = element_blank(),
        axis.text.y.right = element_text(size = 55, color = "black", family = "Arial", face = "bold", hjust = 0.5,
                                         margin = margin(t = 0, r = 0, b = 0, l = 40)),
        axis.text.x = element_text(size = 55, color = "black", family = "Arial"),
        axis.title.x = element_text(size = 58, color = "black", family = "Arial", margin = margin(t = 40)))
p2 = ggplot(shared_loci_ultimate_positive, aes(x = value, y = disease)) +
  geom_col(fill = "#1ABC9C") + 
  labs(title = bquote(atop(bold("Positively"), "correlated loci with breast cancer"))) +
  xlab("Number of loci") +
  ylab("") +
  scale_x_continuous(breaks = seq(0, 40, 5), limits = c(0, 40)) +
  theme_bw() +
  theme(plot.title = element_text(size = 60, color = "black", family = "Arial", margin = margin(t = 40, b = 40)),
        panel.grid.major.x = element_line(color = "lightgrey"),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 55, color = "black", family = "Arial"),
        axis.title.x = element_text(size = 58, color = "black", family = "Arial", margin = margin(t = 40)))
figure_2A = ggarrange(p1, p2, nrow = 1, align = "hv")
figure_2A
ggsave("figures/figure_2A.png", device = "png", width = 36, height = 12)

# 2B: Shared loci between BRCA and clinically associated diseases - figure style like manhattan plot
shared_loci_ultimate = rbind(bc_depression_shared_loci, bc_hdl_shared_loci, bc_ldl_shared_loci, bc_prostate_shared_loci, bc_scz_shared_loci, bc_t2dm_shared_loci)
shared_loci_ultimate = shared_loci_ultimate %>% 
  arrange(chr, begin_pos) %>%
  dplyr::select(chr, begin_pos, stop_pos, stat, qval, direction, disease) %>% 
  distinct()
# length of each chromosome
chrom_length = data.frame(chr = seq(1,22,1),
                          start = c(rep(1, 22)),
                          end = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 
                                  141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753,
                                  81195210, 78077248, 59128983, 63025520, 48129895, 51304566))
chrom_length = reshape2::melt(chrom_length, "chr", colnames(chrom_length)[2:ncol(chrom_length)])
chrom_length = chrom_length %>% 
  dplyr::select(chr, pos = variable, bp = value) %>% 
  arrange(chr, bp)
chr_length_for_plot = chrom_length %>%
  dplyr::group_by(chr) %>%
  dplyr::summarize(chr_len=max(bp), .group = "drop") %>%
  dplyr::mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  dplyr::left_join(chrom_length, ., by=c("chr"="chr")) %>%
  dplyr::arrange(chr, bp) %>%
  dplyr::mutate(BPcum = bp + tot)
# define x-axis position (middle of each chromosome range)
axis.set = chr_length_for_plot %>% 
  dplyr::group_by(chr) %>% 
  dplyr::summarize(center=(max(BPcum) + min(BPcum) ) / 2 , .groups = 'drop')
# define x-axis labels
axis_label = axis.set[order(as.numeric(as.character(axis.set$chr))), ]
# find coordinates of shared loci in the figure
shared_loci_ultimate = left_join(shared_loci_ultimate, 
                                 chr_length_for_plot %>% dplyr::filter(pos == "start") %>% dplyr::select(chr, BPcum), by = "chr")
shared_loci_ultimate$begin_pos_figure = shared_loci_ultimate$begin_pos + shared_loci_ultimate$BPcum - 1
shared_loci_ultimate$stop_pos_figure = shared_loci_ultimate$stop_pos + shared_loci_ultimate$BPcum - 1

shared_loci_ultimate$neg_log10_qval = -log10(shared_loci_ultimate$qval)
shared_loci_ultimate$neg_log10_qval = ifelse(shared_loci_ultimate$direction == "positive", shared_loci_ultimate$neg_log10_qval, shared_loci_ultimate$neg_log10_qval * - 1)
shared_loci_ultimate$disease = factor(shared_loci_ultimate$disease, levels = c("Type 2 diabetes", "Prostate cancer", "high LDL", "Depression", "Schizophrenia", "high HDL"), labels = c("Type 2 diabetes", "Prostate cancer", "high LDL", "Depression", "Schizophrenia", "high HDL"))
figure_2B = ggplot(chr_length_for_plot, aes(x = BPcum)) +
  # geom_rect(xmin = 1, xmax = 249250621, ymin = -Inf,ymax = Inf, fill = "gray95") +
  geom_rect(xmin = 249250622, xmax = 492449994, ymin = -Inf,ymax = Inf, fill = "gray95") +
  # geom_rect(xmin = 492449995, xmax = 690472424, ymin = -Inf,ymax = Inf,  fill = "gray95") +
  geom_rect(xmin = 690472425, xmax = 881626700, ymin = -Inf,ymax = Inf, fill = "gray95") +
  # geom_rect(xmin = 881626701, xmax = 1062541960, ymin = -Inf,ymax = Inf, fill = "gray95") +
  geom_rect(xmin = 1062541961, xmax = 1233657027, ymin = -Inf,ymax = Inf, fill = "gray95") +
  # geom_rect(xmin = 1233657028, xmax = 1392795690, ymin = -Inf,ymax = Inf, fill = "gray95") +
  geom_rect(xmin = 1392795691, xmax = 1539159712, ymin = -Inf,ymax = Inf, fill = "gray95") +
  # geom_rect(xmin = 1539159713, xmax = 1680373143, ymin = -Inf,ymax = Inf, fill = "gray95") +
  geom_rect(xmin = 1680373144, xmax = 1815907890, ymin = -Inf,ymax = Inf, fill = "gray95") +
  # geom_rect(xmin = 1815907891, xmax = 1950914406, ymin = -Inf,ymax = Inf, fill = "gray95") +
  geom_rect(xmin = 1950914407, xmax = 2084766301, ymin = -Inf,ymax = Inf, fill = "gray95") +
  # geom_rect(xmin = 2084766302, xmax = 2199936179, ymin = -Inf,ymax = Inf, fill = "gray95") +
  geom_rect(xmin = 2199936180, xmax = 2307285719, ymin = -Inf,ymax = Inf, fill = "gray95") +
  # geom_rect(xmin = 2307285720, xmax = 2409817111, ymin = -Inf,ymax = Inf, fill = "gray95") +
  geom_rect(xmin = 2409817112, xmax = 2500171864, ymin = -Inf,ymax = Inf, fill = "gray95") +
  # geom_rect(xmin = 2500171865, xmax = 2581367074, ymin = -Inf,ymax = Inf, fill = "gray95") +
  geom_rect(xmin = 2581367075, xmax = 2659444322, ymin = -Inf,ymax = Inf, fill = "gray95") +
  # geom_rect(xmin = 2659444323, xmax = 2718573305, ymin = -Inf,ymax = Inf, fill = "gray95") +
  geom_rect(xmin = 2718573306, xmax = 2781598825, ymin = -Inf,ymax = Inf, fill = "gray95") +
  # geom_rect(xmin = 2781598826, xmax = 2829728720, ymin = -Inf,ymax = Inf, fill = "gray95") +
  geom_rect(xmin = 2829728721, xmax = 2881033286, ymin = -Inf,ymax = Inf, fill = "gray95") +
  geom_point(data = shared_loci_ultimate, 
             aes(x = begin_pos_figure, y = neg_log10_qval, color = disease, shape = disease),
             size = 12, alpha = 1) +
  scale_shape_manual(values = c(16,17,15,18,9,8)) + 
  scale_x_continuous(labels = axis_label$chr, breaks = axis_label$center, expand = c(0, 1), guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(breaks = seq(-3, 3, 1), labels = abs(seq(-3, 3, 1))) +
  scale_color_manual(values = c("#1f78b4", "#e31a1c", "#33a02c", "#6a3d9a", "#ff7f00", "deeppink4")) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  labs(x = "Chromosome", y = bquote(-log[10]("LOGODetect q-value"))) +
  geom_label(aes(x = 0, hjust = -0.02, y = 3, label = "Positively correlated loci "), fill = "#1ABC9C", color = "white", fontface = "bold", family = "Arial", size = 24) +
  geom_label(aes(x = 0, hjust = -0.02, y = -3, label = "Negatively correlated loci "), fill = "#F8766D", color = "white", fontface = "bold", family = "Arial", size = 24) +
  theme_classic() +
  theme(panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(size = 60, family = "Arial", color = "black", angle = 0, vjust = 0.5),
        axis.text.y = element_text(size = 60, family = "Arial", color = "black"),
        axis.title.y = element_text(size = 70, family = "Arial", color = "black", 
                                    margin = margin(t = 0, r = 25, b = 0, l = 0)),
        axis.title.x = element_text(size = 70, family = "Arial", color = "black",
                                    margin = margin(t = 25, r = 0, b = 0, l = 0)),
        legend.text = element_text(size = 60, family = "Arial", color = "black", margin = margin(t = 15, b = 15, r = 0, l = 0)),
        legend.title = element_blank(), legend.key.size = unit(2, units = "cm")) +
  guides(shape = guide_legend(override.aes = list(size = 25)))
figure_2B
ggsave("figures/figure_2B.png", device = "png", width = 54, height = 20, limitsize = FALSE)

## 2C: number of genes in positively correlated regions and number of genes after applying MAGMA/S-MultiXcan filters
# filtered for genes in STRING network
string_high_conf = fread("preprocessed_data/string_high_conf_processed.txt")
genes_in_string_high_conf = unique(c(string_high_conf$from, string_high_conf$to))
# shared genes in positively correlated regions (from FUMA SNP2GENE)
bc_depression_genes_in_pos_reg = fread("fuma_snp2gene/mapped_genes/brca_depression_fuma_snp2gene_genes_logodetect_qvalues.txt") %>% dplyr::filter(entrezID %in% genes_in_string_high_conf)
bc_hdl_genes_in_pos_reg = fread("fuma_snp2gene/mapped_genes/brca_hdl_fuma_snp2gene_genes_logodetect_qvalues.txt") %>% dplyr::filter(entrezID %in% genes_in_string_high_conf)
bc_ldl_genes_in_pos_reg = fread("fuma_snp2gene/mapped_genes/brca_ldl_fuma_snp2gene_genes_logodetect_qvalues.txt") %>% dplyr::filter(entrezID %in% genes_in_string_high_conf)
bc_prostate_genes_in_pos_reg = fread("fuma_snp2gene/mapped_genes/brca_prostate_fuma_snp2gene_genes_logodetect_qvalues.txt") %>% dplyr::filter(entrezID %in% genes_in_string_high_conf)
bc_scz_genes_in_pos_reg = fread("fuma_snp2gene/mapped_genes/brca_schizophrenia_fuma_snp2gene_genes_logodetect_qvalues.txt") %>% dplyr::filter(entrezID %in% genes_in_string_high_conf)
bc_t2dm_genes_in_pos_reg = fread("fuma_snp2gene/mapped_genes/brca_t2dm_fuma_snp2gene_genes_logodetect_qvalues.txt") %>% dplyr::filter(entrezID %in% genes_in_string_high_conf)
# shared genes filtered for MAGMA sig BC-related disease and same S-MultiXcan direction
bc_depression_genes_in_pos_reg_filtered = fread("magma_smultixcan_filtered_shared_genes/bc_depression_filtered_shared_genes.txt")
bc_hdl_genes_in_pos_reg_filtered = fread("magma_smultixcan_filtered_shared_genes/bc_hdl_filtered_shared_genes.txt")
bc_ldl_genes_in_pos_reg_filtered = fread("magma_smultixcan_filtered_shared_genes/bc_ldl_filtered_shared_genes.txt")
bc_prostate_genes_in_pos_reg_filtered = fread("magma_smultixcan_filtered_shared_genes/bc_prostate_filtered_shared_genes.txt")
bc_scz_genes_in_pos_reg_filtered = fread("magma_smultixcan_filtered_shared_genes/bc_schizophrenia_filtered_shared_genes.txt")
bc_t2dm_genes_in_pos_reg_filtered = fread("magma_smultixcan_filtered_shared_genes/bc_t2dm_diamante_filtered_shared_genes.txt")
# create data frame with needed info
fig2c_data = data.frame(disease = c("high HDL", "Schizophrenia", "Depression", "high LDL", "Prostate cancer", "Type 2 diabetes"),
                        all = c(length(unique(bc_hdl_genes_in_pos_reg$entrezID)),
                                length(unique(bc_scz_genes_in_pos_reg$entrezID)),
                                length(unique(bc_depression_genes_in_pos_reg$entrezID)),
                                length(unique(bc_ldl_genes_in_pos_reg$entrezID)),
                                length(unique(bc_prostate_genes_in_pos_reg$entrezID)),
                                length(unique(bc_t2dm_genes_in_pos_reg$entrezID))),
                        magma_smultixcan_filtered = c(length(unique(bc_hdl_genes_in_pos_reg_filtered$entrezID)),
                                                      length(unique(bc_scz_genes_in_pos_reg_filtered$entrezID)),
                                                      length(unique(bc_depression_genes_in_pos_reg_filtered$entrezID)),
                                                      length(unique(bc_ldl_genes_in_pos_reg_filtered$entrezID)),
                                                      length(unique(bc_prostate_genes_in_pos_reg_filtered$entrezID)),
                                                      length(unique(bc_t2dm_genes_in_pos_reg_filtered$entrezID))))
fig2c_data = reshape2::melt(fig2c_data, "disease", colnames(fig2c_data)[2:ncol(fig2c_data)])
fig2c_data$disease = factor(fig2c_data$disease, levels = fig2c_data$disease, labels = fig2c_data$disease)
fig2c_data$variable = factor(fig2c_data$variable, levels = c("all", "magma_smultixcan_filtered"), labels = c("All", "After MAGMA &\nS-MultiXcan filtering"))
figure_2C = ggplot(fig2c_data, aes(x = value, y = disease, fill = variable)) +
  geom_col(position = "identity") +
  scale_fill_manual(values = c("#1ABC9C", "#00798C")) +
  labs(title = "\n\n") + # for aesthetics 
  xlab("\nNumber of genes\nin positively correlated loci with breast cancer") +
  ylab("") +
  scale_x_continuous(breaks = seq(0, 75, 10), limits = c(0, 75)) +
  theme_bw() +
  theme(plot.title = element_text(size = 44, color = "black", family = "Arial"),
        panel.grid.major.x = element_line(color = "lightgrey"),
        panel.grid.major.y = element_blank(),
        axis.text.x = element_text(size = 55, color = "black", family = "Arial"),
        axis.text.y = element_text(size = 55, color = "black", family = "Arial", face = "bold"),
        axis.title.x = element_text(size = 58, color = "black", family = "Arial"),
        legend.justification = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 55, color = "black", family = "Arial"),
        legend.key.size = unit(2,"line"), 
        legend.spacing.y = unit(0.5, 'cm')) +
  guides(fill = guide_legend(byrow = TRUE))
figure_2C
ggsave("figures/figure_2C.png", device = "png", width = 26, height = 14)

## full figure 2
figure_2_top_row = ggarrange(figure_2A, figure_2C, ncol = 2, nrow = 1, heights = c(0.5, 1), widths = c(1.5, 1)) +
  theme(plot.margin = margin(t = 0, r = 0, b = 80, l = 0)) 
figure_2_all = ggarrange(figure_2_top_row, figure_2B, ncol = 1, nrow = 2)
figure_2_all
ggsave("figures/figure_2_all.png", device = "png", width = 60, height = 40, limitsize = FALSE, dpi = 300)

rm(axis_label, axis.set, bc_depression_genes_in_pos_reg, chrom_length, chr_length_for_plot, bc_depression_genes_in_pos_reg_filtered, bc_depression_shared_loci, bc_hdl_genes_in_pos_reg, bc_hdl_genes_in_pos_reg_filtered, bc_hdl_shared_loci, bc_ldl_genes_in_pos_reg, bc_ldl_genes_in_pos_reg_filtered, bc_ldl_shared_loci, bc_prostate_genes_in_pos_reg, bc_prostate_genes_in_pos_reg_filtered, bc_prostate_shared_loci, bc_scz_genes_in_pos_reg, bc_scz_genes_in_pos_reg_filtered, bc_scz_shared_loci, bc_t2dm_genes_in_pos_reg, bc_t2dm_genes_in_pos_reg_filtered, bc_t2dm_shared_loci, fig2c_data, genes_in_string_high_conf, nCHR, p1, p2, shared_loci_ultimate, string_high_conf)

# ------------------------------------------------------------------------------------------------------------ #

### Figure 3 and S1 - Drugs indicated for the breast cancer clinically associated diseases that target the identified shared biology

### S1 - Connecting shared genetics to canonical pathways

## data
# canonical pathways significantly connected to the shared genes
bc_depression_sig_cps = fread("shared_genes_to_canonical_pathways/BC_DEPRESSION_CPs_pvalues_bonferroni.txt")
bc_hdl_sig_cps = fread("shared_genes_to_canonical_pathways/BC_DEPRESSION_CPs_pvalues_bonferroni.txt")
bc_ldl_sig_cps = fread("shared_genes_to_canonical_pathways/BC_LDL_CPs_pvalues_bonferroni.txt")
bc_prostate_sig_cps = fread("shared_genes_to_canonical_pathways/BC_PROSTATE_CPs_pvalues_bonferroni.txt")
bc_schizophrenia_sig_cps = fread("shared_genes_to_canonical_pathways/BC_SCZ_CPs_pvalues_bonferroni.txt")
bc_t2dm_sig_cps = fread("shared_genes_to_canonical_pathways/BC_T2DM_DIAMANTE_CPs_pvalues_bonferroni.txt")

## visualize KEGG only - rest of pathways are available as supplementary table
all_of_cps = unique(c(bc_depression_sig_cps$canonical_pathway, bc_hdl_sig_cps$canonical_pathway, bc_ldl_sig_cps$canonical_pathway,
                      bc_prostate_sig_cps$canonical_pathway, bc_schizophrenia_sig_cps$canonical_pathway, bc_t2dm_sig_cps$canonical_pathway))
## keep KEGG canonical pathways
kegg_cps = grepl("^KEGG_.*$", all_of_cps)
kegg_cps = all_of_cps[kegg_cps]
rm(all_of_cps)

## find which shared genes are connected to canonical pathways for each disease pair
# bc-depression
bc_depression_permuted = readRDS("shared_genes_to_canonical_pathways/permuted_BC_DEPRESSION.RDS")
# canonical pathways significantly connected to the shared genes
bc_depression_sig_cp_genes = bc_depression_permuted
for (i in 1:length(bc_depression_sig_cp_genes)) {
  bc_depression_sig_cp_genes[[i]] = bc_depression_sig_cp_genes[[i]] %>%
    mutate(ppr_pvalue = ppr_greater_than_obs / 1000) %>%
    mutate(ppr_pvalue_adj = p.adjust(ppr_pvalue, "bonferroni", length(ppr_pvalue))) %>%
    dplyr::select(shared_gene, ppr_pvalue_adj) %>%
    distinct()
  cat(i, "\n")
} ; rm(i)
bc_depression_sig_cp_genes = do.call(rbind, bc_depression_sig_cp_genes)
bc_depression_sig_cp_genes$canonical_pathway = rownames(bc_depression_sig_cp_genes) ; rownames(bc_depression_sig_cp_genes) = NULL
bc_depression_sig_cp_genes$canonical_pathway = gsub("\\..*$", "", bc_depression_sig_cp_genes$canonical_pathway)
temp = bc_depression_sig_cp_genes %>%
  dplyr::select(canonical_pathway, ppr_pvalue_adj) %>%
  distinct()
# adjust for the number of canonical pathways tested
temp$p_adjusted_final = p.adjust(temp$ppr_pvalue_adj, method = "bonferroni")
temp = temp %>%
  dplyr::select(canonical_pathway, pvalue_adj = p_adjusted_final) %>%
  filter(pvalue_adj < 0.05) %>%
  distinct()
bc_depression_sig_cp_genes = bc_depression_sig_cp_genes %>% filter(canonical_pathway %in% temp$canonical_pathway, ppr_pvalue_adj < 0.05)
rm(bc_depression_permuted, temp)
# keep KEGG pathways only
bc_depression_sig_cp_genes = bc_depression_sig_cp_genes[which(grepl("KEGG_", bc_depression_sig_cp_genes$canonical_pathway)), ] ; rownames(bc_depression_sig_cp_genes) = NULL
bc_depression_sig_cp_genes = bc_depression_sig_cp_genes %>% arrange(canonical_pathway, shared_gene) %>% dplyr::select(-ppr_pvalue_adj) %>% mutate(disease_pair = "BC-Depression", .before = shared_gene)

# bc-hdl
bc_hdl_permuted = readRDS("shared_genes_to_canonical_pathways/permuted_BC_HDL.RDS")
# canonical pathways significantly connected to the shared genes
bc_hdl_sig_cp_genes = bc_hdl_permuted
for (i in 1:length(bc_hdl_sig_cp_genes)) {
  bc_hdl_sig_cp_genes[[i]] = bc_hdl_sig_cp_genes[[i]] %>%
    mutate(ppr_pvalue = ppr_greater_than_obs / 1000) %>%
    mutate(ppr_pvalue_adj = p.adjust(ppr_pvalue, "bonferroni", length(ppr_pvalue))) %>%
    dplyr::select(shared_gene, ppr_pvalue_adj) %>%
    distinct()
  cat(i, "\n")
} ; rm(i)
bc_hdl_sig_cp_genes = do.call(rbind, bc_hdl_sig_cp_genes)
bc_hdl_sig_cp_genes$canonical_pathway = rownames(bc_hdl_sig_cp_genes) ; rownames(bc_hdl_sig_cp_genes) = NULL
bc_hdl_sig_cp_genes$canonical_pathway = gsub("\\..*$", "", bc_hdl_sig_cp_genes$canonical_pathway)
temp = bc_hdl_sig_cp_genes %>%
  dplyr::select(canonical_pathway, ppr_pvalue_adj) %>%
  distinct()
# adjust for the number of canonical pathways tested
temp$p_adjusted_final = p.adjust(temp$ppr_pvalue_adj, method = "bonferroni")
temp = temp %>%
  dplyr::select(canonical_pathway, pvalue_adj = p_adjusted_final) %>%
  filter(pvalue_adj < 0.05) %>%
  distinct()
bc_hdl_sig_cp_genes = bc_hdl_sig_cp_genes %>% filter(canonical_pathway %in% temp$canonical_pathway, ppr_pvalue_adj < 0.05)
rm(bc_hdl_permuted, temp)
# keep KEGG pathways only
bc_hdl_sig_cp_genes = bc_hdl_sig_cp_genes[which(grepl("KEGG_", bc_hdl_sig_cp_genes$canonical_pathway)), ] ; rownames(bc_hdl_sig_cp_genes) = NULL
bc_hdl_sig_cp_genes = bc_hdl_sig_cp_genes %>% arrange(canonical_pathway, shared_gene) %>% dplyr::select(-ppr_pvalue_adj) %>% mutate(disease_pair = "BC-HDL", .before = shared_gene)

# bc-ldl
bc_ldl_permuted = readRDS("shared_genes_to_canonical_pathways/permuted_BC_LDL.RDS")
# canonical pathways significantly connected to the shared genes
bc_ldl_sig_cp_genes = bc_ldl_permuted
for (i in 1:length(bc_ldl_sig_cp_genes)) {
  bc_ldl_sig_cp_genes[[i]] = bc_ldl_sig_cp_genes[[i]] %>%
    mutate(ppr_pvalue = ppr_greater_than_obs / 1000) %>%
    mutate(ppr_pvalue_adj = p.adjust(ppr_pvalue, "bonferroni", length(ppr_pvalue))) %>%
    dplyr::select(shared_gene, ppr_pvalue_adj) %>%
    distinct()
  cat(i, "\n")
} ; rm(i)
bc_ldl_sig_cp_genes = do.call(rbind, bc_ldl_sig_cp_genes)
bc_ldl_sig_cp_genes$canonical_pathway = rownames(bc_ldl_sig_cp_genes) ; rownames(bc_ldl_sig_cp_genes) = NULL
bc_ldl_sig_cp_genes$canonical_pathway = gsub("\\..*$", "", bc_ldl_sig_cp_genes$canonical_pathway)
temp = bc_ldl_sig_cp_genes %>%
  dplyr::select(canonical_pathway, ppr_pvalue_adj) %>%
  distinct()
# adjust for the number of canonical pathways tested
temp$p_adjusted_final = p.adjust(temp$ppr_pvalue_adj, method = "bonferroni")
temp = temp %>%
  dplyr::select(canonical_pathway, pvalue_adj = p_adjusted_final) %>%
  filter(pvalue_adj < 0.05) %>%
  distinct()
bc_ldl_sig_cp_genes = bc_ldl_sig_cp_genes %>% filter(canonical_pathway %in% temp$canonical_pathway, ppr_pvalue_adj < 0.05)
rm(bc_ldl_permuted, temp)
# keep KEGG pathways only
bc_ldl_sig_cp_genes = bc_ldl_sig_cp_genes[which(grepl("KEGG_", bc_ldl_sig_cp_genes$canonical_pathway)), ] ; rownames(bc_ldl_sig_cp_genes) = NULL
bc_ldl_sig_cp_genes = bc_ldl_sig_cp_genes %>% arrange(canonical_pathway, shared_gene) %>% dplyr::select(-ppr_pvalue_adj) %>% mutate(disease_pair = "BC-LDL", .before = shared_gene)

# bc-prostate
bc_prostate_permuted = readRDS("shared_genes_to_canonical_pathways/permuted_BC_PROSTATE.RDS")
# canonical pathways significantly connected to the shared genes
bc_prostate_sig_cp_genes = bc_prostate_permuted
for (i in 1:length(bc_prostate_sig_cp_genes)) {
  bc_prostate_sig_cp_genes[[i]] = bc_prostate_sig_cp_genes[[i]] %>%
    mutate(ppr_pvalue = ppr_greater_than_obs / 1000) %>%
    mutate(ppr_pvalue_adj = p.adjust(ppr_pvalue, "bonferroni", length(ppr_pvalue))) %>%
    dplyr::select(shared_gene, ppr_pvalue_adj) %>%
    distinct()
  cat(i, "\n")
} ; rm(i)
bc_prostate_sig_cp_genes = do.call(rbind, bc_prostate_sig_cp_genes)
bc_prostate_sig_cp_genes$canonical_pathway = rownames(bc_prostate_sig_cp_genes) ; rownames(bc_prostate_sig_cp_genes) = NULL
bc_prostate_sig_cp_genes$canonical_pathway = gsub("\\..*$", "", bc_prostate_sig_cp_genes$canonical_pathway)
temp = bc_prostate_sig_cp_genes %>%
  dplyr::select(canonical_pathway, ppr_pvalue_adj) %>%
  distinct()
# adjust for the number of canonical pathways tested
temp$p_adjusted_final = p.adjust(temp$ppr_pvalue_adj, method = "bonferroni")
temp = temp %>%
  dplyr::select(canonical_pathway, pvalue_adj = p_adjusted_final) %>%
  filter(pvalue_adj < 0.05) %>%
  distinct()
bc_prostate_sig_cp_genes = bc_prostate_sig_cp_genes %>% filter(canonical_pathway %in% temp$canonical_pathway, ppr_pvalue_adj < 0.05)
rm(bc_prostate_permuted, temp)
# keep KEGG pathways only
bc_prostate_sig_cp_genes = bc_prostate_sig_cp_genes[which(grepl("KEGG_", bc_prostate_sig_cp_genes$canonical_pathway)), ] ; rownames(bc_prostate_sig_cp_genes) = NULL
bc_prostate_sig_cp_genes = bc_prostate_sig_cp_genes %>% arrange(canonical_pathway, shared_gene) %>% dplyr::select(-ppr_pvalue_adj) %>% mutate(disease_pair = "BC-PC", .before = shared_gene)

# bc-scz
bc_scz_permuted = readRDS("shared_genes_to_canonical_pathways/permuted_BC_SCZ.RDS")
# canonical pathways significantly connected to the shared genes
bc_scz_sig_cp_genes = bc_scz_permuted
for (i in 1:length(bc_scz_sig_cp_genes)) {
  bc_scz_sig_cp_genes[[i]] = bc_scz_sig_cp_genes[[i]] %>%
    mutate(ppr_pvalue = ppr_greater_than_obs / 1000) %>%
    mutate(ppr_pvalue_adj = p.adjust(ppr_pvalue, "bonferroni", length(ppr_pvalue))) %>%
    dplyr::select(shared_gene, ppr_pvalue_adj) %>%
    distinct()
  cat(i, "\n")
} ; rm(i)
bc_scz_sig_cp_genes = do.call(rbind, bc_scz_sig_cp_genes)
bc_scz_sig_cp_genes$canonical_pathway = rownames(bc_scz_sig_cp_genes) ; rownames(bc_scz_sig_cp_genes) = NULL
bc_scz_sig_cp_genes$canonical_pathway = gsub("\\..*$", "", bc_scz_sig_cp_genes$canonical_pathway)
temp = bc_scz_sig_cp_genes %>%
  dplyr::select(canonical_pathway, ppr_pvalue_adj) %>%
  distinct()
# adjust for the number of canonical pathways tested
temp$p_adjusted_final = p.adjust(temp$ppr_pvalue_adj, method = "bonferroni")
temp = temp %>%
  dplyr::select(canonical_pathway, pvalue_adj = p_adjusted_final) %>%
  filter(pvalue_adj < 0.05) %>%
  distinct()
bc_scz_sig_cp_genes = bc_scz_sig_cp_genes %>% filter(canonical_pathway %in% temp$canonical_pathway, ppr_pvalue_adj < 0.05)
rm(bc_scz_permuted, temp)
# keep KEGG pathways only
bc_scz_sig_cp_genes = bc_scz_sig_cp_genes[which(grepl("KEGG_", bc_scz_sig_cp_genes$canonical_pathway)), ] ; rownames(bc_scz_sig_cp_genes) = NULL
bc_scz_sig_cp_genes = bc_scz_sig_cp_genes %>% arrange(canonical_pathway, shared_gene) %>% dplyr::select(-ppr_pvalue_adj) %>% mutate(disease_pair = "BC-SCZ", .before = shared_gene)

# bc-t2dm
bc_t2dm_permuted = readRDS("shared_genes_to_canonical_pathways/permuted_BC_T2DM_DIAMANTE.RDS")
# canonical pathways significantly connected to the shared genes
bc_t2dm_sig_cp_genes = bc_t2dm_permuted
for (i in 1:length(bc_t2dm_sig_cp_genes)) {
  bc_t2dm_sig_cp_genes[[i]] = bc_t2dm_sig_cp_genes[[i]] %>%
    mutate(ppr_pvalue = ppr_greater_than_obs / 1000) %>%
    mutate(ppr_pvalue_adj = p.adjust(ppr_pvalue, "bonferroni", length(ppr_pvalue))) %>%
    dplyr::select(shared_gene, ppr_pvalue_adj) %>%
    distinct()
  cat(i, "\n")
} ; rm(i)
bc_t2dm_sig_cp_genes = do.call(rbind, bc_t2dm_sig_cp_genes)
bc_t2dm_sig_cp_genes$canonical_pathway = rownames(bc_t2dm_sig_cp_genes) ; rownames(bc_t2dm_sig_cp_genes) = NULL
bc_t2dm_sig_cp_genes$canonical_pathway = gsub("\\..*$", "", bc_t2dm_sig_cp_genes$canonical_pathway)
temp = bc_t2dm_sig_cp_genes %>%
  dplyr::select(canonical_pathway, ppr_pvalue_adj) %>%
  distinct()
# adjust for the number of canonical pathways tested
temp$p_adjusted_final = p.adjust(temp$ppr_pvalue_adj, method = "bonferroni")
temp = temp %>%
  dplyr::select(canonical_pathway, pvalue_adj = p_adjusted_final) %>%
  filter(pvalue_adj < 0.05) %>%
  distinct()
bc_t2dm_sig_cp_genes = bc_t2dm_sig_cp_genes %>% filter(canonical_pathway %in% temp$canonical_pathway, ppr_pvalue_adj < 0.05)
rm(bc_t2dm_permuted, temp)
# keep KEGG pathways only
bc_t2dm_sig_cp_genes = bc_t2dm_sig_cp_genes[which(grepl("KEGG_", bc_t2dm_sig_cp_genes$canonical_pathway)), ] ; rownames(bc_t2dm_sig_cp_genes) = NULL
bc_t2dm_sig_cp_genes = bc_t2dm_sig_cp_genes %>% arrange(shared_gene, canonical_pathway, shared_gene) %>% dplyr::select(-ppr_pvalue_adj) %>% mutate(disease_pair = "BC-T2DM", .before = shared_gene)

## create matrix input for heatmap
bc_depression_sig_cp_genes
bc_hdl_sig_cp_genes
length(c(unique(bc_depression_sig_cp_genes$shared_gene), 
         unique(bc_hdl_sig_cp_genes$shared_gene),
         unique(bc_ldl_sig_cp_genes$shared_gene),
         unique(bc_prostate_sig_cp_genes$shared_gene),
         unique(bc_scz_sig_cp_genes$shared_gene),
         unique(bc_t2dm_sig_cp_genes$shared_gene))) # 35 --> these must be the columns of the heatmap
## The following genes are shared between more than one BC-disease pair
## BC-LDL / BC-T2DM: 7059, 55974 
## BC-T2DM / BC-Prostate: 6934
## BC-SCZ / BC-Prostate: 1463, 23383, 374887

kegg_cps_heatmap_data = matrix(ncol = 35, nrow = length(kegg_cps)) # these KEGG CPs are the ones connected to at least one of the shared genes for any disease pair
# colnames(kegg_cps_heatmap_data) = c("BRCA-Depression", "BRCA-high HDL", "BRCA-high LDL", "BRCA-PC", "BRCA-SCZ", "BRCA-T2DM")
colnames(kegg_cps_heatmap_data) = c(unique(bc_depression_sig_cp_genes$shared_gene), 
                                    unique(bc_hdl_sig_cp_genes$shared_gene),
                                    unique(bc_ldl_sig_cp_genes$shared_gene),
                                    unique(bc_prostate_sig_cp_genes$shared_gene),
                                    unique(bc_scz_sig_cp_genes$shared_gene),
                                    unique(bc_t2dm_sig_cp_genes$shared_gene))
rownames(kegg_cps_heatmap_data) = kegg_cps
# add one row at the end of the matrix to help populate it
kegg_cps_heatmap_data = rbind(kegg_cps_heatmap_data, 
                              disease_pair = c(rep("BRCA-Depression", length(unique(bc_depression_sig_cp_genes$shared_gene))), 
                                               rep("BRCA-high HDL", length(unique(bc_hdl_sig_cp_genes$shared_gene))),
                                               rep("BRCA-high LDL", length(unique(bc_ldl_sig_cp_genes$shared_gene))),
                                               rep("BRCA-PC", length(unique(bc_prostate_sig_cp_genes$shared_gene))),
                                               rep("BRCA-SCZ", length(unique(bc_scz_sig_cp_genes$shared_gene))),
                                               rep("BRCA-T2DM", length(unique(bc_t2dm_sig_cp_genes$shared_gene)))))

for (row in 1:(nrow(kegg_cps_heatmap_data) - 1)) {
  for (col in 1:ncol(kegg_cps_heatmap_data)) {
    
    cp = rownames(kegg_cps_heatmap_data)[row]
    shared_gene_temp = colnames(kegg_cps_heatmap_data)[col]
    disease_pair = kegg_cps_heatmap_data[nrow(kegg_cps_heatmap_data), col]
    
    if (disease_pair == "BRCA-Depression") {
      if (cp %in% bc_depression_sig_cp_genes$canonical_pathway) {
        temp = bc_depression_sig_cp_genes %>% filter(shared_gene == shared_gene_temp, canonical_pathway == cp)
        if (nrow(temp) > 0) {
          kegg_cps_heatmap_data[row, col] = 1
        } else {
          kegg_cps_heatmap_data[row, col] = 0
        }
      } else {
        kegg_cps_heatmap_data[row, col] = 0
      }
    }
    
    if (disease_pair == "BRCA-high HDL") {
      if (cp %in% bc_hdl_sig_cp_genes$canonical_pathway) {
        temp = bc_hdl_sig_cp_genes %>% filter(shared_gene == shared_gene_temp, canonical_pathway == cp)
        if (nrow(temp) > 0) {
          kegg_cps_heatmap_data[row, col] = 1
        } else {
          kegg_cps_heatmap_data[row, col] = 0
        }
      } else {
        kegg_cps_heatmap_data[row, col] = 0
      }
    }
    
    if (disease_pair == "BRCA-high LDL") {
      if (cp %in% bc_ldl_sig_cp_genes$canonical_pathway) {
        temp = bc_ldl_sig_cp_genes %>% filter(shared_gene == shared_gene_temp, canonical_pathway == cp)
        if (nrow(temp) > 0) {
          kegg_cps_heatmap_data[row, col] = 1
        } else {
          kegg_cps_heatmap_data[row, col] = 0
        }
      } else {
        kegg_cps_heatmap_data[row, col] = 0
      }
    }
    
    if (disease_pair == "BRCA-PC") {
      if (cp %in% bc_prostate_sig_cp_genes$canonical_pathway) {
        temp = bc_prostate_sig_cp_genes %>% filter(shared_gene == shared_gene_temp, canonical_pathway == cp)
        if (nrow(temp) > 0) {
          kegg_cps_heatmap_data[row, col] = 1
        } else {
          kegg_cps_heatmap_data[row, col] = 0
        }
      } else {
        kegg_cps_heatmap_data[row, col] = 0
      }
    }
    
    if (disease_pair == "BRCA-SCZ") {
      if (cp %in% bc_scz_sig_cp_genes$canonical_pathway) {
        temp = bc_scz_sig_cp_genes %>% filter(shared_gene == shared_gene_temp, canonical_pathway == cp)
        if (nrow(temp) > 0) {
          kegg_cps_heatmap_data[row, col] = 1
        } else {
          kegg_cps_heatmap_data[row, col] = 0
        }
      } else {
        kegg_cps_heatmap_data[row, col] = 0
      }
    }
    
    if (disease_pair == "BRCA-T2DM") {
      if (cp %in% bc_t2dm_sig_cp_genes$canonical_pathway) {
        temp = bc_t2dm_sig_cp_genes %>% filter(shared_gene == shared_gene_temp, canonical_pathway == cp)
        if (nrow(temp) > 0) {
          kegg_cps_heatmap_data[row, col] = 1
        } else {
          kegg_cps_heatmap_data[row, col] = 0
        }
      } else {
        kegg_cps_heatmap_data[row, col] = 0
      }
    }
  }
}
# remove last row
kegg_cps_heatmap_data = kegg_cps_heatmap_data[-nrow(kegg_cps_heatmap_data), ]
row_names = rownames(kegg_cps_heatmap_data)
kegg_cps_heatmap_data = apply(kegg_cps_heatmap_data, 2 ,as.numeric)
rownames(kegg_cps_heatmap_data) = row_names ; rm(row_names)

## combine KEGG pathways that are similar and connected to shared genetics of same disease pairs
which(grepl("HEDGEHOG_SIGNALING_PATHWAY", rownames(kegg_cps_heatmap_data)))
kegg_cps_heatmap_data[c(1,2,3,4,111), ] # 5727, 2736
rownames(kegg_cps_heatmap_data)[c(1,2,3,4,111)] = "KEGG_HEDGEHOG_SIGNALING_PATHWAY"
kegg_cps_heatmap_data[c(1,2,3,4,111), which(colnames(kegg_cps_heatmap_data) == "5727")] = 1 ; kegg_cps_heatmap_data[c(1,2,3,4,111), which(colnames(kegg_cps_heatmap_data) == "2736")] = 1
kegg_cps_heatmap_data = kegg_cps_heatmap_data[-c(2,3,4,111), ]

which(grepl("PI3K_SIGNALING_PATHWAY", rownames(kegg_cps_heatmap_data)))
kegg_cps_heatmap_data[c(39,46,47,54,57,58,59,60,62,65,66,69,71,74,75,77,80,81,82,84,85,86,87,88,89,90,91,92,93,96), ] # 4824
rownames(kegg_cps_heatmap_data)[c(39,46,47,54,57,58,59,60,62,65,66,69,71,74,75,77,80,81,82,84,85,86,87,88,89,90,91,92,93,96)] = "KEGG_MEDICUS_PI3K_SIGNALING_PATHWAY"
kegg_cps_heatmap_data = kegg_cps_heatmap_data[-c(46,47,54,57,58,59,60,62,65,66,69,71,74,75,77,80,81,82,84,85,86,87,88,89,90,91,92,93,96), ]

which(grepl("CELL_CYCLE", rownames(kegg_cps_heatmap_data)))
kegg_cps_heatmap_data[c(34,40,44,54,55,56,60,61,62,63,64), ] # 1030
rownames(kegg_cps_heatmap_data)[c(34,40,44,54,55,56,60,61,62,63,64)] = "KEGG_CELL_CYCLE"
kegg_cps_heatmap_data = kegg_cps_heatmap_data[-c(40,44,54,55,56,60,61,62,63,64), ]

which(grepl("WNT_SIGNALING", rownames(kegg_cps_heatmap_data)))
which(grepl("WNT5A_ROR_SIGNALING_PATHWAY", rownames(kegg_cps_heatmap_data)))
kegg_cps_heatmap_data[c(5,6,7,8,42,56,57), ] # 8326, 6934, 6934
rownames(kegg_cps_heatmap_data)[c(5,6,7,8,42,56,57)] = "KEGG_MEDICUS_WNT_SIGNALING_PATHWAY"
kegg_cps_heatmap_data[c(5,6,7,8,42,56,57), which(colnames(kegg_cps_heatmap_data) == "8326")] = 1 ; kegg_cps_heatmap_data[c(5,6,7,8,42,56,57), which(colnames(kegg_cps_heatmap_data) == "6934")] = 1
kegg_cps_heatmap_data = kegg_cps_heatmap_data[-c(6,7,8,42,56,57), ]

which(grepl("AUTOPHAGY_VESICLE_NUCLEATION", rownames(kegg_cps_heatmap_data)))
kegg_cps_heatmap_data[c(20,21,22,30),] # 987, 25989
rownames(kegg_cps_heatmap_data)[c(20,21,22,30)] = "KEGG_MEDICUS_AUTOPHAGY_VESICLE_NUCLEATION_ELONGATION_MATURATION"
kegg_cps_heatmap_data[c(20,21,22,30), which(colnames(kegg_cps_heatmap_data) == "987")] = 1 ; kegg_cps_heatmap_data[c(20,21,22,30), which(colnames(kegg_cps_heatmap_data) == "25989")] = 1
kegg_cps_heatmap_data = kegg_cps_heatmap_data[-c(21,22,30), ]

which(grepl("PRNP_PI3K_NOX2_SIGNALING_PATHWAY", rownames(kegg_cps_heatmap_data)))
kegg_cps_heatmap_data[c(64,65),] # 653361
rownames(kegg_cps_heatmap_data)[c(64,65)] = "KEGG_MEDICUS_PRNP_PI3K_NOX2_SIGNALING_PATHWAY"
kegg_cps_heatmap_data = kegg_cps_heatmap_data[-c(65), ]

which(grepl("BETA_CATENIN_SIGNALING_PATHWAY", rownames(kegg_cps_heatmap_data)))
kegg_cps_heatmap_data[c(36,48),] # 6934, 6934
rownames(kegg_cps_heatmap_data)[c(36,48)] = "KEGG_MEDICUS_BETA_CATENIN_SIGNALING_PATHWAY"
kegg_cps_heatmap_data = kegg_cps_heatmap_data[-c(48), ]

# make pathway names shorter
rownames(kegg_cps_heatmap_data) = gsub("_MEDICUS_REFERENCE_", replacement = "_MEDICUS_", rownames(kegg_cps_heatmap_data))
rownames(kegg_cps_heatmap_data) = gsub("_MEDICUS_PATHOGEN_", replacement = "_MEDICUS_", rownames(kegg_cps_heatmap_data))

# convert entrezIDs to gene names
colnames(kegg_cps_heatmap_data) = c("PTCH1", "FZD9", "CYP7A1", "EFNA3", "THBS3", "MPI", "CSPG4", "CYP1A2", "SLC50A1", "LRBA", "ULK3", "CLK3", "GCNT4", "CDKN2B", "NCAN", "NKX3-1", "MAU2", "TCF7L2",
                                    "ZW10", "DFFA", "YJEFN3", "NCAN", "MAU2", "YJEFN3", "DNMT3A", "FDPS", "GBA1", "GLI2", "TCF7L2", "THBS3", "PLA2G6", "SLC50A1", "BAIAP2L2", "TOM1L2", "NCF1")

## heatmap
colors = structure(c("white", "gray70"), names = c(0, 1))
top_ha = HeatmapAnnotation("BRCA clinically associated disease" = c("Depression", "high HDL", rep("high LDL",11), rep("Prostate cancer",8), rep("Schizophrenia",3), rep("Type 2 diabetes",11)),
                           show_annotation_name = FALSE,
                           col = list("BRCA clinically associated disease" = c("Depression" = "#6a3d9a", "high HDL" = "deeppink4", "high LDL" = "#33a02c",
                                                                               "Prostate cancer" = "#e31a1c", "Schizophrenia" = "#ff7f00", "Type 2 diabetes" = "#1f78b4")), show_legend = FALSE)
kegg_hp = Heatmap(kegg_cps_heatmap_data,
                  name = "\n\t",
                  row_names_side = "left",
                  column_names_side = "top",
                  col = colors,
                  row_dend_reorder = TRUE,
                  cluster_rows = TRUE,
                  show_row_dend = FALSE,
                  column_dend_reorder = FALSE,
                  show_column_dend = FALSE,
                  cluster_columns = FALSE,
                  show_heatmap_legend = FALSE,
                  rect_gp = gpar(col = "black", lwd = 0.4),
                  width = unit(40, "cm"),
                  height = unit(60, "cm"),
                  row_names_gp = gpar(fontsize = 20),
                  column_names_gp = gpar(fontsize = 20),
                  column_names_rot = 45, 
                  top_annotation = top_ha)
lgd_list = list(
  Legend(labels = c("Depression", "high HDL", "high LDL", "Prostate cancer", "Schizophrenia", "Type 2 diabetes"),
         title = "BRCA clinically associated disease",
         background = "white", 
         border = TRUE,
         legend_gp = gpar(fill = c("#6a3d9a", "deeppink4", "#33a02c", "#e31a1c", "#ff7f00", "#1f78b4")),
         labels_gp = gpar(fontsize = 20),
         title_gp = gpar(col = "black", fontsize = 20, fontface = "bold"),
         grid_height = unit(0.7, "cm"),
         grid_width = unit(0.7, "cm"),
         row_gap = unit(0.1, "cm"),
         title_gap = unit(0.4, "cm")),
  Legend(labels = c("Yes", "No"),
         title = "Shared gene connected to canonical pathway",
         background = "white", 
         border = TRUE,
         legend_gp = gpar(fill = c("gray70", "white")),
         labels_gp = gpar(fontsize = 20),
         title_gp = gpar(col = "black", fontsize = 20, fontface = "bold"),
         grid_height = unit(0.7, "cm"),
         grid_width = unit(0.7, "cm"),
         row_gap = unit(0.1, "cm"),
         title_gap = unit(0.4, "cm"))
)
tiff("figures/figure_S1.tiff",
     width = 110, height = 70, units = "cm", 
     res = 300, compression = "lzw")
draw(kegg_hp, annotation_legend_list = lgd_list, annotation_legend_side = "right")
dev.off()

## Figure 3
shared_genes_to_canonical_pathways = function(permuted_values_path, shared_genes, disease) {
  permuted_values = readRDS(permuted_values_path)
  sig_cp = permuted_values
  for (i in 1:length(sig_cp)) {
    sig_cp[[i]] = sig_cp[[i]] %>%
      mutate(ppr_pvalue = ppr_greater_than_obs / 1000) %>%
      # keep the minimum p-value for each canonical pathway, after adjusting for the number of shared genes tested
      mutate(ppr_pvalue_adj = p.adjust(ppr_pvalue, "bonferroni", length(ppr_pvalue))) %>%
      dplyr::select(shared_gene, ppr_pvalue_adj) %>%
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
  sig_cp = sig_cp %>% dplyr::select(shared_gene = gene_name, connected_canonical_pathway = "Connected canonical pathway") %>% arrange(shared_gene)
  return(sig_cp)
}

bc_hdl_sig_cp = shared_genes_to_canonical_pathways(permuted_values_path = "shared_genes_to_canonical_pathways/permuted_BC_HDL.RDS", 
                                                   shared_genes = fread("magma_smultixcan_filtered_shared_genes/bc_hdl_filtered_shared_genes.txt"),
                                                   disease = "high HDL")
bc_hdl_sig_cp$connected_canonical_pathway = gsub("KEGG_MEDICUS_REFERENCE_WNT5A_ROR_SIGNALING_PATHWAY", "KEGG_WNT5A_ROR_SIGNALING_PATHWAY", bc_hdl_sig_cp$connected_canonical_pathway)
bc_hdl_sig_cp$connected_canonical_pathway = gsub("KEGG_MEDICUS_REFERENCE_WNT_SIGNALING_MODULATION_SOST_LRP4", "KEGG_WNT_SIGNALING_MODULATION_SOST_LRP4", bc_hdl_sig_cp$connected_canonical_pathway)
bc_hdl_sig_cp$connected_canonical_pathway = gsub("KEGG_MEDICUS_REFERENCE_WNT_SIGNALING_PATHWAY", "KEGG_WNT_SIGNALING_PATHWAY", bc_hdl_sig_cp$connected_canonical_pathway)
bc_hdl_sig_cp$connected_canonical_pathway = gsub("KEGG_MEDICUS_VARIANT_LRP6_OVEREXPRESSION_TO_WNT_SIGNALING_PATHWAY", "KEGG_LRP6_OVEREXPRESSION_TO_WNT_SIGNALING_PATHWAY", bc_hdl_sig_cp$connected_canonical_pathway)

hdl_drugs = fread("drugs_to_shared_canonical_pathways/BC_HDL_drug_targets_perm_pvalues_bonferroni.txt") %>% dplyr::select(drug) %>% distinct()
bc_hdl_drugs_cp = fread("drugs_to_shared_canonical_pathways/BC_HDL_drug_targets_perm_pvalues_bonferroni.txt") %>% filter(perm_pvalue_adj < 0.05)
bc_hdl_drugs_cp$canonical_pathway = gsub("KEGG_MEDICUS_REFERENCE_WNT5A_ROR_SIGNALING_PATHWAY", "KEGG_WNT5A_ROR_SIGNALING_PATHWAY", bc_hdl_drugs_cp$canonical_pathway)
bc_hdl_drugs_cp$canonical_pathway = gsub("KEGG_MEDICUS_REFERENCE_WNT_SIGNALING_MODULATION_SOST_LRP4", "KEGG_WNT_SIGNALING_MODULATION_SOST_LRP4", bc_hdl_drugs_cp$canonical_pathway)
bc_hdl_drugs_cp$canonical_pathway = gsub("KEGG_MEDICUS_REFERENCE_WNT_SIGNALING_PATHWAY", "KEGG_WNT_SIGNALING_PATHWAY", bc_hdl_drugs_cp$canonical_pathway)
bc_hdl_drugs_cp$canonical_pathway = gsub("KEGG_MEDICUS_VARIANT_LRP6_OVEREXPRESSION_TO_WNT_SIGNALING_PATHWAY", "KEGG_LRP6_OVEREXPRESSION_TO_WNT_SIGNALING_PATHWAY", bc_hdl_drugs_cp$canonical_pathway)

# filter for cps that are also connected to a drug
bc_hdl_sig_cp = bc_hdl_sig_cp %>% filter(connected_canonical_pathway %in% bc_hdl_drugs_cp$canonical_pathway) %>% dplyr::select(edge_start = shared_gene, edge_stop = connected_canonical_pathway) %>% distinct()
bc_hdl_drugs_cp = bc_hdl_drugs_cp %>% dplyr::select(edge_start = drug, edge_stop = canonical_pathway) %>% distinct()
edges = rbind(bc_hdl_sig_cp, bc_hdl_drugs_cp)
graph = graph_from_data_frame(edges, directed = TRUE)  

# create a 3 layer graph
v_layers_df = unique( rbind(
  expand.grid( ID = bc_hdl_sig_cp$edge_start, Layer = 1),
  expand.grid( ID = bc_hdl_sig_cp$edge_stop, Layer = 2),
  expand.grid( ID = bc_hdl_drugs_cp$edge_start, Layer = 3),
  expand.grid( ID = bc_hdl_drugs_cp$edge_stop, Layer = 2)
))
v_layers = setNames( v_layers_df$Layer, v_layers_df$ID)
V(graph)$layer = v_layers[V(graph)$name]
layout.k_partite = function(g) {
  l = layout.sugiyama(g)$layout[,2:1]
  l[,1] = V(g)$layer
  l[,2] = - l[,2] + 1 + max(l[,2])
  l
}
V(graph)$label.cex <- c(1,1,1,1,1,1,1,0.7,0.7,0.7,0.7,0.7)
V(graph)$shape <- c("circle","circle","rectangle","rectangle","rectangle","rectangle","rectangle","raster","raster","raster","raster","raster")
plot(graph, layout = layout.k_partite(graph),
     # vertex.shape = "raster",
     vertex.frame.color = "black",
     vertex.size = 15,
     vertex.color = c("orange", "orange", "green4", "green4", "green4", "green4", "green4", "steelblue", "steelblue", "steelblue" ,"steelblue", "steelblue"),
     
     vertex.label.color = "black",
     vertex.label.family = "Arial",                  # Font family of the label (e.g.“Times”, “Helvetica”)  
     vertex.label.font = c(4,4,2,2,2,2,2,2,2,2,2,2),
     vertex.label.dist = 2,
     vertex.label.degree = 300,
     
     edge.color = "gray60",
     edge.width = 2,
     edge.arrow.size = 0,
     edge.arrow.width = 0,
     edge.curved = 0.3)
# manual save
# 850 x 100

rm(shared_genes_to_canonical_pathways, bc_hdl_sig_cp, hdl_drugs, edges, graph, v_layers_df, v_layers, layout.k_partite)
rm(bc_depression_sig_cp_genes, bc_depression_sig_cps, bc_hdl_drug_to_sig_cp, bc_hdl_sig_cp_genes, bc_hdl_sig_cps, bc_hdl_temp, bc_ldl_sig_cp_genes, bc_ldl_sig_cps, bc_prostate_sig_cp_genes, bc_prostate_sig_cps, bc_schizophrenia_sig_cps, bc_scz_sig_cp_genes, bc_t2dm_sig_cp_genes, bc_t2dm_sig_cps, kegg_cps_heatmap_data, temp, col , colors, cp, disease_pair, drug_temp, kegg_cps, row, shared_gene_temp)
# ------------------------------------------------------------------------------------------------ #

### Figure 4: evaluation of drug repurposing and prioritization of new candidates

## drugs indicated for each disease
depression_drugs = fread("preprocessed_data/depression_drug_targets.txt") %>% dplyr::select(drug) %>% distinct()
hdl_drugs = fread("preprocessed_data/hdl_drug_targets.txt") %>% dplyr::select(drug) %>% distinct()
ldl_drugs = fread("preprocessed_data/ldl_drug_targets.txt") %>% dplyr::select(drug) %>% distinct()
prostate_drugs = fread("preprocessed_data/prostate_drug_targets.txt") %>% dplyr::select(drug) %>% distinct()
scz_drugs = fread("preprocessed_data/schizophrenia_drug_targets.txt") %>% dplyr::select(drug) %>% distinct()
t2dm_drugs = fread("preprocessed_data/t2dm_drug_targets.txt") %>% dplyr::select(drug) %>% distinct()
total_drugs = data.frame(disease = c("Depression", "high HDL", "high LDL", "Prostate cancer", "Schizophrenia", "Type 2 diabetes"),
                         nr_approved_drugs = c(nrow(depression_drugs),nrow(hdl_drugs),nrow(ldl_drugs),nrow(prostate_drugs),nrow(scz_drugs),nrow(t2dm_drugs)))
total_drugs_unique = unique(rbind(depression_drugs, hdl_drugs, ldl_drugs, prostate_drugs, scz_drugs, t2dm_drugs))
rm(depression_drugs, hdl_drugs, ldl_drugs, prostate_drugs, scz_drugs, t2dm_drugs)
total_drugs = total_drugs %>% arrange(nr_approved_drugs)
total_drugs$disease = factor(total_drugs$disease, levels = total_drugs$disease, labels = total_drugs$disease)
ggplot(total_drugs, aes(x = nr_approved_drugs, y = disease)) +
  geom_col(fill = "#EFC1A1", color = "black") + 
  xlab("\nNumber of drugs") +
  ylab("") +
  scale_x_continuous(breaks = seq(0,30,5), limits = c(0, 31), expand = c(0.01,0)) +
  theme_bw() +
  theme(plot.title = element_text(size = 26, color = "black", family = "Arial"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y.right = element_text(size = 30, color = "black", family = "Arial", face = "bold", hjust = 0.5,
                                         margin = margin(t = 0, r = 0, b = 0, l = 40)),
        axis.text = element_text(size = 50, color = "black", family = "Arial"),
        axis.title.x = element_text(size = 60, color = "black", family = "Arial"))
ggsave("figure4A_nr_drugs_other_diseases.png", device = "png", path = "figures/",
       width = 16, height = 8, dpi = 300)

# investigated and approved drugs for BRCA
bc_inv_ind_drugs = fread("preprocessed_data/bc_indicated_investigated_drugs_all.txt")
bc_drugs_per_phase = data.frame(phase = c("Approved", "Phase III",  "Phase II", "Phase I", "Pre-clinical"),
                                nr_drugs = c(table(bc_inv_ind_drugs$max_phase)[5],
                                             table(bc_inv_ind_drugs$max_phase)[4],
                                             table(bc_inv_ind_drugs$max_phase)[3],
                                             table(bc_inv_ind_drugs$max_phase)[2],
                                             table(bc_inv_ind_drugs$max_phase)[1]))
bc_drugs_per_phase$phase = factor(bc_drugs_per_phase$phase, levels = bc_drugs_per_phase$phase, labels = bc_drugs_per_phase$phase) # these drugs are from 4,237 clinical trials (the rest of them are Phase 4 --> already approved)
ggplot(bc_drugs_per_phase, aes(x = nr_drugs, y = phase)) +
  geom_col(fill = "#B7CDA4", color = "black") + 
  xlab("\nNumber of drugs") +
  ylab("") +
  scale_x_continuous(breaks = seq(0,300,50), limits = c(0, 305), expand = c(0.01, 5)) +
  theme_bw() +
  theme(plot.title = element_text(size = 26, color = "black", family = "Arial"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y.right = element_text(size = 30, color = "black", family = "Arial", face = "bold", hjust = 0.5,
                                         margin = margin(t = 0, r = 0, b = 0, l = 40)),
        axis.text = element_text(size = 50, color = "black", family = "Arial"),
        axis.title.x = element_text(size = 60, color = "black", family = "Arial"))
ggsave("figure4A_nr_drugs_bc.png", device = "png", path = "figures/",
       width = 16, height = 8, dpi = 300)

# drugs for the other diseases that are also investigated/indicated for breast cancer
bc_inv_ind_diseases_drugs_overlap = fread("preprocessed_data/bc_indicated_investigated_drugs_of_predisposing_diseases.txt")

# for venn diagram
nrow(total_drugs_unique) # 112 unique drugs approved for the BC clinically associated diseases
nrow(bc_inv_ind_drugs) # 773 unique drugs investigated or indicated for BC
nrow(bc_inv_ind_diseases_drugs_overlap) # 16 unique drugs are approved for the BC clinically associated diseases and are also investigated/indicated for breast cancer
table(bc_inv_ind_diseases_drugs_overlap$max_phase) # the 16 drugs per phase

## Figure 4B
depression_drug_targets = fread("drugs_to_shared_canonical_pathways/BC_DEPRESSION_drug_targets_perm_pvalues_bonferroni.txt")
hdl_drug_targets = fread("drugs_to_shared_canonical_pathways/BC_HDL_drug_targets_perm_pvalues_bonferroni.txt")
ldl_drug_targets = fread("drugs_to_shared_canonical_pathways/BC_LDL_drug_targets_perm_pvalues_bonferroni.txt")
prostate_drug_targets = fread("drugs_to_shared_canonical_pathways/BC_PROSTATE_drug_targets_perm_pvalues_bonferroni.txt")
scz_drug_targets = fread("drugs_to_shared_canonical_pathways/BC_SCZ_drug_targets_perm_pvalues_bonferroni.txt")
t2dm_drug_targets = fread("drugs_to_shared_canonical_pathways/BC_T2DM_drug_targets_perm_pvalues_bonferroni.txt")
# find candidate drugs for repurposing by looking if a drug is significantly connected to any of the CPs that are significantly connected to the shared genes
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
rm(i, drug_sig_cps)
# combine
drugs_PAall = rbind(bc_depression_recommended_drugs, bc_hdl_recommended_drugs, bc_ldl_recommended_drugs, bc_prostate_recommended_drugs, bc_scz_recommended_drugs, bc_t2dm_recommended_drugs)
drugs_PAall = drugs_PAall %>%
  group_by(drug) %>% 
  mutate(recommended = max(recommended)) %>%
  ungroup() %>%
  distinct()
# from thhe above drugs, ones that are investigater/indicated for bc
bc_inv_ind_diseases_drugs_overlap$drug_name = tolower(bc_inv_ind_diseases_drugs_overlap$drug_name)

# data for plot
figure_4b_data = data.frame(drug = drugs_PAall$drug, repurposing_candidate = drugs_PAall$recommended)
figure_4b_data$investigated_indicated_for_bc = ifelse(figure_4b_data$drug %in% bc_inv_ind_diseases_drugs_overlap$drug_name, 1, 0)
figure_4b_data$phase = left_join(figure_4b_data, bc_inv_ind_diseases_drugs_overlap[, c(2,4)], by = c("drug" = "drug_name"))
figure_4b_data$repurposing_candidate = factor(figure_4b_data$repurposing_candidate, levels = c(1, 0), labels = c("Candidate drugs\nfor repurposing", "Not candidate drugs\nfor repurposing"))
figure_4b_data$investigated_indicated_for_bc = factor(figure_4b_data$investigated_indicated_for_bc, levels = c(1, 0), labels = c("Yes", "No"))
ggplot(figure_4b_data, aes(x = repurposing_candidate, fill = investigated_indicated_for_bc)) +
  geom_bar(position = "dodge", color = "white", width = 0.6) +
  xlab("") +
  ylab("Number of drugs\n") +
  scale_y_continuous(breaks = seq(0, 70, 10), limits = c(0,63)) + 
  scale_fill_manual("Investigated/Approved for BRCA", values = c("#1ABC9C", "#F8766D")) +
  geom_label(label = "Novel candidate drugs for repurposing", x=1.2, y=63,
             label.padding = unit(0.4, "lines"), # Rectangle size around label
             size = 10,
             color = "white",
             fill = "red", alpha = 0.05, fontface = "bold") +
  theme_bw() +
  theme(plot.title = element_text(size = 48, color = "black", family = "Arial"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.y.right = element_text(size = 44, color = "black", family = "Arial", face = "bold", hjust = 0.5,
                                         margin = margin(t = 0, r = 0, b = 0, l = 40)),
        axis.text = element_text(size = 40, color = "black", family = "Arial"),
        axis.title.y = element_text(size = 40, color = "black", family = "Arial"),
        legend.text = element_text(size = 40, color = "black", family = "Arial"),
        legend.title = element_text(size = 44, color = "black", family = "Arial"), legend.position = "top")
ggsave("figure4B.png", device = "png", path = "figures/", 
       width = 14, height = 9, dpi = 300)

## drug development phase of drugs we recommend and are investigated/approved for breast cancer
bc_inv_ind_diseases_drugs_overlap$max_phase = factor(bc_inv_ind_diseases_drugs_overlap$max_phase, levels = c(4,3,2,1), labels = c("Approved", "Phase III", "Phase II", "Phase I"))
# keep the ones we recommend
candidate_drugs_repurposing = drugs_PAall %>% filter(recommended == 1)
bc_inv_ind_diseases_drugs_overlap = bc_inv_ind_diseases_drugs_overlap %>% filter(drug_name %in% candidate_drugs_repurposing$drug)
ggplot(bc_inv_ind_diseases_drugs_overlap, aes(y = max_phase, fill = max_phase)) +
  geom_bar(width = 0.4) +
  xlab("Number of drugs") +
  ylab("") +
  scale_x_continuous(breaks = seq(0,7, 1), limits = c(0,7)) +
  scale_fill_grey(start = 0, end = 0.8) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.y.right = element_text(size = 44, color = "black", family = "Arial", face = "bold", hjust = 0.5,
                                         margin = margin(t = 0, r = 0, b = 0, l = 40)),
        axis.text = element_text(size = 40, color = "black", family = "Arial"),
        axis.title = element_text(size = 40, color = "black", family = "Arial"),
        legend.position = "none")
ggsave("figure4B_candidate_drugs_inv_ind_per_phase.png", device = "png", path = "figures/", 
       width = 10, height = 6, dpi = 300)

rm(bc_inv_ind_drugs, depression_drug_targets, hdl_drug_targets, ldl_drug_targets, prostate_drug_targets, scz_drug_targets, t2dm_drug_targets, bc_depression_recommended_drugs, bc_hdl_recommended_drugs, bc_ldl_recommended_drugs, bc_prostate_recommended_drugs, bc_scz_recommended_drugs, bc_t2dm_recommended_drugs, candidate_drugs_repurposing, drugs_PAall, total_drugs_unique)
