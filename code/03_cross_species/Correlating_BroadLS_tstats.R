#Goal: Correlate standardized.logFC from all LS clusters. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(here)


##################
### Prep Human ###
##################

#load the SingleCellExperiment object for human Lateral Septum
load(here("processed-data","sce_human_sub.rda"),verbose = TRUE)
# Loading objects:
#   sce_human_sub

sce_human_sub
# class: SingleCellExperiment 
# dim: 16574 9225 
# metadata(1): Samples
# assays(3): counts binomial_pearson_residuals logcounts
# rownames(16574): ENSG00000187634 ENSG00000188976 ... ENSG00000205542
# ENSG00000092345
# rowData names(9): source type ... hs.entrezIds JAX.geneID
# colnames(9225): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 3_TTTGGTTTCTTCGACC-1 3_TTTGTTGTCCCGATCT-1
# colData names(60): Sample Barcode ... k_20_louvain_1 CellType.Final
# reducedDimNames(4): GLMPCA_approx tSNE_50 mnn tSNE_mnn_50
# mainExpName: NULL
# altExpNames(0):

identical(rownames(colData(sce_human_sub)),colnames(sce_human_sub))
# [1] TRUE

#Load the broad LS DEGs
load(here("processed-data","human_1vALL_broad.rda"),verbose = TRUE)
# Loading objects:
#   human_1vALL_broad

#subset the dataframe for LS
human_1vALL_broad <- subset(human_1vALL_broad, subset = (cellType.target == "LS"))

dim(human_1vALL_broad)
#[1] 36601     9

#subset DEGs for only homologs 
human_1vALL_homol <- merge(human_1vALL_broad,
                           as.data.frame(rowData(sce_human_sub)[,c("gene_id","JAX.geneID")]),
                           by = "gene_id") 

dim(human_1vALL_homol)
#[1] 16574    10

#Change rownames to JAX.geneID
rownames(human_1vALL_homol) <- human_1vALL_homol$JAX.geneID

##################
### Prep Mouse ###
##################

#Load mouse object 
#This object generated in 03_cross_species/mouse_LS_1vALL_DEGs.R after removing problematic clusters. 
load(here("processed-data","sce_mouse_sub.rda"),verbose = TRUE)
# Loading objects:
#   sce_mouse_sub

sce_mouse_sub
# dim: 16574 21884 
# metadata(1): Samples
# assays(3): counts binomial_pearson_residuals logcounts
# rownames(16574): ENSMUSG00000096351 ENSMUSG00000095567 ...
# ENSMUSG00000049775 ENSMUSG00000010592
# rowData names(9): source type ... mm.entrezIds JAX.geneID
# colnames(21884): 1_AAACCCAAGGTACATA-1 1_AAACCCACATCCGAGC-1 ...
# 4_TTTGTTGCATACAGCT-1 4_TTTGTTGGTCAAACGG-1
# colData names(17): Sample Barcode ... cellType.final cellType.broad
# reducedDimNames(4): GLMPCA_approx UMAP TSNE GLMPCA_50
# mainExpName: NULL
# altExpNames(0):

identical(rownames(colData(sce_mouse_sub)),colnames(sce_mouse_sub))
#[1] TRUE

##Load in DEGs
load(here("processed-data","mouse_1vALL_broad.rda"),verbose = TRUE)
# Loading objects:
#   mouse_1vALL_broad

#Subset for only LS clusters. 
mouse_1vALL_broad <- subset(mouse_1vALL_broad,subset=(cellType.target == "LS"))

dim(mouse_1vALL_broad)
#[1] 32285     9

#Add gene name information 
mouse_1vALL_homol <- merge(mouse_1vALL_broad,
                           as.data.frame(rowData(sce_mouse_sub)[,c("gene_id","JAX.geneID")]),
                           by = "gene_id") 

dim(mouse_1vALL_homol)
#[1] 16574    10

#Change rownames to JAX.geneID
rownames(mouse_1vALL_homol) <- mouse_1vALL_homol$JAX.geneID 

#Alter order of mouse DEGs to be that of human DEGs
mouse_1vALL_homol <- mouse_1vALL_homol[match(human_1vALL_homol$JAX.geneID,
                                             mouse_1vALL_homol$JAX.geneID),]

#Sanity check to make sure the rownames are in the corret order
identical(rownames(human_1vALL_homol),rownames(mouse_1vALL_homol))
#[1] TRUE

#Correlate the std.logFC 
cor(human_1vALL_homol$std.logFC,
    mouse_1vALL_homol$std.logFC)
#[1] 0.5954178

#Correlate the two dataframes to identify genes that are shared markers and those that are divergent markers
#To do this merge the two dataframes
colnames(human_1vALL_homol)[c(5,9)] <- paste0(colnames(human_1vALL_homol)[c(5,9)],"_human")
colnames(mouse_1vALL_homol)[c(5,9)] <- paste0(colnames(mouse_1vALL_homol)[c(5,9)],"_mouse")

#Merge the two dataframes
all_DEGs_homol <- merge(x = human_1vALL_homol[,c(5,9,10)],
                        y = mouse_1vALL_homol[,c(5,9,10)],
                        by = "JAX.geneID")

dim(all_DEGs_homol)
#[1] 16574     5

#Save the all_DEGs_homol dataframe
save(all_DEGs_homol,file = here("processed-data","Human_mouse_homologous_DEGs.rda"))

cor.test(all_DEGs_homol[,"std.logFC_human"],
         all_DEGs_homol[,"std.logFC_mouse"])
# Pearson's product-moment correlation
# 
# data:  all_DEGs_homol[, "std.logFC_human"] and all_DEGs_homol[, "std.logFC_mouse"]
# t = 95.404, df = 16572, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.5855009 0.6051566
# sample estimates:
#       cor 
# 0.5954178

#Create a dataframe subset for genes that we want to highlight. 
homol_subset <- subset(all_DEGs_homol,subset=(gene_name_human == "TRPC4" |
                                                gene_name_human == "FREM2" |
                                                gene_name_human == "MYO5B" |
                                                gene_name_human == "DGKG" |
                                                gene_name_human == "ANO1"))

#Make the plot
all_LS_plot <- ggplot(data=all_DEGs_homol,aes(x = std.logFC_human,y = std.logFC_mouse)) + 
  geom_point(alpha=0.5) +
  xlim(c(-3.5,3.5)) +
  ylim(c(-3.5,3.5)) +
  geom_hline(yintercept = 0,lty = 2) +
  geom_vline(xintercept = 0,lty = 2) +
  labs(x = "Human LS\nstandardized log-fold change", 
       y = "Mouse LS\nstandardized log-fold change") +
  theme_bw() +
  geom_text(data = homol_subset,
            label = homol_subset$gene_name_human,
            nudge_y = -.2,
            nudge_x = .2) +
  annotate("text",x = 3,y = 3.5,label = "Human Enriched\nMouse Enriched") +
  annotate("text",x = 3,y = -3.5, label = "Human Enriched\nMouse Depleted") +
  annotate("text",x = -3,y = 3.5,label = "Human Depleted\nMouse Enriched") +
  annotate("text",x = -3,y = -3.5,label = "Human Depleted\nMouse Depleted") 
ggsave(filename = here("plots","Conservation","all_LS_homologs_std.logFC_correlation.pdf"),plot = all_LS_plot)

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# [1] "Reproducibility information:"
# [1] "2024-09-06 14:59:37 EDT"
# user  system elapsed 
# 45.873   4.106 563.577
# ─ Session info ───────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.4.0 Patched (2024-05-22 r86590)
# os       Rocky Linux 9.4 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-09-06
# pandoc   3.1.13 @ /jhpce/shared/community/core/conda_R/4.4/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────
# package              * version date (UTC) lib source
# abind                  1.4-5   2016-07-21 [2] CRAN (R 4.4.0)
# Biobase              * 2.64.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# BiocGenerics         * 0.50.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# cli                    3.6.3   2024-06-21 [2] CRAN (R 4.4.0)
# colorspace             2.1-1   2024-07-26 [2] CRAN (R 4.4.0)
# crayon                 1.5.3   2024-06-20 [2] CRAN (R 4.4.0)
# DelayedArray           0.30.1  2024-05-07 [2] Bioconductor 3.19 (R 4.4.0)
# dplyr                  1.1.4   2023-11-17 [2] CRAN (R 4.4.0)
# fansi                  1.0.6   2023-12-08 [2] CRAN (R 4.4.0)
# farver                 2.1.2   2024-05-13 [2] CRAN (R 4.4.0)
# generics               0.1.3   2022-07-05 [2] CRAN (R 4.4.0)
# GenomeInfoDb         * 1.40.1  2024-05-24 [2] Bioconductor 3.19 (R 4.4.0)
# GenomeInfoDbData       1.2.12  2024-05-23 [2] Bioconductor
# GenomicRanges        * 1.56.1  2024-06-12 [2] Bioconductor 3.19 (R 4.4.0)
# ggplot2              * 3.5.1   2024-04-23 [2] CRAN (R 4.4.0)
# glue                   1.7.0   2024-01-09 [2] CRAN (R 4.4.0)
# gtable                 0.3.5   2024-04-22 [2] CRAN (R 4.4.0)
# here                 * 1.0.1   2020-12-13 [2] CRAN (R 4.4.0)
# httr                   1.4.7   2023-08-15 [2] CRAN (R 4.4.0)
# IRanges              * 2.38.1  2024-07-03 [2] Bioconductor 3.19 (R 4.4.0)
# jsonlite               1.8.8   2023-12-04 [2] CRAN (R 4.4.0)
# labeling               0.4.3   2023-08-29 [2] CRAN (R 4.4.0)
# lattice                0.22-6  2024-03-20 [3] CRAN (R 4.4.0)
# lifecycle              1.0.4   2023-11-07 [2] CRAN (R 4.4.0)
# magrittr               2.0.3   2022-03-30 [2] CRAN (R 4.4.0)
# Matrix                 1.7-0   2024-04-26 [3] CRAN (R 4.4.0)
# MatrixGenerics       * 1.16.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# matrixStats          * 1.3.0   2024-04-11 [2] CRAN (R 4.4.0)
# munsell                0.5.1   2024-04-01 [2] CRAN (R 4.4.0)
# pillar                 1.9.0   2023-03-22 [2] CRAN (R 4.4.0)
# pkgconfig              2.0.3   2019-09-22 [2] CRAN (R 4.4.0)
# R6                     2.5.1   2021-08-19 [2] CRAN (R 4.4.0)
# ragg                   1.3.2   2024-05-15 [2] CRAN (R 4.4.0)
# rlang                  1.1.4   2024-06-04 [2] CRAN (R 4.4.0)
# rprojroot              2.0.4   2023-11-05 [2] CRAN (R 4.4.0)
# S4Arrays               1.4.1   2024-05-20 [2] Bioconductor 3.19 (R 4.4.0)
# S4Vectors            * 0.42.1  2024-07-03 [2] Bioconductor 3.19 (R 4.4.0)
# scales                 1.3.0   2023-11-28 [2] CRAN (R 4.4.0)
# sessioninfo          * 1.2.2   2021-12-06 [2] CRAN (R 4.4.0)
# SingleCellExperiment * 1.26.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# SparseArray            1.4.8   2024-05-24 [2] Bioconductor 3.19 (R 4.4.0)
# SummarizedExperiment * 1.34.0  2024-05-01 [2] Bioconductor 3.19 (R 4.4.0)
# systemfonts            1.1.0   2024-05-15 [2] CRAN (R 4.4.0)
# textshaping            0.4.0   2024-05-24 [2] CRAN (R 4.4.0)
# tibble                 3.2.1   2023-03-20 [2] CRAN (R 4.4.0)
# tidyselect             1.2.1   2024-03-11 [2] CRAN (R 4.4.0)
# UCSC.utils             1.0.0   2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# utf8                   1.2.4   2023-10-22 [2] CRAN (R 4.4.0)
# vctrs                  0.6.5   2023-12-01 [2] CRAN (R 4.4.0)
# withr                  3.0.1   2024-07-31 [2] CRAN (R 4.4.0)
# XVector                0.44.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# zlibbioc               1.50.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# 
# [1] /users/rphillip/R/4.4
# [2] /jhpce/shared/community/core/conda_R/4.4/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.4/R/lib64/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────
