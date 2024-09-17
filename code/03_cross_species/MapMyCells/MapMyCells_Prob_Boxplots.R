#Make a boxplot for class bootstrapping probability for a response to reviewer #3 
library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(here)

#load in the sce object
load(here("processed-data","02_build_sce","sce_celltype.rda"),verbose = TRUE)
# Loading objects:
#   sce

#load the cluster colors
load(here("processed-data","cluster_cols_CellType_Final.rda"),verbose = TRUE)
# Loading objects:
#   cluster_cols


#Read in the map my cells output
Map_Cells_Out_hi <- read.delim(file = here("processed-data",
                                           "MapMyCells_Output",
                                           "h_ls_anndata_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1726506607196",
                                           "h_ls_anndata_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1726506607196.csv"),
                               comment.char = "#",
                               sep = ",")

#Merge with celltype final information from the sce object
hi_output <- merge(x = colData(sce)[,c("key","CellType.Final")],
                   y = Map_Cells_Out_hi,
                   by.x = "key",
                   by.y = "cell_id")


#Boxplot with x axis being the class that each human cell is assigned. 
mouse_class <- ggplot(data = hi_output,aes(x = class_name,
                            y = class_bootstrapping_probability,
                            fill = class_name)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "none")

ggsave(filename = here("plots","Conservation","MapMyCells","Mouse_Class_Bootstrap_Prob_Boxplot.png"),
       plot = mouse_class)

#Boxplot with x axis being the class that each human cell is assigned. 
human_class <- ggplot(data = hi_output,aes(x = CellType.Final, y = class_bootstrapping_probability,fill = CellType.Final)) +
  scale_fill_manual(values = cluster_cols) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "none")

ggsave(filename = here("plots","Conservation","MapMyCells","Human_Class_Bootstrap_Prob_Boxplot.png"),
       plot = human_class)

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# [1] "Reproducibility information:"
# [1] "2024-09-17 17:33:30 EDT"
# user  system elapsed 
# 32.350   2.338  90.532
# ─ Session info ────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.4.0 Patched (2024-05-22 r86590)
# os       Rocky Linux 9.4 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-09-17
# pandoc   3.1.13 @ /jhpce/shared/community/core/conda_R/4.4/bin/pandoc
# 
# ─ Packages ────────────────────────────────────────────────────────────────
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
# ───────────────────────────────────────────────────────────────────────────
# 
# 
# 
