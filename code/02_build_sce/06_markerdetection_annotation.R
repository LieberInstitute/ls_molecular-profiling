#Goal: Identify marker genes for each cluster and annotate the clusters. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

library(SingleCellExperiment)
library(DeconvoBuddies)
library(sessioninfo)
library(ggplot2)
library(scater)
library(scran)
library(here)

#load the clustered object
load(here("processed-data","sce_clustered.rda"))

sce
# class: SingleCellExperiment 
# dim: 36601 9225 
# metadata(1): Samples
# assays(3): counts binomial_deviance_residuals logcounts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(9225): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 3_TTTGGTTTCTTCGACC-1 3_TTTGTTGTCCCGATCT-1
# colData names(59): Sample Barcode ... k_50_walktrap k_75_walktrap
# reducedDimNames(14): GLMPCA_approx UMAP_15 ... UMAP_mnn_25 UMAP_mnn_50
# mainExpName: NULL
# altExpNames(0):

##################################################
############Annotate the clusters.################
##################################################
annotation_df <- data.frame(cluster=c(1:15),
                            celltype = c("Astrocyte_1","Astrocyte_2","LS_GABA_2","Glutamatergic","MS_GABA_1",
                                         "Striosome","Polydendrocyte","LS_GABA_1","Oligo_1","Astrocyte_3",
                                         "Microglia","Oligo_2","Oligo_3","MS_GABA_2","Endothelial"))

sce$CellType <- annotation_df$celltype[match(sce$k_50_walktrap,
                                             annotation_df$cluster)]
sce$CellType <- factor(sce$CellType,
                       levels = c("LS_GABA_1","LS_GABA_2",
                                  "MS_GABA_1","MS_GABA_2",
                                  "Glutamatergic","Striosome",
                                  "Astrocyte_1","Astrocyte_2","Astrocyte_3",
                                  "Oligo_1","Oligo_2","Oligo_3",
                                  "Microglia",
                                  "Polydendrocyte",
                                  "Endothelial"))

###check expression profiles of the clusters
genes <- c("SYT1","SNAP25","GAD1","GAD2","SLC32A1",
           "TRPC4","HOMER2","PTPN3",
           "ELAVL2",
           "SLC17A7", "SLC17A6", "SLC17A8",
           "DRD1","OPRM1","FOXP2",
           "GFAP", "TNC", "AQP4", "SLC1A2",
           "MBP","MOBP",
           "CD74", "CSF1R", "C3",
           "PDGFRA", "VCAN", "CSPG4",
           "CLDN5", "FLT1", "VTN")


Expression_dotplot <- plotDots(object = sce,
                               features = rev(genes),
                               group = "CellType",swap_rownames = "gene_name") +
    scale_color_gradientn(colours = c("lightgrey","orange","red")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = Expression_dotplot,filename = here("plots","Expression_plots",
                                                 "annotated_dotplot.pdf"),
       height = 8)

###Annotate the umap. 
annotated_umap <- plotReducedDim(object = sce,dimred = "UMAP_mnn_15",colour_by = "CellType",text_by = "CellType")
ggsave(plot = annotated_umap,filename = here("plots","Dim_Red",
                                             "annotated_umap.pdf"))

#Another dotplot with other canonical LS genes. 
genes <- c("SYT1","SNAP25","GAD1","GAD2","SLC32A1",
           "TRPC4","HOMER2","PTPN3",
           "ELAVL2",
           "SLC17A7", "SLC17A6", "SLC17A8",
           "DRD1","OPRM1","FOXP2",
           "NTS","CRHR1","CRHR2","OXTR","AVPR1A","DRD3",
           "GFAP", "TNC", "AQP4", "SLC1A2",
           "MBP","MOBP",
           "CD74", "CSF1R", "C3",
           "PDGFRA", "VCAN", "CSPG4",
           "CLDN5", "FLT1", "VTN")


Expression_dotplot <- plotDots(object = sce,
                               features = rev(genes),
                               group = "CellType",swap_rownames = "gene_name") +
    scale_color_gradientn(colours = c("lightgrey","orange","red")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = Expression_dotplot,filename = here("plots","Expression_plots",
                                                 "annotated_dotplot_withadditionalLSmarkers.pdf"),
       height = 8)

violins <- plotExpression(object = sce,
                          features = c("TRPC4","HOMER2","PTPN3",
                                       "NTS","CRHR1","CRHR2",
                                       "OXTR","AVPR1A","DRD3"),
                          swap_rownames = "gene_name",x = "CellType",ncol = 3,
                          colour_by = "CellType") +
    theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")
ggsave(plot = violins,filename = here("plots","Expression_plots",
                                                 "annotated_violin_LSmarkers.pdf"),
       width = 12,
       height = 12)


##################################################
#########Prep sce object for DEG testing##########
##################################################

#Do any genes have 0 counts for every cell. 
table(rowSums(assay(sce, "counts")) == 0)
# FALSE  TRUE 
# 33556  3045 

#Remove 3045 genes that are all 0s 
sce <- sce[!rowSums(assay(sce, "counts")) == 0, ]

##################################################
###############run 1 vs all testing###############
##################################################
markers_1vALL_enrich <- findMarkers_1vAll(sce, 
                                          assay_name = "logcounts", 
                                          cellType_col = "CellType", 
                                          mod = "~Sample")
# LS_GABA_1 - '2023-10-24 16:13:20.146492
# LS_GABA_2 - '2023-10-24 16:13:36.971206
# Glutamatergic - '2023-10-24 16:13:54.641428
# Polydendrocyte - '2023-10-24 16:14:12.438615
# Microglia - '2023-10-24 16:14:30.290416
# Oligo_2 - '2023-10-24 16:14:48.158806
# Astrocyte_1 - '2023-10-24 16:15:06.106722
# Oligo_1 - '2023-10-24 16:15:24.129632
# Striosome - '2023-10-24 16:15:42.010409
# MS_GABA_2 - '2023-10-24 16:15:59.972938
# MS_GABA_1 - '2023-10-24 16:16:17.891058
# Endothelial - '2023-10-24 16:16:35.84846
# Astrocyte_2 - '2023-10-24 16:16:53.755235
# Astrocyte_3 - '2023-10-24 16:17:11.606782
# Oligo_3 - '2023-10-24 16:17:30.754411
# Building Table - 2023-10-24 16:17:48.676136
# ** Done! **

#Add symbol information to the table
#First change the ensembl gene id column to have same name as what is in rowData(sce)
colnames(markers_1vALL_enrich)[1] <- "gene_id"
markers_1vALL_df <- dplyr::left_join(x = as.data.frame(markers_1vALL_enrich),
                                     y = as.data.frame(rowData(sce)[,c("gene_id","gene_name")]),
                                     by = "gene_id")

#save the dataframe. 
save(markers_1vALL_df,file = here("processed-data","markers_1vAll_ttest.rda"))


##################################################
###############run pairwise testing###############
##################################################
mod <- with(colData(sce), model.matrix(~ Sample))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

# Run pairwise t-tests
markers_pairwise <- findMarkers(sce, 
                                groups=sce$CellType,
                                assay.type="logcounts", 
                                design=mod, 
                                test="t",
                                direction="up", 
                                pval.type="all", 
                                full.stats=T)

#How many DEGs for each cluster? 
sapply(markers_pairwise, function(x){table(x$FDR<0.05)})
#       LS_GABA_1 LS_GABA_2 MS_GABA_1 MS_GABA_2 Glutamatergic Striosome
# FALSE     33268     32954     33408     32704         32700     33311
# TRUE        288       602       148       852           856       245
#        Astrocyte_1 Astrocyte_2 Astrocyte_3 Oligo_1 Oligo_2 Oligo_3 Microglia
# FALSE       31437       33483       32991   33555   33484   32492     32410
# TRUE         2119          73         565       1      72    1064      1146
#       Polydendrocyte Endothelial
# FALSE          33056       32349
# TRUE             500        1207

#Add gene info to each list. 
for(i in names(markers_pairwise)){
    markers_pairwise[[i]] <- as.data.frame(markers_pairwise[[i]])
    markers_pairwise[[i]]$gene_id <- row.names(markers_pairwise[[i]])
    markers_pairwise[[i]] <- dplyr::left_join(x  =  markers_pairwise[[i]],
                                              y  =  as.data.frame(rowData(sce)[,c("gene_id","gene_name")]),
                                              by = "gene_id")
}

save(markers_pairwise,file = here("processed-data","markers_pairwise_list.rda"))

#Save the object with celltype
save(sce,file = here("processed-data","sce_with_CellType.rda"))

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# [1] "Reproducibility information:"
# [1] "2023-10-24 16:36:28 EDT"
#    user   system  elapsed 
# 710.000   10.332 1536.325 
# ─ Session info ────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 Patched (2023-07-19 r84711)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2023-10-24
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/pandoc
# 
# ─ Packages ────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [2] CRAN (R 4.3.1)
# beachmat               2.16.0    2023-04-25 [2] Bioconductor
# beeswarm               0.4.0     2021-06-01 [2] CRAN (R 4.3.1)
# Biobase              * 2.60.0    2023-04-25 [2] Bioconductor
# BiocGenerics         * 0.46.0    2023-04-25 [2] Bioconductor
# BiocNeighbors          1.18.0    2023-04-25 [2] Bioconductor
# BiocParallel           1.34.2    2023-05-22 [2] Bioconductor
# BiocSingular           1.16.0    2023-04-25 [2] Bioconductor
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.1)
# bluster                1.10.0    2023-04-25 [2] Bioconductor
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.1)
# cluster                2.1.4     2022-08-22 [3] CRAN (R 4.3.1)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.3.1)
# colorout             * 1.2-2     2023-09-22 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.3.1)
# cowplot                1.1.1     2020-12-30 [1] CRAN (R 4.3.1)
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.1)
# DeconvoBuddies       * 0.99.0    2023-10-23 [1] Github (LieberInstitute/DeconvoBuddies@9ce4a42)
# DelayedArray           0.26.7    2023-07-28 [2] Bioconductor
# DelayedMatrixStats     1.22.6    2023-08-28 [2] Bioconductor
# dplyr                  1.1.3     2023-09-03 [2] CRAN (R 4.3.1)
# dqrng                  0.3.1     2023-08-30 [2] CRAN (R 4.3.1)
# edgeR                  3.42.4    2023-05-31 [2] Bioconductor
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.3.1)
# farver                 2.1.1     2022-07-06 [2] CRAN (R 4.3.1)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.3.1)
# GenomeInfoDb         * 1.36.3    2023-09-07 [2] Bioconductor
# GenomeInfoDbData       1.2.10    2023-07-20 [2] Bioconductor
# GenomicRanges        * 1.52.0    2023-04-25 [2] Bioconductor
# ggbeeswarm             0.7.2     2023-04-29 [2] CRAN (R 4.3.1)
# ggplot2              * 3.4.3     2023-08-14 [2] CRAN (R 4.3.1)
# ggrepel                0.9.3     2023-02-03 [2] CRAN (R 4.3.1)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.3.1)
# gridExtra              2.3       2017-09-09 [2] CRAN (R 4.3.1)
# gtable                 0.3.4     2023-08-21 [2] CRAN (R 4.3.1)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.3.1)
# igraph                 1.5.1     2023-08-10 [2] CRAN (R 4.3.1)
# IRanges              * 2.34.1    2023-06-22 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.3.1)
# labeling               0.4.3     2023-08-29 [2] CRAN (R 4.3.1)
# lattice                0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.3.1)
# limma                  3.56.2    2023-06-04 [2] Bioconductor
# locfit                 1.5-9.8   2023-06-11 [2] CRAN (R 4.3.1)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.3.1)
# Matrix                 1.6-1.1   2023-09-18 [3] CRAN (R 4.3.1)
# MatrixGenerics       * 1.12.3    2023-07-30 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.3.1)
# metapod                1.8.0     2023-04-25 [2] Bioconductor
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.3.1)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.3.1)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.3.1)
# purrr                  1.0.2     2023-08-10 [2] CRAN (R 4.3.1)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.1)
# rafalib                1.0.0     2015-08-09 [1] CRAN (R 4.3.1)
# ragg                   1.2.5     2023-01-12 [2] CRAN (R 4.3.1)
# RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.3.1)
# Rcpp                   1.0.11    2023-07-06 [2] CRAN (R 4.3.1)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.1)
# rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.1)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.3.1)
# S4Arrays               1.0.6     2023-08-30 [2] Bioconductor
# S4Vectors            * 0.38.1    2023-05-02 [2] Bioconductor
# ScaledMatrix           1.8.1     2023-05-03 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.3.1)
# scater               * 1.28.0    2023-04-25 [2] Bioconductor
# scran                * 1.28.2    2023-07-23 [2] Bioconductor
# scuttle              * 1.10.2    2023-08-03 [2] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.1)
# SingleCellExperiment * 1.22.0    2023-04-25 [2] Bioconductor
# sparseMatrixStats      1.12.2    2023-07-02 [2] Bioconductor
# statmod                1.5.0     2023-01-06 [2] CRAN (R 4.3.1)
# stringi                1.7.12    2023-01-11 [2] CRAN (R 4.3.1)
# stringr                1.5.0     2022-12-02 [2] CRAN (R 4.3.1)
# SummarizedExperiment * 1.30.2    2023-06-06 [2] Bioconductor
# systemfonts            1.0.4     2022-02-11 [2] CRAN (R 4.3.1)
# textshaping            0.3.6     2021-10-13 [2] CRAN (R 4.3.1)
# tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.3.1)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.3.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.3.1)
# vctrs                  0.6.3     2023-06-14 [2] CRAN (R 4.3.1)
# vipor                  0.4.5     2017-03-22 [2] CRAN (R 4.3.1)
# viridis                0.6.4     2023-07-22 [2] CRAN (R 4.3.1)
# viridisLite            0.4.2     2023-05-02 [2] CRAN (R 4.3.1)
# withr                  2.5.0     2022-03-03 [2] CRAN (R 4.3.1)
# XVector                0.40.0    2023-04-25 [2] Bioconductor
# zlibbioc               1.46.0    2023-04-25 [2] Bioconductor
# 
# [1] /users/rphillip/R/4.3
# [2] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library
# 
# ───────────────────────────────────────────────────────────────────────────────────────────
