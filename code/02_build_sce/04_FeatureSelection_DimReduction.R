#Goal: compile droplet scores, calculate QC metrics, and detect doublets. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/
#code modified from https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/03_reducedDimensions_clustering.R

library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(scater)
library(scran)
library(scry)
library(here)

## load QCed and cleaned object. 
#load(file = here("processed-data","sce_clean.rda"))
load(file = here("processed-data","sce_clean_numeric_cutoffs.rda"))

#sce
# class: SingleCellExperiment 
# dim: 36601 9883 
# metadata(1): Samples
# assays(1): counts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(6): source type ... gene_name gene_type
# colnames(9883): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 3_TTTGTTGAGGCTCCCA-1 3_TTTGTTGTCCCGATCT-1
# colData names(44): Sample Barcode ... discard_auto doubletScore
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):

sce
# class: SingleCellExperiment 
# dim: 36601 10000 
# metadata(1): Samples
# assays(1): counts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(6): source type ... gene_name gene_type
# colnames(10000): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 3_TTTGGTTTCTTCGACC-1 3_TTTGTTGTCCCGATCT-1
# colData names(54): Sample Barcode ... discared_numeric discard_numeric
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):
#Run deviance feature selection with default parameters. 
# sce <- devianceFeatureSelection(sce,
#                                 assay = "counts",
#                                 fam = "binomial",
#                                 sorted = FALSE,
#                                 batch = as.factor(sce$Sample))

# pdf(here("plots","featureSelxn_binomialDeviance-byGene.pdf"))
# plot(sort(rowData(sce)$binomial_deviance, decreasing = T),
#      type = "l", xlab = "ranked genes",
#      ylab = "binomial deviance"
# )
# abline(v = 2000,lty = 2, col = "red")
# dev.off()

sce <- devianceFeatureSelection(sce,
                                assay = "counts",
                                fam = "binomial",
                                sorted = FALSE,
                                batch = as.factor(sce$Sample))

pdf(here("plots","featureSelxn_binomialDeviance-byGene_numericcutoffs.pdf"))
plot(sort(rowData(sce)$binomial_deviance, decreasing = T),
     type = "l", xlab = "ranked genes",
     ylab = "binomial deviance"
)
abline(v = 2000,lty = 2, col = "red")
dev.off()


#2000 should be good. 
sce <- nullResiduals(sce,
                     assay = "counts", 
                     fam   = "binomial", 
                     type  = "deviance")
# In addition: Warning messages:
#     1: In asMethod(object) :
#     sparse->dense coercion: allocating vector of size 1.0 GiB
# 2: In asMethod(object) :
#     sparse->dense coercion: allocating vector of size 1.0 GiB
# 3: In .sparse2dense(x) :
#     sparse->dense coercion: allocating vector of size 1.0 GiB
# 4: In asMethod(object) :
#     sparse->dense coercion: allocating vector of size 1.0 GiB
# 5: In asMethod(object) :
#     sparse->dense coercion: allocating vector of size 1.0 GiB
# 6: In sqrt(x@x) : NaNs produced
#Not sure about warnings. sparse-->dense coercion doesn't seem to affect function. 
#No binomical deviances are NAs and  #ull residual matrix seems to be alright. 
#Need to figure out why NANs are being produced and why this could be happening? 
#In the meantime, all values are accounted for. Will move forward. 


#Take top 2000 highly deviant genes
hdgs <- rownames(sce)[order(rowData(sce)$binomial_deviance, decreasing = T)][1:2000]
hdgs.symbols <- rowData(sce)$gene_name[match(hdgs, rowData(sce)$gene_id)]

#Run PCA
sce_uncorrected <- runPCA(sce,
                          exprs_values = "binomial_deviance_residuals",
                          subset_row = hdgs, 
                          ncomponents = 100,
                          name = "GLMPCA_approx")

# UMAP
set.seed(1234)
sce_uncorrected <- runUMAP(sce_uncorrected,
                           dimred = "GLMPCA_approx",
                           n_dimred = 50, 
                           name = "UMAP")
# t-SNE
set.seed(1234)
sce_uncorrected <- runTSNE(sce_uncorrected,
                           dimred = "GLMPCA_approx",
                           n_dimred = 50, 
                           name = "TSNE")

#PCA plot of top 10 PCs
PCA_plots <- plotReducedDim(sce_uncorrected,
               dimred = "GLMPCA_approx", 
               colour_by = "Sample",
               ncomponents = 6, 
               point_alpha = 0.3)
#ggsave(PCA_plots,filename = here("plots","Dim_Red","multi_PCAs.png"))
ggsave(PCA_plots,filename = here("plots","Dim_Red","multi_PCAs_numericCutoff.png"))

# UMAPs
#Colored by sample, library size, and doublet score
#Sample
sample_umap <- plotReducedDim(sce_uncorrected,
                              dimred = "UMAP", 
                              colour_by = "Sample",
                              point_alpha = 0.3)
#ggsave(sample_umap,filename = here("plots","Dim_Red","sample_umap.png"))
ggsave(sample_umap,filename = here("plots","Dim_Red","sample_umap_numericCutoffs.png"))

#library size
sum_umap <- plotReducedDim(sce_uncorrected,
               dimred = "UMAP", colour_by = "sum",
               point_alpha = 0.3)
#ggsave(sum_umap,filename = here("plots","Dim_Red","sum_umap.png"))
ggsave(sum_umap,filename = here("plots","Dim_Red","sum_umap_numericCutoffs.png"))

#Doublet score
doublet_umap <- plotReducedDim(sce_uncorrected,
                               dimred = "UMAP", 
                               colour_by = "doubletScore")
#ggsave(doublet_umap,filename = here("plots","Dim_Red","doublet_umap.png"))
ggsave(doublet_umap,filename = here("plots","Dim_Red","doublet_umap_numericCutoffs.png"))

# TSNE by sample
sample_TSNE <- plotReducedDim(sce_uncorrected,
                              dimred = "TSNE", 
                              colour_by = "Sample")
#ggsave(sample_TSNE,filename = here("plots","Dim_Red","sample_TSNE.png"))
ggsave(sample_TSNE,filename = here("plots","Dim_Red","sample_TSNE_numericCutoffs.png"))


#batch effect results in cluster coming from a single sample. 
#Will need to run MNN to fix this. 
#save uncorrected object. 
#save(sce_uncorrected,file = here("processed-data","sce_uncorrected.rda"))
save(sce_uncorrected,file = here("processed-data","sce_uncorrected_numericCutoffs.rda"))
#Stopping point

#Load the uncorrected sce file that contains 
load(file = here("processed-data","sce_uncorrected.rda"))

sce_uncorrected
# class: SingleCellExperiment 
# dim: 36601 9883 
# metadata(1): Samples
# assays(2): counts binomial_deviance_residuals
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(9883): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 3_TTTGTTGAGGCTCCCA-1 3_TTTGTTGTCCCGATCT-1
# colData names(44): Sample Barcode ... discard_auto doubletScore
# reducedDimNames(3): GLMPCA_approx UMAP TSNE
# mainExpName: NULL
# altExpNames(0):

glmpca_mnn <- batchelor::reducedMNN(reducedDim(sce_uncorrected, "GLMPCA_approx"),
                                    batch=as.factor(sce_uncorrected$Sample))

#Add mnn to the object
reducedDim(sce_uncorrected,"mnn") <- glmpca_mnn$corrected

#Rename object
sce <- sce_uncorrected
rm(sce_uncorrected)

#Rerun umpa and tsne
set.seed(1234)

#umap
sce <- runUMAP(sce,
               dimred = "mnn",
               name = "UMAP_mnn")

# pdf(here("plots","UMAP_corrected_by_sample.pdf"))
# plotReducedDim(sce,
#                dimred = "UMAP_mnn", 
#                colour_by = "Sample")
# dev.off()
pdf(here("plots","Dim_Red","UMAP_corrected_by_sample_numericCutoffs.pdf"))
plotReducedDim(sce,
               dimred = "UMAP_mnn", 
               colour_by = "Sample")
dev.off()

# pdf(here("plots","Dim_Red","UMAP_corrected_mito.pdf"))
# plotReducedDim(sce,
#                dimred = "UMAP_mnn", colour_by = "subsets_Mito_percent",
#                point_alpha = 0.3)
# dev.off()
pdf(here("plots","Dim_Red","UMAP_corrected_mito_numericCutoffs.pdf"))
plotReducedDim(sce,
               dimred = "UMAP_mnn", colour_by = "subsets_Mito_percent",
               point_alpha = 0.3)
dev.off()

# pdf(here("plots","Dim_Red","UMAP_corrected_sum.pdf"))
# plotReducedDim(sce,
#                dimred = "UMAP_mnn", colour_by = "sum",
#                point_alpha = 0.3)
# dev.off()
pdf(here("plots","Dim_Red","UMAP_corrected_sum_numericCutoffs.pdf"))
plotReducedDim(sce,
               dimred = "UMAP_mnn", colour_by = "sum",
               point_alpha = 0.3)
dev.off()


pdf(here("plots","Dim_Red","UMAP_corrected_detected_numericCutoffs.pdf"))
plotReducedDim(sce,
               dimred = "UMAP_mnn", colour_by = "detected",
               point_alpha = 0.3)
dev.off()

#One cluster is dominated by Sample 1. Checking expression values to identify if this is a 
#batch correction issue, or if this is due to the 

#Comoute log counts to plot expression
sce <- batchelor::multiBatchNorm(sce, batch = sce$Sample)

genes <- c("SYT1","SNAP25", #pan neuron
           "MBP","MOBP", #OLIGODENDROCYTE
           "CD74", "CSF1R", "C3", #MICROGLIA
           "GFAP", "TNC", "AQP4", "SLC1A2", #ASTROCYTEs
           "GAD1","GAD2","SLC32A1",#Pan GABA
           "SLC17A7", "SLC17A6", "SLC17A8",
           "TRPC4","HOMER2","PTPN3", #Mouse LS markers
           "ELAVL2", #Mouse LS markers
           "CRHR1","CRHR2", 
           "OXTR","AVPR1A", 
           "DRD3")

for(i in genes){
    print(i)
    x <- plotReducedDim(sce,
                        dimred = "UMAP_mnn", 
                        colour_by = i,
                        swap_rownames = "gene_name") +
        scale_color_gradientn(colours = c("lightgrey","red")) +
        ggtitle(i) +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(filename = paste0("plots/Expression_plots/",i,"_expression_umap_numericCutoffs.pdf"),
           plot = x,
           height = 8,width = 8)
}



#Save the object
save(sce,file = here("processed-data","sce_numeric_cutoff_QC.rda"))

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# [1] "Reproducibility information:"
# [1] "2023-10-16 17:16:43 EDT"
#     user   system  elapsed 
# 1024.822   90.870 7503.115 
# ─ Session info ─────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 Patched (2023-07-19 r84711)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2023-10-16
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/pandoc
# 
# ─ Packages ─────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [2] CRAN (R 4.3.1)
# batchelor              1.16.0    2023-04-25 [2] Bioconductor
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
# cowplot                1.1.1     2020-12-30 [2] CRAN (R 4.3.1)
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.1)
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
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.1)
# ragg                   1.2.5     2023-01-12 [2] CRAN (R 4.3.1)
# Rcpp                   1.0.11    2023-07-06 [2] CRAN (R 4.3.1)
# RcppAnnoy              0.0.21    2023-07-02 [2] CRAN (R 4.3.1)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.1)
# ResidualMatrix         1.10.0    2023-04-25 [2] Bioconductor
# rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.1)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.3.1)
# Rtsne                  0.16      2022-04-17 [2] CRAN (R 4.3.1)
# S4Arrays               1.0.6     2023-08-30 [2] Bioconductor
# S4Vectors            * 0.38.1    2023-05-02 [2] Bioconductor
# ScaledMatrix           1.8.1     2023-05-03 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.3.1)
# scater               * 1.28.0    2023-04-25 [2] Bioconductor
# scran                * 1.28.2    2023-07-23 [2] Bioconductor
# scry                 * 1.12.0    2023-04-25 [2] Bioconductor
# scuttle              * 1.10.2    2023-08-03 [2] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.1)
# SingleCellExperiment * 1.22.0    2023-04-25 [2] Bioconductor
# sparseMatrixStats      1.12.2    2023-07-02 [2] Bioconductor
# statmod                1.5.0     2023-01-06 [2] CRAN (R 4.3.1)
# SummarizedExperiment * 1.30.2    2023-06-06 [2] Bioconductor
# systemfonts            1.0.4     2022-02-11 [2] CRAN (R 4.3.1)
# textshaping            0.3.6     2021-10-13 [2] CRAN (R 4.3.1)
# tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.3.1)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.3.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.3.1)
# uwot                   0.1.16    2023-06-29 [2] CRAN (R 4.3.1)
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
# ────────────────────────────────────────────────────────────────────────────────────
# 
# 
