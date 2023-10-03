#Goal: Map expressed genes to human homologs. 
#Secondary goal: Make .gene.loc files for the expressed genes. 
library(here)
library(SingleCellExperiment)
library(scran)
library(scater)
library(sessioninfo)


#Load the medianNon0 list and the DEGs 
load(here("MAGMA",
          "mouse_analysis",
          "markers-stats_LS-n4_findMarkers_33cellTypes.rda"),
     verbose = TRUE)

#Can load express genes like this because genes with 0 counts were previously removed
expressed_genes <- rownames(markers.ls.t.1vAll[["Astro"]][["1"]])
length(unique(expressed_genes))
#[1] 27751

#Pull the two ensembl marts. 
#Sep 2019 is ensembl release 98, which was used for alinging mouse snRNA-seq data. 
hs_mart <- biomaRt::useMart("ensembl", 
                            dataset="hsapiens_gene_ensembl",
                            host = "https://sep2019.archive.ensembl.org")
mm_mart <- biomaRt::useMart("ensembl",
                            dataset="mmusculus_gene_ensembl",
                            host = "https://sep2019.archive.ensembl.org") 
#getLDS() pulls information from two linked datasets
#From the help page: In Ensembl this translates to homology mapping.
expressed_hom <- biomaRt::getLDS(attributes  = "ensembl_gene_id",
                                 mart        = mm_mart,
                                 values      = expressed_genes,
                                 filters     = "ensembl_gene_id",
                                 martL       = hs_mart,
                                 attributesL = c("ensembl_gene_id",
                                                 "external_gene_name",
                                                 "entrezgene_id"))
length(unique(expressed_hom$NCBI.gene.ID))
#[1] 17121

colnames(expressed_hom) <- c("Mouse_ensembl","Human_ensembl",
                             "Gene_Symbol","Entrez_gene_id")
#Ensembl release 98 contains GRCH38/Hg38 gene coordinates. These coordinates are 
#needed for the MDD and OUD GWAS. 
#Downloaded the gene loc file for hg19, now just need to merge the two dataframes. 
hg38_magma <- read.delim(file = here("MAGMA",
                                     "mouse_analysis",
                                     "NCBI38",
                                     "NCBI38.gene.loc"),
                         header = FALSE)
nrow(hg38_magma)
#[1] 20137
#20k annotations in hg38
#Subset hg38_magma for expressed genes. 
hg38_magma_homs <- hg38_magma[which(hg38_magma$V1 %in% expressed_hom$Entrez_gene_id),]
nrow(hg38_magma_homs)
#[1] 16954
#Moving forward with 16,954 genes that are expresseed (human homologs of mouse genes)

#Write out the file which will be the gene loc file for MAGMA analyses.
#Specifically those with Hg38 gene coordinates. 
write.table(x         = hg38_magma_homs,
            file      = here("MAGMA",
                             "mouse_analysis",
                             "mouse_expressing_hg38_homs.gene.loc"),
            col.names = FALSE,
            row.names = FALSE,
            sep       = "\t",
            quote     = FALSE)

#Now do the same thing for hg19
#Downloaded the gene loc file for hg19, now just need to merge the two dataframes. 
hg19_magma <- read.delim(file = here("MAGMA",
                                     "mouse_analysis",
                                     "NCBI37",
                                     "NCBI37.3.gene.loc"),
                         header = FALSE)

nrow(hg19_magma)
# [1] 19427

#Subset hg19_magma for the expressed genes
hg19_magma_homs <- hg19_magma[which(hg19_magma$V1 %in% expressed_hom$Entrez_gene_id),]
nrow(hg19_magma_homs)
#[1] 16892

#Write out the file which will be the gene loc file for MAGMA analyses.
write.table(x         = hg19_magma_homs,
            file      = here("MAGMA",
                             "mouse_analysis",
                             "mouse_expressing_hg19_homs.gene.loc"),
            col.names = FALSE,
            row.names = FALSE,
            sep       = "\t",
            quote     = FALSE)

Sys.time()
# [1] "2023-10-03 13:50:24 EDT"
options(width = 120)
session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 Patched (2023-07-19 r84711)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2023-10-03
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [2] CRAN (R 4.3.1)
# AnnotationDbi          1.62.2    2023-07-02 [2] Bioconductor
# beachmat               2.16.0    2023-04-25 [2] Bioconductor
# beeswarm               0.4.0     2021-06-01 [2] CRAN (R 4.3.1)
# Biobase              * 2.60.0    2023-04-25 [2] Bioconductor
# BiocFileCache          2.8.0     2023-04-25 [2] Bioconductor
# BiocGenerics         * 0.46.0    2023-04-25 [2] Bioconductor
# BiocNeighbors          1.18.0    2023-04-25 [2] Bioconductor
# BiocParallel           1.34.2    2023-05-22 [2] Bioconductor
# BiocSingular           1.16.0    2023-04-25 [2] Bioconductor
# biomaRt                2.56.1    2023-06-09 [2] Bioconductor
# Biostrings             2.68.1    2023-05-16 [2] Bioconductor
# bit                    4.0.5     2022-11-15 [2] CRAN (R 4.3.1)
# bit64                  4.0.5     2020-08-30 [2] CRAN (R 4.3.1)
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.1)
# blob                   1.2.4     2023-03-17 [2] CRAN (R 4.3.1)
# bluster                1.10.0    2023-04-25 [2] Bioconductor
# cachem                 1.0.8     2023-05-01 [2] CRAN (R 4.3.1)
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.1)
# cluster                2.1.4     2022-08-22 [3] CRAN (R 4.3.1)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.3.1)
# colorout             * 1.2-2     2023-09-22 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.3.1)
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.1)
# curl                   5.0.2     2023-08-14 [2] CRAN (R 4.3.1)
# DBI                    1.1.3     2022-06-18 [2] CRAN (R 4.3.1)
# dbplyr                 2.3.3     2023-07-07 [2] CRAN (R 4.3.1)
# DelayedArray           0.26.7    2023-07-28 [2] Bioconductor
# DelayedMatrixStats     1.22.6    2023-08-28 [2] Bioconductor
# digest                 0.6.33    2023-07-07 [2] CRAN (R 4.3.1)
# dplyr                  1.1.3     2023-09-03 [2] CRAN (R 4.3.1)
# dqrng                  0.3.1     2023-08-30 [2] CRAN (R 4.3.1)
# edgeR                  3.42.4    2023-05-31 [2] Bioconductor
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.3.1)
# fastmap                1.1.1     2023-02-24 [2] CRAN (R 4.3.1)
# filelock               1.0.2     2018-10-05 [2] CRAN (R 4.3.1)
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
# hms                    1.1.3     2023-03-21 [2] CRAN (R 4.3.1)
# httr                   1.4.7     2023-08-15 [2] CRAN (R 4.3.1)
# igraph                 1.5.1     2023-08-10 [2] CRAN (R 4.3.1)
# IRanges              * 2.34.1    2023-06-22 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.3.1)
# KEGGREST               1.40.0    2023-04-25 [2] Bioconductor
# lattice                0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.3.1)
# limma                  3.56.2    2023-06-04 [2] Bioconductor
# locfit                 1.5-9.8   2023-06-11 [2] CRAN (R 4.3.1)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.3.1)
# Matrix                 1.6-1.1   2023-09-18 [3] CRAN (R 4.3.1)
# MatrixGenerics       * 1.12.3    2023-07-30 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.3.1)
# memoise                2.0.1     2021-11-26 [2] CRAN (R 4.3.1)
# metapod                1.8.0     2023-04-25 [2] Bioconductor
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.3.1)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.3.1)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.3.1)
# png                    0.1-8     2022-11-29 [2] CRAN (R 4.3.1)
# prettyunits            1.1.1     2020-01-24 [2] CRAN (R 4.3.1)
# progress               1.2.2     2019-05-16 [2] CRAN (R 4.3.1)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.1)
# rappdirs               0.3.3     2021-01-31 [2] CRAN (R 4.3.1)
# Rcpp                   1.0.11    2023-07-06 [2] CRAN (R 4.3.1)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.1)
# rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.1)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
# RSQLite                2.3.1     2023-04-03 [2] CRAN (R 4.3.1)
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
# tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.3.1)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.3.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.3.1)
# vctrs                  0.6.3     2023-06-14 [2] CRAN (R 4.3.1)
# vipor                  0.4.5     2017-03-22 [2] CRAN (R 4.3.1)
# viridis                0.6.4     2023-07-22 [2] CRAN (R 4.3.1)
# viridisLite            0.4.2     2023-05-02 [2] CRAN (R 4.3.1)
# withr                  2.5.0     2022-03-03 [2] CRAN (R 4.3.1)
# XML                    3.99-0.14 2023-03-19 [2] CRAN (R 4.3.1)
# xml2                   1.3.5     2023-07-06 [2] CRAN (R 4.3.1)
# XVector                0.40.0    2023-04-25 [2] Bioconductor
# zlibbioc               1.46.0    2023-04-25 [2] Bioconductor
# 
# [1] /users/rphillip/R/4.3
# [2] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────