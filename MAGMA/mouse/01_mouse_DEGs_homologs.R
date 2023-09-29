#Goal: Map mouse cell type specific DEGs to human. 
library(here)
library(SingleCellExperiment)
library(scran)
library(scater)
library(sessioninfo)

#Load in the DEG list
#Code to calculate DEGs can be found at https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/04_markerDetection.R
load(here("mouse",
          "markers-stats_LS-n4_findMarkers_33cellTypes.rda"),
     verbose = TRUE)

#DEGs generated with the cluster-vs-all method are in the markers.ls.t.1vAll object
#Add the medianNon0 information into the cluster. Code taken from the same github link. 
#Each iteration of the list contains two lists.
#Want the list "1" which includes genes that are enriched in the cluster. 
for(i in names(markers.ls.t.1vAll)){
    print(i)
    markers.ls.t.1vAll[[i]][["1"]] <- cbind(
        markers.ls.t.1vAll[[i]][["1"]],
        medianNon0.ls[[i]][
            match(row.names(markers.ls.t.1vAll[[i]][["1"]]),
                  names(medianNon0.ls[[i]]))
        ]
    )
    colnames(markers.ls.t.1vAll[[i]][["1"]])[5] <- "non0Median"
    markers.ls.t.1vAll[[i]][["1"]]$Ensembl_gene_id <- row.names(markers.ls.t.1vAll[[i]][["1"]])
}

#FDR<1e-2 and median > 0
markers <- lapply(markers.ls.t.1vAll,function(x){
    as.data.frame(subset(x[["1"]],
                         subset=(summary.stats > 0 & FDR < 1e-2 & non0Median == TRUE)))
})

#How many genes make up each set? 
as.matrix(lapply(markers,FUN = nrow))
# [,1]
# Astro         486 
# Chol_Ex.D     608 
# ChP           1232
# Endo          854 
# Ependymal     1400
# IoC_In.E      569 
# LS_In.C       1032
# LS_In.D       1633
# LS_In.M       990 
# LS_In.N       368 
# LS_In.O       515 
# LS_In.P       1386
# LS_In.Q       435 
# LS_In.R       1030
# Micro         437 
# MS_In.J       965 
# MS_In.K       1405
# Mural         239 
# Oligo         385 
# OPC           583 
# OPC_COP       527 
# Sept_In.G     720 
# Sept_In.I     1495
# Str_In.A      980 
# Str_In.F      1576
# Str_In.H      604 
# Str_In.L      400 
# Thal_Ex.B     1779
# TNoS_Ex.A     773 
# TT.IG.SH_Ex.C 1173
# TT.IG.SH_Ex.E 466 
# TT.IG.SH_Ex.F 1615
# Ventr_In.B    585 

##Using sep2019 ensembl (v98) because that was used to align the data. 
hs_mart <- biomaRt::useMart("ensembl", 
                            dataset="hsapiens_gene_ensembl",
                            host = "https://sep2019.archive.ensembl.org")
mm_mart <- biomaRt::useMart("ensembl",
                            dataset="mmusculus_gene_ensembl",
                            host = "https://sep2019.archive.ensembl.org") 

markers_hom <- lapply(markers,function(x){
    biomaRt::getLDS(attributes  = "ensembl_gene_id",
                    mart        = mm_mart,
                    values      = x$Ensembl_gene_id,
                    filters     = "ensembl_gene_id",
                    martL       = hs_mart,
                    attributesL = c("ensembl_gene_id",
                                    "external_gene_name",
                                    "entrezgene_id"))
})

##Need the entrez gene id because I will be using the gene loc files
#provided by MAGMA. 

#Rearrange the dataframes. 
for(i in names(markers_hom)){
    colnames(markers_hom[[i]]) <-  c("Mouse_ensembl","Hg38_Ensembl",
                                     "gene_symbol","entrez_id")
}

#Save this entire list. 
save(markers_hom,
     file = here("mouse","mouse_markers_human_homologs_list.rda"))

#The markers_hom list is set up so that column 1 is the mouse ensembl and column 2 is human ensembl
#For magma we just need to go ahead and have column 1 be the set or cell type.
#Column 2 needs to be human ensembl gene name. 
for(i in names(markers_hom)){
    markers_hom[[i]]$Set <- i
    markers_hom[[i]] <- markers_hom[[i]][,c("Set","entrez_id")]
    colnames(markers_hom[[i]])[2] <- "Gene"
}

#collapse the list into a single file. 
markers_hom_df <- do.call(what = rbind,markers_hom)

table(is.na(markers_hom_df$Gene))
# FALSE  TRUE 
# 27173   288 
#288 genes don't have matching entrez gene ids. 
markers_hom_df <- markers_hom_df[!is.na(markers_hom_df$Gene),]

#write out the table. 
write.table(x = markers_hom_df,
            file = here("mouse",
                        "mouse_markers_human_homologs_1e-2.txt"),
            col.names = TRUE,
            row.names = FALSE,
            quote     = FALSE,
            sep       = "\t")

Sys.time()
#[1] "2023-09-29 13:42:04 EDT"
options(width = 120)
session_info()
# ─ Session info ─────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 Patched (2023-07-19 r84711)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2023-09-29
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/pandoc
# 
# ─ Packages ─────────────────────────────────────────────────────────────────
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
# ────────────────────────────────────────────────────────────────────────────
# 
# 
# 
# 
# 
# 
# 
