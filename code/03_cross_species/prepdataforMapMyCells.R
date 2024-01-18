#Goal: Build anndata file for mapmycells
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling
library(SingleCellExperiment)
library(zellkonverter)
library(sessioninfo)
library(anndata)
library(here)

#Load human LS sce object
load(here("processed-data","sce_with_CellType.rda"),verbose = TRUE)
# Loading objects:
#     sce

######Python environment setup######
#Set up python environment with the reticulate package 
library(reticulate)

# install anndata
virtualenv_install("r-reticulate", "anndata")

# import anndata (it will be automatically discovered in "r-reticulate")
anndata <- import("anndata")
############################

######Convert gene ids to mouse#########
#Pull counts matrix from sce
#Map my cells correlation method expects raw counts. 
counts <- counts(sce)
dim(counts)
#[1] 33556  9225

# Read in a conversion table generated from GeneOrthology
convert <- read.csv("https://github.com/AllenInstitute/GeneOrthology/raw/main/csv/mouse_human_marmoset_macaque_orthologs_20231113.csv")
dim(convert)
#[1] 18381    14

# Remove NAs from conversion table.
#Only keep mouse and human because we will map human to mouse
convert_by_ensembl <- convert[!(is.na(convert$human_EnsemblID)|is.na(convert$mouse_EnsemblID)),
                              c("human_EnsemblID","mouse_EnsemblID")]
dim(convert_by_ensembl)
#[1] 16529     2

# Remove genes not in data matrix
convert_by_ensembl <- convert_by_ensembl[is.element(convert_by_ensembl$human_EnsemblID,rownames(counts)),] 
dim(convert_by_ensembl)
#[1] 16193     2

#There is a row that is duplicated within the dataframe. Remove it. 
convert_by_ensembl <- convert_by_ensembl[!duplicated(convert_by_ensembl),]
dim(convert_by_ensembl)
#[1] 16192     2

# Subset data to include only genes with mouse othologs
new_order <- match(convert_by_ensembl$human_EnsemblID,rownames(counts))
new_order <- new_order[!is.na(new_order)]

counts_out <- counts[new_order,]
dim(counts_out)
#[1] 16192  9225

#Double check that everything is in the same order and no duplicates
all(rownames(counts_out) == convert_by_ensembl$human_EnsemblID)
#[1] TRUE

table(duplicated(rownames(counts_out)))
# FALSE 
# 16192 

#Convert the rownames to the mouse ensembl IDs
rownames(counts_out) <- convert_by_ensembl$mouse_EnsemblID

##Some human genes map to the same mouse ensembl id.
#Check to see if there are any duplicated mouse genes within the convert_by_ensembl table
dups <- convert_by_ensembl$mouse_EnsemblID[duplicated(convert_by_ensembl$mouse_EnsemblID)]
dups
# [1] "ENSMUSG00000087408" "ENSMUSG00000029723" "ENSMUSG00000023156"
# [4] "ENSMUSG00000056629"

#There are 4 genes that are duplicated. Loop through those genes to find the row that is less expressed. 
rows_to_remove <- vector(length = 4)
for(i in 1:4){
    dup_rows <- grep(dups[i],rownames(counts_out))
    rows_to_remove[i] <- dup_rows[which(rowSums(counts_out[dup_rows,]) == min(rowSums(counts_out[dup_rows,])))]
}

#Remove those rows. 
counts_out <- counts_out[-rows_to_remove,]

#Check that all duplicates are gone. 
table(duplicated(rownames(counts_out)))
# FALSE 
# 16188

#Transpose because mapmycells needs genes in column and cells in rows. 
counts_out <- t(counts_out)

# Convert to anndata format
h_LS_anndata <- AnnData(X = counts_out,
                        obs = data.frame(group = rownames(counts_out),
                                         row.names = rownames(counts_out)),
                        var = data.frame(group = colnames(counts_out),
                                         row.names = colnames(counts_out)))


object.size(h_LS_anndata)
#352 bytes

#Write out the file. 
write_h5ad(h_LS_anndata,
           here("processed-data","h_ls_anndata.h5ad"),
           compression="gzip")

# Check file size. File MUST be <500MB to upload for MapMyCells
print(paste("Size in MB:",
            round(file.size(here("processed-data",
                                 "h_ls_anndata.h5ad"))/2^20)))

#[1] "Size in MB: 88"
#File is of the right size and is ready to go. 
#Used 10x genomics whole mouse brain and the both mapping algorithms.


#Reproducibility information. 
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# [1] "Reproducibility information:"
# [1] "2024-01-18 13:36:35 EST"
# user  system elapsed 
# 71.362   4.447 767.001 
# ─ Session info ───────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 Patched (2023-07-19 r84711)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-01-18
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [2] CRAN (R 4.3.1)
# anndata              * 0.7.5.6   2023-03-17 [1] CRAN (R 4.3.1)
# assertthat             0.2.1     2019-03-21 [2] CRAN (R 4.3.1)
# basilisk               1.12.1    2023-06-30 [2] Bioconductor
# basilisk.utils         1.12.1    2023-05-19 [2] Bioconductor
# Biobase              * 2.60.0    2023-04-25 [2] Bioconductor
# BiocGenerics         * 0.46.0    2023-04-25 [2] Bioconductor
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.1)
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.1)
# colorout             * 1.3-0.1   2023-12-01 [1] Github (jalvesaq/colorout@deda341)
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.1)
# DelayedArray           0.26.7    2023-07-28 [2] Bioconductor
# dir.expiry             1.8.0     2023-04-25 [2] Bioconductor
# filelock               1.0.2     2018-10-05 [2] CRAN (R 4.3.1)
# GenomeInfoDb         * 1.36.3    2023-09-07 [2] Bioconductor
# GenomeInfoDbData       1.2.10    2023-07-20 [2] Bioconductor
# GenomicRanges        * 1.52.0    2023-04-25 [2] Bioconductor
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.3.1)
# IRanges              * 2.34.1    2023-06-22 [2] Bioconductor
# jsonlite               1.8.7     2023-06-29 [2] CRAN (R 4.3.1)
# lattice                0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
# Matrix                 1.6-1.1   2023-09-18 [3] CRAN (R 4.3.1)
# MatrixGenerics       * 1.12.3    2023-07-30 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.3.1)
# png                    0.1-8     2022-11-29 [2] CRAN (R 4.3.1)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.1)
# rappdirs               0.3.3     2021-01-31 [2] CRAN (R 4.3.1)
# Rcpp                   1.0.11    2023-07-06 [2] CRAN (R 4.3.1)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.1)
# reticulate           * 1.32.0    2023-09-11 [2] CRAN (R 4.3.1)
# rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.1)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
# S4Arrays               1.0.6     2023-08-30 [2] Bioconductor
# S4Vectors            * 0.38.1    2023-05-02 [2] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.1)
# SingleCellExperiment * 1.22.0    2023-04-25 [2] Bioconductor
# SummarizedExperiment * 1.30.2    2023-06-06 [2] Bioconductor
# withr                  2.5.0     2022-03-03 [2] CRAN (R 4.3.1)
# XVector                0.40.0    2023-04-25 [2] Bioconductor
# zellkonverter        * 1.10.1    2023-05-23 [1] Bioconductor
# zlibbioc               1.46.0    2023-04-25 [2] Bioconductor
# 
# [1] /users/rphillip/R/4.3
# [2] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library
# 
# ─ Python configuration ───────────────────────────────────────────────────────
# python:         /users/rphillip/.virtualenvs/r-reticulate/bin/python
# libpython:      /jhpce/shared/community/core/conda_R/4.3/lib/libpython3.11.so
# pythonhome:     /users/rphillip/.virtualenvs/r-reticulate:/users/rphillip/.virtualenvs/r-reticulate
# version:        3.11.4 | packaged by conda-forge | (main, Jun 10 2023, 18:26:14) [GCC 12.2.0]
# numpy:          /users/rphillip/.virtualenvs/r-reticulate/lib/python3.11/site-packages/numpy
# numpy_version:  1.26.3
# anndata:        /users/rphillip/.virtualenvs/r-reticulate/lib/python3.11/site-packages/anndata
# 
# ──────────────────────────────────────────────────────────────────────────────
# 
