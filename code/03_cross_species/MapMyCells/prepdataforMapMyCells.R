#Goal: Build anndata file for mapmycells
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling
#module load conda_R/4.4

library(SingleCellExperiment)
library(zellkonverter)
library(sessioninfo)
library(anndata)
library(here)

#Load human LS sce object
load(here("processed-data","02_build_sce","sce_celltype.rda"),verbose = TRUE)
# Loading objects:
#   sce

identical(rownames(colData(sce)),colnames(sce))
#[1] TRUE

######Python environment setup######
#Set up python environment with the reticulate package 
library(reticulate)

# create a new environment 
virtualenv_create("r-reticulate")

# install anndata
virtualenv_install("r-reticulate","anndata")

#Use the virtual environment 
use_virtualenv("r-reticulate")

# import anndata (it will be automatically discovered in "r-reticulate")
anndata <- import("anndata")
############################

######Convert gene ids to mouse#########
#Pull counts matrix from sce
#Map my cells correlation method expects raw counts. 
counts <- counts(sce)
dim(counts)
#[1] 36601  9225

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
convert_by_ensembl <- convert_by_ensembl[is.element(convert_by_ensembl$human_EnsemblID,
                                                    rownames(counts)),] 
dim(convert_by_ensembl)
#[1] 16507     2

#There is a row that is duplicated within the dataframe. Remove it. 
convert_by_ensembl <- convert_by_ensembl[!duplicated(convert_by_ensembl),]
dim(convert_by_ensembl)
#[1] 16506     2

# Subset data to include only genes with mouse othologs
new_order <- match(convert_by_ensembl$human_EnsemblID,rownames(counts))
new_order <- new_order[!is.na(new_order)]

counts_out <- counts[new_order,]
dim(counts_out)
#[1] 16506  9225

#Double check that everything is in the same order and no duplicates
identical(rownames(counts_out),convert_by_ensembl$human_EnsemblID)
#[1] TRUE

table(duplicated(rownames(counts_out)))
# FALSE 
# 16506 

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
# 16502 

#Transpose because mapmycells needs genes in column and cells in rows. 
counts_out <- t(counts_out)

# Convert to anndata format
h_LS_anndata <- AnnData(X = counts_out,
                        obs = data.frame(group = rownames(counts_out),
                                         row.names = rownames(counts_out)),
                        var = data.frame(group = colnames(counts_out),
                                         row.names = colnames(counts_out)))
# FutureWarning: The dtype argument is deprecated and will be removed in late 2024.
# warnings.warn(

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
# [1] "2024-09-16 13:02:31 EDT"
# user  system elapsed 
# 61.814   8.575 387.610 
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.4.0 Patched (2024-05-22 r86590)
# os       Rocky Linux 9.4 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-09-16
# pandoc   3.1.13 @ /jhpce/shared/community/core/conda_R/4.4/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version date (UTC) lib source
# abind                  1.4-5   2016-07-21 [2] CRAN (R 4.4.0)
# anndata              * 0.7.5.6 2023-03-17 [1] CRAN (R 4.4.0)
# assertthat             0.2.1   2019-03-21 [2] CRAN (R 4.4.0)
# basilisk               1.16.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# basilisk.utils         1.16.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# Biobase              * 2.64.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# BiocGenerics         * 0.50.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# cli                    3.6.3   2024-06-21 [2] CRAN (R 4.4.0)
# crayon                 1.5.3   2024-06-20 [2] CRAN (R 4.4.0)
# DelayedArray           0.30.1  2024-05-07 [2] Bioconductor 3.19 (R 4.4.0)
# dir.expiry             1.12.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# filelock               1.0.3   2023-12-11 [2] CRAN (R 4.4.0)
# GenomeInfoDb         * 1.40.1  2024-05-24 [2] Bioconductor 3.19 (R 4.4.0)
# GenomeInfoDbData       1.2.12  2024-05-23 [2] Bioconductor
# GenomicRanges        * 1.56.1  2024-06-12 [2] Bioconductor 3.19 (R 4.4.0)
# here                 * 1.0.1   2020-12-13 [2] CRAN (R 4.4.0)
# httr                   1.4.7   2023-08-15 [2] CRAN (R 4.4.0)
# IRanges              * 2.38.1  2024-07-03 [2] Bioconductor 3.19 (R 4.4.0)
# jsonlite               1.8.8   2023-12-04 [2] CRAN (R 4.4.0)
# lattice                0.22-6  2024-03-20 [3] CRAN (R 4.4.0)
# Matrix                 1.7-0   2024-04-26 [3] CRAN (R 4.4.0)
# MatrixGenerics       * 1.16.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# matrixStats          * 1.3.0   2024-04-11 [2] CRAN (R 4.4.0)
# png                    0.1-8   2022-11-29 [2] CRAN (R 4.4.0)
# R6                     2.5.1   2021-08-19 [2] CRAN (R 4.4.0)
# rappdirs               0.3.3   2021-01-31 [2] CRAN (R 4.4.0)
# Rcpp                   1.0.13  2024-07-17 [2] CRAN (R 4.4.0)
# reticulate           * 1.38.0  2024-06-19 [2] CRAN (R 4.4.0)
# rlang                  1.1.4   2024-06-04 [2] CRAN (R 4.4.0)
# rprojroot              2.0.4   2023-11-05 [2] CRAN (R 4.4.0)
# S4Arrays               1.4.1   2024-05-20 [2] Bioconductor 3.19 (R 4.4.0)
# S4Vectors            * 0.42.1  2024-07-03 [2] Bioconductor 3.19 (R 4.4.0)
# sessioninfo          * 1.2.2   2021-12-06 [2] CRAN (R 4.4.0)
# SingleCellExperiment * 1.26.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# SparseArray            1.4.8   2024-05-24 [2] Bioconductor 3.19 (R 4.4.0)
# SummarizedExperiment * 1.34.0  2024-05-01 [2] Bioconductor 3.19 (R 4.4.0)
# UCSC.utils             1.0.0   2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# withr                  3.0.1   2024-07-31 [2] CRAN (R 4.4.0)
# XVector                0.44.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# zellkonverter        * 1.14.1  2024-06-23 [1] Bioconductor 3.19 (R 4.4.0)
# zlibbioc               1.50.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# 
# [1] /users/rphillip/R/4.4
# [2] /jhpce/shared/community/core/conda_R/4.4/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.4/R/lib64/R/library
# 
# ─ Python configuration ───────────────────────────────────────────────────────────────────────────────────────
# python:         /users/rphillip/.virtualenvs/r-reticulate/bin/python
# libpython:      /usr/lib64/libpython3.9.so
# pythonhome:     /users/rphillip/.virtualenvs/r-reticulate:/users/rphillip/.virtualenvs/r-reticulate
# version:        3.9.18 (main, Jan 24 2024, 00:00:00)  [GCC 11.4.1 20231218 (Red Hat 11.4.1-3)]
# numpy:          /users/rphillip/.virtualenvs/r-reticulate/lib/python3.9/site-packages/numpy
# numpy_version:  2.0.2
# anndata:        /users/rphillip/.virtualenvs/r-reticulate/lib64/python3.9/site-packages/anndata
# 
# NOTE: Python version was forced by use_python() function
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────
# 

