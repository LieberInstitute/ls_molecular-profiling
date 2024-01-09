#Goal: Compare gene expression signatures of the mouse and human LS
#This analysis will focus on using MNN to integrate. 
#Code modified from https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/51d15ef9f5f2c4c53f55e22e3fe467de1a724668/10x_NAc-n8_step04_cross-species_rnNAc_MNT.R#L4 and
#https://github.com/LieberInstitute/BLA_crossSpecies/blob/devel/code/costa_rds/05_species_comparisons/crossSpecies_PCA.R
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling

library(SingleCellExperiment)
library(sessioninfo)
library(batchelor)
library(harmony)
library(scater)
library(here)
library(scry)

#Load the human and mouse matched SCE objects.
load(here("processed-data","human_mouse_matched_by_JAX.rda"),verbose = TRUE)
# Loading objects:
#     sce_mouse_sub
#     sce_human_sub

#Make sure the mixed and doublet populations are gone. 
levels(sce_mouse_sub$cellType.final)
# [1] "Astro"         "Chol_Ex.D"     "ChP"           "Endo"         
# [5] "Ependymal"     "IoC_In.E"      "LS_In.C"       "LS_In.D"      
# [9] "LS_In.M"       "LS_In.N"       "LS_In.O"       "LS_In.P"      
# [13] "LS_In.Q"       "LS_In.R"       "Micro"         "MS_In.J"      
# [17] "MS_In.K"       "Mural"         "Neuroblast"    "Oligo"        
# [21] "OPC"           "OPC_COP"       "Sept_In.G"     "Sept_In.I"    
# [25] "Str_In.A"      "Str_In.F"      "Str_In.H"      "Str_In.L"     
# [29] "Thal_Ex.B"     "TNoS_Ex.A"     "TT.IG.SH_Ex.C" "TT.IG.SH_Ex.E"
# [33] "TT.IG.SH_Ex.F"

sce_human_sub
# class: SingleCellExperiment 
# dim: 16562 9225 
# metadata(1): Samples
# assays(3): counts binomial_deviance_residuals logcounts
# rownames(16562): ENSG00000187634 ENSG00000188976 ... ENSG00000214026
# ENSG00000100033
# rowData names(9): source type ... hs.entrezIds JAX.geneID
# colnames(9225): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 3_TTTGGTTTCTTCGACC-1 3_TTTGTTGTCCCGATCT-1
# colData names(61): Sample Barcode ... CellType_k_20_louvain
# CellType.Final
# reducedDimNames(18): GLMPCA_approx UMAP_15 ... tSNE_mnn_25 tSNE_mnn_50
# mainExpName: NULL
# altExpNames(0):

sce_mouse_sub
# class: SingleCellExperiment 
# dim: 16562 21884 
# metadata(1): Samples
# assays(3): counts binomial_pearson_residuals logcounts
# rownames(16562): ENSMUSG00000096351 ENSMUSG00000095567 ...
# ENSMUSG00000037772 ENSMUSG00000003526
# rowData names(9): source type ... mm.entrezIds JAX.geneID
# colnames(21884): 1_AAACCCAAGGTACATA-1 1_AAACCCACATCCGAGC-1 ...
# 4_TTTGTTGCATACAGCT-1 4_TTTGTTGGTCAAACGG-1
# colData names(17): Sample Barcode ... cellType.final cellType.broad
# reducedDimNames(4): GLMPCA_approx UMAP TSNE GLMPCA_50
# mainExpName: NULL
# altExpNames(0):

#Everything should be in order, but check to make sure? 
all(rowData(sce_human_sub)$JAX.geneID == rowData(sce_mouse_sub)$JAX.geneID) 
# [1] TRUE

#Add coldata column that corresponds to the species. 
sce_human_sub$Species <- "Human"
sce_mouse_sub$Species <- "Mouse"

#Clean up celltype column names. 
sce_human_sub$CellType <- sce_human_sub$CellType.Final
sce_mouse_sub$CellType <- sce_mouse_sub$cellType.final

#subset the objects 
colData(sce_human_sub) <- colData(sce_human_sub)[,c("Sample","Barcode","Species","sum","detected","CellType")]
colData(sce_mouse_sub) <- colData(sce_mouse_sub)[,c("Sample","Barcode","Species","sum","detected","CellType")]

#To combine the objects, concatenate the count data + colData
#Then, simply create the object with SingleCellExperiment()
combo_counts <- cbind(assay(sce_human_sub,"counts"),
                      assay(sce_mouse_sub,"counts"))
combo_coldata <- rbind(colData(sce_human_sub),
                       colData(sce_mouse_sub))

#Create the combo object
sce_combo <- SingleCellExperiment(assays=list(counts=combo_counts), colData=combo_coldata)

sce_combo
# class: SingleCellExperiment 
# dim: 16562 31109 
# metadata(0):
#     assays(1): counts
# rownames(16562): ENSG00000187634 ENSG00000188976 ... ENSG00000214026
# ENSG00000100033
# rowData names(0):
#     colnames(31109): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 4_TTTGTTGCATACAGCT-1 4_TTTGTTGGTCAAACGG-1
# colData names(6): Sample Barcode ... detected CellType
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):
    
#compute logcounts
sce_combo <- multiBatchNorm(sce_combo, batch=sce_combo$Sample)

##########Feature selection. 
#Deviance feature selection
sce_combo <- devianceFeatureSelection(sce_combo,
                                      assay = "counts",
                                      fam = "binomial",
                                      sorted = FALSE,
                                      batch = as.factor(sce_combo$Sample))

pdf(here("plots","Conservation","featureSelxn_binomialDeviance-byGene_numericcutoffs.pdf"))
plot(sort(rowData(sce_combo)$binomial_deviance, decreasing = T),
     type = "l", xlab = "ranked genes",
     ylab = "binomial deviance"
)
abline(v = 2000,lty = 2, col = "red")
dev.off()

#2000 should be good. 
sce_combo <- nullResiduals(sce_combo,
                           assay = "counts", 
                           fam   = "binomial", 
                           type  = "deviance")

#Take top 2000 highly deviant genes
hdgs <- rownames(sce_combo)[order(rowData(sce_combo)$binomial_deviance, decreasing = T)][1:2000]

#current object doesn't include gene name info. Pull from the sce_human_sub object. 
hdgs.symbols <- rowData(sce_human_sub)$gene_name[match(hdgs, rowData(sce_human_sub)$gene_id)]

#PCA
sce_combo <- runPCA(sce_combo,
                    exprs_values = "binomial_deviance_residuals",
                    subset_row = hdgs, 
                    ncomponents = 100,
                    name = "GLMPCA_approx")


#PCA plot of top 6 PCs
PCA_plots <- plotReducedDim(sce_combo,
                            dimred = "GLMPCA_approx", 
                            colour_by = "Sample",
                            ncomponents = 6, 
                            point_alpha = 0.3)
ggsave(PCA_plots,filename = here("plots","Conservation","multi_PCAs.png"))

#Run the UMAP + tSNE with 50 dimensions. 
#UMAP
set.seed(1234)
sce_combo <- runUMAP(sce_combo,
                     dimred = "GLMPCA_approx",
                     n_dimred = 50,
                     name = "UMAP_50")
#tSNE
set.seed(1234)
sce_combo <- runTSNE(sce_combo,
                     dimred = "GLMPCA_approx",
                     n_dimred = 50,
                     name = "tSNE_50")

#tSNE
tSNE_cross_species <- plotReducedDim(sce_combo,
                                     dimred = "tSNE_50", 
                                     colour_by = "Species",
                                     point_alpha = 0.3)
ggsave(tSNE_cross_species,filename = here("plots",
                                          "Conservation",
                                          "tSNE_cross_species.png"))


#As expected, tSNE split by species. Next run fast fastMNN
mnn.hold <- fastMNN(sce_combo, 
                    batch=sce_combo$Sample,
                    subset.row=hdgs, 
                    d=50,
                    correct.all=TRUE, 
                    get.variance=TRUE,
                    BSPARAM=BiocSingular::IrlbaParam())

#Add mnn to sce_combo object
reducedDim(sce_combo, "PCA_MNN") <- reducedDim(mnn.hold, "corrected")

#Rerun the UMAP and tSNE
set.seed(1234)
sce_combo <- runUMAP(sce_combo,
                     dimred = "PCA_MNN",
                     n_dimred = 50,
                     name = "UMAP_corrected_50")
#tSNE
set.seed(1234)
sce_combo <- runTSNE(sce_combo,
                     dimred = "PCA_MNN",
                     n_dimred = 50,
                     name = "tSNE_corrected_50")


#Plot the tSNE and UMAP
#UMAP by species
umap_cross_species <- plotReducedDim(sce_combo,
                                     dimred = "UMAP_corrected_50", 
                                     colour_by = "Species",
                                     point_alpha = 0.3)
ggsave(umap_cross_species,filename = here("plots",
                                          "Conservation",
                                          "umap_cross_species_mnn_corrected.png"))

#tSNE
tSNE_cross_species <- plotReducedDim(sce_combo,
                                     dimred = "tSNE_corrected_50", 
                                     colour_by = "Species",
                                     point_alpha = 0.3)
ggsave(tSNE_cross_species,filename = here("plots",
                                          "Conservation",
                                          "tSNE_cross_species_mnn_corrected.png"))

save(sce_combo,file = here("processed-data","sce_combo.rda"))

#By celltype
#First create a new column within sce_combo that is celltype_species
sce_combo$CellType_Species <- paste(sce_combo$CellType,sce_combo$Species,sep = "_")

#Assign colors
cluster_cols <- Polychrome::createPalette(length(unique(sce_combo$CellType_Species)),
                                          c("#D81B60", "#1E88E5","#FFC107","#009E73"))
names(cluster_cols) <- unique(sce_combo$CellType_Species)

#UMAP
umap_cross_species_ct <- plotReducedDim(sce_combo,
                                        dimred = "UMAP_corrected_50", 
                                        colour_by = "CellType_Species",
                                        text_by = "CellType_Species",
                                        point_alpha = 0.3) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_color_manual(values = cluster_cols)
ggsave(umap_cross_species_ct,filename = here("plots",
                                             "Conservation",
                                             "umap_cross_species_mnn_corrected_byCellType.png"),
       height = 10, width = 10)

#tSNE
tSNE_cross_species_ct <- plotReducedDim(sce_combo,
                                        dimred = "tSNE_corrected_50", 
                                        colour_by = "CellType_Species",
                                        text_by = "CellType_Species",
                                        point_alpha = 0.3) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_color_manual(values = cluster_cols)
ggsave(tSNE_cross_species_ct,filename = here("plots",
                                             "Conservation",
                                             "tSNE_cross_species_mnn_corrected_byCellType.png"),
       height = 10, width = 10)

#These look okay. Will also try Harmony. 
#Add a duplicate reducedDim that is only called PCA. This is required for RunHarmony 
reducedDim(sce_combo,"PCA") <- reducedDim(sce_combo,"GLMPCA_approx")
sce_harmony_Species <- RunHarmony(sce_combo, group.by.vars = "Species", verbose = TRUE)
# Transposing data matrix
# Hard k-means centroids initialization
# Harmony 1/10
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# Harmony 2/10
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# Harmony converged after 2 iterations

#Run tSNE adn UMAP with harmony
#UMAP
set.seed(1234)
sce_harmony_Species <- runUMAP(sce_harmony_Species,
                               dimred = "HARMONY",
                               name = "UMAP_HARMONY")
#tSNE
set.seed(1234)
sce_harmony_Species <- runTSNE(sce_harmony_Species,
                               dimred = "HARMONY",
                               name = "tSNE_HARMONY")
#Plot the tSNE and UMAP
#UMAP by species
umap_Harmony <- plotReducedDim(sce_harmony_Species,
                               dimred = "UMAP_HARMONY", 
                               colour_by = "Species",
                               point_alpha = 0.3)
ggsave(umap_Harmony,filename = here("plots",
                                    "Conservation",
                                    "umap_cross_species_Harmony_corrected.png"))

#tSNE
tSNE_Harmony <- plotReducedDim(sce_harmony_Species,
                               dimred = "tSNE_HARMONY", 
                               colour_by = "Species",
                               point_alpha = 0.3)
ggsave(tSNE_Harmony,filename = here("plots",
                                    "Conservation",
                                    "tSNE_cross_species_Harmony.png"))

save(sce_harmony_Species,file = here("processed-data","sce_harmony_species.rda"))

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# [1] "Reproducibility information:"
# [1] "2024-01-09 11:55:36 EST"
# user   system  elapsed 
# 2234.974   91.092 3555.720 
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
# date     2024-01-09
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/pandoc
# 
# ─ Packages ─────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [2] CRAN (R 4.3.1)
# batchelor            * 1.16.0    2023-04-25 [2] Bioconductor
# beachmat               2.16.0    2023-04-25 [2] Bioconductor
# beeswarm               0.4.0     2021-06-01 [2] CRAN (R 4.3.1)
# Biobase              * 2.60.0    2023-04-25 [2] Bioconductor
# BiocGenerics         * 0.46.0    2023-04-25 [2] Bioconductor
# BiocNeighbors          1.18.0    2023-04-25 [2] Bioconductor
# BiocParallel           1.34.2    2023-05-22 [2] Bioconductor
# BiocSingular           1.16.0    2023-04-25 [2] Bioconductor
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.1)
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.1)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.3.1)
# colorout             * 1.3-0.1   2023-12-01 [1] Github (jalvesaq/colorout@deda341)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.3.1)
# cowplot                1.1.1     2020-12-30 [2] CRAN (R 4.3.1)
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.1)
# DelayedArray           0.26.7    2023-07-28 [2] Bioconductor
# DelayedMatrixStats     1.22.6    2023-08-28 [2] Bioconductor
# dplyr                  1.1.3     2023-09-03 [2] CRAN (R 4.3.1)
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
# harmony              * 1.0.1     2023-09-20 [2] CRAN (R 4.3.1)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.3.1)
# igraph                 1.5.1     2023-08-10 [2] CRAN (R 4.3.1)
# IRanges              * 2.34.1    2023-06-22 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.3.1)
# labeling               0.4.3     2023-08-29 [2] CRAN (R 4.3.1)
# lattice                0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.3.1)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.3.1)
# Matrix                 1.6-1.1   2023-09-18 [3] CRAN (R 4.3.1)
# MatrixGenerics       * 1.12.3    2023-07-30 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.3.1)
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.3.1)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.3.1)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.3.1)
# Polychrome             1.5.1     2022-05-03 [1] CRAN (R 4.3.1)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.1)
# ragg                   1.2.5     2023-01-12 [2] CRAN (R 4.3.1)
# Rcpp                 * 1.0.11    2023-07-06 [2] CRAN (R 4.3.1)
# RcppAnnoy              0.0.21    2023-07-02 [2] CRAN (R 4.3.1)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.1)
# ResidualMatrix         1.10.0    2023-04-25 [2] Bioconductor
# RhpcBLASctl            0.23-42   2023-02-11 [2] CRAN (R 4.3.1)
# rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.1)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.3.1)
# Rtsne                  0.16      2022-04-17 [2] CRAN (R 4.3.1)
# S4Arrays               1.0.6     2023-08-30 [2] Bioconductor
# S4Vectors            * 0.38.1    2023-05-02 [2] Bioconductor
# ScaledMatrix           1.8.1     2023-05-03 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.3.1)
# scater               * 1.28.0    2023-04-25 [2] Bioconductor
# scatterplot3d          0.3-44    2023-05-05 [1] CRAN (R 4.3.1)
# scry                 * 1.12.0    2023-04-25 [2] Bioconductor
# scuttle              * 1.10.2    2023-08-03 [2] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.1)
# SingleCellExperiment * 1.22.0    2023-04-25 [2] Bioconductor
# sparseMatrixStats      1.12.2    2023-07-02 [2] Bioconductor
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