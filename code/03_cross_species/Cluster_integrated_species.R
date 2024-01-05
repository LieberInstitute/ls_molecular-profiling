#Goal: Investigate and cluster the integrated object. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling

library(SingleCellExperiment)
library(sessioninfo)
library(scater)
library(here)

#Load the integrated object. 
load(here("processed-data","sce_harmony_species.rda"),verbose = TRUE)
# Loading objects:
#     sce_harmony_Species


sce_harmony_Species
# class: SingleCellExperiment 
# dim: 16562 31109 
# metadata(0):
#     assays(3): counts logcounts binomial_deviance_residuals
# rownames(16562): ENSG00000187634 ENSG00000188976 ... ENSG00000214026
# ENSG00000100033
# rowData names(1): binomial_deviance
# colnames(31109): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 4_TTTGTTGCATACAGCT-1 4_TTTGTTGGTCAAACGG-1
# colData names(8): Sample Barcode ... sizeFactor CellType_Species
# reducedDimNames(10): GLMPCA_approx UMAP_50 ... UMAP_HARMONY
# tSNE_HARMONY
# mainExpName: NULL
# altExpNames(0):

#Would be good to have plots where one species is greyed out and the other species clusters are colored. 
#First make two coldata columns that will contain the information needed to make the plots. 
####Human first
sce_harmony_Species$Human_CellType <- ifelse(sce_harmony_Species$Species == "Human",
                                             sce_harmony_Species$CellType_Species,
                                             "Mouse")

#Generate colors
human_cell_cols <- Polychrome::createPalette(length(unique(sce_harmony_Species$Human_CellType)),
                                             c("#D81B60", "#1E88E5","#FFC107","#009E73"))
names(human_cell_cols) <- unique(sce_harmony_Species$Human_CellType)

#Change the mouse color to grey. 
human_cell_cols["Mouse"] <- "#767574"

#Now Mouse
sce_harmony_Species$Mouse_CellType <- ifelse(sce_harmony_Species$Species == "Mouse",
                                             sce_harmony_Species$CellType_Species,
                                             "Human")

#Generate colors
mouse_cell_cols <- Polychrome::createPalette(length(unique(sce_harmony_Species$Mouse_CellType)),
                                             c("#D81B60", "#1E88E5","#FFC107","#009E73"))
names(mouse_cell_cols) <- unique(sce_harmony_Species$Mouse_CellType)

#Change the mouse color to grey. 
human_cell_cols["Human"] <- "#767574"







