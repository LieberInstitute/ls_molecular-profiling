#Goal: Compare gene expression signatures of the mouse and human LS
#This analysis will focus on using MNN to integrate. 
#Code modified from https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/51d15ef9f5f2c4c53f55e22e3fe467de1a724668/10x_NAc-n8_step04_cross-species_rnNAc_MNT.R#L4 and
#https://github.com/LieberInstitute/BLA_crossSpecies/blob/devel/code/costa_rds/05_species_comparisons/crossSpecies_PCA.R
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling

library(SingleCellExperiment)
library(batchelor)
library(here)

#Load the human and mouse matched SCE objects.
load(here("processed-data","human_mouse_matched_by_JAX.rda"),verbose = TRUE)
# Loading objects:
#     sce_mouse_sub
#     sce_human_sub

sce_mouse_sub
# class: SingleCellExperiment 
# dim: 16572 22860 
# metadata(1): Samples
# assays(3): counts binomial_pearson_residuals logcounts
# rownames(16572): ENSMUSG00000096351 ENSMUSG00000095567 ...
# ENSMUSG00000037772 ENSMUSG00000003526
# rowData names(9): source type ... mm.entrezIds JAX.geneID
# colnames(22860): 1_AAACCCAAGGTACATA-1 1_AAACCCACATCCGAGC-1 ...
# 4_TTTGTTGCATACAGCT-1 4_TTTGTTGGTCAAACGG-1
# colData names(17): Sample Barcode ... cellType.final cellType.broad
# reducedDimNames(4): GLMPCA_approx UMAP TSNE GLMPCA_50
# mainExpName: NULL
# altExpNames(0):

sce_human_sub
# class: SingleCellExperiment 
# dim: 16572 9225 
# metadata(1): Samples
# assays(3): counts binomial_deviance_residuals logcounts
# rownames(16572): ENSG00000187634 ENSG00000188976 ... ENSG00000214026
# ENSG00000100033
# rowData names(9): source type ... hs.entrezIds JAX.geneID
# colnames(9225): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 3_TTTGGTTTCTTCGACC-1 3_TTTGTTGTCCCGATCT-1
# colData names(61): Sample Barcode ... CellType_k_20_louvain
# CellType.Final
# reducedDimNames(18): GLMPCA_approx UMAP_15 ... tSNE_mnn_25 tSNE_mnn_50
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
# dim: 16572 32085 
# metadata(0):
#     assays(1): counts
# rownames(16572): ENSG00000187634 ENSG00000188976 ... ENSG00000214026
# ENSG00000100033
# rowData names(0):
#     colnames(32085): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 4_TTTGTTGCATACAGCT-1 4_TTTGTTGGTCAAACGG-1
# colData names(6): Sample Barcode ... detected CellType
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):

#compute logcounts
sce_combo <- multiBatchNorm(sce_combo, batch=sce_combo$Sample)






