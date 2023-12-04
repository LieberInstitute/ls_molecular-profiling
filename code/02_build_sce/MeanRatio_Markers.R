##Goal: Evaluate the mean ratio method to help identify useful genes for identification of LS broadly
#as well as LS subclusters. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling

library(SingleCellExperiment)
library(DeconvoBuddies)
library(sessioninfo)
library(here)

#load the SingleCellExperiment object
load(here("processed-data","sce_with_CellType.rda"))

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
# colData names(61): Sample Barcode ... CellType_k_20_louvain
# CellType.Final
# reducedDimNames(18): GLMPCA_approx UMAP_15 ... tSNE_mnn_25 tSNE_mnn_50
# mainExpName: NULL
# altExpNames(0):

#load in the colors. 
load(here("processed-data","Final_CellTypes_colors_cb_Friendly.rda"),verbose = TRUE)
# Loading objects:
#     new_cluster_cols

#load the DEG lists. 
load(here("processed-data","markers_1vAll_ttest_CellTypeFinal_20Clusters.rda"),verbose = TRUE)
# Loading objects:
#     markers_1vALL_enrich_Final
load(here("processed-data","markers_pairwise_list_CellTypeFinal_20CellTypes.rda"),verbose = TRUE) 
# Loading objects:
#     markers_pairwise


#Now use the mean ratio method to help identify better markers. Function = get_mean_ratio2() from DeconvoBuddies
#From http://research.libd.org/DeconvoBuddies/articles/DeconvoBuddies.html#using-meanratio-to-find-cell-type-markers: 
# "To select genes specific for each cell type, you can evaluate the mean ratio for each gene x each cell type, where 
# mean ratio = mean(Expression of target cell type)/mean(Expression of highest non-target cell type)"



