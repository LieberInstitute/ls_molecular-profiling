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
mean_ratios_CellTypeFinal <- get_mean_ratio2(sce,
                                             cellType_col = "CellType.Final")

#Check out structure and first couple lines of the dataframe. 
dim(mean_ratios_CellTypeFinal)
# [1] 89642     8

mean_ratios_CellTypeFinal[1:5,1:8]
# A tibble: 5 × 8
# gene    cellType.target mean.target cellType  mean ratio rank_ratio anno_ratio
# <chr>   <chr>                 <dbl> <chr>    <dbl> <dbl>      <int> <chr>     
# 1 ENSG00… LS_Inh_A              0.807 LS_Inh_G 0.372  2.17          1 LS_Inh_A/…
# 2 ENSG00… LS_Inh_A              2.00  LS_Inh_B 1.09   1.83          2 LS_Inh_A/…
# 3 ENSG00… LS_Inh_A              1.06  Excit_A  0.593  1.78          3 LS_Inh_A/…
# 4 ENSG00… LS_Inh_A              4.21  MS_Inh_E 2.38   1.77          4 LS_Inh_A/…
# 5 ENSG00… LS_Inh_A              1.87  MS_Inh_A 1.07   1.75          5 LS_Inh_A/…

#Add the gene symbol. 
#The get_mean_ratio2() function has an add_symbol flag, but using that didn't add the symbol.  
mean_ratios_CellTypeFinal_symbol <- merge(x = as.data.frame(mean_ratios_CellTypeFinal),
                                          y = as.data.frame(rowData(sce)[,c("gene_id","gene_name")]),
                                          by.x = "gene",
                                          by.y = "gene_id")
dim(mean_ratios_CellTypeFinal_symbol)
# [1] 89642     9
