#Goal: Investigate and cluster the integrated object. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling

library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(scater)
library(scran)
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
                                             c("#E69F00", "#009E73","#0072B2","#000000"))
names(human_cell_cols) <- unique(sce_harmony_Species$Human_CellType)

#Change the mouse color to grey. 
human_cell_cols["Mouse"] <- "#c6c6c6"

#Now Mouse
sce_harmony_Species$Mouse_CellType <- ifelse(sce_harmony_Species$Species == "Mouse",
                                             sce_harmony_Species$CellType_Species,
                                             "Human")

#Generate colors
mouse_cell_cols <- Polychrome::createPalette(length(unique(sce_harmony_Species$Mouse_CellType)),
                                             c("#E69F00", "#009E73","#0072B2","#000000"))
names(mouse_cell_cols) <- unique(sce_harmony_Species$Mouse_CellType)

#Change the mouse color to grey. 
mouse_cell_cols["Human"] <- "#c6c6c6"

#Generate the plots. 
###Human 
#Legend
human_tSNE <- plotReducedDim(object = sce_harmony_Species,
                             dimred = "tSNE_HARMONY",
                             colour_by = "Human_CellType",
                             point_alpha = 0.4) +
    scale_color_manual(values = human_cell_cols)
ggsave(plot = human_tSNE,
       filename = here("plots","Conservation","harmony_tSNE_Human_CellType_Colored_withLegend.pdf"))
#No Legend
human_tSNE_noLegend <- plotReducedDim(object = sce_harmony_Species,
                                      dimred = "tSNE_HARMONY",
                                      colour_by = "Human_CellType",
                                      point_alpha = 0.4) +
    scale_color_manual(values = human_cell_cols) +
    theme(legend.position = "none")
ggsave(plot = human_tSNE_noLegend,
       filename = here("plots","Conservation","harmony_tSNE_Human_CellType_Colored_withoutLegend.pdf"),
       height = 10,width = 10)

###Mouse 
#Legend
Mouse_tSNE <- plotReducedDim(object = sce_harmony_Species,
                             dimred = "tSNE_HARMONY",
                             colour_by = "Mouse_CellType",
                             point_alpha = 0.4) +
    scale_color_manual(values = mouse_cell_cols)
ggsave(plot = Mouse_tSNE,
       filename = here("plots","Conservation","harmony_tSNE_Mouse_CellType_Colored_withLegend.pdf"))
#No Legend
mouse_tSNE_noLegend <- plotReducedDim(object = sce_harmony_Species,
                                      dimred = "tSNE_HARMONY",
                                      colour_by = "Mouse_CellType",
                                      point_alpha = 0.4) +
    scale_color_manual(values = mouse_cell_cols) +
    theme(legend.position = "none")
ggsave(plot = mouse_tSNE_noLegend,
       filename = here("plots","Conservation","harmony_tSNE_Mouse_CellType_Colored_withoutLegend.pdf"),
       height = 10,width = 10)


#30,000 cells so will use higher k values. 
#jaccard + louvain is similar to seurat workflow. 
#Resolution=1 is the default
snn_k_50 <- buildSNNGraph(sce_harmony_Species, k = 50, use.dimred = "HARMONY",type="jaccard")
snn_k_75 <- buildSNNGraph(sce_harmony_Species, k = 75, use.dimred = "HARMONY",type="jaccard")
snn_k_100 <- buildSNNGraph(sce_harmony_Species, k = 100, use.dimred = "HARMONY",type="jaccard")

#Louvain clustering
#k=50
set.seed(1234)
clust_50 <- igraph::cluster_louvain(snn_k_50,resolution=1)$membership
table(clust_50)
#clust_50
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 2582 2798 1186  675  388 3452 1010 1434 1934  398 3161  686  321 1509  234  329 
# 17   18   19   20   21   22   23   24   25   26   27 
# 467  150 4302  110  712 1665  494   72  345   57  638 

#k=75
set.seed(1234)
clust_75 <- igraph::cluster_louvain(snn_k_75,resolution=1)$membership
table(clust_75)
# clust_75
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 2673 3017 1181  673  385 3572 1005 1439 1785  382 3298  619 1397  338  528 4634 
# 17   18   19   20   21   22 
# 703 1604  443  346  450  637 

#k=100
set.seed(1234)
clust_100 <- igraph::cluster_louvain(snn_k_100,resolution=1)$membership
table(clust_100)
# clust_100
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 2680 3025 1181  665  383 3573 1000 1453 1819  380 3568  338 1341  336  529  142 
# 17   18   19   20   21   22   23 
# 4513  644  699 1594  447  345  454 


