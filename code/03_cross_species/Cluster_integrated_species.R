#Goal: Investigate and cluster the integrated object. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling

library(SingleCellExperiment)
library(DeconvoBuddies)
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
clust_50
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
# 2673 3034 1181  673  385 3572 1005 1439 1773  380 3590  324 1390  338  528  147 
# 17   18   19   20   21   22   23 
# 4487  703 1604  445  346  450  642 

#k=100
set.seed(1234)
clust_100 <- igraph::cluster_louvain(snn_k_100,resolution=1)$membership
table(clust_100)
# clust_100
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 2680 3025 1181  665  383 3573 1000 1453 1819  380 3568  338 1341  336  529  142 
# 17   18   19   20   21   22   23 
# 4513  644  699 1594  447  345  454 

#Add cluster information to object
sce_harmony_Species$k_50_louvain_1 <- factor(clust_50)
sce_harmony_Species$k_75_louvain_1 <- factor(clust_75)
sce_harmony_Species$k_100_louvain_1 <- factor(clust_100)

#Plot the clusters on the tSNE with 15 dimensions. 
for(i in c("k_50_louvain_1",
           "k_75_louvain_1",
           "k_100_louvain_1")){
    print(i)
    x <- plotReducedDim(object = sce_harmony_Species,
                        dimred = "tSNE_HARMONY",
                        colour_by = i,text_by = i) +
        ggtitle(i) +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(plot = x,filename = here("plots","Conservation","Clustering",paste0(i,".pdf")))
}

#k=75 looks best so far. 
#Now plot expression of marker genes to identify clusters. 
#First need to add in gene name from ensembl gene id. To do this load in the original humnan sce and 
#merge the rowData
load(here("processed-data","sce_with_CellType.rda"),verbose = TRUE)
# Loading objects:
#     sce
gene_names <- rowData(sce)[,c("gene_id","gene_name","gene_type")]

#Add a gene id column to sce_harmony_Species 
rowData(sce_harmony_Species)$gene_id <- row.names(rowData(sce_harmony_Species))
nrow(rowData(sce_harmony_Species))
#[1] 16562
rowData(sce_harmony_Species) <- merge(rowData(sce_harmony_Species),
                                      gene_names,
                                      by = "gene_id")
nrow(rowData(sce_harmony_Species))
#[1] 16562

#Pull ensembl IDs for these genes.
genes <- c("SYT1","SNAP25", #Pan-neuronal
           "GAD1","GAD2","SLC32A1", #GABAergic
           "SLC17A6","SLC17A7","SLC17A8", #Glutamatergic
           "TRPC4","DGKG","HOMER2","PTPN3","TRHDE","CPNE7","NRP1","FREM2", #Lateral Septum markers from mouse
           "ELAVL2","TRPC5", #Medial Septum markers from mouse
           "RARB","BCL11B","PPP1R1B", #Broad striatal + MSN markers 
           "FXYD6", #Broad septal
           "OPRM1","DRD1", #Striosome markers
           "MBP","MOBP", #Oligodendrocyte
           "CD74", "CSF1R", "C3", #Microglia
           "GFAP", "TNC", "AQP4", "SLC1A2", #Astrocytes
           "CLDN5", "FLT1", "VTN", #Endothelial
           "COL1A2", "TBX18", "RBPMS", #Mural
           "SKAP1", "ITK", "CD247", #Tcell
           "CD163", "SIGLEC1", "F13A1", #Macrophage
           "FOXJ1","PIFO","CFAP44", #Ependymal
           "PDGFRA", "VCAN", "CSPG4") #Polydendrocytes

for(i in genes){
    print(i)
    x <- rowData(sce_harmony_Species)[which(rowData(sce_harmony_Species)$gene_name == i),"gene_id"]
    gene_plot <- plotReducedDim(object = sce_harmony_Species,
                                dimred = "tSNE_HARMONY",
                                color_by = x) +
        scale_color_gradientn(colours = c("lightgrey","orange","red")) +
        ggtitle(paste0(i,"\n",
                       rowData(sce_harmony_Species)[which(rowData(sce_harmony_Species)$gene_name == i),"gene_id"])) +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(plot = gene_plot,filename = here("plots","Conservation",
                                            "Clustering","k_75_louvain_1_expression",
                                            paste0(i,"_cross_species_k_75.pdf")))
}


plotReducedDim(object = sce_harmony_Species,
               dimred = "tSNE_HARMONY",
               color_by = "ENSG00000149295") +
    scale_color_gradientn(colours = c("lightgrey","orange","red"))

#Run a quick DEG analysis to help with cluster identification. 
markers_1vALL_enrich <- findMarkers_1vAll(sce_harmony_Species, 
                                          assay_name = "logcounts", 
                                          cellType_col = "k_75_louvain_1", 
                                          mod = "~Sample")
# 1 - '2024-01-08 14:43:56.808779
# 2 - '2024-01-08 14:46:04.164763
# 3 - '2024-01-08 14:48:16.745078
# 4 - '2024-01-08 14:50:30.802523
# 5 - '2024-01-08 14:52:45.587948
# 6 - '2024-01-08 14:55:03.956395
# 7 - '2024-01-08 14:57:09.731182
# 8 - '2024-01-08 14:59:25.046667
# 9 - '2024-01-08 15:01:39.360235
# 10 - '2024-01-08 15:03:55.02235
# 11 - '2024-01-08 15:06:06.442347
# 12 - '2024-01-08 15:08:06.443018
# 13 - '2024-01-08 15:10:04.40993
# 14 - '2024-01-08 15:11:47.439724
# 15 - '2024-01-08 15:13:36.992938
# 16 - '2024-01-08 15:15:30.60645
# 17 - '2024-01-08 15:17:21.666588
# 18 - '2024-01-08 15:19:13.508835
# 19 - '2024-01-08 15:21:06.690871
# 20 - '2024-01-08 15:22:57.691572
# 21 - '2024-01-08 15:24:50.009232
# 22 - '2024-01-08 15:26:40.621126
# 23 - '2024-01-08 15:28:33.457227
# Building Table - 2024-01-08 15:30:06.178589
# ** Done! **

#Merge the two dataframes. 
markers_gene_name <- merge(x    = as.data.frame(markers_1vALL_enrich),
                           y    = as.data.frame(rowData(sce_harmony_Species)[,c("gene_id","gene_name")]),
                           by.x = "gene",
                           by.y = "gene_id")

#Subset for top 50 markers for each cluster
markers_gene_name <- subset(markers_gene_name,subset=(rank_marker %in% 1:50))


#Order by cellType.target
markers_gene_name <- markers_gene_name[order(markers_gene_name$rank_marker,decreasing = FALSE),]

#save markers_df and object
save(markers_gene_name,file = here("processed-data","integrated_markers_df.rda"))
save(sce_harmony_Species,file = here("processed-data","sce_integrated.rda"))

#Would also help to figure out what makes up each cluster. 
for(i in unique(sce_harmony_Species$k_75_louvain_1)){
    print(unique(subset(colData(sce_harmony_Species),subset=(k_75_louvain_1 == i))$CellType_Species))
}
#Need to figure out a better way to do this. 


#Start the annotation. 
# 1 - LS_Inh_A
# 2 - Drd1_MSN_A
# 3 - Excit_A
# 4 - Polydendrocyte
# 5 - Microglia
# 6 - Oligodendrocyte
# 7 - Ependymal
# 8 - Drd1_MSN_Patch
# 9 - Drd2_MSN
# 10 - MS_Inh_A
# 11 - Septal_Inh_A
# 12 - MS_Inh_B
# 13 - Str_Inh_A 
# 14 - Excit_B
# 15 - Mural
# 16 - Astrocyte_1
# 17 - Astrocyte_2
# 18 - 
# 19 - Astrocyte_3
# 20 - 
# 21
# 22
# 23






