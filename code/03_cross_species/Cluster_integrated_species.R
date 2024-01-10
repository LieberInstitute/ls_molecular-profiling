#Goal: Investigate and cluster the integrated object. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling

library(SingleCellExperiment)
library(DeconvoBuddies)
library(sessioninfo)
library(reshape2)
library(ggridges)
library(ggplot2)
library(scater)
library(scran)
library(here)
library(lisi)

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

###################################
############# HARMONY #############
###################################
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

#No Legend, labelled 
human_tSNE_labeled <- plotReducedDim(object = sce_harmony_Species,
                                     dimred = "tSNE_HARMONY",
                                     colour_by = "Human_CellType",
                                     text_by = "Human_CellType",
                                     point_alpha = 0.4) +
    scale_color_manual(values = human_cell_cols) +
    theme(legend.position = "none")
ggsave(plot = human_tSNE_labeled,
       filename = here("plots","Conservation","harmony_tSNE_Human_CellType_Colored_labeled.pdf"),
       height = 13,width = 13)

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
#No Legend, labelled 
mouse_tSNE_labeled <- plotReducedDim(object = sce_harmony_Species,
                                     dimred = "tSNE_HARMONY",
                                     colour_by = "Mouse_CellType",
                                     text_by = "Mouse_CellType",
                                     point_alpha = 0.4) +
    scale_color_manual(values = mouse_cell_cols) +
    theme(legend.position = "none")
ggsave(plot = mouse_tSNE_labeled,
       filename = here("plots","Conservation","harmony_tSNE_mouse_CellType_Colored_labeled.pdf"),
       height = 13,width = 13)

#30,000 cells so will use higher k values. 
#jaccard + louvain is similar to seurat workflow. 
#Resolution=1 is the default
snn_k_10 <- buildSNNGraph(sce_harmony_Species, k = 10, use.dimred = "HARMONY",type="jaccard")
snn_k_25 <- buildSNNGraph(sce_harmony_Species, k = 25, use.dimred = "HARMONY",type="jaccard")
snn_k_50 <- buildSNNGraph(sce_harmony_Species, k = 50, use.dimred = "HARMONY",type="jaccard")
snn_k_75 <- buildSNNGraph(sce_harmony_Species, k = 75, use.dimred = "HARMONY",type="jaccard")
snn_k_100 <- buildSNNGraph(sce_harmony_Species, k = 100, use.dimred = "HARMONY",type="jaccard")


#Louvain clustering
#k=10
set.seed(1234)
clust_10 <- igraph::cluster_louvain(snn_k_10,resolution=1)$membership
table(clust_10)
# clust_10
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 3790 1083 1353  548  927 1521 1037  591 5316 1485  930 3717 3214 5025  276  296 

#k=25
set.seed(1234)
clust_25 <- igraph::cluster_louvain(snn_k_25,resolution=1)$membership
table(clust_25)
# clust_25
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 3891 1050 1340  534 1890  162 1075  571 7015 1433   38 6809 4871  150  280 

#k=50
set.seed(1234)
clust_50 <- igraph::cluster_louvain(snn_k_50,resolution=1)$membership
table(clust_50)
# clust_50
# 1    2    3    4    5    6    7    8    9   10   11   12 
# 3882 1588 1322  523 2129 7235 1062 1400 6448 4916  342  262 

#k=75
set.seed(1234)
clust_75 <- igraph::cluster_louvain(snn_k_75,resolution=1)$membership
table(clust_75)
# clust_75
# 1    2    3    4    5    6    7    8    9   10   11 
# 3864 1552 1310  508  604 4164 1047 5148 2898 4967 5047 

#k=100
set.seed(1234)
clust_100 <- igraph::cluster_louvain(snn_k_100,resolution=1)$membership
table(clust_100)
# clust_100
# 1    2    3    4    5    6    7    8    9   10   11 
# 3889 1516 1309  492  558 3283 1041 6012 5038 3005 4966 

#Add cluster information to object
sce_harmony_Species$k_10_louvain_1 <- factor(clust_10)
sce_harmony_Species$k_25_louvain_1 <- factor(clust_25)
sce_harmony_Species$k_50_louvain_1 <- factor(clust_50)
sce_harmony_Species$k_75_louvain_1 <- factor(clust_75)
sce_harmony_Species$k_100_louvain_1 <- factor(clust_100)

#Plot the clusters on the tSNE
for(i in c("k_10_louvain_1",
           "k_25_louvain_1",
           "k_50_louvain_1",
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
           "DRD2","ADORA2A","PENK", #D2 MARKERS
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
                                            "Clustering","Expression",
                                            paste0(i,"_cross_species.pdf")))
}


#save object
save(sce_harmony_Species,file = here("processed-data","sce_integrated.rda"))

#Would also help to figure out what makes up each cluster for each iteration of clustering. 
for(i in c("k_10_louvain_1","k_25_louvain_1","k_50_louvain_1",
           "k_75_louvain_1","k_100_louvain_1")){
    for(l in unique(levels(colData(sce_harmony_Species)[,i]))){
        rows <- which(colData(sce_harmony_Species)[,i] == l)
        data <- as.data.frame(table(colData(sce_harmony_Species)[rows,"CellType_Species"]))
        data$prop <- data$Freq/sum(data$Freq)*100
        prop_plot <- ggplot(data, aes(x=reorder(Var1,prop), y=prop, fill=Var1)) +
            geom_bar(stat="identity") +
            labs(x = "Cell Type",
                 y = "Proportion") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "none")
        ggsave(plot = prop_plot,
               filename = here("plots","Conservation","Cluster_Proportion_barplots","Harmony",i,
                               paste0(l,"_prop_barplot.pdf")))
 
}}

#unsure about the annotation. Mouse and Human striatal MSNs are not integrating well with HARMONY
#will check MNN integration. 

###################################
############### MNN ###############
###################################
#Cluster the mnn integration. 
#PCA_MNN is the mnn corrected space. 
#tSNE_corrected_50 and umao_corrected_50
#First make same plots as before with Harmony. 
#Generate the plots. 
###Human 
#Legend
human_tSNE <- plotReducedDim(object = sce_harmony_Species,
                             dimred = "tSNE_corrected_50",
                             colour_by = "Human_CellType",
                             point_alpha = 0.4) +
    scale_color_manual(values = human_cell_cols)
ggsave(plot = human_tSNE,
       filename = here("plots","Conservation","mnn_tSNE_Human_CellType_Colored_withLegend.pdf"))
#No Legend
human_tSNE_noLegend <- plotReducedDim(object = sce_harmony_Species,
                                      dimred = "tSNE_corrected_50",
                                      colour_by = "Human_CellType",
                                      point_alpha = 0.4) +
    scale_color_manual(values = human_cell_cols) +
    theme(legend.position = "none")
ggsave(plot = human_tSNE_noLegend,
       filename = here("plots","Conservation","mnn_tSNE_Human_CellType_Colored_withoutLegend.pdf"),
       height = 10,width = 10)

#No Legend, labelled 
human_tSNE_labeled <- plotReducedDim(object = sce_harmony_Species,
                                     dimred = "tSNE_corrected_50",
                                     colour_by = "Human_CellType",
                                     text_by = "Human_CellType",
                                     point_alpha = 0.4) +
    scale_color_manual(values = human_cell_cols) +
    theme(legend.position = "none")
ggsave(plot = human_tSNE_labeled,
       filename = here("plots","Conservation","mnn_tSNE_Human_CellType_Colored_labeled.pdf"),
       height = 13,width = 13)

###Mouse 
#Legend
Mouse_tSNE <- plotReducedDim(object = sce_harmony_Species,
                             dimred = "tSNE_corrected_50",
                             colour_by = "Mouse_CellType",
                             point_alpha = 0.4) +
    scale_color_manual(values = mouse_cell_cols)
ggsave(plot = Mouse_tSNE,
       filename = here("plots","Conservation","mnn_tSNE_Mouse_CellType_Colored_withLegend.pdf"))
#No Legend
mouse_tSNE_noLegend <- plotReducedDim(object = sce_harmony_Species,
                                      dimred = "tSNE_corrected_50",
                                      colour_by = "Mouse_CellType",
                                      point_alpha = 0.4) +
    scale_color_manual(values = mouse_cell_cols) +
    theme(legend.position = "none")
ggsave(plot = mouse_tSNE_noLegend,
       filename = here("plots","Conservation","mnn_tSNE_Mouse_CellType_Colored_withoutLegend.pdf"),
       height = 10,width = 10)
#No Legend, labelled 
mouse_tSNE_labeled <- plotReducedDim(object = sce_harmony_Species,
                                     dimred = "tSNE_corrected_50",
                                     colour_by = "Mouse_CellType",
                                     text_by = "Mouse_CellType",
                                     point_alpha = 0.4) +
    scale_color_manual(values = mouse_cell_cols) +
    theme(legend.position = "none")
ggsave(plot = mouse_tSNE_labeled,
       filename = here("plots","Conservation","mnn_tSNE_mouse_CellType_Colored_labeled.pdf"),
       height = 13,width = 13)

#30,000 cells so will use higher k values. 
#jaccard + louvain is similar to seurat workflow. 
#Resolution=1 is the default
snn_k_10_mnn <- buildSNNGraph(sce_harmony_Species, k = 10, use.dimred = "PCA_MNN",type="jaccard")
snn_k_25_mnn <- buildSNNGraph(sce_harmony_Species, k = 25, use.dimred = "PCA_MNN",type="jaccard")
snn_k_50_mnn <- buildSNNGraph(sce_harmony_Species, k = 50, use.dimred = "PCA_MNN",type="jaccard")
snn_k_75_mnn<- buildSNNGraph(sce_harmony_Species, k = 75, use.dimred = "PCA_MNN",type="jaccard")
snn_k_100_mnn <- buildSNNGraph(sce_harmony_Species, k = 100, use.dimred = "PCA_MNN",type="jaccard")


#Louvain clustering
#k=10
set.seed(1234)
clust_10_mnn <- igraph::cluster_louvain(snn_k_10_mnn,resolution=1)$membership
table(clust_10_mnn)
# clust_10_mnn
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 2262 2493 1463  564  390  589 1794 1032 1225 1715  560  213  632 2009 1703 1631 
# 17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32 
# 187  481  183 1992 3997  155  180   76  429  397  382  508  724   57   59  623 
# 33   34 
# 109  295

#k=25
set.seed(1234)
clust_25_mnn <- igraph::cluster_louvain(snn_k_25_mnn,resolution=1)$membership
table(clust_25_mnn)
# clust_25_mnn
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 248 2769 1427  614  389 2276 1760  980 1250 1865  397  217 2872  177   78 1642 
# 17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32 
# 1676  237  465  177 1834 4208  227  171  310   75  592  487   57  696  634  302 

#k=50
set.seed(1234)
clust_50_mnn <- igraph::cluster_louvain(snn_k_50_mnn,resolution=1)$membership
table(clust_50_mnn)
# clust_50_mnn
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 254 2848 1509  675  387 2274 1804 1010 1397  395 3120 1895  309 1618 1648  231 
# 17   18   19   20   21   22   23   24   25   26   27   28 
# 467  171 1842 4214   73  492  276  474   57  723  650  296 

#k=75
set.seed(1234)
clust_75_mnn <- igraph::cluster_louvain(snn_k_75_mnn,resolution=1)$membership
table(clust_75_mnn)
# clust_75_mnn
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 2614 2769 1524  673  386 1811 1004 1392  385 3043  594 1828 1595 1767  236  526 
# 17   18   19   20   21   22   23   24   25 
# 173 1866 4198  170  429  473  713  645  295 

#k=100
set.seed(1234)
clust_100_mnn <- igraph::cluster_louvain(snn_k_100_mnn,resolution=1)$membership
table(clust_100_mnn)
clust_100_mnn
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 2625 4770 1523  765  386 1816  998 1407  380 3341  328 1575 1760  218  526  168 
# 17   18   19   20   21   22   23 
# 1876 4210  471  299  704  642  321 

#Add cluster information to object
sce_harmony_Species$k_10_louvain_1_mnn <- factor(clust_10_mnn)
sce_harmony_Species$k_25_louvain_1_mnn <- factor(clust_25_mnn)
sce_harmony_Species$k_50_louvain_1_mnn <- factor(clust_50_mnn)
sce_harmony_Species$k_75_louvain_1_mnn <- factor(clust_75_mnn)
sce_harmony_Species$k_100_louvain_1_mnn <- factor(clust_100_mnn)

#Plot the clusters on the tSNE
for(i in c("k_10_louvain_1_mnn",
           "k_25_louvain_1_mnn",
           "k_50_louvain_1_mnn",
           "k_75_louvain_1_mnn",
           "k_100_louvain_1_mnn")){
    print(i)
    x <- plotReducedDim(object = sce_harmony_Species,
                        dimred = "tSNE_corrected_50",
                        colour_by = i,text_by = i) +
        ggtitle(i) +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(plot = x,filename = here("plots","Conservation","Clustering","mnn",paste0(i,".pdf")))
}

#Plot genes from above on top of the mnn corrected tsne
for(i in genes){
    print(i)
    x <- rowData(sce_harmony_Species)[which(rowData(sce_harmony_Species)$gene_name == i),"gene_id"]
    gene_plot <- plotReducedDim(object = sce_harmony_Species,
                                dimred = "tSNE_corrected_50",
                                color_by = x) +
        scale_color_gradientn(colours = c("lightgrey","orange","red")) +
        ggtitle(paste0(i,"\n",
                       rowData(sce_harmony_Species)[which(rowData(sce_harmony_Species)$gene_name == i),"gene_id"])) +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(plot = gene_plot,filename = here("plots","Conservation",
                                            "Clustering","mnn","Expression",
                                            paste0(i,"_cross_species_mnn_corrected.pdf")))
}


#barplot 
for(i in c("k_10_louvain_1_mnn","k_25_louvain_1_mnn","k_50_louvain_1_mnn",
           "k_75_louvain_1_mnn","k_100_louvain_1_mnn")){
    for(l in unique(levels(colData(sce_harmony_Species)[,i]))){
        rows <- which(colData(sce_harmony_Species)[,i] == l)
        data <- as.data.frame(table(colData(sce_harmony_Species)[rows,"CellType_Species"]))
        data$prop <- data$Freq/sum(data$Freq)*100
        prop_plot <- ggplot(data, aes(x=reorder(Var1,prop), y=prop, fill=Var1)) +
            geom_bar(stat="identity") +
            labs(x = "Cell Type",
                 y = "Proportion") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "none")
        ggsave(plot = prop_plot,
               filename = here("plots","Conservation","Cluster_Proportion_barplots","mnn",i,
                               paste0(l,"_prop_barplot_mnn.pdf")))
        
}}

#rename the object
sce_cs <- sce_harmony_Species
save(sce_cs,file = here("processed-data","sce_cs.rda"))

#Will move forward with the mnn object. 
#k=75 looks best. 
#To aid with cluster annotation, I will 1vALL run DEG testing to pull top 50 genes. 


markers_1vALL_enrich <- findMarkers_1vAll(sce_cs,
                                          assay_name = "logcounts",
                                          cellType_col = "k_75_louvain_1_mnn",
                                          mod = "~Sample")
# 1 - '2024-01-10 16:00:15.54975
# 2 - '2024-01-10 16:01:19.797639
# 3 - '2024-01-10 16:02:23.971816
# 4 - '2024-01-10 16:03:28.01513
# 5 - '2024-01-10 16:04:32.727875
# 6 - '2024-01-10 16:05:36.879443
# 7 - '2024-01-10 16:06:40.161971
# 8 - '2024-01-10 16:07:44.528281
# 9 - '2024-01-10 16:08:48.210028
# 10 - '2024-01-10 16:09:52.699159
# 11 - '2024-01-10 16:10:56.554436
# 12 - '2024-01-10 16:12:00.625374
# 13 - '2024-01-10 16:13:04.777749
# 14 - '2024-01-10 16:14:08.734222
# 15 - '2024-01-10 16:15:12.800999
# 16 - '2024-01-10 16:16:17.302194
# 17 - '2024-01-10 16:17:21.370131
# 18 - '2024-01-10 16:18:25.728023
# 19 - '2024-01-10 16:19:25.61338
# 20 - '2024-01-10 16:20:16.207099
# 21 - '2024-01-10 16:21:06.762948
# 22 - '2024-01-10 16:21:56.038117
# 23 - '2024-01-10 16:22:39.692608
# 24 - '2024-01-10 16:23:23.459407
# 25 - '2024-01-10 16:24:07.644038
# Building Table - 2024-01-10 16:24:51.917303
# ** Done! **

nrow(markers_1vALL_enrich)
#[1] 414050

#Add gene name information. 
markers_gene_name <- merge(x    = markers_1vALL_enrich,
                           y    = rowData(sce_cs)[,c("gene_id","gene_name")],
                           by.x = "gene",
                           by.y = "gene_id")
nrow(markers_gene_name)
#[1] 414050

#Subset for top 50 markers for each cluster
markers_gene_name <- subset(markers_gene_name,subset=(rank_marker %in% 1:50))

#Order by cellType.target
markers_gene_name <- markers_gene_name[order(markers_gene_name$cellType.target,decreasing = FALSE),]

#save the dataframe
save(markers_gene_name,file = here("processed-data","integrated_unannotated_1vALL_DEGs.rda"))

#Begin annotation
# 1 - LS_Inh
# 2 - Drd1_MSN
# 3 - Excit_A 
# 4 - Polydendrocyte
# 5 - Microglia
# 6 - Oligo_Human
# 7 - Ependymal
# 8 - Drd1_MSN_Patch
# 9 - MS_Inh
# 10 - Sept_Inh_A
# 11 - Sept_Inh_B
# 12 - Drd2_MSN
# 13 - Sept_Str_Inh
# 14 - Oligo_Mouse
# 15 - MS_Excit_Human
# 16 - Mural
# 17 - Astrocyte_Human
# 18 - Astrocyte_Mouse_A
# 19 - Astrocyte_Mouse_B
# 20 - Drd1_MSN_Matrix_Human
# 21 - Sept_Mouse
# 22 - TNoS_Mouse
# 23 - Neuroblast_Mouse
# 24 - Thal_Mouse
# 25 - IoC_Mouse


