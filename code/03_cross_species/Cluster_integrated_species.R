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


###################################
############# HARMONY #############
#First cluster and evaluate harmony integration. 

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

#Plot the clusters on the tSNE with 15 dimensions. 
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


#Create annotation dataframes
annotation_df <- data.frame(cluster= 1:23,
                            celltype = c("LS_Inh_A","Drd1_MSN_A","Excit_A","Polydendrocyte",
                                         "Microglia","Oligodendrocyte","Ependymal","Drd1_MSN_Patch",
                                         "Drd2_MSN","MS_Inh_A","Septal_Inh_A","Septal_Inh_B",
                                         "Str_Inh_A","Excit_B","Mural","Astrocyte_1",
                                         "Astrocyte_2","Neuroblast","Astrocyte_3","Septal_Inh_C",
                                         "IoC","TNoS","Thal"))

#add celltype info
sce_harmony_Species$CellType_k_75_louvain <- annotation_df$celltype[match(sce_harmony_Species$k_75_louvain_1,
                                                                          annotation_df$cluster)]


#Make celltype a factor
sce_harmony_Species$CellType_k_75_louvain <- factor(sce_harmony_Species$CellType_k_75_louvain,
                                                    levels = c("LS_Inh_A","MS_Inh_A","Septal_Inh_A","Septal_Inh_B",
                                                               "Septal_Inh_C","Drd1_MSN_A","Drd1_MSN_Patch","Drd2_MSN",
                                                               "Str_Inh_A","Excit_A","Excit_B","IoC","TNoS","Thal",
                                                               "Astrocyte_1","Astrocyte_2","Astrocyte_3","Ependymal",
                                                               "Oligodendrocyte","Polydendrocyte","Microglia",
                                                               "Mural","Neuroblast"))

#Annotate the tSNE
cluster_cols <- Polychrome::createPalette(length(unique(sce_harmony_Species$CellType_k_75_louvain)),
                                          c("#D81B60", "#1E88E5","#FFC107","#009E73"))
names(cluster_cols) <- unique(sce_harmony_Species$CellType_k_75_louvain)

#Save the cluster cols
save(cluster_cols,file = here("plots","Conservation","harmony_integration_cluster_cols.rda"))

#make plot
annotated_tSNE <- plotReducedDim(object = sce_harmony_Species,
                                 dimred = "tSNE_HARMONY",
                                 colour_by = "CellType_k_75_louvain",
                                 text_by = "CellType_k_75_louvain") +
    scale_color_manual(values = cluster_cols) +
    ggtitle("HARMONY Integration") +
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5)) 
ggsave(filename = here("plots","Conservation","tSNE_Harmony_annotatedCellTypes.pdf"),
       plot = annotated_tSNE)

#Calculate LISI to assess integration. 
#Pull tSNE embeddings from HARMONY reduced dim
harmony_tSNE_embeds <- as.data.frame(reducedDim(sce_harmony_Species,"tSNE_HARMONY"))

#DF of metadata
harmony_metadata <- as.data.frame(colData(sce_harmony_Species)[,c("Species","CellType_Species","CellType_k_75_louvain")])

#Run lisi
res <- compute_lisi(harmony_tSNE_embeds,harmony_metadata,"Species")
colnames(res)[1] <- "Species_LISI"

#Add celltype to res
res$CellID <- rownames(res)
harmony_metadata$CellID <- rownames(harmony_metadata)
res_all <- merge(x = res,
                 y = harmony_metadata,
                 by = "CellID")

#ridge plot of LISI by cluster. Peaks closer to 2 are more well integrated. 
harmony_lisi_ridge <- ggplot(res_all,aes(x = Species_LISI,y = CellType_k_75_louvain,fill = CellType_k_75_louvain)) +
    geom_density_ridges() +
    scale_fill_manual(values = cluster_cols) +
    theme(legend.position = "none")

ggsave(plot = harmony_lisi_ridge,filename = here("plots","Conservation","harmony_lisi_ridge.pdf"))

#Many of the distributions are centered on 1. Most likely due to the fact that there are double 
#the number of mouse to human cells so cells within neighborhood are mostly going to be mouse. 

#Calculate the proportion of each species that are within each cluster. 
#human
human_only <- subset(colData(sce_harmony_Species),subset=(Species == "Human"))
human_props <- as.data.frame((table(human_only$CellType_k_75_louvain)/nrow(human_only))*100)
colnames(human_props) <- c("CellType","Human")
#mouse
mouse_only <- subset(colData(sce_harmony_Species),subset=(Species == "Mouse"))
mouse_props <- as.data.frame((table(mouse_only$CellType_k_75_louvain)/nrow(mouse_only))*100)
colnames(mouse_props)  <- c("CellType","Mouse")

#Combine human and mouse props
merged_props <- merge(x = human_props,
                      y = mouse_props,
                      by = "CellType")
#Check that everything was calculated correctly. 
sum(merged_props$Human)

sum(merged_props$Mouse)


#Melt merged_props
merged_props_melt <- reshape2::melt(merged_props)

#barplot with proportions
prop_barplot <- ggplot(data = merged_props_melt,aes(x = reorder(CellType,value),y = value,fill = variable)) +
    geom_bar(stat = "identity",position="dodge") +
    theme_bw() +
    labs(x    = "Cell Type",
         y    = "Percentage of Species Total Cells in Cluster",
         fill = "Species",
         title = "HARMONY Integration") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)) 
ggsave(filename = here("plots","Conservation","proportion_barplot.pdf"),plot = prop_barplot)

#Finalize DEG analysis with cluster names. 
Harmony_cellType_1vALL <- findMarkers_1vAll(sce_harmony_Species, 
                                            assay_name = "logcounts", 
                                            cellType_col = "CellType_k_75_louvain", 
                                            mod = "~Sample")


save(Harmony_cellType_1vALL,file = here("processed-data","Harmony_cellType_1vALL.rda"))

###################################
############### MNN ###############
###################################
#Cluster the mnn integration. 
#PCA_Corrected is the 




