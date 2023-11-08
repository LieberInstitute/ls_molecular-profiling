#Goal: Identify marker genes for each cluster and annotate the clusters. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

library(SingleCellExperiment)
library(DeconvoBuddies)
library(sessioninfo)
library(Polychrome)
library(ggplot2)
library(scater)
library(scran)
library(here)

#load the clustered object
load(here("processed-data","sce_clustered.rda"))

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
# colData names(59): Sample Barcode ... k_50_walktrap k_75_walktrap
# reducedDimNames(18): GLMPCA_approx UMAP_15 ... tSNE_mnn_25 tSNE_mnn_50
# mainExpName: NULL
# altExpNames(0):

##################################################
############ Explore the clusters ################
##################################################
#Dotplot with major markers of LS, MS, Septal, Str, + major oligo clusters. 
genes <- c("SYT1","SNAP25", #Pan-neuronal
           "GAD1","GAD2","SLC32A1", #GABAergic
           "SLC17A6","SLC17A7","SLC17A8", #Glutamatergic
           "TRPC4","DGKG","HOMER2","PTPN3","TRHDE","CPNE7","NRP1", #Lateral Septum markers from mouse
           "ELAVL2","TRPC5", #Medial Septum markers from mouse
           "RARB","BCL11B","PPP1R1B", #Broad striatal + MSN markers 
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

#Louvain clustering
Expression_dotplot <- plotDots(object = sce,
                               features = rev(genes),
                               group = "k_20_louvain_1",
                               swap_rownames = "gene_name") +
    scale_color_gradientn(colours = c("lightgrey","orange","red"))
ggsave(plot = Expression_dotplot,filename = here("plots","Expression_plots",
                                                 "post_k_20_louvain_1_clustering_general_dotplot_v2.pdf"),
       height = 8)


#Plot these genes on top of the tSNE
#Add several additional markers of MSNs. 
for(i in c(genes,"FOXP2","PDE1B","KIAA1211L","PDE2A","SLIT3","NGEF")){
    print(i)
    x <- plotReducedDim(object = sce,
                        dimred = "tSNE_mnn_15",
                        colour_by = i,
                        swap_rownames = "gene_name") +
        scale_color_gradientn(colours = c("lightgrey","orange","red"))
    ggsave(plot = x,
           filename = here("plots","Expression_plots","post_k_20_louvain_clustering",paste0(i,"_tSNE_mnn_15.pdf")))
}

##################################################
############   Annotate clusters      ############    
##################################################
#Create annotation dataframe
annotation_df <- data.frame(cluster= 1:24,
                            celltype = c("Sept_Inh_A","Str_Inh_A","Excit_A","Polydendrocyte",
                                         "Microglia","Sept_Inh_B","Oligo_A","Ependymal",
                                         "Oligo_B","Str_Inh_B","Sept_Inh_C","Sept_Inh_D",
                                         "Sept_Inh_E","Str_Inh_C","Sept_Inh_F","Sept_Excit_A",
                                         "Excit_B","Sept_Inh_G","Mural","Astrocyte",
                                         "Sept_Inh_H","Sept_Inh_I","Oligo_C","Str_Inh_D"))

#add celltype info
sce$CellType_k_20_louvain <- annotation_df$celltype[match(sce$k_20_louvain_1,
                                                          annotation_df$cluster)]
#factorize. 
sce$CellType_k_20_louvain <- factor(sce$CellType_k_20_louvain,
                                    levels = c("Sept_Inh_A","Sept_Inh_B","Sept_Inh_C","Sept_Inh_D",
                                               "Sept_Inh_E","Sept_Inh_F","Sept_Inh_G","Sept_Inh_H",
                                               "Sept_Inh_I","Sept_Excit_A","Str_Inh_A","Str_Inh_B",
                                               "Str_Inh_C","Str_Inh_D","Excit_A","Excit_B",
                                               "Oligo_A","Oligo_B","Oligo_C","Polydendrocyte",
                                               "Astrocyte","Ependymal","Microglia","Mural"))

#Annotate the tSNE
cluster_cols <- Polychrome::createPalette(length(unique(sce$CellType_k_20_louvain)),c("#FF0000", "#00FF00", "#0000FF"))
names(cluster_cols) <- unique(sce$CellType_k_20_louvain)

annotated_tSNE <- plotReducedDim(object = sce,
                                 dimred = "tSNE_mnn_15",
                                 colour_by = "CellType_k_20_louvain",
                                 text_by = "CellType_k_20_louvain") +
    scale_color_manual(values = cluster_cols)
ggsave(filename = here("plots","Dim_Red","tSNE_mnn_15_annotated_CellType_k20_louvain.pdf"),
       plot = annotated_tSNE)

########Calculate modularity scores.
set.seed(1234)
#Make the graph. 
snn_k_20 <- buildSNNGraph(sce, k = 20, use.dimred = "mnn",type="jaccard")

#Calcualte modularity scores. 
k_20_modularity <- bluster::pairwiseModularity(graph = snn_k_20,
                                               clusters = sce$CellType_k_20_louvain,
                                               as.ratio = TRUE)

#make a heatmap. 
library(pheatmap)
pdf(file = here("plots","CellType_k_20_louvain_pairwise_modularity.pdf"))
pheatmap(log2(k_20_modularity+1), 
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         display_numbers=TRUE, 
         number_format="%.2f", 
         fontsize_number=6.5,
         main = "Modularity ratio for 24 clusters in human LS (n=3)",
         color=colorRampPalette(c("white","orange","red"))(100))
dev.off()

##################################################
#########Prep sce object for DEG testing##########
##################################################

#Do any genes have 0 counts for every cell. 
table(rowSums(assay(sce, "counts")) == 0)
# FALSE  TRUE 
# 33556  3045 

#Remove 3045 genes that are all 0s 
sce <- sce[!rowSums(assay(sce, "counts")) == 0, ]

##################################################
###############run 1 vs all testing###############
##################################################
markers_1vALL_enrich <- findMarkers_1vAll(sce, 
                                          assay_name = "logcounts", 
                                          cellType_col = "CellType_k_20_louvain", 
                                          mod = "~Sample")
# Sept_Inh_A - '2023-11-08 09:41:12.810285
# Str_Inh_A - '2023-11-08 09:41:29.786506
# Excit_A - '2023-11-08 09:41:46.657498
# Polydendrocyte - '2023-11-08 09:42:03.591051
# Microglia - '2023-11-08 09:42:20.535435
# Sept_Inh_B - '2023-11-08 09:42:37.437445
# Oligo_A - '2023-11-08 09:42:54.335511
# Ependymal - '2023-11-08 09:43:11.236241
# Oligo_B - '2023-11-08 09:43:28.417719
# Str_Inh_B - '2023-11-08 09:43:45.380726
# Sept_Inh_C - '2023-11-08 09:44:02.623924
# Sept_Inh_D - '2023-11-08 09:44:19.597667
# Sept_Inh_E - '2023-11-08 09:44:36.583435
# Str_Inh_C - '2023-11-08 09:44:53.842077
# Sept_Inh_F - '2023-11-08 09:45:11.948174
# Sept_Excit_A - '2023-11-08 09:45:29.059849
# Excit_B - '2023-11-08 09:45:45.996569
# Sept_Inh_G - '2023-11-08 09:46:02.548759
# Mural - '2023-11-08 09:46:19.526886
# Astrocyte - '2023-11-08 09:46:36.846822
# Sept_Inh_H - '2023-11-08 09:46:53.798798
# Sept_Inh_I - '2023-11-08 09:47:10.814651
# Oligo_C - '2023-11-08 09:47:27.890011
# Str_Inh_D - '2023-11-08 09:47:44.854139
# Building Table - 2023-11-08 09:48:02.068503
# ** Done! **

#Add symbol information to the table
#First change the ensembl gene id column to have same name as what is in rowData(sce)
colnames(markers_1vALL_enrich)[1] <- "gene_id"
markers_1vALL_df <- dplyr::left_join(x = as.data.frame(markers_1vALL_enrich),
                                     y = as.data.frame(rowData(sce)[,c("gene_id","gene_name")]),
                                     by = "gene_id")

#save the dataframe. 
save(markers_1vALL_df,file = here("processed-data","markers_1vAll_ttest_k_20_louvain_24Clusters.rda"))
###############################

#Modularity scores suggest Str_Inh_A, Str_Inh_B, and Str_Inh_C are related. 
#Check out their gene markers and see if their is significant overlap. 

#DEGs: logFC>=0.5 & FDR<=0.001
str_a <- subset(markers_1vALL_df,subset=(cellType.target == "Str_Inh_A" & logFC >= 0.5 & log.FDR < 0.001))
str_c <- subset(markers_1vALL_df,subset=(cellType.target == "Str_Inh_C" & logFC >= 0.5 & log.FDR < 0.001))
str_d <- subset(markers_1vALL_df,subset=(cellType.target == "Str_Inh_D" & logFC >= 0.5 & log.FDR < 0.001))

dim(str_a)
# [1] 964   9
dim(str_c)
# [1] 1011    9
dim(str_d)
# [1] 1321    9

#Check how many genes intersect
length(intersect(str_a$gene_name,str_c$gene_name))
# [1] 763
length(intersect(str_a$gene_name,str_d$gene_name))
# [1] 812
length(intersect(str_c$gene_name,str_d$gene_name))
# [1] 932
#61-92% of DEGs are shared between these populations. Will merge these. 
###############################

##################################################
############  Begin final annotation  ############    
##################################################
#Create final celltype 
sce$CellType.Final <- as.character(sce$CellType_k_20_louvain)

#Combine str_a,c,d and then combine the three oligo populations. 
sce$CellType.Final[sce$CellType.Final %in% c("Str_Inh_A", "Str_Inh_C","Str_Inh_D")] <- "Str_Inh_A"
sce$CellType.Final[sce$CellType.Final %in% c("Oligo_A", "Oligo_B","Oligo_C")] <- "Oligo"

#Check out the tSNE
cluster_cols <- Polychrome::createPalette(length(unique(sce$CellType.Final)),c("#FF0000", "#00FF00", "#0000FF"))
names(cluster_cols) <- unique(sce$CellType.Final)
tSNE_celltype_final <- plotReducedDim(object = sce,
                                      dimred = "tSNE_mnn_15",
                                      colour_by = "CellType.Final",
                                      text_by = "CellType.Final")
ggsave(plot = tSNE_celltype_final,
       filename = here("plots","Dim_Red","tSNE_mnn_15_20Clusters.pdf"))

###Begin annotation of Septal clusters. 
#subset sce for just the septal clusters to investigate. 
sce_sept <- sce[,grep("Sept",sce$CellType.Final)]

#Plot tSNE to ensure the subset worked. 
tSNE_sept_only <- plotReducedDim(object = sce_sept,
                                 dimred = "tSNE_mnn_1f5",
                                 colour_by = "CellType.Final",
                                 text_by = "CellType.Final")
ggsave(plot = tSNE_sept_only,
       filename = here("plots","Dim_Red","tSNE_mnn_15_SeptOnly.pdf"))


#Septal violin plots
Sept_violin <- plotExpression(object = sce_sept,features = c("TRPC4","HOMER2","PTPN3","TRHDE", #LS markers from mouse
                                                             "DGKG", #boad septal mouse from mouse
                                                             "ELAVL2","TRPC5",#MS markers from mouse
                                                             "SLC17A6","CHAT",#MS subsets (Excit + Chat-int.) 
                                                             "CRHR2", "OPRM1","CCK","VIP","KIT",#Expressed in LS NEURONS
                                                             "FRMD3","UNC5D","ABLIM3",#BROAD SEPTAL MARKERS from mouse
                                                             "NOS1","SST",
                                                             "BCL11B","PPP1R1B","ISL1","RARB","DRD1"), 
                              swap_rownames = "gene_name",x = "CellType.Final",
                              colour_by = "CellType.Final",ncol = 4) +
    theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")
ggsave(plot = Sept_violin,
       filename = here("plots","Expression_plots","Sept_Str_markers_SeptOnly.pdf"),
       width = 10,height =12)

#Plot RARB and LHX6
RARB_violin <- plotExpression(object = sce,features = c("RARB"), 
                              swap_rownames = "gene_name",x = "CellType.Final",
                              colour_by = "CellType.Final",ncol = 4) +
    theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")
ggsave(plot = RARB_violin,
       filename = here("plots","Expression_plots","RARB_VIOLIN.pdf"))

LHX6_violin <- plotExpression(object = sce,features = c("LHX6"), 
                              swap_rownames = "gene_name",x = "CellType.Final",
                              colour_by = "CellType.Final",ncol = 4) +
    theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")
ggsave(plot = LHX6_violin,
       filename = here("plots","Expression_plots","LHX6_VIOLIN.pdf"))

#Make final Celltype designations. 



#Save the object with celltype
save(sce,file = here("processed-data","sce_with_CellType.rda"))

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
