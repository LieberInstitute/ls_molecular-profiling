#Goal: Identify marker genes for each cluster and annotate the clusters. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

library(SingleCellExperiment)
library(DeconvoBuddies)
library(sessioninfo)
library(scater)
library(ggplot2)
library(scran)
library(here)

#load the clustered object
load(here("processed-data","sce_clustered.rda"))

#Going to use k=50 walktrap clusters for cluster identification. 
table(sce$k_50_walktrap)
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 544  344 1403  361 1357  643  197 1399  630  272  254 1093  355  224  149 

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
                                          cellType_col = "k_50_walktrap", 
                                          mod = "~Sample")
# 8 - '2023-10-23 12:52:49.400386
# 3 - '2023-10-23 12:53:09.2398
# 4 - '2023-10-23 12:53:28.766212
# 7 - '2023-10-23 12:53:48.195631
# 11 - '2023-10-23 12:54:07.929927
# 12 - '2023-10-23 12:54:27.729303
# 1 - '2023-10-23 12:54:47.466659
# 9 - '2023-10-23 12:55:07.404354
# 6 - '2023-10-23 12:55:27.113037
# 14 - '2023-10-23 12:55:46.896635
# 5 - '2023-10-23 12:56:06.688447
# 15 - '2023-10-23 12:56:26.534018
# 2 - '2023-10-23 12:56:46.312353
# 10 - '2023-10-23 12:57:06.135446
# 13 - '2023-10-23 12:57:25.947295
# Building Table - 2023-10-23 12:57:45.824007
# ** Done! **

#Add symbol information to the table
#First change the ensemble gene id column to have same name as what is in rowData(sce)
colnames(markers_1vALL_enrich)[1] <- "gene_id"
markers_1vALL_df <- dplyr::left_join(x = as.data.frame(markers_1vALL_enrich),
                                     y = as.data.frame(rowData(sce)[,c("gene_id","gene_name")]),
                                     by = "gene_id")

#save the dataframe. 
save(markers_1vALL_df,file = here("processed-data","markers_1vAll_ttest.rda"))


##################################################
############Annotate the clusters.################
##################################################
annotation_df <- data.frame(cluster=c(1:15),
                            celltype = c("Astrocyte_1","Astrocyte_2","LS_GABA_2","Glutamatergic","MS_GABA_1",
                                         "Striosome","Polydendrocyte","LS_GABA_1","Oligo_1","Astrocyte_3",
                                         "Microglia","Oligo_2","Oligo_3","MS_GABA_2","Endothelial"))



sce$CellType <- annotation_df$celltype[match(sce$k_50_walktrap,
                                             annotation_df$cluster)]
sce$CellType <- factor(sce$CellType,
                       levels = c("LS_GABA_1","LS_GABA_2",
                                  "MS_GABA_1","MS_GABA_2",
                                  "Glutamatergic","Striosome",
                                  "Astrocyte_1","Astrocyte_2","Astrocyte_3",
                                  "Oligo_1","Oligo_2","Oligo_3",
                                  "Microglia",
                                  "Polydendrocyte",
                                  "Endothelial"))

###check expression profiles of the clusters to convince youself that these are correct. 
genes <- c("SYT1","SNAP25","GAD1","GAD2","SLC32A1",
           "TRPC4","HOMER2","PTPN3",
           "ELAVL2",
           "SLC17A7", "SLC17A6", "SLC17A8",
           "DRD1","OPRM1","FOXP2",
           "GFAP", "TNC", "AQP4", "SLC1A2",
           "MBP","MOBP",
           "CD74", "CSF1R", "C3",
           "PDGFRA", "VCAN", "CSPG4",
           "CLDN5", "FLT1", "VTN")


Expression_dotplot <- plotDots(object = sce,
                               features = rev(genes),
                               group = "CellType",swap_rownames = "gene_name") +
    scale_color_gradientn(colours = c("lightgrey","orange","red")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = Expression_dotplot,filename = here("plots","Expression_plots",
                                                 "annotated_dotplot.pdf"),
       height = 8)

###Annotate the umap. 
annotated_umap <- plotReducedDim(object = sce,dimred = "UMAP_mnn_15",colour_by = "CellType",text_by = "CellType")
ggsave(plot = annotated_umap,filename = here("plots","Dim_Red",
                                                 "annotated_umap.pdf"))

