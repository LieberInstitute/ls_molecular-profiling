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
for(i in c(genes,
           "RXFP1","EBF1","PDYN", "CHST9","SEMA5B","TAC1","STXBP6",#Additional D1 markers including D1 islands. 
           "DRD2","PENK","ADORA2A", #D2 markers. 
           "CRYM",#Medial dorsal Striatum marker
           "DLK1",#ventral medial dorsal striaum 
           "FOXP2","PDE1B","KIAA1211L","PDE2A","SLIT3","NGEF")){
    print(i)
    x <- plotReducedDim(object = sce,
                        dimred = "tSNE_mnn_50",
                        colour_by = i,
                        swap_rownames = "gene_name") +
        scale_color_gradientn(colours = c("lightgrey","orange","red"))
    ggsave(plot = x,
           filename = here("plots","Expression_plots","post_k_20_louvain_clustering",paste0(i,"_tSNE_mnn_50.pdf")))
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
cluster_cols <- Polychrome::createPalette(length(unique(sce$CellType_k_20_louvain)),
                                          c("#D81B60", "#1E88E5","#FFC107","#009E73"))
names(cluster_cols) <- unique(sce$CellType_k_20_louvain)

annotated_tSNE <- plotReducedDim(object = sce,
                                 dimred = "tSNE_mnn_50",
                                 colour_by = "CellType_k_20_louvain",
                                 text_by = "CellType_k_20_louvain") +
    scale_color_manual(values = cluster_cols) +
    theme(legend.position = "none")
ggsave(filename = here("plots","Dim_Red","tSNE_mnn_50_annotated_CellType_k20_louvain.pdf"),
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
# Sept_Inh_A - '2023-12-08 09:34:03.357655
# Str_Inh_A - '2023-12-08 09:34:31.410573
# Excit_A - '2023-12-08 09:34:58.293269
# Polydendrocyte - '2023-12-08 09:35:25.292636
# Microglia - '2023-12-08 09:35:52.505842
# Sept_Inh_B - '2023-12-08 09:36:19.828381
# Oligo_A - '2023-12-08 09:36:46.716865
# Ependymal - '2023-12-08 09:37:13.470004
# Oligo_B - '2023-12-08 09:37:39.993648
# Str_Inh_B - '2023-12-08 09:38:07.057986
# Sept_Inh_C - '2023-12-08 09:38:33.341049
# Sept_Inh_D - '2023-12-08 09:38:59.867286
# Sept_Inh_E - '2023-12-08 09:39:26.45466
# Str_Inh_C - '2023-12-08 09:39:53.230684
# Sept_Inh_F - '2023-12-08 09:40:20.98131
# Sept_Excit_A - '2023-12-08 09:40:47.121283
# Excit_B - '2023-12-08 09:41:13.534617
# Sept_Inh_G - '2023-12-08 09:41:39.272678
# Mural - '2023-12-08 09:42:05.578006
# Astrocyte - '2023-12-08 09:42:32.089838
# Sept_Inh_H - '2023-12-08 09:42:59.001963
# Sept_Inh_I - '2023-12-08 09:43:25.420513
# Oligo_C - '2023-12-08 09:43:51.70494
# Str_Inh_D - '2023-12-08 09:44:18.039134
# Building Table - 2023-12-08 09:44:44.270452
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

#subset for just striatal populations and plot some Drd1 and Drd2 marker genes. 
#str_a,c, and d might not be separated well because they represent D1/D2 MSNs. 
#These cells do have a ton of CRYM which is a known marker of the medial striatum. 
str_vln <- plotExpression(object = sce[,sce$CellType_k_20_louvain %in% c("Str_Inh_A","Str_Inh_B","Str_Inh_C","Str_Inh_D")],
               features = c("DRD1","TAC1","PDYN",#D1-MSN markers. 
                            "OPRM1","SEMA5B", #Patch
                            "CRYM", #MEDIAL DORSAL STRIATUM
                            "DLK1", #VENTRAL MEDIAL DORSAL STRIAUM
                            "STXBP6","RELN",
                            "DRD2","ADORA2A","PENK",#Drd2 markers
                            "RARB","PPP1R1B","BCL11B", #Striatum/MSN
                            "TRPC4","DGKG","FREM2"), 
               x = "CellType_k_20_louvain",
               colour_by = "CellType_k_20_louvain",
               ncol = 3,
               swap_rownames = "gene_name") +
    theme(axis.text.x = element_text(angle = 90,vjust = 1),legend.position = "none")
ggsave(plot = str_vln,
       filename = here("plots","Expression_plots","Striatal_markers_violin_plot.pdf"),
       height = 12,width = 8)

#Str_A Str_A-Drd1-MSN
#Str_B Str_B_Drd1-Patch
#Str_C Str_C_Drd2-MSN
#Str_D Str_D_Drd1-Matrix
###############################

##################################################
############  Begin final annotation  ############    
##################################################
#Create final celltype 
sce$CellType.Final <- as.character(sce$CellType_k_20_louvain)

#Combine str_a,c,d and then combine the three oligo populations. 
sce$CellType.Final[sce$CellType.Final %in% c("Oligo_A", "Oligo_B","Oligo_C")] <- "Oligo"

#Check out the tSNE
cluster_cols <- Polychrome::createPalette(length(unique(sce$CellType.Final)),
                                          c("#D81B60", "#1E88E5","#FFC107","#009E73"))
names(cluster_cols) <- unique(sce$CellType.Final)
tSNE_celltype_final <- plotReducedDim(object = sce,
                                      dimred = "tSNE_mnn_50",
                                      colour_by = "CellType.Final",
                                      text_by = "CellType.Final") +
    scale_color_manual(values = cluster_cols) +
    theme(legend.position = "none")
ggsave(plot = tSNE_celltype_final,
       filename = here("plots","Dim_Red","tSNE_mnn_50_22Clusters.pdf"))

###Begin annotation of Septal clusters. 
#subset sce for just the septal clusters to investigate. 
sce_sept <- sce[,grep("Sept",sce$CellType.Final)]

#Plot tSNE to ensure the subset worked. 
tSNE_sept_only <- plotReducedDim(object = sce_sept,
                                 dimred = "tSNE_mnn_50",
                                 colour_by = "CellType.Final",
                                 text_by = "CellType.Final")
ggsave(plot = tSNE_sept_only,
       filename = here("plots","Dim_Red","tSNE_mnn_50_SeptOnly.pdf"))


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
sce$CellType.Final[grep("Sept_Excit_A", sce$CellType.Final)] <- "MS_Excit_A"
sce$CellType.Final[grep("Sept_Inh_A", sce$CellType.Final)] <- "LS_Inh_A"
sce$CellType.Final[grep("Sept_Inh_B", sce$CellType.Final)] <- "LS_Inh_B"
sce$CellType.Final[grep("Sept_Inh_C", sce$CellType.Final)] <- "MS_Inh_A"
sce$CellType.Final[grep("Sept_Inh_E", sce$CellType.Final)] <- "MS_Inh_E"
sce$CellType.Final[grep("Sept_Inh_G", sce$CellType.Final)] <- "LS_Inh_G"
sce$CellType.Final[grep("Sept_Inh_H", sce$CellType.Final)] <- "MS_Inh_H"
sce$CellType.Final[grep("Sept_Inh_I", sce$CellType.Final)] <- "LS_Inh_I"
sce$CellType.Final[grep("Str_Inh_A", sce$CellType.Final)] <- "Str_Drd1-MSN"
sce$CellType.Final[grep("Str_Inh_B", sce$CellType.Final)] <- "Str_Drd1-Patch"
sce$CellType.Final[grep("Str_Inh_C", sce$CellType.Final)] <- "Str_Drd2-MSN"
sce$CellType.Final[grep("Str_Inh_D", sce$CellType.Final)] <- "Str_Drd1-Matrix"

#Plot the tSNE with final celltypes annotated. 
cluster_cols <- Polychrome::createPalette(length(unique(sce$CellType.Final)),
                                          c("#D81B60", "#1E88E5","#FFC107","#009E73"))
names(cluster_cols) <- unique(sce$CellType.Final)
tSNE_final <- plotReducedDim(object = sce,
                             dimred = "tSNE_mnn_50",
                             colour_by = "CellType.Final",
                             text_by = "CellType.Final") +
    scale_color_manual(values = cluster_cols) +
    theme(legend.position = "none")
ggsave(plot = tSNE_final,
       filename = here("plots","Dim_Red","tSNE_mnn_50_Final_CellTypes.pdf"))

#Save the cluster_cols vector as well (
save(cluster_cols,file = here("processed-data","Final_CellTypes_colors.rda"))

#Save the object with the final celltypes. 
save(sce,file = here("processed-data","sce_with_CellType.rda"))


#Rerun DEG testing with new celltype designations. 
##################################################
##############rerun 1 vs all testing##############
##################################################
markers_1vALL_enrich_Final <- findMarkers_1vAll(sce, 
                                                assay_name = "logcounts", 
                                                cellType_col = "CellType.Final", 
                                                mod = "~Sample")
# LS_Inh_A - '2023-12-11 14:25:42.758302
# Str_Drd1-MSN - '2023-12-11 14:25:53.132007
# Excit_A - '2023-12-11 14:26:03.387258
# Polydendrocyte - '2023-12-11 14:26:13.318529
# Microglia - '2023-12-11 14:26:23.600446
# LS_Inh_B - '2023-12-11 14:26:33.773973
# Oligo - '2023-12-11 14:26:43.772765
# Ependymal - '2023-12-11 14:26:54.009164
# Str_Drd1-Patch - '2023-12-11 14:27:04.306281
# MS_Inh_A - '2023-12-11 14:27:14.263405
# Sept_Inh_D - '2023-12-11 14:27:24.40818
# MS_Inh_E - '2023-12-11 14:27:34.763483
# Str_Drd2-MSN - '2023-12-11 14:27:45.118414
# Sept_Inh_F - '2023-12-11 14:27:55.49254
# MS_Excit_A - '2023-12-11 14:28:05.908915
# Excit_B - '2023-12-11 14:28:16.334723
# LS_Inh_G - '2023-12-11 14:28:26.436024
# Mural - '2023-12-11 14:28:36.47635
# Astrocyte - '2023-12-11 14:28:46.535918
# MS_Inh_H - '2023-12-11 14:28:56.790901
# LS_Inh_I - '2023-12-11 14:29:07.011881
# Str_Drd1-Matrix - '2023-12-11 14:29:17.386279
# Building Table - 2023-12-11 14:29:27.993341
# ** Done! **


#Add symbol information to the table
#First change the ensembl gene id column to have same name as what is in rowData(sce)
colnames(markers_1vALL_enrich_Final)[1] <- "gene_id"
markers_1vALL_enrich_Final <- dplyr::left_join(x = as.data.frame(markers_1vALL_enrich_Final),
                                               y = as.data.frame(rowData(sce)[,c("gene_id","gene_name")]),
                                               by = "gene_id")

#save the dataframe. 
save(markers_1vALL_enrich_Final,file = here("processed-data","markers_1vAll_ttest_CellTypeFinal_22Clusters.rda"))

##################################################
###############run pairwise testing###############
##################################################
mod <- with(colData(sce), model.matrix(~ Sample))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

# Run pairwise t-tests
markers_pairwise <- findMarkers(sce, 
                                groups=sce$CellType.Final,
                                assay.type="logcounts", 
                                design=mod, 
                                test="t",
                                direction="up", 
                                pval.type="all", 
                                full.stats=T)

#How many DEGs for each cluster? 
sapply(markers_pairwise, function(x){table(x$FDR<0.05)})
# Astrocyte Ependymal Excit_A Excit_B LS_Inh_A LS_Inh_B LS_Inh_G LS_Inh_I
# FALSE     36298     34926   36352   36277    36557    36590    36432    36569
# TRUE        303      1675     249     324       44       11      169       32
# Microglia MS_Excit_A MS_Inh_A MS_Inh_E MS_Inh_H Mural Oligo
# FALSE     35684      36460    36431    36414    36456 35617 36046
# TRUE        917        141      170      187      145   984   555
# Polydendrocyte Sept_Inh_D Sept_Inh_F Str_Drd1-Matrix Str_Drd1-MSN
# FALSE          36254      36597      36451           36335        36548
# TRUE             347          4        150             266           53
# Str_Drd1-Patch Str_Drd2-MSN
# FALSE          36534        36561
# TRUE              67           40

#Add gene info to each list.
for(i in names(markers_pairwise)){
    markers_pairwise[[i]] <- as.data.frame(markers_pairwise[[i]])
    markers_pairwise[[i]]$gene_id <- row.names(markers_pairwise[[i]])
    markers_pairwise[[i]] <- dplyr::left_join(x  =  markers_pairwise[[i]],
                                              y  =  as.data.frame(rowData(sce)[,c("gene_id","gene_name")]),
                                              by = "gene_id")
}

save(markers_pairwise,file = here("processed-data","markers_pairwise_list_CellTypeFinal_22CellTypes.rda"))


print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# [1] "Reproducibility information:"
# [1] "2023-12-11 14:37:41 EST"  
# user   system  elapsed 
# 428.098    2.584 1281.674 
# ─ Session info ─────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 Patched (2023-07-19 r84711)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2023-12-11
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/pandoc
# 
# ─ Packages ─────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [2] CRAN (R 4.3.1)
# beachmat               2.16.0    2023-04-25 [2] Bioconductor
# beeswarm               0.4.0     2021-06-01 [2] CRAN (R 4.3.1)
# Biobase              * 2.60.0    2023-04-25 [2] Bioconductor
# BiocGenerics         * 0.46.0    2023-04-25 [2] Bioconductor
# BiocNeighbors          1.18.0    2023-04-25 [2] Bioconductor
# BiocParallel           1.34.2    2023-05-22 [2] Bioconductor
# BiocSingular           1.16.0    2023-04-25 [2] Bioconductor
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.1)
# bluster                1.10.0    2023-04-25 [2] Bioconductor
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.1)
# cluster                2.1.4     2022-08-22 [3] CRAN (R 4.3.1)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.3.1)
# colorout             * 1.3-0.1   2023-12-01 [1] Github (jalvesaq/colorout@deda341)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.3.1)
# cowplot                1.1.1     2020-12-30 [2] CRAN (R 4.3.1)
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.1)
# DeconvoBuddies       * 0.99.0    2023-12-04 [1] Github (LieberInstitute/DeconvoBuddies@9ce4a42)
# DelayedArray           0.26.7    2023-07-28 [2] Bioconductor
# DelayedMatrixStats     1.22.6    2023-08-28 [2] Bioconductor
# dplyr                  1.1.3     2023-09-03 [2] CRAN (R 4.3.1)
# dqrng                  0.3.1     2023-08-30 [2] CRAN (R 4.3.1)
# edgeR                  3.42.4    2023-05-31 [2] Bioconductor
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.3.1)
# farver                 2.1.1     2022-07-06 [2] CRAN (R 4.3.1)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.3.1)
# GenomeInfoDb         * 1.36.3    2023-09-07 [2] Bioconductor
# GenomeInfoDbData       1.2.10    2023-07-20 [2] Bioconductor
# GenomicRanges        * 1.52.0    2023-04-25 [2] Bioconductor
# ggbeeswarm             0.7.2     2023-04-29 [2] CRAN (R 4.3.1)
# ggplot2              * 3.4.3     2023-08-14 [2] CRAN (R 4.3.1)
# ggrepel                0.9.3     2023-02-03 [2] CRAN (R 4.3.1)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.3.1)
# gridExtra              2.3       2017-09-09 [2] CRAN (R 4.3.1)
# gtable                 0.3.4     2023-08-21 [2] CRAN (R 4.3.1)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.3.1)
# igraph                 1.5.1     2023-08-10 [2] CRAN (R 4.3.1)
# IRanges              * 2.34.1    2023-06-22 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.3.1)
# labeling               0.4.3     2023-08-29 [2] CRAN (R 4.3.1)
# lattice                0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.3.1)
# limma                  3.56.2    2023-06-04 [2] Bioconductor
# locfit                 1.5-9.8   2023-06-11 [2] CRAN (R 4.3.1)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.3.1)
# Matrix                 1.6-1.1   2023-09-18 [3] CRAN (R 4.3.1)
# MatrixGenerics       * 1.12.3    2023-07-30 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.3.1)
# metapod                1.8.0     2023-04-25 [2] Bioconductor
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.3.1)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.3.1)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.3.1)
# Polychrome           * 1.5.1     2022-05-03 [1] CRAN (R 4.3.1)
# purrr                  1.0.2     2023-08-10 [2] CRAN (R 4.3.1)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.1)
# rafalib                1.0.0     2015-08-09 [1] CRAN (R 4.3.1)
# ragg                   1.2.5     2023-01-12 [2] CRAN (R 4.3.1)
# RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.3.1)
# Rcpp                   1.0.11    2023-07-06 [2] CRAN (R 4.3.1)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.1)
# rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.1)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.3.1)
# S4Arrays               1.0.6     2023-08-30 [2] Bioconductor
# S4Vectors            * 0.38.1    2023-05-02 [2] Bioconductor
# ScaledMatrix           1.8.1     2023-05-03 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.3.1)
# scater               * 1.28.0    2023-04-25 [2] Bioconductor
# scatterplot3d          0.3-44    2023-05-05 [1] CRAN (R 4.3.1)
# scran                * 1.28.2    2023-07-23 [2] Bioconductor
# scuttle              * 1.10.2    2023-08-03 [2] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.1)
# SingleCellExperiment * 1.22.0    2023-04-25 [2] Bioconductor
# sparseMatrixStats      1.12.2    2023-07-02 [2] Bioconductor
# statmod                1.5.0     2023-01-06 [2] CRAN (R 4.3.1)
# stringi                1.7.12    2023-01-11 [2] CRAN (R 4.3.1)
# stringr                1.5.0     2022-12-02 [2] CRAN (R 4.3.1)
# SummarizedExperiment * 1.30.2    2023-06-06 [2] Bioconductor
# systemfonts            1.0.4     2022-02-11 [2] CRAN (R 4.3.1)
# textshaping            0.3.6     2021-10-13 [2] CRAN (R 4.3.1)
# tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.3.1)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.3.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.3.1)
# vctrs                  0.6.3     2023-06-14 [2] CRAN (R 4.3.1)
# vipor                  0.4.5     2017-03-22 [2] CRAN (R 4.3.1)
# viridis                0.6.4     2023-07-22 [2] CRAN (R 4.3.1)
# viridisLite            0.4.2     2023-05-02 [2] CRAN (R 4.3.1)
# withr                  2.5.0     2022-03-03 [2] CRAN (R 4.3.1)
# XVector                0.40.0    2023-04-25 [2] Bioconductor
# zlibbioc               1.46.0    2023-04-25 [2] Bioconductor
# 
# [1] /users/rphillip/R/4.3
# [2] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library
# 
# ────────────────────────────────────────────────────────────────────────────────
