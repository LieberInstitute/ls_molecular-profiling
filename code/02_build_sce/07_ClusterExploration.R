#Goal: Explore clusters and identified marker genes + plot expression of known LS territories.
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

library(SingleCellExperiment)
library(ComplexHeatmap)
library(ggplot2)
library(scater)
library(dplyr)
library(scran)
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

#bargraph of number of cells in each population. 
celltype_bar <- as.data.frame(table(sce$CellType.Final)) %>% 
    rename(CellType = Var1, Number_Nuclei = Freq) %>%
    ggplot(aes(x = reorder(CellType,-Number_Nuclei), y = Number_Nuclei,fill = CellType)) +
    scale_fill_manual(values = new_cluster_cols) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(x = "Cell Type",
         y = "Number of Nuclei")

ggsave(celltype_bar,filename = here("plots","No_Nuclei_per_CellType_bar.png"))

#load the DEG lists. 
load(here("processed-data","markers_1vAll_ttest_CellTypeFinal_20Clusters.rda"),verbose = TRUE)
# Loading objects:
#     markers_1vALL_enrich_Final
load(here("processed-data","markers_pairwise_list_CellTypeFinal_20CellTypes.rda"),verbose = TRUE) 
# Loading objects:
#     markers_pairwise

#Gene list from figure 2 in Besnard + Leroy, Mol. Psy. 2022
#Fig. 2a "Gene expression domains defining non-overlapping LS territories.
domain_genes <- c("CBLN2","CBLN4","PAX6","SEMA3A", #Rostrocaudal Deep
                  "COL15A1","MATN2","WFS1", #Rostrocaudal Superficial
                  "DRD3","GALR1","IGFBP5","POU6F2", #Rostral Deep
                  "FOXP2","NDST4","GDPD2","CDH7","CRHR2","DACH2","FST","LHX2", #Rostral Superficial
                  "ASB4","OTX2","PDE1C","TAC1R","TRPC6", #Caudal Deep
                  "CD24A","IGFBP4") #Caudal Superficial 

#Some gene names might be differnt in the human genome.
#Check if all genes are in the object
table(domain_genes %in% rowData(sce)$gene_name)
# FALSE  TRUE 
#.    2    24

domain_genes[!(domain_genes %in% rowData(sce)$gene_name)]
#[1] "TAC1R" "CD24A"
#TAC1R name is same in human genome. 
#CD24a is CD24 in human genome. 

#Dot plot for domain markers
domain_dot <- plotDots(object = sce,features = rev(c("CBLN2","CBLN4","PAX6","SEMA3A", #Rostrocaudal Deep
                                                     "COL15A1","MATN2","WFS1", #Rostrocaudal Superficial
                                                     "DRD3","GALR1","IGFBP5","POU6F2", #Rostral Deep
                                                     "FOXP2","NDST4","GDPD2","CDH7","CRHR2","DACH2","FST","LHX2", #Rostral Superficial
                                                     "ASB4","OTX2","PDE1C","TRPC6", #Caudal Deep
                                                     "CD24","IGFBP4")), #Caudal Superficial 
                       swap_rownames = "gene_name",group = "CellType.Final") +
    scale_color_gradientn(colours = c("lightgrey","orange","red")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = here("plots","Expression_plots","LS_Domains_dotplot.png"),plot = domain_dot)


#Make volcano plots for these genes
domain_violin <- plotExpression(object = sce,features = c("CBLN2","CBLN4","PAX6","SEMA3A", #Rostrocaudal Deep
                                                          "COL15A1","MATN2","WFS1", #Rostrocaudal Superficial
                                                          "DRD3","GALR1","IGFBP5","POU6F2", #Rostral Deep
                                                          "FOXP2","NDST4","GDPD2","CDH7","CRHR2","DACH2","FST","LHX2", #Rostral Superficial
                                                          "ASB4","OTX2","PDE1C","TRPC6", #Caudal Deep
                                                          "CD24","IGFBP4"),
                                swap_rownames = "gene_name",x = "CellType.Final",
                                colour_by = "CellType.Final",ncol = 5) +
    theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")
ggsave(filename = here("plots","Expression_plots","LS_Domains_violin.png"),
       plot = domain_violin,width = 14,height = 16)

#tSNEs
domain_markers <- c("CBLN2","CBLN4","PAX6","SEMA3A", #Rostrocaudal Deep
                    "COL15A1","MATN2","WFS1", #Rostrocaudal Superficial
                    "DRD3","GALR1","IGFBP5","POU6F2", #Rostral Deep
                    "FOXP2","NDST4","GDPD2","CDH7","CRHR2","DACH2","FST","LHX2", #Rostral Superficial
                    "ASB4","OTX2","PDE1C","TRPC6", #Caudal Deep
                    "CD24","IGFBP4")

for(i in domain_markers){
    tSNE_domain <- plotReducedDim(object = sce,
                                  dimred = "tSNE_mnn_15",
                                  colour_by = i,
                                  swap_rownames = "gene_name") +
        scale_color_gradientn(colours = c("lightgrey","red")) +
        ggtitle(i) +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(filename = paste0("plots/Expression_plots/domain_markers_tSNE_mnn_15/",i,"_expression_tSNE_15_dims.pdf"),
           plot = tSNE_domain,
           height = 8,width = 8)
}

# #Additional investigation of marker genes. 
# LS_GABA_1_markers <- subset(markers_1vALL_df,subset=(cellType.target == "LS_GABA_1"))[1:50,]
# LS_GABA_2_markers <- subset(markers_1vALL_df,subset=(cellType.target == "LS_GABA_2"))[1:50,]
# for(i in 1:nrow(LS_GABA_2_markers)){
#     print(i)
#     x <- plotReducedDim(object = sce,
#                         dimred = "UMAP_mnn_15",
#                         colour_by = LS_GABA_2_markers[i,"gene_name"],
#                         swap_rownames = "gene_name") +
#         scale_color_gradientn(colours = c("lightgrey","orange","red"))
#     ggsave(plot = x,
#            filename = here("plots","Expression_plots","LS_GABA_2_Markers",
#                            paste0(i,"_",LS_GABA_2_markers[i,"gene_name"],"_umap.pdf")))
# }
# p1 <- plotReducedDim(object = sce,
#                dimred = "UMAP_mnn_15",
#                colour_by = "DGKG",
#                swap_rownames = "gene_name") +
#     scale_color_gradientn(colours = c("lightgrey","orange","red")) +
#     ggtitle("DGKG (LS Marker)") +
#     theme(plot.title = element_text(hjust = 0.5))
# 
# p2 <- plotReducedDim(object = sce,
#                      dimred = "UMAP_mnn_15",
#                      colour_by = "TRPC4",
#                      swap_rownames = "gene_name") +
#     scale_color_gradientn(colours = c("lightgrey","orange","red")) +
#     ggtitle("TRPC4 (LS Marker)")
# 
# 
# p3 <- plotReducedDim(object = sce,
#                      dimred = "UMAP_mnn_15",
#                      colour_by = "RARB",
#                      swap_rownames = "gene_name") +
#     scale_color_gradientn(colours = c("lightgrey","orange","red")) +
#     ggtitle("RARB (Striatal Marker)") +
#     theme(plot.title = element_text(hjust = 0.5))
# 
# p4 <- plotReducedDim(object = sce,
#                      dimred = "UMAP_mnn_15",
#                      colour_by = "BCL11B",
#                      swap_rownames = "gene_name") +
#     scale_color_gradientn(colours = c("lightgrey","orange","red")) +
#     ggtitle("BCL11B (Striatal Marker)") +
#     theme(plot.title = element_text(hjust = 0.5))
# 
# x <- cowplot::plot_grid(plotlist = list(p1,p2,p3,p4),ncol = 2)
# ggsave(x, filename = here("plots","Expression_plots","LS_Str_Markers.pdf"))



#####Cluster modularity
#Build the graph again. 
#Used k=50 + walktrap clustering for celltype designations. 
snn_k_20 <- buildSNNGraph(sce, k = 20, use.dimred = "mnn",type="jaccard")

# "Rather, we use the pairwiseModularity() function from bluster with as.ratio=TRUE, 
# which returns the ratio of the observed to expected sum of weights between each pair of clusters. 
# We use the ratio instead of the difference as the former is less dependent on the number of cells 
# in each cluster" - OSCA advanced, Chatper 5.2.5
k_20_modularity <- bluster::pairwiseModularity(graph = snn_k_20,
                                               clusters = sce$CellType.Final,
                                               as.ratio = TRUE)


library(pheatmap)
pdf(file = here("plots","k_20_pairwise_modularity_final_celltypes.pdf"))
pheatmap(log2(k_20_modularity+1), 
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         display_numbers=TRUE, 
         number_format="%.2f", 
         fontsize_number=6.5,
         main = "Modularity ratio for 20 graph-based clusters in human LS (n=3)",
         color=colorRampPalette(c("white","orange","red"))(100))
dev.off()
# "Indeed, concentration of the weight on the diagonal of (Figure 5.5) indicates 
# that most of the clusters are well-separated, while some modest off-diagonal entries 
# represent closely related clusters with more inter-connecting edges." - OSCA advanced, chatper 5.2.5

cluster.gr <- igraph::graph_from_adjacency_matrix(log2(k_20_modularity+1),
                                                  mode="upper", weighted=TRUE, diag=FALSE)

# Increasing the weight to increase the visibility of the lines.
set.seed(11001010)
pdf(file = here("plots","k_20_relationship_between_clusters_graph.pdf"))
plot(cluster.gr, edge.width=igraph::E(cluster.gr)$weight*10,
     layout=igraph::layout_with_lgl)
dev.off()
# Force-based layout showing the relationships between clusters based on the log-ratio of observed to expected total weights 
# between nodes in different clusters. The thickness of the edge between a pair of clusters is proportional to the corresponding 
# log-ratio. - OSCA advanced, chapter 5.2.5, figure 5.6
                
sce$CellType.Final <- factor(x = sce$CellType.Final,
                             levels = c("LS_Inh_A","LS_Inh_B","LS_Inh_G","LS_Inh_I",
                                        "MS_Inh_A","MS_Inh_E","MS_Inh_H","Sept_Inh_D",
                                        "Sept_Inh_F","Str_Inh_A","Str_Inh_B","MS_Excit_A",
                                        "Excit_A","Excit_B","Oligo","Polydendrocyte",
                                        "Astrocyte","Ependymal","Microglia","Mural"))

x <- plotExpression(object = sce,features = c("SYT1","GAD1","SLC17A6","MOBP"),
                    x = "CellType.Final",
                    swap_rownames = "gene_name",
                    ncol = 2,colour_by = "CellType.Final") +
    scale_color_manual(values = new_cluster_cols) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    stat_summary(fun = median, 
                 fun.min = median, 
                 fun.max = median,
                 geom = "crossbar", 
                 width = 0.3)

ggsave(file = "plots/Expression_plots/violin_General_markers.pdf",plot = x,height = 8,width = 8)

#############################################
##Heatmap of genes that I want to highlight##
#############################################
#Code from https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/Complex%20Heatmap.R
splitit <- function(x) split(seq(along = x), x)

cell_idx <- splitit(sce$CellType.Final)
dat <- as.matrix(logcounts(sce))
rownames(dat) <- rowData(sce)$gene_name
dim(dat)
# [1] 36601  9225

############set up columns for heatmaps. 
#Set marker genes to be included on the heatmap.
markers_all <- c("RBFOX3","SNAP25","SYT1",#Pan-neuronal
                 "GAD1","GAD2", "SLC32A1", #Inhibitory
                 "SLC17A6","SLC17A7", #Excitatory
                 "MOBP","MBP", #Oligodendrocyte
                 "PDGFRA","CSPG4", #Polydendrocyte
                 "GFAP","SLC1A2", #Astrocyte
                 "CFAP44","FOXJ1",#Ependymal
                 "TMEM119","C3", #Microglia
                 "RGS5","CLDN5") #Mural

#marker labels
marker_labels <- c(rep("neuronal",3),
                   rep("Inhibitory",3),
                   rep("Excitatory",2),
                   rep("Oligodendrocyte",2),
                   rep("Polydendrocyte",2),
                   rep("Astrocyte",2),
                   rep("Ependymal",2),
                   rep("Microglia",2),
                   rep("Mural",2))

marker_labels <- factor(x = marker_labels,
                        levels = c("neuronal","Inhibitory","Excitatory",
                                   "Oligodendrocyte","Polydendrocyte","Astrocyte",
                                   "Ependymal","Microglia","Mural"))

colors_markers <- list(marker = c(neuronal = "black",
                                  Inhibitory = "#D62728",
                                  Excitatory = "#98DF8A",
                                  Oligodendrocyte = "#FF9E4A",
                                  Polydendrocyte = "#9467BD",
                                  Astrocyte = "#8C564B",
                                  Ependymal = "#9EDAE5",
                                  Microglia = "#17BECF",
                                  Mural = "#1F77B4"))

col_ha <- ComplexHeatmap::columnAnnotation(marker = marker_labels,
                                           show_annotation_name = FALSE,
                                           show_legend = FALSE,
                                           col = colors_markers)

###########set up rows for heatmap. 
# cluster labels
cluster_pops <- list(Inhibitory = c("LS_Inh_A","LS_Inh_B","LS_Inh_G","LS_Inh_I",
                                    "MS_Inh_A","MS_Inh_E","MS_Inh_H","Sept_Inh_D",
                                    "Sept_Inh_F","Str_Inh_A","Str_Inh_B"),
                     Excitatory = c("MS_Excit_A","Excit_A","Excit_B"),
                     Oligodendrocyte = "Oligodendrocyte",
                     Polydendrocyte = "Polydendrocyte",
                     Astrocyte = "Astrocyte",
                     Ependymal = "Ependymal",
                     Microglia = "Microglia",
                     Mural = "Mural")

# cluster labels order
cluster_pops_order <- unname(unlist(cluster_pops))

# swap values and names of list
cluster_pops_rev <- rep(names(cluster_pops), 
                        times = sapply(cluster_pops, length))
names(cluster_pops_rev) <- unname(unlist(cluster_pops))
cluster_pops_rev <- cluster_pops_rev[as.character(sort(cluster_pops_order))]
cluster_pops_rev <- factor(cluster_pops_rev, levels = names(cluster_pops))

# second set of cluster labels
neuron_pops <- ifelse(cluster_pops_rev %in% c("Excitatory", "Inhibitory"),
                      "Neuronal",
                      "Non-neuronal")
neuron_pops <- factor(x = neuron_pops,levels = c("Neuronal","Non-neuronal"))


colors_neurons <- list(class = c(Neuronal = "black",
                                 `Non-neuronal` = "gray90"))

n <- table(sce$CellType.Final)

#row annotation dataframe. 
# row annotation
library(ComplexHeatmap)
pop_markers <- list(population = c(Inhibitory = "#D62728",
                                   Excitatory = "#98DF8A",
                                   Oligodendrocyte = "#FF9E4A",
                                   Polydendrocyte = "#9467BD",
                                   Astrocyte = "#8C564B",
                                   Ependymal = "#9EDAE5",
                                   Microglia = "#17BECF",
                                   Mural = "#1F77B4"))


row_ha <- rowAnnotation(n = anno_barplot(as.numeric(n), 
                                         gp = gpar(fill = "navy"), 
                                         border = FALSE),
                        class = neuron_pops,
                        population = cluster_pops_rev,
                        show_annotation_name = FALSE,
                        col = c(pop_markers,colors_neurons))


hm_mat <- scale(t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers_all, i])))))


hm <- ComplexHeatmap::Heatmap(matrix = hm_mat,
                              name = "centered,scaled",
                              column_title = "General cell class marker \ngene expression across clusters",
                              column_title_gp = gpar(fontface = "bold"),
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              bottom_annotation = col_ha,
                              right_annotation = row_ha,
                              column_split = marker_labels,
                              row_split = cluster_pops_rev,
                              row_title = NULL,
                              rect_gp = gpar(col = "gray50", lwd = 0.5))


pdf(here("plots","ComplexHeatmap_General_cell_class_markers.pdf"))
hm
dev.off()

###################################################
##HEATMAP HIGHLIGHTING  NEURONAL POPULATIONS ONLY##
###################################################
#Code from https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/Complex%20Heatmap.R
splitit <- function(x) split(seq(along = x), x)



sce_neuronal <- sce[,sce$CellType.Final %in% c("LS_Inh_A","LS_Inh_B","LS_Inh_G","LS_Inh_I",
                                               "MS_Inh_A","MS_Inh_E","MS_Inh_H","MS_Excit_A",
                                               "Sept_Inh_D","Sept_Inh_F","Str_Inh_A","Str_Inh_B",
                                               "Excit_A","Excit_B")]


sce_neuronal
# class: SingleCellExperiment 
# dim: 36601 5387 
# metadata(1): Samples
# assays(3): counts binomial_deviance_residuals logcounts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(5387): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 3_TTTGGTTGTGTTAGCT-1 3_TTTGGTTTCTTCGACC-1
# colData names(61): Sample Barcode ... CellType_k_20_louvain
# CellType.Final
# reducedDimNames(18): GLMPCA_approx UMAP_15 ... tSNE_mnn_25 tSNE_mnn_50
# mainExpName: NULL
# altExpNames(0):


cell_idx <- splitit(sce_neuronal$CellType.Final)
dat <- as.matrix(logcounts(sce_neuronal))
rownames(dat) <- rowData(sce_neuronal)$gene_name
dim(dat)
# [1] 36601  5387

############set up columns for heatmaps. 
#Set marker genes to be included on the heatmap.
markers_all <- c("TRPC4","DGKG","CRHR2",#Broad LS
                 "GPR26","COL25A1", #LS A
                 "MYO5B","FREM2", #LS B
                 "HDC","OPRM1", #LS G
                 "SCML4","TMEM215", #LS I 
                 "ELAVL2","ELAVL4",#Broad MS
                 "TACR1", #MS SUB
                 "FXYD6","GRIN2D", "KCNC2",#Septal
                 "CHAT","KIT","SST","CCK","VIP",#Interneuronal. 
                 "RARB","BCL11B", #Striatal 
                 "ISL1","FOXP2", "DRD1","EBF1","DRD2","PENK",#MSN markers
                 "SLC17A7","SLC17A6") #Excitatory

#marker labels
marker_labels <- c(rep("LS-Broad",3),
                   rep("LS-subclusters",8),
                   rep("MS-Broad",3),
                   rep("Septal",3),
                   rep("Interneuron",5),
                   rep("Striatal",2),
                   rep("MSN",6),
                   rep("Excitatory",2))

marker_labels <- factor(x = marker_labels,
                        levels = c("LS-Broad","LS-subclusters","MS-Broad","Septal",
                                   "Interneuron","Striatal","MSN","Excitatory"))

colors_markers <- list(marker = c(`LS-Broad` = "#FF9E4A",
                                  `LS-subclusters`= "#D62728",
                                  `MS-Broad` = "#9EDAE5",
                                  Septal = "#17BECF",
                                  Interneuron = "#9467BD",
                                  Excitatory = "#98DF8A",
                                  Striatal = "#8C564B",
                                  MSN = "#1F77B4"))

col_ha <- ComplexHeatmap::columnAnnotation(marker = marker_labels,
                                           show_annotation_name = FALSE,
                                           show_legend = FALSE,
                                           col = colors_markers)

col_ha_legend <- ComplexHeatmap::columnAnnotation(marker = marker_labels,
                                                  show_annotation_name = TRUE,
                                                  show_legend = TRUE,
                                                  col = colors_markers)


# ###########set up rows for heatmap. 
# # cluster labels
cluster_pops <- list(LS = c("LS_Inh_A","LS_Inh_B","LS_Inh_G","LS_Inh_I"),
                     MS = c("MS_Inh_A","MS_Inh_E","MS_Inh_H"),
                     Sept = c("Sept_Inh_D","Sept_Inh_F"),
                     Str = c("Str_Inh_A","Str_Inh_B"),
                     Excit = c("MS_Excit_A","Excit_A","Excit_B"))
# # cluster labels order
cluster_pops_order <- unname(unlist(cluster_pops))

# swap values and names of list
cluster_pops_rev <- rep(names(cluster_pops),
                        times = sapply(cluster_pops, length))
names(cluster_pops_rev) <- unname(unlist(cluster_pops))
cluster_pops_rev <- cluster_pops_rev[as.character(sort(cluster_pops_order))]
cluster_pops_rev <- factor(cluster_pops_rev, levels = names(cluster_pops))

n <- table(sce_neuronal$CellType.Final)

#row annotation dataframe.
# row annotation
pop_markers <- list(population = c(LS = "#FF9E4A",
                                   MS = "#9EDAE5",
                                   Sept = "#9467BD",
                                   Str = "#98DF8A",
                                   Excit = "#1F77B4"))


row_ha <- rowAnnotation(n = anno_barplot(as.numeric(n),
                                         gp = gpar(fill = "navy"),
                                         border = FALSE),
                        population = cluster_pops_rev,
                        show_annotation_name = FALSE,
                        col = pop_markers)


hm_mat <- scale(t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers_all, i])))))


hm <- ComplexHeatmap::Heatmap(matrix = hm_mat,
                              name = "centered,scaled",
                              column_title = "Marker gene expression\n across neuronal clusters",
                              column_title_gp = gpar(fontface = "bold"),
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              bottom_annotation = col_ha,
                              right_annotation = row_ha,
                              column_split = marker_labels,
                              row_split = cluster_pops_rev,
                              row_title = NULL,
                              rect_gp = gpar(col = "gray50", lwd = 0.5))


pdf(here("plots","ComplexHeatmap_Neuronal_cell_class_markers.pdf"),height = 10,width = 10)
hm
dev.off()

#####Heatmap with legend. 
hm <- ComplexHeatmap::Heatmap(matrix = hm_mat,
                              name = "centered,scaled",
                              column_title = "Marker gene expression\n across neuronal clusters",
                              column_title_gp = gpar(fontface = "bold"),
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              bottom_annotation = col_ha_legend,
                              right_annotation = row_ha,
                              column_split = marker_labels,
                              row_split = cluster_pops_rev,
                              row_title = NULL,
                              rect_gp = gpar(col = "gray50", lwd = 0.5))


pdf(here("plots","ComplexHeatmap_Neuronal_cell_class_markers_with_Legend.pdf"))
hm
dev.off()



#####Violin plots for figure 2. 
load(here("processed-data","Final_CellTypes_colors.rda"),verbose = TRUE)


sce_neuronal$CellType.Final <- factor(sce_neuronal$CellType.Final,
                                      levels = c("LS_Inh_A","LS_Inh_B","LS_Inh_G","LS_Inh_I",
                                                 "MS_Inh_A","MS_Inh_E","MS_Inh_H","MS_Excit_A",
                                                 "Sept_Inh_D","Sept_Inh_F","Excit_A","Excit_B",
                                                 "Str_Inh_A","Str_Inh_B"))

vln_fig2 <- plotExpression(object = sce_neuronal,
                           features = c("TRPC4","DGKG",
                                        "CRHR2","FREM2",
                                        "OPRM1","GPR26",
                                        "ELAVL2","ELAVL4",
                                        "FXYD6","RARB"),
                           color_by = "CellType.Final",
                           x = "CellType.Final",
                           swap_rownames = "gene_name",
                           ncol = 2) +
    stat_summary(fun = median, 
                 fun.min = median, 
                 fun.max = median,
                 geom = "crossbar", 
                 width = 0.3) +
    scale_color_manual(values  = cluster_cols) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90))

ggsave(plot = vln_fig2,
       filename = here("plots","Expression_plots","Violin_neuronal_markers.pdf"),height = 12,width = 8)


#Make individual violins for 4 genes
trpc4_vln <- plotExpression(object = sce_neuronal,
                            features = "TRPC4",
                            colour_by = "CellType.Final",
                            x = "CellType.Final",
                            swap_rownames = "gene_name") +
    stat_summary(fun = median, 
                 fun.min = median, 
                 fun.max = median,
                 geom = "crossbar", 
                 width = 0.3) +
    scale_color_manual(values  = cluster_cols) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90))
ggsave(plot = trpc4_vln,
       filename = here("plots","Expression_plots","trpc4_violin.pdf"),height = 5, width = 5)

fxyd6_vln <- plotExpression(object = sce_neuronal,
                            features = "FXYD6",
                            colour_by = "CellType.Final",
                            x = "CellType.Final",
                            swap_rownames = "gene_name") +
    stat_summary(fun = median, 
                 fun.min = median, 
                 fun.max = median,
                 geom = "crossbar", 
                 width = 0.3) +
    scale_color_manual(values  = cluster_cols) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90))
ggsave(plot = fxyd6_vln,
       filename = here("plots","Expression_plots","fxyd6_violin.pdf"),height = 5, width = 5)


OPRM1_vln <- plotExpression(object = sce_neuronal,
                            features = "OPRM1",
                            colour_by = "CellType.Final",
                            x = "CellType.Final",
                            swap_rownames = "gene_name") +
    stat_summary(fun = median, 
                 fun.min = median, 
                 fun.max = median,
                 geom = "crossbar", 
                 width = 0.3) +
    scale_color_manual(values  = cluster_cols) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90))
ggsave(plot = OPRM1_vln,
       filename = here("plots","Expression_plots","oprm1_violin.pdf"),height = 5, width = 5)

FREM2_vln <- plotExpression(object = sce_neuronal,
                            features = "FREM2",
                            colour_by = "CellType.Final",
                            x = "CellType.Final",
                            swap_rownames = "gene_name") +
    stat_summary(fun = median, 
                 fun.min = median, 
                 fun.max = median,
                 geom = "crossbar", 
                 width = 0.3) +
    scale_color_manual(values  = cluster_cols) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90))
ggsave(plot = FREM2_vln,
       filename = here("plots","Expression_plots","FREM2_violin.pdf"),height = 5, width = 5)


MYO5B_vln <- plotExpression(object = sce_neuronal,
                            features = "MYO5B",
                            colour_by = "CellType.Final",
                            x = "CellType.Final",
                            swap_rownames = "gene_name") +
    stat_summary(fun = median, 
                 fun.min = median, 
                 fun.max = median,
                 geom = "crossbar", 
                 width = 0.3) +
    scale_color_manual(values  = cluster_cols) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90))
ggsave(plot = MYO5B_vln,
       filename = here("plots","Expression_plots","MYO5B_violin.pdf"),height = 5, width = 5)


#FREM2 feature plot
frem2_featureplot <- plotReducedDim(object = sce,dimred = "tSNE_mnn_15",
                                    color_by = "FREM2",swap_rownames = "gene_name") +
    scale_color_gradientn(colours = c("lightgrey","orange","red"))
ggsave(plot = frem2_featureplot,
       filename = here("plots","Expression_plots","FREM2_featureplot.pdf"))

#MYO5B feature plot
MYO5B_featureplot <- plotReducedDim(object = sce,dimred = "tSNE_mnn_15",
                                    color_by = "MYO5B",swap_rownames = "gene_name") +
    scale_color_gradientn(colours = c("lightgrey","orange","red"))
ggsave(plot = MYO5B_featureplot,
       filename = here("plots","Expression_plots","myo5b_featureplot.pdf"))




