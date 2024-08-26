#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

library(SingleCellExperiment)
library(ComplexHeatmap)
library(sessioninfo)
library(scater)
library(here)

load(here("processed-data","02_build_sce","sce_celltype_082624.rda"),verbose = TRUE)

sce

#Load the cluster colors
load(here("processed-data","cluster_cols_CellType_Final_082524.rda"),verbose = TRUE)

###Feature plots for SYT1,SNAP25, GAD1, GAD2, SLC17A6, SLC17A7, MOBP, GFAP
#Make each plot with and without the figure legend.  
#With axes --> Save as pdf. 
#Without axes --> Save as png
#Theme_void removes the axes. 
######Without axes 
for(i in c("SYT1","SNAP25","GAD1","GAD2","SLC17A6","SLC17A7","MOBP","GFAP")){
  x <- plotReducedDim(sce,
                      dimred = "tSNE_mnn_50",
                      colour_by = i,
                      swap_rownames = "gene_name") +
    scale_color_gradientn(colours = c("lightgrey","orange","red")) +
    theme_void() +
    theme(legend.position = "none")
  ggsave(x,filename = here("plots","Figure1_Plots","FeaturePlot_NoLegend",
                           paste0(i,"_FeaturePlot_NoLegend.png")),
         height = 6, width = 6)
}

######With axes
for(i in c("SYT1","SNAP25","GAD1","GAD2","SLC17A6","SLC17A7","MOBP","GFAP")){
  x <- plotReducedDim(sce,
                      dimred = "tSNE_mnn_50",
                      colour_by = i,
                      swap_rownames = "gene_name") +
    scale_color_gradientn(colours = c("lightgrey","orange","red")) +
    theme_void() 
  ggsave(x,filename = here("plots","Figure1_Plots","FeaturePlot_Legend",
                           paste0(i,"_FeaturePlot_Legend.png")))
}


#############################################
##Complex Heatmap##
#############################################
#Code from https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/Complex%20Heatmap.R
splitit <- function(x) split(seq(along = x), x)

cell_idx <- splitit(sce$CellType.Final)

#Pull logcounts matrix
dat <- as.matrix(logcounts(sce))

#Make sure rownames of dat are in same order as rowData(sce)
identical(rownames(dat),rowData(sce)$gene_id)
stopifnot(identical(rownames(dat),rowData(sce)$gene_id))

#Change the rownames to the gene symbol (here called gene_name)
rownames(dat) <- rowData(sce)$gene_name
dim(dat)

############set up columns for heatmaps. 
#Set marker genes to be included on the heatmap.
markers_all <- c("RBFOX3","SNAP25","SYT1",#Pan-neuronal 3
                 "GAD1","GAD2", "SLC32A1", #Inhibitory 3
                 "SLC17A6","SLC17A7", #Excitatory 2
                 "MOBP","MBP", #Oligodendrocyte 2
                 "PDGFRA","CSPG4", #Polydendrocyte 2
                 "GFAP","SLC1A2", #Astrocyte 2
                 "CFAP44","FOXJ1",#Ependymal 2
                 "TMEM119","C3", #Microglia 2
                 "RGS5","CLDN5") #Mural 2

#marker labels
marker_labels <- c(rep("neuronal",3),
                   rep("Inhibitory",3),
                   rep("Excitatory",2),
                   rep("Oligodendrocyte",2),
                   rep("OPC",2),
                   rep("Astrocyte",2),
                   rep("Ependymal",2),
                   rep("Microglia",2),
                   rep("Mural",2))

marker_labels <- factor(x = marker_labels,
                        levels = c("neuronal","Inhibitory","Excitatory",
                                   "Oligodendrocyte","OPC","Astrocyte",
                                   "Ependymal","Microglia","Mural"))

colors_markers <- list(marker = c(neuronal = "black",
                                  Inhibitory = "#D62728",
                                  Excitatory = "#0a99c0",
                                  Oligodendrocyte = "#32FF0D",
                                  OPC = "#0D996A",
                                  Astrocyte = "#B9C4FB",
                                  Ependymal = "#800D91",
                                  Microglia = "#FC00FC",
                                  Mural = "#32FEA1"))

col_ha <- ComplexHeatmap::columnAnnotation(marker = marker_labels,
                                           show_annotation_name = FALSE,
                                           show_legend = FALSE,
                                           col = colors_markers)

###########set up rows for heatmap. 
# cluster labels
cluster_pops <- list(Inhibitory = c("LS_Inh_A","LS_Inh_B","LS_Inh_C","LS_Inh_D","LS_Inh_E",
                                    "MS_Inh_A","MS_Inh_B","MS_Inh_C","MS_Inh_D",
                                    "Sept_Inh_A","Sept_Inh_B","Sept_Inh_C",
                                    "Str_DRD1_MSN_A","Str_DRD1_MSN_B","Str_DRD1_Patch","Str_DRD2_MSN"),
                     Excitatory = c("MS_Excit_A","Excitatory_A","Excitatory_B"),
                     Oligodendrocyte = "Oligo",
                     OPC = "OPC",
                     Astrocyte = "Astrocyte",
                     Ependymal = "Ependymal",
                     Microglia = "Microglia",
                     Mural = "Mural")

# cluster labels order
# # cluster labels order
cluster_pops_order <- unname(unlist(cluster_pops))

# swap values and names of list
cluster_pops_rev <- rep(names(cluster_pops),
                        times = sapply(cluster_pops, length))
names(cluster_pops_rev) <- unname(unlist(cluster_pops))
#cluster_pops_rev <- cluster_pops_rev[as.character(sort(cluster_pops_order))]
cluster_pops_rev <- factor(cluster_pops_rev, levels = names(cluster_pops))

# second set of cluster labels
neuron_pops <- ifelse(cluster_pops_rev %in% c("Inhibitory","Excitatory"),
                      "Neuronal",
                      "Non-neuronal")

neuron_pops <- factor(x = neuron_pops,levels = c("Neuronal","Non-neuronal"))


colors_neurons <- list(class = c(Neuronal = "black",
                                 `Non-neuronal` = "gray90"))

n <- table(sce$CellType.Final)

#row annotation dataframe. 
# row annotation
pop_markers <- list(population = c(Inhibitory = "#D62728",
                                   Excitatory = "#0a99c0",
                                   Oligodendrocyte = "#32FF0D",
                                   OPC = "#0D996A",
                                   Astrocyte = "#B9C4FB",
                                   Ependymal = "#800D91",
                                   Microglia = "#FC00FC",
                                   Mural = "#32FEA1"))


row_ha <- rowAnnotation(n = anno_barplot(as.numeric(n), 
                                         gp = gpar(fill = "navy"), 
                                         border = FALSE),
                        class = neuron_pops,
                        population = cluster_pops_rev,
                        show_annotation_name = FALSE,
                        col = c(pop_markers,colors_neurons))


hm_mat <- scale(t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers_all, i])))),
                center = TRUE, scale = TRUE)


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

pdf(here("plots","Figure1_Plots","ComplexHeatmap_General_cell_class_markers.pdf"))
hm
dev.off()

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
