#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

library(SingleCellExperiment)
library(ComplexHeatmap)
library(sessioninfo)
library(scater)
library(here)

#load sce object
load(here("processed-data","02_build_sce","sce_celltype_082724.rda"),verbose = TRUE)

sce

#load 1vALL DEGs
load(here("processed-data","markers_1vAll_ttest_k_20_louvain_CellType_Final.rda"),verbose = TRUE)

#Create a column for neurons and non-neurons
sce$Neuronal_NonNeuronal <- ifelse(sce$CellType.Final %in% c("LS_Inh_A","LS_Inh_B","LS_Inh_C","LS_Inh_D",
                                                             "LS_Inh_E","MS_Inh_A","MS_Inh_B","MS_Inh_C",
                                                             "MS_Inh_D","Sept_Inh_A","Sept_Inh_B","Sept_Inh_C",
                                                             "Str_DRD1_MSN_A","Str_DRD1_MSN_B","Str_DRD1_Patch","Str_DRD2_MSN",
                                                             "MS_Excit_A","Excit_A","Excit_B"),
                                   "Neuronal",
                                   "Non-neuronal")

#Susbet for neurons only
sce_sub <- sce[,sce$Neuronal_NonNeuronal == "Neuronal"]

###Violin plots for FXYD6, TRPC4, and OPRM1
#Match colors above. 
neuronal_colors <- c(rep("#FF9E4A",5),
                     rep("#9EDAE5",4),
                     rep("#9467BD",3),
                     rep("#8C564B",4),
                     rep("#98DF8A",3))
names(neuronal_colors) <- c("LS_Inh_A","LS_Inh_B","LS_Inh_C","LS_Inh_D","LS_Inh_E",
                            "MS_Inh_A","MS_Inh_B","MS_Inh_C","MS_Inh_D",
                            "Sept_Inh_A","Sept_Inh_B","Sept_Inh_C",
                            "Str_DRD1_MSN_A","Str_DRD1_MSN_B","Str_DRD1_Patch","Str_DRD2_MSN",
                            "MS_Excit_A","Excit_A","Excit_B")

#Now make the plots
##FXYD6
fxyd6 <- plotExpression(sce_sub,features = "FXYD6",
                        x = "CellType.Final",colour_by = "CellType.Final",
                        swap_rownames = "gene_name") +
  scale_color_manual(values = neuronal_colors) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3)
ggsave(plot = fxyd6,
       filename = here("plots","Figure2_Plots","FXYD6_violin.pdf"),height = 5, width = 5)

##TRPC4
trpc4 <- plotExpression(sce_sub,features = "TRPC4",
                        x = "CellType.Final",colour_by = "CellType.Final",
                        swap_rownames = "gene_name") +
  scale_color_manual(values = neuronal_colors) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3)
ggsave(plot = trpc4,
       filename = here("plots","Figure2_Plots","TRPC4_violin.pdf"),height = 5, width = 5)

##OPRM1
oprm1 <- plotExpression(sce_sub,features = "OPRM1",
                        x = "CellType.Final",colour_by = "CellType.Final",
                        swap_rownames = "gene_name") +
  scale_color_manual(values = neuronal_colors) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3)
ggsave(plot = oprm1,
       filename = here("plots","Figure2_Plots","OPRM1_violin.pdf"),height = 5, width = 5)


##########Feature plots
######Without axes 
for(i in c("FXYD6","TRPC4","OPRM1")){
  x <- plotReducedDim(sce,
                      dimred = "tSNE_mnn_50",
                      colour_by = i,
                      swap_rownames = "gene_name") +
    scale_color_gradientn(colours = c("lightgrey","orange","red")) +
    theme_void() +
    theme(legend.position = "none")
  ggsave(x,filename = here("plots","Figure2_Plots","FeaturePlot_NoLegend",
                           paste0(i,"_FeaturePlot_NoLegend.png")),
         height = 6, width = 6)
}

######With axes
for(i in c("FXYD6","TRPC4","OPRM1")){
  x <- plotReducedDim(sce,
                      dimred = "tSNE_mnn_50",
                      colour_by = i,
                      swap_rownames = "gene_name") +
    scale_color_gradientn(colours = c("lightgrey","orange","red")) +
    theme_void() 
  ggsave(x,filename = here("plots","Figure2_Plots","FeaturePlot_Legend",
                           paste0(i,"_FeaturePlot_Legend.png")))
}

###################################################
##HEATMAP HIGHLIGHTING  NEURONAL POPULATIONS ONLY##
###################################################
#Code from https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/Complex%20Heatmap.R
splitit <- function(x) split(seq(along = x), x)

sce_sub$CellType.Final <- as.character(sce_sub$CellType.Final)

cell_idx <- splitit(sce_sub$CellType.Final)

#Pull count matrix
dat <- as.matrix(logcounts(sce_sub))

#Make sure rowdata in same order as dat
identical(rownames(dat),rownames(rowData(sce_sub)))
stopifnot(identical(rownames(dat),rownames(rowData(sce_sub))))

rownames(dat) <- rowData(sce_sub)$gene_name
dim(dat)


############set up columns for heatmaps. 
#Set marker genes to be included on the heatmap.
markers_all <- c("TRPC4","DGKG","CRHR2",#Broad LS
                 "GPR26","PAX6", #LS A
                 "FGF19","FREM2", #LS B
                 "HS3ST2","MYO5B", #LS C
                 "HDC","SP8", #LS D
                 "SCML4","TMEM215", #LS E
                 "ADARB2","SOX6","ELAVL2", #Broad MS
                 "GBX1","RAB3B",#MS_A
                 "CHAT","SLC5A7",#MS_B
                 "CCK","VIP",#MS_C
                 "KIT","LHX6", #MS_D
                 "FXYD6","DNER","SST",#Septal
                 "RARB","BCL11B","PPP1R1B", #Striatal 
                 "ISL1","FOXP2", "DRD1","DRD2",#MSN markers
                 "SLC17A6","SLC17A7") #Excitatory

#marker labels
marker_labels <- c(rep("LS-Broad",3),
                   rep("LS-subclusters",10),
                   rep("MS-Broad",3),
                   rep("MS-subclusters",8),
                   rep("Septal",3),
                   rep("Striatal",3),
                   rep("MSN",4),
                   rep("Excitatory",2))

marker_labels <- factor(x = marker_labels,
                        levels = c("LS-Broad","LS-subclusters","MS-Broad","MS-subclusters",
                                   "Septal","Striatal","MSN","Excitatory"))

colors_markers <- list(marker = c(`LS-Broad` = "#FF9E4A",
                                  `LS-subclusters`= "#D62728",
                                  `MS-Broad` = "#9EDAE5",
                                  `MS-subclusters` = "#3b8314", 
                                  Septal = "#17BECF",
                                  Interneuron = "#9467BD",
                                  Striatal = "#8C564B",
                                  MSN = "#1F77B4",
                                  Excitatory = "#98DF8A"))

col_ha <- ComplexHeatmap::columnAnnotation(marker = marker_labels,
                                           show_annotation_name = FALSE,
                                           show_legend = TRUE,
                                           col = colors_markers)

col_ha_legend <- ComplexHeatmap::columnAnnotation(marker = marker_labels,
                                                  show_annotation_name = FALSE,
                                                  show_legend = TRUE,
                                                  col = colors_markers)


# ###########set up rows for heatmap. 
# # cluster labels
cluster_pops <- list(LS = c("LS_Inh_A","LS_Inh_B","LS_Inh_C",
                            "LS_Inh_D","LS_Inh_E"),
                     MS = c("MS_Inh_A","MS_Inh_B","MS_Inh_C",
                            "MS_Inh_D"),
                     Sept = c("Sept_Inh_A","Sept_Inh_B","Sept_Inh_C"),
                     Str = c("Str_DRD1_MSN_A","Str_DRD1_MSN_B",
                             "Str_DRD1_Patch","Str_DRD2_MSN"),
                     Excit = c("MS_Excit_A","Excit_A","Excit_B"))
# # cluster labels order
cluster_pops_order <- unname(unlist(cluster_pops))

# swap values and names of list
cluster_pops_rev <- rep(names(cluster_pops),
                        times = sapply(cluster_pops, length))

names(cluster_pops_rev) <- cluster_pops_order

cluster_pops_rev <- factor(cluster_pops_rev, levels = names(cluster_pops))

n <- table(sce_sub$CellType.Final)

#row annotation dataframe.
# row annotation
pop_markers <- list(population = c(LS = "#FF9E4A",
                                   MS = "#9EDAE5",
                                   Sept = "#9467BD",
                                   Str = "#8C564B",
                                   Excit = "#98DF8A"))


row_ha <- rowAnnotation(population = cluster_pops_rev,
                        show_annotation_name = FALSE,
                        show_legend = TRUE,
                        col = pop_markers)


hm_mat <- scale(t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers_all, i])))),
                center = TRUE,scale = TRUE)

hm_mat <- hm_mat[names(cluster_pops_rev),]


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


pdf(here("plots","Figure2_Plots","ComplexHeatmap_Neuronal_cell_class_markers.pdf"),height = 10,width = 10)
hm
dev.off()

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()




