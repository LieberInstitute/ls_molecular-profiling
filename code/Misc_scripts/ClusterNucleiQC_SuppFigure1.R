#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling
#module load conda_R/4.4

library(SingleCellExperiment)
library(sessioninfo)
library(scater)
library(here)

#load the object
load(here("processed-data","02_build_sce","sce_celltype.rda"),verbose = TRUE)

#Load the cluster colors
load(here("processed-data","cluster_cols_CellType_Final.rda"),verbose = TRUE)

#Reverse the levels because it is opposite order in the violin plot.
sce$CellType.Final <- factor(x = sce$CellType.Final,
                             levels = rev(levels(sce$CellType.Final)))


#Violin plots for # genes, log10(# reads), % mito, and doubletscore
detected_Vln <-  plotColData(object = sce,
                             y = "CellType.Final",
                             x = "detected",
                             colour_by = "CellType.Final") +
  scale_color_manual(values = cluster_cols) +
  theme(legend.position = "none") +
  labs(x = "Cell Type",y = "Number of Genes Detected")

ggsave(filename = here("plots","Plots_for_Supp","Detected_Features_Violin.pdf"),plot = detected_Vln)

#Violin plot for number of total reads
sum_Vln <- plotColData(object = sce,
                       x = "sum",
                       y = "CellType.Final",
                       colour_by = "CellType.Final") +
  scale_color_manual(values = cluster_cols) +
  scale_y_log10() +
  theme(legend.position = "none") +
  labs(x = "Cell Type",y = "log10(Number of Reads)")

ggsave(filename = here("plots","Plots_for_Supp","Number_of_Reads_Violin.pdf"),plot = sum_Vln)

#Violin plot for detected mitochondrial percentage
mito_Vln <- plotColData(object = sce,
                        x = "subsets_Mito_percent",
                        y = "CellType.Final",
                        colour_by = "CellType.Final") +
  scale_color_manual(values = cluster_cols) +
  theme(legend.position = "none") +
  labs(x = "Cell Type",y = "% Mitochondria")

ggsave(filename = here("plots","Plots_for_Supp","Percent_Mito_Violin.pdf"),plot = mito_Vln)

#Violin plot for detected mitochondrial percentage
doublet_Vln <- plotColData(object = sce,
                           x = "doubletScore",
                           y = "CellType.Final",
                           colour_by = "CellType.Final") +
  scale_color_manual(values = cluster_cols) +
  theme(legend.position = "none") +
  labs(x = "Cell Type",y = "doubletScore")

ggsave(filename = here("plots","Plots_for_Supp","doublet_score_Violin.pdf"),plot = doublet_Vln)


##Now make a tSNE and color by sample
sample_cols <- c("#1E88E5","#D81B60","#004D40")
names(sample_cols) <- unique(sce$Sample)
sample_tSNE <- plotReducedDim(object = sce,dimred = "tSNE_mnn_50",color_by = "Sample") +
  scale_color_manual(values = sample_cols)
ggsave(plot = sample_tSNE,filename = here("plots","Plots_for_Supp","Sample_tSNE.pdf"))

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()



