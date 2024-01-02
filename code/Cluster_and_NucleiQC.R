#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

###Make plots for supp figures. 
library(SingleCellExperiment)
library(ggplot2)
library(scater)
library(here)



#load the SingleCellExperiment object for human Lateral Septum
load(here("processed-data","sce_with_CellType.rda"))


sce
# class: SingleCellExperiment 
# dim: 33556 9225 
# metadata(1): Samples
# assays(3): counts binomial_deviance_residuals logcounts
# rownames(33556): ENSG00000243485 ENSG00000238009 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(9225): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 3_TTTGGTTTCTTCGACC-1 3_TTTGTTGTCCCGATCT-1
# colData names(61): Sample Barcode ... CellType_k_20_louvain
# CellType.Final
# reducedDimNames(18): GLMPCA_approx UMAP_15 ... tSNE_mnn_25 tSNE_mnn_50
# mainExpName: NULL
# altExpNames(0):


#Load colors
load(here("processed-data","Final_CellTypes_colors.rda"),verbose = TRUE)


#Make CellType.Final a factor
sce$CellType.Final <- factor(x = sce$CellType.Final,
                             levels = rev(c("LS_Inh_A","LS_Inh_B","LS_Inh_G","LS_Inh_I",
                                            "MS_Inh_A","MS_Inh_E","MS_Inh_H",
                                            "Sept_Inh_D","Sept_Inh_F",
                                            "MS_Excit_A","Excit_A","Excit_B",
                                            "Str_Drd1-MSN","Str_Drd1-Patch","Str_Drd1-Matrix","Str_Drd2-MSN",
                                            "Astrocyte","Ependymal","Microglia",
                                            "Mural","Oligo","Polydendrocyte")))

#Violin plot for number of detected genes. 
detected_Vln <- plotColData(object = sce,
                            x = "detected",y = "CellType.Final",colour_by = "CellType.Final") +
    scale_color_manual(values = cluster_cols) +
    theme(legend.position = "none") +
    labs(y = "Cell Type",x = "Number of Genes Detected")

ggsave(filename = here("plots","Plots_for_Supp","Detected_Features_Violin.pdf"),plot = detected_Vln)

#Violin plot for number of total reads
sum_Vln <- plotColData(object = sce,
                       x = "sum",y = "CellType.Final",colour_by = "CellType.Final") +
    scale_color_manual(values = cluster_cols) +
    scale_y_log10() +
    theme(legend.position = "none") +
    labs(x = "Cell Type",y = "log10(Number of Reads)")

ggsave(filename = here("plots","Plots_for_Supp","Number_of_Reads_Violin.pdf"),plot = sum_Vln)

#Violin plot for detected mitochondrial percentage
mito_Vln <- plotColData(object = sce,
                        x = "subsets_Mito_percent",y = "CellType.Final",colour_by = "CellType.Final") +
    scale_color_manual(values = cluster_cols) +
    theme(legend.position = "none") +
    labs(x = "Cell Type",y = "% Mitochondria")

ggsave(filename = here("plots","Plots_for_Supp","Number_of_Reads_Violin.pdf"),plot = mito_Vln)

#Violin plot for detected mitochondrial percentage
doublet_Vln <- plotColData(object = sce,
                           x = "doubletScore",y = "CellType.Final",colour_by = "CellType.Final") +
    scale_color_manual(values = cluster_cols) +
    theme(legend.position = "none") +
    labs(x = "Cell Type",y = "doubletScore")

ggsave(filename = here("plots","Plots_for_Supp","doublet_score_Violin.pdf"),plot = mito_Vln)


#Make umaps where points are colored by donor. 
Br8331 <- plotReducedDim(object = sce[,sce$Brain == "Br8331"],dimred = "tSNE_mnn_50",color_by = "Brain") +
    scale_color_manual(values = "#1E88E5")
Br8354 <- plotReducedDim(object = sce[,sce$Brain == "Br8354"],dimred = "tSNE_mnn_50",color_by = "Brain") +
    scale_color_manual(values = "#D81B60")
Br9103 <- plotReducedDim(object = sce[,sce$Brain == "Br8354"],dimred = "tSNE_mnn_50",color_by = "Brain") +
    scale_color_manual(values = "#004D40")

ggsave(filename = here("plots","Plots_for_Supp","8331_tSNE.pdf"),plot = Br8331)
ggsave(filename = here("plots","Plots_for_Supp","8354_tSNE.pdf"),plot = Br8354)
ggsave(filename = here("plots","Plots_for_Supp","9103_tSNE.pdf"),plot = Br9103)




