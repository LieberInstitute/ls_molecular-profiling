##Goal: Evaluate the mean ratio method to help identify useful genes for identification of LS broadly
#as well as LS subclusters. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling

library(SingleCellExperiment)
library(DeconvoBuddies)
library(sessioninfo)
library(ggplot2)
library(scater)
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

#load the DEG lists. 
load(here("processed-data","markers_1vAll_ttest_CellTypeFinal_20Clusters.rda"),verbose = TRUE)
# Loading objects:
#     markers_1vALL_enrich_Final
load(here("processed-data","markers_pairwise_list_CellTypeFinal_20CellTypes.rda"),verbose = TRUE) 
# Loading objects:
#     markers_pairwise


#Now use the mean ratio method to help identify better markers. Function = get_mean_ratio2() from DeconvoBuddies
#From http://research.libd.org/DeconvoBuddies/articles/DeconvoBuddies.html#using-meanratio-to-find-cell-type-markers: 
# "To select genes specific for each cell type, you can evaluate the mean ratio for each gene x each cell type, where 
# mean ratio = mean(Expression of target cell type)/mean(Expression of highest non-target cell type)"

#Source code for the get_mean_ratio2 function shows that the gene symbol column needs to be called "Symbol" 
#Currently called "gene_name" so change it. 
colnames(rowData(sce))[5] <- "Symbol"

#Calculate the mean ratio. 
mean_ratios_CellTypeFinal <- get_mean_ratio2(sce,
                                             cellType_col = "CellType.Final",
                                             add_symbol = TRUE)

#Check out structure and first couple lines of the dataframe. 
dim(mean_ratios_CellTypeFinal)
# [1] 89642     9

mean_ratios_CellTypeFinal[1:5,]
# # A tibble: 5 × 9
# gene        cellType.target mean.target cellType  mean ratio rank_ratio Symbol
# <chr>       <chr>                 <dbl> <chr>    <dbl> <dbl>      <int> <chr> 
#     1 ENSG000000… LS_Inh_A              0.807 LS_Inh_G 0.372  2.17          1 SLC12…
# 2 ENSG000001… LS_Inh_A              2.00  LS_Inh_B 1.09   1.83          2 SLC27…
# 3 ENSG000000… LS_Inh_A              1.06  Excit_A  0.593  1.78          3 CROT  
# 4 ENSG000001… LS_Inh_A              4.21  MS_Inh_E 2.38   1.77          4 COL25…
# 5 ENSG000002… LS_Inh_A              1.87  MS_Inh_A 1.07   1.75          5 EPHA5…
# # ℹ 1 more variable: anno_ratio <chr>

#Plot the top 10 markers of each cluster. 
for(i in unique(mean_ratios_CellTypeFinal$cellType.target)){
    print(i)
    top_plot <- plot_marker_express(sce       = sce,
                                    stats     = mean_ratios_CellTypeFinal,
                                    cell_type = i,
                                    n_genes   = 10,
                                    rank_col  = "rank_ratio",
                                    anno_col  = "anno_ratio",
                                    cellType_col = "CellType.Final",
                                    color_pal = new_cluster_cols)
    ggsave(plot = top_plot,
           filename = here("plots",
                           "mean_ratio_plots",
                           "CellType_Final_plots",
                           paste0(i,"_Top10_meanratio.pdf")
                           ),
           height = 12,
           width = 8)
}

#plot the top 10 marker genes for each cluster. 
for(i in unique(mean_ratios_CellTypeFinal$cellType.target)){
    x <- subset(mean_ratios_CellTypeFinal,subset=(cellType.target == i))
    for(l in 1:10){
        y <- plotReducedDim(object = sce,
                            dimred = "tSNE_mnn_15",
                            colour_by = as.character(x[l,"Symbol"]),
                            swap_rownames = "Symbol") +
            scale_color_gradientn(colours = c("lightgrey","orange","red"))
        ggsave(filename = here("plots",
                               "mean_ratio_plots",
                               "CellType_Final_plots",
                               "Feature_Plots",
                               paste(i,x[l,"Symbol"],"FeaturePlot_tSNE.pdf",sep = "_")),
               plot = y,
               height = 8,
               width = 8)
    }
}

#Subset the top 100 by each cluster. 
mean_ratios_CTF_top100 <- subset(mean_ratios_CellTypeFinal,subset=(rank_ratio %in% 1:100))

#Write out the file. 
write.csv(x = mean_ratios_CTF_top100,
          file = here("processed-data","mean_ratios_top100_CellTypeFinal.csv"))



#Use the mean ratio method to identify markers of region specific neuronal populations. 
#To do this first subset to neuronal clusters only. 
#Going to leave out Excit_A and Excit_B because I am not totally sure where they are coming from. 
sce_neuronal <- sce[,sce$CellType.Final %in% c("LS_Inh_A","LS_Inh_B","LS_Inh_G","LS_Inh_I",
                                               "MS_Inh_A","MS_Inh_E","MS_Inh_H","MS_Excit_A",
                                               "Sept_Inh_D","Sept_Inh_F","Str_Inh_A","Str_Inh_B",
                                               "Excit_A","Excit_B")]

#Make several columns to facilitate region specific vs all others. 
#Lateral Septum
sce_neuronal$Lateral_septum <- ifelse(sce_neuronal$CellType.Final %in% c("LS_Inh_A","LS_Inh_B","LS_Inh_G","LS_Inh_I"),
                                      "LS",
                                      "Other")

#Calculate the mean ratio. Focus on lateral septum 
mean_ratios_LS <- get_mean_ratio2(sce_neuronal,
                                  cellType_col = "Lateral_septum",
                                  add_symbol = TRUE)

#Top 100 LS genes
LS_broad_top100 <- subset(mean_ratios_LS,subset=(cellType.target == "LS" & rank_ratio %in% 1:100))

#Write out the file. 
write.csv(x = LS_broad_top100,
          file = here("processed-data","mean_ratios_top100_LS-broad.csv"))

#Plot top 10 broad LS
top_LS_plot <- plot_marker_express(sce          = sce_neuronal,
                                   stats        = LS_broad_top100,
                                   cell_type    = "LS",
                                   n_genes      = 10,
                                   rank_col     = "rank_ratio",
                                   anno_col     = "anno_ratio",
                                   cellType_col = "Lateral_septum")

ggsave(plot = top_LS_plot,
       filename = here("plots",
                       "mean_ratio_plots",
                       "Lateral_Septum_Broad",
                       "LS-Broad_Top10_meanratio.pdf"),
       height = 12,
       width = 8)


#plot the top 10 marker genes for each cluster. 
for(i in unique(mean_ratios_CellTypeFinal$cellType.target)){
    x <- subset(mean_ratios_CellTypeFinal,subset=(cellType.target == i))
    for(l in 1:10){
        y <- plotReducedDim(object = sce,
                            dimred = "tSNE_mnn_15",
                            colour_by = as.character(x[l,"Symbol"]),
                            swap_rownames = "Symbol") +
            scale_color_gradientn(colours = c("lightgrey","orange","red"))
        ggsave(filename = here("plots",
                               "mean_ratio_plots",
                               "CellType_Final_plots",
                               "Feature_Plots",
                               paste(i,x[l,"Symbol"],"FeaturePlot_tSNE.pdf",sep = "_")),
               plot = y,
               height = 8,
               width = 8)
    }
}


