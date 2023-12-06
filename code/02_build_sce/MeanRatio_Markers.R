##Goal: Evaluate the mean ratio method to help identify useful genes for identification of LS broadly
#as well as LS subclusters. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling

library(SingleCellExperiment)
library(DeconvoBuddies)
library(sessioninfo)
library(ggplot2)
library(scater)
library(dplyr)
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

#plot the top 25 marker genes for each cluster. 
for(i in unique(mean_ratios_CellTypeFinal$cellType.target)){
    print(i)
    x <- subset(mean_ratios_CellTypeFinal,subset=(cellType.target == i))
    for(l in 1:25){
        print(l)
        y <- plotReducedDim(object = sce,
                            dimred = "tSNE_mnn_15",
                            colour_by = as.character(x[l,"Symbol"]),
                            swap_rownames = "Symbol") +
            scale_color_gradientn(colours = c("lightgrey","orange","red"))
        ggsave(filename = here("plots",
                               "mean_ratio_plots",
                               "CellType_Final_plots",
                               "Feature_Plots",
                               paste(i,x[l,"Symbol"],l,"FeaturePlot_tSNE.pdf",sep = "_")),
               plot = y,
               height = 8,
               width = 8)
    }
}

#Subset the top 100 by each cluster. 
mean_ratios_CTF_top100 <- subset(mean_ratios_CellTypeFinal,subset=(rank_ratio %in% 1:100))

#Write out the files. 
write.csv(x = mean_ratios_CTF_top100,
          file = here("processed-data","mean_ratios_top100_CellTypeFinal.csv"))
write.csv(x = mean_ratios_CellTypeFinal,
          file = here("processed-data","mean_ratios_all_CellTypeFinal.csv"))

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

write.csv(x = mean_ratios_LS,
          file = here("processed-data","mean_ratios_all_LS-broad.csv"))

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


#Plot the top 25 marker genes by mean ratio method.  
for(i in 1:25){
    print(i)
    y <- plotReducedDim(object = sce,
                        dimred = "tSNE_mnn_15",
                        colour_by = as.character(mean_ratios_LS[i,"Symbol"]),
                        swap_rownames = "Symbol") +
        scale_color_gradientn(colours = c("lightgrey","orange","red"))
    ggsave(filename = here("plots","mean_ratio_plots",
                           "Lateral_Septum_Broad","Feature_Plots",
                           paste0(as.character(mean_ratios_LS[i,"Symbol"]),
                                  "_LS-broad_marker_",
                                  i,
                                  ".pdf")),
           plot = y,
           height = 8,
           width = 8)
}

#Cell ids for LS vs other
cell_ids_ls <- lapply(X = split(x = colData(sce_neuronal)[,c("unique_rowname","Lateral_septum")],
                             f = colData(sce_neuronal)$Lateral_septum),
                   FUN = function(x){x[["unique_rowname"]]})

#Identify the pct of nuclei in which a transcript for each gene is detected. 
pct_detected_ls <- lapply(cell_ids_ls,FUN = function(x){
    as.data.frame(apply(logcounts(sce_neuronal)[,x]>0,MARGIN = 1,FUN = sum)/length(x)*100)
})

mean_ratios_LS <- as.data.frame(mean_ratios_LS)
colnames(mean_ratios_LS)[1] <- "gene_id"

for(i in names(pct_detected_ls)){
    colnames(pct_detected_ls[[i]])[1] <- "Percent_detected"
    pct_detected_ls[[i]]$CellType <- i
    pct_detected_ls[[i]]$gene_id <- row.names(pct_detected_ls[[i]])
    #Add gene symbol info
    pct_detected_ls[[i]] <- left_join(x  = as.data.frame(pct_detected_ls[[i]]),
                                      y  = as.data.frame(rowData(sce)[,c("gene_id","Symbol")]),
                                      by = "gene_id")
    #Add mean ratio information 
    pct_detected_ls[[i]] <- left_join(x  = pct_detected_ls[[i]],
                                      y = subset(mean_ratios_LS,
                                                 subset=(cellType.target == i))[,c("gene_id","mean.target","mean",
                                                                                   "ratio","rank_ratio")],
                                      by = "gene_id")
}

#Save 
save(pct_detected_ls,file = here("processed-data","LSBroad_pctDetected_meanratio.rda"))


#On a per cluster basis, calculate what percentage of cells each gene is detected
#Get cell ids (unique_rowname from colData(sce))
cell_ids <- lapply(X = split(x = colData(sce)[,c("unique_rowname","CellType.Final")],f = colData(sce)$CellType.Final),
                   FUN = function(x){x[["unique_rowname"]]})

#Identify the pct of nuclei in which a transcript for each gene is detected. 
pct_detected <- lapply(cell_ids,FUN = function(x){
   as.data.frame(apply(logcounts(sce)[,x]>0,MARGIN = 1,FUN = sum)/length(x)*100)
})

#Manipulate the list to make a dataframe. 
mean_ratios_CellTypeFinal_df <- as.data.frame(mean_ratios_CellTypeFinal)
colnames(mean_ratios_CellTypeFinal_df)[1] <- "gene_id"

#Add gene symbol information, as well as 1vALL DEG and mean ratio info
for(i in names(pct_detected)){
    colnames(pct_detected[[i]])[1] <- "Percent_detected"
    pct_detected[[i]]$CellType <- i
    pct_detected[[i]]$gene_id <- row.names(pct_detected[[i]])
    #Add gene symbol info
    pct_detected[[i]] <- left_join(x  = as.data.frame(pct_detected[[i]]),
                                   y  = as.data.frame(rowData(sce)[,c("gene_id","Symbol")]),
                                   by = "gene_id")
    #Add mean ratio information 
    pct_detected[[i]] <- left_join(x  = pct_detected[[i]],
                                   y = subset(mean_ratios_CellTypeFinal_df,
                                              subset=(cellType.target == i))[,c("gene_id","mean.target","mean",
                                                                                "ratio","rank_ratio")],
                                   by = "gene_id")
    # #Add 1vALL DEG info
    pct_detected[[i]] <- left_join(x =  pct_detected[[i]],
                                   y = subset(markers_1vALL_enrich_Final,
                                              subset=(cellType.target == i))[,c("gene_id","logFC","std.logFC",
                                                                                "log.FDR","rank_marker","gene_name")],
                                   by = "gene_id")
    
}

#Save the list 
save(pct_detected,file = here("processed-data","Clusterspecific_DEG_pctDetected_meanratio.rda"))

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# [1] "Reproducibility information:"
# [1] "2023-12-06 15:47:04 EST"
# user   system  elapsed 
# 1049.993   15.377 1799.725 
# ─ Session info ────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 Patched (2023-07-19 r84711)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2023-12-06
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/pandoc
# 
# ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────
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
# dplyr                * 1.1.3     2023-09-03 [2] CRAN (R 4.3.1)
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
# plyr                   1.8.8     2022-11-11 [2] CRAN (R 4.3.1)
# purrr                  1.0.2     2023-08-10 [2] CRAN (R 4.3.1)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.1)
# rafalib                1.0.0     2015-08-09 [1] CRAN (R 4.3.1)
# ragg                   1.2.5     2023-01-12 [2] CRAN (R 4.3.1)
# RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.3.1)
# Rcpp                   1.0.11    2023-07-06 [2] CRAN (R 4.3.1)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.1)
# reshape2               1.4.4     2020-04-09 [2] CRAN (R 4.3.1)
# rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.1)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.3.1)
# S4Arrays               1.0.6     2023-08-30 [2] Bioconductor
# S4Vectors            * 0.38.1    2023-05-02 [2] Bioconductor
# ScaledMatrix           1.8.1     2023-05-03 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.3.1)
# scater               * 1.28.0    2023-04-25 [2] Bioconductor
# scran                  1.28.2    2023-07-23 [2] Bioconductor
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
# ───────────────────────────────────────────────────────────────────────────────────────────────────────
