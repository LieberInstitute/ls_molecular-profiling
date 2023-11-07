#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

#Load libraries
library(SingleCellExperiment)
library(sessioninfo)
library(rafalib)
library(ggplot2)
library(scater)
library(dplyr)
library(scran)
library(here)

#Load the object 
load(here("processed-data","sce_postMNN.rda"))

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
# colData names(53): Sample Barcode ... doubletScore sizeFactor
# reducedDimNames(18): GLMPCA_approx UMAP_15 ... tSNE_mnn_25 tSNE_mnn_50
# mainExpName: NULL
# altExpNames(0):

#Lower K= fewer, larger clusters
#Smaller= more, small clusters.
#jaccard + louvain is similar to seurat workflow. 
#Resolution=1 is the default
snn_k_20 <- buildSNNGraph(sce, k = 20, use.dimred = "mnn",type="jaccard")
snn_k_50 <- buildSNNGraph(sce, k = 50, use.dimred = "mnn",type="jaccard")
snn_k_75 <- buildSNNGraph(sce, k = 75, use.dimred = "mnn",type="jaccard")

#Louvain clustering
set.seed(1234)
clust_20 <- igraph::cluster_louvain(snn_k_20,resolution=1)$membership
table(clust_20)
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 264  588  291  197  254 1078  635  544  628  643  228  739  133  357  101  253 
# 17   18   19   20   21   22   23   24 
# 71   49  149  630  121  296  801  175 

set.seed(1234)
clust_50 <- igraph::cluster_louvain(snn_k_50,resolution=1)$membership
table(clust_50)
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 241 1407  361  197  254 1128  649  545  627  643  974  168  211  252  149  615 
# 17 
# 804

set.seed(1234)
clust_75 <- igraph::cluster_louvain(snn_k_75,resolution=1)$membership
table(clust_75)
clust_75
#    1    2    3    4    5    6    7    8    9   10   11   12 
# 1352 1395  347  196  403 1456  545  627  645 1483  162  614 


#Add cluster information to object
sce$k_20_louvain_1 <- factor(clust_20)
sce$k_50_louvain_1 <- factor(clust_50)
sce$k_75_louvain_1 <- factor(clust_75)

#Also run walktrap clustering. 
set.seed(1234)
wt_clusters_k_20 <- igraph::cluster_walktrap(snn_k_20)$membership
table(wt_clusters_k_20)
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 436  243  189  255  645  698  228  420 1018  618  774  253  161  347  381   61 
#  17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32 
# 264  149  172  269  114  297  175   83  206   65  107   58  119   51   73   56 
#  33   34   35   36   37   38 
# 63   19   29   40   42   47 

set.seed(1234)
wt_clusters_k_50 <- igraph::cluster_walktrap(snn_k_50)$membership
table(wt_clusters_k_50)
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 544  344 1420  361 1357  643  197 1399  630  272  254 1093  355  207  149 

set.seed(1234)
wt_clusters_k_75 <- igraph::cluster_walktrap(snn_k_75)$membership
table(wt_clusters_k_75)
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14 
# 545 1390  363 1571  643  625  359 1093 1419  254  251  196  365  151

#Add cluster information to object
sce$k_20_walktrap <- factor(wt_clusters_k_20)
sce$k_50_walktrap <- factor(wt_clusters_k_50)
sce$k_75_walktrap <- factor(wt_clusters_k_75)


#Plot the clusters on the tSNE with 15 dimensions. 
for(i in c("k_20_louvain_1","k_50_louvain_1","k_75_louvain_1",
           "k_20_walktrap","k_50_walktrap","k_75_walktrap")){
    print(i)
    x <- plotReducedDim(object = sce,
                        dimred = "tSNE_mnn_15",
                        colour_by = i,text_by = i) +
        ggtitle(i) +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(plot = x,filename = here("plots","Dim_Red",paste0(i,"_tSNE_mnn_15.pdf")))
}

#Tables by cluster and sample
table(sce$k_20_louvain_1,sce$Sample)
#    1c_LS_SCP 2c_LS_SCP 3c_LS_SCP
# 1        162        92        10
# 2        437        34       117
# 3        291         0         0
# 4         55        79        63
# 5        151        77        26
# 6        544       484        50
# 7         39       376       220
# 8        322       197        25
# 9        270       268        90
# 10       470       102        71
# 11       161        28        39
# 12       224       218       297
# 13        53        30        50
# 14       183         0       174
# 15        59        37         5
# 16        83        58       112
# 17        71         0         0
# 18        39        10         0
# 19        91        26        32
# 20       220       281       129
# 21        91         4        26
# 22        21        94       181
# 23         3       211       587
# 24        23         0       152

table(sce$k_50_louvain_1,sce$Sample)
#    1c_LS_SCP 2c_LS_SCP 3c_LS_SCP
# 1        153        82         6
# 2        663       125       619
# 3        361         0         0
# 4         55        79        63
# 5        151        77        26
# 6        585       493        50
# 7         37       385       227
# 8        323       197        25
# 9        268       269        90
# 10       470       102        71
# 11       337       268       369
# 12       109        34        25
# 13       152        27        32
# 14        82        59       111
# 15        91        26        32
# 16       219       271       125
# 17         7       212       585

table(sce$k_75_louvain_1,sce$Sample)
#    1c_LS_SCP 2c_LS_SCP 3c_LS_SCP
# 1        725       568        59
# 2        663       121       611
# 3        347         0         0
# 4         54        79        63
# 5        242       103        58
# 6         45       599       812
# 7        323       197        25
# 8        268       269        90
# 9        471       103        71
# 10       600       367       516
# 11       106        30        26
# 12       219       270       125

table(sce$k_20_walktrap,sce$Sample)
#     1c_LS_SCP 2c_LS_SCP 3c_LS_SCP
# 1        321        94        21
# 2         75       122        46
# 3        135        37        17
# 4         18        89       148
# 5        268       283        94
# 6        186       213       299
# 7        161        28        39
# 8          1        81       338
# 9         43       506       469
# 10       444        39       135
# 11       372       360        42
# 12        83        58       112
# 13        75        51        35
# 14       176         0       171
# 15       244        87        50
# 16        36        25         0
# 17       227        15        22
# 18        91        26        32
# 19       172         0         0
# 20       167        92        10
# 21         5        54        55
# 22       165       125         7
# 23        24         0       151
# 24        50        25         8
# 25        70        93        43
# 26        16        40         9
# 27         0       103         4
# 28        54         3         1
# 29       119         0         0
# 30        40        11         0
# 31        13        15        45
# 32        38        14         4
# 33        37         1        25
# 34         0         0        19
# 35        29         0         0
# 36        23        12         5
# 37        42         0         0
# 38        43         4         0

table(sce$k_50_walktrap,sce$Sample)
#    1c_LS_SCP 2c_LS_SCP 3c_LS_SCP
# 1        322       197        25
# 2        136       133        75
# 3        665       131       624
# 4        361         0         0
# 5        519       342       496
# 6        470       102        71
# 7         55        79        63
# 8        748       589        62
# 9        271       269        90
# 10        84       138        50
# 11       151        77        26
# 12        40       532       521
# 13         1        64       290
# 14       149        27        31
# 15        91        26        32

table(sce$k_75_walktrap,sce$Sample)
#    1c_LS_SCP 2c_LS_SCP 3c_LS_SCP
# 1        323       197        25
# 2        663       121       606
# 3        142       143        78
# 4        668       367       536
# 5        470       102        71
# 6        267       270        88
# 7        359         0         0
# 8         45       530       518
# 9        751       598        70
# 10       151        77        26
# 11        77       127        47
# 12        54        79        63
# 13         1        68       296
# 14        92        27        32

#Preliminary explorration of distribution of expression of marker genes found that one cluster 
#was dominated by Slc17a7 expression and only contained sample 1. This is not driven by batch effect
#but rather biological (due to anatomy of section).

#save the object
save(sce,file = here("processed-data","sce_clustered.rda"))

#check doublet score per cluster. 
#will move forward with k=20 louvain
doublet_violin <- plotColData(object = sce,
                              x = "k_20_louvain_1",
                              y = "doubletScore",
                              colour_by = "k_20_louvain_1") +
    labs(x = "Cluster",
         y = "Doublet Score",
         title = "Doublet Score by Cluster") +
    theme(plot.title = element_text(hjust=0.5),legend.position = "none") +
    geom_hline(yintercept = 5)
ggsave(doublet_violin,filename = here("plots","doublet_score_by_cluster_k_20_louvain_1_violin.pdf"))

#number of genes per cluster
genes_violin <- plotColData(object = sce,
                            x = "k_20_louvain_1",
                            y = "detected",
                            colour_by = "k_20_louvain_1") +
    labs(x = "Cluster",
         y = "Number of genes/cell",
         title = "Number of Genes/Cell by Cluster") +
    theme(plot.title = element_text(hjust=0.5),legend.position = "none")
ggsave(genes_violin,filename = here("plots","Genes_by_cluster_k_20_louvain_1_violin.pdf"))

#library size per cluster
lib_violin <- plotColData(object = sce,
                          x = "k_20_louvain_1",
                          y = "sum",
                          colour_by = "k_20_louvain_1") +
    scale_y_log10() +
    labs(x = "Cluster",
         title = "Total UMIs") +
    theme(plot.title = element_text(hjust=0.5),legend.position = "none")
ggsave(lib_violin,filename = here("plots","lib_size_by_cluster_k_20_louvain_1_violin.pdf"))

#Dfine some genes that are good markers. 
#Primary goal of making this plot is to make sure that there are no low quality clusters.   
genes <- c("SYT1","SNAP25", #pan neuron
           "MBP","MOBP", #OLIGODENDROCYTE
           "CD74", "CSF1R", "C3", #MICROGLIA
           "GFAP", "TNC", "AQP4", "SLC1A2", #ASTROCYTEs
           "GAD1","GAD2","SLC32A1",#Pan GABA
           "SLC17A7", "SLC17A6", "SLC17A8", #Glutamatergic
           "TRPC4","HOMER2","PTPN3", #Mouse LS markers
           "ELAVL2", #Mouse LS markers
           "CRHR1","CRHR2", 
           "OXTR","AVPR1A", 
           "DRD3",
           "CLDN5", "FLT1", "VTN",#endothelial
           "COL1A2", "TBX18", "RBPMS",#Mural
           "SKAP1", "ITK", "CD247", #Tcell
           "CD163", "SIGLEC1", "F13A1",#Macrophage
           "PDGFRA", "VCAN", "CSPG4")#Polydendrocytes

#Walktrap
Expression_dotplot <- plotDots(object = sce,
                               features = rev(genes),
                               group = "k_20_louvain_1",swap_rownames = "gene_name") +
    scale_color_gradientn(colours = c("lightgrey","orange","red"))
ggsave(plot = Expression_dotplot,filename = here("plots","Expression_plots",
                                                 "post_k_20_louvain_1_clustering_general_dotplot.pdf"),
       height = 8)

##########This commented code from 10/20/23#########
# #The louvain algorithm identifies a cluster (cluster 11) that exhibits low # of genes/cell, but expresses
# #only neuronal genes. I am pretty sure these are low quality nuclei. The walktrap algorithm just lumps 
# #these cells into another cluster. We are going to remove them. 
# #First what is the sample make up of this cluster. 
# table(sce$Sample,sce$k_50_louvain_1 == 11)
# #           FALSE TRUE
# # 1c_LS_SCP  4063  269
# # 2c_LS_SCP  2706  218
# # 3c_LS_SCP  2456   95
# #Primarily coming from samples 1 and 2, but all samples have cells within that cluster. 
# #Will need to remove the cluster and rerun dimensionality reduction steps. 
# 
# #Identify the cell IDs that make up cluster 11 so that they can be removed on the front end of the analysis. 
# low_quality_nuclei <- colnames(sce[,sce$k_50_louvain_1 == 11])
# save(low_quality_nuclei,file = here("processed-data","cluster_11_low_quality_IDs.rda"))
##########################################
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# [1] "Reproducibility information:"
# [1] "2023-11-07 12:18:40 EST"
# user   system  elapsed 
# 870.180   10.914 1909.545 
# ─ Session info ──────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 Patched (2023-07-19 r84711)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2023-11-07
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/pandoc
# 
# ─ Packages ──────────────────────────────────────────────────────────────────
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
# colorout             * 1.2-2     2023-09-22 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.3.1)
# cowplot                1.1.1     2020-12-30 [1] CRAN (R 4.3.1)
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.1)
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
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.1)
# rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 4.3.1)
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
# scran                * 1.28.2    2023-07-23 [2] Bioconductor
# scuttle              * 1.10.2    2023-08-03 [2] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.1)
# SingleCellExperiment * 1.22.0    2023-04-25 [2] Bioconductor
# sparseMatrixStats      1.12.2    2023-07-02 [2] Bioconductor
# statmod                1.5.0     2023-01-06 [2] CRAN (R 4.3.1)
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
# ─────────────────────────────────────────────────────────────────────────────
