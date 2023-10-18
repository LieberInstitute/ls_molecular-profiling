#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

#Load libraries
library(SingleCellExperiment)
library(sessioninfo)
library(rafalib)
library(ggplot2)
library(scater)
library(scran)
library(here)

#Load the object 
load(here("processed-data","sce_numeric_cutoff_QC.rda"))

sce
# class: SingleCellExperiment 
# dim: 36601 10000 
# metadata(1): Samples
# assays(3): counts binomial_deviance_residuals logcounts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(10000): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 3_TTTGGTTTCTTCGACC-1 3_TTTGTTGTCCCGATCT-1
# colData names(55): Sample Barcode ... discard_numeric sizeFactor
# reducedDimNames(14): GLMPCA_approx UMAP_15 ... UMAP_mnn_25 UMAP_mnn_50
# mainExpName: NULL
# altExpNames(0):

#Build SNN graph with k=20 and k=50. Can change k value later but will start with this. 
#Lower K= fewer, larger clusters
#Smaller= more, small clusters.
#Run louvain clustering. 
#jaccard + louvain is similar to seurat workflow. 
#Resolution=1 is the default method
snn_k_20 <- buildSNNGraph(sce, k = 20, use.dimred = "mnn",type="jaccard")
snn_k_50 <- buildSNNGraph(sce, k = 50, use.dimred = "mnn",type="jaccard")
snn_k_75 <- buildSNNGraph(sce, k = 75, use.dimred = "mnn",type="jaccard")

set.seed(1234)
clust_20 <- igraph::cluster_louvain(snn_k_20,resolution=1)$membership
table(clust_20)
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 276 1125  307  200  258 1096  475  542  674  645  225  838  143  106  517  252 
# 17   18   19   20   21   22   23 
# 72   49  168  637  126  310  959

set.seed(1234)
clust_50 <- igraph::cluster_louvain(snn_k_50,resolution=1)$membership
table(clust_50)
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 1430 1125  362  197  257 1446  545  657  643  226  981  350  694  154  637  296 

set.seed(1234)
clust_75 <- igraph::cluster_louvain(snn_k_75,resolution=1)$membership
table(clust_75)
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 1365 1224  355  193  256 1456  545  650  645 1382  750  214  143  636  186 

#Add cluster information to object
sce$k_20_louvain_1 <- factor(clust_20)
sce$k_50_louvain_1 <- factor(clust_50)
sce$k_75_louvain_1 <- factor(clust_75)

#Also run walktrap clustering. 
set.seed(1234)
wt_clusters_k_20 <- igraph::cluster_walktrap(snn_k_20)$membership
table(wt_clusters_k_20)
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 258 1087  721  436  375  192  218 1215  807  854  253  338  157   61  379  286 
# 17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32 
# 307  265  164   86  114  179  170  143  144   58  107   95   77  119   48   54 
# 33   34   35   36   37 
# 56   63   29   43   42 

set.seed(1234)
wt_clusters_k_50 <- igraph::cluster_walktrap(snn_k_50)$membership
table(wt_clusters_k_50)
wt_clusters_k_50
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 289  163 1866 1228  643 1417  257  197 1381  373  316  353  215  332  154  149 
#  17   18   19   20 
# 228  271   98   70 

set.seed(1234)
wt_clusters_k_75 <- igraph::cluster_walktrap(snn_k_75)$membership
table(wt_clusters_k_75)
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14 
# 625  545 1213 2344  646  642  343 1360 1081  255  197  377  230  142 

#Add cluster information to object
sce$k_20_walktrap <- factor(wt_clusters_k_20)
sce$k_50_walktrap <- factor(wt_clusters_k_50)
sce$k_75_walktrap <- factor(wt_clusters_k_75)

#Plot the clusters on the umap
#k=20 louvain
k_20_umap <- plotReducedDim(object = sce,dimred = "UMAP_mnn_20",
                            colour_by = "k_20_louvain_1",text_by = "k_20_louvain_1") #Color and label by cluster 
ggsave(plot = k_20_umap,filename = here("plots","Dim_Red","k_20_louvain_umap_20components.pdf"))

#k=50 louvain
k_50_umap <- plotReducedDim(object = sce,dimred = "UMAP_mnn_20",
                            colour_by = "k_50_louvain_1",text_by = "k_50_louvain_1") #Color and label by cluster 
ggsave(plot = k_50_umap,filename = here("plots","Dim_Red","k_50_louvain_umap_20components.pdf"))

#k=75 louvain
k_75_umap <- plotReducedDim(object = sce,dimred = "UMAP_mnn_20",
                            colour_by = "k_75_louvain_1",text_by = "k_75_louvain_1") #Color and label by cluster 
ggsave(plot = k_50_umap,filename = here("plots","Dim_Red","k_75_louvain_umap_20components.pdf"))

#k=20 walktrap
k_20_wt_umap <- plotReducedDim(object = sce,dimred = "UMAP_mnn_20",
                               colour_by = "k_20_walktrap",text_by = "k_20_walktrap") #Color and label by cluster 
ggsave(plot = k_20_wt_umap,filename = here("plots","Dim_Red","k_20_walktrap_umap_20components.pdf"))

#k=50 walktrap
k_50_wt_umap <- plotReducedDim(object = sce,dimred = "UMAP_mnn_20",
                               colour_by = "k_50_walktrap",text_by = "k_50_walktrap") #Color and label by cluster 
ggsave(plot = k_50_wt_umap,filename = here("plots","Dim_Red","k_50_walktrap_umap_20components.pdf"))

#k=50 walktrap
k_50_wt_umap <- plotReducedDim(object = sce,dimred = "UMAP_mnn_20",
                               colour_by = "k_50_walktrap",text_by = "k_50_walktrap") #Color and label by cluster 
ggsave(plot = k_50_wt_umap,filename = here("plots","Dim_Red","k_50_walktrap_umap_20components.pdf"))

#k=75 walktrap
k_75_wt_umap <- plotReducedDim(object = sce,dimred = "UMAP_mnn_20",
                               colour_by = "k_75_walktrap",text_by = "k_75_walktrap") #Color and label by cluster 
ggsave(plot = k_75_wt_umap,filename = here("plots","Dim_Red","k_75_wt_umap_umap_20components.pdf"))

#Tables by cluster and sample
table(sce$k_20_louvain_1,sce$Sample)
#     1c_LS_SCP 2c_LS_SCP 3c_LS_SCP
# 1        167        96        13
# 2        649        33       443
# 3        307         0         0
# 4         57        80        63
# 5        152        75        31
# 6        550       493        53
# 7         38       284       153
# 8        321       196        25
# 9        295       281        98
# 10       472       102        71
# 11       158        29        38
# 12       259       235       344
# 13        59        31        53
# 14        62        39         5
# 15       257       131       129
# 16        83        57       112
# 17        72         0         0
# 18        39        10         0
# 19       109        27        32
# 20       234       263       140
# 21        96         4        26
# 22        20       103       187
# 23         5       301       653

table(sce$k_50_louvain_1,sce$Sample)
#    1c_LS_SCP 2c_LS_SCP 3c_LS_SCP
# 1        753       608        69
# 2        646        35       444
# 3        362         0         0
# 4         55        79        63
# 5        151        75        31
# 6         41       595       810
# 7        323       197        25
# 8        297       267        93
# 9        470       102        71
# 10       161        28        37
# 11       310       256       415
# 12       198        71        81
# 13       350       168       176
# 14        96        26        32
# 15       229       268       140
# 16        19        95       182

table(sce$k_75_louvain_1,sce$Sample)
#    1c_LS_SCP 2c_LS_SCP 3c_LS_SCP
# 1        738       571        56
# 2        660        82       482
# 3        355         0         0
# 4         52        78        63
# 5        150        75        31
# 6         47       597       812
# 7        323       197        25
# 8        292       265        93
# 9        471       103        71
# 10       526       347       509
# 11       375       190       185
# 12       155        27        32
# 13        86        26        31
# 14       228       268       140
# 15         3        44       139

table(sce$k_20_walktrap,sce$Sample)
#     1c_LS_SCP 2c_LS_SCP 3c_LS_SCP
# 1        152        75        31
# 2        544       492        51
# 3        468        88       165
# 4        321        94        21
# 5        132       171        72
# 6          2        45       145
# 7        154        27        37
# 8        465       347       403
# 9          5       288       514
# 10       151       476       227
# 11        83        57       113
# 12       175         0       163
# 13        75        47        35
# 14        36        25         0
# 15       242        88        49
# 16       174        99        13
# 17       183        95        29
# 18       229        14        22
# 19       106        26        32
# 20        52        26         8
# 21         5        54        55
# 22        24         0       155
# 23       170         0         0
# 24         0         9       134
# 25        48        19        77
# 26        55         3         0
# 27         0       103         4
# 28        20        43        32
# 29        13        16        48
# 30       119         0         0
# 31        38        10         0
# 32        50         4         0
# 33        38        14         4
# 34        37         1        25
# 35        29         0         0
# 36        24        14         5 

table(sce$k_50_walktrap,sce$Sample)
#     1c_LS_SCP 2c_LS_SCP 3c_LS_SCP
# 1        289         0         0
# 2         95        38        30
# 3        719       466       681
# 4        661        85       482
# 5        470       102        71
# 6        137       692       588
# 7        151        75        31
# 8         55        79        63
# 9        738       585        58
# 10       124       176        73
# 11       126       170        20
# 12       201       115        37
# 13       155        27        33
# 14         0        54       278
# 15        73        47        34
# 16        92        26        31
# 17       196        27         5
# 18        88        62       121
# 19        21        44        33
# 20        70         0         0

table(sce$k_75_walktrap,sce$Sample)
#     1c_LS_SCP 2c_LS_SCP 3c_LS_SCP
# 1        220       267       138
# 2        323       197        25
# 3        658        79       476
# 4        922       583       839
# 5        472       103        71
# 6        289       263        90
# 7        343         0         0
# 8        736       572        52
# 9         47       534       500
# 10       149        75        31
# 11        55        79        63
# 12         0        64       313
# 13       162        28        40
# 14        85        26        31 

#Preliminary explorration of distribution of expression of marker genes found that one cluster 
#was dominated by Slc17a7 expression and only contained sample 1. This is not driven by batch effect
#but rather biological (due to anatomy of section?). k=20 seems to contain several clusters dominated by one 
#sample, while k=50 only has one cluster that is dominated by sample 1 (and is the cluster identified yesterday).
#k=75 looks best here. 

#save the object
save(sce,file = here("processed-data","sce_clustered.rda"))

#check doublet score per cluster. 
#violin
doublet_violin <- plotColData(object = sce,
                              x = "k_50_louvain_1",
                              y = "doubletScore",
                              colour_by = "k_50_louvain_1") +
    labs(x = "Cluster \n(jaccard + louvain k=50 + res=1)",
         y = "Doublet Score",
         title = "Doublet Score by Cluster") +
    theme(plot.title = element_text(hjust=0.5),legend.position = "none") +
    geom_hline(yintercept = 5)
ggsave(doublet_violin,filename = here("plots","doublet_score_by_cluster_k50_louvain_violin.pdf"))

#number of genes per cluster
genes_violin <- plotColData(object = sce,
                              x = "k_50_louvain_1",
                              y = "detected",
                              colour_by = "k_50_louvain_1") +
    labs(x = "Cluster \n(jaccard + louvain k=50 + res=1)",
         y = "Number of genes/cell",
         title = "Number of Genes/Cell by Cluster") +
    theme(plot.title = element_text(hjust=0.5),legend.position = "none")
ggsave(genes_violin,filename = here("plots","Genes_by_cluster_k50_louvain_violin.pdf"))


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


Expression_dotplot <- plotDots(object = sce,
                               features = rev(genes),
                               group = "k_50_louvain_1",swap_rownames = "gene_name") +
    scale_color_gradientn(colours = c("lightgrey","orange","red"))
ggsave(plot = Expression_dotplot,filename = here("plots","Expression_plots",
                                                 "post_k_50_louvain_1_clustering","general_dotplot.pdf"),
       height = 8)

#cluster 13 expresses Gad1/2 but not much else. 
#Will run some basic DEG testing to try and figure out what it is.
colLabels(sce) <- sce$k_50_louvain_1
markers_1vALL <- scran::findMarkers(sce,pval.type = "all", direction = "up",test.type="t")
#Pull out cluster 13 
cluster_13 <- as.data.frame(subset(markers_1vALL[[13]],subset=(FDR <= 0.05)))
cluster_13$gene_id <- row.names(cluster_13)
#add in common gene name info
cluster_13 <- dplyr::left_join(x = cluster_13,y = as.data.frame(rowData(sce)),by = "gene_id")
nrow(cluster_13)
#[1] 4
#4 DEGs. 
#DEGs are SNHG25, UCHL1, KIF2C, LENG8

#Annotate clusters and remake the Expression dotplot. 
annotation_df <- data.frame(cluster = c(1:16))
annotation_df$CellType <- c("LS_GABA_1","LS_GABA_Drd3","Glutamatergic","Polydendrocyte",
                            "Microglia","Oligodendrocyte_1","Astrocyte_1","Oligodendrocyte_2",
                            "LS_GABA_2","MS_GABA_1","MS_GABA_2","MS_GABA_3",
                            "undefined","Mural","Astrocyte_2","LS_GABA_3")

sce$CellType <- annotation_df$CellType[match(sce$k_50_louvain_1,annotation_df$cluster)]

sce$CellType <- factor(x = sce$CellType,
                       levels = c("LS_GABA_1","LS_GABA_2","LS_GABA_3","LS_GABA_Drd3",
                                  "MS_GABA_1","MS_GABA_2","MS_GABA_3","Polydendrocyte",
                                  "Glutamatergic","Astrocyte_1","Astrocyte_2","Oligodendrocyte_1",
                                  "Oligodendrocyte_2","Microglia","Mural","undefined"))

Expression_dotplot <- plotDots(object = sce,
                               features = rev(c("SYT1","SNAP25",
                                                "GAD1","GAD2","SLC32A1",
                                                "TRPC4","HOMER2","PTPN3","DRD3",
                                                "CRHR1","CRHR2","OXTR","AVPR1A",
                                                "ELAVL2",
                                                "SLC17A7", "SLC17A6", "SLC17A8",
                                                "GFAP", "TNC", "AQP4", "SLC1A2",
                                                "MBP","MOBP",
                                                "CD74", "CSF1R", "C3",
                                                "COL1A2", "TBX18", "RBPMS",
                                                "PDGFRA", "VCAN", "CSPG4",
                                                "SKAP1", "ITK", "CD247",
                                                "CD163", "SIGLEC1", "F13A1")),
                               group = "CellType",swap_rownames = "gene_name") +
    theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
    scale_color_gradientn(colours = c("lightgrey","orange","red")) 
ggsave(plot = Expression_dotplot,filename = here("plots","Expression_plots",
                                                 "post_k_50_louvain_1_clustering","general_dotplot_annotated.pdf"),height = 8)


#Make umap annotated
#k=50 louvain
annotated_umap <- plotReducedDim(object = sce,dimred = "UMAP_mnn",
                                 colour_by = "CellType",text_by = "CellType")
ggsave(plot = annotated_umap,filename = here("plots","Dim_Red","Annotated_k_50_louvain_umap_50components.pdf"))

#Remake violin plots with cluster information on the x axis
#violin
doublet_violin <- plotColData(object = sce,
                              x = "CellType",
                              y = "doubletScore",
                              colour_by = "CellType") +
    labs(x = "CellType",
         y = "Doublet Score",
         title = "Doublet Score by Cluster") +
    theme(plot.title = element_text(hjust=0.5),
          legend.position = "none",
          axis.text.x = element_text(angle=45,hjust=1)) +
    geom_hline(yintercept = 5)
ggsave(doublet_violin,filename = here("plots","Annotated_doublet_score_by_cluster_k50_louvain_violin.pdf"))

#number of genes per cluster
genes_violin <- plotColData(object = sce,
                            x = "CellType",
                            y = "detected",
                            colour_by = "CellType") +
    labs(x = "CellType",
         y = "Number of genes/cell",
         title = "Number of Genes/Cell by Cluster") +
    theme(plot.title = element_text(hjust=0.5),
          legend.position = "none",
          axis.text.x = element_text(angle=45,hjust=1))
ggsave(genes_violin,filename = here("plots","Annotated_Genes_by_cluster_k50_louvain_violin.pdf"))

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# [1] "Reproducibility information:"
# [1] "2023-10-18 09:05:14 EDT"
# user   system  elapsed 
# 390.223    4.851 1847.954 
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
# date     2023-10-18
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
# colorout             * 1.2-2     2023-09-22 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.3.1)
# cowplot                1.1.1     2020-12-30 [2] CRAN (R 4.3.1)
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.1)
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
# ───────────────────────────────────────────────────────────────────────────────────────────────────────
# 
