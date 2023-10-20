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
# dim: 36601 9807 
# metadata(1): Samples
# assays(3): counts binomial_deviance_residuals logcounts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(9807): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 3_TTTGGTTTCTTCGACC-1 3_TTTGTTGTCCCGATCT-1
# colData names(53): Sample Barcode ... doubletScore sizeFactor
# reducedDimNames(14): GLMPCA_approx UMAP_15 ... UMAP_mnn_25 UMAP_mnn_50
# mainExpName: NULL
# altExpNames(0):

#Lower K= fewer, larger clusters
#Smaller= more, small clusters.
#Run louvain clustering. 
#jaccard + louvain is similar to seurat workflow. 
#Resolution=1 is the default
snn_k_20 <- buildSNNGraph(sce, k = 20, use.dimred = "mnn",type="jaccard")
snn_k_50 <- buildSNNGraph(sce, k = 50, use.dimred = "mnn",type="jaccard")
snn_k_75 <- buildSNNGraph(sce, k = 75, use.dimred = "mnn",type="jaccard")

set.seed(1234)
clust_20 <- igraph::cluster_louvain(snn_k_20,resolution=1)$membership
table(clust_20)
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 276  773  300  200  256 1092  694  543  645  351  524  284  168  416  139  105 
# 17   18   19   20   21   22   23   24   25 
# 478  255   71   49  162  635  126  313  952 

set.seed(1234)
clust_50 <- igraph::cluster_louvain(snn_k_50,resolution=1)$membership
table(clust_50)
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 1398 1124  361  197  254 1174  545  644  226  991  582  349  149  615  294  904 

set.seed(1234)
clust_75 <- igraph::cluster_louvain(snn_k_75,resolution=1)$membership
table(clust_75)
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14 
# 1389 1224  355  193  392 1455  545  626  644  218 1291  675  614  186

#Add cluster information to object
sce$k_20_louvain_1 <- factor(clust_20)
sce$k_50_louvain_1 <- factor(clust_50)
sce$k_75_louvain_1 <- factor(clust_75)

#Also run walktrap clustering. 
set.seed(1234)
wt_clusters_k_20 <- igraph::cluster_walktrap(snn_k_20)$membership
table(wt_clusters_k_20)
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 436  211  236  237 1063  668  777  750  251  338  162  734   61  375  314  275 
# 17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32 
# 269   95  178  303  160  176  105  284   63  210  144  188  107   45   77  113 
# 33   34   35   36   37   38   39   40   41 
# 49   54   63   51   51   29   44   19   42 

set.seed(1234)
wt_clusters_k_50 <- igraph::cluster_walktrap(snn_k_50)$membership
table(wt_clusters_k_50)
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 643  197  255 1387  662 1495  634  763  328  316 1056  389  227  158  100  127 
# 17   18   19   20   21   22   23   24 
# 228  174  144  232   70   44  115   63 

set.seed(1234)
wt_clusters_k_75 <- igraph::cluster_walktrap(snn_k_75)$membership
table(wt_clusters_k_75)
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14 
# 612  545 1214 2201  730  646  287 1411 1019  251  197  331  226  137 

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
# 1        168        95        13
# 2        469        33       271
# 3        300         0         0
# 4         57        80        63
# 5        151        79        26
# 6        551       491        50
# 7         98       397       199
# 8        321       197        25
# 9        472       102        71
# 10       179         0       172
# 11       218       111       195
# 12        72        87       125
# 13        98        47        23
# 14       208       163        45
# 15        55        31        53
# 16        61        39         5
# 17       209       192        77
# 18        83        58       114
# 19        71         0         0
# 20        39        10         0
# 21       104        26        32
# 22       224       282       129
# 23        96         4        26
# 24        21       105       187
# 25         7       295       650

table(sce$k_50_louvain_1,sce$Sample)
#    1c_LS_SCP 2c_LS_SCP 3c_LS_SCP
# 1        747       590        61
# 2        644        35       445
# 3        361         0         0
# 4         55        79        63
# 5        151        77        26
# 6        309       596       269
# 7        323       197        25
# 8        470       103        71
# 9        161        28        37
# 10       312       268       411
# 11       269       218        95
# 12       197        71        81
# 13        91        26        32
# 14       219       271       125
# 15        20        95       179
# 16         3       270       631

table(sce$k_75_louvain_1,sce$Sample)
#    1c_LS_SCP 2c_LS_SCP 3c_LS_SCP
# 1        740       588        61
# 2        660        82       482
# 3        355         0         0
# 4         52        78        63
# 5        233       102        57
# 6         46       597       812
# 7        323       197        25
# 8        267       270        89
# 9        470       103        71
# 10       158        27        33
# 11       498       317       476
# 12       307       249       119
# 13       220       269       125
# 14         3        45       138

table(sce$k_20_walktrap,sce$Sample)
#    1c_LS_SCP 2c_LS_SCP 3c_LS_SCP
# 1        321        94        21
# 2        144        48        19
# 3         73       118        45
# 4        165        30        42
# 5        374       385       304
# 6        469        50       149
# 7        374       361        42
# 8        133       423       194
# 9         83        57       111
# 10       173         0       165
# 11        75        51        36
# 12        21       305       408
# 13        36        25         0
# 14       241        86        48
# 15         0        42       272
# 16       166        96        13
# 17       230        16        23
# 18        54        31        10
# 19        24         0       154
# 20       164       130         9
# 21       102        26        32
# 22       176         0         0
# 23         3        49        53
# 24       159       100        25
# 25        59         3         1
# 26        71        96        43
# 27        48        19        77
# 28         4        53       131
# 29         0       103         4
# 30         7        31         7
# 31        13        16        48
# 32       113         0         0
# 33        39        10         0
# 34        50         4         0
# 35        37         1        25
# 36        35        14         2
# 37         0        37        14
# 38        29         0         0
# 39        25        14         5
# 40         0         0        19
# 41        42         0         0

table(sce$k_50_walktrap,sce$Sample)
#    1c_LS_SCP 2c_LS_SCP 3c_LS_SCP
# 1        470       102        71
# 2         55        79        63
# 3        151        78        26
# 4        742       587        58
# 5        294       134       234
# 6        590       454       451
# 7        272       272        90
# 8        371         0       392
# 9        114       159        55
# 10       126       170        20
# 11        40       527       489
# 12         1        67       321
# 13       161        28        38
# 14        75        49        34
# 15        58        37         5
# 16        30        61        36
# 17       196        27         5
# 18       174         0         0
# 19        87        26        31
# 20        74        53       105
# 21        70         0         0
# 22        29        13         2
# 23       115         0         0
# 24        37         1        25

table(sce$k_75_walktrap,sce$Sample)
#     1c_LS_SCP 2c_LS_SCP 3c_LS_SCP
# 1        218       269       125
# 2        323       197        25
# 3        659        79       476
# 4        870       602       729
# 5        276       338       116
# 6        472       103        71
# 7        287         0         0 #probably glutamatergic cluster
# 8        746       598        67
# 9         36       476       507
# 10       149        76        26
# 11        55        79        63
# 12         0        53       278
# 13       161        28        37
# 14        80        26        31

#Preliminary explorration of distribution of expression of marker genes found that one cluster 
#was dominated by Slc17a7 expression and only contained sample 1. This is not driven by batch effect
#but rather biological (due to anatomy of section?).

#save the object
save(sce,file = here("processed-data","sce_clustered.rda"))

#check doublet score per cluster. 
#Going to make these plots and compare k=50 louvain and k=75 walktrap. 
#violin
#louvain
doublet_violin <- plotColData(object = sce,
                              x = "k_50_louvain_1",
                              y = "doubletScore",
                              colour_by = "k_50_louvain_1") +
    labs(x = "Cluster",
         y = "Doublet Score",
         title = "Doublet Score by Cluster") +
    theme(plot.title = element_text(hjust=0.5),legend.position = "none") +
    geom_hline(yintercept = 5)
ggsave(doublet_violin,filename = here("plots","doublet_score_by_cluster_k_50_louvain_violin.pdf"))

#walktrap
doublet_violin <- plotColData(object = sce,
                              x = "k_75_walktrap",
                              y = "doubletScore",
                              colour_by = "k_75_walktrap") +
    labs(x = "Cluster",
         y = "Doublet Score",
         title = "Doublet Score by Cluster") +
    theme(plot.title = element_text(hjust=0.5),legend.position = "none") +
    geom_hline(yintercept = 5)
ggsave(doublet_violin,filename = here("plots","doublet_score_by_cluster_k_75_walktrap_violin.pdf"))

#number of genes per cluster
#louvain
genes_violin <- plotColData(object = sce,
                              x = "k_50_louvain_1",
                              y = "detected",
                              colour_by = "k_50_louvain_1") +
    labs(x = "Cluster",
         y = "Number of genes/cell",
         title = "Number of Genes/Cell by Cluster") +
    theme(plot.title = element_text(hjust=0.5),legend.position = "none")
ggsave(genes_violin,filename = here("plots","Genes_by_cluster_k_50_louvain_violin.pdf"))

#walktrap
genes_violin <- plotColData(object = sce,
                            x = "k_75_walktrap",
                            y = "detected",
                            colour_by = "k_75_walktrap") +
    labs(x = "Cluster",
         y = "Number of genes/cell",
         title = "Number of Genes/Cell by Cluster") +
    theme(plot.title = element_text(hjust=0.5),legend.position = "none")
ggsave(genes_violin,filename = here("plots","Genes_by_cluster_k_75_walktrap_violin.pdf"))

#library size per cluster
#Louvain
lib_violin <- plotColData(object = sce,
                          x = "k_50_louvain_1",
                          y = "sum",
                          colour_by = "k_50_louvain_1") +
    scale_y_log10() +
    labs(x = "Cluster",
         title = "Total UMIs") +
    theme(plot.title = element_text(hjust=0.5),legend.position = "none")
ggsave(lib_violin,filename = here("plots","lib_size_by_cluster_k50_louvain_violin.pdf"))

#Walktrap
lib_violin <- plotColData(object = sce,
                          x = "k_75_walktrap",
                          y = "sum",
                          colour_by = "k_75_walktrap") +
    scale_y_log10() +
    labs(x = "Cluster",
         title = "Total UMIs") +
    theme(plot.title = element_text(hjust=0.5),legend.position = "none")
ggsave(lib_violin,filename = here("plots","lib_size_by_cluster_k75_walktrap_violin.pdf"))

#Check expression profiles of clusters with both algorithms. 
#DEfine some genes that are good markers. 
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

#louvain
Expression_dotplot <- plotDots(object = sce,
                               features = rev(genes),
                               group = "k_50_louvain_1",swap_rownames = "gene_name") +
    scale_color_gradientn(colours = c("lightgrey","orange","red"))
ggsave(plot = Expression_dotplot,filename = here("plots","Expression_plots",
                                                 "post_k_50_louvain_1_clustering","general_dotplot.pdf"),
       height = 8)

#Walktrap
Expression_dotplot <- plotDots(object = sce,
                               features = rev(genes),
                               group = "k_75_walktrap",swap_rownames = "gene_name") +
    scale_color_gradientn(colours = c("lightgrey","orange","red"))
ggsave(plot = Expression_dotplot,filename = here("plots","Expression_plots",
                                                 "post_k_75_walktrap_clustering","general_dotplot.pdf"),
       height = 8)

#The louvain algorithm identifies a cluster (cluster 11) that exhibits low # of genes/cell, but expresses
#only neuronal genes. I am pretty sure these are low quality nuclei. The walktrap algorithm just lumps 
#these cells into another cluster. We are going to remove them. 
#First what is the sample make up of this cluster. 
table(sce$Sample,sce$k_50_louvain_1 == 11)
#           FALSE TRUE
# 1c_LS_SCP  4063  269
# 2c_LS_SCP  2706  218
# 3c_LS_SCP  2456   95
#Primarily coming from samples 1 and 2, but all samples have cells within that cluster. 
#Will need to remove the cluster and rerun dimensionality reduction steps. 

