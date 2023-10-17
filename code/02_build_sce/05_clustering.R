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
# reducedDimNames(5): GLMPCA_approx UMAP TSNE mnn UMAP_mnn
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
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 276 1125  307  200  258 1095  475  542  674  645  225  147  692  143  106  517 
#  17   18   19   20   21   22   23   24 
# 252   72   49  168  637  126  310  959

set.seed(1234)
clust_50 <- igraph::cluster_louvain(snn_k_50,resolution=1)$membership
table(clust_50)
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 1430 1125  362  197  257 1446  545  657  643  226  981  350  694  154  637  296 

set.seed(1234)
clust_75 <- igraph::cluster_louvain(snn_k_75,resolution=1)$membership
table(clust_75)
clust_75
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 1365 1224  355  193  256 1456  545  650  645 1382  750  214  143  636  18

#Add cluster information to object
sce$k_20_louvain_1 <- factor(clust_20)
sce$k_50_louvain_1 <- factor(clust_50)
sce$k_75_louvain_1 <- factor(clust_75)

#Also run walktrap clustering. 
set.seed(1234)
wt_clusters_k_20 <- igraph::cluster_walktrap(snn_k_20)$membership
table(wt_clusters_k_20)
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 258 1087  301  953  436  379  218 1214  807  851  253  157   61  379  286  307 
# 17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32 
# 265  164  113  176   87  170  143  144   58  107   95   77  119   48   54   56 
# 33   34   35   36 
# 63   29   43   42 

set.seed(1234)
wt_clusters_k_50 <- igraph::cluster_walktrap(snn_k_50)$membership
table(wt_clusters_k_50)
# wt_clusters_k_50
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 289  163 1879 1469 1228  643  257  197 1381  283  316  215  332  301  154  149 
#  17   18   19   20 
# 228  258  188   70 

#Add cluster information to object
sce$k_20_walktrap <- factor(wt_clusters_k_20)
sce$k_50_walktrap <- factor(wt_clusters_k_50)

#Plot the clusters on the umap
#k=20 louvain
k_20_umap <- plotReducedDim(object = sce,dimred = "UMAP_mnn",
                            colour_by = "k_20_louvain_1",text_by = "k_20_louvain_1") #Color and label by cluster 
ggsave(plot = k_20_umap,filename = here("plots","Dim_Red","k_20_louvain_umap_50components.pdf"))

#k=50 louvain
k_50_umap <- plotReducedDim(object = sce,dimred = "UMAP_mnn",
                            colour_by = "k_50_louvain_1",text_by = "k_50_louvain_1") #Color and label by cluster 
ggsave(plot = k_50_umap,filename = here("plots","Dim_Red","k_50_louvain_umap_50components.pdf"))

#k=75 louvain
k_75_umap <- plotReducedDim(object = sce,dimred = "UMAP_mnn",
                            colour_by = "k_75_louvain_1",text_by = "k_75_louvain_1") #Color and label by cluster 
ggsave(plot = k_50_umap,filename = here("plots","Dim_Red","k_75_louvain_umap_50components.pdf"))

#k=20 walktrap
k_20_wt_umap <- plotReducedDim(object = sce,dimred = "UMAP_mnn",
                               colour_by = "k_20_walktrap",text_by = "k_20_walktrap") #Color and label by cluster 
ggsave(plot = k_20_wt_umap,filename = here("plots","Dim_Red","k_20_walktrap_umap_50components.pdf"))

#k=50 walktrap
k_50_wt_umap <- plotReducedDim(object = sce,dimred = "UMAP_mnn",
                               colour_by = "k_50_walktrap",text_by = "k_50_walktrap") #Color and label by cluster 
ggsave(plot = k_50_wt_umap,filename = here("plots","Dim_Red","k_50_walktrap_umap_50components.pdf"))

#Tables by cluster and sample
table(sce$k_20_louvain_1,sce$Sample)
#    1c_LS_SCP 2c_LS_SCP 3c_LS_SCP
# 1        167        96        13
# 2        649        33       443
# 3        307         0         0
# 4         57        80        63
# 5        152        75        31
# 6        550       493        52
# 7         38       284       153
# 8        321       196        25
# 9        295       281        98
# 10       472       102        71
# 11       158        29        38
# 12        48        21        78
# 13       211       214       267
# 14        59        31        53
# 15        62        39         5
# 16       257       131       129
# 17        83        57       112
# 18        72         0         0
# 19        39        10         0
# 20       109        27        32
# 21       234       263       140
# 22        96         4        26
# 23        20       103       187
# 24         5       301       653

table(sce$k_50_louvain_1,sce$Sample)
#      1c_LS_SCP 2c_LS_SCP 3c_LS_SCP
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
#     1c_LS_SCP 2c_LS_SCP 3c_LS_SCP
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
# 3         19        99       183
# 4        626        34       293
# 5        321        94        21
# 6        134       173        72
# 7        154        27        37
# 8        464       347       403
# 9          5       288       514
# 10       150       474       227
# 11        83        57       113
# 12        75        47        35
# 13        36        25         0
# 14       242        88        49
# 15       174        99        13
# 16       183        95        29
# 17       229        14        22
# 18       106        26        32
# 19         5        53        55
# 20        24         0       152
# 21        52        27         8
# 22       170         0         0
# 23         0         9       134
# 24        48        19        77
# 25        55         3         0
# 26         0       103         4
# 27        20        43        32
# 28        13        16        48
# 29       119         0         0
# 30        38        10         0
# 31        50         4         0
# 32        38        14         4
# 33        37         1        25
# 34        29         0         0
# 35        24        14         5
# 36        42         0         0

table(sce$k_50_walktrap,sce$Sample)
#     1c_LS_SCP 2c_LS_SCP 3c_LS_SCP
# 1        289         0         0
# 2         95        38        30
# 3        724       468       687
# 4        158       715       596
# 5        661        85       482
# 6        470       102        71
# 7        151        75        31
# 8         55        79        63
# 9        738       585        58
# 10        86       143        54
# 11       126       170        20
# 12       155        27        33
# 13         0        54       278
# 14       180        92        29
# 15        73        47        34
# 16        92        26        31
# 17       196        27         5
# 18        83        60       115
# 19        59        77        52
# 20        70         0         0
#Preliminary explorration of distribution of expression of marker genes found that one cluster 
#was dominated by Slc17a7 expression and only contained sample 1. This is not driven by batch effect
#but rather biological (due to anatomy of section?). k=20 seems to contain several clusters dominated by one 
#sample, while k=50 only has one cluster that is dominated by sample 1 (and is the cluster identified yesterday)
#going to check expression profiles of k=50 louvain. 

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
           "COL1A2", "TBX18", "RBPMS") #Mural

Expression_dotplot <- plotDots(object = sce,
                               features = rev(genes),
                               group = "k_50_louvain_1",swap_rownames = "gene_name") +
    scale_color_gradientn(colours = c("lightgrey","red"))
ggsave(plot = Expression_dotplot,filename = here("plots","Expression_plots",
                                                 "post_k_50_louvain_1_clustering","general_dotplot.pdf"))

#Annotate clusters and remake the Expression dotplot. 
annotation_df <- data.frame(cluster = c(1:16))
annotation_df$CellType <- c("LS_GABA_1","LS_GABA_Drd3","Glutamatergic","GABA_undefined_1",
                            "Microglia","Oligodendrocyte_1","Astrocyte_1","Oligodendrocyte_2",
                            "LS_GABA_2","MS_GABA_1","MS_GABA_2","MS_GABA_3",
                            "GABA_undefined_2","Mural","Astrocyte_2","LS_GABA_3")

sce$CellType <- annotation_df$CellType[match(sce$k_50_louvain_1,annotation_df$cluster)]

sce$CellType <- factor(x = sce$CellType,
                       levels = c("LS_GABA_1","LS_GABA_2","LS_GABA_3","LS_GABA_Drd3",
                                  "MS_GABA_1","MS_GABA_2","MS_GABA_3","GABA_undefined_1",
                                  "GABA_undefined_2","Glutamatergic","Astrocyte_1","Astrocyte_2",
                                  "Oligodendrocyte_1","Oligodendrocyte_2","Microglia","Mural"))

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
                                                "COL1A2", "TBX18", "RBPMS")),
                               group = "CellType",swap_rownames = "gene_name") +
    theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
    scale_color_gradientn(colours = c("lightgrey","orange","red")) 
ggsave(plot = Expression_dotplot,filename = here("plots","Expression_plots",
                                                 "post_k_50_louvain_1_clustering","general_dotplot_annotated.pdf"))


#Make umap annotated
#k=50 louvain
annotated_umap <- plotReducedDim(object = sce,dimred = "UMAP_mnn",
                                 colour_by = "CellType",text_by = "CellType")
ggsave(plot = annotated_umap,filename = here("plots","Dim_Red","Annotated_k_50_louvain_umap_50components.pdf"))

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


