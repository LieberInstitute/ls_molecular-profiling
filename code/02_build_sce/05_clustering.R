#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

#Load libraries
library(SingleCellExperiment)
library(sessioninfo)
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
#Smaller= more, small clusters
snn_k_20 <- buildSNNGraph(sce, k = 20, use.dimred = "mnn",type="jaccard")
snn_k_50 <- buildSNNGraph(sce, k = 50, use.dimred = "mnn",type="jaccard")

#Run louvain clustering. 
#jaccard + louvain is similar to seurat workflow. 
#Resolution=1 is the default method
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

#Add cluster information to object
sce$k_20_louvain_1 <- factor(clust_20)
sce$k_50_louvain_1 <- factor(clust_50)

#Plot the clusters on the umap
#k=20
k_20_umap <- plotReducedDim(object = sce,dimred = "UMAP_mnn",
                            colour_by = "k_20_louvain_1",text_by = "k_20_louvain_1_50components") #Color and label by cluster 
ggsave(plot = k_20_umap,filename = here("plots","Dim_Red","k_20_louvain_umap.pdf"))

#k=50
k_50_umap <- plotReducedDim(object = sce,dimred = "UMAP_mnn",
                            colour_by = "k_50_louvain_1",text_by = "k_50_louvain_1_50components") #Color and label by cluster 
ggsave(plot = k_50_umap,filename = here("plots","Dim_Red","k_50_louvain_umap.pdf"))


