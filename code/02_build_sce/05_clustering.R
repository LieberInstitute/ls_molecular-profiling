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
load(here("processed-data","02_build_sce","sce_postMNN_081724.rda"),verbose = TRUE)

sce

identical(rownames(colData(sce)),colnames(sce))

#Run clustering with k=10, k=15 and k=20 at 2 different resolution values (0.75 and 1.0)
##k=10##
set.seed(10)
snn_k_10 <- buildSNNGraph(sce, k = 10, use.dimred = "mnn",type="jaccard")

#Louvain clustering
set.seed(10)
clust_10_pt75 <- igraph::cluster_louvain(snn_k_10,resolution=0.75)$membership
table(clust_10_pt75)

#res = 1
set.seed(10)
clust_10_1 <- igraph::cluster_louvain(snn_k_10,resolution=1)$membership
table(clust_10_1)

#Add the cluster information to the colData
sce$k_10_louvain_pt75 <- factor(clust_10_pt75)
table(sce$k_10_louvain_pt75,sce$Sample)

sce$k_10_louvain_1 <- factor(clust_10_1)
table(sce$k_10_louvain_1,sce$Sample)

#Plot tSNE with k=10, res=0.75
x <- plotReducedDim(object = sce,
                    dimred = "tSNE_mnn_50",
                    colour_by = "k_10_louvain_pt75",
                    text_by = "k_10_louvain_pt75") +
  ggtitle("k=10 louvain res=0.75 clustering") +
  theme(plot.title = element_text(hjust = 0.5))
#Save as png and pdf
ggsave(plot = x,filename = here("plots","Dim_Red","k_10_louvain_pt75_tSNE_mnn_50.pdf"),height = 8, width = 8)
ggsave(plot = x,filename = here("plots","Dim_Red","k_10_louvain_pt75_tSNE_mnn_50.png"),height = 8, width = 8)

#Plot tSNE with k=10, res=1
x <- plotReducedDim(object = sce,
                    dimred = "tSNE_mnn_50",
                    colour_by = "k_10_louvain_1",
                    text_by = "k_10_louvain_1") +
  ggtitle("k=10 louvain 1 clustering") +
  theme(plot.title = element_text(hjust = 0.5))
#Save as png and pdf
ggsave(plot = x,filename = here("plots","Dim_Red","k_10_louvain_1_tSNE_mnn_50.pdf"),height = 8, width = 8)
ggsave(plot = x,filename = here("plots","Dim_Red","k_10_louvain_1_tSNE_mnn_50.png"),height = 8, width = 8)


##k=15##
set.seed(15)
snn_k_15 <- buildSNNGraph(sce, k = 15, use.dimred = "mnn",type="jaccard")

#Louvain clustering
set.seed(15)
clust_15_pt75 <- igraph::cluster_louvain(snn_k_15,resolution=0.75)$membership
table(clust_15_pt75)

#res = 1
set.seed(15)
clust_15_1 <- igraph::cluster_louvain(snn_k_15,resolution=1)$membership
table(clust_15_1)

#Add the cluster information to the colData
sce$k_15_louvain_pt75 <- factor(clust_15_pt75)
table(sce$k_15_louvain_pt75,sce$Sample)

sce$k_15_louvain_1 <- factor(clust_15_1)
table(sce$k_15_louvain_1,sce$Sample)

#Plot tSNE with k=15, res=0.75
x <- plotReducedDim(object = sce,
                    dimred = "tSNE_mnn_50",
                    colour_by = "k_15_louvain_pt75",
                    text_by = "k_15_louvain_pt75") +
  ggtitle("k=15 louvain res=0.75 clustering") +
  theme(plot.title = element_text(hjust = 0.5))
#Save as png and pdf
ggsave(plot = x,filename = here("plots","Dim_Red","k_15_louvain_pt75_tSNE_mnn_50.pdf"),height = 8, width = 8)
ggsave(plot = x,filename = here("plots","Dim_Red","k_15_louvain_pt75_tSNE_mnn_50.png"),height = 8, width = 8)

#Plot tSNE with k=15, res=1
x <- plotReducedDim(object = sce,
                    dimred = "tSNE_mnn_50",
                    colour_by = "k_15_louvain_1",
                    text_by = "k_15_louvain_1") +
  ggtitle("k=15 louvain 1 clustering") +
  theme(plot.title = element_text(hjust = 0.5))
#Save as png and pdf
ggsave(plot = x,filename = here("plots","Dim_Red","k_15_louvain_1_tSNE_mnn_50.pdf"),height = 8, width = 8)
ggsave(plot = x,filename = here("plots","Dim_Red","k_15_louvain_1_tSNE_mnn_50.png"),height = 8, width = 8)

##k=20##
set.seed(20)
snn_k_20 <- buildSNNGraph(sce, k = 20, use.dimred = "mnn",type="jaccard")

#Louvain clustering
set.seed(20)
clust_20_pt75 <- igraph::cluster_louvain(snn_k_20,resolution=0.75)$membership
table(clust_20_pt75)

#res = 1
set.seed(20)
clust_20_1 <- igraph::cluster_louvain(snn_k_20,resolution=1)$membership
table(clust_20_1)

#Add the cluster information to the colData
sce$k_20_louvain_pt75 <- factor(clust_20_pt75)
table(sce$k_20_louvain_pt75,sce$Sample)

sce$k_20_louvain_1 <- factor(clust_20_1)
table(sce$k_20_louvain_1,sce$Sample)

#Plot tSNE with k=20, res=0.75
x <- plotReducedDim(object = sce,
                    dimred = "tSNE_mnn_50",
                    colour_by = "k_20_louvain_pt75",
                    text_by = "k_20_louvain_pt75") +
  ggtitle("k=20 louvain res=0.75 clustering") +
  theme(plot.title = element_text(hjust = 0.5))
#Save as png and pdf
ggsave(plot = x,filename = here("plots","Dim_Red","k_20_louvain_pt75_tSNE_mnn_50.pdf"),height = 8, width = 8)
ggsave(plot = x,filename = here("plots","Dim_Red","k_20_louvain_pt75_tSNE_mnn_50.png"),height = 8, width = 8)

#Plot tSNE with k=20, res=1
x <- plotReducedDim(object = sce,
                    dimred = "tSNE_mnn_50",
                    colour_by = "k_20_louvain_1",
                    text_by = "k_20_louvain_1") +
  ggtitle("k=20 louvain 1 clustering") +
  theme(plot.title = element_text(hjust = 0.5))
#Save as png and pdf
ggsave(plot = x,filename = here("plots","Dim_Red","k_20_louvain_1_tSNE_mnn_50.pdf"),height = 8, width = 8)
ggsave(plot = x,filename = here("plots","Dim_Red","k_20_louvain_1_tSNE_mnn_50.png"),height = 8, width = 8)

#Preliminary explorration of distribution of expression of marker genes found that one cluster 
#was dominated by Slc17a7 expression and only contained sample 1. This is not driven by batch effect
#but rather biological (due to anatomy of section).

#Now that we have clusters, calculate sum factors and compute logcounts
logcounts(sce) <- NULL
sce <- computeSumFactors(sce,cluster = sce$k_20_louvain_1,min.mean = 0.1)
sce <- logNormCounts(sce)

#save the object
save(sce,file = here("processed-data","02_build_sce","sce_clustered_082024.rda"))

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
#No cluster dominated by high doublet score

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

#Define some genes that are good markers. 
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
           "PDGFRA", "VCAN", "CSPG4", #Polydendrocytes
           "DRD1","OPRM1","RXFP1","EBF1","PDYN", "CHST9","SEMA5B","TAC1","STXBP6",#Additional D1 markers including D1 islands.
           "DRD2","PENK","ADORA2A", #D2 markers.
           "CRYM",#Medial dorsal Striatum marker
           "DLK1",#ventral medial dorsal striaum
           "FOXP2","PDE1B","KIAA1211L","PDE2A","SLIT3","NGEF")

#Primary goal of making this plot is to make sure that there are no low quality clusters.  
Expression_dotplot <- plotDots(object = sce,
                               features = rev(genes),
                               group = "k_20_louvain_1",swap_rownames = "gene_name") +
  scale_color_gradientn(colours = c("lightgrey","orange","red"))
ggsave(plot = Expression_dotplot,filename = here("plots","Expression_plots",
                                                 "post_k_20_louvain_1_clustering_general_dotplot.pdf"),
       height = 8)

#Plot expression on tSNE
for(i in genes){
  x <- plotReducedDim(object = sce,
                      dimred = "tSNE_mnn_50",
                      colour_by = i,
                      swap_rownames = "gene_name") +
    scale_color_gradientn(colours = c("lightgrey","orange","red"))
  ggsave(plot = x,
         filename = here("plots","Expression_plots","post_k_20_louvain_clustering","tSNE",paste0(i,"_tSNE_mnn_50.png")))
  }

#Violin plots
for(i in genes){
  x <- plotExpression(sce,features = i,x = "k_20_louvain_1",colour_by = "k_20_louvain_1",swap_rowname = "gene_name") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3) 
  ggsave(plot = x,
         filename = here("plots","Expression_plots","post_k_20_louvain_clustering","Violin",paste0(i,"_k_20_louvain_1_violin.png")),
         height = 8, width = 12)
         }

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
