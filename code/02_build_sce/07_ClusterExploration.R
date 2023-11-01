#Goal: Explore clusters and identified marker genes + plot expression of known LS territories.
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

library(SingleCellExperiment)
library(ggplot2)
library(scater)
library(dplyr)
library(scran)
library(here)

#load the SingleCellExperiment object
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
# colData names(60): Sample Barcode ... k_75_walktrap CellType
# reducedDimNames(14): GLMPCA_approx UMAP_15 ... UMAP_mnn_25 UMAP_mnn_50
# mainExpName: NULL
# altExpNames(0):

#bargraph of number of cells in each population. 
celltype_bar <- as.data.frame(table(sce$CellType)) %>% 
    rename(CellType = Var1, Number_Nuclei = Freq) %>%
    ggplot(aes(x = reorder(CellType,-Number_Nuclei), y = Number_Nuclei,fill = CellType)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(x = "Cell Type",
         y = "Number of Nuclei")

ggsave(celltype_bar,filename = here("plots","No_Nuclei_per_CellType_bar.png"))

#load the DEG lists. 
load(here("processed-data","markers_1vAll_ttest.rda")) #markers_1vALL_df
load(here("processed-data","markers_pairwise_list.rda")) #markers_pairwise

#Gene list from figure 2 in Besnard + Leroy, Mol. Psy. 2022
#Fig. 2a "Gene expression domains defining non-overlapping LS territories.
domain_genes <- c("CBLN2","CBLN4","PAX6","SEMA3A", #Rostrocaudal Deep
                  "COL15A1","MATN2","WFS1", #Rostrocaudal Superficial
                  "DRD3","GALR1","IGFBP5","POU6F2", #Rostral Deep
                  "FOXP2","NDST4","GDPD2","CDH7","CRHR2","DACH2","FST","LHX2", #Rostral Superficial
                  "ASB4","OTX2","PDE1C","TAC1R","TRPC6", #Caudal Deep
                  "CD24A","IGFBP4") #Caudal Superficial 

#Some gene names might be differnt in the human genome.
#Check if all genes are in the object
table(domain_genes %in% rowData(sce)$gene_name)
# FALSE  TRUE 
#.    2    24

domain_genes[!(domain_genes %in% rowData(sce)$gene_name)]
#[1] "TAC1R" "CD24A"
#TAC1R name is same in human genome. 
#CD24a is CD24 in human genome. 

#Dot plot for domain markers
domain_dot <- plotDots(object = sce,features = rev(c("CBLN2","CBLN4","PAX6","SEMA3A", #Rostrocaudal Deep
                                                     "COL15A1","MATN2","WFS1", #Rostrocaudal Superficial
                                                     "DRD3","GALR1","IGFBP5","POU6F2", #Rostral Deep
                                                     "FOXP2","NDST4","GDPD2","CDH7","CRHR2","DACH2","FST","LHX2", #Rostral Superficial
                                                     "ASB4","OTX2","PDE1C","TRPC6", #Caudal Deep
                                                     "CD24","IGFBP4")), #Caudal Superficial 
                       swap_rownames = "gene_name",group = "CellType") +
    scale_color_gradientn(colours = c("lightgrey","orange","red")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = here("plots","Expression_plots","LS_Domains_dotplot.png"),plot = domain_dot)


#Make volcano plots for these genes
domain_violin <- plotExpression(object = sce,features = c("CBLN2","CBLN4","PAX6","SEMA3A", #Rostrocaudal Deep
                                                          "COL15A1","MATN2","WFS1", #Rostrocaudal Superficial
                                                          "DRD3","GALR1","IGFBP5","POU6F2", #Rostral Deep
                                                          "FOXP2","NDST4","GDPD2","CDH7","CRHR2","DACH2","FST","LHX2", #Rostral Superficial
                                                          "ASB4","OTX2","PDE1C","TRPC6", #Caudal Deep
                                                          "CD24","IGFBP4"),
                                swap_rownames = "gene_name",x = "CellType",
                                colour_by = "CellType",ncol = 5) +
    theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")
ggsave(filename = here("plots","Expression_plots","LS_Domains_violin.png"),
       plot = domain_violin,width = 14,height = 16)

#UMAPs
domain_markers <- c("CBLN2","CBLN4","PAX6","SEMA3A", #Rostrocaudal Deep
                    "COL15A1","MATN2","WFS1", #Rostrocaudal Superficial
                    "DRD3","GALR1","IGFBP5","POU6F2", #Rostral Deep
                    "FOXP2","NDST4","GDPD2","CDH7","CRHR2","DACH2","FST","LHX2", #Rostral Superficial
                    "ASB4","OTX2","PDE1C","TRPC6", #Caudal Deep
                    "CD24","IGFBP4")

for(i in domain_markers){
    umap_domain <- plotReducedDim(object = sce,
                                  dimred = "UMAP_mnn_15",
                                  colour_by = i,
                                  swap_rownames = "gene_name") +
        scale_color_gradientn(colours = c("lightgrey","red")) +
        ggtitle(i) +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(filename = paste0("plots/Expression_plots/domain_markers_UMAP_mnn_15/",i,"_expression_umap_15_dims.pdf"),
           plot = umap_domain,
           height = 8,width = 8)
}

#Additional investigation of marker genes. 
LS_GABA_1_markers <- subset(markers_1vALL_df,subset=(cellType.target == "LS_GABA_1"))[1:50,]
LS_GABA_2_markers <- subset(markers_1vALL_df,subset=(cellType.target == "LS_GABA_2"))[1:50,]
for(i in 1:nrow(LS_GABA_2_markers)){
    print(i)
    x <- plotReducedDim(object = sce,
                        dimred = "UMAP_mnn_15",
                        colour_by = LS_GABA_2_markers[i,"gene_name"],
                        swap_rownames = "gene_name") +
        scale_color_gradientn(colours = c("lightgrey","orange","red"))
    ggsave(plot = x,
           filename = here("plots","Expression_plots","LS_GABA_2_Markers",
                           paste0(i,"_",LS_GABA_2_markers[i,"gene_name"],"_umap.pdf")))
}
p1 <- plotReducedDim(object = sce,
               dimred = "UMAP_mnn_15",
               colour_by = "DGKG",
               swap_rownames = "gene_name") +
    scale_color_gradientn(colours = c("lightgrey","orange","red")) +
    ggtitle("DGKG (LS Marker)") +
    theme(plot.title = element_text(hjust = 0.5))

p2 <- plotReducedDim(object = sce,
                     dimred = "UMAP_mnn_15",
                     colour_by = "TRPC4",
                     swap_rownames = "gene_name") +
    scale_color_gradientn(colours = c("lightgrey","orange","red")) +
    ggtitle("TRPC4 (LS Marker)")


p3 <- plotReducedDim(object = sce,
                     dimred = "UMAP_mnn_15",
                     colour_by = "RARB",
                     swap_rownames = "gene_name") +
    scale_color_gradientn(colours = c("lightgrey","orange","red")) +
    ggtitle("RARB (Striatal Marker)") +
    theme(plot.title = element_text(hjust = 0.5))

p4 <- plotReducedDim(object = sce,
                     dimred = "UMAP_mnn_15",
                     colour_by = "BCL11B",
                     swap_rownames = "gene_name") +
    scale_color_gradientn(colours = c("lightgrey","orange","red")) +
    ggtitle("BCL11B (Striatal Marker)") +
    theme(plot.title = element_text(hjust = 0.5))

x <- cowplot::plot_grid(plotlist = list(p1,p2,p3,p4),ncol = 2)
ggsave(x, filename = here("plots","Expression_plots","LS_Str_Markers.pdf"))



#####Cluster modularity
#Build the graph again. 
#Used k=50 + walktrap clustering for celltype designations. 
snn_k_50 <- buildSNNGraph(sce, k = 50, use.dimred = "mnn",type="jaccard")

# "Rather, we use the pairwiseModularity() function from bluster with as.ratio=TRUE, 
# which returns the ratio of the observed to expected sum of weights between each pair of clusters. 
# We use the ratio instead of the difference as the former is less dependent on the number of cells 
# in each cluster" - OSCA advanced, Chatper 5.2.5
k_50_modularity <- bluster::pairwiseModularity(graph = snn_k_50,
                                               clusters = sce$CellType,
                                               as.ratio = TRUE)


library(pheatmap)
pdf(file = here("plots","k_50_pairwise_modularity_celltypes.pdf"))
pheatmap(log2(k_50_modularity+1), cluster_rows=FALSE, cluster_cols=FALSE,
         color=colorRampPalette(c("white","orange","red"))(100))
dev.off()
# "Indeed, concentration of the weight on the diagonal of (Figure 5.5) indicates 
# that most of the clusters are well-separated, while some modest off-diagonal entries 
# represent closely related clusters with more inter-connecting edges." - OSCA advanced, chatper 5.2.5

cluster.gr <- igraph::graph_from_adjacency_matrix(log2(k_50_modularity+1),
                                                  mode="upper", weighted=TRUE, diag=FALSE)

# Increasing the weight to increase the visibility of the lines.
set.seed(11001010)
pdf(file = here("plots","k_50_relationship_between_clusters_graph.pdf"))
plot(cluster.gr, edge.width=igraph::E(cluster.gr)$weight*10,
     layout=igraph::layout_with_lgl)
dev.off()

# Force-based layout showing the relationships between clusters based on the log-ratio of observed to expected total weights 
# between nodes in different clusters. The thickness of the edge between a pair of clusters is proportional to the corresponding 
# log-ratio. - OSCA advanced, chapter 5.2.5, figure 5.6
table(sce$CellType,sce$Sample)





