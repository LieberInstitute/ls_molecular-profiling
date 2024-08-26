#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

library(SingleCellExperiment)
library(DeconvoBuddies)
library(sessioninfo)
library(Polychrome)
library(ggplot2)
library(scater)
library(scran)
library(here)

#load the clustered object
load(here("processed-data","02_build_sce","sce_clustered_082024.rda"),verbose = TRUE)

sce

identical(rownames(colData(sce)),colnames(sce))

stopifnot(identical(rownames(colData(sce)),colnames(sce)))

##################################################
############   Annotate clusters      ############    
##################################################
#Read in annotation dataframe
annotation_df <- read.csv(file = here("processed-data","k_20_louvain_1_Annotation.csv"))

#add celltype info
sce$CellType.Final <- annotation_df$CellType[match(sce$k_20_louvain_1,
                                                   annotation_df$k_20_louvain_1)]

#factorize. 
sce$CellType.Final <- factor(sce$CellType.Final,
                             levels = c("LS_Inh_A","LS_Inh_B","LS_Inh_C","LS_Inh_D","LS_Inh_E",
                                        "MS_Inh_A","MS_Inh_B","MS_Inh_C","MS_Inh_D",
                                        "Sept_Inh_A","Sept_Inh_B","Sept_Inh_C",
                                        "Str_DRD1_MSN_A","Str_DRD1_MSN_B","Str_DRD1_Patch","Str_DRD2_MSN",
                                        "MS_Excit_A","Excitatory_A","Excitatory_B",
                                        "Oligo","OPC","Astrocyte","Ependymal","Microglia","Mural"))

#Annotate the tSNE
cluster_cols <- Polychrome::createPalette(length(unique(sce$CellType.Final)),
                                          c("#D81B60", "#1E88E5","#FFC107","#009E73"))
names(cluster_cols) <- unique(sce$CellType.Final)
save(cluster_cols,file = here("processed-data","cluster_cols_CellType_Final_082524.rda"))

annotated_tSNE <- plotReducedDim(object = sce,
                                 dimred = "tSNE_mnn_50",
                                 colour_by = "CellType.Final",
                                 text_by = "CellType.Final") +
  scale_color_manual(values = cluster_cols) +
  theme(legend.position = "none")
ggsave(filename = here("plots","Dim_Red","tSNE_mnn_50_annotated_CellType_Final.pdf"),
       plot = annotated_tSNE)


#Save the object post cell type annotation
save(sce,file = here("processed-data","02_build_sce","sce_celltype_082624.rda"))

########Calculate modularity scores.
set.seed(20)
#Make the graph. 
snn_k_20 <- buildSNNGraph(sce, k = 20, use.dimred = "mnn",type="jaccard")

#Calcualte modularity scores. 
k_20_modularity <- bluster::pairwiseModularity(graph = snn_k_20,
                                               clusters = sce$CellType.Final,
                                               as.ratio = TRUE)

#make a heatmap. 
library(pheatmap)
pdf(file = here("plots","CellType_k_20_louvain_pairwise_modularity.pdf"))
pheatmap(log2(k_20_modularity+1), 
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         display_numbers=TRUE, 
         number_format="%.2f", 
         fontsize_number=6.5,
         main = "Modularity ratio for 25 clusters in human LS (n=3)",
         color=colorRampPalette(c("white","orange","red"))(100))
dev.off()

##################################################
#########Prep sce object for DEG testing##########
##################################################

#Do any genes have 0 counts for every cell. 
table(rowSums(assay(sce, "counts")) == 0)

#Remove genes that are all 0s 
sce <- sce[!rowSums(assay(sce, "counts")) == 0, ]

##################################################
###############run 1 vs all testing###############
##################################################
markers_1vALL_enrich <- findMarkers_1vAll(sce, 
                                          assay_name = "logcounts", 
                                          cellType_col = "CellType.Final", 
                                          direction = "up",
                                          mod = "~Sample")
#From DeconvoBuddies help page: If "up" genes with logFC < 0 will have p.value = 1.

#Add symbol information to the table
#First change the ensembl gene id column to have same name as what is in rowData(sce)
colnames(markers_1vALL_enrich)[1] <- "gene_id"
markers_1vALL_df <- dplyr::left_join(x = as.data.frame(markers_1vALL_enrich),
                                     y = as.data.frame(rowData(sce)[,c("gene_id","gene_name")]),
                                     by = "gene_id")

#save the dataframe. 
save(markers_1vALL_df,file = here("processed-data","markers_1vAll_ttest_k_20_louvain_CellType_Final.rda"))
###############################


##################################################
###############run pairwise testing###############
##################################################
mod <- with(colData(sce), model.matrix(~ Sample))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

# Run pairwise t-tests
markers_pairwise <- findMarkers(sce, 
                                groups=sce$CellType.Final,
                                assay.type="logcounts", 
                                design=mod, 
                                test="t",
                                direction="up", 
                                pval.type="all", 
                                full.stats=T)

#How many DEGs for each cluster? 
sapply(markers_pairwise, function(x){table(x$FDR<0.05)})

#Add gene info to each list.
for(i in names(markers_pairwise)){
  markers_pairwise[[i]] <- as.data.frame(markers_pairwise[[i]])
  markers_pairwise[[i]]$gene_id <- row.names(markers_pairwise[[i]])
  markers_pairwise[[i]] <- dplyr::left_join(x  =  markers_pairwise[[i]],
                                            y  =  as.data.frame(rowData(sce)[,c("gene_id","gene_name")]),
                                            by = "gene_id")
}

save(markers_pairwise,file = here("processed-data","markers_pairwise_list_k_20_louvain_CellTypeFinal.rda"))


print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
