#Goal: Compare gene expression signatures of the mouse and human LS
#Code modified from https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/51d15ef9f5f2c4c53f55e22e3fe467de1a724668/10x_NAc-n8_step04_cross-species_rnNAc_MNT.R#L4
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/
#module load conda_R/4.4

library(SingleCellExperiment)
library(sparseMatrixStats)
library(DeconvoBuddies)
library(RColorBrewer)
library(org.Hs.eg.db) #human
library(org.Mm.eg.db) #mouse
library(sessioninfo)
library(pheatmap)
library(rafalib)
library(here)

#####load and prep data for the analysis######
######HUMAN prep
#load the SingleCellExperiment object for human Lateral Septum
load(here("processed-data","02_build_sce","sce_celltype.rda"),verbose = TRUE)

sce
# class: SingleCellExperiment 
# dim: 36601 9225 
# metadata(1): Samples
# assays(3): counts binomial_pearson_residuals logcounts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(9225): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 3_TTTGGTTTCTTCGACC-1 3_TTTGTTGTCCCGATCT-1
# colData names(60): Sample Barcode ... k_20_louvain_1 CellType.Final
# reducedDimNames(4): GLMPCA_approx tSNE_50 mnn tSNE_mnn_50
# mainExpName: NULL
# altExpNames(0):

#Make sure everything is in correct order. 
identical(rownames(colData(sce)),colnames(sce))
#[1] TRUE

#load human DEG list from human
load(here("processed-data","markers_1vAll_ttest_k_20_louvain_CellType_Final.rda"),verbose = TRUE) 
# Loading objects:
#   markers_1vALL_df

dim(markers_1vALL_df)
#[1] 838900      9

#Split the markers_1vALL_df into a list. 
#f = splits the list by the cell type. 
markers_1vALL_list_human <- split(markers_1vALL_df,
                                  f = markers_1vALL_df$cellType.target)

#Make the rownames of each element of the list the ensemble gene id. 
#Change the rownames of the list to be the gene_id
markers_1vALL_list_human <- lapply(markers_1vALL_list_human, function(x){
  rownames(x) <- x$gene_id
  return(x)
})

#Save the markers_1vALL_list 
save(markers_1vALL_list_human,file = here("processed-data","markers_1vAll_human_list.rda"))

###### MOUSE prep
#load the SingleCellExperiment object for mouse Lateral Septum
#This object generated in 03_cross_species/mouse_LS_1vALL_DEGs.R after removing problematic clusters. 
load(file = here("processed-data","02_build_sce","mouse_sce.rda"),verbose = TRUE)
# Loading objects:
#   sce.ls

sce.ls
# class: SingleCellExperiment 
# dim: 32285 21884 
# metadata(1): Samples
# assays(3): counts binomial_pearson_residuals logcounts
# rownames(32285): ENSMUSG00000051951 ENSMUSG00000089699 ...
# ENSMUSG00000095019 ENSMUSG00000095041
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(21884): 1_AAACCCAAGGTACATA-1 1_AAACCCACATCCGAGC-1 ...
# 4_TTTGTTGCATACAGCT-1 4_TTTGTTGGTCAAACGG-1
# colData names(17): Sample Barcode ... cellType.final cellType.broad
# reducedDimNames(4): GLMPCA_approx UMAP TSNE GLMPCA_50
# mainExpName: NULL
# altExpNames(0):

identical(rownames(colData(sce.ls)),colnames(sce.ls))
#[1] TRUE


#load the 1vALL mouse DEGs
load(here("processed-data","mouse_markers_1vAll.rda"),verbose = TRUE)
# Loading objects:
#   markers_1vALL_mouse

#Split the markers_1vALL_mouse into a list
#Split list by celltype with f=
markers_1vALL_mouse_list <- split(markers_1vALL_mouse,
                                  f = markers_1vALL_mouse$cellType.target)

#Change the rownames of the list to be the gene_id
markers_1vALL_mouse_list <- lapply(markers_1vALL_mouse_list, function(x){
  rownames(x) <- x$gene_id
  return(x)
})

#save the list. 
save(markers_1vALL_mouse_list,
     file = here("processed-data","mouse_markers_1vAll_list.rda"))

#rename the objects to make more sense. 
sce_human_ls <- sce 
sce_mouse_ls <- sce.ls
rm(sce,sce.ls)
