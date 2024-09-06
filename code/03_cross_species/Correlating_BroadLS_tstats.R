#Goal: Correlate standardized.logFC from all LS clusters. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(here)


##################
### Prep Human ###
##################

#load the SingleCellExperiment object for human Lateral Septum
load(here("processed-data","sce_human_sub.rda"),verbose = TRUE)
# Loading objects:
#   sce_human_sub

sce_human_sub
# class: SingleCellExperiment 
# dim: 16574 9225 
# metadata(1): Samples
# assays(3): counts binomial_pearson_residuals logcounts
# rownames(16574): ENSG00000187634 ENSG00000188976 ... ENSG00000205542
# ENSG00000092345
# rowData names(9): source type ... hs.entrezIds JAX.geneID
# colnames(9225): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 3_TTTGGTTTCTTCGACC-1 3_TTTGTTGTCCCGATCT-1
# colData names(60): Sample Barcode ... k_20_louvain_1 CellType.Final
# reducedDimNames(4): GLMPCA_approx tSNE_50 mnn tSNE_mnn_50
# mainExpName: NULL
# altExpNames(0):

identical(rownames(colData(sce_human_sub)),colnames(sce_human_sub))
# [1] TRUE

#Load the broad LS DEGs
load(here("processed-data","human_1vALL_broad.rda"),verbose = TRUE)
# Loading objects:
#   human_1vALL_broad

#subset the dataframe for LS
human_1vALL_broad <- subset(human_1vALL_broad, subset = (cellType.target == "LS"))

dim(human_1vALL_broad)
#[1] 36601     9

#subset DEGs for only homologs 
human_1vALL_homol <- merge(human_1vALL_broad,
                           as.data.frame(rowData(sce_human_sub)[,c("gene_id","JAX.geneID")]),
                           by = "gene_id") 

dim(human_1vALL_homol)
#[1] 16574    10

#Change rownames to JAX.geneID
rownames(human_1vALL_homol) <- human_1vALL_homol$JAX.geneID

##################
### Prep Mouse ###
##################

#Load mouse object 
#This object generated in 03_cross_species/mouse_LS_1vALL_DEGs.R after removing problematic clusters. 
load(here("processed-data","sce_mouse_sub.rda"),verbose = TRUE)
# Loading objects:
#   sce_mouse_sub

sce_mouse_sub
# dim: 16574 21884 
# metadata(1): Samples
# assays(3): counts binomial_pearson_residuals logcounts
# rownames(16574): ENSMUSG00000096351 ENSMUSG00000095567 ...
# ENSMUSG00000049775 ENSMUSG00000010592
# rowData names(9): source type ... mm.entrezIds JAX.geneID
# colnames(21884): 1_AAACCCAAGGTACATA-1 1_AAACCCACATCCGAGC-1 ...
# 4_TTTGTTGCATACAGCT-1 4_TTTGTTGGTCAAACGG-1
# colData names(17): Sample Barcode ... cellType.final cellType.broad
# reducedDimNames(4): GLMPCA_approx UMAP TSNE GLMPCA_50
# mainExpName: NULL
# altExpNames(0):

identical(rownames(colData(sce_mouse_sub)),colnames(sce_mouse_sub))
#[1] TRUE

##Load in DEGs
load(here("processed-data","mouse_1vALL_broad.rda"),verbose = TRUE)
# Loading objects:
#   mouse_1vALL_broad

#Subset for only LS clusters. 
mouse_1vALL_broad <- subset(mouse_1vALL_broad,subset=(cellType.target == "LS"))

dim(mouse_1vALL_broad)
#[1] 32285     9

#Add gene name information 
mouse_1vALL_homol <- merge(mouse_1vALL_broad,
                           as.data.frame(rowData(sce_mouse_sub)[,c("gene_id","JAX.geneID")]),
                           by = "gene_id") 

dim(mouse_1vALL_homol)
#[1] 16574    10

#Change rownames to JAX.geneID
rownames(mouse_1vALL_homol) <- mouse_1vALL_homol$JAX.geneID 

#Alter order of mouse DEGs to be that of human DEGs
mouse_1vALL_homol <- mouse_1vALL_homol[match(human_1vALL_homol$JAX.geneID,
                                             mouse_1vALL_homol$JAX.geneID),]

#Sanity check to make sure the rownames are in the corret order
identical(rownames(human_1vALL_homol),rownames(mouse_1vALL_homol))
#[1] TRUE

#Correlate the std.logFC 
cor(human_1vALL_homol$std.logFC,
    mouse_1vALL_homol$std.logFC)
#[1] 0.5954178
