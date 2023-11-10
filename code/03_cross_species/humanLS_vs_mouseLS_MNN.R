#Goal: Compare gene expression signatures of the mouse and human LS
#This analysis will focus on using MNN to integrate. 
#Code modified from https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/51d15ef9f5f2c4c53f55e22e3fe467de1a724668/10x_NAc-n8_step04_cross-species_rnNAc_MNT.R#L4 and
#https://github.com/LieberInstitute/BLA_crossSpecies/blob/devel/code/costa_rds/05_species_comparisons/crossSpecies_PCA.R
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling

library(SingleCellExperiment)
library(batchelor)
library(harmony)
library(scater)
library(here)
library(scry)

#Load the human and mouse matched SCE objects.
load(here("processed-data","human_mouse_matched_by_JAX.rda"),verbose = TRUE)
# Loading objects:
#     sce_mouse_sub
#     sce_human_sub

sce_mouse_sub
# class: SingleCellExperiment 
# dim: 16572 22860 
# metadata(1): Samples
# assays(3): counts binomial_pearson_residuals logcounts
# rownames(16572): ENSMUSG00000096351 ENSMUSG00000095567 ...
# ENSMUSG00000037772 ENSMUSG00000003526
# rowData names(9): source type ... mm.entrezIds JAX.geneID
# colnames(22860): 1_AAACCCAAGGTACATA-1 1_AAACCCACATCCGAGC-1 ...
# 4_TTTGTTGCATACAGCT-1 4_TTTGTTGGTCAAACGG-1
# colData names(17): Sample Barcode ... cellType.final cellType.broad
# reducedDimNames(4): GLMPCA_approx UMAP TSNE GLMPCA_50
# mainExpName: NULL
# altExpNames(0):

sce_human_sub
# class: SingleCellExperiment 
# dim: 16572 9225 
# metadata(1): Samples
# assays(3): counts binomial_deviance_residuals logcounts
# rownames(16572): ENSG00000187634 ENSG00000188976 ... ENSG00000214026
# ENSG00000100033
# rowData names(9): source type ... hs.entrezIds JAX.geneID
# colnames(9225): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 3_TTTGGTTTCTTCGACC-1 3_TTTGTTGTCCCGATCT-1
# colData names(61): Sample Barcode ... CellType_k_20_louvain
# CellType.Final
# reducedDimNames(18): GLMPCA_approx UMAP_15 ... tSNE_mnn_25 tSNE_mnn_50
# mainExpName: NULL
# altExpNames(0):

#Everything should be in order, but check to make sure? 
all(rowData(sce_human_sub)$JAX.geneID == rowData(sce_mouse_sub)$JAX.geneID) 
# [1] TRUE

#Add coldata column that corresponds to the species. 
sce_human_sub$Species <- "Human"
sce_mouse_sub$Species <- "Mouse"

#Clean up celltype column names. 
sce_human_sub$CellType <- sce_human_sub$CellType.Final
sce_mouse_sub$CellType <- sce_mouse_sub$cellType.final

#subset the objects 
colData(sce_human_sub) <- colData(sce_human_sub)[,c("Sample","Barcode","Species","sum","detected","CellType")]
colData(sce_mouse_sub) <- colData(sce_mouse_sub)[,c("Sample","Barcode","Species","sum","detected","CellType")]

#To combine the objects, concatenate the count data + colData
#Then, simply create the object with SingleCellExperiment()
combo_counts <- cbind(assay(sce_human_sub,"counts"),
                      assay(sce_mouse_sub,"counts"))
combo_coldata <- rbind(colData(sce_human_sub),
                       colData(sce_mouse_sub))

#Create the combo object
sce_combo <- SingleCellExperiment(assays=list(counts=combo_counts), colData=combo_coldata)

sce_combo
# class: SingleCellExperiment 
# dim: 16572 32085 
# metadata(0):
#     assays(1): counts
# rownames(16572): ENSG00000187634 ENSG00000188976 ... ENSG00000214026
# ENSG00000100033
# rowData names(0):
#     colnames(32085): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 4_TTTGTTGCATACAGCT-1 4_TTTGTTGGTCAAACGG-1
# colData names(6): Sample Barcode ... detected CellType
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):
    
#compute logcounts
sce_combo <- multiBatchNorm(sce_combo, batch=sce_combo$Sample)

##########Feature selection. 
#Deviance feature selection
sce_combo <- devianceFeatureSelection(sce_combo,
                                      assay = "counts",
                                      fam = "binomial",
                                      sorted = FALSE,
                                      batch = as.factor(sce_combo$Sample))

pdf(here("plots","Conservation","featureSelxn_binomialDeviance-byGene_numericcutoffs.pdf"))
plot(sort(rowData(sce_combo)$binomial_deviance, decreasing = T),
     type = "l", xlab = "ranked genes",
     ylab = "binomial deviance"
)
abline(v = 2000,lty = 2, col = "red")
dev.off()

#2000 should be good. 
sce_combo <- nullResiduals(sce_combo,
                           assay = "counts", 
                           fam   = "binomial", 
                           type  = "deviance")

#Take top 2000 highly deviant genes
hdgs <- rownames(sce_combo)[order(rowData(sce_combo)$binomial_deviance, decreasing = T)][1:2000]


#current object doesn't include gene name info. Pull from the sce_human_sub object. 
hdgs.symbols <- rowData(sce_human_sub)$gene_name[match(hdgs, rowData(sce_human_sub)$gene_id)]

#PCA
sce_combo <- runPCA(sce_combo,
                    exprs_values = "binomial_deviance_residuals",
                    subset_row = hdgs, 
                    ncomponents = 100,
                    name = "GLMPCA_approx")


#PCA plot of top 6 PCs
PCA_plots <- plotReducedDim(sce_combo,
                            dimred = "GLMPCA_approx", 
                            colour_by = "Sample",
                            ncomponents = 6, 
                            point_alpha = 0.3)
ggsave(PCA_plots,filename = here("plots","Conservation","multi_PCAs.png"))

#Run the UMAP + tSNE with 50 dimensions. 
#UMAP
set.seed(1234)
sce_combo <- runUMAP(sce_combo,
                     dimred = "GLMPCA_approx",
                     n_dimred = 50,
                     name = "UMAP_50")
#tSNE
set.seed(1234)
sce_combo <- runTSNE(sce_combo,
                     dimred = "GLMPCA_approx",
                     n_dimred = 50,
                     name = "tSNE_50")

#Plot the tSNE and UMAP
#UMAP  
umap_cross_species <- plotReducedDim(sce_combo,
                                     dimred = "UMAP_50", 
                                     colour_by = "Species",
                                     point_alpha = 0.3)
ggsave(umap_cross_species,filename = here("plots",
                                          "Conservation",
                                          "umap_cross_species.png"))

#tSNE
tSNE_cross_species <- plotReducedDim(sce_combo,
                                     dimred = "tSNE_50", 
                                     colour_by = "Species",
                                     point_alpha = 0.3)
ggsave(tSNE_cross_species,filename = here("plots",
                                          "Conservation",
                                          "tSNE_cross_species.png"))


#As expected, tSNE and UMAP are split by species. Next run fast fastMNN
mnn.hold <- fastMNN(sce_combo, 
                    batch=sce_combo$Sample,
                    subset.row=hdgs, d=50,
                    correct.all=TRUE, get.variance=TRUE,
                    BSPARAM=BiocSingular::IrlbaParam())

#Add mnn to sce_combo object
reducedDim(sce_combo, "PCA_corrected") <- reducedDim(mnn.hold, "corrected")

#Rerun the UMAP and tSNE
set.seed(1234)
sce_combo <- runUMAP(sce_combo,
                     dimred = "PCA_corrected",
                     n_dimred = 50,
                     name = "UMAP_corrected_50")
#tSNE
set.seed(1234)
sce_combo <- runTSNE(sce_combo,
                     dimred = "PCA_corrected",
                     n_dimred = 50,
                     name = "tSNE_corrected_50")


#Plot the tSNE and UMAP
#UMAP by species
umap_cross_species <- plotReducedDim(sce_combo,
                                     dimred = "UMAP_corrected_50", 
                                     colour_by = "Species",
                                     point_alpha = 0.3)
ggsave(umap_cross_species,filename = here("plots",
                                          "Conservation",
                                          "umap_cross_species_mnn_corrected.png"))

#tSNE
tSNE_cross_species <- plotReducedDim(sce_combo,
                                     dimred = "tSNE_corrected_50", 
                                     colour_by = "Species",
                                     point_alpha = 0.3)
ggsave(tSNE_cross_species,filename = here("plots",
                                          "Conservation",
                                          "tSNE_cross_species_mnn_corrected.png"))

save(sce_combo,file = here("processed-data","sce_combo.rda"))

#By celltype
#First create a new column within sce_combo that is celltype_species
sce_combo$CellType_Species <- paste(sce_combo$CellType,sce_combo$Species,sep = "_")

#Assign colors
cluster_cols <- Polychrome::createPalette(length(unique(sce_combo$CellType_Species)),c("#FF0000", "#00FF00", "#0000FF"))
names(cluster_cols) <- unique(sce_combo$CellType_Species)

#UMAP
umap_cross_species_ct <- plotReducedDim(sce_combo,
                                        dimred = "UMAP_corrected_50", 
                                        colour_by = "CellType_Species",
                                        text_by = "CellType_Species",
                                        point_alpha = 0.3) +
    ggplot2::theme(legend.position = "none")
ggsave(umap_cross_species_ct,filename = here("plots",
                                             "Conservation",
                                             "umap_cross_species_mnn_corrected_byCellType.png"),
       height = 10, width = 10)

#tSNE
tSNE_cross_species_ct <- plotReducedDim(sce_combo,
                                        dimred = "tSNE_corrected_50", 
                                        colour_by = "CellType_Species",
                                        text_by = "CellType_Species",
                                        point_alpha = 0.3) +
    ggplot2::theme(legend.position = "none")
ggsave(tSNE_cross_species_ct,filename = here("plots",
                                             "Conservation",
                                             "tSNE_cross_species_mnn_corrected_byCellType.png"),
       height = 10, width = 10)

#These look okay. Will also try Harmony. 
#Add a duplicate reducedDim that is only called PCA. This is required for RunHarmony 
reducedDim(sce_combo,"PCA") <- reducedDim(sce_combo,"PCA_corrected")
sce_harmony_Species <- RunHarmony(sce_combo, group.by.vars = "Species", verbose = TRUE)
# Transposing data matrix
# Hard k-means centroids initialization
# Harmony 1/10
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# Harmony 2/10
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# Harmony 3/10
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# Harmony converged after 3 iterations

#Run tSNE adn UMAP with harmony
#UMAP
set.seed(1234)
sce_harmony_Species <- runUMAP(sce_harmony_Species,
                               dimred = "HARMONY",
                               name = "UMAP_HARMONY")
#tSNE
set.seed(1234)
sce_harmony_Species <- runTSNE(sce_harmony_Species,
                               dimred = "HARMONY",
                               name = "tSNE_HARMONY")
#Plot the tSNE and UMAP
#UMAP by species
umap_Harmony <- plotReducedDim(sce_harmony_Species,
                               dimred = "UMAP_HARMONY", 
                               colour_by = "Species",
                               point_alpha = 0.3)
ggsave(umap_Harmony,filename = here("plots",
                                    "Conservation",
                                    "umap_cross_species_Harmony_corrected.png"))

#tSNE
tSNE_Harmony <- plotReducedDim(sce_harmony_Species,
                               dimred = "tSNE_HARMONY", 
                               colour_by = "Species",
                               point_alpha = 0.3)
ggsave(tSNE_Harmony,filename = here("plots",
                                    "Conservation",
                                    "tSNE_cross_species_Harmony.png"))





