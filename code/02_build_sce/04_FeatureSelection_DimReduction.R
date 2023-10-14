#Goal: compile droplet scores, calculate QC metrics, and detect doublets. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/
#code modified from https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/03_reducedDimensions_clustering.R

library(SingleCellExperiment)
library(ggplot2)
library(scater)
library(scran)
library(scry)
library(here)

## load QCed and cleaned object. 
load(file = here("processed-data","sce_clean.rda"))

sce
# class: SingleCellExperiment 
# dim: 36601 9883 
# metadata(1): Samples
# assays(1): counts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(6): source type ... gene_name gene_type
# colnames(9883): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 3_TTTGTTGAGGCTCCCA-1 3_TTTGTTGTCCCGATCT-1
# colData names(44): Sample Barcode ... discard_auto doubletScore
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):

#Run deviance feature selection with default parameters. 
sce <- devianceFeatureSelection(sce,assay = "counts",fam = "binomial",sorted = FALSE,batch = as.factor(sce$Sample))

pdf(here("plots","featureSelxn_binomialDeviance-byGene.pdf"))
plot(sort(rowData(sce)$binomial_deviance, decreasing = T),
     type = "l", xlab = "ranked genes",
     ylab = "binomial deviance"
)
abline(v = 2000,lty = 2, col = "red")
dev.off()

#2000 should be good. 
sce <- nullResiduals(sce,
                     assay = "counts", 
                     fam   = "binomial", 
                     type  = "deviance",
                     batch = as.factor(sce$Sample))
# In addition: Warning messages:
#     1: In asMethod(object) :
#     sparse->dense coercion: allocating vector of size 1.0 GiB
# 2: In asMethod(object) :
#     sparse->dense coercion: allocating vector of size 1.0 GiB
# 3: In .sparse2dense(x) :
#     sparse->dense coercion: allocating vector of size 1.0 GiB
# 4: In asMethod(object) :
#     sparse->dense coercion: allocating vector of size 1.0 GiB
# 5: In asMethod(object) :
#     sparse->dense coercion: allocating vector of size 1.0 GiB
# 6: In sqrt(x@x) : NaNs produced
#Not sure about warnings. sparse-->dense coercion doesn't seem to affect function. 
#No binomical deviances are NAs and  #ull residual matrix seems to be alright. 
#Need to figure out why NANs are being produced and why this could be happening? 
#In the meantime, all values are accounted for. Will move forward. 


#Take top 2000 highly deviant genes
hdgs <- rownames(sce)[order(rowData(sce)$binomial_deviance, decreasing = T)][1:2000]
hdgs.symbols <- rowData(sce)$gene_name[match(hdgs, rowData(sce)$gene_id)]

#Run PCA
sce_uncorrected <- runPCA(sce,
                          exprs_values = "binomial_deviance_residuals",
                          subset_row = hdgs, 
                          ncomponents = 100,
                          name = "GLMPCA_approx")

# UMAP
set.seed(1234)
sce_uncorrected <- runUMAP(sce_uncorrected,
                           dimred = "GLMPCA_approx",
                           n_dimred = 50, 
                           name = "UMAP")
# t-SNE
set.seed(1234)
sce_uncorrected <- runTSNE(sce_uncorrected,
                           dimred = "GLMPCA_approx",
                           n_dimred = 50, 
                           name = "TSNE")

#PCA plot of top 10 PCs
PCA_plots <- plotReducedDim(sce_uncorrected,
               dimred = "GLMPCA_approx", 
               colour_by = "Sample",
               ncomponents = 6, 
               point_alpha = 0.3)
ggsave(PCA_plots,filename = here("plots","Dim_Red","multi_PCAs.png"))

# UMAPs
#Colored by sample, library size, and doublet score
#Sample
sample_umap <- plotReducedDim(sce_uncorrected,
                              dimred = "UMAP", 
                              colour_by = "Sample",
                              point_alpha = 0.3)
ggsave(sample_umap,filename = here("plots","Dim_Red","sample_umap.png"))

#library size
sum_umap <- plotReducedDim(sce_uncorrected,
               dimred = "UMAP", colour_by = "sum",
               point_alpha = 0.3)
ggsave(sum_umap,filename = here("plots","Dim_Red","sum_umap.png"))

#Doublet score
doublet_umap <- plotReducedDim(sce_uncorrected,
                               dimred = "UMAP", 
                               colour_by = "doubletScore")
ggsave(doublet_umap,filename = here("plots","Dim_Red","doublet_umap.png"))

# TSNE by sample
sample_TSNE <- plotReducedDim(sce_uncorrected,
                              dimred = "TSNE", 
                              colour_by = "Sample")
ggsave(sample_TSNE,filename = here("plots","Dim_Red","sample_TSNE.png"))

#batch effect results in cluster coming from a single sample. 
#Will need to run MNN to fix this. 
#save uncorrected object. 
save(sce_uncorrected,file = here("processed-data","sce_uncorrected.rda"))


