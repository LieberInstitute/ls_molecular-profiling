#Goal: compile droplet scores, calculate QC metrics, and detect doublets. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/
#code modified from https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/03_reducedDimensions_clustering.R

library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(scater)
library(scran)
library(scry)
library(here)

## load QCed and cleaned object. 
load(file = here("processed-data","02_build_sce","sce_clean.rda"),verbose = TRUE)

dim(sce)

identical(rownames(colData(sce)),colnames(sce))

identical(rownames(colData(sce)),colnames(counts(sce)))

sce

#Deviance feature selection
set.seed(1234)
sce <- devianceFeatureSelection(sce,
                                assay = "counts",
                                fam = "binomial",
                                sorted = FALSE,
                                batch = as.factor(sce$Sample))

pdf(here("plots","featureSelxn_binomialDeviance-byGene.pdf"))
plot(sort(rowData(sce)$binomial_deviance, decreasing = T),
     type = "l", xlab = "ranked genes",
     ylab = "binomial deviance"
)
abline(v = 2000,lty = 2, col = "red")
dev.off()


#2000 should be good. 
message("Running null residuals")
Sys.time()
sce <- nullResiduals(sce,
                     assay = "counts", 
                     fam   = "binomial", 
                     type  = "pearson")

#Take top 2000 highly deviant genes
hdgs <- rownames(sce)[order(rowData(sce)$binomial_deviance, decreasing = T)][1:2000]
hdgs.symbols <- rowData(sce)$gene_name[match(hdgs, rowData(sce)$gene_id)]

#Run PCA
message("Running PCA")
Sys.time()
set.seed(1234)
sce_uncorrected <- runPCA(sce,
                          exprs_values = "binomial_pearson_residuals",
                          subset_row = hdgs, 
                          ncomponents = 100,
                          name = "GLMPCA_approx")

#PCA plot of top 6 PCs
PCA_plots <- plotReducedDim(sce_uncorrected,
                            dimred = "GLMPCA_approx", 
                            colour_by = "Sample",
                            ncomponents = 6, 
                            point_alpha = 0.3)
ggsave(PCA_plots,filename = here("plots","Dim_Red","multi_PCAs.png"))


# tSNE
#50 dimensions as in Tran, Maynard et al Neuron
message("Running tSNE-uncorrected")
Sys.time()
set.seed(1234)
sce_uncorrected <- runTSNE(sce_uncorrected,
                           dimred = "GLMPCA_approx",
                           n_dimred = 50,
                           name = "tSNE_50")

##Plot the TSNE by Sample, library size and Doublet score
#Sample 
sample_tSNE_uncorrected <- plotReducedDim(sce_uncorrected,
                                          dimred = "tSNE_50", 
                                          colour_by = "Sample",
                                          point_alpha = 0.3) +
  ggtitle("50 Dimensions") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(sample_tSNE_uncorrected,
       filename = here("plots",
                       "Dim_Red",
                       "Sample_TSNE_50dimensions.png"))

#Library Size
libsize_tSNE_uncorrected <- plotReducedDim(sce_uncorrected,
                                           dimred = "tSNE_50", 
                                           colour_by = "sum",
                                           point_alpha = 0.3) +
  ggtitle("50 Dimensions") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(libsize_tSNE_uncorrected,
       filename = here("plots",
                       "Dim_Red",
                       "libsize_TSNE_50dimensions.png"))

#Doublet Score
dubscore_tSNE_uncorrected <- plotReducedDim(sce_uncorrected,
                                            dimred = "tSNE_50", 
                                            colour_by = "doubletScore",
                                            point_alpha = 0.3) +
  ggtitle("50 Dimensions") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(dubscore_tSNE_uncorrected,
       filename = here("plots",
                       "Dim_Red",
                       "DoubletScore_TSNE_50dimensions.png"))

#batch effect results in cluster coming from a single sample. 
#Will need to run MNN to fix this. 
#save uncorrected object. 
save(sce_uncorrected,file = here("processed-data","02_build_sce","sce_uncorrected.rda"))

#Run batch correction with mutual nearest neighbors. 
message("Running mnn")
Sys.time()
set.seed(1234)
glmpca_mnn <- batchelor::reducedMNN(reducedDim(sce_uncorrected, "GLMPCA_approx"),
                                    k=20,
                                    batch=as.factor(sce_uncorrected$Sample))

identical(rownames(colData(sce_uncorrected)),rownames(glmpca_mnn))

#Add mnn to the object
reducedDim(sce_uncorrected,"mnn") <- glmpca_mnn$corrected

#Rename object
sce <- sce_uncorrected
rm(sce_uncorrected)

#tSNE in MNN space
message("Running tSNE post-mnn")
Sys.time()
set.seed(1234)
sce <- runTSNE(sce,
               dimred = "mnn",
               n_dimred = 50,
               name = "tSNE_mnn_50")

##Plot the TSNE by Sample, library size and Doublet score
#Sample 
sample_tSNE_corrected <- plotReducedDim(sce,
                                        dimred = "tSNE_mnn_50", 
                                        colour_by = "Sample",
                                        point_alpha = 0.3) + 
  ggtitle("50 Dimensions\nmnn corrected") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(sample_tSNE_corrected,
       filename = here("plots",
                       "Dim_Red",
                       "Sample_mnn_corrected_tSNE_50dimensions.png"))

#Library Size
libsize_tSNE_corrected <- plotReducedDim(sce,
                                         dimred = "tSNE_mnn_50", 
                                         colour_by = "sum",
                                         point_alpha = 0.3) +
  ggtitle("50 Dimensions\nmnn corrected") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(libsize_tSNE_corrected,
       filename = here("plots",
                       "Dim_Red",
                       "libsize_mnn_corrected_tSNE_50dimensions.png"))

#Doublet Score
dubscore_tSNE_corrected <- plotReducedDim(sce,
                                          dimred = "tSNE_mnn_50", 
                                          colour_by = "doubletScore",
                                          point_alpha = 0.3) +
  ggtitle("50 Dimensions\nmnn corrected") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(dubscore_tSNE_corrected,
       filename = here("plots",
                       "Dim_Red",
                       "DoubletScore_mnn_corrected_tSNE_50dimensions.png"))


#detected (number of genes)
detected_plot <- plotReducedDim(sce,
                                dimred = "tSNE_mnn_50", 
                                colour_by = "detected",
                                point_alpha = 0.3) +
  ggtitle("50 Dimensions\nmnn corrected") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(detected_plot,filename = here("plots",
                                     "Dim_Red",
                                     "detected_mnn_corrected_tSNE_50dimensions.png"))

#subsets mito percent
mito_plot <- plotReducedDim(sce,
                            dimred = "tSNE_mnn_50", 
                            colour_by = "subsets_Mito_percent",
                            point_alpha = 0.3) +
  ggtitle("50 Dimensions\nmnn corrected") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(mito_plot,filename = here("plots",
                                 "Dim_Red",
                                 "mito_percentage_mnn_corrected_tSNE_50dimensions.png"))

########
#One cluster is dominated by Sample 1. Checking expression values to identify if this is a 
#batch correction issue, or if this is due to something about the anatomy/biology of the sample

#Compute log counts to plot expression.
#Solely for the purpose of exploration at this point. 
sce <- batchelor::multiBatchNorm(sce, batch = sce$Sample)

genes <- c("SYT1","SNAP25", #pan neuron
           "MBP","MOBP", #OLIGODENDROCYTE
           "CD74", "CSF1R", "C3", #MICROGLIA
           "GFAP", "TNC", "AQP4", "SLC1A2", #ASTROCYTEs
           "GAD1","GAD2","SLC32A1",#Pan GABA
           "SLC17A7", "SLC17A6", "SLC17A8", #GLut markers
           "TRPC4","HOMER2","PTPN3", #Mouse LS markers
           "OPRM1","DRD1","DRD2","CRYM", #Striatal markers for good measure
           "ELAVL2", #Mouse LS markers
           "CRHR1","CRHR2", 
           "OXTR","AVPR1A", 
           "DRD3")

#Check cluster express
for(i in genes){
  print(i)
  x <- plotReducedDim(sce,
                      dimred = "tSNE_mnn_50", 
                      colour_by = i,
                      swap_rownames = "gene_name") +
    scale_color_gradientn(colours = c("lightgrey","red")) +
    ggtitle(i) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename = here("plots","Expression_plots","FeaturePlots_tSNE_50dims",paste0(i,"_expression_tSNE_50_dims.png")),
         plot = x,
         height = 8,width = 8)
}

#This cluster is the result of glutamatergic cells present in sample 1 only. 

#Save the object
save(sce,file = here("processed-data","02_build_sce","sce_postMNN.rda"))

sce

#session info
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
