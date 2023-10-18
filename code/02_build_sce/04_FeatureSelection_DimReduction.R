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
#load(file = here("processed-data","sce_clean.rda"))
load(file = here("processed-data","sce_clean_numeric_cutoffs.rda"))

sce
# class: SingleCellExperiment 
# dim: 36601 10000 
# metadata(1): Samples
# assays(1): counts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(6): source type ... gene_name gene_type
# colnames(10000): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 3_TTTGGTTTCTTCGACC-1 3_TTTGTTGTCCCGATCT-1
# colData names(54): Sample Barcode ... discared_numeric discard_numeric
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):

sce <- devianceFeatureSelection(sce,
                                assay = "counts",
                                fam = "binomial",
                                sorted = FALSE,
                                batch = as.factor(sce$Sample))

pdf(here("plots","featureSelxn_binomialDeviance-byGene_numericcutoffs.pdf"))
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
                     type  = "deviance")
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
#NaNs are produced when all count values for a specific gene are 0. 
#This can be seen at https://github.com/kstreet13/scry/blob/master/R/nullResiduals.R line 87
#To prove it 
table(rowSums(counts(sce)) == 0)
# FALSE  TRUE 
# 33564  3037 
table(rowSums(assay(sce,"binomial_deviance_residuals")) == 0)
# FALSE  TRUE 
# 33564  3037
#Same number of true and false

#Take top 2000 highly deviant genes
hdgs <- rownames(sce)[order(rowData(sce)$binomial_deviance, decreasing = T)][1:2000]
hdgs.symbols <- rowData(sce)$gene_name[match(hdgs, rowData(sce)$gene_id)]

#Run PCA
sce_uncorrected <- runPCA(sce,
                          exprs_values = "binomial_deviance_residuals",
                          subset_row = hdgs, 
                          ncomponents = 100,
                          name = "GLMPCA_approx")

#PCA plot of top 6 PCs
PCA_plots <- plotReducedDim(sce_uncorrected,
                            dimred = "GLMPCA_approx", 
                            colour_by = "Sample",
                            ncomponents = 6, 
                            point_alpha = 0.3)
#ggsave(PCA_plots,filename = here("plots","Dim_Red","multi_PCAs.png"))
ggsave(PCA_plots,filename = here("plots","Dim_Red","multi_PCAs_numericCutoff.png"))

# UMAP
#With testing 15,20,25,and 50
#50 dimensions as in Tran, Maynard et al Neuron
for(i in c(15,20,25,50)){
    print(i)
    set.seed(1234)
    sce_uncorrected <- runUMAP(sce_uncorrected,
                               dimred = "GLMPCA_approx",
                               n_dimred = i,
                               name = paste0("UMAP_",i))
}

# t-SNE
#with 15 dimensions
for(i in c(15,20,25,50)){
    print(i)
    set.seed(1234)
    sce_uncorrected <- runTSNE(sce_uncorrected,
                               dimred = "GLMPCA_approx",
                               n_dimred = i,
                               name = paste0("TSNE_",i))
}

# UMAPs
#Colored by sample, library size, and doublet score
#Sample
for(i in c(15,20,25,50)){
    x <- plotReducedDim(sce_uncorrected,
                        dimred = paste0("UMAP_",i), 
                        colour_by = "Sample",
                        point_alpha = 0.3) +
        ggtitle(paste(i,"Dimensions")) +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(x,filename = here("plots",
                             "Dim_Red",
                             paste0("sample_umap_numericCutoffs_",i,"dims.png")))
}

#####
#library size
#Sample
for(i in c(15,20,25,50)){
    x <- plotReducedDim(sce_uncorrected,
                        dimred = paste0("UMAP_",i), 
                        colour_by = "sum",
                        point_alpha = 0.3) +
        ggtitle(paste(i,"Dimensions")) +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(x,filename = here("plots",
                             "Dim_Red",
                             paste0("sum_umap_numericCutoffs_",i,"dims.png")))
}


#####
#Doublet score
for(i in c(15,20,25,50)){
    x <- plotReducedDim(sce_uncorrected,
                        dimred = paste0("UMAP_",i), 
                        colour_by = "doubletScore",
                        point_alpha = 0.3) +
        ggtitle(paste(i,"Dimensions")) +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(x,filename = here("plots",
                             "Dim_Red",
                             paste0("doublet_umap_numericCutoffs_",i,"dims.png")))
}

######
# TSNE by sample
for(i in c(15,20,25,50)){
    x <- plotReducedDim(sce_uncorrected,
                        dimred = paste0("TSNE_",i), 
                        colour_by = "Sample")
    ggsave(x,filename = here("plots",
                             "Dim_Red",
                             paste0("sample_tSNE_numericCutoffs_",i,"dims.png")))
    
}

#batch effect results in cluster coming from a single sample. 
#Will need to run MNN to fix this. 
#save uncorrected object. 
#save(sce_uncorrected,file = here("processed-data","sce_uncorrected.rda"))
save(sce_uncorrected,file = here("processed-data","sce_uncorrected_numericCutoffs.rda"))


#Run batch correction with mutual nearest neighbors. 
glmpca_mnn <- batchelor::reducedMNN(reducedDim(sce_uncorrected, "GLMPCA_approx"),
                                    batch=as.factor(sce_uncorrected$Sample))

#Add mnn to the object
reducedDim(sce_uncorrected,"mnn") <- glmpca_mnn$corrected

#Rename object
sce <- sce_uncorrected
rm(sce_uncorrected)

#umap with 15,20,25,50 dimensions
#Only going to run umap at this point. 
for(i in c(15,20,25,50)){
    set.seed(1234)
    sce <- runUMAP(sce,
                   dimred = "mnn",
                   name = paste0("UMAP_mnn_",i))
}


#Plot umaps by sample, sum, detected, and subsets_mito_percent
#Sample
for(i in c(15,20,25,50)){
    x <- plotReducedDim(sce_uncorrected,
                        dimred = paste0("UMAP_",i), 
                        colour_by = "Sample",
                        point_alpha = 0.3) +
        ggtitle(paste(i,"Dimensions")) +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(x,filename = here("plots",
                             "Dim_Red",
                             paste0("sample_umap_numericCutoffs_",i,"dims_postMNN.png")))
}

#sum
for(i in c(15,20,25,50)){
    x <- plotReducedDim(sce_uncorrected,
                        dimred = paste0("UMAP_",i), 
                        colour_by = "sum",
                        point_alpha = 0.3) +
        ggtitle(paste(i,"Dimensions")) +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(x,filename = here("plots",
                             "Dim_Red",
                             paste0("sum_umap_numericCutoffs_",i,"dims_postMNN.png")))
}

#detected
for(i in c(15,20,25,50)){
    x <- plotReducedDim(sce_uncorrected,
                        dimred = paste0("UMAP_",i), 
                        colour_by = "detected",
                        point_alpha = 0.3) +
        ggtitle(paste(i,"Dimensions")) +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(x,filename = here("plots",
                             "Dim_Red",
                             paste0("detected_umap_numericCutoffs_",i,"dims_postMNN.png")))
}

#subsets mito percent
for(i in c(15,20,25,50)){
    x <- plotReducedDim(sce_uncorrected,
                        dimred = paste0("UMAP_",i), 
                        colour_by = "detected",
                        point_alpha = 0.3) +
        ggtitle(paste(i,"Dimensions")) +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(x,filename = here("plots",
                             "Dim_Red",
                             paste0("detected_umap_numericCutoffs_",i,"dims_postMNN.png")))
}




#One cluster is dominated by Sample 1. Checking expression values to identify if this is a 
#batch correction issue, or if this is due to the 

#Comoute log counts to plot expression
sce <- batchelor::multiBatchNorm(sce, batch = sce$Sample)

genes <- c("SYT1","SNAP25", #pan neuron
           "MBP","MOBP", #OLIGODENDROCYTE
           "CD74", "CSF1R", "C3", #MICROGLIA
           "GFAP", "TNC", "AQP4", "SLC1A2", #ASTROCYTEs
           "GAD1","GAD2","SLC32A1",#Pan GABA
           "SLC17A7", "SLC17A6", "SLC17A8",
           "TRPC4","HOMER2","PTPN3", #Mouse LS markers
           "ELAVL2", #Mouse LS markers
           "CRHR1","CRHR2", 
           "OXTR","AVPR1A", 
           "DRD3")

for(i in genes){
    print(i)
    x <- plotReducedDim(sce,
                        dimred = "UMAP_mnn", 
                        colour_by = i,
                        swap_rownames = "gene_name") +
        scale_color_gradientn(colours = c("lightgrey","red")) +
        ggtitle(i) +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(filename = paste0("plots/Expression_plots/",i,"_expression_umap_numericCutoffs.pdf"),
           plot = x,
           height = 8,width = 8)
}



#Save the object
save(sce,file = here("processed-data","sce_numeric_cutoff_QC.rda"))

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
