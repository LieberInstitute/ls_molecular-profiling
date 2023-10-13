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


Sys.time()
#[1] "2023-10-13 17:06:44 EDT"
sce.ls <- nullResiduals(sce,
                        assay = "counts", 
                        fam   = "binomial", 
                        type  = "deviance")
Sys.time()

