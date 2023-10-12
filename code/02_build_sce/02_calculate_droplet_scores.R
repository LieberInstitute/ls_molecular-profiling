#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling
#Code modified from https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/02_get_droplet_scores.R
#and https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/master/10x_all-FACS-n10_2021rev_step01_processing-QC_MNT.R

library(SingleCellExperiment)
library(DropletUtils)
library(scuttle)
library(tidyverse)
library(here)
library(sessioninfo)

#Load the sce object 
load(here("processed-data","sce_raw.rda"),verbose = TRUE)

for(i in unique(sce$Sample)){
    #First subset the sample
    print(paste("Subsetting for sample",i))
    sce_sub <- sce[,sce$Sample == i]
    
    #Rank barcodes
    print(paste("Ranking barcodes for sample",i))
    bcRanks <- barcodeRanks(sce_sub, fit.bounds = c(10,1e3))
    
    knee_highest <- metadata(bcRanks)$knee - 200
    message(
        "'First knee point' = ", metadata(bcRanks)$knee, "\n",
        "knee_highest =", knee_highest
    )
    
    knee_higher <- metadata(bcRanks)$knee - 100
    message(
        "'Second knee point' = ", metadata(bcRanks)$knee, "\n",
        "knee_higher =", knee_higher
    )
    
    message(
        "'Third knee point' = ", metadata(bcRanks)$knee, "\n",
        "knee =", metadata(bcRanks)$knee
    )
    
    knee_lower <- metadata(bcRanks)$knee + 100
    message(
        "'Fourth knee point' = ", metadata(bcRanks)$knee, "\n",
        "knee_lower =", knee_lower
    )
    
    knee_lowest <- metadata(bcRanks)$knee + 200
    message(
        "'Fifth knee point' = ", metadata(bcRanks)$knee, "\n",
        "knee_lowest =", knee_lowest
    )
    
    #Run emptyDrops with knee + 100 first
    set.seed(1234)
    print("Starting emptyDrops")
    Sys.time()
    e.out <- DropletUtils::emptyDrops(sce_sample,
                                      niters = 30000,
                                      lower = knee_lower)
    print("Done - saving data")
    Sys.time()
    
}




