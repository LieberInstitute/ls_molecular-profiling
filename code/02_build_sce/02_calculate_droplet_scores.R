#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling
#Code modified from https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/02_get_droplet_scores.R
#and https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/master/10x_all-FACS-n10_2021rev_step01_processing-QC_MNT.R

library(SingleCellExperiment)
library(DropletUtils)
library(scuttle)
library(tidyverse)
library(here)
library(sessioninfo)
library(ggplot2)

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
    e.out <- DropletUtils::emptyDrops(sce_sub,
                                      niters = 30000,
                                      lower = knee_lower)
    print("Done - saving data")
    Sys.time()
    #save the droplet data. 
    save(e.out,here("processed-data","02_build_sce","droplet_scores",paste0(i,"_droplet_scores",".Rdata")))
    
    #Generate plots to check the cutoff value. 
    FDR_cutoff <- 0.001
    addmargins(table(Signif = e.out$FDR <= FDR_cutoff, Limited = e.out$Limited, useNA = "ifany"))
    
    n_cell_anno <- paste("Non-empty:", sum(e.out$FDR < FDR_cutoff, na.rm = TRUE))
    message(n_cell_anno)
    
    my_theme <- theme_bw() +
        theme(text = element_text(size = 15))
    
    droplet_elbow_plot <- as.data.frame(bcRanks) %>%
        add_column(FDR = e.out$FDR) %>%
        ggplot(aes(x = rank, y = total, color = FDR < FDR_cutoff)) +
        geom_point(alpha = 0.5, size = 1) +
        geom_hline(yintercept = metadata(bcRanks)$knee, linetype = "dotted", color = "gray") +
        annotate("text", x = 10, y = metadata(bcRanks)$knee, label = "Knee", vjust = -1, color = "gray") +
        geom_hline(yintercept = knee_highest, linetype = "dashed") +
        annotate("text", x = 10, y = knee_highest, label = "Knee est 'highest'") +
        geom_hline(yintercept = knee_higher, linetype = "dashed") +
        annotate("text", x = 10, y = knee_higher, label = "Knee est 'higher'") +
        geom_hline(yintercept = knee_lower, linetype = "dashed") +
        annotate("text", x = 10, y = knee_lower, label = "Knee est 'lower'") +
        geom_hline(yintercept = knee_lowest, linetype = "dashed") +
        annotate("text", x = 10, y = knee_lowest, label = "Knee est 'lowest'") +
        scale_x_continuous(trans = "log10") +
        scale_y_continuous(trans = "log10") +
        labs(
            x = "Barcode Rank",
            y = "Total UMIs",
            title = paste("Sample", sample_run),
            subtitle = n_cell_anno,
            color = paste("FDR <", FDR_cutoff)
        ) +
        my_theme +
        theme(legend.position = "bottom")
    
    ggsave(droplet_elbow_plot, 
           filename = here("plots","02_build_sce", "droplet_scores", paste0(i,"droplet_qc_",i,".png")))
}


#session info
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

