#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling
#module load conda_R/4.4

library(SingleCellExperiment)
library(sessioninfo)
library(scater)
library(here)


#load the object
load(here("processed-data","02_build_sce","sce_celltype.rda"),verbose = TRUE)


for(i in c("DGKG","CRHR2","ELAVL2",
           "RARB","BCL11B","DRD1",
           "DRD2","OPRM1","CRYM")){
  x <- plotReducedDim(sce,
                      dimred = "tSNE_mnn_50",
                      colour_by = i,
                      swap_rownames = "gene_name") +
    scale_color_gradientn(colours = c("lightgrey","orange","red")) +
    theme_void() 
  ggsave(x,filename = here("plots","Plots_for_Supp",
                           paste0(i,"_FeaturePlot_Legend.pdf")))
}



print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
