#Downstream analysis of MapMyCells output
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling
#processed-data/MapMyCells_Output/'h_ls_anndata_10xWholeMouseBrain(CCN20230722)_CorrelationMapping_UTC_1705602688838'

library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(here)

#Load human LS sce object
load(here("processed-data","sce_with_CellType.rda"),verbose = TRUE)
# Loading objects:
#     sce

#Read in the MapMyCells output
#First hierarchial mapping. 
Map_Cells_Out_hi <- read.delim(file = here("processed-data",
                                           "MapMyCells_Output",
                                           'h_ls_anndata_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1705605534631',
                                           'h_ls_anndata_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1705605534631.csv'),
                               comment.char = "#",
                               sep = ",")

hi_output <- merge(x = colData(sce)[,c("unique_rowname","CellType.Final")],
                   y = Map_Cells_Out_hi,
                   by.x = "unique_rowname",
                   by.y = "cell_id")








