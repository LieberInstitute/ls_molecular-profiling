#Goal: Investigate and cluster the integrated object. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling

library(SingleCellExperiment)
library(sessioninfo)
library(scater)
library(here)

#Load the integrated object. 
load(here("processed-data","sce_harmony_species.rda"))