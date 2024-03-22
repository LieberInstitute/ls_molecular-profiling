library("SingleCellExperiment")
library("iSEE")
library("shiny")
library("paletteer")
library("scuttle")
library("SpatialExperiment")

load("sce_clean.rda", verbose = TRUE)

## Don't run this on app.R since we don't want to run this every single time
#lobstr::obj_size(spe_pseudo)




source("initial.R", print.eval = TRUE)


#rse_gene <- registerAppOptions(rse_gene, color.maxlevels = length(Sample_ID)
iSEE(
  spe_pseudo,
  appTitle = "LS snRNAseq data",
  initial = initial,
  colormap = ExperimentColorMap(colData = list(
    domain = function(n) {
      cols <- spatial.palette
    },
    broad.domain = function(n) {
      cols <- spatial.palette2
    }
    
  ))
)
