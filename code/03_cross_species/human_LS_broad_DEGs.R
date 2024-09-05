#Goal: Perform 1vALL DEG testing with DeconvoBuddies. Human LS neurons vs all other human cells. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

library(SingleCellExperiment)
library(DeconvoBuddies)
library(sessioninfo)
library(scater)
library(here)

#load the SingleCellExperiment object for human Lateral Septum
load(here("processed-data","02_build_sce","sce_celltype.rda"),verbose = TRUE)

sce

identical(rownames(colData(sce)),colnames(sce))

stopifnot(identical(rownames(colData(sce)),colnames(sce)))

############
#Label cell types as broad LS or not. 
sce$Broad_LS <- ifelse(sce$CellType.Final %in% c("LS_Inh_A","LS_Inh_B","LS_Inh_C","LS_Inh_D","LS_Inh_E"),
                       "LS",
                       "Other")

table(sce$Broad_LS)

############
##Make a colored tsne by LS vs other
LS_cols <- c("red","gray44")
names(LS_cols) <- c("LS","Other")

#Make the tSNE
LS_tSNE <- plotReducedDim(sce,
                          dimred      = "tSNE_mnn_50",
                          colour_by   = "Broad_LS",
                          point_alpha = 0.3) +
  scale_color_manual(values = LS_cols)
ggsave(filename = here("plots","Conservation","human_tSNE_LSvsother.pdf"),plot = LS_tSNE)

############
##DEG Testing
#Run 1vALL DEG testing with DeconvoBuddies function as was done for human data. 
human_1vALL_broad <- findMarkers_1vAll(sce, 
                                       assay_name = "logcounts", 
                                       cellType_col = "Broad_LS", 
                                       direction = "up",
                                       mod = "~Sample")

#Add the gene_name information to the markers dataframe. 
colnames(human_1vALL_broad)[1] <- "gene_id"
human_1vALL_broad <- dplyr::left_join(x = as.data.frame(human_1vALL_broad),
                                      y = as.data.frame(rowData(sce)[,c("gene_id","gene_name")]),
                                      by = "gene_id")

#Save the DEG dataframe 
save(human_1vALL_broad,
     file = here("processed-data","human_1vALL_broad.rda"))
#######################################
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

