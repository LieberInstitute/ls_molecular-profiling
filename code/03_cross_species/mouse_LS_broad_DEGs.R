#Goal: Perform 1vALL DEG testing with DeconvoBuddies. LS cells vs all others. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

library(SingleCellExperiment)
library(DeconvoBuddies)
library(sessioninfo)
library(scater)
library(here)


#load the SingleCellExperiment object for mouse Lateral Septum
load(file = here("processed-data","02_build_sce","mouse_sce.rda"),verbose = TRUE)

sce.ls 

stopifnot(identical(rownames(colData(sce.ls)),colnames(sce.ls)))

################
#Need to remove drop.doublet, drop.likelyDoublet, drop.lowNTx, and Neuron.mixed from cell types
#Keep only true cell types. 
sce.ls <- sce.ls[,sce.ls$cellType.final %in% c("Astro","Chol_Ex.D","ChP",
                                               "Endo","Ependymal","IoC_In.E",
                                               "LS_In.C","LS_In.D","LS_In.M",
                                               "LS_In.N","LS_In.O","LS_In.P",
                                               "LS_In.Q","LS_In.R","Micro",
                                               "MS_In.J","MS_In.K","Mural",
                                               "Neuroblast","Oligo","OPC",
                                               "OPC_COP","Sept_In.G","Sept_In.I",
                                               "Str_In.A","Str_In.F","Str_In.H","Str_In.L",
                                               "Thal_Ex.B","TNoS_Ex.A","TT.IG.SH_Ex.C",
                                               "TT.IG.SH_Ex.E","TT.IG.SH_Ex.F")]

sce.ls$Broad_LS <- ifelse(sce.ls$cellType.final %in% c("LS_In.C","LS_In.D","LS_In.M","LS_In.N","LS_In.O","LS_In.P","LS_In.Q","LS_In.R"),"LS","Other")

table(sce.ls$Broad_LS)

stopifnot(identical(rownames(colData(sce.ls)),colnames(sce.ls)))

####################
##Make a colored tsne by LS vs other
LS_cols <- c("red","gray44")
names(LS_cols) <- c("LS","Other")

#Make the tSNE
LS_tSNE_mouse <- plotReducedDim(sce.ls,
                                dimred      = "TSNE",
                                colour_by   = "Broad_LS",
                                point_alpha = 0.3) +
    scale_color_manual(values = LS_cols)
ggsave(filename = here("plots","Conservation","mouse_tSNE_LSvsother.pdf"),plot = LS_tSNE_mouse)


####################
#Run 1vALL DEG testing with DeconvoBuddies function as was done for human data. 
mouse_1vALL_broad <- findMarkers_1vAll(sce.ls, 
                                         assay_name = "logcounts", 
                                         cellType_col = "Broad_LS", 
                                         direction = "up",
                                         mod = "~Sample")

#Add the gene_name information to the markers dataframe. 
colnames(mouse_1vALL_broad)[1] <- "gene_id"
mouse_1vALL_broad <- dplyr::left_join(x = as.data.frame(mouse_1vALL_broad),
                                        y = as.data.frame(rowData(sce.ls)[,c("gene_id","gene_name")]),
                                        by = "gene_id")

#Save the DEG dataframe 
save(mouse_1vALL_broad,
     file = here("processed-data","mouse_1vALL_broad.rda"))
#######################################

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
