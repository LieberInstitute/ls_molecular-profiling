#Goal: Perform 1vALL DEG testing with DeconvoBuddies as was done with human dataset.
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

library(SingleCellExperiment)
library(DeconvoBuddies)
library(here)


#load the SingleCellExperiment object for mouse Lateral Septum
load(file = "/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/processed_data/SCE/sce_updated_LS.rda",verbose = TRUE)


sce.ls 

stopifnot(identical(rownames(colData(sce)),colnames(sce)))

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


#Refactor the celltypes 
sce.ls$cellType.final <- factor(sce.ls$cellType.final,
                                levels = c("Astro","Chol_Ex.D","ChP",
                                           "Endo","Ependymal","IoC_In.E",
                                           "LS_In.C","LS_In.D","LS_In.M",
                                           "LS_In.N","LS_In.O","LS_In.P",
                                           "LS_In.Q","LS_In.R","Micro",
                                           "MS_In.J","MS_In.K","Mural",
                                           "Neuroblast","Oligo","OPC",
                                           "OPC_COP","Sept_In.G","Sept_In.I",
                                           "Str_In.A","Str_In.F","Str_In.H","Str_In.L",
                                           "Thal_Ex.B","TNoS_Ex.A","TT.IG.SH_Ex.C",
                                           "TT.IG.SH_Ex.E","TT.IG.SH_Ex.F"))


sce.ls 

stopifnot(identical(rownames(colData(sce)),colnames(sce)))

#Run 1vALL DEG testing with DeconvoBuddies function as was done for human data. 
markers_1vALL_mouse <- findMarkers_1vAll(sce.ls, 
                                         assay_name = "logcounts", 
                                         cellType_col = "cellType.final", 
                                         direction = "up",
                                         mod = "~Sample")

#Add the gene_name information to the markers dataframe. 
colnames(markers_1vALL_mouse)[1] <- "gene_id"
markers_1vALL_mouse <- dplyr::left_join(x = as.data.frame(markers_1vALL_mouse),
                                        y = as.data.frame(rowData(sce.ls)[,c("gene_id","gene_name")]),
                                        by = "gene_id")

#Save the DEG dataframe 
save(markers_1vALL_mouse,
     file = here("processed-data","mouse_markers_1vAll.rda"))
#######################################

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

