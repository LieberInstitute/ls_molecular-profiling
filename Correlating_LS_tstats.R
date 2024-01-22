#Goal: Correlate t-statistics from all LS clusters. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

library(SingleCellExperiment)
library(DeconvoBuddies)
library(sessioninfo)
library(here)

#Load the SingleCellExperiment object for human Lateral Septum
load(here("processed-data","sce_with_CellType.rda"),verbose = TRUE)
# Loading objects:
#     sce

sce
# class: SingleCellExperiment 
# dim: 33556 9225 
# metadata(1): Samples
# assays(3): counts binomial_deviance_residuals logcounts
# rownames(33556): ENSG00000243485 ENSG00000238009 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(9225): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 3_TTTGGTTTCTTCGACC-1 3_TTTGTTGTCCCGATCT-1
# colData names(61): Sample Barcode ... CellType_k_20_louvain
# CellType.Final
# reducedDimNames(18): GLMPCA_approx UMAP_15 ... tSNE_mnn_25 tSNE_mnn_50
# mainExpName: NULL
# altExpNames(0):

#Mark LS clusters as just "LS" without subtype designation 
sce$LS_vs_other <- ifelse(sce$CellType.Final %in% c("LS_Inh_A","LS_Inh_B","LS_Inh_G","LS_Inh_I"),
                          "LS",
                          "Other")

#Find markers for the merged LS cluster
h_LS_clusters <- findMarkers_1vAll(sce,
                                   assay_name   = "logcounts",
                                   cellType_col = "LS_vs_other",
                                   mod          = "~Sample")
# LS - '2024-01-22 13:44:51.101302
# Other - '2024-01-22 13:45:07.978748
# Building Table - 2024-01-22 13:45:24.865171
# ** Done! **

#Keep only the LS
h_LS_DEGs <- subset(h_LS_clusters,subset=(cellType.target == "LS"))

#Add symbol information to the table
#First change the ensembl gene id column to have same name as what is in rowData(sce)
colnames(h_LS_DEGs)[1] <- "gene_id"
h_LS_DEGs <- dplyr::left_join(x = as.data.frame(h_LS_DEGs),
                              y = as.data.frame(rowData(sce)[,c("gene_id","gene_name")]),
                              by = "gene_id")

#Calcualte t statistics for human
h_LS_DEGs$t.stat <- h_LS_DEGs$std.logFC * sqrt(ncol(sce))

##load the SingleCellExperiment object for mouse Lateral Septum
load(file = "/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/processed_data/SCE/sce_updated_LS.rda",verbose = TRUE)
# Loading objects:
#     sce.ls
#     annotationTab.ls
#     cell_colors.ls

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

#Mark LS clusters as just "LS" without subtype designation 
sce.ls$LS_vs_other <- ifelse(sce.ls$cellType.final %in% c("LS_In.C","LS_In.D","LS_In.M",
                                                          "LS_In.N","LS_In.O","LS_In.P",
                                                          "LS_In.Q","LS_In.R"),
                             "LS",
                             "Other")

#Run 1vALL DEG for mouse LS 
m_LS_clusters <- findMarkers_1vAll(sce.ls,
                                   assay_name   = "logcounts",
                                   cellType_col = "LS_vs_other",
                                   mod          = "~Sample")
# LS - '2024-01-22 14:36:22.898585
# Other - '2024-01-22 14:37:47.913619
# Building Table - 2024-01-22 14:39:11.516143
# ** Done! **

#Subset for only LS clusters. 
m_LS_DEGs <- subset(m_LS_clusters,subset=(cellType.target == "LS"))

#Add gene name information 
colnames(m_LS_DEGs)[1] <- "gene_id"
m_LS_DEGs <- dplyr::left_join(x = as.data.frame(m_LS_DEGs),
                              y = as.data.frame(rowData(sce.ls)[,c("gene_id","gene_name")]),
                              by = "gene_id")

#Calculate t statistic for mouse data
m_LS_DEGs$t.stat <- m_LS_DEGs$std.logFC * sqrt(ncol(sce.ls))


#Load the human and mouse matched SCE objects.
load(here("processed-data","human_mouse_matched_by_JAX.rda"),verbose = TRUE)
# Loading objects:
#     sce_mouse_sub
#     sce_human_sub
#The objects above are subsetted to contain only genes that are homologous between the two species

#Add Jax.GeneID info to the human LS DEGs. 
#This also removes any genes from the DEG table that are not homologous
h_LS_DEGs_homol <- merge(h_LS_DEGs,
                         rowData(sce_human_sub)[,c("gene_id","gene_name","JAX.geneID")],
                         by = "gene_id") 

#Add Jax.GeneID info to the mouse LS DEGs and also remove any genes that are not homologus
m_LS_DEGs_homol <- merge(m_LS_DEGs,
                         rowData(sce_mouse_sub)[,c("gene_id","gene_name","JAX.geneID")],
                         by = "gene_id") 



