#Several goals for this analysis: 
#1. Investigate Astrocyte_1 and identify if those are in fact Ependymal cells. 
#2. Collapse glial populations that are currently split. 
#3. Subcluster LS_GABA_1 + LS_GABA_2 to assess if any additional heterogeneity within human LS neurons exists.
#Rationale for subclustering:
#Rationale for subclustering LS_GABA_1: Cross-species analysis demonstrated that LS_GABA_1 was highly correlated
#with several of the additional LS clusters found in mouse. 
#Rationale for subclustering LS_GABA_2: Large majority of LS_GABA_2 seems to be striatal cells due to the expression of 
#RARB + Bcl11b. However, one group of cells within that population does not express much RARB and instead 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

library(SingleCellExperiment)
library(ggplot2)
library(scater)
library(scran)
library(here)

#load the SingleCellExperiment object
load(here("processed-data","sce_with_CellType.rda"),verbose = TRUE)

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
# colData names(60): Sample Barcode ... k_75_walktrap CellType
# reducedDimNames(14): GLMPCA_approx UMAP_15 ... UMAP_mnn_25 UMAP_mnn_50
# mainExpName: NULL
# altExpNames(0):

########Goal 1: Investigate Astrocyte_1 and identify if those are in fact Ependymal cells. 
#load in the mouse DEG list to get the top ependymal markers. 
load(here("processed-data","mouse_markers_1vAll_ttest_withnon0median.rda"),verbose = TRUE) #markers.ls.t.1vAll

ependymal_marker_plot <- plotExpression(object = sce,
                                        features = c("FOXJ1","PIFO","CFAP44",#Literature markers
                                                     "DNAH12","WDR49","FRMPD2"), #Top 3 markers from mouse paper
                                        x = "CellType",
                                        swap_rownames = "gene_name",
                                        ncol = 3) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
ggsave(plot = ependymal_marker_plot,filename = here("plots","Expression_plots","ependymal_markers.png"))

#Astrocyte_1 expresses many canonical ependymal markers. Will rename this cluster Ependymal. 
sce$CellType <- as.character(sce$CellType)
sce$CellType[sce$CellType == "Astrocyte_1"] <- "Ependymal"


########Goal 2: Collapse glial populations that are currently split. 
#Merge 3 oligodendrocyte populations.
sce$CellType[sce$CellType %in% c("Oligo_1","Oligo_2","Oligo_3")] <- "Oligo"
#Merge 2 Astrocyte populations. 
sce$CellType[sce$CellType %in% c("Astrocyte_2","Astrocyte_3")] <- "Astrocyte"

######Goal 3: Subcluster LS_GABA_1 and LS_GABA_2
#Subcluster LS_GABA_1 population using k=20 + walktrap
#Force CellType column to a character
sce$CellType_sub <- as.character(sce$CellType)
sce$CellType_sub[sce$CellType_sub == "LS_GABA_1"] <- paste0("LS_1_",as.character(
    clusterCells(x = sce[,sce$CellType == "LS_GABA_1"],
                 use.dimred = "mnn",
                 BLUSPARAM = bluster::NNGraphParam(k=20))
))
#6 subpopulations from LS_GABA_1
table(sce[,sce$CellType == "LS_GABA_1"]$CellType_sub)
# LS_1_1 LS_1_2 LS_1_3 LS_1_4 LS_1_5 LS_1_6 
# 277    253    295    403     48    123 


#Subcluster LS_GABA_1 population using k=20 + walktrap
sce$CellType_sub[sce$CellType_sub == "LS_GABA_2"] <- paste0("LS_2_",as.character(
    clusterCells(x = sce[,sce$CellType == "LS_GABA_2"],
                 use.dimred = "mnn",
                 BLUSPARAM = bluster::NNGraphParam(k=20))
))
#8 subpopulations from LS_GABA_2
table(sce[,sce$CellType == "LS_GABA_2"]$CellType_sub)
# LS_2_1 LS_2_2 LS_2_3 LS_2_4 LS_2_5 LS_2_6 LS_2_7 LS_2_8 
# 177    205    159    175    282    333     56     16 


#Run 1vALL Marker testing to evaluate whether we are oversplitting these clusters. 
library(DeconvoBuddies)
OnevAll_enriched_sub <- findMarkers_1vAll(sce, 
                                          assay_name = "logcounts", 
                                          cellType_col = "CellType_sub", 
                                          mod = "~Sample")
# LS_1_2 - '2023-11-02 15:36:20.294838
# LS_2_5 - '2023-11-02 15:36:36.370232
# Glutamatergic - '2023-11-02 15:36:52.564961
# Polydendrocyte - '2023-11-02 15:37:08.601832
# LS_2_6 - '2023-11-02 15:37:24.664279
# Microglia - '2023-11-02 15:37:40.803494
# LS_1_1 - '2023-11-02 15:37:56.949976
# Oligo - '2023-11-02 15:38:13.119035
# Ependymal - '2023-11-02 15:38:29.240786
# Striosome - '2023-11-02 15:38:46.39432
# LS_1_4 - '2023-11-02 15:39:02.650818
# MS_GABA_2 - '2023-11-02 15:39:18.864675
# LS_1_3 - '2023-11-02 15:39:35.022745
# MS_GABA_1 - '2023-11-02 15:39:51.268955
# LS_2_4 - '2023-11-02 15:40:07.387008
# LS_1_5 - '2023-11-02 15:40:23.518827
# Endothelial - '2023-11-02 15:40:39.707098
# Astrocyte - '2023-11-02 15:40:55.823201
# LS_1_6 - '2023-11-02 15:41:11.89402
# LS_2_3 - '2023-11-02 15:41:28.159031
# LS_2_2 - '2023-11-02 15:41:44.206101
# LS_2_1 - '2023-11-02 15:42:00.238676
# LS_2_7 - '2023-11-02 15:42:16.253879
# LS_2_8 - '2023-11-02 15:42:32.275946
# Building Table - 2023-11-02 15:42:48.276712
# ** Done! **

#Add gene symbols. 
colnames(OnevAll_enriched_sub)[1] <- "gene_id"
OnevAll_enriched_sub <- dplyr::left_join(x = as.data.frame(OnevAll_enriched_sub),
                                         y = as.data.frame(rowData(sce)[,c("gene_id","gene_name")]),
                                         by = "gene_id")

#Add an FDR column 
OnevAll_enriched_sub$FDR <- exp(OnevAll_enriched_sub$log.FDR)

#Subset for significant FDR only
OnevAll_enriched_sig <- subset(OnevAll_enriched_sub,subset=(FDR <= 0.001))

#One major question about the current clustering is whether the LS_GABA_2
#cluster contains a mixture of striatal and Lateral Septum cell types.
#TRPC4 + DGKG are good markers of LS, while RARB and BCL11B are good markers of striatum. 
#Best way is to determine which populations contain sets of DEGs that map to each region. 
#For example, classifying a population that contains GAD1, GAD2, TRPC4, DGKG as DEGs as LS is much more 
#convincing than just TRPC4 itself. One gene does not define a cell type. 
#GAD1, GAD2, TRPC4, DGKG, HOMER2, PTPN3, NRP1, TRHDE
#FOXP2, RARB, BCL11B, OPRM1, DRD1, DRD2, PENK, ADORA2A
DEG_mat_list <- vector(mode = "list",length = length(unique(sce$CellType_sub)))
names(DEG_mat_list) <- unique(sce$CellType_sub)
for(i in names(DEG_mat_list)){
    DEG_mat_list[[i]] <- matrix(ncol = 22,nrow = 22)
    colnames(DEG_mat_list[[i]]) <- c("GAD1", "GAD2", "TRPC4", "DGKG", "HOMER2", "PTPN3", "NRP1", "TRHDE", #LS
                                     "FOXP2", "RARB", "BCL11B","PPP1R1B","OPRM1", "DRD1", "DRD2", "PENK", #Str
                                     "VIPR2", "CPA6", "LGR5", "BAIAP3", "PDCD7", "TRPC5")
    rownames(DEG_mat_list[[i]]) <- colnames(DEG_mat_list[[i]])
}

#Fill the matrix
for(i in names(DEG_mat_list)){
    print(i)
    x <- subset(OnevAll_enriched_sig,subset=(cellType.target == i))
    for(j in colnames(DEG_mat_list[[i]])){
        for(l in rownames(DEG_mat_list[[i]])){
            if(nrow(subset(x,subset=(gene_name == j | gene_name == l))) > 1){
                DEG_mat_list[[i]][l,j] <- 1
            }else{
                DEG_mat_list[[i]][l,j] <- 0
            }
        }
    }
}






