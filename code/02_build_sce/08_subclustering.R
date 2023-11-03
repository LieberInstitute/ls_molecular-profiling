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
library(ComplexHeatmap)
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
set.seed(1234)
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
set.seed(1234)
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

# LS_1_2 - '2023-11-03 11:57:47.398536
# LS_2_5 - '2023-11-03 11:58:32.625525
# Glutamatergic - '2023-11-03 11:59:16.967022
# Polydendrocyte - '2023-11-03 12:00:01.764599
# LS_2_6 - '2023-11-03 12:00:46.003714
# Microglia - '2023-11-03 12:01:31.023445
# LS_1_1 - '2023-11-03 12:02:15.986136
# Oligo - '2023-11-03 12:03:00.906465
# Ependymal - '2023-11-03 12:03:45.591706
# Striosome - '2023-11-03 12:04:29.166364
# LS_1_4 - '2023-11-03 12:05:13.30039
# MS_GABA_2 - '2023-11-03 12:05:57.406989
# LS_1_3 - '2023-11-03 12:06:41.428329
# MS_GABA_1 - '2023-11-03 12:07:25.143667
# LS_2_4 - '2023-11-03 12:08:08.894753
# LS_1_5 - '2023-11-03 12:08:52.729599
# Endothelial - '2023-11-03 12:09:36.202296
# Astrocyte - '2023-11-03 12:10:20.199761
# LS_1_6 - '2023-11-03 12:11:03.626639
# LS_2_3 - '2023-11-03 12:11:47.699489
# LS_2_2 - '2023-11-03 12:12:30.810959
# LS_2_1 - '2023-11-03 12:13:14.611372
# LS_2_7 - '2023-11-03 12:13:58.633656
# LS_2_8 - '2023-11-03 12:14:42.805443
# Building Table - 2023-11-03 12:15:24.505502
# ** Done! **

#Add gene symbols. 
colnames(OnevAll_enriched_sub)[1] <- "gene_id"
OnevAll_enriched_sub <- dplyr::left_join(x = as.data.frame(OnevAll_enriched_sub),
                                         y = as.data.frame(rowData(sce)[,c("gene_id","gene_name")]),
                                         by = "gene_id")

save(OnevAll_enriched_sub,file = here("processed-data","OnevAll_enriched_subclustered.rda"))

#Add an FDR column 
OnevAll_enriched_sub$FDR <- exp(OnevAll_enriched_sub$log.FDR)

#Subset for significant FDR + logFC >= 0.25 only
OnevAll_enriched_sig <- subset(OnevAll_enriched_sub,subset=(log.FDR <= log(0.001) & logFC >= 0.25))


#One major question about the current clustering is whether the LS_GABA_2
#cluster contains a mixture of striatal and Lateral Septum cell types.
#TRPC4 + DGKG are good markers of LS, while RARB and BCL11B are good markers of striatum. 
#Best way is to determine which populations contain sets of DEGs that map to each region. 
#For example, classifying a population that contains GAD1, GAD2, TRPC4, DGKG as DEGs as LS is much more 
#convincing than just TRPC4 itself. One gene does not define a cell type. 
#Will use Vipr2, Cpa6, Lgr5, Baiap3, Pdcd7, Trpc5 as broad septal markers.
#To define MS markers, find shared DEGs within the two MS clusters. 
Top100_MS_1 <- subset(OnevAll_enriched_sub,subset=(cellType.target == "MS_GABA_1"))[1:100,"gene_name"]
Top100_MS_2 <- subset(OnevAll_enriched_sub,subset=(cellType.target == "MS_GABA_2"))[1:100,"gene_name"]
MS_genes <- intersect(Top100_MS_1,Top100_MS_2)
MS_genes
# [1] "KCNC2"      "ELAVL2"     "ANK1"       "GRIN2D"     "AL033504.1"
# [6] "SAMD5"      "KLHL1"      "AL445623.2" "ELAVL4"     "AL445123.2"
# [11] "NHS"        "RANBP17"    "AFF2"       "RAB3B"      "FLT3"      
# [16] "LHX6"       "SRRM4" 
#AL445623.2 is a lncRNA for ELAVL2, interesting. 

#First try to make binary matrices that tell whether each gene 
DEG_mat_list <- vector(mode = "list",length = 18)
names(DEG_mat_list) <- c(paste0("LS_1_",1:6),
                         paste0("LS_2_",1:8),
                         "Glutamatergic",
                         "MS_GABA_1",
                         "MS_GABA_2",
                         "Striosome")
for(i in names(DEG_mat_list)){
    DEG_mat_list[[i]] <- matrix(ncol = 39,nrow = 39)
    colnames(DEG_mat_list[[i]]) <- c("GAD1", "GAD2", "TRPC4", "DGKG", "HOMER2", "PTPN3", "NRP1", "TRHDE", #LS
                                     "FOXP2", "RARB", "BCL11B","PPP1R1B","OPRM1", "DRD1", "DRD2", "PENK", #Str
                                     "VIPR2", "CPA6", "LGR5", "BAIAP3", "PDCD7", "TRPC5", #Broad septal
                                     MS_genes) 
    rownames(DEG_mat_list[[i]]) <- colnames(DEG_mat_list[[i]])
}

#Fill the matrix
for(i in names(DEG_mat_list)){
    print(i)
    x <- subset(OnevAll_enriched_sig,subset=(cellType.target == i))
    for(j in colnames(DEG_mat_list[[i]])){
        for(l in rownames(DEG_mat_list[[i]])){
            DEG_mat_list[[i]][l,j] <- ifelse(nrow(subset(x,subset=(gene_name == j | gene_name == l))) == 2,
                                             1,
                                             0)
        }
    }
}

#Create a row annotation dataframe for the heatmaps.
gene_celltype_markers <- c(rep("GABA",2),
                           rep("LS",6),
                           rep("Str",8),
                           rep("Sept",6),
                           rep("MS",17))
row_anno <- rowAnnotation(CellType = gene_celltype_markers)

#Generate the heatmaps. 
for(i in names(DEG_mat_list)){
    hm <- Heatmap(DEG_mat_list[[i]],
                  cluster_rows = FALSE,cluster_columns = FALSE,col = c("white","red"),
                  rect_gp = gpar(col = "black", lwd = 1),
                  right_annotation = row_anno)
    pdf(file = here("plots","DEG_matrices",paste0(i,"_DEG_Matrix.pdf")),
        height = 8,
        width = 8)
    draw(hm,column_title=i)
    dev.off()
}




