#Goal: Correlate t-statistics from all LS clusters. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

library(SingleCellExperiment)
library(DeconvoBuddies)
library(sessioninfo)
library(ggplot2)
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
                         rowData(sce_human_sub)[,c("gene_id","JAX.geneID")],
                         by = "gene_id") 

#Add Jax.GeneID info to the mouse LS DEGs and also remove any genes that are not homologus
m_LS_DEGs_homol <- merge(m_LS_DEGs,
                         rowData(sce_mouse_sub)[,c("gene_id","JAX.geneID")],
                         by = "gene_id") 

#To correlate, each dataframe needs to be in the same order. 
#First, make the JAX.geneID the rownames for each dataframe
rownames(h_LS_DEGs_homol) <- h_LS_DEGs_homol$JAX.geneID
rownames(m_LS_DEGs_homol) <- m_LS_DEGs_homol$JAX.geneID

#Alter order of mouse DEGs to be that of human DEGs
fixTo <- rownames(h_LS_DEGs_homol)
m_LS_DEGs_homol <- m_LS_DEGs_homol[fixTo,]

#Sanity check to make sure the rownames are in the corret order
all(rownames(h_LS_DEGs_homol) == rownames(m_LS_DEGs_homol))
#[1] TRUE

#Correlate 
cor(h_LS_DEGs_homol$t.stat,
    m_LS_DEGs_homol$t.stat)
#[1] 0.598573

#Correlate the two dataframes to identify genes that are shared markers and those that are divergent markers
#To do this merge the two dataframes
colnames(h_LS_DEGs_homol)[c(2,6,7,9,10)] <- paste0(colnames(h_LS_DEGs_homol)[c(2,6,7,9,10)],"_human")
colnames(m_LS_DEGs_homol)[c(2,6,7,9,10)] <- paste0(colnames(m_LS_DEGs_homol)[c(2,6,7,9,10)],"_mouse")

#Merge the two dataframes
all_DEGs_homol <- merge(x = h_LS_DEGs_homol[,c(2,6,7,9:11)],
                        y = m_LS_DEGs_homol[,c(2,6,7,9:11)],
                        by = "JAX.geneID")


#Save the all_DEGs_homol dataframe
save(all_DEGs_homol,file = here("processed-data","Human_mouse_homologous_DEGs_and_tstats.rda"))

#Make the plot
all_LS_plot <- ggplot(data=as.data.frame(all_DEGs_homol),aes(x = t.stat_human,y = t.stat_mouse)) + 
    geom_point(alpha=0.5) +
    xlim(c(-350,350)) +
    ylim(c(-600,600)) +
    geom_hline(yintercept = 0,lty = 2) +
    geom_vline(xintercept = 0,lty = 2) +
    labs(x = "Human LS t-statistic", 
         y = "Mouse LS t-statistic") +
    theme_bw() +
    geom_text(data = as.data.frame(subset(all_DEGs_homol,subset=(t.stat_mouse >= 300 & t.stat_human >= 100))),
              label = as.data.frame(subset(all_DEGs_homol,subset=(t.stat_mouse >= 300 & t.stat_human >= 100)))$gene_name_human,
              nudge_y = 15,
              nudge_x = -15) +
    geom_text(data = as.data.frame(subset(all_DEGs_homol,subset=(gene_name_human == "FREM2"))),
              label = as.data.frame(subset(all_DEGs_homol,subset=(gene_name_human == "FREM2")))$gene_name_human,
              nudge_y = 15,
              nudge_x = -15) +
    geom_text(data = as.data.frame(subset(all_DEGs_homol,subset=(t.stat_mouse <= (-150) & t.stat_human <= (-100)))),
              label = as.data.frame(subset(all_DEGs_homol,subset=(t.stat_mouse <= (-150) & t.stat_human <= (-100))))$gene_name_human,
              nudge_y = 15,
              nudge_x = -15) +
    annotate("text",x = 275,y = 600,label = "Human Enriched\nMouse Enriched") +
    annotate("text",x = 275,y = -600, label = "Human Enriched\nMouse Depleted") +
    annotate("text",x = -275,y = 600,label = "Human Depleted\nMouse Enriched") +
    annotate("text",x = -275,y = -600,label = "Human Depleted\nMouse Depleted") 
ggsave(filename = here("plots","Conservation","all_LS_homologs_tstat_correlation.pdf"),plot = all_LS_plot)


#Get the top 100 identifiers for the LS in each species
h_top250 <- subset(h_LS_DEGs_homol,subset=(rank_marker %in% 1:250))
m_top250 <- subset(m_LS_DEGs_homol,subset=(rank_marker %in% 1:250))

#Find the shared identifiers
shared_identifiers <- intersect(rownames(h_top250),
                                rownames(m_top250))
length(shared_identifiers)
# [1] 99

#correlate
cor(h_top250[shared_identifiers,"t.stat"],
    m_top250[shared_identifiers,"t.stat"])
#[1] 0.6182991

#Combine the two dataframes to generate a plot. 
#To do this, it's critical to have column names that tell us something about the speices
colnames(h_top250)[c(2,6,7,9,10)] <- paste0(colnames(h_top250)[c(2,6,7,9,10)],"_human")
colnames(m_top250)[c(2,6,7,9,10)] <- paste0(colnames(m_top250)[c(2,6,7,9,10)],"_mouse")

#subset both dataframes
h_top250_sub <- h_top250[,c(2,6,7,9:11)]
m_top250_sub <- m_top250[,c(2,6,7,9:11)]


#merge
top250_sub <- merge(x  = h_top250_sub,
                    y  = m_top250_sub,
                    by = "JAX.geneID")







