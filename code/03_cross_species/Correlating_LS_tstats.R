#Goal: Correlate t-statistics from all LS clusters. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

library(SingleCellExperiment)
library(DeconvoBuddies)
library(sessioninfo)
library(ggplot2)
library(scater)
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

#Generate a tSNE where LS clusters are colored red and everything else is gray
#color coding for LS clusters to be red and everything else to be gray
LS_cols <- c("red","gray44")
names(LS_cols) <- c("LS","Other")

#Make the tSNE
LS_tSNE <- plotReducedDim(sce,
                          dimred      = "tSNE_mnn_50",
                          colour_by   = "LS_vs_other",
                          point_alpha = 0.3) +
    scale_color_manual(values = LS_cols)
ggsave(filename = here("plots","Conservation","human_tSNE_LSvsother.pdf"),plot = LS_tSNE)

#Find markers for the merged LS cluster
h_LS_clusters <- findMarkers_1vAll(sce,
                                   assay_name   = "logcounts",
                                   cellType_col = "LS_vs_other",
                                   mod          = "~Sample")
# LS - '2024-01-25 14:45:57.152908
# Other - '2024-01-25 14:46:09.434334
# Building Table - 2024-01-25 14:46:21.040722
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


#Generate a tSNE where LS clusters are colored red and everything else is gray
#color coding for LS clusters to be red and everything else to be gray
LS_cols <- c("red","gray44")
names(LS_cols) <- c("LS","Other")

#Make the tSNE
LS_tSNE_mouse <- plotReducedDim(sce.ls,
                                dimred      = "TSNE",
                                colour_by   = "LS_vs_other",
                                point_alpha = 0.3) +
    scale_color_manual(values = LS_cols)
ggsave(filename = here("plots","Conservation","mouse_tSNE_LSvsother.pdf"),plot = LS_tSNE_mouse)


#Run 1vALL DEG for mouse LS 
m_LS_clusters <- findMarkers_1vAll(sce.ls,
                                   assay_name   = "logcounts",
                                   cellType_col = "LS_vs_other",
                                   mod          = "~Sample")
# LS - '2024-01-25 14:50:40.16095
# Other - '2024-01-25 14:51:09.498538
# Building Table - 2024-01-25 14:51:38.845341
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

cor.test(all_DEGs_homol[,"t.stat_human"],
         all_DEGs_homol[,"t.stat_mouse"])
# Pearson's product-moment correlation
# data:  all_DEGs_homol[, "t.stat_human"] and all_DEGs_homol[, "t.stat_mouse"]
# t = 96.156, df = 16560, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.5887099 0.6082579
# sample estimates:
#      cor 
# 0.598573 

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

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# [1] "Reproducibility information:"
# [1] "2024-01-25 14:54:53 EST"
#    user  system elapsed 
# 233.528   8.031 829.548
# ─ Session info ────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 Patched (2023-07-19 r84711)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-01-25
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/pandoc
# 
# ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [2] CRAN (R 4.3.1)
# beachmat               2.16.0    2023-04-25 [2] Bioconductor
# beeswarm               0.4.0     2021-06-01 [2] CRAN (R 4.3.1)
# Biobase              * 2.60.0    2023-04-25 [2] Bioconductor
# BiocGenerics         * 0.46.0    2023-04-25 [2] Bioconductor
# BiocNeighbors          1.18.0    2023-04-25 [2] Bioconductor
# BiocParallel           1.34.2    2023-05-22 [2] Bioconductor
# BiocSingular           1.16.0    2023-04-25 [2] Bioconductor
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.1)
# bluster                1.10.0    2023-04-25 [2] Bioconductor
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.1)
# cluster                2.1.4     2022-08-22 [3] CRAN (R 4.3.1)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.3.1)
# colorout             * 1.3-0.1   2023-12-01 [1] Github (jalvesaq/colorout@deda341)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.3.1)
# cowplot                1.1.1     2020-12-30 [2] CRAN (R 4.3.1)
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.1)
# DeconvoBuddies       * 0.99.0    2023-12-04 [1] Github (LieberInstitute/DeconvoBuddies@9ce4a42)
# DelayedArray           0.26.7    2023-07-28 [2] Bioconductor
# DelayedMatrixStats     1.22.6    2023-08-28 [2] Bioconductor
# dplyr                  1.1.3     2023-09-03 [2] CRAN (R 4.3.1)
# dqrng                  0.3.1     2023-08-30 [2] CRAN (R 4.3.1)
# edgeR                  3.42.4    2023-05-31 [2] Bioconductor
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.3.1)
# farver                 2.1.1     2022-07-06 [2] CRAN (R 4.3.1)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.3.1)
# GenomeInfoDb         * 1.36.3    2023-09-07 [2] Bioconductor
# GenomeInfoDbData       1.2.10    2023-07-20 [2] Bioconductor
# GenomicRanges        * 1.52.0    2023-04-25 [2] Bioconductor
# ggbeeswarm             0.7.2     2023-04-29 [2] CRAN (R 4.3.1)
# ggplot2              * 3.4.3     2023-08-14 [2] CRAN (R 4.3.1)
# ggrepel                0.9.3     2023-02-03 [2] CRAN (R 4.3.1)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.3.1)
# gridExtra              2.3       2017-09-09 [2] CRAN (R 4.3.1)
# gtable                 0.3.4     2023-08-21 [2] CRAN (R 4.3.1)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.3.1)
# igraph                 1.5.1     2023-08-10 [2] CRAN (R 4.3.1)
# IRanges              * 2.34.1    2023-06-22 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.3.1)
# labeling               0.4.3     2023-08-29 [2] CRAN (R 4.3.1)
# lattice                0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.3.1)
# limma                  3.56.2    2023-06-04 [2] Bioconductor
# locfit                 1.5-9.8   2023-06-11 [2] CRAN (R 4.3.1)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.3.1)
# Matrix                 1.6-1.1   2023-09-18 [3] CRAN (R 4.3.1)
# MatrixGenerics       * 1.12.3    2023-07-30 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.3.1)
# metapod                1.8.0     2023-04-25 [2] Bioconductor
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.3.1)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.3.1)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.3.1)
# purrr                  1.0.2     2023-08-10 [2] CRAN (R 4.3.1)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.1)
# rafalib                1.0.0     2015-08-09 [1] CRAN (R 4.3.1)
# ragg                   1.2.5     2023-01-12 [2] CRAN (R 4.3.1)
# RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.3.1)
# Rcpp                   1.0.11    2023-07-06 [2] CRAN (R 4.3.1)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.1)
# rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.1)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.3.1)
# S4Arrays               1.0.6     2023-08-30 [2] Bioconductor
# S4Vectors            * 0.38.1    2023-05-02 [2] Bioconductor
# ScaledMatrix           1.8.1     2023-05-03 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.3.1)
# scater               * 1.28.0    2023-04-25 [2] Bioconductor
# scran                  1.28.2    2023-07-23 [2] Bioconductor
# scuttle              * 1.10.2    2023-08-03 [2] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.1)
# SingleCellExperiment * 1.22.0    2023-04-25 [2] Bioconductor
# sparseMatrixStats      1.12.2    2023-07-02 [2] Bioconductor
# statmod                1.5.0     2023-01-06 [2] CRAN (R 4.3.1)
# stringi                1.7.12    2023-01-11 [2] CRAN (R 4.3.1)
# stringr                1.5.0     2022-12-02 [2] CRAN (R 4.3.1)
# SummarizedExperiment * 1.30.2    2023-06-06 [2] Bioconductor
# systemfonts            1.0.4     2022-02-11 [2] CRAN (R 4.3.1)
# textshaping            0.3.6     2021-10-13 [2] CRAN (R 4.3.1)
# tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.3.1)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.3.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.3.1)
# vctrs                  0.6.3     2023-06-14 [2] CRAN (R 4.3.1)
# vipor                  0.4.5     2017-03-22 [2] CRAN (R 4.3.1)
# viridis                0.6.4     2023-07-22 [2] CRAN (R 4.3.1)
# viridisLite            0.4.2     2023-05-02 [2] CRAN (R 4.3.1)
# withr                  2.5.0     2022-03-03 [2] CRAN (R 4.3.1)
# XVector                0.40.0    2023-04-25 [2] Bioconductor
# zlibbioc               1.46.0    2023-04-25 [2] Bioconductor
# 
# [1] /users/rphillip/R/4.3
# [2] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library
# 
# ───────────────────────────────────────────────────────────────────────────────────────────────────────