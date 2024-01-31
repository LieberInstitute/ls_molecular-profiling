#Downstream analysis of MapMyCells output
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling

library(SingleCellExperiment)
library(sessioninfo)
library(pheatmap)
library(ggplot2)
library(here)

#Load human LS sce object
load(here("processed-data","sce_with_CellType.rda"),verbose = TRUE)
# Loading objects:
#     sce

#Read in the MapMyCells output
#First hierarchial mapping. 
Map_Cells_Out_hi <- read.delim(file = here("processed-data",
                                           "MapMyCells_Output",
                                           'h_ls_anndata_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1705605534631',
                                           'h_ls_anndata_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1705605534631.csv'),
                               comment.char = "#",
                               sep = ",")

hi_output <- merge(x = colData(sce)[,c("unique_rowname","CellType.Final")],
                   y = Map_Cells_Out_hi,
                   by.x = "unique_rowname",
                   by.y = "cell_id")

#Bargraph to see where each of the LS clusters are mapping to
LS_cells <- hi_output[grep("LS",hi_output$CellType.Final),]

#Create an empty dataframe and add proportion data to it via for loop
LS_cells_props <- data.frame(CellType.Final  = NA,
                             Mapped_CellType = NA,
                             Freq = NA,
                             Prop = NA)

for(i in unique(LS_cells$CellType.Final)){
    prop_df <- as.data.frame(table(subset(LS_cells,subset=(CellType.Final == i))$class_name))
    prop_df$CellType.Final <- i
    prop_df$Prop <- prop_df$Freq/sum(prop_df$Freq)*100
    colnames(prop_df)[1] <- "Mapped_CellType"
    LS_cells_props <- rbind(LS_cells_props,prop_df[,c("CellType.Final","Mapped_CellType","Freq","Prop")])
}

#remove the NA row 
LS_cells_props <- LS_cells_props[!is.na(LS_cells_props$CellType.Final),]

#Make cell colors. 
cluster_cols <- Polychrome::createPalette(length(unique(LS_cells_props$Mapped_CellType)),
c("#D81B60", "#1E88E5","#FFC107","#009E73"))
names(cluster_cols) <- unique(LS_cells_props$Mapped_CellType)

#Make the plot. 
LS_cluster_mapping <- ggplot(data = LS_cells_props,aes(x = CellType.Final,y = Prop, fill = Mapped_CellType)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cluster_cols) +
    theme_bw() +
    labs(x = "Human LS Cell Type",
         y = "Proportion of Human Cell Type",
         fill = "Mapped Mouse\nCell Type")

ggsave(filename = here("plots","Conservation","MapMyCells","Human_LS_neurons_mapped.pdf"),
       plot = LS_cluster_mapping)

#Plot a heatmap that will determine the proportion of each human cell type represented by each 
#mapped mouse celltypes. 
#Build a matrix where the columns are the human cell types and the rows are the mapped mouse classes
mapped_mat <- matrix(ncol = length(unique(hi_output$CellType.Final)),
                     nrow = length(unique(hi_output$class_name)))

dim(mapped_mat)
#[1] 24 22

#add col and rownames
colnames(mapped_mat) <- unique(hi_output$CellType.Final)
rownames(mapped_mat) <- unique(hi_output$class_name)

#Build the matrix
for(i in colnames(mapped_mat)){
    x <- as.data.frame(table(subset(hi_output,subset=(CellType.Final == i))$class_name))
    x$Prop <- x$Freq/sum(x$Freq)*100
    rownames(x) <- x$Var1
    for(l in rownames(mapped_mat)){
        if(l %in% rownames(x)){
            mapped_mat[l,i] <- x[l,"Prop"]
        }else{
            mapped_mat[l,i] <- 0
        }
    }
}

#Plot the heatmap and save 
pdf(here("plots","Conservation","MapMyCells","Human_to_mouse_classname_heatmap.pdf"))
pheatmap(mat = mapped_mat,
         color = colorRampPalette(c("gray86", "#FFC107", "#D81B60"))(200))
dev.off()
       
#Make a heatmap of the subclass as well. 
mapped_subclass <- matrix(ncol = length(unique(LS_cells$CellType.Final)),
                          nrow = length(unique(LS_cells$subclass_name)))

dim(mapped_subclass)
#[1] 34  4

#add col and rownames
colnames(mapped_subclass) <- unique(LS_cells$CellType.Final)
rownames(mapped_subclass) <- unique(LS_cells$subclass_name)

#Build the matrix
for(i in colnames(mapped_subclass)){
    x <- as.data.frame(table(subset(LS_cells,subset=(CellType.Final == i))$subclass_name))
    x$Prop <- x$Freq/sum(x$Freq)*100
    rownames(x) <- x$Var1
    for(l in rownames(mapped_subclass)){
        if(l %in% rownames(x)){
            mapped_subclass[l,i] <- x[l,"Prop"]
        }else{
            mapped_subclass[l,i] <- 0
        }
    }
}

#Plot the heatmap and save 
pdf(here("plots","Conservation","MapMyCells","LS_Only_Subclass_name.pdf"))
pheatmap(mat = t(mapped_subclass),
         color = colorRampPalette(c("gray86", "#FFC107", "#D81B60"))(200),
         cluster_rows = FALSE)
dev.off()


print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# [1] "Reproducibility information:"
# [1] "2024-01-31 15:42:58 EST"
#    user   system  elapsed 
# 110.561    5.586 1422.261 
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 Patched (2023-07-19 r84711)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-01-31
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [2] CRAN (R 4.3.1)
# Biobase              * 2.60.0    2023-04-25 [2] Bioconductor
# BiocGenerics         * 0.46.0    2023-04-25 [2] Bioconductor
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.1)
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.1)
# colorout             * 1.3-0.1   2023-12-01 [1] Github (jalvesaq/colorout@deda341)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.3.1)
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.1)
# DelayedArray           0.26.7    2023-07-28 [2] Bioconductor
# dplyr                  1.1.3     2023-09-03 [2] CRAN (R 4.3.1)
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.3.1)
# farver                 2.1.1     2022-07-06 [2] CRAN (R 4.3.1)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.3.1)
# GenomeInfoDb         * 1.36.3    2023-09-07 [2] Bioconductor
# GenomeInfoDbData       1.2.10    2023-07-20 [2] Bioconductor
# GenomicRanges        * 1.52.0    2023-04-25 [2] Bioconductor
# ggplot2              * 3.4.3     2023-08-14 [2] CRAN (R 4.3.1)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.3.1)
# gtable                 0.3.4     2023-08-21 [2] CRAN (R 4.3.1)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.3.1)
# IRanges              * 2.34.1    2023-06-22 [2] Bioconductor
# labeling               0.4.3     2023-08-29 [2] CRAN (R 4.3.1)
# lattice                0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.3.1)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.3.1)
# Matrix                 1.6-1.1   2023-09-18 [3] CRAN (R 4.3.1)
# MatrixGenerics       * 1.12.3    2023-07-30 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.3.1)
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.3.1)
# pheatmap             * 1.0.12    2019-01-04 [2] CRAN (R 4.3.1)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.3.1)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.3.1)
# Polychrome             1.5.1     2022-05-03 [1] CRAN (R 4.3.1)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.1)
# ragg                   1.2.5     2023-01-12 [2] CRAN (R 4.3.1)
# RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.3.1)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.1)
# rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.1)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
# S4Arrays               1.0.6     2023-08-30 [2] Bioconductor
# S4Vectors            * 0.38.1    2023-05-02 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.3.1)
# scatterplot3d          0.3-44    2023-05-05 [1] CRAN (R 4.3.1)
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.1)
# SingleCellExperiment * 1.22.0    2023-04-25 [2] Bioconductor
# SummarizedExperiment * 1.30.2    2023-06-06 [2] Bioconductor
# systemfonts            1.0.4     2022-02-11 [2] CRAN (R 4.3.1)
# textshaping            0.3.6     2021-10-13 [2] CRAN (R 4.3.1)
# tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.3.1)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.3.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.3.1)
# vctrs                  0.6.3     2023-06-14 [2] CRAN (R 4.3.1)
# withr                  2.5.0     2022-03-03 [2] CRAN (R 4.3.1)
# XVector                0.40.0    2023-04-25 [2] Bioconductor
# zlibbioc               1.46.0    2023-04-25 [2] Bioconductor
# 
# [1] /users/rphillip/R/4.3
# [2] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────
