#Downstream analysis of MapMyCells output
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling
#module load conda_R/4.4

library(SingleCellExperiment)
library(sessioninfo)
library(pheatmap)
library(ggplot2)
library(here)

#Load human LS sce object
load(here("processed-data","02_build_sce","sce_celltype.rda"),verbose = TRUE)
# Loading objects:
#   sce


#Read in the MapMyCells output
#First hierarchial mapping. 
Map_Cells_Out_hi <- read.delim(file = here("processed-data",
                                           "MapMyCells_Output",
                                           "h_ls_anndata_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1726506607196",
                                           "h_ls_anndata_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1726506607196.csv"),
                               comment.char = "#",
                               sep = ",")

hi_output <- merge(x = colData(sce)[,c("key","CellType.Final")],
                   y = Map_Cells_Out_hi,
                   by.x = "key",
                   by.y = "cell_id")

####Write out hi_ouput as a csv to serve as a supplementary table
write.csv(x = hi_output,file = here("processed-data","MapMyCells_Output","MapMyCells_Output.csv"))

#Bargraph to see where each of the LS clusters are mapping to
LS_cells <- hi_output[grep("LS",hi_output$CellType.Final),]

#Create an empty dataframe and add proportion data to it via for loop
LS_cells_props <- data.frame(CellType.Final  = NA,
                             Mapped_CellType = NA,
                             Freq = NA,
                             Prop = NA)

for(i in unique(LS_cells$CellType.Final)){
  print(i)
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
#[1] 24 25

#add col and rownames
colnames(mapped_mat) <- unique(hi_output$CellType.Final)
rownames(mapped_mat) <- unique(hi_output$class_name)

#Build the matrix
for(i in colnames(mapped_mat)){
  x <- as.data.frame(table(subset(hi_output,subset=(CellType.Final == i))$class_name))
  x$Prop <- (x$Freq/sum(x$Freq))*100
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
#[1] 31  5

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
# [1] "2024-09-16 14:20:46 EDT"
# user   system  elapsed 
# 31.455    2.368 3785.283
# ─ Session info ─────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.4.0 Patched (2024-05-22 r86590)
# os       Rocky Linux 9.4 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-09-16
# pandoc   3.1.13 @ /jhpce/shared/community/core/conda_R/4.4/bin/pandoc
# 
# ─ Packages ─────────────────────────────────────────────────────────────────
# package              * version date (UTC) lib source
# abind                  1.4-5   2016-07-21 [2] CRAN (R 4.4.0)
# Biobase              * 2.64.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# BiocGenerics         * 0.50.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# cli                    3.6.3   2024-06-21 [2] CRAN (R 4.4.0)
# colorspace             2.1-1   2024-07-26 [2] CRAN (R 4.4.0)
# crayon                 1.5.3   2024-06-20 [2] CRAN (R 4.4.0)
# DelayedArray           0.30.1  2024-05-07 [2] Bioconductor 3.19 (R 4.4.0)
# dplyr                  1.1.4   2023-11-17 [2] CRAN (R 4.4.0)
# fansi                  1.0.6   2023-12-08 [2] CRAN (R 4.4.0)
# farver                 2.1.2   2024-05-13 [2] CRAN (R 4.4.0)
# generics               0.1.3   2022-07-05 [2] CRAN (R 4.4.0)
# GenomeInfoDb         * 1.40.1  2024-05-24 [2] Bioconductor 3.19 (R 4.4.0)
# GenomeInfoDbData       1.2.12  2024-05-23 [2] Bioconductor
# GenomicRanges        * 1.56.1  2024-06-12 [2] Bioconductor 3.19 (R 4.4.0)
# ggplot2              * 3.5.1   2024-04-23 [2] CRAN (R 4.4.0)
# glue                   1.7.0   2024-01-09 [2] CRAN (R 4.4.0)
# gtable                 0.3.5   2024-04-22 [2] CRAN (R 4.4.0)
# here                 * 1.0.1   2020-12-13 [2] CRAN (R 4.4.0)
# httr                   1.4.7   2023-08-15 [2] CRAN (R 4.4.0)
# IRanges              * 2.38.1  2024-07-03 [2] Bioconductor 3.19 (R 4.4.0)
# jsonlite               1.8.8   2023-12-04 [2] CRAN (R 4.4.0)
# labeling               0.4.3   2023-08-29 [2] CRAN (R 4.4.0)
# lattice                0.22-6  2024-03-20 [3] CRAN (R 4.4.0)
# lifecycle              1.0.4   2023-11-07 [2] CRAN (R 4.4.0)
# magrittr               2.0.3   2022-03-30 [2] CRAN (R 4.4.0)
# Matrix                 1.7-0   2024-04-26 [3] CRAN (R 4.4.0)
# MatrixGenerics       * 1.16.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# matrixStats          * 1.3.0   2024-04-11 [2] CRAN (R 4.4.0)
# munsell                0.5.1   2024-04-01 [2] CRAN (R 4.4.0)
# pheatmap             * 1.0.12  2019-01-04 [2] CRAN (R 4.4.0)
# pillar                 1.9.0   2023-03-22 [2] CRAN (R 4.4.0)
# pkgconfig              2.0.3   2019-09-22 [2] CRAN (R 4.4.0)
# Polychrome             1.5.1   2022-05-03 [1] CRAN (R 4.4.0)
# R6                     2.5.1   2021-08-19 [2] CRAN (R 4.4.0)
# ragg                   1.3.2   2024-05-15 [2] CRAN (R 4.4.0)
# RColorBrewer           1.1-3   2022-04-03 [2] CRAN (R 4.4.0)
# rlang                  1.1.4   2024-06-04 [2] CRAN (R 4.4.0)
# rprojroot              2.0.4   2023-11-05 [2] CRAN (R 4.4.0)
# S4Arrays               1.4.1   2024-05-20 [2] Bioconductor 3.19 (R 4.4.0)
# S4Vectors            * 0.42.1  2024-07-03 [2] Bioconductor 3.19 (R 4.4.0)
# scales                 1.3.0   2023-11-28 [2] CRAN (R 4.4.0)
# scatterplot3d          0.3-44  2023-05-05 [1] CRAN (R 4.4.0)
# sessioninfo          * 1.2.2   2021-12-06 [2] CRAN (R 4.4.0)
# SingleCellExperiment * 1.26.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# SparseArray            1.4.8   2024-05-24 [2] Bioconductor 3.19 (R 4.4.0)
# SummarizedExperiment * 1.34.0  2024-05-01 [2] Bioconductor 3.19 (R 4.4.0)
# systemfonts            1.1.0   2024-05-15 [2] CRAN (R 4.4.0)
# textshaping            0.4.0   2024-05-24 [2] CRAN (R 4.4.0)
# tibble                 3.2.1   2023-03-20 [2] CRAN (R 4.4.0)
# tidyselect             1.2.1   2024-03-11 [2] CRAN (R 4.4.0)
# UCSC.utils             1.0.0   2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# utf8                   1.2.4   2023-10-22 [2] CRAN (R 4.4.0)
# vctrs                  0.6.5   2023-12-01 [2] CRAN (R 4.4.0)
# withr                  3.0.1   2024-07-31 [2] CRAN (R 4.4.0)
# XVector                0.44.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# zlibbioc               1.50.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# 
# [1] /users/rphillip/R/4.4
# [2] /jhpce/shared/community/core/conda_R/4.4/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.4/R/lib64/R/library
# 
# ────────────────────────────────────────────────────────────────────────────
# 
