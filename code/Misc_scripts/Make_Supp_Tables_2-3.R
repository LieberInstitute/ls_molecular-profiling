#Make supplementary tables 2 and 3
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling

#load librarires
library(here)
library(sessioninfo)

#Load in the sce object
load(here("processed-data","02_build_sce","sce_celltype.rda"),verbose = TRUE)
# Loading objects:
#   sce

x <- table(sce$CellType.Final,sce$Sample)
write.csv(x,file = here("processed-data","Supp_table_2.csv"))


####Supplememtary Table 3
#Load pairwise
load(here("processed-data","markers_pairwise_list_k_20_louvain_CellTypeFinal.rda"),verbose = TRUE)
# Loading objects:
#   markers_pairwise

for(i in names(markers_pairwise)){
  markers_pairwise[[i]]$CellType <- i
  markers_pairwise[[i]] <- subset(markers_pairwise[[i]],subset=(FDR <= 0.05))
  markers_pairwise[[i]] <- markers_pairwise[[i]][,c("p.value","FDR","summary.stats",
                                                    "gene_id","gene_name","CellType")]
}

#Combine into a single data.frame
markers_pairwise_collapse <- do.call(what = rbind,markers_pairwise)

rownames(markers_pairwise_collapse) <- NULL


write.csv()

#Load 1vALL 
load(here("processed-data","markers_1vAll_ttest_k_20_louvain_CellType_Final.rda"),verbose = TRUE)
# Loading objects:
#   markers_1vALL_df

dim(markers_1vALL_df)
#[1] 838900      9

#subset for only significant DEGs
markers_1vALL_df_sig <- subset(markers_1vALL_df,subset=(log.FDR <= log(0.05)))

dim(markers_1vALL_df_sig)
#[1] 138846      9

#Load in the LS broad DEGs
load(here("processed-data","human_1vALL_broad.rda"),verbose = TRUE)
# Loading objects:
#   human_1vALL_broad

#Subset for the LS genes only
human_1vALL_broad <- subset(human_1vALL_broad,subset=(cellType.target == "LS" & 
                                                        log.FDR <= log(0.05)))


#Write out an excel sheet that has all of the DEGs
library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

#Create worksheets within the excel workbook 
addWorksheet(wb, "1vALL")
addWorksheet(wb, "Pairwise")
addWorksheet(wb, "LS_vs_All_Other")

#Write out the data
writeData(wb, sheet = "1vALL", markers_1vALL_df_sig)
writeData(wb, sheet = "Pairwise", markers_pairwise_collapse)
writeData(wb, sheet = "LS_vs_All_Other", human_1vALL_broad)

# Save the workbook as an Excel file
saveWorkbook(wb, file = here("processed-data","Supp_Table_3.xlsx"), overwrite = TRUE)



print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# [1] "Reproducibility information:"
# [1] "2024-09-16 15:41:14 EDT"
# user   system  elapsed 
# 106.912    6.057 4030.303 
# ─ Session info ───────────────────────────────────────────────────────────────────────────────
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
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────
# package              * version date (UTC) lib source
# abind                  1.4-5   2016-07-21 [2] CRAN (R 4.4.0)
# Biobase              * 2.64.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# BiocGenerics         * 0.50.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# cli                    3.6.3   2024-06-21 [2] CRAN (R 4.4.0)
# crayon                 1.5.3   2024-06-20 [2] CRAN (R 4.4.0)
# DelayedArray           0.30.1  2024-05-07 [2] Bioconductor 3.19 (R 4.4.0)
# GenomeInfoDb         * 1.40.1  2024-05-24 [2] Bioconductor 3.19 (R 4.4.0)
# GenomeInfoDbData       1.2.12  2024-05-23 [2] Bioconductor
# GenomicRanges        * 1.56.1  2024-06-12 [2] Bioconductor 3.19 (R 4.4.0)
# here                 * 1.0.1   2020-12-13 [2] CRAN (R 4.4.0)
# httr                   1.4.7   2023-08-15 [2] CRAN (R 4.4.0)
# IRanges              * 2.38.1  2024-07-03 [2] Bioconductor 3.19 (R 4.4.0)
# jsonlite               1.8.8   2023-12-04 [2] CRAN (R 4.4.0)
# lattice                0.22-6  2024-03-20 [3] CRAN (R 4.4.0)
# Matrix                 1.7-0   2024-04-26 [3] CRAN (R 4.4.0)
# MatrixGenerics       * 1.16.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# matrixStats          * 1.3.0   2024-04-11 [2] CRAN (R 4.4.0)
# openxlsx             * 4.2.7   2024-08-30 [1] CRAN (R 4.4.0)
# R6                     2.5.1   2021-08-19 [2] CRAN (R 4.4.0)
# Rcpp                   1.0.13  2024-07-17 [2] CRAN (R 4.4.0)
# rprojroot              2.0.4   2023-11-05 [2] CRAN (R 4.4.0)
# S4Arrays               1.4.1   2024-05-20 [2] Bioconductor 3.19 (R 4.4.0)
# S4Vectors            * 0.42.1  2024-07-03 [2] Bioconductor 3.19 (R 4.4.0)
# sessioninfo          * 1.2.2   2021-12-06 [2] CRAN (R 4.4.0)
# SingleCellExperiment * 1.26.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# SparseArray            1.4.8   2024-05-24 [2] Bioconductor 3.19 (R 4.4.0)
# stringi                1.8.4   2024-05-06 [2] CRAN (R 4.4.0)
# SummarizedExperiment * 1.34.0  2024-05-01 [2] Bioconductor 3.19 (R 4.4.0)
# UCSC.utils             1.0.0   2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# XVector                0.44.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# zip                    2.3.1   2024-01-27 [2] CRAN (R 4.4.0)
# zlibbioc               1.50.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# 
# [1] /users/rphillip/R/4.4
# [2] /jhpce/shared/community/core/conda_R/4.4/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.4/R/lib64/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────
# 
# 
# 




