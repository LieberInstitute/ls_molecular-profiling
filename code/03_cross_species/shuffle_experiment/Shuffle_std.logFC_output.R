#Make plots for correlation of random shuffling
#Make distributions for correlation coefficients 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/
library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(here)

#Load output dataframes
load(here("processed-data","randomized_std.logFC_human_500.rda"),verbose = TRUE)
# Loading objects:
#   out_h_dataframe

dim(out_h_dataframe)
#[1] 16574   500

load(here("processed-data","randomized_std.logFC_mouse_500.rda"),verbose = TRUE)
# Loading objects:
#   out_m_dataframe

dim(out_m_dataframe)
# [1] 16574   500

#Load the human and mouse matched SCE objects.
load(here("processed-data","sce_human_sub.rda"),verbose = TRUE)
# Loading objects:
#   sce_human_sub

load(here("processed-data","sce_mouse_sub.rda"),verbose = TRUE)
# Loading objects:
#   sce_mouse_sub

#Add JAX.gene ID info
####Human
out_h_dataframe$gene_id <- rownames(out_h_dataframe)
out_h_dataframe <- merge(x  = out_h_dataframe,
                         y  = rowData(sce_human_sub)[,c("gene_id","JAX.geneID")],
                         by = "gene_id")
dim(out_h_dataframe)
#[1] 16574   502
#2 extra columns from ensemble and JAX gene ID

#Make the rownames the JAX.geneID 
rownames(out_h_dataframe) <- out_h_dataframe$JAX.geneID

#Get rid of the geneID columns and just keep the numeric values. 
out_h_dataframe <- out_h_dataframe[,2:501]

####Mouse
out_m_dataframe$gene_id <- rownames(out_m_dataframe)
out_m_dataframe <- merge(x  = out_m_dataframe,
                         y  = rowData(sce_mouse_sub)[,c("gene_id","JAX.geneID")],
                         by = "gene_id")
dim(out_m_dataframe)
#[1] 16574   502
#2 extra column names are ensembl and jax gene id

#Make the rownames the jax gene id
rownames(out_m_dataframe) <- out_m_dataframe$JAX.geneID

#Subset the dataframe
out_m_dataframe <- out_m_dataframe[,2:501]

#make sure the mouse dataframe is in the same order as the human dataframe
out_m_dataframe <- out_m_dataframe[rownames(out_h_dataframe),]

#Sanity check
identical(rownames(out_m_dataframe),rownames(out_h_dataframe))
#[1] TRUE

#Make a dataframe that contains the correlation coefficients and p-values for every set. 
results_df <- as.data.frame(matrix(nrow = 500,ncol = 3))
colnames(results_df) <- c("correlation_coefficient","p_value","test_statistic")

#Make plots for every paired column
for(i in 1:500){
  print(i)
  #Cbind for correlation
  x <- as.data.frame(cbind(out_m_dataframe[,i],
                           out_h_dataframe[,i]))
  colnames(x) <- c("mouse","human")
  
  #Add statistics to dataframe
  results_df[i,"correlation_coefficient"] <- as.numeric(cor.test(x[,"mouse"],x[,"human"])$estimate)
  results_df[i,"p_value"] <- cor.test(x[,"mouse"],x[,"human"])$p.value
  results_df[i,"test_statistic"] <- as.numeric(cor.test(x[,"mouse"],x[,"human"])$statistic)
  
  #Make plot
  cor_plot <- ggplot(data = x,aes(x = human,y = mouse)) +
    geom_point() +
    ggtitle(paste0("r=",
                   round(cor(x[,"mouse"],x[,"human"]),digits = 3),
                   "\np=",
                   round(cor.test(x[,"mouse"],x[,"human"])$p.val,digits = 3))) + 
    theme_bw() +
    geom_hline(yintercept = 0,lty = 2) +
    geom_vline(xintercept = 0,lty = 2) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(plot = cor_plot,filename = here("plots","Conservation","Shuffled_cells_cor_plots",
                                         paste0(i,"_shuffled_plot.pdf")))
}

#save the results dataframe
save(results_df,file = here("processed-data","shuffle_stats_results_df.rda"))

#correlation coefficient distribution
cor_coef_distr <- ggplot(results_df,aes(x = correlation_coefficient)) +
  geom_histogram(aes(y=after_stat(density)),binwidth = 0.01) +
  geom_density(alpha=.2, fill="#FF6666") +
  labs(x = "Pearson's correlation coefficient",
       y = "Density") +
  theme_bw() +
  ggtitle("Correlation coefficient distribution") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = here("plots","Conservation","Correlation_Coefficient_distribution_random_shuffling.pdf"),
       plot = cor_coef_distr)

#p-value distribution
p_distr <- ggplot(results_df,aes(x = p_value)) +
  geom_histogram(aes(y=after_stat(density)),binwidth = 0.01) +
  geom_density(alpha=.2, fill="#FF6666") +
  labs(x = "p-value",
       y = "Density") +
  theme_bw() +
  ggtitle("p-value distribution") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = here("plots","Conservation","p_value_distribution_random_shuffling.pdf"),
       plot = p_distr)

#test statistic distribution
test_stat_distr <- ggplot(results_df,aes(x = test_statistic)) +
  geom_histogram(aes(y=after_stat(density))) +
  geom_density(alpha=.2, fill="#FF6666") +
  xlim(c(-20,20)) +
  labs(x = "Test statistic",
       y = "Density") +
  theme_bw() +
  ggtitle("Test statistic distribution") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = here("plots","Conservation","Test_statistic_distribution_random_shuffling.pdf"),
       plot = test_stat_distr)

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# [1] "Reproducibility information:"
# [1] "2024-09-08 15:42:22 EDT"
# user   system  elapsed 
# 281.904    5.913 1534.350 
# ─ Session info ───────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.4.0 Patched (2024-05-22 r86590)
# os       Rocky Linux 9.4 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-09-08
# pandoc   3.1.13 @ /jhpce/shared/community/core/conda_R/4.4/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────
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
# pillar                 1.9.0   2023-03-22 [2] CRAN (R 4.4.0)
# pkgconfig              2.0.3   2019-09-22 [2] CRAN (R 4.4.0)
# R6                     2.5.1   2021-08-19 [2] CRAN (R 4.4.0)
# ragg                   1.3.2   2024-05-15 [2] CRAN (R 4.4.0)
# rlang                  1.1.4   2024-06-04 [2] CRAN (R 4.4.0)
# rprojroot              2.0.4   2023-11-05 [2] CRAN (R 4.4.0)
# S4Arrays               1.4.1   2024-05-20 [2] Bioconductor 3.19 (R 4.4.0)
# S4Vectors            * 0.42.1  2024-07-03 [2] Bioconductor 3.19 (R 4.4.0)
# scales                 1.3.0   2023-11-28 [2] CRAN (R 4.4.0)
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
# ──────────────────────────────────────────────────────────────────────────────────
