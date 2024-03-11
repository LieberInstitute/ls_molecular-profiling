#Make plots for correlation of random shuffling
#Make distributions for correlation coefficients 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/
library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(here)

#Load output dataframes
load(here("processed-data","randomized_t_stats_human_500.rda"),verbose = TRUE)
# Loading objects:
#     out_h_dataframe

load(here("processed-data","randomized_t_stats_mouse_500.rda"),verbose = TRUE)
# Loading objects:
#     out_m_dataframe

#Load the human and mouse matched SCE objects.
load(here("processed-data","human_mouse_matched_by_JAX.rda"),verbose = TRUE)
# Loading objects:
#     sce_mouse_sub
#     sce_human_sub

#Add JAX.gene ID info
####Human
out_h_dataframe$gene_id <- rownames(out_h_dataframe)
out_h_dataframe <- merge(x  = out_h_dataframe,
                         y  = rowData(sce_human_sub)[,c("gene_id","JAX.geneID")],
                         by = "gene_id")
rownames(out_h_dataframe) <- out_h_dataframe$JAX.geneID
out_h_dataframe <- out_h_dataframe[,2:501]

####Mouse
out_m_dataframe$gene_id <- rownames(out_m_dataframe)
out_m_dataframe <- merge(x  = out_m_dataframe,
                         y  = rowData(sce_mouse_sub)[,c("gene_id","JAX.geneID")],
                         by = "gene_id")
rownames(out_m_dataframe) <- out_m_dataframe$JAX.geneID
out_m_dataframe <- out_m_dataframe[,2:501]

#make sure the mouse dataframe is in the same order as the human dataframe
out_m_dataframe <- out_m_dataframe[rownames(out_h_dataframe),]

#Sanity check
all(rownames(out_m_dataframe) == rownames(out_h_dataframe))
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
# [1] "2024-01-26 11:43:13 EST"
#     user   system  elapsed 
# 835.899   43.211 4722.918 
# ─ Session info ────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 Patched (2023-07-19 r84711)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-01-26
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/pandoc
# 
# ─ Packages ────────────────────────────────────────────────────────────────────────
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
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.3.1)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.3.1)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.1)
# ragg                   1.2.5     2023-01-12 [2] CRAN (R 4.3.1)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.1)
# rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.1)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
# S4Arrays               1.0.6     2023-08-30 [2] Bioconductor
# S4Vectors            * 0.38.1    2023-05-02 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.3.1)
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
# ───────────────────────────────────────────────────────────────────────────────────


