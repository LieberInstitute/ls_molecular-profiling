#Goal: compile droplet scores, calculate QC metrics, and detect doublets. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/
#code modified from https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/03_droplet_qc.R

library(SingleCellExperiment)
library(DropletUtils)
library(scuttle)
library(ggplot2)
library(scater)
library(here)
library(purrr)
library(dplyr)
library(tidyr)
library(rafalib)
library(scran)
library(scDblFinder)
library(sessioninfo)

#Where are the droplet QC files? 
droplet_paths <- list.files(path = here("processed-data","02_build_sce","droplet_scores"),
                            full.names = TRUE)

names(droplet_paths) <- gsub(x = basename(droplet_paths),
                             pattern = "_droplet_scores.Rdata",
                             replacement = "")

#Read in the droplet scores
e.out <- lapply(droplet_paths, function(x) get(load(x)))

#To make sure we aren't throwing out any cells check if Limited=TRUE and SIG==FALSE
#If both are true, then we could be throwing out non-empty droplets. 
lapply(e.out,function(x){
    table(x$Limited == TRUE & x$FDR>0.001)
})
# $`1c_LS_SCP`
# 
# FALSE 
# 15259 
# 
# $`2c_LS_SCP`
# 
# FALSE 
# 7157 
# 
# $`3c_LS_SCP`
# 
# FALSE 
# 6923 

#Another way to look at this
map(e.out, ~ addmargins(table(Signif = .x$FDR <= 0.001, Limited = .x$Limited)))
# $`1c_LS_SCP`
# Limited
# Signif  FALSE  TRUE   Sum
# FALSE 10707     0 10707
# TRUE    151  4401  4552
# Sum   10858  4401 15259
# 
# $`2c_LS_SCP`
# Limited
# Signif  FALSE TRUE  Sum
# FALSE  3457    0 3457
# TRUE    268 3432 3700
# Sum    3725 3432 7157
# 
# $`3c_LS_SCP`
# Limited
# Signif  FALSE TRUE  Sum
# FALSE  3708    0 3708
# TRUE    320 2895 3215
# Sum    4028 2895 6923

#Not losing any droplets due to number of iterations. Can move on to QC plots and metrics. 

#Pull knee lower values
std_out <- readLines(here("code","02_build_sce","emptyDrops_dropletQC.err"))
knee_lowers <- as.numeric(lapply(strsplit(std_out[grep("knee_lower",std_out)],split="="),"[",2))
names(knee_lowers) <- names(e.out)
knee_lowers
# 1c_LS_SCP 2c_LS_SCP 3c_LS_SCP 
# 307       207       207

#Create droplet summary table
droplet_summary <- stack(map_int(e.out,nrow)) %>% 
    rename(total_drops=values) %>% 
    left_join(stack(map_int(e.out, ~ sum(.x$FDR < 0.001, na.rm = TRUE)))) %>%
    rename(non_empty=values) %>%
    left_join(stack(knee_lowers)) %>%
    rename(Sample=ind) %>%
    select(Sample,total_drops,non_empty,knee_lower=values)
droplet_summary
# Sample total_drops non_empty knee_lower
# 1 1c_LS_SCP     2093387      4552        307
# 2 2c_LS_SCP     1392235      3700        207
# 3 3c_LS_SCP     1423592      3215        207

#Write out file. 
write.csv(x = droplet_summary,
          file = here("processed-data","02_build_sce","droplet_summary.csv"),
          row.names = FALSE,
          quote = FALSE)
    

#Make a barplot summarizing the number of empty and non-empty droplets. 
droplet_barplot <- droplet_summary %>%
    mutate(empty = total_drops - non_empty) %>%
    select(-total_drops) %>%
    select(-knee_lower) %>%
    pivot_longer(!Sample,names_to = "drop_type",values_to = "number_drops") %>%
    ggplot(aes(x = Sample,y=number_drops,fill = drop_type)) +
    geom_col() +
    scale_y_continuous(trans = "log10") +
    labs(x = "Sample",
         y = "Number of Droplets",
         fill = "Droplet Status")

ggsave(plot = droplet_barplot,filename = here("plots","droplet_barplot_per_sample.png"))

#Load in the sce object
load(here("processed-data","sce_raw.rda"),verbose = TRUE)

dim(sce)
# [1]   36601 4909214

#### Eliminate empty droplets ####
e.out.all <- do.call("rbind", e.out)[colnames(sce), ]
sce <- sce[, which(e.out.all$FDR <= 0.001)]

dim(sce)
# [1] 36601 11467
#11467 droplets containing cells a this point. 

#Save object
save(sce,file = here("processed-data","sce_emptyDrops_removed.rda"))

####Begin QC
sce <- scuttle::addPerCellQC(sce,subsets = list(Mito=which(seqnames(sce) == "chrM")))

## High mito
##Comment in https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/51d15ef9f5f2c4c53f55e22e3fe467de1a724668/10x_all-FACS-n10_2021rev_step01_processing-QC_MNT.R#L4
#suggests that MAD approach may unneccisarily throw out cells from samples in which the mito percentage  distribution is centered around 0
sce$high_mito <- isOutlier(sce$subsets_Mito_percent, nmads = 3, type = "higher", batch = sce$Sample)
table(sce$high_mito)
# FALSE  TRUE 
# 10580   887 
table(sce$high_mito, sce$Sample)
#       1c_LS_SCP 2c_LS_SCP 3c_LS_SCP
# FALSE      4114      3469      2997
# TRUE        438       231       218

#Check to see if the cells being dropped have large mito percentages. 
mito_violin <- plotColData(sce, x = "Sample", y = "subsets_Mito_percent", colour_by = "high_mito") +
    ggtitle("Mito Precent") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = 5,lty = 2) +
    annotate(geom="text",label = "5% Mito",x = 0.75,y=6)

ggsave(mito_violin,file=here("plots","mito_percentage_violin.png"))

#I think the MAD approach is appropriate here. 
## low library size
sce$low_lib <- isOutlier(sce$sum, log = TRUE, type = "lower", batch = sce$Sample)
table(sce$low_lib)
# FALSE  TRUE 
# 11188   279 

lib_size_violin <- plotColData(sce, x = "Sample", y = "sum", colour_by = "low_lib") +
    scale_y_log10() +
    ggtitle("Total UMIs")

ggsave(lib_size_violin,file=here("plots","lib_size_violin.png"))

## low detected features
# sce$qc.detected
sce$low_genes <- isOutlier(sce$detected, log = TRUE, type = "lower", batch = sce$Sample)
table(sce$low_genes)
# FALSE  TRUE 
# 10548   919 

low_genes_violin <- plotColData(sce, x = "Sample", y = "detected", colour_by = "low_genes") +
    scale_y_log10() +
    ggtitle("Total Detected Features")

ggsave(low_genes_violin,file=here("plots","low_genes_violin.png"))

table(sce$low_lib, sce$low_genes)
#       FALSE  TRUE
# FALSE 10548   640
# TRUE      0   279
#Low library size also have low genes. 

## Annotate nuclei to drop
sce$discard_auto <- sce$high_mito | sce$low_lib | sce$low_genes

table(sce$discard_auto)
# FALSE  TRUE 
# 9883  1584 

## discard 13% of nuc
100 * sum(sce$discard_auto) / ncol(sce)
# [1] 13.81355

(qc_t <- addmargins(table(sce$Sample, sce$discard_auto)))
#            FALSE  TRUE   Sum
# 1c_LS_SCP  3717   835  4552
# 2c_LS_SCP  3469   231  3700
# 3c_LS_SCP  2697   518  3215
# Sum        9883  1584 11467

round(100 * sweep(qc_t, 1, qc_t[, 3], "/"), 1)
# FALSE  TRUE   Sum
# 1c_LS_SCP  81.7  18.3 100.0
# 2c_LS_SCP  93.8   6.2 100.0
# 3c_LS_SCP  83.9  16.1 100.0
# Sum        86.2  13.8 100.0

#### Doublet detection ####
## To speed up, run on sample-level top-HVGs - just take top 1000
set.seed(1234)

colData(sce)$doubletScore <- NA

for (i in splitit(sce$Sample)) {
    sce_temp <- sce[, i]
    ## To speed up, run on sample-level top-HVGs - just take top 1000
    normd <- logNormCounts(sce_temp)
    geneVar <- modelGeneVar(normd)
    topHVGs <- getTopHVGs(geneVar, n = 1000)
    
    dbl_dens <- computeDoubletDensity(normd, subset.row = topHVGs)
    colData(sce)$doubletScore[i] <- dbl_dens
}

summary(sce$doubletScore)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.1730  0.3996  0.5438  0.6660 17.7661 


## Visualize doublet scores ##
dbl_df <- colData(sce) %>%
    as.data.frame() %>%
    select(Sample, doubletScore)

dbl_box_plot <- dbl_df %>%
    ggplot(aes(x = Sample, y = doubletScore, fill = Sample)) +
    geom_boxplot() +
    labs(x = "Sample") +
    geom_hline(yintercept = 5, color = "red", linetype = "dashed") +
    coord_flip() +
    theme_bw()

ggsave(dbl_box_plot, filename = here("plots", "doublet_scores_boxplot.png"))

dbl_density_plot <- dbl_df %>%
    ggplot(aes(x = doubletScore,fill = Sample)) +
    geom_density() +
    labs(x = "doublet score") +
    theme_bw()

ggsave(dbl_density_plot, filename = here("plots", "doublet_scores_desnity.png"))

dbl_df %>%
    group_by(Sample) %>%
    summarize(
        median = median(doubletScore),
        q95 = quantile(doubletScore, .95),
        drop = sum(doubletScore >= 5),
        drop_precent = 100 * drop / n()
    )

# A tibble: 3 × 5
#     Sample    median   q95  drop drop_precent
#.    <chr>      <dbl> <dbl> <int>        <dbl>
# 1 1c_LS_SCP  0.282 0.992    31        0.681
# 2 2c_LS_SCP  0.496 1.37     18        0.486
# 3 3c_LS_SCP  0.521 1.20     41        1.28 

table(sce$discard_auto, sce$doubletScore >= 5)
#        FALSE TRUE
# FALSE  9806   77
# TRUE   1571   13

#Save object
save(sce,file = here("processed-data","sce_emptyDrops_removed.rda"))

#Remove problematic cells. 
sce <- sce[, !sce$discard_auto]
dim(sce)
#[1] 36601  9883

## save QCed and cleaned object. 
save(sce,file = here("processed-data","sce_clean.rda"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# [1] "Reproducibility information:"
# [1] "2023-10-13 15:32:46 EDT"
# # user   system  elapsed 
# 342.463    6.905 1656.542 
# [1] "202─ Session info ────────────────────────────────────────────
#  setting  value
#  version  R version 4.3.1 Patched (2023-07-19 r84711)
#  os       Rocky Linux 9.2 (Blue Onyx)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2023-10-13
#  pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/pandoc
# 
# ─ Packages ────────────────────────────────────────────────
#  package              * version   date (UTC) lib source
#  abind                  1.4-5     2016-07-21 [2] CRAN (R 4.3.1)
#  beachmat               2.16.0    2023-04-25 [2] Bioconductor
#  beeswarm               0.4.0     2021-06-01 [2] CRAN (R 4.3.1)
#  Biobase              * 2.60.0    2023-04-25 [2] Bioconductor
#  BiocGenerics         * 0.46.0    2023-04-25 [2] Bioconductor
#  BiocIO                 1.10.0    2023-04-25 [2] Bioconductor
#  BiocNeighbors          1.18.0    2023-04-25 [2] Bioconductor
#  BiocParallel           1.34.2    2023-05-22 [2] Bioconductor
#  BiocSingular           1.16.0    2023-04-25 [2] Bioconductor
#  Biostrings             2.68.1    2023-05-16 [2] Bioconductor
#  bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.1)
#  bluster                1.10.0    2023-04-25 [2] Bioconductor
#  cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.1)
#  cluster                2.1.4     2022-08-22 [3] CRAN (R 4.3.1)
#  codetools              0.2-19    2023-02-01 [3] CRAN (R 4.3.1)
#  colorout             * 1.2-2     2023-09-22 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.3.1)
#  cowplot                1.1.1     2020-12-30 [2] CRAN (R 4.3.1)
#  crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.1)
#  data.table             1.14.8    2023-02-17 [2] CRAN (R 4.3.1)
#  DelayedArray           0.26.7    2023-07-28 [2] Bioconductor
#  DelayedMatrixStats     1.22.6    2023-08-28 [2] Bioconductor
#  dplyr                * 1.1.3     2023-09-03 [2] CRAN (R 4.3.1)
#  dqrng                  0.3.1     2023-08-30 [2] CRAN (R 4.3.1)
#  DropletUtils         * 1.20.0    2023-04-25 [2] Bioconductor
#  edgeR                  3.42.4    2023-05-31 [2] Bioconductor
#  fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.3.1)
#  farver                 2.1.1     2022-07-06 [2] CRAN (R 4.3.1)
#  generics               0.1.3     2022-07-05 [2] CRAN (R 4.3.1)
#  GenomeInfoDb         * 1.36.3    2023-09-07 [2] Bioconductor
#  GenomeInfoDbData       1.2.10    2023-07-20 [2] Bioconductor
#  GenomicAlignments      1.36.0    2023-04-25 [2] Bioconductor
#  GenomicRanges        * 1.52.0    2023-04-25 [2] Bioconductor
#  ggbeeswarm             0.7.2     2023-04-29 [2] CRAN (R 4.3.1)
#  ggplot2              * 3.4.3     2023-08-14 [2] CRAN (R 4.3.1)
#  ggrepel                0.9.3     2023-02-03 [2] CRAN (R 4.3.1)
#  glue                   1.6.2     2022-02-24 [2] CRAN (R 4.3.1)
#  gridExtra              2.3       2017-09-09 [2] CRAN (R 4.3.1)
#  gtable                 0.3.4     2023-08-21 [2] CRAN (R 4.3.1)
#  HDF5Array              1.28.1    2023-05-01 [2] Bioconductor
#  here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.3.1)
#  igraph                 1.5.1     2023-08-10 [2] CRAN (R 4.3.1)
#  IRanges              * 2.34.1    2023-06-22 [2] Bioconductor
#  irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.3.1)
#  jsonlite               1.8.7     2023-06-29 [2] CRAN (R 4.3.1)
#  labeling               0.4.3     2023-08-29 [2] CRAN (R 4.3.1)
#  lattice                0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
#  lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.3.1)
#  limma                  3.56.2    2023-06-04 [2] Bioconductor
#  locfit                 1.5-9.8   2023-06-11 [2] CRAN (R 4.3.1)
#  magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.3.1)
#  MASS                   7.3-60    2023-05-04 [3] CRAN (R 4.3.1)
#  Matrix                 1.6-1.1   2023-09-18 [3] CRAN (R 4.3.1)
#  MatrixGenerics       * 1.12.3    2023-07-30 [2] Bioconductor
#  matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.3.1)
#  metapod                1.8.0     2023-04-25 [2] Bioconductor
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 4.3.1)
#  pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.3.1)
#  pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.3.1)
#  purrr                * 1.0.2     2023-08-10 [2] CRAN (R 4.3.1)
#  R.methodsS3            1.8.2     2022-06-13 [2] CRAN (R 4.3.1)
#  R.oo                   1.25.0    2022-06-12 [2] CRAN (R 4.3.1)
#  R.utils                2.12.2    2022-11-11 [2] CRAN (R 4.3.1)
#  R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.1)
#  rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 4.3.1)
#  ragg                   1.2.5     2023-01-12 [2] CRAN (R 4.3.1)
#  RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.3.1)
#  Rcpp                   1.0.11    2023-07-06 [2] CRAN (R 4.3.1)
#  RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.1)
#  restfulr               0.0.15    2022-06-16 [2] CRAN (R 4.3.1)
#  rhdf5                  2.44.0    2023-04-25 [2] Bioconductor
#  rhdf5filters           1.12.1    2023-04-30 [2] Bioconductor
#  Rhdf5lib               1.22.1    2023-09-10 [2] Bioconductor
#  rjson                  0.2.21    2022-01-09 [2] CRAN (R 4.3.1)
#  rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.1)
#  rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
#  Rsamtools              2.16.0    2023-04-25 [2] Bioconductor
#  rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.3.1)
#  rtracklayer            1.60.1    2023-08-15 [2] Bioconductor
#  S4Arrays               1.0.6     2023-08-30 [2] Bioconductor
#  S4Vectors            * 0.38.1    2023-05-02 [2] Bioconductor
#  ScaledMatrix           1.8.1     2023-05-03 [2] Bioconductor
#  scales                 1.2.1     2022-08-20 [2] CRAN (R 4.3.1)
#  scater               * 1.28.0    2023-04-25 [2] Bioconductor
#  scDblFinder          * 1.14.0    2023-04-25 [1] Bioconductor
#  scran                * 1.28.2    2023-07-23 [2] Bioconductor
#  scuttle              * 1.10.2    2023-08-03 [2] Bioconductor
#  sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.1)
#  SingleCellExperiment * 1.22.0    2023-04-25 [2] Bioconductor
#  sparseMatrixStats      1.12.2    2023-07-02 [2] Bioconductor
#  statmod                1.5.0     2023-01-06 [2] CRAN (R 4.3.1)
#  SummarizedExperiment * 1.30.2    2023-06-06 [2] Bioconductor
#  systemfonts            1.0.4     2022-02-11 [2] CRAN (R 4.3.1)
#  textshaping            0.3.6     2021-10-13 [2] CRAN (R 4.3.1)
#  tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.3.1)
#  tidyr                * 1.3.0     2023-01-24 [2] CRAN (R 4.3.1)
#  tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.3.1)
#  utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.3.1)
#  vctrs                  0.6.3     2023-06-14 [2] CRAN (R 4.3.1)
#  vipor                  0.4.5     2017-03-22 [2] CRAN (R 4.3.1)
#  viridis                0.6.4     2023-07-22 [2] CRAN (R 4.3.1)
#  viridisLite            0.4.2     2023-05-02 [2] CRAN (R 4.3.1)
#  withr                  2.5.0     2022-03-03 [2] CRAN (R 4.3.1)
#  xgboost                1.7.5.1   2023-03-30 [2] CRAN (R 4.3.1)
#  XML                    3.99-0.14 2023-03-19 [2] CRAN (R 4.3.1)
#  XVector                0.40.0    2023-04-25 [2] Bioconductor
#  yaml                   2.3.7     2023-01-23 [2] CRAN (R 4.3.1)
#  zlibbioc               1.46.0    2023-04-25 [2] Bioconductor
# 
#  [1] /users/rphillip/R/4.3
#  [2] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/site-library
#  [3] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library
# 
# ───────────────────────────────────────────────────────────
# 
