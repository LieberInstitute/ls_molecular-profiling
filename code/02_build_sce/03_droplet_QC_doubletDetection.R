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
#       Sample total_drops non_empty knee_lower
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

#Plot mitochondria vs detected
#All samples together
mito_vs_detected <- plotColData(object = sce,
                                y = "subsets_Mito_percent",
                                x = "detected",
                                colour_by = "Sample")
ggsave(filename = here("plots","mito_vs_detected_bySample.png"),plot = mito_vs_detected)

#Now samples separately.
for(i in unique(sce$Sample)){
    x <- plotColData(object = sce[,sce$Sample == i],
                     y = "subsets_Mito_percent",
                     x = "detected",
                     colour_by = "Sample")
    ggsave(filename = here("plots",paste0("mito_vs_detected_",i,"_only.png")),
           plot = x)
}

#########################################
############### High mito ###############
#########################################
##Comment in https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/51d15ef9f5f2c4c53f55e22e3fe467de1a724668/10x_all-FACS-n10_2021rev_step01_processing-QC_MNT.R#L4
#suggests that MAD approach may unneccisarily throw out cells from samples in which the mito percentage  distribution is centered around 0
#Will try several MAD values + a numeric cutoff of 5% to see how many cells are getting thrown out 
sce$high_mito_1 <- isOutlier(sce$subsets_Mito_percent, nmads = 1, type = "higher", batch = sce$Sample)
table(sce$Sample,sce$high_mito_1)
#           FALSE TRUE
# 1c_LS_SCP  3380 1172
# 2c_LS_SCP  2656 1044
# 3c_LS_SCP  2263  952
sce$high_mito_2 <- isOutlier(sce$subsets_Mito_percent, nmads = 2, type = "higher", batch = sce$Sample)
table(sce$Sample,sce$high_mito_2)
#           FALSE TRUE
# 1c_LS_SCP  3870  682
# 2c_LS_SCP  3221  479
# 3c_LS_SCP  2745  470
sce$high_mito_3 <- isOutlier(sce$subsets_Mito_percent, nmads = 3, type = "higher", batch = sce$Sample)
table(sce$Sample,sce$high_mito_3)
#           FALSE TRUE
# 1c_LS_SCP  4114  438
# 2c_LS_SCP  3469  231
# 3c_LS_SCP  2997  218

for(i in c("high_mito_1","high_mito_2","high_mito_3")){
    x <- plotColData(sce,x = "Sample",y= "subsets_Mito_percent",colour_by = i) +
        ggtitle(paste0("Mito Percent\n",i)) +
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_hline(yintercept = 5,lty = 2) #Line at 5% which is a cutoff that is widely used. 
    ggsave(x,file=here("plots",paste0(i,"_violin.png")))
}

####Numeric cutoff
table(sce$Sample,sce$subsets_Mito_percent > 5.0)
#           FALSE TRUE
# 1c_LS_SCP  4548    4
# 2c_LS_SCP  3411  289
# 3c_LS_SCP  3105  110
sce$high_mito_numeric <- ifelse(sce$subsets_Mito_percent > 5.0,
                                TRUE,
                                FALSE)
table(sce$high_mito_numeric)
# FALSE  TRUE 
# 11064   403 

#Check to see if the cells being dropped have large mito percentages.
mito_violin <- plotColData(sce, x = "Sample", y = "subsets_Mito_percent", colour_by = "high_mito_numeric") +
    ggtitle("Mito Precent") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = 5,lty = 2) +
    annotate(geom="text",label = "5% Mito",x = 0.75,y=6)

ggsave(mito_violin,file=here("plots","mito_percentage_violin_numericCutoff.png"))


#########################################
########## Low Library Size #############
#########################################
# ## low library size
#Plot library size per sample to look at distributions. 
lib_size_violin <- plotColData(sce, x = "Sample", y = "sum",colour_by = "Sample") +
    scale_y_log10() +
    ggtitle("Total UMIs")

ggsave(lib_size_violin,file=here("plots","lib_size_violin.png"))

#Sample 1 has a unimodal distribution while samples 2 and 3 have bimodal distributions. 
#Going to check calculate low library size using the MAD approach with differing numbers of MAD 1-3. 
##1 
sce$low_lib_1 <- isOutlier(sce$sum, log = TRUE, type = "lower", batch = sce$Sample,nmads = 1)
table(sce$Sample,sce$low_lib_1)
#            FALSE TRUE
# 1c_LS_SCP  3181 1371
# 2c_LS_SCP  2849  851
# 3c_LS_SCP  2457  758

##2
sce$low_lib_2 <- isOutlier(sce$sum, log = TRUE, type = "lower", batch = sce$Sample,nmads = 2)
table(sce$Sample,sce$low_lib_2)
#           FALSE TRUE
# 1c_LS_SCP  3858  694
# 2c_LS_SCP  3412  288
# 3c_LS_SCP  2646  569

##3
sce$low_lib_3 <- isOutlier(sce$sum, log = TRUE, type = "lower", batch = sce$Sample,nmads = 3)
table(sce$Sample,sce$low_lib_3)
#           FALSE TRUE
# 1c_LS_SCP  4336  216
# 2c_LS_SCP  3700    0
# 3c_LS_SCP  3152   63

#Plot each 
for(i in c(1:3)){
    sum_violon <- plotColData(object = sce,y = "sum",x = "Sample",colour_by = paste0("low_lib_",i)) +
        scale_y_log10() +
        ggtitle(paste0("Total UMIs"))
    ggsave(here("plots",paste0("lib_size_violin_nmad_",i,".png"))) 
}


#########################################
####### Low Detected Features ###########
#########################################
#Plot number of detected features per sample to look at distributions. 
detected_features_violin <- plotColData(sce, x = "Sample", y = "detected",colour_by = "Sample") +
    scale_y_log10()+
    ggtitle("Detected Features")

ggsave(detected_features_violin,file=here("plots","Detected_Features_violin.png"))

#Similar bimodal distributions in samples 2+3 and unimodal in sample 1. 
#Similar to library size, will investigate different MAD values for number of genes 
##1 
sce$low_genes_1 <- isOutlier(sce$detected, log = TRUE, type = "lower", batch = sce$Sample,nmads = 1)
table(sce$Sample,sce$low_genes_1)
#           FALSE TRUE
# 1c_LS_SCP  3084 1468
# 2c_LS_SCP  2818  882
# 3c_LS_SCP  2439  776

##2
sce$low_genes_2 <- isOutlier(sce$detected, log = TRUE, type = "lower", batch = sce$Sample,nmads = 2)
table(sce$Sample,sce$low_genes_2)
#           FALSE TRUE
# 1c_LS_SCP  3534 1018
# 2c_LS_SCP  3183  517
# 3c_LS_SCP  2608  607

##3
sce$low_genes_3 <- isOutlier(sce$detected, log = TRUE, type = "lower", batch = sce$Sample,nmads = 3)
table(sce$Sample,sce$low_genes_3)
#           FALSE TRUE
# 1c_LS_SCP  3961  591
# 2c_LS_SCP  3700    0
# 3c_LS_SCP  2887  328

#Plot each 
for(i in c(1:3)){
    detected_violin <- plotColData(object = sce,y = "detected",x = "Sample",colour_by = paste0("low_genes_",i)) +
        scale_y_log10() +
        geom_hline(yintercept = 500) +
        ggtitle(paste0("Total Detected features"))
    ggsave(here("plots",paste0("low_genes_violin_nmad_",i,".png"))) 
}

#Check what the distribution looks like when not splitting by sample
x <- plotColData(object = sce,y = "detected")
ggsave(here("plots","detected_features.png"),plot = x)

######################################################################
#Another option is to identify outliers in high-dimensional space based 
#on the QC metrics for each cell" (OSCA ch. 1.3.3)
df <- colData(sce)[,c("sum","detected","subsets_Mito_percent")]
df$sum <- log10(df$sum)
df$detected <- log10(df$detected)

library(robustbase)
outlying <- adjOutlyingness(df, only.outlyingness = TRUE)
multi.outlier <- isOutlier(outlying, type = "higher")
table(multi.outlier)
# multi.outlier
# FALSE  TRUE 
# 10349  1118 

#With this method, interpreting why a cell is dropped is difficult. 
#Plot genes detected, library size, and mito percentage and color by it being defined as an outlier
all(rownames(colData(sce)) == names(multi.outlier))
#[1] TRUE
#Everything is in the same order here, so can just add right to the object. 
sce$robustbase_outlier <- multi.outlier

#Now plot
#Genes
genes_robustbase <- plotColData(object = sce,y = "detected",x = "Sample",colour_by = "robustbase_outlier") +
    scale_y_log10() +
    ggtitle("Total Detected Features")
ggsave(filename = here("plots","genes_violin_robustbaseoutliers.png"))
#library size
libsize_robustbase <- plotColData(object = sce,y = "sum",x = "Sample",colour_by = "robustbase_outlier") +
    scale_y_log10() +
    ggtitle("Total UMIs")
ggsave(filename = here("plots","libsize_violin_robustbaseoutliers.png"))
#mitochondrial genes
mito_robustbase <- plotColData(object = sce,y = "subsets_Mito_percent",x = "Sample",colour_by = "robustbase_outlier") +
    ggtitle("Reads mapping to mitochondrial genome")
ggsave(filename = here("plots","mito_violin_robustbaseoutliers.png"))

#This approach keys in on the mitochondrial percentage and not much else. Still keepign cells with low library 
#and low number of genes. 

############################################
######### Annotate nuclei to drop ##########
############################################
#The MAd approach for mitochondria unnecessarily punishes sample 1 because the distribution is centered 
#around 0. To prove: 
summary(sce[,sce$Sample == "1c_LS_SCP"]$subsets_Mito_percent)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.06194 0.12597 0.21432 0.24457 6.66667
#Using the same MAD values for each sample is also problematic for library size and number of genes. 
#Sample 1 is unimodal and using MAD approach for number of genes throws out cells with >1000 genes, which
#I believe are most likely high quality cells. However, samples 2 and 3 are bimodal. This is a particular
#issue for sample 2 because MAD=3 for library size removes 0 cells and includes those with <1000 reads +
#<500 genes. With this information, I am going to use sample specific MAD values for library size and a numeric cutoff
#of mitochondrial reads of 5%.

#Using sample specific cutoffs. 
####Sample 1
sample_1_drop <- sce[,sce$Sample == "1c_LS_SCP"]$low_lib_3 | sce[,sce$Sample == "1c_LS_SCP"]$high_mito_numeric
table(sample_1_drop)
# FALSE  TRUE 
# 4332   220 
#Get cell IDs that need to be dropped
sample_1_drop_names <- rownames(colData(sce[,sce$Sample == "1c_LS_SCP"])[which(sample_1_drop),])

####Sample 2
sample_2_drop <- sce[,sce$Sample == "2c_LS_SCP"]$low_lib_2 | sce[,sce$Sample == "2c_LS_SCP"]$low_genes_2 | sce[,sce$Sample == "2c_LS_SCP"]$high_mito_numeric
table(sample_2_drop)
# sample_2_drop
# FALSE  TRUE 
# 2924   776 
table(sce[,sce$Sample == "2c_LS_SCP"]$low_lib_2, sce[,sce$Sample == "2c_LS_SCP"]$low_genes_2)
#        FALSE TRUE
# FALSE  3183  229
# TRUE      0  288
#All that have low library size are also low genes.  
#Some that are low genes, do pass library size. 

#Get cell IDs that need to be dropped
sample_2_drop_names <- rownames(colData(sce[,sce$Sample == "2c_LS_SCP"])[which(sample_2_drop),])

####Sample 3
sample_3_drop <- sce[,sce$Sample == "3c_LS_SCP"]$low_lib_2 | sce[,sce$Sample == "3c_LS_SCP"]$low_genes_3 | sce[,sce$Sample == "3c_LS_SCP"]$high_mito_numeric
table(sample_3_drop)
# sample_3_drop
# FALSE  TRUE 
# 2551   664 

table(sce[,sce$Sample == "3c_LS_SCP"]$low_lib_2, sce[,sce$Sample == "3c_LS_SCP"]$low_genes_3)
#       FALSE TRUE
# FALSE  2646    0
# TRUE    241  328
#Some that have low library size, do pass gene cutoffs. 
#However, everything that has low genes is also low libary size. 

#Get cell IDs that need to be dropped
sample_3_drop_names <- rownames(colData(sce[,sce$Sample == "3c_LS_SCP"])[which(sample_3_drop),])

#Concatenate all
cells_to_drop <- c(sample_1_drop_names,sample_2_drop_names,sample_3_drop_names)

#add information to the object
sce$discard_sample_specific <- ifelse(sce$unique_rowname %in% cells_to_drop,
                                      TRUE,
                                      FALSE)

table(sce$discard_sample_specific)
# FALSE  TRUE 
# 9807  1660 

# table(sce$discard_auto)
# FALSE  TRUE 
# 9883  1584 

#table(sce$discard_numeric)
# FALSE  TRUE 
# 10000  1467 

100 * sum(sce$discard_sample_specific) / ncol(sce)
#[1] 14.47632

(qc_t <- addmargins(table(sce$Sample, sce$discard_sample_specific)))
#           FALSE  TRUE   Sum
# 1c_LS_SCP  4332   220  4552
# 2c_LS_SCP  2924   776  3700
# 3c_LS_SCP  2551   664  3215
# Sum        9807  1660 11467

round(100 * sweep(qc_t, 1, qc_t[, 3], "/"), 1)
#           FALSE  TRUE   Sum
# 1c_LS_SCP  95.2   4.8 100.0
# 2c_LS_SCP  79.0  21.0 100.0
# 3c_LS_SCP  79.3  20.7 100.0
# Sum        85.5  14.5 100.0

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
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
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
# Sample    median   q95  drop drop_precent
# <chr>      <dbl> <dbl> <int>        <dbl>
# 1 1c_LS_SCP  0.282 0.992    31        0.681
# 2 2c_LS_SCP  0.496 1.37     18        0.486
# 3 3c_LS_SCP  0.521 1.20     41        1.28 

table(sce$discard_sample_specific, sce$doubletScore >= 5)
#       FALSE TRUE
# FALSE  9726   81
# TRUE   1651    9

#Save object
save(sce,file = here("processed-data","sce_emptyDrops_removed_withQC.rda"))

sce <- sce[,!sce$discard_sample_specific]
dim(sce)
# [1] 36601  9807

## save QCed and cleaned object. 
save(sce,file=here("processed-data","sce_clean.rda"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# [1] "Reproducibility information:"
# [1] "2023-10-20 14:48:40 EDT"
# user   system  elapsed 
# 548.662   22.111 6976.418 
# ─ Session info ───────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 Patched (2023-07-19 r84711)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2023-10-20
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [2] CRAN (R 4.3.1)
# beachmat               2.16.0    2023-04-25 [2] Bioconductor
# beeswarm               0.4.0     2021-06-01 [2] CRAN (R 4.3.1)
# Biobase              * 2.60.0    2023-04-25 [2] Bioconductor
# BiocGenerics         * 0.46.0    2023-04-25 [2] Bioconductor
# BiocIO                 1.10.0    2023-04-25 [2] Bioconductor
# BiocNeighbors          1.18.0    2023-04-25 [2] Bioconductor
# BiocParallel           1.34.2    2023-05-22 [2] Bioconductor
# BiocSingular           1.16.0    2023-04-25 [2] Bioconductor
# Biostrings             2.68.1    2023-05-16 [2] Bioconductor
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.1)
# bluster                1.10.0    2023-04-25 [2] Bioconductor
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.1)
# cluster                2.1.4     2022-08-22 [3] CRAN (R 4.3.1)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.3.1)
# colorout             * 1.2-2     2023-09-22 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.3.1)
# cowplot                1.1.1     2020-12-30 [1] CRAN (R 4.3.1)
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.1)
# data.table             1.14.8    2023-02-17 [2] CRAN (R 4.3.1)
# DelayedArray           0.26.7    2023-07-28 [2] Bioconductor
# DelayedMatrixStats     1.22.6    2023-08-28 [2] Bioconductor
# DEoptimR               1.1-2     2023-08-28 [2] CRAN (R 4.3.1)
# dplyr                * 1.1.3     2023-09-03 [2] CRAN (R 4.3.1)
# dqrng                  0.3.1     2023-08-30 [2] CRAN (R 4.3.1)
# DropletUtils         * 1.20.0    2023-04-25 [2] Bioconductor
# edgeR                  3.42.4    2023-05-31 [2] Bioconductor
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.3.1)
# farver                 2.1.1     2022-07-06 [2] CRAN (R 4.3.1)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.3.1)
# GenomeInfoDb         * 1.36.3    2023-09-07 [2] Bioconductor
# GenomeInfoDbData       1.2.10    2023-07-20 [2] Bioconductor
# GenomicAlignments      1.36.0    2023-04-25 [2] Bioconductor
# GenomicRanges        * 1.52.0    2023-04-25 [2] Bioconductor
# ggbeeswarm             0.7.2     2023-04-29 [2] CRAN (R 4.3.1)
# ggplot2              * 3.4.3     2023-08-14 [2] CRAN (R 4.3.1)
# ggrepel                0.9.3     2023-02-03 [2] CRAN (R 4.3.1)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.3.1)
# gridExtra              2.3       2017-09-09 [2] CRAN (R 4.3.1)
# gtable                 0.3.4     2023-08-21 [2] CRAN (R 4.3.1)
# HDF5Array              1.28.1    2023-05-01 [2] Bioconductor
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.3.1)
# igraph                 1.5.1     2023-08-10 [2] CRAN (R 4.3.1)
# IRanges              * 2.34.1    2023-06-22 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.3.1)
# jsonlite               1.8.7     2023-06-29 [2] CRAN (R 4.3.1)
# labeling               0.4.3     2023-08-29 [2] CRAN (R 4.3.1)
# lattice                0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.3.1)
# limma                  3.56.2    2023-06-04 [2] Bioconductor
# locfit                 1.5-9.8   2023-06-11 [2] CRAN (R 4.3.1)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.3.1)
# MASS                   7.3-60    2023-05-04 [3] CRAN (R 4.3.1)
# Matrix                 1.6-1.1   2023-09-18 [3] CRAN (R 4.3.1)
# MatrixGenerics       * 1.12.3    2023-07-30 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.3.1)
# metapod                1.8.0     2023-04-25 [2] Bioconductor
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.3.1)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.3.1)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.3.1)
# purrr                * 1.0.2     2023-08-10 [2] CRAN (R 4.3.1)
# R.methodsS3            1.8.2     2022-06-13 [2] CRAN (R 4.3.1)
# R.oo                   1.25.0    2022-06-12 [2] CRAN (R 4.3.1)
# R.utils                2.12.2    2022-11-11 [2] CRAN (R 4.3.1)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.1)
# rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 4.3.1)
# ragg                   1.2.5     2023-01-12 [2] CRAN (R 4.3.1)
# RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.3.1)
# Rcpp                   1.0.11    2023-07-06 [2] CRAN (R 4.3.1)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.1)
# restfulr               0.0.15    2022-06-16 [2] CRAN (R 4.3.1)
# rhdf5                  2.44.0    2023-04-25 [2] Bioconductor
# rhdf5filters           1.12.1    2023-04-30 [2] Bioconductor
# Rhdf5lib               1.22.1    2023-09-10 [2] Bioconductor
# rjson                  0.2.21    2022-01-09 [2] CRAN (R 4.3.1)
# rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.1)
# robustbase           * 0.99-0    2023-06-16 [2] CRAN (R 4.3.1)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
# Rsamtools              2.16.0    2023-04-25 [2] Bioconductor
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.3.1)
# rtracklayer            1.60.1    2023-08-15 [2] Bioconductor
# S4Arrays               1.0.6     2023-08-30 [2] Bioconductor
# S4Vectors            * 0.38.1    2023-05-02 [2] Bioconductor
# ScaledMatrix           1.8.1     2023-05-03 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.3.1)
# scater               * 1.28.0    2023-04-25 [2] Bioconductor
# scDblFinder          * 1.14.0    2023-04-25 [1] Bioconductor
# scran                * 1.28.2    2023-07-23 [2] Bioconductor
# scuttle              * 1.10.2    2023-08-03 [2] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.1)
# SingleCellExperiment * 1.22.0    2023-04-25 [2] Bioconductor
# sparseMatrixStats      1.12.2    2023-07-02 [2] Bioconductor
# statmod                1.5.0     2023-01-06 [2] CRAN (R 4.3.1)
# SummarizedExperiment * 1.30.2    2023-06-06 [2] Bioconductor
# systemfonts            1.0.4     2022-02-11 [2] CRAN (R 4.3.1)
# textshaping            0.3.6     2021-10-13 [2] CRAN (R 4.3.1)
# tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.3.1)
# tidyr                * 1.3.0     2023-01-24 [2] CRAN (R 4.3.1)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.3.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.3.1)
# vctrs                  0.6.3     2023-06-14 [2] CRAN (R 4.3.1)
# vipor                  0.4.5     2017-03-22 [2] CRAN (R 4.3.1)
# viridis                0.6.4     2023-07-22 [2] CRAN (R 4.3.1)
# viridisLite            0.4.2     2023-05-02 [2] CRAN (R 4.3.1)
# withr                  2.5.0     2022-03-03 [2] CRAN (R 4.3.1)
