#goal: Color clusters by brain region of origin + glia (LS,MS,Sept,Str,Excit,Glia)
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/
#module load conda_R/4.4

library(SingleCellExperiment)
library(Polychrome)
library(ggplot2)
library(scater)
library(scran)
library(here)

#Load the sce object
load(here("processed-data","02_build_sce","sce_celltype.rda"),verbose = TRUE)
# Loading objects:
#   sce

#What are the cell types again? 
levels(sce$CellType.Final)
# [1] "LS_Inh_A"       "LS_Inh_B"       "LS_Inh_C"       "LS_Inh_D"      
# [5] "LS_Inh_E"       "MS_Inh_A"       "MS_Inh_B"       "MS_Inh_C"      
# [9] "MS_Inh_D"       "Sept_Inh_A"     "Sept_Inh_B"     "Sept_Inh_C"    
# [13] "Str_DRD1_MSN_A" "Str_DRD1_MSN_B" "Str_DRD1_Patch" "Str_DRD2_MSN"  
# [17] "MS_Excit_A"     "Excit_A"        "Excit_B"        "Oligo"         
# [21] "OPC"            "Astrocyte"      "Ependymal"      "Microglia"     
# [25] "Mural"  

#create a dataframe to add to the coldata
Origin_df <- data.frame(CellType = levels(sce$CellType.Final),
                        Origin    = c(rep("LS",5),
                                      rep("MS",4),
                                      rep("Septal",3),
                                      rep("Str",4),
                                      "MS",
                                      rep("TT",2),
                                      rep("Non-Neuronal",6)
                                      )
                        )
Origin_df
#          CellType       Origin
# 1        LS_Inh_A           LS
# 2        LS_Inh_B           LS
# 3        LS_Inh_C           LS
# 4        LS_Inh_D           LS
# 5        LS_Inh_E           LS
# 6        MS_Inh_A           MS
# 7        MS_Inh_B           MS
# 8        MS_Inh_C           MS
# 9        MS_Inh_D           MS
# 10     Sept_Inh_A       Septal
# 11     Sept_Inh_B       Septal
# 12     Sept_Inh_C       Septal
# 13 Str_DRD1_MSN_A          Str
# 14 Str_DRD1_MSN_B          Str
# 15 Str_DRD1_Patch          Str
# 16   Str_DRD2_MSN          Str
# 17     MS_Excit_A           MS
# 18        Excit_A           TT
# 19        Excit_B           TT
# 20          Oligo Non-Neuronal
# 21            OPC Non-Neuronal
# 22      Astrocyte Non-Neuronal
# 23      Ependymal Non-Neuronal
# 24      Microglia Non-Neuronal
# 25          Mural Non-Neuronal

#Create a colData column that is the broad region where the cells originate
sce$Origin <- Origin_df$Origin[match(sce$CellType.Final,Origin_df$CellType)]

table(sce$Origin)
#   LS           MS Non-Neuronal       Septal          Str           TT 
# 1680          784         3837          798         1763          363 

#Create a vector of colors 
#LS,MS,Sept
origin_cols <- c("#1E88E5","#D81B60","#FFC107","#004D40","#000000","#B3B3B3")

#name the colors
names(origin_cols) <- c("LS","MS","Str","Septal","TT","Non-Neuronal")  

#plot tSNE colored by region of origin
x <- plotReducedDim(object = sce,dimred = "tSNE_mnn_50",color_by = "Origin") +
  scale_color_manual(values = origin_cols) +
  theme_void()

ggsave(plot = x,filename = here("plots","Dim_Red","tSNE_Brain_RegionofOrigin.pdf"))


sessionInfo()  
# R version 4.4.0 Patched (2024-05-22 r86590)
# Platform: x86_64-conda-linux-gnu
# Running under: Rocky Linux 9.4 (Blue Onyx)
# 
# Matrix products: default
# BLAS:   /jhpce/shared/community/core/conda_R/4.4/R/lib64/R/lib/libRblas.so 
# LAPACK: /jhpce/shared/community/core/conda_R/4.4/R/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: US/Eastern
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices datasets  utils     methods  
# [8] base     
# 
# other attached packages:
#   [1] here_1.0.1                  scran_1.32.0               
# [3] scater_1.32.1               scuttle_1.14.0             
# [5] ggplot2_3.5.1               Polychrome_1.5.1           
# [7] SingleCellExperiment_1.26.0 SummarizedExperiment_1.34.0
# [9] Biobase_2.64.0              GenomicRanges_1.56.1       
# [11] GenomeInfoDb_1.40.1         IRanges_2.38.1             
# [13] S4Vectors_0.42.1            BiocGenerics_0.50.0        
# [15] MatrixGenerics_1.16.0       matrixStats_1.3.0          
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1          viridisLite_0.4.2        
# [3] dplyr_1.1.4               vipor_0.4.7              
# [5] farver_2.1.2              viridis_0.6.5            
# [7] bluster_1.14.0            rsvd_1.0.5               
# [9] lifecycle_1.0.4           cluster_2.1.6            
# [11] statmod_1.5.0             magrittr_2.0.3           
# [13] compiler_4.4.0            rlang_1.1.4              
# [15] tools_4.4.0               igraph_2.0.3             
# [17] utf8_1.2.4                S4Arrays_1.4.1           
# [19] labeling_0.4.3            dqrng_0.4.1              
# [21] scatterplot3d_0.3-44      DelayedArray_0.30.1      
# [23] abind_1.4-5               BiocParallel_1.38.0      
# [25] withr_3.0.1               grid_4.4.0               
# [27] fansi_1.0.6               beachmat_2.20.0          
# [29] colorspace_2.1-1          edgeR_4.2.1              
# [31] scales_1.3.0              cli_3.6.3                
# [33] crayon_1.5.3              ragg_1.3.2               
# [35] generics_0.1.3            metapod_1.12.0           
# [37] httr_1.4.7                DelayedMatrixStats_1.26.0
# [39] ggbeeswarm_0.7.2          zlibbioc_1.50.0          
# [41] parallel_4.4.0            XVector_0.44.0           
# [43] vctrs_0.6.5               Matrix_1.7-0             
# [45] jsonlite_1.8.8            BiocSingular_1.20.0      
# [47] BiocNeighbors_1.22.0      ggrepel_0.9.5            
# [49] irlba_2.3.5.1             beeswarm_0.4.0           
# [51] systemfonts_1.1.0         locfit_1.5-9.10          
# [53] limma_3.60.4              glue_1.7.0               
# [55] codetools_0.2-20          cowplot_1.1.3            
# [57] gtable_0.3.5              UCSC.utils_1.0.0         
# [59] ScaledMatrix_1.12.0       munsell_0.5.1            
# [61] tibble_3.2.1              pillar_1.9.0             
# [63] GenomeInfoDbData_1.2.12   R6_2.5.1                 
# [65] textshaping_0.4.0         sparseMatrixStats_1.16.0 
# [67] rprojroot_2.0.4           lattice_0.22-6           
# [69] Rcpp_1.0.13               gridExtra_2.3            
# [71] SparseArray_1.4.8         pkgconfig_2.0.3      
#   
#   
#   
