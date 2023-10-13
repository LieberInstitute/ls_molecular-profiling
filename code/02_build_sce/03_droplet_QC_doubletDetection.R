#Goal: compile droplet scores, calculate QC metrics, and detect doublets. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/
#code modified from https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/03_droplet_qc.R

library(SingleCellExperiment)
library(DropletUtils)
library(scuttle)
library(ggplot2)
library(here)
library(purrr)
library(dplyr)
library(tidyr)

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
#All FALSE, so don't need to rerun with increased iterations. 

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

sce
# class: SingleCellExperiment 
# dim: 36601 4909214 
# metadata(1): Samples
# assays(1): counts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(6): source type ... gene_name gene_type
# colnames(4909214): 1_AAACCCAAGAAACCAT-1 1_AAACCCAAGAAACCCA-1 ...
# 3_TTTGTTGTCTTTGCGC-1 3_TTTGTTGTCTTTGGAG-1
# colData names(33): Sample Barcode ... Mean_Reads_per_Cell
# Median_Genes_per_Cell
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):

#### Eliminate empty droplets ####
e.out.all <- do.call("rbind", e.out)[colnames(sce), ]
sce <- sce[, which(e.out.all$FDR <= 0.001)]

sce
# class: SingleCellExperiment 
# dim: 36601 11467 
# metadata(1): Samples
# assays(1): counts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(6): source type ... gene_name gene_type
# colnames(11467): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 3_TTTGTTGAGGCTCCCA-1 3_TTTGTTGTCCCGATCT-1
# colData names(33): Sample Barcode ... Mean_Reads_per_Cell
# Median_Genes_per_Cell
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):

#11467 droplets containing cells a this point. 
#Save object
save()

