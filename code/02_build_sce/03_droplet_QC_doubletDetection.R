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

#Another way to look at this
map(e.out, ~ addmargins(table(Signif = .x$FDR <= 0.001, Limited = .x$Limited,useNA = "ifany")))

#Pull knee lower values
std_out <- readLines(here("code","02_build_sce","logs","02_emptyDrops.log"))
knee_lowers <- as.numeric(lapply(strsplit(std_out[grep("knee_lower",std_out)],split="="),"[",2))
names(knee_lowers) <- names(e.out)
knee_lowers

#Create droplet summary table
droplet_summary <- stack(map_int(e.out,nrow)) %>% 
  rename(total_drops=values) %>% 
  left_join(stack(map_int(e.out, ~ sum(.x$FDR < 0.001, na.rm = TRUE)))) %>%
  rename(non_empty=values) %>%
  left_join(stack(knee_lowers)) %>%
  rename(Sample=ind) %>%
  select(Sample,total_drops,non_empty,knee_lower=values)
droplet_summary

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
load(here("processed-data","02_build_sce","sce_raw_081724.rda"),verbose = TRUE)

sce

dim(sce)

stopifnot(identical(rownames(colData(sce)),colnames(sce)))

#### Eliminate empty droplets ####
e.out.all <- do.call("rbind", e.out)[colnames(sce), ]

#Double check that e.out.all is in same order as sce
stopifnot(identical(rownames(e.out.all),colnames(sce)))

sce <- sce[, which(e.out.all$FDR <= 0.001)]

dim(sce)

#Save object
save(sce,file = here("processed-data","02_build_sce","sce_emptyDrops_removed_081724.rda"))

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
message("High mito 1")
sce$high_mito_1 <- isOutlier(sce$subsets_Mito_percent, nmads = 1, type = "higher", batch = sce$Sample)
table(sce$Sample,sce$high_mito_1)

message("High mito 2")
sce$high_mito_2 <- isOutlier(sce$subsets_Mito_percent, nmads = 2, type = "higher", batch = sce$Sample)
table(sce$Sample,sce$high_mito_2)

message("High mito 3")
sce$high_mito_3 <- isOutlier(sce$subsets_Mito_percent, nmads = 3, type = "higher", batch = sce$Sample)
table(sce$Sample,sce$high_mito_3)


for(i in c("high_mito_1","high_mito_2","high_mito_3")){
  x <- plotColData(sce,x = "Sample",y= "subsets_Mito_percent",colour_by = i) +
    ggtitle(paste0("Mito Percent\n",i)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = 5,lty = 2) #Line at 5% which is a cutoff that is widely used. 
  ggsave(x,file=here("plots",paste0(i,"_violin.png")))
}

####Numeric cutoff
message("Summary of Mito percentage QC")
summary(sce$subsets_Mito_percent)

message("How many nuclei thrown out with percentage > 5?")

sce$high_mito_numeric <- ifelse(sce$subsets_Mito_percent > 5.0,
                                TRUE,
                                FALSE)

#Check to see if the cells being dropped have large mito percentages.
mito_violin <- plotColData(sce, x = "Sample", 
                           y = "subsets_Mito_percent", 
                           colour_by = "high_mito_numeric") +
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
message("Low lib size  1")
sce$low_lib_1 <- isOutlier(sce$sum, log = TRUE, type = "lower", batch = sce$Sample,nmads = 1)
table(sce$Sample,sce$low_lib_1)

##2
message("Low lib size 2")
sce$low_lib_2 <- isOutlier(sce$sum, log = TRUE, type = "lower", batch = sce$Sample,nmads = 2)
table(sce$Sample,sce$low_lib_2)

##3
message("Low lib size  3")
sce$low_lib_3 <- isOutlier(sce$sum, log = TRUE, type = "lower", batch = sce$Sample,nmads = 3)
table(sce$Sample,sce$low_lib_3)

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
message("Low genes  1")
sce$low_genes_1 <- isOutlier(sce$detected, log = TRUE, type = "lower", batch = sce$Sample,nmads = 1)
table(sce$Sample,sce$low_genes_1)

##2
message("Low genes 2")
sce$low_genes_2 <- isOutlier(sce$detected, log = TRUE, type = "lower", batch = sce$Sample,nmads = 2)
table(sce$Sample,sce$low_genes_2)

##3
message("Low genes 3")
sce$low_genes_3 <- isOutlier(sce$detected, log = TRUE, type = "lower", batch = sce$Sample,nmads = 3)
table(sce$Sample,sce$low_genes_3)

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

############################################
######### Annotate nuclei to drop ##########
############################################
#The MAd approach for mitochondria unnecessarily punishes sample 1 because the distribution is centered 
#around 0. To prove: 
summary(sce[,sce$Sample == "1c_LS_SCP"]$subsets_Mito_percent)
#Using the same MAD values for each sample is also problematic for library size and number of genes. 
#Sample 1 is unimodal and using MAD approach for number of genes throws out cells with >1000 genes, which
#I believe are most likely high quality cells. However, samples 2 and 3 are bimodal. This is a particular
#issue for sample 2 because MAD=3 for library size removes 0 cells and includes those with <1000 reads +
#<500 genes. With this information, I am going to use sample specific MAD values for library size and a numeric cutoff
#of mitochondrial reads of 5%.

#Using sample specific cutoffs. 
####Sample 1
message("Sample 1 drop")
sample_1_drop <- sce[,sce$Sample == "1c_LS_SCP"]$low_lib_3 | sce[,sce$Sample == "1c_LS_SCP"]$high_mito_numeric
table(sample_1_drop)
sample_1_drop_names <- rownames(colData(sce[,sce$Sample == "1c_LS_SCP"])[which(sample_1_drop),])

####Sample 2
message("Sample 2 drop")
sample_2_drop <- sce[,sce$Sample == "2c_LS_SCP"]$low_lib_2 | sce[,sce$Sample == "2c_LS_SCP"]$low_genes_2 | sce[,sce$Sample == "2c_LS_SCP"]$high_mito_numeric
table(sample_2_drop)
table(sce[,sce$Sample == "2c_LS_SCP"]$low_lib_2, sce[,sce$Sample == "2c_LS_SCP"]$low_genes_2)

#Get cell IDs that need to be dropped
sample_2_drop_names <- rownames(colData(sce[,sce$Sample == "2c_LS_SCP"])[which(sample_2_drop),])

####Sample 3
message("Sample 3 drop")
sample_3_drop <- sce[,sce$Sample == "3c_LS_SCP"]$low_lib_2 | sce[,sce$Sample == "3c_LS_SCP"]$low_genes_3 | sce[,sce$Sample == "3c_LS_SCP"]$high_mito_numeric
table(sample_3_drop)

table(sce[,sce$Sample == "3c_LS_SCP"]$low_lib_2, sce[,sce$Sample == "3c_LS_SCP"]$low_genes_3)

#Get cell IDs that need to be dropped
sample_3_drop_names <- rownames(colData(sce[,sce$Sample == "3c_LS_SCP"])[which(sample_3_drop),])

#Concatenate all
cells_to_drop <- c(sample_1_drop_names,sample_2_drop_names,sample_3_drop_names)

#add information to the object
#sce$key is the rownames of the colData, which is the cell_id
message("identical(sce$key,rownames(colData(sce)))")
identical(sce$key,rownames(colData(sce)))

#Label cells to be removed. 
sce$discard_sample_specific <- ifelse(sce$key %in% cells_to_drop,
                                      TRUE,
                                      FALSE)

message("Number of total cells discarded")
table(sce$discard_sample_specific)

message("Percentage of non-empty droplets removed")
100 * sum(sce$discard_sample_specific) / ncol(sce)

message("Sample x discard table")
(qc_t <- addmargins(table(sce$Sample, sce$discard_sample_specific)))

message("Percentage of each sample dropped")
round(100 * sweep(qc_t, 1, qc_t[, 3], "/"), 1)

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


message("doublet score statistics")
summary(sce$doubletScore)


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

table(sce$discard_sample_specific, sce$doubletScore >= 5)

#Save object
save(sce,file = here("processed-data","02_build_sce","sce_emptyDrops_removed_withQC_081724.rda"))

# #10/23/23
# #load in the empty drops removed with QC object
# load(here("processed-data","sce_emptyDrops_removed_withQC.rda"))

#Also load in vector of cell IDs that need to be removed. 
load(here("processed-data","cluster_11_low_quality_IDs.rda"),verbose = TRUE)

message("How many low quality nuclei identified?")
length(low_quality_nuclei)

#Add the low quality nuclei to the discard 
sce[,low_quality_nuclei]$discard_sample_specific <- TRUE

message("How many total nuclei removed?")
table(sce$discard_sample_specific)

sce <- sce[,!sce$discard_sample_specific]
message("Dimensions of QCed object")
dim(sce)

sce

## save QCed and cleaned object. 
save(sce,file=here("processed-data","02_build_sce","sce_clean_081724.rda"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
