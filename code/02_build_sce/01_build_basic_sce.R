#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling
#Code modified from https://github.com/LieberInstitute/spatial_hpc/blob/main/snRNAseq_hpc/code/build_sce/01_build_basic_sce.R
#RP

library(SingleCellExperiment)
library(DropletUtils)
library(here)
library(rtracklayer)
library(dplyr)
library(scuttle)
library(sessioninfo)

#Read in a dataframe consisting of identifying information for the samples
sample_data <- read.delim(here("tables",
                               "ls_molecular_profiling_sample_info.csv"),
                          header = TRUE,sep = ",")

#Check that all files exist. Also, stop if they do not. 
all(file.exists(sample_data$Raw_data_path))

stopifnot(all(file.exists(sample_data$Raw_data_path)))

## Build basic SCE
message("Read 10x data and create sce - ", Sys.time())

sce <- read10xCounts(samples = sample_data$Raw_data_path,
                     sample.names = sample_data$Sample_ID,
                     type = "sparse",
                     col.names = TRUE)

message("RDone - ", Sys.time())

sce

#Add information about the study design to the colData
#merging removes the rownames that are unique to each sample/cell. 
#create a column that will remain after merging.
colData(sce)$key <- rownames(colData(sce))

#Then merge and restore the rownames. 
new_column_data<- merge(x = colData(sce),
                        y = sample_data,
                        by.x = "Sample",
                        by.y = "Sample_ID")

#reorder the colData so that it is in the original order. 
new_column_data <- new_column_data[match(sce$key, new_column_data$key), ]

#Check that the key column within the new column data dataframe and sce are in the same order
#Stop if not
identical(sce$key,new_column_data$key)

stopifnot(identical(sce$key, new_column_data$key))

#Everything is in order. 
#Change the rownames back
rownames(new_column_data) <- new_column_data$key

#Update the column data. 
colData(sce) <- new_column_data

#rownames of the coldata needs to be in same order as column names of the count matrix
identical(rownames(colData(sce)),colnames(sce))

stopifnot(identical(rownames(colData(sce)),colnames(sce)))

#now update the rowData
gtf <- rtracklayer::import("/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf")
gtf <- gtf[gtf$type == "gene"]
names(gtf) <- gtf$gene_id

#match the genes
match_genes <- match(rownames(sce),gtf$gene_id)
stopifnot(all(!is.na(match_genes)))

#Keep only specific columns from the gtf
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_version", "gene_name", "gene_type")]

#Add gene info
rowRanges(sce) <- gtf[match_genes]

#Print the object again. 
message("Printing raw sce object - ", Sys.time())
sce

#Save object. 
save(sce,
     file = here("processed-data","02_build_sce","sce_raw.rda"))

#Empty droplets have not been removed. Will be the next step of the analysis. 

#sessionInfo
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
