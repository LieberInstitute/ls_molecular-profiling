#srun --cpus-per-task=2  --mem=25G --pty bash
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling
#Code modified from https://github.com/LieberInstitute/spatial_hpc/blob/main/snRNAseq_hpc/code/build_sce/01_build_basic_sce.R

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

## Build basic SCE
message("Read 10x data and create sce - ", Sys.time())
#Read 10x data and create sce - 2023-10-12 10:07:33.021133
sce <- read10xCounts(samples = sample_data$Raw_data_path,
                     sample.names = sample_data$Sample_ID,
                     type = "sparse",
                     col.names = TRUE)
message("RDone - ", Sys.time())
#RDone - 2023-10-12 10:10:00.252362

sce
# class: SingleCellExperiment 
# dim: 36601 4909214 
# metadata(1): Samples
# assays(1): counts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(3): ID Symbol Type
# colnames(4909214): 1_AAACCCAAGAAACCAT-1 1_AAACCCAAGAAACCCA-1 ...
# 3_TTTGTTGTCTTTGCGC-1 3_TTTGTTGTCTTTGGAG-1
# colData names(2): Sample Barcode
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):

save(sce,
     file = here("processed-data","sce_raw.rda"))

#Add information about the study design to the colData
#merging removes the rownames that are unique to each sample/cell. 
#create a column that will remain after merging. 
colData(sce)$unique_rowname <- rownames(colData(sce))

#Then merge and restore the rownames. 
new_column_data<- merge(x = colData(sce),
                        y = sample_data[,-which(colnames(sample_data) == "Raw_data_path")],
                        by.x = "Sample",
                        by.y = "Sample_ID")
rownames(new_column_data) <- new_column_data$unique_rowname

#Update the column data. 
colData(sce) <- new_column_data

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

#Add QC information about each feature;

#Overwrite previous save. 
save(sce,
     file = here("processed-data","sce_raw.rda"))

#Empty droplets

