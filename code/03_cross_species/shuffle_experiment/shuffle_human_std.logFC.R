library(SingleCellExperiment)
library(DeconvoBuddies)
library(sessioninfo)
library(scater)
library(scran)
library(here)

load(here("processed-data","sce_human_sub.rda"),verbose = TRUE)

#Number of LS cells in the dataset 
size_h <- 1680

#Make a dataframe to output the std.logFC into
out_h_dataframe <- as.data.frame(matrix(nrow = nrow(sce_human_sub),ncol = 500))
rownames(out_h_dataframe) <- rownames(sce_human_sub)
colnames(out_h_dataframe) <- 1:500

#Randomly sample the cells 500x
set.seed(001110001)
random_cells <- replicate(500,sample(x = 1:9225,size = size_h,replace = TRUE))
random_cellnames <- as.data.frame(matrix(nrow = nrow(random_cells),ncol = 500))
colnames(random_cellnames) <- 1:500
for(i in 1:500){
  random_cellnames[,i] <- colData(sce_human_sub)[random_cells[,i],"key"]
}

#Calculate the std.logFC
for(i in 1:500){
  #Set cell names
  sce_human_sub$random_designation <- ifelse(colData(sce_human_sub)$key %in% random_cellnames[,i],
                                             "random",
                                             "other")
  #DEG testing to get std.logFC 
  DEGs <- findMarkers_1vAll(sce_human_sub,
                            assay_name = "logcounts",
                            cellType_col = "random_designation",
                            direction = "up",
                            mod = "~Sample")

  DEGs <- as.data.frame(DEGs)
    
  #subet for random DEGs
  DEGs <- subset(DEGs,subset=(cellType.target == "random"))
   
  #Make the rownames the gene id
  rownames(DEGs) <- DEGs$gene
  
  print(head(DEGs))
  
  #Force the rows to be in the same order. 
  DEGs <- DEGs[match(rownames(out_h_dataframe),rownames(DEGs)),]
  
  stopifnot(identical(rownames(out_h_dataframe),rownames(DEGs)))

#Input the std.logFC  into the dataframe
out_h_dataframe[,i] <- DEGs$std.logFC
print(i)
}

save(out_h_dataframe,file = here("processed-data","randomized_std.logFC_human_500.rda"))

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
