library(SingleCellExperiment)
library(DeconvoBuddies)
library(sessioninfo)
library(scater)
library(scran)
library(here)


load(here("processed-data","sce_mouse_sub.rda"),verbose = TRUE)

sce_mouse_sub$key <- rownames(colData(sce_mouse_sub))

#Number of LS cells in the dataset 
size_m <- 1841

#Make a dataframe to output the std.logFC into
out_m_dataframe <- as.data.frame(matrix(nrow = nrow(sce_mouse_sub),ncol = 500))
rownames(out_m_dataframe) <- rownames(sce_mouse_sub)
colnames(out_m_dataframe) <- 1:500

#Randomly sample the cells 500x
set.seed(15001110)
random_cells <- replicate(500,sample(x = 1:21884,size = size_m,replace = TRUE))
random_cellnames <- as.data.frame(matrix(nrow = nrow(random_cells),ncol = 500))
colnames(random_cellnames) <- 1:500

for(i in 1:500){
  random_cellnames[,i] <- colData(sce_mouse_sub)[random_cells[,i],"key"]
}

#Calculate the std.logFC
for(i in 1:500){
  #Set cell names
  sce_mouse_sub$random_designation <- ifelse(colData(sce_mouse_sub)$key %in% random_cellnames[,i],
                                             "random",
                                             "other")
  
  #DEG testing to get std.logFC 
  DEGs <- findMarkers_1vAll(sce_mouse_sub,
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
  DEGs <- DEGs[match(rownames(out_m_dataframe),rownames(DEGs)),]
   
  stopifnot(identical(rownames(out_m_dataframe),rownames(DEGs)))
  
  #Input the t-statistic into the dataframe
  out_m_dataframe[,i] <- DEGs$std.logFC
  print(i)
}

save(out_m_dataframe,file = here("processed-data","randomized_std.logFC_mouse_500.rda"))

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
