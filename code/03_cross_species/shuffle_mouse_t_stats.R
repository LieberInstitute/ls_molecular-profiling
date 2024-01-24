library(SingleCellExperiment)
library(sessioninfo)
library(scater)
library(scran)
library(here)


load(here("processed-data","human_mouse_matched_by_JAX.rda"),verbose = TRUE)

#Randomly sample 18.3% of the human LS dataset, 500x
size_m <- round(.0841*ncol(sce_mouse_sub))

#Make a dataframe to output the t-statistics into
out_m_dataframe <- as.data.frame(matrix(nrow = nrow(sce_mouse_sub),ncol = 500))
rownames(out_m_dataframe) <- rownames(sce_mouse_sub)
colnames(out_m_dataframe) <- 1:500

#Randomly sample the cells 500x
set.seed(150001110)
random_cells <- replicate(500,sample(x = 1:9225,size = size_m,replace = FALSE))
random_cellnames <- as.data.frame(matrix(nrow = nrow(random_cells),ncol = 500))
colnames(random_cellnames) <- 1:500
colData(sce_mouse_sub)$unique_rowname <- rownames(colData(sce_mouse_sub))
for(i in 1:500){
    random_cellnames[,i] <- colData(sce_mouse_sub)$unique_rowname[random_cells[,i]]
}

#Calculate the t-statistics
for(i in 1:500){
    #Set cell names
    sce_mouse_sub$random_designation <- ifelse(colData(sce_mouse_sub)$unique_rowname %in% random_cellnames[,i],
                                               "random",
                                               "other")
    
    #DEG testing to get std.logFC 
    pd <- as.data.frame(SummarizedExperiment::colData(sce_mouse_sub))
    mod <- with(pd, stats::model.matrix(as.formula("~Sample")))
    mod <- mod[, -1, drop = F]
    DEGs <- scran::findMarkers(sce_mouse_sub,
                               groups = sce_mouse_sub$random_designation,
                               assay.type = "logcounts", 
                               design = mod, 
                               test.type = "t",
                               log.p = TRUE,
                               std.lfc = TRUE,
                               direction = "up", 
                               pval.type = "all", 
                               full.stats = T)
    
    #Pull stats for random cells
    DEGs <- DEGs[["random"]]
    
    #Calculate t-statistic
    DEGs$std.logFC <- DEGs$summary.stats 
    DEGs$t.stat <-  DEGs$std.logFC * sqrt(ncol(sce_mouse_sub))
    
    #Input the t-statistic into the dataframe
    out_m_dataframe[,i] <- DEGs[rownames(out_m_dataframe),"t.stat"]
    print(i)
}


save(out_m_dataframe,file = here("processed-data","randomized_t_stats_mouse_500.rda"))

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()