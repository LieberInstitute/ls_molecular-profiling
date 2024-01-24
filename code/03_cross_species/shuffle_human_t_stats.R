library(SingleCellExperiment)
library(sessioninfo)
library(scater)
library(scran)
library(here)


load(here("processed-data","human_mouse_matched_by_JAX.rda"),verbose = TRUE)


size_h <- round(.183*ncol(sce_human_sub))


#Make a dataframe to output the t-statistics into
out_h_dataframe <- as.data.frame(matrix(nrow = nrow(sce_human_sub),ncol = 100))
rownames(out_h_dataframe) <- rownames(sce_human_sub)
colnames(out_h_dataframe) <- 1:500

#Randomly sample the cells 100x
set.seed(001110001)
random_cells <- replicate(500,sample(x = 1:9225,size = size_h,replace = FALSE))
random_cellnames <- as.data.frame(matrix(nrow = nrow(random_cells),ncol = 500))
colnames(random_cellnames) <- 1:500
for(i in 1:500){
    random_cellnames[,i] <- colData(sce_human_sub)$unique_rowname[random_cells[,i]]
}

#Calculate the t-statistics
for(i in 1:500){
    #Set cell names
    sce_human_sub$random_designation <- ifelse(colData(sce_human_sub)$unique_rowname %in% random_cellnames[,i],
                                               "random",
                                               "other")
    #DEG testing to get std.logFC 
    pd <- as.data.frame(SummarizedExperiment::colData(sce_human_sub))
    mod <- with(pd, stats::model.matrix(as.formula("~Sample")))
    mod <- mod[, -1, drop = F]
    DEGs <- scran::findMarkers(sce_human_sub,
                               groups = sce_human_sub$random_designation,
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
    DEGs$t.stat <-  DEGs$std.logFC * sqrt(ncol(sce_human_sub))
    
    #Input the t-statistic into the dataframe
    out_h_dataframe[,i] <- DEGs[rownames(out_h_dataframe),"t.stat"]
    print(i)
}

save(out_h_dataframe,file = here("processed-data","randomized_t_stats_human_500.rda"))

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()