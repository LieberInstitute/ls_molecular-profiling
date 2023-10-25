#Goal: Compare gene expression signatures of the 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

library(SingleCellExperiment)
library(rafalib)
library(here)

#####load and prep data for the analysis######

#load the SingleCellExperiment object for human Lateral Septum
load(here("processed-data","sce_with_CellType.rda"))

sce
# class: SingleCellExperiment 
# dim: 33556 9225 
# metadata(1): Samples
# assays(3): counts binomial_deviance_residuals logcounts
# rownames(33556): ENSG00000243485 ENSG00000238009 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(9225): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 3_TTTGGTTTCTTCGACC-1 3_TTTGTTGTCCCGATCT-1
# colData names(60): Sample Barcode ... k_75_walktrap CellType
# reducedDimNames(14): GLMPCA_approx UMAP_15 ... UMAP_mnn_25 UMAP_mnn_50
# mainExpName: NULL
# altExpNames(0):

#load human DEG list from
load(here("processed-data","markers_1vAll_ttest.rda")) #markers_1vALL_df

#Figure out which genes have non0 means.
cellClust.idx <- splitit(sce$CellType)
medianNon0.human_ls <- lapply(cellClust.idx, function(x) {
    apply(as.matrix(assay(sce, "logcounts")), 1, function(y) {
        median(y[x]) > 0
    })
})

#load the SingleCellExperiment object for mouse Lateral Septum
load(here("MAGMA","mouse_analysis","sce_updated_LS.rda"))

sce.ls
# class: SingleCellExperiment 
# dim: 32285 22860 
# metadata(1): Samples
# assays(3): counts binomial_pearson_residuals logcounts
# rownames(32285): ENSMUSG00000051951 ENSMUSG00000089699 ...
# ENSMUSG00000095019 ENSMUSG00000095041
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(22860): 1_AAACCCAAGGTACATA-1 1_AAACCCACATCCGAGC-1 ...
# 4_TTTGTTGCATACAGCT-1 4_TTTGTTGGTCAAACGG-1
# colData names(17): Sample Barcode ... cellType.final cellType.broad
# reducedDimNames(4): GLMPCA_approx UMAP TSNE GLMPCA_50
# mainExpName: NULL
# altExpNames(0):

#Run mouse DEGs 
load(here("MAGMA","mouse_analysis","markers-stats_LS-n4_findMarkers_33cellTypes.rda"),verbose = TRUE)
# Loading objects:
#     markers.ls.t.pw
#     markers.ls.t.1vAll
#     medianNon0.ls

#DEGs generated with the cluster-vs-all method are in the markers.ls.t.1vAll object
#Add the medianNon0 information into the cluster.  
#Each iteration of the list contains two lists.
#Want the list "1" which includes genes that are enriched in the cluster. 
for(i in names(markers.ls.t.1vAll)){
    print(i)
    print(nrow(markers.ls.t.1vAll[[i]][["1"]]))
    markers.ls.t.1vAll[[i]][["1"]] <- cbind(
        markers.ls.t.1vAll[[i]][["1"]],
        medianNon0.ls[[i]][
            match(row.names(markers.ls.t.1vAll[[i]][["1"]]),
                  names(medianNon0.ls[[i]]))
        ]
    )
    colnames(markers.ls.t.1vAll[[i]][["1"]])[5] <- "non0Median"
    markers.ls.t.1vAll[[i]][["1"]] <- as.data.frame(markers.ls.t.1vAll[[i]][["1"]])
    markers.ls.t.1vAll[[i]][["1"]]$gene_id <- row.names(markers.ls.t.1vAll[[i]][["1"]])
    markers.ls.t.1vAll[[i]][["1"]] <- dplyr::left_join(x = markers.ls.t.1vAll[[i]][["1"]],
                                                       y = as.data.frame(rowData(sce.ls)[,c("gene_id","gene_name")]),
                                                       by = "gene_id")
    print(nrow(markers.ls.t.1vAll[[i]][["1"]]))
}











