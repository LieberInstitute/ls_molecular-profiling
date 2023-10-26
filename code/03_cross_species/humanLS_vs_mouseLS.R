#Goal: Compare gene expression signatures of the mouse and human LS
#Will use two different approaches: 
#1. Correlation of t-statistics for markers genes. Code modified from https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/51d15ef9f5f2c4c53f55e22e3fe467de1a724668/10x_NAc-n8_step04_cross-species_rnNAc_MNT.R#L4
#2. Integrate human and mouse data with mnn. 
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

library(SingleCellExperiment)
library(sparseMatrixStats)
library(org.Hs.eg.db) #human
library(org.Mm.eg.db) #mouse
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
dim(markers_1vALL_df)
# [1] 503340      9

#Split the markers_1vALL_df into a lsit
markers_1vALL_list <- split(markers_1vALL_df,
                            f = markers_1vALL_df$cellType.target)

#Figure out which genes have non0 medians
cellClust.idx <- splitit(sce$CellType)
non0median_human_ls <- vector(mode = "list",length = 15)
names(non0median_human_ls) <- names(cellClust.idx)
for(i in names(cellClust.idx)){
    print(i)
    #Create dataframe that consists of cell type, ensembl gene id, and whether gene has a non0median or not. 
    non0median_human_ls[[i]] <- data.frame(celltype = i,
                                           gene_id = rownames(sce[,cellClust.idx[[i]]]), 
                                           non0median = rowMedians(assay(sce[,cellClust.idx[[i]]],"logcounts")) > 0)
    #add non0median information oto the list of dataframes. 
    markers_1vALL_list[[i]] <- dplyr::left_join(x  = markers_1vALL_list[[i]],
                                                y  = non0median_human_ls[[i]],
                                                by = "gene_id")
    markers_1vALL_list[[i]]$FDR <- exp(markers_1vALL_list[[i]]$log.FDR)
    
}


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
}

#rename the objects to make more sense. 
sce_human_ls <- sce 
sce_mouse_ls <- sce.ls
rm(sce,sce.ls)


############Correlation of t-statistics for markers genes###############

#Do the objects contain any genes with no counts? 
table(rowSums(assay(sce_human_ls, "counts"))==0)
# FALSE 
# 33556 

table(rowSums(assay(sce_mouse_ls, "counts"))==0)
# FALSE  TRUE 
# 27797  4488 

#Remove genes with all 0s from the mouse object. 
sce_mouse_ls <- sce_mouse_ls[!rowSums(assay(sce_mouse_ls, "counts"))==0, ]

#sanity check
table(rowSums(assay(sce_mouse_ls, "counts"))==0)
# FALSE 
# 27797 

##Add entrez gene ids and jax ids to the human object
#Entrez ids 
hs.entrezIds <- mapIds(org.Hs.eg.db, 
                       keys=rowData(sce_human_ls)$gene_id, 
                       column="ENTREZID", 
                       keytype="ENSEMBL")
# 'select()' returned 1:many mapping between keys and columns

table(is.na(hs.entrezIds))
# FALSE  TRUE 
# 22561 10995 
#22561 valid entries. 

#Which genes do not have entrez ids
withoutEntrez <- names(hs.entrezIds)[is.na(hs.entrezIds)]
names(withoutEntrez) <- rowData(sce_human_ls)[rowData(sce_human_ls)$gene_id %in% withoutEntrez, ]$gene_name

#Add entrez ids to the object
rowData(sce_human_ls) <- cbind(rowData(sce_human_ls), hs.entrezIds)

# JAX annotation info
hom <-  read.delim("http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt",
                   as.is=TRUE)

#Save dataframe with date. In case we need to use later and it gets updated. 
write.table(x = hom,
            file = here("processed-data","HOM_AllOrganism_JAX_102623.csv"),
            sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE)

#Subset for human 
hom_hs <- hom[hom$Common.Organism.Name == "human", ]

table(rowData(sce_human_ls)$hs.entrezIds %in% hom_hs$EntrezGene.ID)
# FALSE  TRUE 
# 16025 17531 

#Add the IDs to the sce_human_ls object. 
rowData(sce_human_ls)$JAX.geneID <- hom_hs$DB.Class.Key[match(rowData(sce_human_ls)$hs.entrezIds,hom_hs$EntrezGene.ID)]

##Add entrez gene ids and jax ids to the mouse object
#Entrez ids 
mm.entrezIds <- mapIds(org.Mm.eg.db, 
                       keys=rowData(sce_mouse_ls)$gene_id, 
                       column="ENTREZID", 
                       keytype="ENSEMBL")
# 'select()' returned 1:many mapping between keys and columns

table(is.na(mm.entrezIds))
# FALSE  TRUE 
# 21576  6221 

#Which genes do not have entrez ids
withoutEntrez_mouse <- names(mm.entrezIds)[is.na(mm.entrezIds)]
names(withoutEntrez_mouse) <- rowData(sce_mouse_ls)[rowData(sce_mouse_ls)$gene_id %in% withoutEntrez_mouse, ]$gene_name

#Add entrez ids to the object
rowData(sce_mouse_ls) <- cbind(rowData(sce_mouse_ls), mm.entrezIds)






