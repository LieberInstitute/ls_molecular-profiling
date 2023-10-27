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
    rownames(markers_1vALL_list[[i]]) <- markers_1vALL_list[[i]]$gene_id
    
}

#Save the markers_1vALL_list 
save(markers_1vALL_list,file = here("processed-data","markers_1vAll_ttest_withnon0median.rda"))

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

#Lines 88-122 are straight from https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/51d15ef9f5f2c4c53f55e22e3fe467de1a724668/10x_NAc-n8_step04_cross-species_rnNAc_MNT.R#L4
#This cleans up the list nicely. 
markers.ls.t.1vAll <- lapply(markers.ls.t.1vAll, function(x) {
    # Basically take the 'stats.[1 or 0]' since is redundant with the 'summary'-level stats
    lapply(x, function(y) {
        y[, 4]
    })
})

# Re-name std.lfc column and the entries; add non-0-median info
for (i in names(markers.ls.t.1vAll)) {
    colnames(markers.ls.t.1vAll[[i]][["0"]])[1] <- "std.logFC"
    colnames(markers.ls.t.1vAll[[i]][["1"]])[1] <- "std.logFC"
    # Add non0median Boolean - might be informative for both sets of stats
    markers.ls.t.1vAll[[i]][["0"]] <- cbind(
        markers.ls.t.1vAll[[i]][["0"]],
        medianNon0.ls[[i]][match(
            rownames(markers.ls.t.1vAll[[i]][["0"]]),
            names(medianNon0.ls[[i]])
        )]
    )
    colnames(markers.ls.t.1vAll[[i]][["0"]])[4] <- "non0median"
    
    # "1" aka 'enriched'
    markers.ls.t.1vAll[[i]][["1"]] <- cbind(
        markers.ls.t.1vAll[[i]][["1"]],
        medianNon0.ls[[i]][match(
            rownames(markers.ls.t.1vAll[[i]][["1"]]),
            names(medianNon0.ls[[i]])
        )]
    )
    colnames(markers.ls.t.1vAll[[i]][["1"]])[4] <- "non0median"
    #Add additional columns
    markers.ls.t.1vAll[[i]][["1"]]$gene_id <- rownames(markers.ls.t.1vAll[[i]][["1"]]) 
    markers.ls.t.1vAll[[i]][["1"]] <- dplyr::left_join(x = as.data.frame(markers.ls.t.1vAll[[i]][["1"]]),
                                                       y = as.data.frame(rowData(sce.ls)[,c("gene_id","gene_name")]),
                                                       by = "gene_id")
    rownames(markers.ls.t.1vAll[[i]][["1"]]) <- markers.ls.t.1vAll[[i]][["1"]]$gene_id
    # Then re-name the entries to more interpretable, because we'll keeping both contrasts
    names(markers.ls.t.1vAll[[i]]) <- paste0(i, c("_depleted", "_enriched"))
}

#save the list. 
save(markers.ls.t.1vAll,file = here("processed-data","mouse_markers_1vAll_ttest_withnon0median.rda"))

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
            file = here("processed-data","HOM_AllOrganism_JAX_102723.csv"),
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

#Subset for mouse
hom_mm <- hom[hom$Common.Organism.Name == "mouse, laboratory", ]

table(rowData(sce_mouse_ls)$mm.entrezIds %in% hom_mm$EntrezGene.ID)
# FALSE  TRUE 
# 8990 18807 

#Add the IDs to the sce_human_ls object. 
rowData(sce_mouse_ls)$JAX.geneID <- hom_mm$DB.Class.Key[match(rowData(sce_mouse_ls)$mm.entrezIds,hom_mm$EntrezGene.ID)]

##Identify shared genes
length(intersect(rowData(sce_human_ls)$JAX.geneID,
                 rowData(sce_mouse_ls)$JAX.geneID)) 
# [1] 16573
#16,573 shared genes. 

shared_homologs <- intersect(rowData(sce_human_ls)$JAX.geneID,
                             rowData(sce_mouse_ls)$JAX.geneID)
shared_homologs <- shared_homologs[-1] #first is na

# Human not in mouse
length(setdiff(rowData(sce_human_ls)$JAX.geneID,
               rowData(sce_mouse_ls)$JAX.geneID)) 
# [1] 460

# Mouse not in human
length(setdiff(rowData(sce_mouse_ls)$JAX.geneID,
               rowData(sce_human_ls)$JAX.geneID)) 
#[1] 2213

# Subset for shared homologs
sce_human_sub <- sce_human_ls[rowData(sce_human_ls)$JAX.geneID %in% shared_homologs, ]
dim(sce_human_sub)
# [1] 16996  9225
#16996 genes and 9225 cells

sce_mouse_sub <- sce_mouse_ls[rowData(sce_mouse_ls)$JAX.geneID %in% shared_homologs, ]
dim(sce_mouse_sub)
# [1] 16588 22860
#16588 genes adn 22860 cells

#Are any of the JAX IDs duplicated
length(rowData(sce_human_sub)$gene_name[duplicated(rowData(sce_human_sub)$JAX.geneID)])
#[1] 424

length(rowData(sce_mouse_sub)$gene_name[duplicated(rowData(sce_mouse_sub)$JAX.geneID)])
# [1] 16

######Need to identify genes that are duplicated and keep the higher expressing. 
###Mouse first. 
mouse_dups <- rowData(sce_mouse_sub)[which(duplicated(rowData(sce_mouse_sub)$JAX.geneID)),"JAX.geneID"]
mouse_genes_to_compare <- list()
mouse_genes_to_keep <- character()
for(i in 1:length(mouse_dups)){
    print(i)
    mouse_genes_to_compare[[i]] <- rownames(sce_mouse_sub)[rowData(sce_mouse_sub)$JAX.geneID == mouse_dups[i]]
    rowmeans_dups <- rowMeans(assay(sce_mouse_sub[mouse_genes_to_compare[[i]], ], "logcounts"))
    mouse_genes_to_keep[i] <- names(rowmeans_dups[order(rowmeans_dups, decreasing=TRUE)])[1]
}

#Get the genes that were not duplicated. 
non_dups_mouse <- rownames(sce_mouse_sub)[!(rownames(sce_mouse_sub) %in% unlist(mouse_genes_to_compare))]

# Finally combine and subset
sce_mouse_sub <- sce_mouse_sub[c(non_dups_mouse, unique(mouse_genes_to_keep)), ]

table(rowData(sce_mouse_sub)$JAX.geneID %in% shared_homologs)
# TRUE 
# 16572 

table(duplicated(rowData(sce_mouse_sub)$JAX.geneID))
# FALSE 
# 16572 

#Nothing is duplicated so can move forward. 

###human first. 
human_dups <- rowData(sce_human_sub)[which(duplicated(rowData(sce_human_sub)$JAX.geneID)),"JAX.geneID"]
human_genes_to_compare <- list()
human_genes_to_keep <- character()
for(i in 1:length(human_dups)){
    print(i)
    human_genes_to_compare[[i]] <- rownames(sce_human_sub)[rowData(sce_human_sub)$JAX.geneID == human_dups[i]]
    rowmeans_dups <- rowMeans(assay(sce_human_sub[human_genes_to_compare[[i]], ], "logcounts"))
    human_genes_to_keep[i] <- names(rowmeans_dups[order(rowmeans_dups, decreasing=TRUE)])[1]
}

#Get the genes that were not duplicated. 
non_dups_human <- rownames(sce_human_sub)[!(rownames(sce_human_sub) %in% unlist(human_genes_to_compare))]

# Finally combine and subset
sce_human_sub <- sce_human_sub[c(non_dups_human, unique(human_genes_to_keep)), ]

table(rowData(sce_human_sub)$JAX.geneID %in% shared_homologs)
# TRUE 
# 16572 

table(duplicated(rowData(sce_human_sub)$JAX.geneID))
# FALSE 
# 16572 

#Nothing is duplicated so can move forward. 

## Match order and save
sce_mouse_sub <- sce_mouse_sub[match(rowData(sce_human_sub)$JAX.geneID,
                                     rowData(sce_mouse_sub)$JAX.geneID), ]

#sanity_check
all(rowData(sce_mouse_sub)$JAX.geneID == rowData(sce_human_sub)$JAX.geneID)
# [1] TRUE

##Calculate the t statistic
#Mouse first. 
#Take the _enriched dataframe within this list. 
mouse_enriched <- lapply(markers.ls.t.1vAll, function(x){x[[2]]})

#Now calculate
fixTo <- rownames(mouse_enriched[["Astro"]])
for(x in names(mouse_enriched)){
    print(x)
    mouse_enriched[[x]]$t.stat <- mouse_enriched[[x]]$std.logFC * sqrt(ncol(sce_mouse_sub))
    mouse_enriched[[x]] <- mouse_enriched[[x]][fixTo, ]
}

# Pull out the t's
mouse_t_stats <- sapply(mouse_enriched, function(x){x$t.stat})
rownames(mouse_t_stats) <- fixTo

table(rownames(sce_mouse_sub) %in% rownames(mouse_t_stats))
# FALSE  TRUE 
# 10     16562 
#There are 10 genes within the sce_mouse_sub object that are not within the mouse_t_statistics. 
#Pull those genes. 
m_missing <- rownames(sce_mouse_sub)[!(rownames(sce_mouse_sub) %in% rownames(mouse_t_stats))]
#Are they found in the DEG list? 
lapply(mouse_enriched,function(x){ table(m_missing %in% rownames(x)) }) #All FALSE

# Subset for those with homologous genes in human
mouse_t_stats_hom <- mouse_t_stats[rownames(mouse_t_stats) %in% rownames(sce_mouse_sub), ]
dim(mouse_t_stats_hom)
# [1] 16562    33

# Need to change the rownames to JAX.geneID so we can assure everything is in the same order. 
mouse_t_stats_df <- as.data.frame(mouse_t_stats_hom)
mouse_t_stats_df$gene_id <- rownames(mouse_t_stats_df)
#Add JAx ID to the dataframe
mouse_t_stats_df <- dplyr::left_join(x = mouse_t_stats_df,
                                     y = as.data.frame(rowData(sce_mouse_sub)[,c("gene_id","JAX.geneID")]),
                                     by = "gene_id")
#Make rownames the JAx ID
rownames(mouse_t_stats_df) <- mouse_t_stats_df$JAX.geneID
#Convert it to a matrix. 
mouse_t_stats_mat <- as.matrix(mouse_t_stats_df[,c(1:33)]) #Just keeping the celltype columns and removing both gene columns


#Human now 
#Now calculate
fixTo <- rownames(markers_1vALL_list[["LS_GABA_1"]])
for(x in names(markers_1vALL_list)){
    print(x)
    markers_1vALL_list[[x]]$t.stat <- markers_1vALL_list[[x]]$std.logFC * sqrt(ncol(sce_human_sub))
    markers_1vALL_list[[x]] <- markers_1vALL_list[[x]][fixTo, ]
}


# Pull out the t's
human_t_stats <- sapply(markers_1vALL_list, function(x){x$t.stat})
rownames(human_t_stats) <- fixTo

table(rownames(sce_human_sub) %in% rownames(human_t_stats))
# TRUE 
# 16572 

human_t_stats_hom <- human_t_stats[rownames(sce_human_sub),]

# Need to change the rownames to JAX.geneID so we can assure everything is in the same order. 
human_t_stats_df <- as.data.frame(human_t_stats_hom)
human_t_stats_df$gene_id <- rownames(human_t_stats_df)
#Add JAx ID to the dataframe
human_t_stats_df <- dplyr::left_join(x = human_t_stats_df,
                                     y = as.data.frame(rowData(sce_human_sub)[,c("gene_id","JAX.geneID")]),
                                     by = "gene_id")
#Make rownames the JAx ID
rownames(human_t_stats_df) <- human_t_stats_df$JAX.geneID
#Convert it to a matrix. 
human_t_stats_mat <- as.matrix(human_t_stats_df[,c(1:15)]) #Just keeping the celltype columns and removing both gene columns


#Save the subset objects and t stat matrices. 
save(human_t_stats_mat,mouse_t_stats_mat,file = here("processed-data","t_stats_mats.rda"))
save(sce_mouse_sub,sce_human_sub,file = here("processed-data","human_mouse_matched_by_JAX.rda"))

#the mouse t stat matrix has 10 less genes. Need to subset the human mat to have the same genes
human_t_stats_mat <- human_t_stats_mat[rownames(mouse_t_stats_mat),]

all(rownames(human_t_stats_mat) == rownames(mouse_t_stats_mat))
# [1] TRUE
#everything in the same order. 

#Correlate all of the homologs. 
cor_t_all <- cor(human_t_stats_mat, mouse_t_stats_mat)
rownames(cor_t_all) <- paste0(rownames(cor_t_all),"_Human")
colnames(cor_t_all) <- paste0(colnames(cor_t_all),"_Mouse")
range(cor_t_all) 
# [1] -0.3748624  0.6449819


#Get top 100 genes for each human LS cluster. Top chosen by t-statistic. 
human_top_100 <- mapply(as.data.frame(human_t_stats_mat), FUN = function(t) {
    o <- order(t, decreasing = TRUE)[1:100]
})

#Now top 100 for each mouse. 
mouse_top_100 <- mapply(as.data.frame(mouse_t_stats_mat), FUN = function(t) {
    o <- order(t, decreasing = TRUE)[1:100]
})

#get the unique identifiers for each species plus the shared. 
human_unique <- unique(as.numeric(human_top_100))
length(human_unique)
#[1] 1235

mouse_unique <- unique(as.numeric(mouse_top_100))
length(mouse_unique)
#[1] 2000

shared_identifiers <- intersect(rownames(human_t_stats_mat)[human_unique], 
                                rownames(mouse_t_stats_mat)[mouse_unique])
length(shared_identifiers)
#706


#Correlate with just the human identifiers. 
cor_t_human_unique <- cor(human_t_stats_mat[human_unique, ],
                          mouse_t_stats_mat[human_unique, ])
rownames(cor_t_human_unique) <- paste0(rownames(cor_t_human_unique),"_Human")
colnames(cor_t_human_unique) <- paste0(colnames(cor_t_human_unique),"_Mouse")
range(cor_t_human_unique)
# [1] -0.4211243  0.7728237


#Correlate with just the mouse identifiers. 
cor_t_mouse_unique <- cor(human_t_stats_mat[mouse_unique, ],
                          mouse_t_stats_mat[mouse_unique, ])
rownames(cor_t_mouse_unique) <- paste0(rownames(cor_t_mouse_unique),"_Human")
colnames(cor_t_mouse_unique) <- paste0(colnames(cor_t_mouse_unique),"_Mouse")
range(cor_t_mouse_unique)
# [1] -0.4143752  0.7664956

#Correlate with just the shared identifiers. 
cor_t_shared <- cor(human_t_stats_mat[shared_identifiers, ],
                    mouse_t_stats_mat[shared_identifiers, ])
rownames(cor_t_shared) <- paste0(rownames(cor_t_shared),"_Human")
colnames(cor_t_shared) <- paste0(colnames(cor_t_shared),"_Mouse")
range(cor_t_shared)
# [1] -0.4683984  0.8702079

#Save all of the correlation matrices. 
save(cor_t_all,cor_t_human_unique,cor_t_mouse_unique,cor_t_shared,
     file = here("processed-data","correlation_matrices_conservation_analysis.rda"))


