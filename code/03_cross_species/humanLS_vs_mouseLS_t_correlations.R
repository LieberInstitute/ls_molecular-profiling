#Goal: Compare gene expression signatures of the mouse and human LS
#This analysis will focus on correlation of t-statistics for markers genes. 
#Code modified from https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/51d15ef9f5f2c4c53f55e22e3fe467de1a724668/10x_NAc-n8_step04_cross-species_rnNAc_MNT.R#L4
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/

library(SingleCellExperiment)
library(sparseMatrixStats)
library(RColorBrewer)
library(org.Hs.eg.db) #human
library(org.Mm.eg.db) #mouse
library(sessioninfo)
library(pheatmap)
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
# colData names(61): Sample Barcode ... CellType_k_20_louvain
# CellType.Final
# reducedDimNames(18): GLMPCA_approx UMAP_15 ... tSNE_mnn_25 tSNE_mnn_50
# mainExpName: NULL
# altExpNames(0):

#load human DEG list from human
load(here("processed-data","markers_1vAll_ttest_CellTypeFinal_20Clusters.rda"),verbose = TRUE) #markers_1vALL_enrich_Final
dim(markers_1vALL_enrich_Final)
# [1] 671120      9

#Split the markers_1vALL_df into a lsit
markers_1vALL_list <- split(markers_1vALL_enrich_Final,
                            f = markers_1vALL_enrich_Final$cellType.target)

#Figure out which genes have non0 medians
cellClust.idx <- splitit(sce$CellType.Final)
non0median_human_ls <- vector(mode = "list",length = length(unique(sce$CellType.Final)))
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
load(file = "/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/processed_data/SCE/sce_updated_LS.rda")

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
load(file = "/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/processed_data/SCE/markers-stats_LS-n4_findMarkers_33cellTypes.rda",
     verbose = TRUE)
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
            file = here("processed-data","HOM_AllOrganism_JAX_110923.csv"),
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
# [1] 459

# Mouse not in human
length(setdiff(rowData(sce_mouse_ls)$JAX.geneID,
               rowData(sce_human_ls)$JAX.geneID)) 
#[1] 2213

# Subset for shared homologs
sce_human_sub <- sce_human_ls[rowData(sce_human_ls)$JAX.geneID %in% shared_homologs, ]
dim(sce_human_sub)
# [1] 16997  9225
#16997 genes and 9225 cells

sce_mouse_sub <- sce_mouse_ls[rowData(sce_mouse_ls)$JAX.geneID %in% shared_homologs, ]
dim(sce_mouse_sub)
#[1] 16588 22860
#16588 genes and 22860 cells

#Are any of the JAX IDs duplicated
length(rowData(sce_human_sub)$gene_name[duplicated(rowData(sce_human_sub)$JAX.geneID)])
#[1] 425

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
fixTo <- rownames(markers_1vALL_list[["LS_Inh_A"]])
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
human_t_stats_mat <- as.matrix(human_t_stats_df[,c(1:20)]) #Just keeping the celltype columns and removing both gene columns


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
#[1] -0.4335106  0.6448713


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
#[1] 1546

mouse_unique <- unique(as.numeric(mouse_top_100))
length(mouse_unique)
#[1] 2000

shared_identifiers <- intersect(rownames(human_t_stats_mat)[human_unique], 
                                rownames(mouse_t_stats_mat)[mouse_unique])
length(shared_identifiers)
#[1] 865


#Correlate with just the human identifiers. 
cor_t_human_unique <- cor(human_t_stats_mat[human_unique, ],
                          mouse_t_stats_mat[human_unique, ])
rownames(cor_t_human_unique) <- paste0(rownames(cor_t_human_unique),"_Human")
colnames(cor_t_human_unique) <- paste0(colnames(cor_t_human_unique),"_Mouse")
range(cor_t_human_unique)
# [1] -0.5175540  0.7730306


#Correlate with just the mouse identifiers. 
cor_t_mouse_unique <- cor(human_t_stats_mat[mouse_unique, ],
                          mouse_t_stats_mat[mouse_unique, ])
rownames(cor_t_mouse_unique) <- paste0(rownames(cor_t_mouse_unique),"_Human")
colnames(cor_t_mouse_unique) <- paste0(colnames(cor_t_mouse_unique),"_Mouse")
range(cor_t_mouse_unique)
# [1] -0.5078417  0.7664956

#Correlate with just the shared identifiers. 
cor_t_shared <- cor(human_t_stats_mat[shared_identifiers, ],
                    mouse_t_stats_mat[shared_identifiers, ])
rownames(cor_t_shared) <- paste0(rownames(cor_t_shared),"_Human")
colnames(cor_t_shared) <- paste0(colnames(cor_t_shared),"_Mouse")
range(cor_t_shared)
# [1] -0.5626473  0.8641498

#Save all of the correlation matrices. 
save(cor_t_all,cor_t_human_unique,cor_t_mouse_unique,cor_t_shared,
     file = here("processed-data","correlation_matrices_conservation_analysis.rda"))

#######################################
############ Make heatmaps ############
#######################################
#First for correlations between all 16000+ genes. 
colrange <-  seq(-.65,.65, by = 0.01)
colorpal <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(colrange))

#Rearrange the matrix
cor_t_all <- cor_t_all[rownames(cor_t_all)[order(rownames(cor_t_all))],
                       colnames(cor_t_all)[order(colnames(cor_t_all))]]
pdf(file = here("plots","Conservation","Human_Mouse_allHomologs_t_correlation_HM.pdf"),
    height = 12,width = 12)
pheatmap(cor_t_all,
         color=colorpal,
         cluster_cols=FALSE, 
         cluster_rows=FALSE,
         breaks=colrange,
         fontsize=11, 
         fontsize_row=11.5, 
         fontsize_col=12,
         display_numbers=T, 
         number_format="%.2f", 
         fontsize_number=6.5,
         legend_breaks=c(seq(-.65,.65, by = 0.325)),
         main = "Using All Homologs")
dev.off()

#Just human identifiers. 
colrange <-  seq(-.8,.8, by = 0.01)
colorpal <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(colrange))

#Rearrange the matrix
cor_t_human_unique <- cor_t_human_unique[rownames(cor_t_human_unique)[order(rownames(cor_t_human_unique))],
                                         colnames(cor_t_human_unique)[order(colnames(cor_t_human_unique))]]
pdf(file = here("plots","Conservation","Human_Mouse_TopHumanOnly_t_correlation_HM.pdf"),
    height = 12,width = 12)
pheatmap(cor_t_human_unique,
         color=colorpal,
         cluster_cols=FALSE, 
         cluster_rows=FALSE,
         breaks=colrange,
         fontsize=11, 
         fontsize_row=11.5, 
         fontsize_col=12,
         display_numbers=T, 
         number_format="%.2f", 
         fontsize_number=6.5,
         legend_breaks=c(seq(-.8,.8, by = 0.2)),
         main = "Top 100 Genes for Human Clusters Only")
dev.off()


#Just mouse identifiers. 
colrange <-  seq(-.8,.8, by = 0.01)
colorpal <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(colrange))

#Rearrange the matrix
cor_t_mouse_unique <- cor_t_mouse_unique[rownames(cor_t_mouse_unique)[order(rownames(cor_t_mouse_unique))],
                                         colnames(cor_t_mouse_unique)[order(colnames(cor_t_mouse_unique))]]
pdf(file = here("plots","Conservation","Human_Mouse_TopMouseOnly_t_correlation_HM.pdf"),
    height = 12,width = 12)
pheatmap(cor_t_mouse_unique,
         color=colorpal,
         cluster_cols=FALSE, 
         cluster_rows=FALSE,
         breaks=colrange,
         fontsize=11, 
         fontsize_row=11.5, 
         fontsize_col=12,
         display_numbers=T, 
         number_format="%.2f", 
         fontsize_number=6.5,
         legend_breaks=c(seq(-.8,.8, by = 0.2)),
         main = "Top 100 Genes for Human Clusters Only")
dev.off()


#Shared identifiers. 
colrange <-  seq(-.9,.9, by = 0.01)
colorpal <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(colrange))

#Rearrange the matrix
cor_t_shared <- cor_t_shared[rownames(cor_t_shared)[order(rownames(cor_t_shared))],
                             colnames(cor_t_shared)[order(colnames(cor_t_shared))]]
pdf(file = here("plots","Conservation","Human_Mouse_SharedMarkers_t_correlation_HM.pdf"),
    height = 12,width = 12)
pheatmap(cor_t_shared,
         color=colorpal,
         cluster_cols=FALSE, 
         cluster_rows=FALSE,
         breaks=colrange,
         fontsize=11, 
         fontsize_row=11.5, 
         fontsize_col=12,
         display_numbers=T, 
         number_format="%.2f", 
         fontsize_number=6.5,
         legend_breaks=c(seq(-.9,.9, by = 0.45)),
         main = "Shared Markers")
dev.off()
#######################################

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# [1] "Reproducibility information:"
# [1] "2023-11-09 13:57:43 EST"
#    user   system  elapsed 
# 675.350   36.125 2598.221 
# ─ Session info ────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 Patched (2023-07-19 r84711)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2023-11-09
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/pandoc
# 
# ─ Packages ────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [2] CRAN (R 4.3.1)
# AnnotationDbi        * 1.62.2    2023-07-02 [2] Bioconductor
# Biobase              * 2.60.0    2023-04-25 [2] Bioconductor
# BiocGenerics         * 0.46.0    2023-04-25 [2] Bioconductor
# Biostrings             2.68.1    2023-05-16 [2] Bioconductor
# bit                    4.0.5     2022-11-15 [2] CRAN (R 4.3.1)
# bit64                  4.0.5     2020-08-30 [2] CRAN (R 4.3.1)
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.1)
# blob                   1.2.4     2023-03-17 [2] CRAN (R 4.3.1)
# cachem                 1.0.8     2023-05-01 [2] CRAN (R 4.3.1)
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.1)
# colorout             * 1.2-2     2023-09-22 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.3.1)
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.1)
# DBI                    1.1.3     2022-06-18 [2] CRAN (R 4.3.1)
# DelayedArray           0.26.7    2023-07-28 [2] Bioconductor
# dplyr                  1.1.3     2023-09-03 [2] CRAN (R 4.3.1)
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.3.1)
# fastmap                1.1.1     2023-02-24 [2] CRAN (R 4.3.1)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.3.1)
# GenomeInfoDb         * 1.36.3    2023-09-07 [2] Bioconductor
# GenomeInfoDbData       1.2.10    2023-07-20 [2] Bioconductor
# GenomicRanges        * 1.52.0    2023-04-25 [2] Bioconductor
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.3.1)
# gtable                 0.3.4     2023-08-21 [2] CRAN (R 4.3.1)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.3.1)
# httr                   1.4.7     2023-08-15 [2] CRAN (R 4.3.1)
# IRanges              * 2.34.1    2023-06-22 [2] Bioconductor
# KEGGREST               1.40.0    2023-04-25 [2] Bioconductor
# lattice                0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.3.1)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.3.1)
# Matrix                 1.6-1.1   2023-09-18 [3] CRAN (R 4.3.1)
# MatrixGenerics       * 1.12.3    2023-07-30 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.3.1)
# memoise                2.0.1     2021-11-26 [2] CRAN (R 4.3.1)
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.3.1)
# org.Hs.eg.db         * 3.17.0    2023-07-20 [2] Bioconductor
# org.Mm.eg.db         * 3.17.0    2023-10-26 [1] Bioconductor
# pheatmap             * 1.0.12    2019-01-04 [2] CRAN (R 4.3.1)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.3.1)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.3.1)
# png                    0.1-8     2022-11-29 [2] CRAN (R 4.3.1)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.1)
# rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 4.3.1)
# RColorBrewer         * 1.1-3     2022-04-03 [2] CRAN (R 4.3.1)
# Rcpp                   1.0.11    2023-07-06 [2] CRAN (R 4.3.1)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.1)
# rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.1)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
# RSQLite                2.3.1     2023-04-03 [2] CRAN (R 4.3.1)
# S4Arrays               1.0.6     2023-08-30 [2] Bioconductor
# S4Vectors            * 0.38.1    2023-05-02 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.3.1)
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.1)
# SingleCellExperiment * 1.22.0    2023-04-25 [2] Bioconductor
# sparseMatrixStats    * 1.12.2    2023-07-02 [2] Bioconductor
# SummarizedExperiment * 1.30.2    2023-06-06 [2] Bioconductor
# tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.3.1)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.3.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.3.1)
# vctrs                  0.6.3     2023-06-14 [2] CRAN (R 4.3.1)
# XVector                0.40.0    2023-04-25 [2] Bioconductor
# zlibbioc               1.46.0    2023-04-25 [2] Bioconductor
# 
# [1] /users/rphillip/R/4.3
# [2] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library
# 
# ───────────────────────────────────────────────────────────────────────────────