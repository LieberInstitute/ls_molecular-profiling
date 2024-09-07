#Goal: Compare gene expression signatures of the mouse and human LS
#Code modified from https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/51d15ef9f5f2c4c53f55e22e3fe467de1a724668/10x_NAc-n8_step04_cross-species_rnNAc_MNT.R#L4
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling/
#module load conda_R/4.4

library(SingleCellExperiment)
library(sparseMatrixStats)
library(DeconvoBuddies)
library(RColorBrewer)
library(org.Hs.eg.db) #human
library(org.Mm.eg.db) #mouse
library(sessioninfo)
library(pheatmap)
library(rafalib)
library(here)

#####load and prep data for the analysis######
######HUMAN prep
#load the SingleCellExperiment object for human Lateral Septum
load(here("processed-data","02_build_sce","sce_celltype.rda"),verbose = TRUE)

sce
# class: SingleCellExperiment 
# dim: 36601 9225 
# metadata(1): Samples
# assays(3): counts binomial_pearson_residuals logcounts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(9225): 1_AAACCCACAGCGTTGC-1 1_AAACCCACATGGCGCT-1 ...
# 3_TTTGGTTTCTTCGACC-1 3_TTTGTTGTCCCGATCT-1
# colData names(60): Sample Barcode ... k_20_louvain_1 CellType.Final
# reducedDimNames(4): GLMPCA_approx tSNE_50 mnn tSNE_mnn_50
# mainExpName: NULL
# altExpNames(0):

#Make sure everything is in correct order. 
identical(rownames(colData(sce)),colnames(sce))
#[1] TRUE

#load human DEG list from human
load(here("processed-data","markers_1vAll_ttest_k_20_louvain_CellType_Final.rda"),verbose = TRUE) 
# Loading objects:
#   markers_1vALL_df

dim(markers_1vALL_df)
#[1] 838900      9

#Split the markers_1vALL_df into a list. 
#f = splits the list by the cell type. 
markers_1vALL_list_human <- split(markers_1vALL_df,
                                  f = markers_1vALL_df$cellType.target)

#Make the rownames of each element of the list the ensemble gene id. 
#Change the rownames of the list to be the gene_id
markers_1vALL_list_human <- lapply(markers_1vALL_list_human, function(x){
  rownames(x) <- x$gene_id
  return(x)
})

#Save the markers_1vALL_list 
save(markers_1vALL_list_human,file = here("processed-data","markers_1vAll_human_list.rda"))

###### MOUSE prep
#load the SingleCellExperiment object for mouse Lateral Septum
#This object generated in 03_cross_species/mouse_LS_1vALL_DEGs.R after removing problematic clusters. 
load(file = here("processed-data","02_build_sce","mouse_sce.rda"),verbose = TRUE)
# Loading objects:
#   sce.ls

sce.ls
# class: SingleCellExperiment 
# dim: 32285 21884 
# metadata(1): Samples
# assays(3): counts binomial_pearson_residuals logcounts
# rownames(32285): ENSMUSG00000051951 ENSMUSG00000089699 ...
# ENSMUSG00000095019 ENSMUSG00000095041
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(21884): 1_AAACCCAAGGTACATA-1 1_AAACCCACATCCGAGC-1 ...
# 4_TTTGTTGCATACAGCT-1 4_TTTGTTGGTCAAACGG-1
# colData names(17): Sample Barcode ... cellType.final cellType.broad
# reducedDimNames(4): GLMPCA_approx UMAP TSNE GLMPCA_50
# mainExpName: NULL
# altExpNames(0):

identical(rownames(colData(sce.ls)),colnames(sce.ls))
#[1] TRUE


#load the 1vALL mouse DEGs
load(here("processed-data","mouse_markers_1vAll.rda"),verbose = TRUE)
# Loading objects:
#   markers_1vALL_mouse

#Split the markers_1vALL_mouse into a list
#Split list by celltype with f=
markers_1vALL_mouse_list <- split(markers_1vALL_mouse,
                                  f = markers_1vALL_mouse$cellType.target)

#Change the rownames of the list to be the gene_id
markers_1vALL_mouse_list <- lapply(markers_1vALL_mouse_list, function(x){
  rownames(x) <- x$gene_id
  return(x)
})

#save the list. 
save(markers_1vALL_mouse_list,
     file = here("processed-data","mouse_markers_1vAll_list.rda"))

#rename the objects to make more sense. 
sce_human_ls <- sce 
sce_mouse_ls <- sce.ls
rm(sce,sce.ls)


#######

#Do the objects contain any genes with no counts? 
table(rowSums(assay(sce_human_ls, "counts"))==0)
# FALSE  TRUE 
# 33556  3045 

table(rowSums(assay(sce_mouse_ls, "counts"))==0)
# FALSE  TRUE 
# 27751  4534 

#Remove genes with all 0s from the human and mouse object. 
sce_human_ls <- sce_human_ls[!rowSums(assay(sce_human_ls, "counts"))==0, ]
sce_mouse_ls <- sce_mouse_ls[!rowSums(assay(sce_mouse_ls, "counts"))==0, ]

#sanity check
table(rowSums(assay(sce_human_ls, "counts"))==0)
# FALSE 
# 33556 

table(rowSums(assay(sce_mouse_ls, "counts"))==0)
# FALSE 
# 27751 

##Add entrez gene ids and jax ids to the human object
#Entrez ids 
hs.entrezIds <- mapIds(org.Hs.eg.db, 
                       keys=rowData(sce_human_ls)$gene_id, 
                       column="ENTREZID", 
                       keytype="ENSEMBL")
#'select()' returned 1:many mapping between keys and columns

table(is.na(hs.entrezIds))
# FALSE  TRUE 
# 22855 10701 

#Which genes do not have entrez ids
withoutEntrez <- names(hs.entrezIds)[is.na(hs.entrezIds)]
names(withoutEntrez) <- rowData(sce_human_ls)[rowData(sce_human_ls)$gene_id %in% withoutEntrez, ]$gene_name
#Most of these are genes that start with AC. 

#Add entrez ids to the object
#Make sure they are in the same order. 
identical(rowData(sce_human_ls)$gene_id,names(hs.entrezIds))
#[1] TRUE

rowData(sce_human_ls) <- cbind(rowData(sce_human_ls), hs.entrezIds)

# JAX annotation info
hom <-  read.delim("http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt",
                   as.is=TRUE)

#Save dataframe with date. In case we need to use later and it gets updated. 
write.table(x = hom,
            file = here("processed-data","HOM_AllOrganism_JAX_090524.csv"),
            sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE)

#Subset for human 
hom_hs <- hom[hom$Common.Organism.Name == "human", ]

table(rowData(sce_human_ls)$hs.entrezIds %in% hom_hs$EntrezGene.ID)
# FALSE  TRUE 
# 16037 17519 

#Add the IDs to the sce_human_ls object. 
rowData(sce_human_ls)$JAX.geneID <- hom_hs$DB.Class.Key[match(rowData(sce_human_ls)$hs.entrezIds,
                                                              hom_hs$EntrezGene.ID)]

##Add entrez gene ids and jax ids to the mouse object
#Entrez ids 
mm.entrezIds <- mapIds(org.Mm.eg.db, 
                       keys=rowData(sce_mouse_ls)$gene_id, 
                       column="ENTREZID", 
                       keytype="ENSEMBL")
#'select()' returned 1:many mapping between keys and columns


table(is.na(mm.entrezIds))
# FALSE  TRUE 
# 21522  6229

#Which genes do not have entrez ids
withoutEntrez_mouse <- names(mm.entrezIds)[is.na(mm.entrezIds)]
names(withoutEntrez_mouse) <- rowData(sce_mouse_ls)[rowData(sce_mouse_ls)$gene_id %in% withoutEntrez_mouse, ]$gene_name
#Many of these are non-coding RNAs or predicted genes. 

#Add entrez ids to the object
identical(rowData(sce_mouse_ls)$gene_id,names(mm.entrezIds))
#[1] TRUE

#Add the entrez gene ids to the rowData
rowData(sce_mouse_ls) <- cbind(rowData(sce_mouse_ls), mm.entrezIds)

#Subset for mouse
hom_mm <- hom[hom$Common.Organism.Name == "mouse, laboratory", ]

#How many entrez gene IDs are also in the jax lab database. 
table(rowData(sce_mouse_ls)$mm.entrezIds %in% hom_mm$EntrezGene.ID)
# FALSE  TRUE 
# 8985 18766 

#Add the IDs to the sce_human_ls object. 
rowData(sce_mouse_ls)$JAX.geneID <- hom_mm$DB.Class.Key[match(rowData(sce_mouse_ls)$mm.entrezIds,
                                                              hom_mm$EntrezGene.ID)]

##Identify shared genes
length(intersect(rowData(sce_human_ls)$JAX.geneID,
                 rowData(sce_mouse_ls)$JAX.geneID)) 
#[1] 16575

shared_homologs <- intersect(rowData(sce_human_ls)$JAX.geneID,
                             rowData(sce_mouse_ls)$JAX.geneID)
shared_homologs <- shared_homologs[-1] #first is na

length(shared_homologs) 
#[1] 16574

# Human not in mouse
length(setdiff(rowData(sce_human_ls)$JAX.geneID,
               rowData(sce_mouse_ls)$JAX.geneID)) 
#[1] 445

# Mouse not in human
length(setdiff(rowData(sce_mouse_ls)$JAX.geneID,
               rowData(sce_human_ls)$JAX.geneID)) 
#[1] 2173

# Subset for shared homologs
sce_human_sub <- sce_human_ls[rowData(sce_human_ls)$JAX.geneID %in% shared_homologs, ]
dim(sce_human_sub)
#[1] 16999  9225

sce_mouse_sub <- sce_mouse_ls[rowData(sce_mouse_ls)$JAX.geneID %in% shared_homologs, ]
dim(sce_mouse_sub)
#[1] 16588 21884

#Are any of the JAX IDs duplicated
table(duplicated(rowData(sce_human_sub)$JAX.geneID))
# FALSE  TRUE 
# 16574   425 

table(duplicated(rowData(sce_mouse_sub)$JAX.geneID))
# FALSE  TRUE 
# 16574    14 

######Need to identify genes that are duplicated and keep the higher expressing. 
###Mouse first. 
m_dup_rows <- which(duplicated(rowData(sce_mouse_sub)$JAX.geneID))
mouse_dups <- rowData(sce_mouse_sub)[m_dup_rows,"JAX.geneID"]
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
# 16574

table(duplicated(rowData(sce_mouse_sub)$JAX.geneID))
# FALSE 
# 16574 

###human
h_dup_rows <- which(duplicated(rowData(sce_human_sub)$JAX.geneID))
human_dups <- rowData(sce_human_sub)[h_dup_rows,"JAX.geneID"]
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
# 16574 

table(duplicated(rowData(sce_human_sub)$JAX.geneID))
# FALSE 
# 16574 

#Nothing is duplicated so can move forward. 
## Match order
sce_mouse_sub <- sce_mouse_sub[match(rowData(sce_human_sub)$JAX.geneID,
                                     rowData(sce_mouse_sub)$JAX.geneID), ]

#sanity_check
identical(rowData(sce_mouse_sub)$JAX.geneID,rowData(sce_human_sub)$JAX.geneID)
#[1] TRUE

##Now make a matrix for mouse and human that contains the std.logFC for all shared homologs. 
#####MOUSE
#Get all of the std.logFC values from the list of DEGs calculated above
mouse_logFC_values <- lapply(markers_1vALL_mouse_list,function(x){
  x[rowData(sce_mouse_sub)$gene_id,"std.logFC"]
  })

#Combine all of the logFC values to make a matrix
mouse_logFC_mat <- do.call(cbind,mouse_logFC_values)
dim(mouse_logFC_mat)
#[1] 16574    33

#Add the rownames which are the ensembl geneIDs. 
rownames(mouse_logFC_mat) <- rowData(sce_mouse_sub)$gene_id

#Add the JAX.GeneID
#Make a column of ensembl gene IDs which are currently the rownames
mouse_logFC_mat <- as.data.frame(mouse_logFC_mat)
mouse_logFC_mat$gene_id <- rownames(mouse_logFC_mat)
mouse_logFC_mat <- dplyr::left_join(x = mouse_logFC_mat,
                                    y = as.data.frame(rowData(sce_mouse_sub)[,c("gene_id","JAX.geneID")]),
                                    by = "gene_id")
#Make the rownames the JAX gene ID
rownames(mouse_logFC_mat) <- mouse_logFC_mat$JAX.geneID

######HUMAN
#Get all of the std.logFC values from the list of DEGs calculated above
human_logFC_values <- lapply(markers_1vALL_list_human,function(x){
  x[rowData(sce_human_sub)$gene_id,"std.logFC"]
})

#Combine all of the logFC values to make a matrix
human_logFC_mat <- do.call(cbind,human_logFC_values)
dim(human_logFC_mat)
#[1] 16574    25

#Add the rownames which are the ensembl geneIDs. 
rownames(human_logFC_mat) <- rowData(sce_human_sub)$gene_id

#Add the JAX.GeneID
#Make a column of ensembl gene IDs which are currently the rownames
human_logFC_mat <- as.data.frame(human_logFC_mat)
human_logFC_mat$gene_id <- rownames(human_logFC_mat)
human_logFC_mat <- dplyr::left_join(x = human_logFC_mat,
                                    y = as.data.frame(rowData(sce_human_sub)[,c("gene_id","JAX.geneID")]),
                                    by = "gene_id")
rownames(human_logFC_mat) <- human_logFC_mat$JAX.geneID

#Force order of the mouse matrix to be the same as the human matrix. 
mouse_logFC_mat <- mouse_logFC_mat[match(human_logFC_mat$JAX.geneID,
                                         mouse_logFC_mat$JAX.geneID),]

#Sanity check that everything is in the same ordre. 
identical(rownames(mouse_logFC_mat),rownames(human_logFC_mat))
#[1] TRUE

identical(mouse_logFC_mat$JAX.geneID,human_logFC_mat$JAX.geneID)
#[1] TRUE

#Reorder the mouse matrix 
mouse_logFC_mat <- mouse_logFC_mat[,c("LS_In.C","LS_In.D","LS_In.M","LS_In.N","LS_In.O","LS_In.P","LS_In.Q","LS_In.R",
                                      "MS_In.J","MS_In.K","Sept_In.G","Sept_In.I","Str_In.A","Str_In.F","Str_In.H","Str_In.L",
                                      "Thal_Ex.B","TNoS_Ex.A","TT.IG.SH_Ex.C","TT.IG.SH_Ex.E","TT.IG.SH_Ex.F","Chol_Ex.D","IoC_In.E",
                                      "Astro","Ependymal","Micro","Oligo","OPC","OPC_COP","ChP","Endo","Mural","Neuroblast")]

cor_t_all <- cor(human_logFC_mat[,1:25], mouse_logFC_mat)
rownames(cor_t_all) <- paste0(rownames(cor_t_all),"_Human")
colnames(cor_t_all) <- paste0(colnames(cor_t_all),"_Mouse")
range(cor_t_all) 
#[1] -0.4235199  0.6445771

#Get top 100 genes for each human LS cluster. Top chosen by std.logFC. 
human_top_100 <- mapply(human_logFC_mat[,1:25], FUN = function(t) {
  o <- order(t, decreasing = TRUE)[1:100]
})

#Now top 100 for each mouse. 
mouse_top_100 <- mapply(mouse_logFC_mat, FUN = function(t) {
  o <- order(t, decreasing = TRUE)[1:100]
})

#get the unique identifiers for each species plus the shared. 
human_unique <- unique(as.numeric(human_top_100))
length(human_unique)
#[1] 1708

mouse_unique <- unique(as.numeric(mouse_top_100))
length(mouse_unique)
#[1] 1996

#Now find intersection of the unique identifiers for each species identified in lines 404-410
shared_identifiers <- intersect(rownames(human_logFC_mat)[human_unique], 
                                rownames(mouse_logFC_mat)[mouse_unique])
length(shared_identifiers)
#[1] 917


#Correlate with just the human identifiers. 
cor_t_human_unique <- cor(human_logFC_mat[human_unique, 1:25],
                          mouse_logFC_mat[human_unique,])
rownames(cor_t_human_unique) <- paste0(rownames(cor_t_human_unique),"_Human")
colnames(cor_t_human_unique) <- paste0(colnames(cor_t_human_unique),"_Mouse")
range(cor_t_human_unique)
#[1] -0.5074776  0.7729077

#Correlate with just the mouse identifiers. 
cor_t_mouse_unique <- cor(human_logFC_mat[mouse_unique, 1:25],
                          mouse_logFC_mat[mouse_unique,])
rownames(cor_t_mouse_unique) <- paste0(rownames(cor_t_mouse_unique),"_Human")
colnames(cor_t_mouse_unique) <- paste0(colnames(cor_t_mouse_unique),"_Mouse")
range(cor_t_mouse_unique)
#[1] -0.4961936  0.7693381

#Correlate with just the shared identifiers. 
cor_t_shared <- cor(human_logFC_mat[shared_identifiers, 1:25],
                    mouse_logFC_mat[shared_identifiers,])
rownames(cor_t_shared) <- paste0(rownames(cor_t_shared),"_Human")
colnames(cor_t_shared) <- paste0(colnames(cor_t_shared),"_Mouse")
range(cor_t_shared)
#[1] -0.5558170  0.8572425

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
#Although, mouse mat was reorder above, jsut go ahead and put everything in alphabetical order. 
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

#######
#Save the human and mouse sub sce objects. 
save(sce_mouse_sub,file = here("processed-data","sce_mouse_sub.rda"))
save(sce_human_sub,file = here("processed-data","sce_human_sub.rda"))

#######################################

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# [1] "Reproducibility information:"
# [1] "2024-09-05 14:28:45 EDT"
# ─ Session info ──────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.4.0 Patched (2024-05-22 r86590)
# os       Rocky Linux 9.4 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-09-05
# pandoc   3.1.13 @ /jhpce/shared/community/core/conda_R/4.4/bin/pandoc
# 
# ─ Packages ──────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [2] CRAN (R 4.4.0)
# AnnotationDbi        * 1.66.0    2024-05-01 [2] Bioconductor 3.19 (R 4.4.0)
# AnnotationHub          3.12.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# attempt                0.3.1     2020-05-03 [2] CRAN (R 4.4.0)
# beachmat               2.20.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# beeswarm               0.4.0     2021-06-01 [2] CRAN (R 4.4.0)
# benchmarkme            1.0.8     2022-06-12 [2] CRAN (R 4.4.0)
# benchmarkmeData        1.0.4     2020-04-23 [2] CRAN (R 4.4.0)
# Biobase              * 2.64.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# BiocFileCache          2.12.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# BiocGenerics         * 0.50.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# BiocIO                 1.14.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# BiocManager            1.30.23   2024-05-04 [2] CRAN (R 4.4.0)
# BiocNeighbors          1.22.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# BiocParallel           1.38.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# BiocSingular           1.20.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# BiocVersion            3.19.1    2024-04-17 [2] Bioconductor 3.19 (R 4.4.0)
# Biostrings             2.72.1    2024-06-02 [2] Bioconductor 3.19 (R 4.4.0)
# bit                    4.0.5     2022-11-15 [2] CRAN (R 4.4.0)
# bit64                  4.0.5     2020-08-30 [2] CRAN (R 4.4.0)
# bitops                 1.0-8     2024-07-29 [2] CRAN (R 4.4.0)
# blob                   1.2.4     2023-03-17 [2] CRAN (R 4.4.0)
# bluster                1.14.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# bslib                  0.8.0     2024-07-29 [2] CRAN (R 4.4.0)
# cachem                 1.1.0     2024-05-16 [2] CRAN (R 4.4.0)
# cli                    3.6.3     2024-06-21 [2] CRAN (R 4.4.0)
# cluster                2.1.6     2023-12-01 [3] CRAN (R 4.4.0)
# codetools              0.2-20    2024-03-31 [3] CRAN (R 4.4.0)
# colorspace             2.1-1     2024-07-26 [2] CRAN (R 4.4.0)
# config                 0.3.2     2023-08-30 [2] CRAN (R 4.4.0)
# cowplot                1.1.3     2024-01-22 [2] CRAN (R 4.4.0)
# crayon                 1.5.3     2024-06-20 [2] CRAN (R 4.4.0)
# curl                   5.2.1     2024-03-01 [2] CRAN (R 4.4.0)
# data.table             1.15.4    2024-03-30 [2] CRAN (R 4.4.0)
# DBI                    1.2.3     2024-06-02 [2] CRAN (R 4.4.0)
# dbplyr                 2.5.0     2024-03-19 [2] CRAN (R 4.4.0)
# DeconvoBuddies       * 0.99.0    2024-08-25 [1] Github (LieberInstitute/DeconvoBuddies@8ed2d3d)
# DelayedArray           0.30.1    2024-05-07 [2] Bioconductor 3.19 (R 4.4.0)
# DelayedMatrixStats     1.26.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# digest                 0.6.36    2024-06-23 [2] CRAN (R 4.4.0)
# doParallel             1.0.17    2022-02-07 [2] CRAN (R 4.4.0)
# dotCall64              1.1-1     2023-11-28 [2] CRAN (R 4.4.0)
# dplyr                  1.1.4     2023-11-17 [2] CRAN (R 4.4.0)
# dqrng                  0.4.1     2024-05-28 [2] CRAN (R 4.4.0)
# DT                     0.33      2024-04-04 [2] CRAN (R 4.4.0)
# edgeR                  4.2.1     2024-07-14 [2] Bioconductor 3.19 (R 4.4.0)
# ExperimentHub          2.12.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# fansi                  1.0.6     2023-12-08 [2] CRAN (R 4.4.0)
# fastmap                1.2.0     2024-05-15 [2] CRAN (R 4.4.0)
# fields                 16.2      2024-06-27 [2] CRAN (R 4.4.0)
# filelock               1.0.3     2023-12-11 [2] CRAN (R 4.4.0)
# foreach                1.5.2     2022-02-02 [2] CRAN (R 4.4.0)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.4.0)
# GenomeInfoDb         * 1.40.1    2024-05-24 [2] Bioconductor 3.19 (R 4.4.0)
# GenomeInfoDbData       1.2.12    2024-05-23 [2] Bioconductor
# GenomicAlignments      1.40.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# GenomicRanges        * 1.56.1    2024-06-12 [2] Bioconductor 3.19 (R 4.4.0)
# ggbeeswarm             0.7.2     2023-04-29 [2] CRAN (R 4.4.0)
# ggplot2                3.5.1     2024-04-23 [2] CRAN (R 4.4.0)
# ggrepel                0.9.5     2024-01-10 [2] CRAN (R 4.4.0)
# glue                   1.7.0     2024-01-09 [2] CRAN (R 4.4.0)
# golem                  0.4.1     2023-06-05 [2] CRAN (R 4.4.0)
# gridExtra              2.3       2017-09-09 [2] CRAN (R 4.4.0)
# gtable                 0.3.5     2024-04-22 [2] CRAN (R 4.4.0)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.4.0)
# htmltools              0.5.8.1   2024-04-04 [2] CRAN (R 4.4.0)
# htmlwidgets            1.6.4     2023-12-06 [2] CRAN (R 4.4.0)
# httpuv                 1.6.15    2024-03-26 [2] CRAN (R 4.4.0)
# httr                   1.4.7     2023-08-15 [2] CRAN (R 4.4.0)
# igraph                 2.0.3     2024-03-13 [2] CRAN (R 4.4.0)
# IRanges              * 2.38.1    2024-07-03 [2] Bioconductor 3.19 (R 4.4.0)
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.4.0)
# iterators              1.0.14    2022-02-05 [2] CRAN (R 4.4.0)
# jquerylib              0.1.4     2021-04-26 [2] CRAN (R 4.4.0)
# jsonlite               1.8.8     2023-12-04 [2] CRAN (R 4.4.0)
# KEGGREST               1.44.1    2024-06-19 [2] Bioconductor 3.19 (R 4.4.0)
# later                  1.3.2     2023-12-06 [2] CRAN (R 4.4.0)
# lattice                0.22-6    2024-03-20 [3] CRAN (R 4.4.0)
# lazyeval               0.2.2     2019-03-15 [2] CRAN (R 4.4.0)
# lifecycle              1.0.4     2023-11-07 [2] CRAN (R 4.4.0)
# limma                  3.60.4    2024-07-17 [2] Bioconductor 3.19 (R 4.4.0)
# locfit                 1.5-9.10  2024-06-24 [2] CRAN (R 4.4.0)
# magick                 2.8.4     2024-07-14 [2] CRAN (R 4.4.0)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.4.0)
# maps                   3.4.2     2023-12-15 [2] CRAN (R 4.4.0)
# Matrix                 1.7-0     2024-04-26 [3] CRAN (R 4.4.0)
# MatrixGenerics       * 1.16.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# matrixStats          * 1.3.0     2024-04-11 [2] CRAN (R 4.4.0)
# memoise                2.0.1     2021-11-26 [2] CRAN (R 4.4.0)
# metapod                1.12.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# mime                   0.12      2021-09-28 [2] CRAN (R 4.4.0)
# munsell                0.5.1     2024-04-01 [2] CRAN (R 4.4.0)
# org.Hs.eg.db         * 3.19.1    2024-05-23 [2] Bioconductor
# org.Mm.eg.db         * 3.19.1    2024-05-23 [2] Bioconductor
# paletteer              1.6.0     2024-01-21 [2] CRAN (R 4.4.0)
# pheatmap             * 1.0.12    2019-01-04 [2] CRAN (R 4.4.0)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.4.0)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.4.0)
# plotly                 4.10.4    2024-01-13 [2] CRAN (R 4.4.0)
# png                    0.1-8     2022-11-29 [2] CRAN (R 4.4.0)
# promises               1.3.0     2024-04-05 [2] CRAN (R 4.4.0)
# purrr                  1.0.2     2023-08-10 [2] CRAN (R 4.4.0)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.4.0)
# rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 4.4.0)
# rappdirs               0.3.3     2021-01-31 [2] CRAN (R 4.4.0)
# RColorBrewer         * 1.1-3     2022-04-03 [2] CRAN (R 4.4.0)
# Rcpp                   1.0.13    2024-07-17 [2] CRAN (R 4.4.0)
# RCurl                  1.98-1.16 2024-07-11 [2] CRAN (R 4.4.0)
# rematch2               2.1.2     2020-05-01 [2] CRAN (R 4.4.0)
# restfulr               0.0.15    2022-06-16 [2] CRAN (R 4.4.0)
# rjson                  0.2.21    2022-01-09 [2] CRAN (R 4.4.0)
# rlang                  1.1.4     2024-06-04 [2] CRAN (R 4.4.0)
# rprojroot              2.0.4     2023-11-05 [2] CRAN (R 4.4.0)
# Rsamtools              2.20.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# RSQLite                2.3.7     2024-05-27 [2] CRAN (R 4.4.0)
# rstudioapi             0.16.0    2024-03-24 [2] CRAN (R 4.4.0)
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.4.0)
# rtracklayer            1.64.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# S4Arrays               1.4.1     2024-05-20 [2] Bioconductor 3.19 (R 4.4.0)
# S4Vectors            * 0.42.1    2024-07-03 [2] Bioconductor 3.19 (R 4.4.0)
# sass                   0.4.9     2024-03-15 [2] CRAN (R 4.4.0)
# ScaledMatrix           1.12.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# scales                 1.3.0     2023-11-28 [2] CRAN (R 4.4.0)
# scater                 1.32.1    2024-07-21 [2] Bioconductor 3.19 (R 4.4.0)
# scran                  1.32.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# scuttle                1.14.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.4.0)
# shiny                  1.9.1     2024-08-01 [2] CRAN (R 4.4.0)
# shinyWidgets           0.8.6     2024-04-24 [2] CRAN (R 4.4.0)
# SingleCellExperiment * 1.26.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# spam                   2.10-0    2023-10-23 [2] CRAN (R 4.4.0)
# SparseArray            1.4.8     2024-05-24 [2] Bioconductor 3.19 (R 4.4.0)
# sparseMatrixStats    * 1.16.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# SpatialExperiment      1.14.0    2024-05-01 [2] Bioconductor 3.19 (R 4.4.0)
# spatialLIBD            1.16.2    2024-05-28 [2] Bioconductor 3.19 (R 4.4.0)
# statmod                1.5.0     2023-01-06 [2] CRAN (R 4.4.0)
# stringi                1.8.4     2024-05-06 [2] CRAN (R 4.4.0)
# stringr                1.5.1     2023-11-14 [2] CRAN (R 4.4.0)
# SummarizedExperiment * 1.34.0    2024-05-01 [2] Bioconductor 3.19 (R 4.4.0)
# tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.4.0)
# tidyr                  1.3.1     2024-01-24 [2] CRAN (R 4.4.0)
# tidyselect             1.2.1     2024-03-11 [2] CRAN (R 4.4.0)
# UCSC.utils             1.0.0     2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# utf8                   1.2.4     2023-10-22 [2] CRAN (R 4.4.0)
# vctrs                  0.6.5     2023-12-01 [2] CRAN (R 4.4.0)
# vipor                  0.4.7     2023-12-18 [2] CRAN (R 4.4.0)
# viridis                0.6.5     2024-01-29 [2] CRAN (R 4.4.0)
# viridisLite            0.4.2     2023-05-02 [2] CRAN (R 4.4.0)
# XML                    3.99-0.17 2024-06-25 [2] CRAN (R 4.4.0)
# xtable                 1.8-4     2019-04-21 [2] CRAN (R 4.4.0)
# XVector                0.44.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# yaml                   2.3.10    2024-07-26 [2] CRAN (R 4.4.0)
# zlibbioc               1.50.0    2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# 
# [1] /users/rphillip/R/4.4
# [2] /jhpce/shared/community/core/conda_R/4.4/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.4/R/lib64/R/library
# 
# ─────────────────────────────────────────────────────────────────────────────────────────────────
