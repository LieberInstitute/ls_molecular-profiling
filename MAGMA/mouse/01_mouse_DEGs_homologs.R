#Goal: Map mouse cell type specific DEGs to human. 
library(here)
library(SingleCellExperiment)
library(scran)
library(scater)
library(sessioninfo)

#Load in the DEG list
#Code to calculate DEGs can be found at https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/04_markerDetection.R
load(here("mouse",
          "markers-stats_LS-n4_findMarkers_33cellTypes.rda"),
     verbose = TRUE)

#DEGs generated with the cluster-vs-all method are in the markers.ls.t.1vAll object
#Add the medianNon0 information into the cluster. Code taken from the same github link. 
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
    markers.ls.t.1vAll[[i]][["1"]]$Ensembl_gene_id <- row.names(markers.ls.t.1vAll[[i]][["1"]])
}

#FDR<1e-2 and median > 0
markers <- lapply(markers.ls.t.1vAll,function(x){
    as.data.frame(subset(x[["1"]],
                         subset=(summary.stats > 0 & FDR < 1e-2 & non0Median == TRUE)))
})

#How many genes make up each set? 
as.matrix(lapply(markers,FUN = nrow))
# [,1]
# Astro         486 
# Chol_Ex.D     608 
# ChP           1232
# Endo          854 
# Ependymal     1400
# IoC_In.E      569 
# LS_In.C       1032
# LS_In.D       1633
# LS_In.M       990 
# LS_In.N       368 
# LS_In.O       515 
# LS_In.P       1386
# LS_In.Q       435 
# LS_In.R       1030
# Micro         437 
# MS_In.J       965 
# MS_In.K       1405
# Mural         239 
# Oligo         385 
# OPC           583 
# OPC_COP       527 
# Sept_In.G     720 
# Sept_In.I     1495
# Str_In.A      980 
# Str_In.F      1576
# Str_In.H      604 
# Str_In.L      400 
# Thal_Ex.B     1779
# TNoS_Ex.A     773 
# TT.IG.SH_Ex.C 1173
# TT.IG.SH_Ex.E 466 
# TT.IG.SH_Ex.F 1615
# Ventr_In.B    585 

##Using sep2019 ensembl (v98) because that was used to align the data. 
hs_mart <- biomaRt::useMart("ensembl", 
                            dataset="hsapiens_gene_ensembl",
                            host = "https://sep2019.archive.ensembl.org")
mm_mart <- biomaRt::useMart("ensembl",
                            dataset="mmusculus_gene_ensembl",
                            host = "https://sep2019.archive.ensembl.org") 

markers_hom <- lapply(markers,function(x){
    biomaRt::getLDS(attributes  = "ensembl_gene_id",
                    mart        = mm_mart,
                    values      = x$Ensembl_gene_id,
                    filters     = "ensembl_gene_id",
                    martL       = hs_mart,
                    attributesL = c("ensembl_gene_id",
                                    "external_gene_name",
                                    "entrezgene_id"))
})

##Need the entrez gene id because I will be using the gene loc files
#provided by MAGMA. 

#Rearrange the dataframes. 
for(i in names(markers_hom)){
    colnames(markers_hom[[i]]) <-  c("Mouse_ensembl","Hg38_Ensembl",
                                     "gene_symbol","entrez_id")
}

#Save this entire list. 
save(markers_hom,
     file = here("mouse","mouse_markers_human_homologs_list.rda"))

#The markers_hom list is set up so that column 1 is the mouse ensembl and column 2 is human ensembl
#For magma we just need to go ahead and have column 1 be the set or cell type.
#Column 2 needs to be human ensembl gene name. 
for(i in names(markers_hom)){
    markers_hom[[i]]$Set <- i
    markers_hom[[i]] <- markers_hom[[i]][,c("Set","entrez_id")]
    colnames(markers_hom[[i]])[2] <- "Gene"
}

#collapse the list into a single file. 
markers_hom_df <- do.call(what = rbind,markers_hom)

table(is.na(markers_hom_df$Gene))
# FALSE  TRUE 
# 27173   288 
#288 genes don't have matching entrez gene ids. 
markers_hom_df <- markers_hom_df[!is.na(markers_hom_df$Gene),]

#write out the table. 
write.table(x = markers_hom_df,
            file = here("mouse",
                        "mouse_markers_human_homologs_1e-2.txt"),
            col.names = TRUE,
            row.names = FALSE,
            quote     = FALSE,
            sep       = "\t")

Sys.time()
# [1] "2023-09-28 17:09:41 EDT"
options(width = 120)
session_info()








