#Goal: Map mouse cell type specific DEGs to human. 
library(here)
library(SingleCellExperiment)
library(scran)
library(scater)


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
    subset(x[["1"]],
           subset=(summary.stats > 0 & FDR < 1e-2 & non0Median == TRUE))
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

#Now mouse genes need to be mapped to their homologs in the human genome

