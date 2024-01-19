#Downstream analysis of MapMyCells output
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling

library(SingleCellExperiment)
library(sessioninfo)
library(pheatmap)
library(ggplot2)
library(here)

#Load human LS sce object
load(here("processed-data","sce_with_CellType.rda"),verbose = TRUE)
# Loading objects:
#     sce

#Read in the MapMyCells output
#First hierarchial mapping. 
Map_Cells_Out_hi <- read.delim(file = here("processed-data",
                                           "MapMyCells_Output",
                                           'h_ls_anndata_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1705605534631',
                                           'h_ls_anndata_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1705605534631.csv'),
                               comment.char = "#",
                               sep = ",")

hi_output <- merge(x = colData(sce)[,c("unique_rowname","CellType.Final")],
                   y = Map_Cells_Out_hi,
                   by.x = "unique_rowname",
                   by.y = "cell_id")

#Bargraph to see where each of the LS clusters are mapping to
LS_cells <- hi_output[grep("LS",hi_output$CellType.Final),]

#Create an empty dataframe and add proportion data to it via for loop
LS_cells_props <- data.frame(CellType.Final  = NA,
                             Mapped_CellType = NA,
                             Freq = NA,
                             Prop = NA)

for(i in unique(LS_cells$CellType.Final)){
    prop_df <- as.data.frame(table(subset(LS_cells,subset=(CellType.Final == i))$class_name))
    prop_df$CellType.Final <- i
    prop_df$Prop <- prop_df$Freq/sum(prop_df$Freq)*100
    colnames(prop_df)[1] <- "Mapped_CellType"
    LS_cells_props <- rbind(LS_cells_props,prop_df[,c("CellType.Final","Mapped_CellType","Freq","Prop")])
}

#remove the NA row 
LS_cells_props <- LS_cells_props[!is.na(LS_cells_props$CellType.Final),]

#Make cell colors. 
cluster_cols <- Polychrome::createPalette(length(unique(LS_cells_props$Mapped_CellType)),
c("#D81B60", "#1E88E5","#FFC107","#009E73"))
names(cluster_cols) <- unique(LS_cells_props$Mapped_CellType)

#Make the plot. 
LS_cluster_mapping <- ggplot(data = LS_cells_props,aes(x = CellType.Final,y = Prop, fill = Mapped_CellType)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cluster_cols) +
    theme_bw() +
    labs(x = "Human LS Cell Type",
         y = "Proportion of Human Cell Type",
         fill = "Mapped Mouse\nCell Type")

ggsave(filename = here("plots","Conservation","MapMyCells","Human_LS_neurons_mapped.pdf"),
       plot = LS_cluster_mapping)

#Plot a heatmap that will determine the proportion of each human cell type represnted by each 
#mapped mouse celltypes. 
#Build a matrix where the columns are the human cell types and the rows are the mapped mouse classes
mapped_mat <- matrix(ncol = length(unique(hi_output$CellType.Final)),
                     nrow = length(unique(hi_output$class_name)))

dim(mapped_mat)
#[1] 24 22

#add col and rownames
colnames(mapped_mat) <- unique(hi_output$CellType.Final)
rownames(mapped_mat) <- unique(hi_output$class_name)

#Build the matrix
for(i in colnames(mapped_mat)){
    x <- as.data.frame(table(subset(hi_output,subset=(CellType.Final == i))$class_name))
    x$Prop <- x$Freq/sum(x$Freq)*100
    rownames(x) <- x$Var1
    for(l in rownames(mapped_mat)){
        if(l %in% rownames(x)){
            mapped_mat[l,i] <- x[l,"Prop"]
        }else{
            mapped_mat[l,i] <- 0
        }
    }
}

#Plot the heatmap and save 
pdf(here("plots","Conservation","MapMyCells","Human_to_mouse_classname_heatmap.pdf"))
pheatmap(mat = mapped_mat,
         color = colorRampPalette(c("gray86", "#FFC107", "#D81B60"))(200))
dev.off()


