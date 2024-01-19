#Downstream analysis of MapMyCells output
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling
#processed-data/MapMyCells_Output/'h_ls_anndata_10xWholeMouseBrain(CCN20230722)_CorrelationMapping_UTC_1705602688838'

library(SingleCellExperiment)
library(sessioninfo)
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

#Create an empty dataframe
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

LS_cells_props <- LS_cells_props[!is.na(LS_cells_props$CellType.Final),]

ggplot(data = LS_cells_props,aes(x = CellType.Final,y = Prop, fill = Mapped_CellType)) +
    geom_bar(stat = "identity")






