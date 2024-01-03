#Short script to get top 50 marker genes for each cluster. pw
#cd /dcs04/lieber/marmaypag/ls_molecular-profiling_LIBD1070/ls_molecular-profiling

#;oad librarires
library(here)
library(sessioninfo)

#Load pairwise
load(here("processed-data","markers_pairwise_list_CellTypeFinal_22CellTypes.rda"),verbose = TRUE)
# Loading objects:
#     markers_pairwise


#Create a dataframe of the top 50 genes 
pw_DEGs_df <- as.data.frame(lapply(markers_pairwise,FUN = function(x){x[1:50,"gene_name"]}))

#add pairwise to the column names
colnames(pw_DEGs_df) <- paste0(colnames(pw_DEGs_df),"_pairwise")

#Remove the periods
colnames(pw_DEGs_df) <- gsub(x = colnames(pw_DEGs_df),pattern = "[.]",replacement = "-")

#Load 1vALL 
load(here("processed-data","markers_1vAll_ttest_CellTypeFinal_22Clusters.rda"),verbose = TRUE)
# Loading objects:
#     markers_1vALL_enrich_Final

#Create an empty dataframe
OnevAll_DEGs_df <- matrix(ncol = length(unique(markers_1vALL_enrich_Final$cellType.target)),
                          nrow = 50)

#Make the column names the cell types
colnames(OnevAll_DEGs_df) <- unique(markers_1vALL_enrich_Final$cellType.target)
    
for(i in unique(markers_1vALL_enrich_Final$cellType.target)){
    print(i)
    OnevAll_DEGs_df[1:50,i] <- subset(markers_1vALL_enrich_Final,
                                      subset=(cellType.target == i & rank_marker %in% 1:50))$gene_name
}

#Add 1vall to column names
colnames(OnevAll_DEGs_df) <- paste0(colnames(OnevAll_DEGs_df),"_1vALL")

#combine the dataframes. 
DEGs <- cbind(OnevAll_DEGs_df,pw_DEGs_df)

#Order the columns by the cell names
DEGs <- DEGs[,order(colnames(DEGs))]

#Write out the DEGss. 
write.csv(x = DEGs,file = here("processed-data","Supp_Table_2.csv"))




