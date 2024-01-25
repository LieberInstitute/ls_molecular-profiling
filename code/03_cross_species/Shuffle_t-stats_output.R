#Make plots for correlation of random shuffling
#Make distributions for correlation coefficients 
library(SingleCellExperiment)
library(ggplot2)
library(here)

#Load output dataframes
load(here("processed-data","randomized_t_stats_human_500.rda"),verbose = TRUE)
# Loading objects:
#     out_h_dataframe
load(here("processed-data","randomized_t_stats_mouse_500.rda"),verbose = TRUE)
# Loading objects:
#     out_m_dataframe

#Load the human and mouse matched SCE objects.
load(here("processed-data","human_mouse_matched_by_JAX.rda"),verbose = TRUE)
# Loading objects:
#     sce_mouse_sub
#     sce_human_sub

#Add JAX.gene ID info
####Human
out_h_dataframe$gene_id <- rownames(out_h_dataframe)
out_h_dataframe <- merge(x  = out_h_dataframe,
                         y  = rowData(sce_human_sub)[,c("gene_id","JAX.geneID")],
                         by = "gene_id")
rownames(out_h_dataframe) <- out_h_dataframe$JAX.geneID
out_h_dataframe <- out_h_dataframe[,2:501]

####Mouse
out_m_dataframe$gene_id <- rownames(out_m_dataframe)
out_m_dataframe <- merge(x  = out_m_dataframe,
                         y  = rowData(sce_mouse_sub)[,c("gene_id","JAX.geneID")],
                         by = "gene_id")
rownames(out_m_dataframe) <- out_m_dataframe$JAX.geneID
out_m_dataframe <- out_m_dataframe[,2:501]

#make sure the mouse dataframe is in the same order as the human dataframe
out_m_dataframe <- out_m_dataframe[rownames(out_h_dataframe),]

#Sanity check
all(rownames(out_m_dataframe) == rownames(out_h_dataframe))
#[1] TRUE

#Make a dataframe that contains the correlation coefficients and p-values for every set. 
results_df <- as.data.frame(matrix(nrow = 500,ncol = 2))
colnames(results_df) <- c("correlation_coefficient","p-value")

#Make plots for every paired column
for(i in 1:500){
    print(i)
    #Cbind for correlation
    x <- as.data.frame(cbind(out_m_dataframe[,i],
                             out_h_dataframe[,i]))
    colnames(x) <- c("mouse","human")
    
    #Add statistics to dataframe
    results_df[i,"correlation_coefficient"] <- cor(x[,"mouse"],x[,"human"])
    results_df[i,"p-value"] <- cor.test(x[,"mouse"],x[,"human"])$p.val
    
    #Make plot
    cor_plot <- ggplot(data = x,aes(x = human,y = mouse)) +
        geom_point() +
        ggtitle(paste0("r=",
                       round(cor(x[,"mouse"],x[,"human"]),digits = 3),
                       "\np=",
                       round(cor.test(x[,"mouse"],x[,"human"])$p.val,digits = 3))) + 
        theme_bw() +
        geom_hline(yintercept = 0,lty = 2) +
        geom_vline(xintercept = 0,lty = 2) +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(plot = cor_plot,filename = here("plots","Conservation","Shuffled_cells_cor_plots",
                                           paste0(i,"_shuffled_plot.pdf")))
}

#save the results dataframe
save(results_df,file = here("processed-data","shuffle_stats_results_df.rda"))

#correlation coefficient distribution
cor_coef_distr <- ggplot(results_df,aes(x = correlation_coefficient)) +
    geom_histogram(aes(y=after_stat(density)),binwidth = 0.01) +
    geom_density(alpha=.2, fill="#FF6666") +
    labs(x = "Pearson's correlation coefficient",
         y = "Count") +
    theme_bw() +
    ggtitle("Correlation coefficient distribution") +
    theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = here("plots","Conservation","Correlation_Coefficient_distribution_random_shuffling.pdf"),plot = cor_coef_distr)


