#Goal: Create heatmap from MAGMA analysis 
library(ggplot2)
library(here)
library(dplyr)
library(tibble)

#Get a table to pull cell types from
ASD <- read.table(file = here("MAGMA",
                              "mouse_analysis",
                              "magma_output",
                              "ASD.gsa.out"),
                  header = TRUE,
                  comment.char = "#")

#List to add files to
gwas_df<- as.data.frame(matrix(ncol = 8,
                               nrow = 33))
colnames(gwas_df) <- c("AD.gsa.out","ADHD.gsa.out","ASD.gsa.out","BIP.gsa.out",
                       "CUD.gsa.out","MDD.gsa.out","OUD.gsa.out","SCZ.gsa.out")
rownames(gwas_df) <- ASD$VARIABLE

for(i in colnames(gwas_df)){
    x <- read.table(file = here("MAGMA",
                                "mouse_analysis",
                                "magma_output",
                                i),
                    comment.char = "#",
                    header = TRUE)
    rownames(x) <- x$VARIABLE
    for(l in rownames(gwas_df)){
        gwas_df[l,i] <- x[l,"P"]
    }
}

#Change the dataframe so that the columns are celltype, GWAS Phenotype, and empirical p-value
dt <- gwas_df %>%
    rownames_to_column() %>%
    gather(colname, value, -rowname)


#Create a new column containing the -log10(empirical p-value)
dt$log10p <- -log10(dt$value)

#make a new column to input FDR values and stars
dt$FDR <- NA
dt$stars <- NA
#generate a column for stars and a p.adjusted column
for(i in 1:nrow(dt)){
    dt[i,"FDR"] <- p.adjust(p = dt$value[i],method = "fdr",n = 33) #33 for 33 celltypes
    if(dt[i,"FDR"]  <= 0.0001){
        dt[i,"stars"] <- "****"
    }else{
        if( dt[i,"FDR"] <=0.001){
            dt[i,"stars"]  <- "***"
        }else{
            if( dt[i,"FDR"] <=0.01){
                dt[i,"stars"]  <- "**"
            }else{
                if( dt[i,"FDR"] <=0.05){
                    dt[i,"stars"]  <- "*"
                }else{
                    dt[i,"stars"]  <- ""
                }
            }
        }
    }
}


x <- ggplot(dt,aes(x = colname, y = rowname, fill = log10p)) +
    geom_tile(color = "black",size = 0.5) +
    geom_text(aes(label = stars)) +
    scale_fill_gradient(low="white",high="red") +
    labs(x = "GWAS Phenotype",
         y = "Cell Type",
         fill = "-log10(Empircial p)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45,hjust = 1)) 

ggsave(x,filename = here("MAGMA",
                         "mouse_analysis",
                         "magma_output","mouse_hm_new.pdf"))

