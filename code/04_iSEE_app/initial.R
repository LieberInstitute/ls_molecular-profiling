## The following list of commands will generate the plots created in iSEE
## Copy them into a script or an R session containing your SingleCellExperiment.
## All commands below refer to your SingleCellExperiment object as `se`.

se <- sce
colormap <- ExperimentColorMap()
se <- iSEE::cleanDataset(se)
colormap <- synchronizeAssays(colormap, se)
all_contents <- list()

################################################################################
# Defining brushes
################################################################################

all_active <- list()
all_active[['ReducedDimensionPlot1']] <- list()
all_active[['FeatureAssayPlot1']] <- list()

################################################################################
## Reduced dimension plot 1
################################################################################


red.dim <- reducedDim(se, "tSNE");
plot.data <- data.frame(X=red.dim[, 1], Y=red.dim[, 2], row.names=colnames(se));

plot.data$ColorBy <- colData(se)[, "CellType.Final"];
plot.data[["ColorBy"]] <- factor(plot.data[["ColorBy"]]);

# Avoid visual biases from default ordering by shuffling the points
set.seed(9225);
plot.data <- plot.data[sample(nrow(plot.data)),,drop=FALSE];

dot.plot <- ggplot() +
    geom_point(aes(x=X, y=Y, color=ColorBy), alpha=1, plot.data, size=2) +
    labs(x="Dimension 1", y="Dimension 2", color="CellType.Final", title="tSNE") +
    coord_cartesian(xlim=range(plot.data$X, na.rm=TRUE),
        ylim=range(plot.data$Y, na.rm=TRUE), expand=TRUE) +
    scale_color_manual(values=colDataColorMap(colormap, "CellType.Final", discrete=TRUE)(22), na.value='grey50', drop=FALSE) +
    scale_fill_manual(values=colDataColorMap(colormap, "CellType.Final", discrete=TRUE)(22), na.value='grey50', drop=FALSE) +
    guides(colour = guide_legend(override.aes = list(size=1)), fill = guide_legend(override.aes = list(size=1))) +
    theme_bw() +
    theme(legend.position='bottom', legend.box='vertical', legend.text=element_text(size=9), legend.title=element_text(size=11),
            axis.text=element_text(size=10), axis.title=element_text(size=12), title=element_text(size=12))

# Saving data for transmission
all_contents[['ReducedDimensionPlot1']] <- plot.data

################################################################################
## Complex heatmap 1
################################################################################


.chosen.rows <- c("AQP4", "SLC17A7", "TRPC4", "DGKG", "FREM2", "CRHR2", "OPRM1", "SCML4", "TMEM215", 
    "TMEM119", "C3", "SLC17A6", "ELAVL2", "CLDN5", "MBP", "FXYD6", "RARB", "BCL11B", 
    "DRD1", "DRD2");
.chosen.columns <- colnames(se);
plot.data <- assay(se, "logcounts")[.chosen.rows, .chosen.columns, drop=FALSE]
plot.data <- as.matrix(plot.data);

plot.data <- plot.data - rowMeans(plot.data)
plot.data <- plot.data / apply(plot.data, 1, sd)

.assay_colors <- c("blue", "white", "red")
.assay_colors <- circlize::colorRamp2(breaks = c(-3, 0, 8), colors = .assay_colors)

# Keep all samples to compute the full range of continuous annotations
.column_data <- colData(se)[, "CellType.Final", drop=FALSE]

.column_col <- list()

.column_data[["CellType.Final"]] <- factor(.column_data[["CellType.Final"]]);
.color_values <- .column_data[["CellType.Final"]]
.color_values <- setdiff(unique(.color_values), NA)
.col_colors <- colDataColorMap(colormap, "CellType.Final", discrete=TRUE)(length(.color_values))
if (is.null(names(.col_colors))) names(.col_colors) <- levels(factor(.color_values))
.column_col[["CellType.Final"]] <- .col_colors

.column_data <- .column_data[colnames(plot.data), , drop=FALSE]
.column_data <- as.data.frame(.column_data, optional=TRUE)
.column_annot_order <- order(.column_data[["CellType.Final"]])
.column_data <- .column_data[.column_annot_order, , drop=FALSE]
plot.data <- plot.data[, .column_annot_order, drop=FALSE]
.column_annot <- ComplexHeatmap::columnAnnotation(df=.column_data, col=.column_col, annotation_legend_param=list(direction="vertical"))

hm <- ComplexHeatmap::Heatmap(matrix=plot.data, col=.assay_colors,
    top_annotation=.column_annot, cluster_rows=FALSE, cluster_columns=FALSE,
    name="logcounts\n(centered, scaled)", show_row_names=TRUE,
    show_column_names=FALSE, row_names_gp=grid::gpar(fontsize=10),
    column_names_gp=grid::gpar(fontsize=10),
    heatmap_legend_param=list(direction="vertical"))

ComplexHeatmap::draw(hm, heatmap_legend_side="right", annotation_legend_side="right")

################################################################################
## Feature assay plot 1
################################################################################


plot.data <- data.frame(Y=assay(se, "logcounts")["TRPC4", ], row.names=colnames(se))
plot.data$X <- colData(se)[, "CellType.Final"];

plot.data[["X"]] <- factor(plot.data[["X"]]);

plot.data$ColorBy <- colData(se)[, "CellType.Final"];
plot.data[["ColorBy"]] <- factor(plot.data[["ColorBy"]]);

plot.data$GroupBy <- plot.data$X;
set.seed(100);
plot.data$jitteredX <- iSEE::jitterViolinPoints(plot.data$X, plot.data$Y, 
    width=0.4, varwidth=FALSE, adjust=1,
    method='quasirandom', nbins=NULL);

# Avoid visual biases from default ordering by shuffling the points
set.seed(9225);
plot.data <- plot.data[sample(nrow(plot.data)),,drop=FALSE];

dot.plot <- ggplot() +
    geom_violin(aes(x=X, y=Y, group=GroupBy), alpha=0.2, data=plot.data, scale='width', width=0.8) +
    geom_point(aes(y=Y, color=ColorBy, x=jitteredX), alpha=1, plot.data, size=1) +
    labs(x="CellType.Final", y="TRPC4 (logcounts)", color="CellType.Final", title="TRPC4 vs CellType.Final") +
    coord_cartesian(ylim=range(plot.data$Y, na.rm=TRUE), expand=TRUE) +
    scale_color_manual(values=colDataColorMap(colormap, "CellType.Final", discrete=TRUE)(22), na.value='grey50', drop=FALSE) +
    scale_fill_manual(values=colDataColorMap(colormap, "CellType.Final", discrete=TRUE)(22), na.value='grey50', drop=FALSE) +
    scale_x_discrete(drop=FALSE) +
    theme_bw() +
    theme(legend.position='bottom', legend.text=element_text(size=9),
            legend.title=element_text(size=11), legend.box='vertical',
            axis.text.x=element_text(angle=90, size=10, hjust=1, vjust=0.5),
            axis.text.y=element_text(size=10),
            axis.title=element_text(size=12), title=element_text(size=12))

# Saving data for transmission
all_contents[['FeatureAssayPlot1']] <- plot.data

################################################################################
## Row data table 1
################################################################################


tab <- as.data.frame(rowData(se));

# Saving data for transmission
all_contents[['RowDataTable1']] <- tab

################################################################################
## To guarantee the reproducibility of your code, you should also
## record the output of sessionInfo()
sessionInfo()
