initial <- list()

################################################################################
# Settings for Reduced dimension plot 1
################################################################################

initial[["ReducedDimensionPlot1"]] <- new("ReducedDimensionPlot", Type = "tSNE", XAxis = 1L, YAxis = 2L, 
    FacetRowByColData = "Sample", FacetColumnByColData = "Sample", 
    ColorByColumnData = "CellType.Final", ColorByFeatureNameAssay = "logcounts", 
    ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "Sample", 
    SizeByColumnData = "doubletScore", TooltipColumnData = character(0), 
    FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "Column data", 
    ColorByDefaultColor = "#000000", ColorByFeatureName = "LINC01409", 
    ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE, 
    ColorBySampleName = "1_AAACCCACAGCGTTGC-1", ColorBySampleSource = "---", 
    ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None", 
    SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(), 
    VisualBoxOpen = FALSE, VisualChoices = c("Color", "Size"), 
    ContourAdd = FALSE, ContourColor = "#0000FF", PointSize = 2, 
    PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200, 
    CustomLabels = FALSE, CustomLabelsText = "1_AAACCCACAGCGTTGC-1", 
    FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom", 
    HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "Sample", 
    LabelCentersColor = "#000000", VersionInfo = list(iSEE = structure(list(
        c(2L, 14L, 0L)), class = c("package_version", "numeric_version"
    ))), PanelId = 1L, PanelHeight = 600L, PanelWidth = 4L, SelectionBoxOpen = FALSE, 
    RowSelectionSource = "---", ColumnSelectionSource = "---", 
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
    SelectionHistory = list())

################################################################################
# Settings for Complex heatmap 1
################################################################################

initial[["ComplexHeatmapPlot1"]] <- new("ComplexHeatmapPlot", Assay = "logcounts", CustomRows = TRUE, 
    CustomRowsText = "AQP4\nSLC17A7\nTRPC4\nDGKG\nFREM2\nCRHR2\nOPRM1\nSCML4\nTMEM215 \nTMEM119\nC3\nSLC17A6\nELAVL2\nCLDN5\nMBP\nFXYD6\nRARB\nBCL11B \nDRD1\nDRD2", 
    ClusterRows = FALSE, ClusterRowsDistance = "spearman", ClusterRowsMethod = "ward.D2", 
    DataBoxOpen = FALSE, VisualChoices = "Annotations", ColumnData = "CellType.Final", 
    RowData = character(0), CustomBounds = FALSE, LowerBound = NA_real_, 
    UpperBound = NA_real_, AssayCenterRows = TRUE, AssayScaleRows = TRUE, 
    DivergentColormap = "blue < white < red", ShowDimNames = "Rows", 
    LegendPosition = "Right", LegendDirection = "Vertical", VisualBoxOpen = FALSE, 
    NamesRowFontSize = 10, NamesColumnFontSize = 10, ShowColumnSelection = FALSE, 
    OrderColumnSelection = FALSE, VersionInfo = list(iSEE = structure(list(
        c(2L, 14L, 0L)), class = c("package_version", "numeric_version"
    ))), PanelId = 1L, PanelHeight = 600L, PanelWidth = 8L, SelectionBoxOpen = FALSE, 
    RowSelectionSource = "---", ColumnSelectionSource = "---", 
    RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
    SelectionHistory = list())

################################################################################
# Settings for Row data table 1
################################################################################

initial[["RowDataTable1"]] <- new("RowDataTable", Selected = "MIR1302-2HG", Search = "", SearchColumns = c("", 
"", "", "", "", "", ""), HiddenColumns = character(0), VersionInfo = list(
    iSEE = structure(list(c(2L, 14L, 0L)), class = c("package_version", 
    "numeric_version"))), PanelId = 1L, PanelHeight = 500L, PanelWidth = 4L, 
    SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---", 
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
    SelectionHistory = list())

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new("FeatureAssayPlot", Assay = "logcounts", XAxis = "Column data", 
    XAxisColumnData = "CellType.Final", XAxisFeatureName = "LINC01409", 
    XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE, 
    YAxisFeatureName = "TRPC4", YAxisFeatureSource = "---", YAxisFeatureDynamicSource = FALSE, 
    FacetRowByColData = "Sample", FacetColumnByColData = "Sample", 
    ColorByColumnData = "CellType.Final", ColorByFeatureNameAssay = "logcounts", 
    ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "Sample", 
    SizeByColumnData = "doubletScore", TooltipColumnData = character(0), 
    FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "Column data", 
    ColorByDefaultColor = "#000000", ColorByFeatureName = "LINC01409", 
    ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE, 
    ColorBySampleName = "1_AAACCCACAGCGTTGC-1", ColorBySampleSource = "---", 
    ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None", 
    SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(), 
    VisualBoxOpen = FALSE, VisualChoices = "Color", ContourAdd = FALSE, 
    ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1, 
    Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE, 
    CustomLabelsText = "1_AAACCCACAGCGTTGC-1", FontSize = 1, 
    LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE, 
    LabelCenters = FALSE, LabelCentersBy = "Sample", LabelCentersColor = "#000000", 
    VersionInfo = list(iSEE = structure(list(c(2L, 14L, 0L)), class = c("package_version", 
    "numeric_version"))), PanelId = 1L, PanelHeight = 600L, PanelWidth = 8L, 
    SelectionBoxOpen = FALSE, RowSelectionSource = "---", ColumnSelectionSource = "---", 
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE, 
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE, 
    SelectionHistory = list())
