## =============================================================================
#'
#' The code in reanalyzes cancer and immune cells, and performs the differential
#'  expression analyses presented. Generates panels for Figure 3, and the
#'  Supplementary Tables
#'
#' Code contributors: Evan Seffar, Rodrigo Gularte Merida
#'
#' Main Figures
#' 
#' * Figure 3:
#' 
#' Supplementary Tables
#' 
#' * Supplementary Table S1
#'
#' * Supplementary Table S2
#' 
#'
## =============================================================================
Purposely placed this message unquoted to prevent an accidental .ess.source


## ----setup--------------------------------------------------------------------
source("vignettes/main.R")
## script specific libraries
library(janitor)
library(harmony)
library(future)
library(EnhancedVolcano)
library(tidyr)
library(DESeq2)
library(readr)

## ----output_directory_setup---------------------------------------------------
sample.name <- "03.p11044.Cancer_Cells"
tabDir <- file.path("tables", sample.name)
tmpDir <- file.path("tmp", sample.name)
sapply(c(tabDir, tmpDir), usethis::use_directory)

## ----load_data----------------------------------------------------------------
data(p11044x)

## ----cancer cells subset
( Cancer_Cells <- subset(p11044x, curated.celltype.l1 == "Cancer Cells") )

## ====Figure_3=================================================================
## pull data for figure
trt.response <- FetchData(
    Cancer_Cells,
    vars = c("patient.id", "cycle", "trt",
             "curated.celltype.l1",
             "bor_type1", "CDK4", "MDM2",
             "GLEASON_SASP_Score", "REACTOME_CELL_CYCLE",
             "REACTOME_CELL_CYCLE_CHECKPOINTS"
             )) %>%
    filter(curated.celltype.l1 == "Cancer Cells") %>%
    mutate(trt = recode(cycle,
                        "C1" = "Pre IO",
                        "C3" = "Post IO")) %>%
    pivot_longer(cols = c("CDK4", "MDM2", "GLEASON_SASP_Score", "REACTOME_CELL_CYCLE",
                          "REACTOME_CELL_CYCLE_CHECKPOINTS"),
                 names_to = "gene", values_to = "expression") %>%
    filter(!is.na(bor_type1)) %>%
    mutate(trt = factor(trt, c("Pre IO", "Post IO")))

## color for publication treatment groups
trt.colors <- c("Pre IO" = "#8CD2C7",
                "Post IO" = "#C51B7D")

## plot element dogdging
pd <- position_dodge(0.6)

## ----Figure 3A
## violin plot with embeded boxplot for MDM2 and CDK4
violin.cdk4.mdm2.trt.response <- trt.response %>%
    filter(gene %in% c( "CDK4", "MDM2")) %>%
    ggplot(aes(trt, expression, fill = bor_type1)) +
    geom_violin(position = pd, alpha = 0.8,
                linewidth = 0.1) +
    geom_boxplot(width = 0.1, linewidth = 0.1,
                 position = pd,
                 outlier.shape = NA) +
    scale_fill_manual(values = bor.colors) +
    facet_grid(~ gene, scales = "free_y") +
    xlab("Biopsy Timepoint") + ylab("Signature Score") +
    theme_palbo() +
    theme(strip.text = element_text(hjust = 0))

## ----Figure 3B
## violin plot with embeded boxplot for SASP, Cell Cycle, and Cell Cycle Checkpoints
violin.signatures.trt.response <- trt.response %>%
    filter(!gene %in% c( "CDK4", "MDM2")) %>%
    ggplot(aes(trt, expression, fill = bor_type1)) +
    geom_violin(position = pd, alpha = 0.8,
                linewidth = 0.1) +
    geom_boxplot(width = 0.1, linewidth = 0.1,
                 position = pd,
                 outlier.shape = NA) +
    coord_cartesian(ylim = c(-.2, 0.4)) +
    scale_fill_manual(values = bor.colors) +
    facet_grid(~ gene, scales = "free_y") +
    xlab("Biopsy Timepoint") + ylab("Signature Score") +
    theme_palbo() +
    theme(strip.text = element_text(hjust = 0))

(violin.cdk4.mdm2.trt.response |
 violin.signatures.trt.response ) +
    plot_layout(widths = c(2, 3))

## Response Summary
( summary.trt.response <- trt.response %>%
      summarySE(
          measurevar = "expression",
          groupvars = c("bor_type1", "trt", "gene")) %>%
      arrange(gene, trt) )


#' =============================================================================
#' P03 and P12 Cancer Cell Subset Reanalysis
#' P12 has one cluster with both increasing SASP and increasing Cell Cycle. We
#'  designated this as cluster 8 from our original analysis.
#' All iterations of the analysis with different seeds produce such cluster,
#'  however, the cluster identity may vary, depending on verion of R, operating
#'  system, etc., as is typical between Seurat runs.
#' We opted to manually re-label as 8 to match our text and downstream code by
#'  swaping the ID e.g. 6<>8 | 7<>8
#'
#' Both analysis code, and the `harmonized` object are provided for
#'  reproducibility


## -----------------------------------------------------------------------------
## Harmony Analysis (P03 and P12)
## -----------------------------------------------------------------------------
## Run harmony integration
set.seed(1996)
harmonized <- subset(x = Cancer_Cells, patient.id %in% c("P03", "P12")) %>%
    RunHarmony("cycle", plot_convergence = TRUE) %>%
    RunUMAP(reduction = "harmony",
            dims = 1:35,
            seed.use = 1996) %>%
    FindNeighbors(reduction = "harmony", dims = 1:35) %>%
    FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
                                0.8, 0.9, 1.0, 1.2)) %>%
    identity()
harmonized
Idents(object = harmonized) <- "SCT_snn_res.0.4"

DimPlot(harmonized)

# Swap cluster 6 with cluster 8
current_idents <- as.character(Idents(harmonized))
new_idents <- current_idents

# Perform the swap // confirm cluster ID before running    
temp_x <- current_idents == "6"
temp_8 <- current_idents == "8"
new_idents[temp_x] <- "8"
new_idents[temp_8] <- "6"

# Convert back to factor and assign
new_idents <- factor(new_idents, levels = levels(Idents(harmonized)))
Idents(harmonized) <- new_idents
# Also store in metadata
harmonized$SCT_snn_res.0.4_swapped <- new_idents
#' usethis::use_harmonized)

## =============================================================================

## ---- load harmonized data
data(harmonized)

## Light blue, Blue, Light green, Green, Pink,
## Light orange (now has cells that were originally cluster 8),
## Orange, Light purple (now has cells that were originally cluster 6),
## Purple, Cyan
color_vector <- c(
    "0" = "#A6CEE3", "1" = "#1F78B4", "2" = "#B2DF8A", "3" = "#33A02C", 
    "4" = "#FB9A99", "5" = "#E31A1C", "6" = "#FDBF6F", "7" = "#FF7F00", 
    "8" = "#CAB2D6", "9" = "#6A3D9A", "10" = "#00C3FC")

## UMAP visualization
UMAP_Visualization <- DimPlot(
    harmonized,
    cols = color_vector,
    seed = 1996,
    pt.size = 0.5,
    repel = TRUE,
    label = TRUE,
    label.box = TRUE,
    label.color = "white") +
    theme_linedraw() +
    theme(
        axis.text = element_blank(),
        axis.title = element_text(hjust = 0),
        axis.ticks = element_blank(),
        panel.grid = element_blank()
    ) +
    xlab(label = "UMAP 1") +
    ylab(label = "UMAP 2") +
    NoLegend()

UMAP_Visualization


## -----------------------------------------------------------------------------
## Cluster Filtering and Analysis
## -----------------------------------------------------------------------------

## Filter clusters with sufficient cells
cluster_cycle_df <- as.data.frame.matrix(table(
    harmonized@meta.data$SCT_snn_res.0.4_swapped,
    harmonized@meta.data$cycle
))
filtered_clusters <- cluster_cycle_df[cluster_cycle_df$C1 >= 50 &
                                      cluster_cycle_df$C3 >= 50, ]
print(filtered_clusters)

## Prepare filtered metadata
metadata_filtered <- harmonized@meta.data %>%
    filter(SCT_snn_res.0.4_swapped %in% rownames(filtered_clusters)) %>%
    dplyr::select(SCT_snn_res.0.4_swapped, bioID, cycle,
                  REACTOME_CELL_CYCLE, GLEASON_SASP_Score) %>%
    pivot_longer(!c("SCT_snn_res.0.4_swapped", "bioID", "cycle"),
                 names_to = "Pathway",
                 values_to = "Pathway_Score")

## -----------------------------------------------------------------------------
## Pathway Analysis Functions
## -----------------------------------------------------------------------------

## Generate pathway expression plot
Pathway_Expression <- plot_median_values(
    data = metadata_filtered,
    variable_name = "Pathway_Score",
    color_scheme = cycle.colors,
    title = "Median Pathway Scores per Cycle",
    ylab = "Median Pathway Score"
)

Pathway_Expression

## -----------------------------------------------------------------------------
## P12 Specific Analysis
## -----------------------------------------------------------------------------

## Subset P12 data
P12_Cancer_Cells <- subset(
    harmonized,
    patient.id == "P12" &
    SCT_snn_res.0.4 %in% rownames(filtered_clusters) &
    cycle == "C3"
)
P12_Cancer_Cells <- SCTransform(P12_Cancer_Cells,
                                do.center = TRUE, do.scale = TRUE)
                                
P12_Cancer_Cells@meta.data <- P12_Cancer_Cells@meta.data %>%
    mutate(ClusterGroup <- ifelse(SCT_snn_res.0.4_swapped == 8,
                                  "C8", "OtherClusters"))                                  

print(P12_Cancer_Cells)

P12_Cancer_Cells@meta.data %>%
    group_by(SCT_snn_res.0.4_swapped) %>%
    tally

## DESeq2 analysis for P12
P12_Cancer_Cells.counts <- AggregateExpression(
    object = P12_Cancer_Cells,
    assays = "RNA",
    group.by = "SCT_snn_res.0.4_swapped") %>%
    as.data.frame()

colData <- data.frame(
    row.names = colnames(P12_Cancer_Cells.counts),
    cluster = gsub("RNA.g", "Cluster_", colnames(P12_Cancer_Cells.counts))
)

dds <- DESeqDataSetFromMatrix(countData = P12_Cancer_Cells.counts,
                              colData = colData,
                              design = ~ cluster)

dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

## Define genes of interest
genes_of_interest <- c(
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "B2M",
    "PSMB2",
    "PSMB4",
    "PSMB8",
    "PSMB10",
    "PSMA7",
    "PSMD9",
    "PSMC1",
    "PSMC3",
    "CALR"
)

## -----------------------------------------------------------------------------
## Differential Expression Analysis
## -----------------------------------------------------------------------------

## Find markers
Markers <- FindMarkers(
    object = P12_Cancer_Cells,
    ident.1 = "8",
    assay = "SCT",
    logfc.threshold = 0.25,
    min.pct = 0.1
)

Markers_to_export <- Markers %>%
    tibble::rownames_to_column(var = "Gene")

write.table(
    Markers_to_export,
    file = file.path(tabDir, "_Cluster8_vs_Other_Panel_G.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

## Define comprehensive gene list
Genes_To_Plot <- c(
    "B2M",
    "BTK",
    "CALR",
    "CD14",
    "CD207",
    "CD36",
    "CHUK",
    "CTSL",
    "CTSS",
    "CTSV",
    "CYBA",
    "CYBB",
    "FCGR1A",
    "FCGR1BP",
    "FGA",
    "FGB",
    "FGG",
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "HLA-E",
    "HLA-F",
    "HLA-G",
    "HMGB1",
    "IKBKB",
    "IKBKG",
    "ITGAV",
    "ITGB5",
    "LNPEP",
    "LY96",
    "MRC1",
    "MRC2",
    "MYD88",
    "NCF1",
    "NCF2",
    "NCF4",
    "PDIA3",
    "PSMA1",
    "PSMA2",
    "PSMA3",
    "PSMA4",
    "PSMA5",
    "PSMA6",
    "PSMA7",
    "PSMA8",
    "PSMB1",
    "PSMB10",
    "PSMB11",
    "PSMB2",
    "PSMB3",
    "PSMB4",
    "PSMB5",
    "PSMB6",
    "PSMB7",
    "PSMB8",
    "PSMB9",
    "PSMC1",
    "PSMC2",
    "PSMC3",
    "PSMC4",
    "PSMC5",
    "PSMC6",
    "PSMD1",
    "PSMD10",
    "PSMD11",
    "PSMD12",
    "PSMD13",
    "PSMD14",
    "PSMD2",
    "PSMD3",
    "PSMD4",
    "PSMD5",
    "PSMD6",
    "PSMD7",
    "PSMD8",
    "PSMD9",
    "PSME1",
    "PSME2",
    "PSME3",
    "PSME4",
    "PSMF1",
    "RPS27A",
    "S100A1",
    "S100A8",
    "S100A9",
    "SEC22B",
    "SEC61A1",
    "SEC61A2",
    "SEC61B",
    "SEC61G",
    "SEM1",
    "SNAP23",
    "STX4",
    "TAP1",
    "TAP2",
    "TAPBP",
    "TIRAP",
    "TLR1",
    "TLR2",
    "TLR4",
    "TLR6",
    "UBA52",
    "UBB",
    "UBC",
    "VAMP3",
    "VAMP8"
)

Markers_HLA <- Markers %>%
    filter(p_val_adj <= 0.05, rownames(.) %in% Genes_To_Plot)

## -----------------------------------------------------------------------------
## Volcano Plot
## -----------------------------------------------------------------------------

## Define colors for volcano plot
keyvals <- ifelse(
    Markers$avg_log2FC <= -0.5 &
    Markers$p_val_adj <= 0.05,
    'dodgerblue',
           ifelse(
               Markers$avg_log2FC >= 0.5 &
               Markers$p_val_adj <= 0.05,
               'firebrick1',
               'azure4'
           )
)

keyvals[is.na(keyvals)] <- 'azure4'
names(keyvals)[keyvals == 'firebrick1'] <- 'UP-REGULATED'
names(keyvals)[keyvals == 'azure4'] <- 'NS'
names(keyvals)[keyvals == 'dodgerblue'] <- 'DOWN-REGULATED'

table(keyvals)

## Create volcano plot
VolcanoPlot_C8_vs_Other <- EnhancedVolcano(
    Markers,
    lab = rownames(Markers),
    x = 'avg_log2FC',
    y = 'p_val_adj',
    selectLab = c(
        "B2M",
        "HLA-A",
        "HLA-B",
        "HLA-C",
        "PSMB8",
        "PSMB4",
        "PSMB2",
        "PSMB10",
        "PSMC1",
        "PSMC3",
        "PSMA7",
        "PSMD9",
        "CALR"
    ),
    max.overlaps = 100,
    boxedLabels = TRUE,
    colCustom = keyvals,
    pointSize = 1,
    drawConnectors = TRUE,
    maxoverlapsConnectors = 100,
    title = "cluster 8 vs other clusters",
    subtitle = "",
    pCutoff = 0.05,
    FCcutoff = 0.5,
    labSize = 2,
    caption = "",
    widthConnectors = TRUE,
    legendPosition = 'bottom'
) +
    theme_linedraw(base_size = 6) +
    theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        text = element_text(size = 6)
    )

VolcanoPlot_C8_vs_Other

## -----------------------------------------------------------------------------
## Pathway Signatures Analysis
## -----------------------------------------------------------------------------
## Load Sarah Watson signatures
data(gruel.liposarcoma.signatures)

## Add module scores one by one for all to have 1 at the end
harmonized <- AddModuleScore(
    object = harmonized,
    features = gruel.liposarcoma.signatures['GRUEL_Adipo_diff'],
    name = "Adipocyte_Differenciation",
)
harmonized <- AddModuleScore(
    object = harmonized,
    features = gruel.liposarcoma.signatures['GRUEL_Angiogenesis'],
    name = "Angiogenesis"
)
harmonized <- AddModuleScore(
    object = harmonized,
    features = gruel.liposarcoma.signatures['GRUEL_ECM_remodeling'],
    name = "ECM remodeling"
)
harmonized <- AddModuleScore(
    object = harmonized,
    features = gruel.liposarcoma.signatures['GRUEL_Hypoxia'],
    name = "Hypoxia"
)
harmonized <- AddModuleScore(
    object = harmonized,
    features = gruel.liposarcoma.signatures['GRUEL_Invasion'],
    name = "Invasion"
)
harmonized <- AddModuleScore(
    object = harmonized,
    features = gruel.liposarcoma.signatures['GRUEL_Stemness'],
    name = "Stemness"
)

## Update filtered metadata
metadata_filtered <- harmonized@meta.data %>%
    filter(SCT_snn_res.0.4_swapped %in% rownames(filtered_clusters))

#' Script specific function to plot pathways // Contains an embedded data object
#' @param cluster_column metadata clster column
plot_pathways <- function(cluster_column) {
    
    pathways <- c(
        "GLEASON_SASP_Score",
        "REACTOME_CELL_CYCLE",
        "REACTOME_DNA_REPAIR",
        "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",
        "Adipocyte_Differenciation1",
        "Angiogenesis1",
        "ECM remodeling1",
        "Hypoxia1",
        "Invasion1",
        "Stemness1"
    )

    ## Process data for all pathways
    data_list <- lapply(pathways, function(pathway) {
        ## Calculate median values per cluster, patient, and cycle
        Median_per_Cycle_Patient_Cluster <- aggregate(
            x = metadata_filtered[[pathway]],
            by = list(
                metadata_filtered$SCT_snn_res.0.4_swapped,
                metadata_filtered$bioID,
                metadata_filtered$cycle
            ),
            FUN = median
        ) %>%
            ## Rename columns for clarity
            dplyr::rename(
                       Cluster = Group.1,
                       PatientID = Group.2,
                       Cycle = Group.3,
                       Value = x
                   ) %>%
            ## Reshape data so C1 and C3 are separate columns
            tidyr::pivot_wider(
                       id_cols = c("Cluster", "PatientID"),
                       names_from = Cycle,
                       values_from = "Value"
                   ) %>%
            ## Remove rows with missing C1 or C3 values
            filter(!is.na(C1) & !is.na(C3)) %>%
            mutate(
                ## Replace zeros with small value to avoid log2(0) errors
                C1 = round(ifelse(C1 == 0, 1e-6, C1), 4),
                C3 = round(ifelse(C3 == 0, 1e-6, C3), 4),
                ## Calculate log2 fold change (C3 vs C1)
                FoldChange = round(log2(C3 / C1), 4),
                ## Mark significant changes (|FC| > 1, i.e., 2-fold change)
                Significance = abs(FoldChange) > 1,
                ## Add pathway name
                Pathway = pathway,
                ## Assign color based on fold change direction/magnitude
                Color = sapply(FoldChange, assign_color),
                ## Create unique label combining cluster and patient
                Cluster_Patient = paste0(Cluster, "_", PatientID)
            )
        return(Median_per_Cycle_Patient_Cluster)
    })
    
    ## Combine all pathway data into one dataframe
    combined_data <- bind_rows(data_list) %>%
        ## Assign patient-specific colors (red for P03, blue for P12)
        mutate(ClusterColor = ifelse(PatientID == "P03", "#E41A1C", "#377EB8")) %>%
        ## Sort by cluster first, then patient
        arrange(Cluster, PatientID)

    ## Convert to factors with specific ordering
    combined_data$Cluster_Patient <- factor(combined_data$Cluster_Patient,
                                            levels = unique(combined_data$Cluster_Patient))
    combined_data$Pathway <- factor(combined_data$Pathway, levels = pathways)

    ## Extract colors for y-axis labels (one per cluster-patient combination)
    label_colors <- combined_data %>%
        dplyr::select(Cluster_Patient, ClusterColor) %>%
        distinct()

    ## Create dotplot
    p <- ggplot(combined_data, aes(x = Pathway, y = Cluster_Patient)) +
        ## Draw points: size = magnitude of FC, fill = color by FC direction
        geom_point(aes(
            size = abs(FoldChange),
            stroke = ifelse(Significance, 1, 0.5),
            fill = Color
        ), shape = 21) +
        ## Scale point size
        scale_size_continuous(name = "Fold Change", range = c(1, 6)) +
        ## Use the color values directly
        scale_fill_identity(name = "Fold Change") +
        ## Apply clean theme
        theme_linedraw(base_size = 6) +
        theme(
            ## Color y-axis labels by patient
            axis.text.y = element_text(color = label_colors$ClusterColor),
            legend.position = "bottom",
            legend.title = element_blank(),
            text = element_text(size = 6),
            ## Rotate x-axis labels for readability
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        xlab("Pathway") + ylab("Cluster_Patient") +
        ggtitle("Dotplot of Fold Change Between C1 and C3")
    
    return(p)
}

## Generate pathway plot
DotPlot_Pathway <- plot_pathways("SCT_snn_res.0.4_swapped")
DotPlot_Pathway

## Export pathway data
DotPlot_Pathway_Data <- DotPlot_Pathway$data %>%
    dplyr::select(Cluster, PatientID, C1, C3, FoldChange, Pathway)
head(DotPlot_Pathway_Data)

#' to export tables
write.table(
    DotPlot_Pathway_Data,
    file = file.path(tabDir, "_DotPlot_Pathway_Data_Panel_F.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

## Save pathway plot
## pdf(file = "DotPlot_Pathway.pdf",
##     width = 3, height = 5.25)    
## DotPlot_Pathway
## dev.off()

## -----------------------------------------------------------------------------
## Additional Pathway Visualizations
## -----------------------------------------------------------------------------

## Dot plot for Sarah Watson signatures
#' DotPlot_scCustom(
#'     harmonized,
#'     features = c(
#'         "Adipocyte_Differenciation1",
#'         "Angiogenesis1",
#'         "ECM remodeling1",
#'         "Hypoxia1",
#'         "Invasion1",
#'         "Stemness1"
#'     ),
#'     x_lab_rotate = TRUE
#' )

## -----------------------------------------------------------------------------
## Final Combined Visualizations
## -----------------------------------------------------------------------------

## Create combined plot
figure_3_top <- ( violin.cdk4.mdm2.trt.response | violin.signatures.trt.response | UMAP_Visualization ) +
    plot_layout(widths = c(2, 3, 2))

figure_3_bottom <- ( ( Pathway_Expression | DotPlot_Pathway | VolcanoPlot_C8_vs_Other ) )

figure_3_top / figure_3_bottom



## Save combined plot draft //
pdf(file = file.path(tmpDir, "Figure_3_draft.pdf"),
    width = 298/25.4, height = 298/25.4)
(figure_3_top / figure_3_bottom) +
    plot_layout(heights = c(1, 1.3))
dev.off()

## -----------------------------------------------------------------------------
## Memory Management
## -----------------------------------------------------------------------------

## Clean up environment (keep only essential objects)
keep_variable <- c("p11044x", "Cancer_Cells", "harmonized")
all_vars <- ls()
vars_to_remove <- setdiff(all_vars, keep_variable)
rm(list = vars_to_remove)


## -----------------------------------------------------------------------------
## Session Info
## -----------------------------------------------------------------------------
sessionInfo()

## __EOF__
