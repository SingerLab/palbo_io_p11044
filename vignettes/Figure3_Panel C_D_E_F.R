# =============================================================================
# PALBO_ANTIPD1 Analysis Script
# =============================================================================

# -----------------------------------------------------------------------------
# Setup and Configuration
# -----------------------------------------------------------------------------
# Load required libraries
library(Seurat)
library(janitor)
library(ggplot2)
library(dplyr)
library(harmony)
library(future)
library(patchwork)
library(EnhancedVolcano)
library(ggpubr)
library(tidyr)
library(DESeq2)
library(readr)
library(scCustomize)

set.seed(1996)

# Load data
data(p11044x) #load("data/processed/Cohort_A/p11044x.rda")

Cancer_Cells <- subset(p11044x, curated.celltype.l1 == "Cancer Cells")

# Define color schemes
cycle.colors <- c("C1" = "#8CD2C7", "C3" = "#C51B7D")

# -----------------------------------------------------------------------------
# Harmony Analysis (P03 and P12)
# -----------------------------------------------------------------------------
set.seed(1996)
# Run harmony integration
harmonized <- subset(x = Cancer_Cells, patient.id %in% c("P03", "P12")) %>%
  RunHarmony("cycle", plot_convergence = TRUE) %>%
  RunUMAP(reduction = "harmony",
          dims = 1:35,
          seed.use = 1996) %>%
  FindNeighbors(reduction = "harmony", dims = 1:35) %>%
  FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2)) %>%
  identity()
harmonized

Idents(object = harmonized) <- "SCT_snn_res.0.4"

# Swap cluster 6 with cluster 8
current_idents <- as.character(Idents(harmonized))
new_idents <- current_idents

# Perform the swap
temp_6 <- current_idents == "6"
temp_8 <- current_idents == "8"
new_idents[temp_6] <- "8"
new_idents[temp_8] <- "6"

# Convert back to factor and assign
new_idents <- factor(new_idents, levels = levels(Idents(harmonized)))
Idents(harmonized) <- new_idents

# Also store in metadata
harmonized$SCT_snn_res.0.4_swapped <- new_idents

# Define colors based on your image (matching cluster numbers to colors)
color_vector <- c(
  "0" = "#A6CEE3",
  # Light blue
  "1" = "#1F78B4",
  # Blue
  "2" = "#B2DF8A",
  # Light green
  "3" = "#33A02C",
  # Green
  "4" = "#FB9A99",
  # Pink
  "5" = "#E31A1C",
  # Red
  "6" = "#FDBF6F",
  # Light orange (now has cells that were originally cluster 8)
  "7" = "#FF7F00",
  # Orange
  "8" = "#CAB2D6",
  # Light purple (now has cells that were originally cluster 6)
  "9" = "#6A3D9A",
  # Purple
  "10" = "#00C3FC"    # Cyan
)

# UMAP visualization
UMAP_Visualization <- DimPlot(
  harmonized,
  group.by = "SCT_snn_res.0.4_swapped",
  cols = color_vector,
  seed = 1996,
  pt.size = 0.5,
  repel = TRUE,
  label = TRUE,
  label.box = TRUE,
  label.color = "white"
) +
  theme_linedraw() +
  theme(
    axis.text = element_blank(),
    axis.title = element_text(hjust = 0),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  xlab(label = "UMAP 1") +
  ylab(label = "UMAP 2") +
  NoLegend() + scale_x_reverse()

UMAP_Visualization

# -----------------------------------------------------------------------------
# Cluster Filtering and Analysis
# -----------------------------------------------------------------------------

# Filter clusters with sufficient cells
cluster_cycle_df <- as.data.frame.matrix(table(
  harmonized@meta.data$SCT_snn_res.0.4_swapped,
  harmonized@meta.data$cycle
))
filtered_clusters <- cluster_cycle_df[cluster_cycle_df$C1 >= 50 &
                                        cluster_cycle_df$C3 >= 50, ]
print(filtered_clusters)

# Prepare filtered metadata
metadata_filtered <- harmonized@meta.data %>%
  filter(SCT_snn_res.0.4_swapped %in% rownames(filtered_clusters)) %>%
  dplyr::select(SCT_snn_res.0.4_swapped,
                bioID,
                cycle,
                REACTOME_CELL_CYCLE,
                SASP_Score) %>%
  pivot_longer(
    !c("SCT_snn_res.0.4_swapped", "bioID", "cycle"),
    names_to = "Pathway",
    values_to = "Pathway_Score"
  )

# -----------------------------------------------------------------------------
# Pathway Analysis Functions
# -----------------------------------------------------------------------------

# Function to plot median values
plot_median_values <- function(data,
                               variable_name,
                               color_scheme = cycle.colors,
                               title,
                               ylab) {
  # Aggregating median values
  Median_per_Cycle_Patient_Cluster <- data %>%
    group_by(SCT_snn_res.0.4_swapped, bioID, cycle, Pathway) %>%
    summarise(Value = median(.data[[variable_name]], na.rm = TRUE)) %>%
    pivot_wider(names_from = cycle, values_from = Value)
  
  # Plotting
  plot <- ggplot(Median_per_Cycle_Patient_Cluster) +
    geom_segment(
      aes(
        x = SCT_snn_res.0.4_swapped,
        xend = SCT_snn_res.0.4_swapped,
        y = C1,
        yend = C3
      ),
      color = "grey"
    ) +
    geom_point(
      aes(x = SCT_snn_res.0.4_swapped, y = C1),
      color = color_scheme["C1"],
      fill = color_scheme["C1"],
      size = 3,
      shape = 16
    ) +
    geom_point(
      aes(x = SCT_snn_res.0.4_swapped, y = C3),
      color = color_scheme["C3"],
      fill = color_scheme["C3"],
      size = 3,
      shape = 24
    ) +
    coord_flip() +
    theme_linedraw(base_size = 6) +
    theme(
      strip.background = element_blank(),
      axis.text.y = element_text(face = "bold"),
      strip.text = element_text(colour = "black", face = "bold"),
      legend.position = "none"
    ) +
    xlab("") +
    ggtitle(title) +
    facet_grid(bioID ~ Pathway, scales = "free_y") +
    ylab(ylab)
  
  return(plot)
}

# Generate pathway expression plot
Pathway_Expression <- plot_median_values(
  data = metadata_filtered,
  variable_name = "Pathway_Score",
  color_scheme = cycle.colors,
  title = "Median Pathway Scores per Cycle",
  ylab = "Median Pathway Score"
)

Pathway_Expression

# -----------------------------------------------------------------------------
# P12 Specific Analysis
# -----------------------------------------------------------------------------

# Subset P12 data
P12_Cancer_Cells <- subset(
  harmonized,
  patient.id == "P12" &
    SCT_snn_res.0.4_swapped %in% rownames(filtered_clusters) &
    cycle == "C3"
)

# P12_Cancer_Cells <- SCTransform(P12_Cancer_Cells, do.center = TRUE, do.scale = TRUE)
P12_Cancer_Cells@meta.data$ClusterGroup <- ifelse(P12_Cancer_Cells@meta.data$SCT_snn_res.0.4_swapped == 8,
                                                  "C8",
                                                  "OtherClusters")

print(P12_Cancer_Cells)
table(P12_Cancer_Cells$SCT_snn_res.0.4_swapped)

# DESeq2 analysis for P12
P12_Cancer_Cells.counts <- Seurat::AggregateExpression(
  object = P12_Cancer_Cells,
  assays = "RNA",
  group.by = c("SCT_snn_res.0.4_swapped")
) %>%
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

# Differential Expression Analysis

# Find markers
Markers <- FindMarkers(
  object = P12_Cancer_Cells,
  ident.1 = "8",
  assay = "SCT",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

Markers_to_export <- Markers %>%
  tibble::rownames_to_column(var = "Gene")

# write.table(
#   Markers_to_export,
#   file = "Cluster8_vs_Other.txt",
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE
# )

# Immune Genes Volcano Plot
# Define comprehensive gene list
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

# -----------------------------------------------------------------------------
# Volcano Plot
# -----------------------------------------------------------------------------

# Define colors for volcano plot
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

# Create volcano plot
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
  ylim = c(0, 65),
  xlim = c(-3, 3),
  colCustom = keyvals,
  pointSize = 1,
  drawConnectors = TRUE,
  maxoverlapsConnectors = 100,
  title = "cluster 7 vs other clusters",
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

# -----------------------------------------------------------------------------
# Pathway Signatures Analysis
# -----------------------------------------------------------------------------

# Function to assign colors based on fold change
assign_color <- function(fold_change) {
  if (is.na(fold_change)) {
    return(NA)
  } else if (fold_change >= 2) {
    return("#9E0101")
  } else if (fold_change >= 1) {
    return("#E63535")
  } else if (fold_change >= 0) {
    return("#F78A8A")
  } else if (fold_change >= -1) {
    return("#86D3FA")
  } else if (fold_change >= -2) {
    return("#07A9FA")
  } else {
    return("#1E6284")
  }
}

# Load Sarah Watson signatures
Sarah_Watson_Signatures <- read_delim(
  "DATA/Signatures/Sarah_Watson_Signatures.txt",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# Add module scores
harmonized <- AddModuleScore(
  object = harmonized,
  features = list(Sarah_Watson_Signatures$`Adipo diff`),
  name = "Adipocyte_Differenciation"
)
harmonized <- AddModuleScore(
  object = harmonized,
  features = list(Sarah_Watson_Signatures$Angiogenesis),
  name = "Angiogenesis"
)
harmonized <- AddModuleScore(
  object = harmonized,
  features = list(Sarah_Watson_Signatures$`ECM remodeling`),
  name = "ECM remodeling"
)
harmonized <- AddModuleScore(
  object = harmonized,
  features = list(Sarah_Watson_Signatures$Hypoxia),
  name = "Hypoxia"
)
harmonized <- AddModuleScore(
  object = harmonized,
  features = list(Sarah_Watson_Signatures$Invasion),
  name = "Invasion"
)
harmonized <- AddModuleScore(
  object = harmonized,
  features = list(Sarah_Watson_Signatures$Stemness),
  name = "Stemness"
)

# Update filtered metadata
metadata_filtered <- harmonized@meta.data %>%
    filter(SCT_snn_res.0.4_swapped %in% rownames(filtered_clusters))

# Define pathways to analyze
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

# Process data for all pathways
data_list <- lapply(pathways, function(pathway) {
  # Calculate median values per cluster, patient, and cycle
  aggregate(
    x = metadata_filtered[[pathway]],
    by = list(
      metadata_filtered$SCT_snn_res.0.4_swapped,
      metadata_filtered$bioID,
      metadata_filtered$cycle
    ),
    FUN = median
  ) %>%
    # Rename columns for clarity
    dplyr::rename(
      Cluster = Group.1,
      PatientID = Group.2,
      Cycle = Group.3,
      Value = x
    ) %>%
    # Reshape data so C1 and C3 are separate columns
    tidyr::pivot_wider(
      id_cols = c("Cluster", "PatientID"),
      names_from = Cycle,
      values_from = "Value"
    ) %>%
    # Remove rows with missing C1 or C3 values
    filter(!is.na(C1) & !is.na(C3)) %>%
    mutate(
      # Replace zeros with small value to avoid log2(0) errors
      C1 = round(ifelse(C1 == 0, 1e-6, C1), 4),
      C3 = round(ifelse(C3 == 0, 1e-6, C3), 4),
      # Calculate log2 fold change (C3 vs C1)
      FoldChange = round(log2(C3 / C1), 4),
      # Mark significant changes (|FC| > 1, i.e., 2-fold change)
      Significance = abs(FoldChange) > 1,
      # Add pathway name
      Pathway = pathway,
      # Assign color based on fold change direction/magnitude
      Color = sapply(FoldChange, assign_color),
      # Create unique label combining cluster and patient
      Cluster_Patient = paste0(Cluster, "_", PatientID)
    )
})

# Combine all pathway data into one dataframe
combined_data <- bind_rows(data_list) %>%
  # Assign patient-specific colors (red for P03, blue for P12)
  mutate(ClusterColor = ifelse(PatientID == "P03", "#E41A1C", "#377EB8")) %>%
  # Sort by cluster first, then patient
  arrange(Cluster, PatientID)

# Convert to factors with specific ordering
combined_data$Cluster_Patient <- factor(combined_data$Cluster_Patient,
                                        levels = unique(combined_data$Cluster_Patient))
combined_data$Pathway <- factor(combined_data$Pathway, levels = pathways)

# Extract colors for y-axis labels (one per cluster-patient combination)
label_colors <- combined_data %>%
  dplyr::select(Cluster_Patient, ClusterColor) %>%
  distinct()

# Create dotplot
DotPlot_Pathway <- ggplot(combined_data, aes(x = Pathway, y = Cluster_Patient)) +
  # Draw points: size = magnitude of FC, fill = color by FC direction
  geom_point(aes(
    size = abs(FoldChange),
    stroke = ifelse(Significance, 1, 0.5),
    fill = Color
  ), shape = 21) +
  # Scale point size
  scale_size_continuous(name = "Fold Change", range = c(1, 6)) +
  # Use the color values directly
  scale_fill_identity(name = "Fold Change") +
  # Apply clean theme
  theme_linedraw(base_size = 6) +
  theme(
    # Color y-axis labels by patient
    axis.text.y = element_text(color = label_colors$ClusterColor),
    legend.position = "bottom",
    legend.title = element_blank(),
    text = element_text(size = 6),
    # Rotate x-axis labels for readability
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  xlab("Pathway") + ylab("Cluster_Patient") +
  ggtitle("Dotplot of Fold Change Between C1 and C3")

DotPlot_Pathway

# Export pathway data
DotPlot_Pathway_Data <- DotPlot_Pathway$data %>%
  dplyr::select(Cluster, PatientID, C1, C3, FoldChange, Pathway)
head(DotPlot_Pathway_Data)

# write.table(
#   DotPlot_Pathway_Data,
#   file = "DotPlot_Pathway_Data.txt",
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE
# )

# -----------------------------------------------------------------------------
# Final Combined Visualizations
# -----------------------------------------------------------------------------

# Create combined plot
combined_plot <- UMAP_Visualization + Pathway_Expression + DotPlot_Pathway + VolcanoPlot_C8_vs_Other +
  plot_layout(widths = c(5, 4, 3.5, 4))
combined_plot
# 
# # Save combined plot
# pdf(file = "FIGURES/UMAP_Pathway.pdf",
#     width = 12,
#     height = 5.5)
# combined_plot
# dev.off()

