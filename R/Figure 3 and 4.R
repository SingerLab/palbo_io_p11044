# =============================================================================
# PALBO_ANTIPD1 Analysis Script
# =============================================================================

# -----------------------------------------------------------------------------
# Setup and Configuration
# -----------------------------------------------------------------------------
setwd("~/Downloads/PALBO_ANTIPD1/")

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

# Load data
load(
  "DATA/Seurat_Objects_Rodrigo/p11044xL.rda"
)
Cancer_Cells <- p11044xL[["Cancer Cells"]]

# Define color schemes
cycle.colors <- c("C1" = "#8CD2C7", "C3" = "#C51B7D")

color_vector <- c(
  "#A6CEE3",
  "#1F78B4",
  "#B2DF8A",
  "#33A02C",
  "#FB9A99",
  "#E31A1C",
  "#FDBF6F",
  "#FF7F00",
  "#CAB2D6",
  "#6A3D9A",
  "#00C3FC",
  "#B15928",
  "grey"
)

# -----------------------------------------------------------------------------
# Harmony Analysis (P03 and P12)
# -----------------------------------------------------------------------------

# Run harmony integration
harmonized <- subset(x = Cancer_Cells, patient.id == c("P03", "P12")) %>%
  RunHarmony("cycle", plot_convergence = TRUE) %>%
  RunUMAP(reduction = "harmony", dims = 1:35) %>%
  FindNeighbors(reduction = "harmony", dims = 1:35) %>%
  FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2)) %>%
  identity()

harmonized
Idents(object = harmonized) <- "SCT_snn_res.0.5"

# UMAP visualization
UMAP_Visualization <- DimPlot(
  harmonized,
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
  NoLegend()

UMAP_Visualization

# -----------------------------------------------------------------------------
# Cluster Filtering and Analysis
# -----------------------------------------------------------------------------

# Filter clusters with sufficient cells
cluster_cycle_df <- as.data.frame.matrix(table(
  harmonized@meta.data$SCT_snn_res.0.5,
  harmonized@meta.data$cycle
))
filtered_clusters <- cluster_cycle_df[cluster_cycle_df$C1 >= 50 &
                                        cluster_cycle_df$C3 >= 50, ]
print(filtered_clusters)

# Prepare filtered metadata
metadata_filtered <- harmonized@meta.data %>%
  filter(SCT_snn_res.0.5 %in% rownames(filtered_clusters)) %>%
  dplyr::select(SCT_snn_res.0.5, bioID, cycle, REACTOME_CELL_CYCLE, SASP_Score) %>%
  pivot_longer(!c("SCT_snn_res.0.5", "bioID", "cycle"),
               names_to = "Pathway",
               values_to = "Pathway_Score")

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
    group_by(SCT_snn_res.0.5, bioID, cycle, Pathway) %>%
    summarise(Value = median(.data[[variable_name]], na.rm = TRUE)) %>%
    pivot_wider(names_from = cycle, values_from = Value)
  
  # Plotting
  plot <- ggplot(Median_per_Cycle_Patient_Cluster) +
    geom_segment(aes(
      x = SCT_snn_res.0.5,
      xend = SCT_snn_res.0.5,
      y = C1,
      yend = C3
    ),
    color = "grey") +
    geom_point(
      aes(x = SCT_snn_res.0.5, y = C1),
      color = color_scheme["C1"],
      fill = color_scheme["C1"],
      size = 3,
      shape = 16
    ) +
    geom_point(
      aes(x = SCT_snn_res.0.5, y = C3),
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
    SCT_snn_res.0.5 %in% rownames(filtered_clusters) &
    cycle == "C3"
)
P12_Cancer_Cells <- SCTransform(P12_Cancer_Cells, do.center = TRUE, do.scale = TRUE)
P12_Cancer_Cells@meta.data$ClusterGroup <- ifelse(P12_Cancer_Cells@meta.data$SCT_snn_res.0.5 == 8,
                                                  "C8",
                                                  "OtherClusters")

print(P12_Cancer_Cells)
table(P12_Cancer_Cells$SCT_snn_res.0.5)

# DESeq2 analysis for P12
P12_Cancer_Cells.counts <- Seurat::AggregateExpression(
  object = P12_Cancer_Cells,
  assays = "RNA",
  group.by = c("SCT_snn_res.0.5")
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

# Define genes of interest
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

# -----------------------------------------------------------------------------
# Differential Expression Analysis
# -----------------------------------------------------------------------------

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

write.table(
  Markers_to_export,
  file = "DATA/DEG_Analysis/Cluster8_vs_Other_Panel_G.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

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
  filter(SCT_snn_res.0.5 %in% rownames(filtered_clusters))

# Function to plot pathways
plot_pathways <- function(cluster_column) {
  pathways <- c(
    "SASP_Score",
    "REACTOME_CELL_CYCLE",
    "HALLMARK_DNA_REPAIR",
    "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",
    "Adipocyte_Differenciation1",
    "Angiogenesis1",
    "ECM remodeling1",
    "Hypoxia1",
    "Invasion1",
    "Stemness1"
  )
  
  data_list <- lapply(pathways, function(pathway) {
    Median_per_Cycle_Patient_Cluster <- aggregate(
      x = metadata_filtered[[pathway]],
      by = list(
        metadata_filtered[[cluster_column]],
        metadata_filtered$bioID,
        metadata_filtered$cycle
      ),
      FUN = median
    ) %>%
      dplyr::rename(
        Cluster = Group.1,
        PatientID = Group.2,
        Cycle = Group.3,
        Value = x
      ) %>%
      tidyr::pivot_wider(
        id_cols = c("Cluster", "PatientID"),
        names_from = Cycle,
        values_from = "Value"
      ) %>%
      mutate(
        C1 = round(ifelse(C1 == 0, 1e-6, C1), 4),
        C3 = round(ifelse(C3 == 0, 1e-6, C3), 4),
        FoldChange = round(log2(C3 / C1), 4),
        Significance = abs(FoldChange) > 1,
        Pathway = pathway,
        Color = sapply(FoldChange, assign_color)
      )
    
    return(Median_per_Cycle_Patient_Cluster)
  })
  
  combined_data <- bind_rows(data_list)
  patient_colors <- setNames(rainbow(length(unique(
    combined_data$PatientID
  ))), unique(combined_data$PatientID))
  cluster_patient_map <- combined_data %>% dplyr::select(Cluster, PatientID) %>% unique()
  
  combined_data$Cluster <- factor(combined_data$Cluster, levels = cluster_patient_map$Cluster)
  cluster_colors <- sapply(levels(combined_data$Cluster), function(cluster) {
    patient_colors[cluster_patient_map$PatientID[cluster_patient_map$Cluster == cluster]]
  })
  
  combined_data$Pathway <- factor(x = combined_data$Pathway, levels = pathways)
  
  p <- ggplot(combined_data, aes(x = Pathway, y = Cluster)) +
    geom_point(aes(
      size = abs(FoldChange),
      stroke = Significance,
      fill = Color
    ), shape = 21) +
    scale_size_continuous(name = "Fold Change") +
    scale_fill_identity(name = "Fold Change") +
    theme_linedraw(base_size = 6) +
    theme(
      axis.text.y = element_text(color = cluster_colors[combined_data$Cluster]),
      legend.position = "bottom",
      legend.title = element_blank(),
      text = element_text(size = 6),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    xlab("Pathway") + ylab("Cluster") +
    ggtitle("Dotplot of Fold Change Between C1 and C3")
  
  return(p)
}

# Generate pathway plot
DotPlot_Pathway <- plot_pathways("SCT_snn_res.0.5")
DotPlot_Pathway

# Export pathway data
DotPlot_Pathway_Data <- DotPlot_Pathway$data %>%
  dplyr::select(Cluster, PatientID, C1, C3, FoldChange, Pathway)
head(DotPlot_Pathway_Data)

write.table(
  DotPlot_Pathway_Data,
  file = "DATA/DEG_Analysis/DotPlot_Pathway_Data_Panel_F.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Save pathway plot
pdf(file = "FIGURES/DotPlot_Pathway.pdf",
    width = 3,
    height = 5.25)
DotPlot_Pathway
dev.off()

# -----------------------------------------------------------------------------
# Additional Pathway Visualizations
# -----------------------------------------------------------------------------

# Dot plot for Sarah Watson signatures
DotPlot_scCustom(
  harmonized,
  features = c(
    "Adipocyte_Differenciation1",
    "Angiogenesis1",
    "ECM remodeling1",
    "Hypoxia1",
    "Invasion1",
    "Stemness1"
  ),
  x_lab_rotate = TRUE
)

# -----------------------------------------------------------------------------
# Final Combined Visualizations
# -----------------------------------------------------------------------------

# Create combined plot
combined_plot <- UMAP_Visualization + Pathway_Expression + DotPlot_Pathway + VolcanoPlot_C8_vs_Other +
  plot_layout(widths = c(3.75, 4, 2, 3.5))
combined_plot

# Save combined plot
pdf(file = "FIGURES/UMAP_Pathway.pdf",
    width = 11.5,
    height = 5.25)
combined_plot
dev.off()

# -----------------------------------------------------------------------------
# Memory Management
# -----------------------------------------------------------------------------

# Clean up environment (keep only essential objects)
keep_variable <- c("p11044xL", "TCells_Harmonized")
all_vars <- ls()
vars_to_remove <- setdiff(all_vars, keep_variable)
rm(list = vars_to_remove)

# Check remaining variables
print(ls())

