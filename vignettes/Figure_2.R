## =============================================================================
#'
#' The code in this scripts generates the figures for the single-cell analysis
#'  of 13 patients with palbociclib lead-in, followed by retifanlimab.
#'
#' Code contributors: Rodrigo Gularte Merida, Evan Seffar, Li-Xuan Qin, and
#'                    Jasme Lee
#'
#' Main Figures
#' 
#' * Figure 1:
#' * Figure 2:
#' * Figure 3:
#' * Figure 4:
#'
#' Supplementary Figures
#'
#' * Supplementary Figure S1
#'
#' * Supplementary Figure S2
#'
#' Supplementary Tables
#' 
#' * Supplementary Table S1
#'
#' * Supplementary Table S2
#' 
#'
## =============================================================================
## Purposely placed this message unquoted to prevent an accidental .ess.source


## ----setup--------------------------------------------------------------------
source("vignettes/main.R")


## max amount of memmory
options(future.globals.maxSize = 24 * 1024^3)

## number of cores for multicore processing
library(future)
plan(multicore, workers = 4)

## ----output_directory_setup---------------------------------------------------
sample.name <- "01.p11044.clinical"
tabDir <- file.path("tables", sample.name)
figDir <- file.path("figures", sample.name)
tmpDir <- file.path("tmp", sample.name)
tmpImage <- "tXX.now.rda"
sapply(c(tabDir, figDir, tmpDir), usethis::use_directory)

## ----load_data----------------------------------------------------------------
## p11044x.metadata <- p11044x[[]]
## usethis::use_data(p11044x.metadata)

data(clinical.pheno.11044)

## alternatively use: data(p11044x.metadata)
if(! exists("p11044x")) {
    data(p11044x) 
}

## Median PFS and OS were estimated by Li-Xuan in separate analysis
median.pfs <- 7.1 * 30
median.os <- 26.8 * 30


## ==== Figure 2 ===============================================================
## ---- Figure 2A
## -- biopsy ID pairs : C1 / C3 distinction
umap.biopsy <- p11044x %>%
    DimPlot(group.by = "biopsy.id") +
    scale_colour_manual(name = "Cluster",
                        values = biopsy.colors) +
    theme_palbo()

## ---- Figure 2B
global.counts <- table(p11044x$curated.celltype.l1)
( global.counts <- setNames(
    paste(names(global.counts), " (n=", format(global.counts, big.mark = ","), ")", sep = ""),
    names(global.counts)) )
          

p11044x@meta.data <- p11044x@meta.data %>%
    mutate(curated.celltype.l1.globals = as.character(
               global.counts[p11044x$curated.celltype.l1]),
           curated.celltype.l1.globals = factor(curated.celltype.l1.globals,
                                                global.counts]))
           
    
dotplot.celltype.marker.genes <- DotPlot_scCustom(p11044x,
                 group.by = "curated.celltype.l1.globals",
                 features =  c("MDM2", "CDK4", "HMGA2", "CD2D",
                               "MS4A1", "CD14", "NKG7", "VWG",
                               "PDGFRB", "IGHG2", "CPA3")) +
    theme_palbo()
                                 

## ---- Figure 2C
barplot.celltype.patient.pct <- p11044x[[]] %>%
    mutate(curated.celltype.l1 = factor(curated.celltype.l1,
           c("Cancer Cells",
             "Mast",
             "Plasma Cells",
             "Macrophage",
             "B","T",
             "Endothelial",
             "NK",
             "Fibroblasts"))) %>%
    ggplot(aes(y = fct_rev(patient.id), fill = fct_rev(curated.celltype.l1))) +
    geom_bar(position = "fill") +
    scale_fill_manual(name = "Cell Type", values = cell.type.colors) +
    ylab("Cell Fraction (%)") +
    coord_cartesian(expand = FALSE) +
    theme_palbo()
barplot.celltype.patient.pct

## ---- Figure 2D.1
barplot.celltype.cycle.patient.pct <- p11044x[[]] %>%
    filter(patient.id %in% c("P03", "P12", "P13")) %>%
    mutate(curated.celltype.l1 = factor(curated.celltype.l1,
           c("Cancer Cells", "Fibroblasts", "Endothelial",
             "Plasma Cells", "Mast", 
             "NK", "B", "Macrophage", "T"))) %>%
    ggplot(aes(y = cycle, fill = fct_rev(curated.celltype.l1))) +
    geom_bar(position = "fill") +
    facet_grid(patient.id ~. , scales = "free", space = "free", switch = "y") +
    scale_fill_manual(name = "Cell Type", values = cell.type.colors) +
    xlab("Proportion of Cells (%)") +
    coord_cartesian(expand = FALSE) +
    theme_palbo() +
    theme(strip.placement = "outside",
          strip.text = element_text(face = "bold", size = 8),
          axis.title.y = element_blank())
barplot.celltype.cycle.patient.pct




## ---- Figure 2D.2
lollipop.delta <- function(df.pct, cell.type) {
    p <- df.pct %>%
        ungroup() %>%
        filter(curated.celltype.l1 %in% cell.type) %>%
        select(patient.id, cycle, cell.pct) %>%
        pivot_wider(names_from = cycle, values_from = "cell.pct") %>%
        mutate(deltaC3C1 = C3 - C1) %>%
        filter(!is.na(deltaC3C1)) %>%
        ggplot(aes(deltaC3C1, patient.id)) +
        geom_vline(xintercept = 0, lwd =  0.2) +
        geom_segment(aes(x = 0, xend = deltaC3C1,
                         y = patient.id, yend = patient.id),
                     colour = "#BDBDBD") +
        geom_point(stat = "identity", size = 1,
                   colour = cell.type.colors[cell.type]) +
        coord_cartesian(xlim = c(-14, 14),
                        ylim = c(0.5, 3.5),
                        expand = FALSE) +
        theme_palbo()

    return(p)
}

## Figure 2D.2
lollipop.delta.comprehensive.pct.c3.c1 <- df.celltype.cycle.patient.counts.pct %>%
    ungroup() %>%
    filter(curated.celltype.l1 %in% c("Cancer Cells", "Macrophage",
                                      "T", "B", "NK", "Smooth Muscle",
                                      "Endothelial")) %>%
    select(patient.id, curated.celltype.l1, cycle, cell.pct) %>%
    pivot_wider(names_from = cycle, values_from = "cell.pct") %>%
    mutate(deltaC3C1 = C3 - C1) %>%
    filter(!is.na(deltaC3C1)) %>%
    ungroup() %>%
    ggplot(aes(deltaC3C1, forcats::fct_rev(curated.celltype.l1),
               colour = curated.celltype.l1)) +
    geom_vline(xintercept = 0, lwd =  0.2) +
    geom_segment(aes(x = 0, xend = deltaC3C1,
                     y = forcats::fct_rev(curated.celltype.l1),
                     yend = forcats::fct_rev(curated.celltype.l1)),
                 colour = "#BDBDBD") +
    geom_point(stat = "identity", size = 1.6) +
    facet_grid(patient.id ~ .) +
    ylab("Cell Type") +
    xlab("Relative Change (%)") +
    coord_cartesian(xlim = c(-14, 14),
                    ylim = c(0.5, 7.5),
                    expand = FALSE)  +
    scale_colour_manual(values = cell.type.colors) +
    theme_palbo() 


pdf(file.path(figDir, "Figure_2.pdf"),
    width = 298/25.4, height = 110/25.4)
( umap.biopsy +
  theme(axis.title = element_text(hjust = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank()) |
    dotplot.celltype.marker.genes +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) |
    barplot.celltype.patient.pct |
( barplot.celltype.cycle.patient.pct | (
    lollipop.delta.comprehensive.pct.c3.c1 +
    NoLegend() +
    theme(axis.text.y =  element_blank(),          
          axis.title.x = element_text(face = "italic")) ) ) ) +
    plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")&
    theme(legend.position = "bottom")
 dev.off()




## ==== Figure 3 ===============================================================
## ----Figure_3A
## Cancer Only
t11044x <- p11044x %>%
    subset(curated.celltype.l1 == "Cancer Cells")

trt.response <- FetchData(
    t11044x,
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


trt.colors <- c("Pre IO" = "#8CD2C7",
                "Post IO" = "#C51B7D")

pd <- position_dodge(0.6)

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
    theme_palbo()

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
    theme_palbo()

summary.trt.response <- trt.response %>% summarySE(
                             measurevar = "expression",
                             groupvars = c("bor_type1", "trt", "gene"))

effect.cdk4.mdm2.p3.p12.trt.response <- summary.trt.response %>%
    filter(gene %in% c( "CDK4", "MDM2")) %>%
    ggplot(aes(trt, expression, colour = bor_type1)) +
    geom_errorbar(aes(ymin = expression-sd, ymax = expression+sd),
                  width = 0.2,
                  position = pd) +
    geom_point(position = pd) +
    geom_point(aes(y = median),
               shape = 3,
               position = pd) +
    facet_grid(~gene) +
    scale_colour_manual(values = bor.colors) +
    xlab("Biopsy Timepoint") + ylab("Expression") +
    theme_palbo()
## 
effect.signatures.p3.p12.trt.response <- summary.trt.response %>%
    filter(!gene %in% c( "CDK4", "MDM2")) %>%
    ggplot(aes(trt, expression, colour = bor_type1)) +
    geom_errorbar(aes(ymin = expression-sd, ymax = expression+sd),
                  width = 0.2,
                  position = pd) +
    geom_point(position = pd) +
    geom_point(aes(y = median),
               shape = 3,
               position = pd) +
    facet_grid(~gene) +
    scale_colour_manual(values = bor.colors) +
    xlab("Biopsy Timepoint") + ylab("Signature Score") +
    theme_palbo() 


## Figure 3A and 3B
( effect.cdk4.mdm2.p3.p12.trt.response |
    effect.signatures.p3.p12.trt.response ) +
    plot_layout(widths = c(2, 3))
(violin.cdk4.mdm2.trt.response |
    violin.signatures.trt.response ) +
    plot_layout(widths = c(2, 3)) &
    plot_annotation(tag_levels = "A")


## ---- Figure 3C
## umap P03 and P12 //

## ---- Figure 3D
paired.id.bor.colors <- c(
    P03 = "#A3D276",
    P12 = "#4C86C6",
    P13 = "#4C86C6",
    Unpaired = "#888888")
    
effect.cdk4.cycle <- gene_cycle_effect_plot_paired(t11044x, gene = "CDK4") +
    coord_cartesian(xlim = c(0.9, 2.1), expand = FALSE) +
    scale_colour_manual(values = paired.id.bor.colors) +
    theme_palbo()

effect.mdm2.cycle <- gene_cycle_effect_plot_paired(t11044x, gene = "MDM2") +
    coord_cartesian(xlim = c(0.9, 2.1), expand = FALSE) +
    scale_colour_manual(values = paired.id.bor.colors) +
    theme_palbo()

effect.cellcycle.cycle <- gene_cycle_effect_plot_paired(
    t11044x,
    gene = "REACTOME_CELL_CYCLE") +
    coord_cartesian(xlim = c(0.9, 2.1), expand = FALSE) +
    scale_colour_manual(values = paired.id.bor.colors) +
    theme_palbo()
    
effect.sasp.cycle <- gene_cycle_effect_plot_paired(
    t11044x,
    gene = "GLEASON_SASP_Score") +
    coord_cartesian(xlim = c(0.9, 2.1), expand = FALSE) +
    scale_colour_manual(values = paired.id.bor.colors) +
    theme_palbo()


pdf(file.path(figDir, "Figure_3X.pdf"),
    width = 30/25.4, height = 40/25.4)
effect.cdk4.cycle + NoLegend()
effect.cellcycle.cycle  + NoLegend()
effect.sasp.cycle + NoLegend()
effect.cdk4.cycle + theme(legend.position = "left")
dev.off()

dg <- 0.8
violin.cdk4.cycle <- gene_cycle_violin_plot_paired(t11044x, gene = "CDK4",
                                                   dodge = dg) +
    coord_cartesian(expand = FALSE) +
    scale_fill_manual(values = paired.id.bor.colors,
                      na.value = "#E5E5E5") +
    theme_palbo()
violin.cellcyle.cycle <- gene_cycle_violin_plot_paired(t11044x, gene = "REACTOME_CELL_CYCLE",
                                                   dodge = dg) +
    coord_cartesian(expand = FALSE) +
    scale_fill_manual(values = paired.id.bor.colors,
                      na.value = "#E5E5E5") +
    theme_palbo()
violin.sasp.cycle <- gene_cycle_violin_plot_paired(t11044x, gene = "GLEASON_SASP_Score",
                                                   dodge = dg) +
    coord_cartesian(expand = FALSE) +
    scale_fill_manual(values = paired.id.bor.colors,
                      na.value = "#E5E5E5") +
    theme_palbo()

##
pdf(file.path(figDir, "Figure_3X(b).pdf"),
    width = 30/25.4, height = 40/25.4)
violin.cdk4.cycle + NoLegend()
violin.cellcyle.cycle + NoLegend()
violin.sasp.cycle + NoLegend()
violin.apoptosis.cycle + NoLegend()
violin.cdk4.cycle + theme(legend.position = "left")
dev.off()

## ----Figure 3E
## Dotplot of pathway expression

## ----Figure 3F
## Volcano plot cluster 8 vs all others


## ==== Figure 4 ===============================================================
## ----Figure 4A
## T-cell barplot

## ----Figure 4B
## Flow Cytometry boxplot waterfall

## ----Figure 4C
## stacked barplot of CD4+ cells

## ----Figure 4D
## stacked barplot of CD8+ cells


