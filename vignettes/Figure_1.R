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
library(ggalluvial)
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

clinical.pheno.11044 <- clinical.pheno.11044 %>%
    mutate(binary.pfs.90d = ifelse(tt_pfs >= 90, ">90 d", "<90 d"),
           binary.pfs.180d = ifelse(tt_pfs >= 180, ">180 d", "<180 d"),
           binary.pfs.median = ifelse(tt_pfs >= median.pfs, ">213 d", "<213 d"),
           pfs.z = z.value(tt_pfs)) 
           

## ---- Figure 1A skeleton
design.skeleton <- function() {
    ## base plot
    par(mar = c(1,3,1,1)) 
    plot(-1:38, rep(0, 40),
         type = "n", pch = "|", xaxs = "i", xpd = TRUE,
         xaxt = "n", yaxt = "n",         
         bty = "n", xlab = "", ylab = "",
         col = "gray50")
    ## line at 0
    ## abline(h = 0, lwd = 2)
    arrows(x0 = 0, y0 = 0, x1 = 37, y1 = 0, lwd =  2,
           length =  0.2, angle = 20)
    ## vertical arrows
    time.points <- c(0, 14, 15, 28)
    arrows(x0 = time.points, y0 = rep(0.2, length(time.points)),
           x1 = time.points, y1 = rep(0.05, length(time.points)),
           length = 0.1, angle = 20)
    points(x = time.points,
           y = rep(0, 4),
           pch = rep(21, 4),
           bg = c(rep("#E6E6E6", 2), rep("#000000", 2)),
           xpd = TRUE)
    ## text
    text.cex <- 1
    text(x = time.points[c(1,4)], y = 0.2,
         labels = c("Palbociclib",
                    "Biopsy\n(Post IO)"),
         pos = 3, xpd = TRUE,
         cex = text.cex,
         srt =  0)
    text(x = c(time.points[2]-0.8,
               time.points[3]-0.4), y = 0.23,
         labels = c("Biopsy (Pre IO)",
                    "Retifanlimab"),
         pos = 4, xpd = TRUE,
         cex = text.cex,
         srt =  45)
    text(x = time.points[-3], y = -0.2,
         labels = c(time.points[1:2], 56), xpd =  TRUE)
    text(x = time.points[1], y = -0.5,
         pos = 4,
         labels =  "Time (days)")
}
#design.skeleton()

## edited post in inkscape
pdf(file.path(figDir, "Figure_1A.pdf"),
    width = 80/25.4, height = 80/25.4)
design.skeleton()
dev.off()

## ----Figure 1B
## Draw Alluvial plot instead of consort / Trial Diagram
cohort.p11044 <- read.delim("data-raw/p11044_alluvial_data_NO_PHI.txt")

alluvial.p11044.trt.fu <- cohort.p11044 %>%
    ggplot(aes(
        axis1 = Cohort,
        axis2 = Treated.Evaluable,
        axis3 = Follow.Up,
        fill = Follow.Up)) +
    geom_alluvium(aes(fill = Follow.Up),
                  colour = "white") +
    geom_stratum() +
    scale_x_discrete(limits = c("Screened\n(n=15)",
                                "Treated & \n Evaluable\n (n=12)",
                                "Follow-Up")) +
        scale_fill_manual(values = c(
                              "POD" = "#A3D276",
                              "Surgery" = "#9970AB",
                              "Toxicity" = "#0570B0")) +
    coord_cartesian(xlim = c(0.80, 3.20), ylim = c(-0.2, 15.20),
                    expand = FALSE) +
    theme_minimal(base_size = 5) +
    theme(panel.grid.major.x = element_blank())
## alluvial.p11044.trt.fu


## ---- Figure 1C
## Waterfall plot | % change
bdx <- clinical.pheno.11044 %>%
    select(patient.id, study_id_1, bor_type1,
           recist_percent_baseline, irAE) %>%
    arrange(desc(recist_percent_baseline)) %>%
    mutate(patient.id = factor(patient.id, levels = patient.id))
##
sc.data <- p11044x[[]] %>% group_by(patient.id, cycle) %>%
    tally %>%
    pivot_wider(names_from = "cycle",
                          values_from = "n") %>%
    mutate(C1 = ifelse(is.na(C1), 0, 1),
           C3 = ifelse(is.na(C3), 0, 1),
           sc.data = paste(C1, C3, sep =  ".")) %>%
    mutate(sc.data = recode(sc.data,
                            "0.0" = "0",
                            "1.0" = "C1",
                            "0.1" = "C3",
                            "1.1" = "C1.C3"))%>%
    select(patient.id, sc.data)
##
bdx <- merge(bdx, sc.data, by = "patient.id", all.x = TRUE) %>%
    mutate(c3.loc = +85,
           c1.loc = +90,
           irae.loc = +95,
           sc.data = ifelse(is.na(sc.data), "0.0", sc.data))

##    
ylims <- ceiling(max(abs(bdx$recist_percent_baseline), na.rm = TRUE))
ylims <- c(ylims *-1, ylims)

##
waterfall.bar.p11044x <- bdx %>%
    filter(!is.na(recist_percent_baseline)) %>%           
    arrange(desc(recist_percent_baseline)) %>%
    mutate(patient.id = factor(patient.id, .$patient.id)) %>%
    ggplot(aes(patient.id, recist_percent_baseline, fill = bor_type1)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = c(-30, 30), col = "gray", lty = 3) +
    coord_cartesian(ylim = c(-100, 100), expand = FALSE) +
    scale_fill_manual(values = bor.colors, na.value = "#E6E6E6") +
    geom_point(aes(patient.id, irae.loc, colour = irAE),
               size = 1.4) +
    geom_point(aes(patient.id, c1.loc, colour = sc.data),
               size = 1.4) +
    ##    geom_point(aes(patient.id, c3.loc, colour = sc.data),
    ##               size = 1.4) +
    scale_colour_manual(values = c("0" = "#FFFFFF",
                                   "C1" = "#C6C6C6",
                                   "C3" = "#737373",
                                   "C1.C3" = "#000000",
                                   irae.colors)) +
    scale_shape_manual(values = c("C1" = "\u25D6",
                                  "C3" = "\u25D7")) +
    theme_palbo() +
    theme(axis.text.x = element_text(angle = 45, hjust =1),
          axis.title.x = element_blank(),
          legend.position = "right")


## ---- Figure 1D
## Swimers plot of Progression Free Survival
## annotated for irAE and reasion for comming off trial
tmp.swimmers.bar.pfs.out.p11044x <- clinical.pheno.11044 %>%
    arrange(desc(tt_pfs)) %>%
    mutate(patient.id = factor(patient.id, .$patient.id)) %>%
    ggplot(aes(tt_pfs, patient.id, fill = bor_type1)) +
    geom_bar(stat = "identity") +
    geom_point(aes(x = irAE.onset.from.pd1.days, shape = irAE), size = 2) +
    geom_point(aes(shape = reason_discontinuation), size = 2) +
    annotate("text", x = median.pfs, y = 12,
             label = paste0(median.pfs, "days"),
             hjust = 0) +
    coord_cartesian(xlim = c(0, 400), expand = FALSE) +
    geom_vline(xintercept = median.pfs, lwd =  0.1, lty = 3) +
    scale_fill_manual(values = bor.colors, na.value = "#E6E6E6") +
    theme_palbo() +
    xlab("PFS (days)") +
    theme(panel.grid.minor = element_blank(),
          legend.position = "right",
          axis.title.y = element_blank())

pdf(file.path(figDir, "tmp.Figure_1D.pdf"),
    width = 140/25.4, height = 50/25.4)
tmp.swimmers.bar.pfs.out.p11044x
dev.off()

## ---- Figure 1E
## Spider Plot
spd <- read.delim("data-raw/p11044_spider_plot_RECIST_NO_PHI.txt") %>%
    mutate(best.response.short = ifelse(is.na(best.response.short), "NE", best.response.short),
           best.response.short = factor(best.response.short, levels = c("PD", "SD", "PR", "NE")))

spd.patient.colors.df <- spd %>%
    distinct(patient.id, best.response.short) %>%
    mutate(bor.color = bor.colors[best.response.short]) %>%
    select(-best.response.short)


spd.patient.colors <- setNames(spd.patient.colors.df$bor.color,
                               spd.patient.colors.df$patient.id)


labels.df <- spd %>% group_by(patient.id) %>% top_n(1, fu.time.point.)

tmp.spider.line.pfs.delta.bor.p11044x <- spd %>%
    ggplot(aes(fu.days, delta.b.pct, colour = patient.id,
               shape = best.response.short)) +
    geom_hline(yintercept = c(-30, 30),
               colour = "gray", lty = 3) +
    geom_point() +
    geom_line() +
    xlab("Folloup Scan (days)") +
    ylab("Change From Baseline (%)") +
    ggrepel::geom_text_repel(data = labels.df,
              aes(x = fu.days, y = delta.b.pct, label = patient.id),
              size = 2) +
    coord_cartesian(ylim = c(-81, 56), expand = FALSE) +
    scale_colour_manual(values = spd.patient.colors) +
    theme_palbo() +
    theme(legend.position = "none")

pdf(file.path(figDir, "tmp.Figure_1E.spiderplot.pdf"),
    width = 50/25.4, height = 50/25.4)
tmp.spider.line.pfs.delta.bor.p11044x
dev.off()



## Final Figure Assembly
pdf(file.path(figDir, "tmp.Figure_1.new2.pdf"),
    width =  210/25.4, height = 60/25.4)
( ( plot_spacer()  / alluvial.p11044.trt.fu ) |
    waterfall.bar.p11044x |
    tmp.swimmers.bar.pfs.out.p11044x |
    tmp.spider.line.pfs.delta.bor.p11044x ) +
    plot_layout(widths = c(0.6, 1, 1, 1)) +
    plot_annotation(tag_levels = 'A') &
    theme(legend.position = "none")
dev.off()


## ==== Figure 2 ===============================================================
## ---- Figure 2A
## -- biopsy ID pairs : C1 / C3 distinction
umap.biopsy <- p11044x %>%
    DimPlot(group.by = "biopsy.id") +
    scale_colour_manual(name = "Cluster",
                        values = biopsy.colors) +
    theme_palbo()

pdf(file.path(figDir, "Figure_2A.pdf"),
    width = 65/25.4, height = 65/25.4)
umap.biopsy + NoLegend()
dev.off()

## ----Figure 2B
barplot.patient.x.cycle <- p11044x[[]] %>%
    ggplot(aes(y = fct_rev(cycle), group = patient.id, 
               fill = patient.id)) +
    geom_bar()  +
    stat_count(geom = "text", size = 1.4,
               aes(label = ..count..),
               position = position_stack(vjust=1),
               hjust = 0) +
    scale_fill_manual(values = patient.colors) +
    facet_grid(patient.id ~.)+
    theme_palbo() +
    scale_x_continuous(position = "top") +
    xlab("n. cells") +
    coord_cartesian(expand = FALSE) +
    theme(panel.spacing = unit(0.3, "mm"),
          axis.title.y = element_blank())

pdf(file.path(figDir, "Figure_2B.pdf"),
    width = 30/25.4, height = 50/25.4)
barplot.patient.x.cycle + NoLegend()
dev.off()

## ---- Figure 2B.1
barplot.celltype.x.cycle <- p11044x[[]] %>%
    ggplot(aes(y = fct_rev(cycle),
               group = curated.celltype.l1,
               fill = curated.celltype.l1)) +
    geom_bar()  +
    stat_count(geom = "text", size = 1.4,
               aes(label = ..count..),
               position = position_stack(vjust = 1),
               hjust = 0) +
    facet_grid(curated.celltype.l1 ~.)+
    scale_fill_manual(values = cell.type.colors) +
    theme_bw(base_size = 5) +
    scale_x_continuous(position = "top") +
    xlab("n. cells") +
    coord_cartesian(expand = FALSE) +
    theme(panel.spacing = unit(0.3, "mm"),
          axis.title.y = element_blank())

pdf(file.path(figDir, "Figure_1F.pdf"),
    width = 30/25.4, height = 50/25.4)
barplot.celltype.x.cycle + NoLegend()
dev.off()


## ---- Figure 2B.2
barplot.celltype.cycle.patient.counts <- p11044x[[]] %>%
    ggplot(aes(x = cycle, fill = curated.celltype.l1)) +
    geom_bar(position = "dodge", preserve = "single") +
    facet_grid(~ patient.id, scales = "free_x", space = "free") +
    scale_y_continuous(breaks = c(1, 10, 100, 1000, 10000),
                       transform = "log10") +
    coord_cartesian(expand = FALSE) +
    scale_fill_manual(name = "Cell Type", values = cell.type.colors) +
    ylab("Cell Counts") +
    theme_palbo()

## ---- Figure 2B.3
df.celltype.cycle.patient.counts.pct <- p11044x[[]] %>%
    group_by(sample.id, curated.celltype.l1) %>%
    tally %>%
    ungroup() %>%
    pivot_wider(names_from = "curated.celltype.l1", values_from = "n")  %>%
    mutate(across(everything(), ~replace_na(., 0))) %>%
    pivot_longer(cols = levels(p11044x$curated.celltype.l1),
                 names_to = "curated.celltype.l1",
                 values_to = "n") %>%
    group_by(sample.id) %>%
    mutate(cell.pct = n/sum(n)*100,
           patient.id = gsub("_C[0-9]", "", sample.id),
           cycle = gsub("P[0-9]+_", "", sample.id),
           cell.pct1 = ifelse(cell.pct <= 1, 0, cell.pct),
           curated.celltype.l1 = factor(curated.celltype.l1,
                                        levels(p11044x$curated.celltype.l1)))


barplot.celltype.cycle.patient.counts.pct <- df.celltype.cycle.patient.counts.pct %>%
    ggplot(aes(x = cycle, y = cell.pct1, fill = curated.celltype.l1,
               label = round(cell.pct, 2))) +    
    geom_bar(colour = "black", linewidth =  0.1,
             position = "dodge", stat = "identity") +
    facet_grid(~ patient.id, scales = "free_x", space = "free") +
    coord_cartesian(ylim = c(1, 100), expand = FALSE) +
    scale_y_continuous(breaks = c(0, 2, 4, 8, 16, 32, 64, 100),
                       transform = "log2") +
    scale_fill_manual(name = "Cell Type", values = cell.type.colors) +
    ylab("Cell Fraction (%)") + xlab("Cell Type") +
    theme_palbo() 

pdf(file.path(tmpDir, "Tumor_Cell_Composition_Side_by_Side.pdf"),
    width = 210/25.4, height = 60/25.4)
barplot.celltype.cycle.patient.counts
barplot.celltype.cycle.patient.counts.pct
dev.off()

## ---- Figure 2B.4
barplot.celltype.cycle.patient.pct <- p11044x[[]] %>%
    mutate(curated.celltype.l1 = factor(curated.celltype.l1,
           c("Cancer Cells", "Smooth Muscle", "Endothelial",
             "Plasma Cells", "Mast", 
             "NK", "B", "Macrophage", "T"))) %>%
    ggplot(aes(x = cycle, fill = fct_rev(curated.celltype.l1))) +
    geom_bar(position = "fill") +
    facet_grid(~ patient.id, scales = "free_x", space = "free") +
    scale_fill_manual(name = "Cell Type", values = cell.type.colors) +
    ylab("Cell Fraction (%)") +
    coord_cartesian(expand = FALSE) +
    theme_palbo()

pdf(file.path(tmpDir, "Tumor_Cell_Composition_Stacked.pdf"),
    width = 105/25.4, height = 105/25.4)
barplot.celltype.cycle.patient.pct
dev.off()

barplot.selected.celltype.cycle.patient.pct <- p11044x[[]] %>%
    mutate(
        curated.celltype.l1 = ifelse(
            as.character(curated.celltype.l1) %in% c("Smooth Muscle", "Endothelial",
                                       "Plasma Cells"), "Other",
            curated.celltype.l1),
        curated.celltype.l1 = factor(curated.celltype.l1,
           c("Cancer Cells", "Other",
             "Mast", 
             "NK", "B", "Macrophage", "T"))) %>%           
    ggplot(aes(x = cycle, fill = fct_rev(curated.celltype.l1))) +
    geom_bar(position = "fill") +
    facet_grid(~ patient.id, scales = "free_x", space = "free") +
    scale_fill_manual(name = "Cell Type", values = cell.type.colors) +
    ylab("Cell Fraction (%)") +
    coord_cartesian(expand = FALSE) +
    theme_palbo()

## ---- Figure 2C
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

lollipop.delta.cancer.pct.c3.c1 <- df.celltype.cycle.patient.counts.pct %>%
    lollipop.delta(cell.type = "Cancer Cells")

lollipop.delta.macrophage.pct.c3.c1 <- df.celltype.cycle.patient.counts.pct %>%
    lollipop.delta(cell.type = "Macrophage")

lollipop.delta.tcell.pct.c3.c1 <- df.celltype.cycle.patient.counts.pct %>%
    lollipop.delta(cell.type = "T")

lollipop.delta.bcell.pct.c3.c1 <- df.celltype.cycle.patient.counts.pct %>%
    lollipop.delta(cell.type = "B")

lollipop.delta.fb.pct.c3.c1 <- df.celltype.cycle.patient.counts.pct %>%
    lollipop.delta(cell.type = "Smooth Muscle")

lollipop.delta.endo.pct.c3.c1 <- df.celltype.cycle.patient.counts.pct %>%
    lollipop.delta(cell.type = "Endothelial")


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

pdf(file.path(figDir, "Figure_1C_deltaC1C3.pdf"),
    width = 30/25.4, height = 45/25.4)
lollipop.delta.comprehensive.pct.c3.c1 +
    NoLegend() +
    theme(axis.text.y =  element_blank(),
          axis.title.x = element_text(face = "italic"))
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
pdf(file.path(figDir, "Figure_3A_3B.pdf"),
    width = 120/25.4, 50/25.4)
( effect.cdk4.mdm2.p3.p12.trt.response |
    effect.signatures.p3.p12.trt.response ) +
    plot_layout(widths = c(2, 3))
(violin.cdk4.mdm2.trt.response |
    violin.signatures.trt.response ) +
    plot_layout(widths = c(2, 3))
dev.off()

1

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


## =============================================================================
##
## SUPPLEMENTARY FIGURES
##
## =============================================================================


## ----Supplementary Figure 1
## Note: Survival analysis and manuscript figure were provided by Biostats
## the code here is an example in an effort to replicate their analysis.
## Some deviations in censoring were noticed.

## Kaplan-Meier curve
library(ggfortify)
library(survival)
library(survminer)

fit.pfs <- survfit(Surv((tt_pfs), event = pt_dod) ~ active_fu,
                   data = clinical.pheno.11044)
km.pfs.plot.trial <- autoplot(fit.pfs) +
    coord_cartesian(xlim = c(0, 392), expand = FALSE) +
    xlab("Time (days)") + ylab("Survival (%)") +
    geom_vline(xintercept = median.pfs, lwd = 0.1, lty = 3) +
    annotate("text", x = median.pfs, y =  .95, label = paste(median.pfs, "days"),
             size = 2, hjust = 0) +
    theme_palbo()
km.pfs.plot.trial
## 
fit.os <- survfit(Surv((tt_os), event = pt_dod) ~ active_fu, data = clinical.pheno.11044)
km.os.plot.trial <- autoplot(fit.os) +
    coord_cartesian(expand = FALSE) +
    xlab("Time (days)") + ylab("Survival (%)") +
    geom_vline(xintercept = median.os, lwd = 0.1, lty = 3) +
    annotate("text", x = median.os, y =  .95, label = paste(median.os, "days"),
             size = 2, hjust = 0) +
    theme_palbo()
km.os.plot.trial


## 
pdf(file.path(figDir, "tmp.km.supplemental.pdf"),
    width = 120/25.4, height = 55/25.4)
km.pfs.plot.trial | km.os.plot.trial
dev.off()




## ---- Supplemenary Figure S2 --
## Gruel Signatures
## data(t11044x)
data(gruel.liposarcoma.signatures)

p11044x %>%
    FeaturePlot_scCustom(
        features = names(gruel.liposarcoma.signatures),
        na_cutoff = 0.3)


t11044x$GRUEL_Signature_Classification <- grep(
    "GRUEL", names(t11044x[[]]),
    value = TRUE)[
    apply(t11044x[[]][, grepl("GRUEL", names(t11044x[[]]))], 1, which.max)]


cluster.colors <- pals::brewer.paired(16)
names(cluster.colors) <- 0:15
##
umap.tumor.clusters <- t11044x %>%
    DimPlot(group.by = "seurat_clusters",
            pt.size = 2,
            raster = TRUE,
            label = TRUE, repel = TRUE,
            label.size = 2, label.box = FALSE) +            
    scale_colour_manual(values = cluster.colors) +
    scale_fill_manual(values = cluster.colors) +
    guides(colour = guide_legend(nrow = 2,
                                 override.aes = list(size = 1))) +
    labs(title = "Seurat Clusters") +
    theme_palbo() +
    theme(axis.text = element_blank(),
          axis.title = element_text(hjust = 0))

umap.tumor.biopsy <- t11044x %>%    
    DimPlot(group.by = "biopsy.id",
            pt.size = 2,
            raster = TRUE,
            label = FALSE, repel = TRUE,
            label.size = 2, label.box = FALSE) +    
    scale_colour_manual(values = biopsy.colors) +
    guides(colour = guide_legend(nrow = 2,
                                 override.aes = list(size = 2))) +
    labs(title = "Biopsy ID") +
    theme_palbo() +
    theme(axis.text = element_blank(),
          axis.title = element_text(hjust = 0))

umap.tumor.gruel.signature.class <- t11044x %>%
    DimPlot(group.by = "GRUEL_Signature_Classification",
            pt.size = 2,
            raster = TRUE) +            
    scale_colour_manual(values =  gruel.signature.colors) +
    labs(title = "Signature Classification") +
    theme_palbo() +
    theme(axis.text = element_blank(),
          axis.title = element_text(hjust = 0))

umap.tumor.gruel.signatures <- t11044x %>%
    FeaturePlot_scCustom(
        pt.size = 2,
        features =names(gruel.liposarcoma.signatures),
        na_cutoff = 0.4,
        raster = TRUE) &
    theme_palbo() +
    theme(axis.title = element_text(hjust = 0),
          axis.text = element_blank())



gruel.df <- t11044x %>%
    FetchData(vars = c("paired.id", "patient.id", "cycle",
                       "biopsy.id", names(gruel.liposarcoma.signatures)))

pd <- position_dodge(width = 0.8) # move them .05 to the left and right

violin.c1c3.tumor.signatures.significance <- gruel.df %>%
    pivot_longer(cols = names(gruel.liposarcoma.signatures),
                 names_to = "signature", values_to = "score") %>%
    ggplot(aes(cycle, score)) +
    geom_violin(aes(fill = cycle),
                position=pd,
                linewidth = 0.1) +
    geom_boxplot(position = pd,
                 linewidth =  0.1,
                 width =  0.08,
                 outlier.shape = NA,
                 outlier.color = "#00000088") +
    ggh4x::facet_grid2(. ~ signature + paired.id, scales = "free") +
    scale_fill_manual(values = cycle.colors) +
    theme_palbo() +
    stat_compare_means(comparisons = list(c("C1", "C3")),
                       method = "wilcox.test",
                       vjust = 1,
                       label = "p.signif",                       
                       bracket.size = .1) 

violin.c1c3.tumor.signatures <- gruel.df %>%
    pivot_longer(cols = names(gruel.liposarcoma.signatures),
                 names_to = "signature", values_to = "score") %>%
    ggplot(aes(paired.id, score,
               fill = cycle)) +
    geom_violin(position=pd,
                linewidth = 0.1) +
    geom_boxplot(position = pd,
                 linewidth =  0.1,
                 width =  0.08,
                 outlier.shape = NA,
                 outlier.color = "#00000088") +
    ggh4x::facet_grid2(~ signature, scales = "free", independent = "y") +
    scale_fill_manual(values = cycle.colors) +
    ylab("Signature Expression") + xlab("Patient ID") +
    theme_palbo() + NoLegend()

1

pdf("figures/Supplementary_Figure_S1C.pdf",
    width =  210/25.4, height = 41/25.4)
violin.c1c3.tumor.signatures
violin.c1c3.tumor.signatures.significance
dev.off()



gruel.signature.colors <- pals::kelly(8)[-(1:2)]
names(gruel.signature.colors) <- names(gruel.liposarcoma.signatures)
##
bar.signature.classification <- t11044x[[]] %>%
    ggplot(aes(y = fct_rev(sample.id), fill = GRUEL_Signature_Classification)) +
    geom_bar(position = "fill") +
    coord_cartesian(expand = FALSE) +
    facet_grid(paired.id ~ ., scales = "free", space = "free") +
    scale_fill_manual(name = "Signature Class",
                      values = gruel.signature.colors) +
    theme_palbo()
           


pdf("figures/Supplementary_Figure_S1A_B_C.pdf",
    width =  210/25.4, height = 123/25.4)
## 
layout <- "
AAB#CD
AAEFGH
IIIIII"
##
( (umap.tumor.gruel.signature.class +
   theme(legend.key.size = unit(3, "mm")) ) +
  (bar.signature.classification + NoLegend() ) +
  (umap.tumor.gruel.signatures[[1]] + NoLegend() ) +
  (umap.tumor.gruel.signatures[[2]] + NoLegend() ) +
  (umap.tumor.gruel.signatures[[3]] + NoLegend() ) +
  (umap.tumor.gruel.signatures[[4]] + NoLegend() ) +
  (umap.tumor.gruel.signatures[[5]] + NoLegend() ) +
  (umap.tumor.gruel.signatures[[6]] + NoLegend() ) +
  violin.c1c3.tumor.signatures ) +
    plot_layout(design =  layout) 
##
( (umap.tumor.biopsy + NoLegend() |
   umap.tumor.clusters + NoLegend() |
    umap.tumor.gruel.signature.class + NoLegend() | 
    bar.signature.classification + NoLegend() |
    plot_spacer() | plot_spacer() ) /
  ( umap.tumor.gruel.signatures[[1]] |
    umap.tumor.gruel.signatures[[2]] |
    umap.tumor.gruel.signatures[[3]] |
    umap.tumor.gruel.signatures[[4]] |
    umap.tumor.gruel.signatures[[5]] |
    umap.tumor.gruel.signatures[[6]] ) /
  violin.c1c3.tumor.signatures ) &
    NoLegend() 
##
umap.tumor.biopsy
umap.tumor.clusters
bar.signature.classification
umap.tumor.gruel.signatures
violin.c1c3.tumor.signatures.significance +
    theme(legend.key.size = unit(3, "mm"))
dev.off()

## __EOF__

