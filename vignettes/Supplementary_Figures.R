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


