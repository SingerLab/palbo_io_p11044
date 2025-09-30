# Libraries
library(Seurat)
library(dbplyr)
library(janitor)
library(ggplot2)
library(dplyr)
library(harmony)
library(scCustomize)
library(RColorBrewer)
library(org.Hs.eg.db)
library(ggalluvial)

setwd("~/Downloads/PALBO_ANTIPD1/")

P12_P13 <- subset(x = p11044xL$T, (bioID %in% c("P12", "P13")))
P12_P13@meta.data <- P12_P13@meta.data %>% dplyr::select(
  bioID,
  cycle,
  predicted.celltype.l3.pbmcref,
  predicted.celltype.l2.pbmcref,
  patient.id
)

# TCells <- p11044x$T

TCells_Harmonized <- P12_P13 %>%
  RunHarmony(c("cycle", "patient.id"), plot_convergence = TRUE) %>%
  RunUMAP(reduction = "harmony", dims = 1:35) %>%
  FindNeighbors(reduction = "harmony", dims = 1:35) %>%
  FindClusters(resolution =  c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2)) %>%
  identity()
TCells_Harmonized

Frequency_Alluvial <- TCells_Harmonized@meta.data %>%
  dplyr::select(bioID,
                cycle,
                predicted.celltype.l3.pbmcref,
                predicted.celltype.l2.pbmcref) %>%
  mutate(
    predicted.celltype.l3.pbmcref = unlist(predicted.celltype.l3.pbmcref),
    predicted.celltype.l2.pbmcref = unlist(predicted.celltype.l2.pbmcref)
  ) %>%
  mutate(
    predicted.celltype.l2.pbmcref = ifelse(
      predicted.celltype.l3.pbmcref == "CD4 TCM_2",
      "CD4 TCM 2",
      predicted.celltype.l2.pbmcref
    )
  ) %>%
  dplyr::group_by(bioID, cycle, predicted.celltype.l2.pbmcref) %>%
  dplyr::summarise(n = n(), .groups = "drop_last") %>%
  mutate(Proportion = round(prop.table(n), 3)) %>%
  dplyr::select(-n)

Alluvial_Plot <- ggplot(
  data = Frequency_Alluvial,
  aes(
    x = cycle,
    stratum = predicted.celltype.l2.pbmcref,
    alluvium = predicted.celltype.l2.pbmcref,
    y = Proportion,
    fill = predicted.celltype.l2.pbmcref
  )
) +
  geom_alluvium(aes(fill = predicted.celltype.l2.pbmcref), width = 0.1) +
  geom_stratum(width = 0.4,
               colour = "white",
               linewidth = 0.4) +
  ylab("Proportion") +
  xlab("Time") + coord_cartesian(expand = FALSE) +
  facet_wrap( ~ bioID, nrow = 1) +
  theme() +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(size = 0.1)
    ),
    shape = guide_legend(override.aes = list(size = 0.5)),
    color = guide_legend(override.aes = list(size = 0.5))
  ) +
  scale_fill_manual(
    values = c(
      "#A6CEE3",
      "#1F78B4",
      "#08519C",
      "#B2DF8A",
      "#33A02C",
      "#FB9A99",
      "#E31A1C",
      "#FDBF6F",
      "#FF7F00",
      "#CAB2D6",
      "#6A3D9A",
      "#FFFF99",
      "#B15928"
    )
  ) +
  theme_linedraw(base_size = 6) +
  theme(
    strip.background = element_blank(),
    legend.key.width = unit(6, 'mm'),
    axis.text.y = element_text(face = "plain"),
    strip.text = element_text(colour = "black", face = "bold"),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.title = element_blank()
  )
Alluvial_Plot

pdf(file = "FIGURES/AlluvialPlot_P12_P13.pdf",
    width = 4,
    height = 2.5)
Alluvial_Plot
dev.off()
