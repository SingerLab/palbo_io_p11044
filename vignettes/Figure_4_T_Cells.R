## =============================================================================
#' 
#' The code in reanalyzes T-cell immune components, and performs the differential
#'  abundance  presented. Generates panels for Figure 4
#'
#' Code contributors: Jasme Lee, Evan Seffar
#'
#' Main Figures
#' 
#' * Figure 4 -- Panel A
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
library(dbplyr)
library(janitor)
library(harmony)
library(scCustomize)
library(RColorBrewer)
library(org.Hs.eg.db)
library(ggalluvial)


## ----output_directory_setup---------------------------------------------------
sample.name <- "04.p11044.T_Cells"
tabDir <- file.path("tables", sample.name)
tmpDir <- file.path("tmp", sample.name)
sapply(c(tabDir, tmpDir), usethis::use_directory)

## ----load_data----------------------------------------------------------------
data(p11044x)

## ---- T-cell Selection 
## patient subset 
( P12_P13 <- p11044x %>%
      subset(
          curated.celltype.l1 == "T" &
          patient.id %in% c("P12", "P13")) )

## metadata thinning
P12_P13@meta.data <- P12_P13@meta.data %>% dplyr::select(
                                                      bioID,
                                                      cycle,
                                                      predicted.celltype.l3.pbmcref,
                                                      predicted.celltype.l2.pbmcref,
                                                      patient.id
                                                  )

## ---- T-cell processing
TCells_Harmonized <- P12_P13 %>%
    RunHarmony(c("cycle", "patient.id"), plot_convergence = TRUE) %>%
    RunUMAP(reduction = "harmony", dims = 1:35) %>%
    FindNeighbors(reduction = "harmony", dims = 1:35) %>%
    FindClusters(resolution =  seq(0.1, 1.2, by =  0.1)) %>%
    identity()
TCells_Harmonized


## ---- get cell type frequencies
Frequency_Alluvial <- TCells_Harmonized@meta.data %>%
    dplyr::select(patient.id,
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
    dplyr::group_by(patient.id, cycle, predicted.celltype.l2.pbmcref) %>%
    dplyr::summarise(n = n(), .groups = "drop_last") %>%
    mutate(Proportion = round(prop.table(n), 3)) %>%
    dplyr::select(-n) %>%
    filter(Proportion > 0)

t_cell.colors <- c(
    "CD4 Proliferating" = "#ABD3E6",
    "CD4 TCM" = "#284EA1",
    "CD4 TEM" = "#ADDE75",
    "CD8 Naive" = "#3AA53D",
    "CD8 Proliferating" = "#F68494",
    "CD8 TCM" = "#E41D25",
    "CD8 TEM" = "#FAB756",
    "dnT" = "#F46E16",
    "gdT" = "#CDAABD",
    "MAIT" = "#F6339A",
    "NK Proliferating" = "#F3E587",
    "NK" = "#AF4C21",
    "Treg" = "#AF4C21")
    "ILC" = "#B15928")

## ---- alluvial plot
Alluvial_Plot <- ggplot(
    data = Frequency_Alluvial,
    aes(
        x = cycle,
        stratum = predicted.celltype.l2.pbmcref,
        alluvium = predicted.celltype.l2.pbmcref,
        y = Proportion,
        fill = predicted.celltype.l2.pbmcref
    )) +
    geom_alluvium(aes(fill = predicted.celltype.l2.pbmcref),
                  width = 0.1) +                  
    geom_stratum(width = 0.4,
                 colour = "white",
                 linewidth = 0.4) +                 
    ylab("Proportion") +
    xlab("Time") + coord_cartesian(expand = FALSE) +
    facet_wrap( ~ patient.id, nrow = 1) +
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
        values = t_cell.colors
    ) +
    theme_linedraw(base_size = 6) +
    theme(
        strip.background = element_blank(),
        legend.key.size = unit(2, 'mm'),
        axis.text.y = element_text(face = "plain"),
        strip.text = element_text(colour = "black", face = "bold"),
        legend.position = "bottom",
        Legend.direction = "horizontal",
        legend.title = element_blank()
    )
Alluvial_Plot

## pdf(file = file.path(tmpDir, "AlluvialPlot_P12_P13.pdf"),
##     width = 4,
##     height = 2.5)
## Alluvial_Plot
## dev.off()

