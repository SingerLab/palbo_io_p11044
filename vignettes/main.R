## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----setup--------------------------------------------------------------------
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(tidyverse)
library(patchwork)
library(ggpubr)
## library(scCustomize)
## library(harmony)

## load custom functions: devtools::load_all() takes a very long time to load
sapply(list.files("R", pattern = "R$", full.names = TRUE), source)


## ----color_maps----------------------------------------------------------------
patient.colors <- c(
    "P01" = "#AECDE1",
    "P02" = "#3C76AF",
    "P03" = "#549E3F",
    "P06" = "#ED9F9B",
    "P07" = "#D0352B",
    "P08" = "#F3C17B",
    "P11" = "#EF8632",
    "P12" = "#643F95",
    "P13" = "#A65E34")

biopsy.colors <- c(
    "P01_C3_D1" = "#AECDE1",
    "P02_C3_D1" = "#3C76AF",
    "P03_C1_D14" = "#BBDE93",
    "P03_C3_D1" = "#549E3F",
    "P06_C3_D1" = "#ED9F9B",
    "P07_C1_D14" = "#D0352B",
    "P08_C1_D14" = "#F3C17B",
    "P11_C3_D1" = "#EF8632",
    "P12_C1_D14" = "#C6B3D3",
    "P12_C3_D1" = "#643F95",
    "P13_C1_D14" = "#FFFFA6",
    "P13_C3_D1"  = "#A65E34")

paired.id.colors = c("P03" = "#000000",
                     "P12" = "#000000",
                     "P13" = "#000000",
                     "Unpaired" = "#888888")

phase.colors <- c("G1" = "#4DAC26", "S" = "#F1B6DA", "G2M" = "#D01C8B")

cell.type.colors <- c(
    "Cancer Cells" = "#333333",
    "Plasma Cells" = "#FFABD0",
    "Mast" = "#8DD3C7",
    "Macrophage" = "#D9A0AB",
    "B" = "#A9A0B2",
    "T" = "#EFB36E",
    "Endothelial" = "#BBDC76",
    "Myeloid" = "#F0D1E1",
    "NK" = "#CCCCFF",
    "Smooth Muscle" = "#FFED6F",
    "Fibroblasts" = "#FFED6F"
)

## change -- one dark one light / helps w/ colorblind safe
## For Cohort A:  C0 / C1 / C3 refer to baseine, Pre-IO, and Post+IO, respectively
cycle.colors <- c( 
    "C0" = "#BABABA",
    "C1" = "#8CD2C7",
    "C3" = "#C51B7D")


## Best overall rseponse colors 
## RECIST response colors 
bor.colors <- c(
    "PD" = "#A3D276",
    "SD" = "#4C86C6",
    "PR" = "#00BFC4",
    "NE" = "#E6E6E6"
)

pfs90.colors <- c(
    "<3 mo" = "#AF79B5",
    ">3 mo" = "#66C2A4" ## "#81A975"
)

irae.colors <- c(
    "TRUE" = "#0570B0", ## "#B0DEED"
    "FALSE" = "#FFD1A9"
)


tissue.colors <- c(
    "Biopsy" = "#C51B7D",
    "PBMC" = "#A50026",
    "Resection" = "#9970AB"
)
    
sample.type.colors <- c(
    "Tumor" = "#4D4B4C",
    "PBMC" = "#A50026",
    "Normal" = "#FEE08B"
)
    
bor.shapes <- c("PD" = 16, "SD" = 17, "PR" = 15)

reason.disc.shapes <- c(
    "Progression of Disease (clinical or radiographic)" = 3,
    "Surgery" = 7,
    "Toxicity" = 8,
    "Remains on treatment" = 17)
    


#' Palbo + IO Clinical Trial ggplot theme
#'
#' The project theme is a modification of \code{\link[ggplot2]{theme_linedraw}}.
#' 
#' We modified the strip background to white, and strip text to black. In
#'  addition panel grid lines were removed, and font size was set to 6
#'  by default. 
#' 
#' @param base_size base fontsize for the plot, default 6
#'
#' @param ... additional parameters passed to \code{\link[ggplot2]{theme_lindraw}
#'            and \code{\link[ggplot2]{theme}}
#'
#' @examples
#'
#' df <- data.frame(patient.id = LETTERS[1:9],
#'                  bor = factor(rep(c("PD", "SD", "PR"), each = 3),
#'                              c("PD", "SD", "PR")),
#'                  pct.change = c(seq(60, 50, by = -5),
#'                                seq(10, -10, by = -10),
#'                                seq(-40, -50, by = -5)))
#'
#' p <- ggplot(df, aes(patient.id, pct.change, fill = bor)) +
#'      geom_bar(stat = "identity") +
#'      scale_fill_manual(values = bor.colors) +
#'      coord_cartesian(expand = FALSE) +
#'      theme_palbo()
#'
#' @export
theme_palbo <- function(base_size = 6, ...) {
    ## use linedraw as base
    theme_linedraw(base_size = base_size, ...) %+replace%
    theme(panel.grid = element_blank(),
          legend.position = "bottom",
          strip.background = element_rect(colour = "#FFFFFF00",
                                          fill = "#FFFFFF00"),
          strip.text = element_text(colour = "#000000"),
          ...)
}


## ====Cell_Type_Markers========================================================
curated.cell.type.markers <- c(
    ## Cancer Cell ? (n=215,813)
    "MDM2", "CDK4", "HMGA2",
    ## T (n=180,214) 
    ## NK (n=70,603)
    "CD3D", "CD4", "CD8A", "NKG7",    
    ## Macrophage (n=54,722)
    ## Monocyte (n=37,838) 
    "CD14", "CD163", "CD68",
    ## Dendritic (n=17,361)
    "ITGAX",
    ## B (n=10,301)
    "CD19", "MS4A1",
    ## Mast (n=900)
    "CPA3", "KIT",
    ## ASPC (n=1,295)
    "APOD",
    ## Adipocyte (n=353)
    "ADIPOQ", "FABP4", "PLIN1",
    ## Endothelial (n=19,405)
    "VWF",
    ## Mesothelial (n=832)
    "SUSD2",
    ## Fibroblasts (n=9,829)
    ## Pericyte (n=1,870) 
    "ACTA2", "PDGFRB",
    ## Lymphatic Endothelial (n=831) 
    "PROX1", "FLT4",
    ## Other (n=2,192)
    ## Unknown (n=93,313)
    "TOP2A", "MKI67"
)

extended.cell.type.markers <- c(
    ## cancer cells
    "MDM2", "CDK4", "HMGA2",
    ## T/NK
    "CD3D", "CD4", "CD8A", "CD96", "NKG7", "IL7R",
    "NCAM1",  "FCGR3A", "CD160",  "ITGA2",  "KLRB1",
    ## Macrohages
    "CD163", "CD14", "CD33", "CD68", "MRC1", "CSF1R",
    ## B cell
    "CD19", "MS4A1", "CD79A", "CD74", "HLA-DRA",
    ## MAST
    "KIT", "CTSG", "TPSAB1", "HPGD",
    ## Plasma
    "CD38", "HLA-DQA1",
    ## ASPC
    "APOD", "CD34", "FBLN5", "PDGFRL",
    ## adipocyte
    "ADIPOQ", "PLIN1", "PNPLA2", "LIPE", "FABP4",
    ## Smooth Muscle / Fibroblast
    "ACTA2", "TAGLN", "MYH11", "RGS5", "GJA4",
    ## Neutral CAF
    "S100A4", 
    ## CAF Acomplice
    "FAP", "VIM", "PDGFRB", "CD70", "PDPN", "ITGA5", "MME",
    ## CAF Defender
    "CIITA", "MCAM", "CAV1", "PDGFRA",
    ## Endothelial
    "CDH5", "VWF", "PECAM1", "ESAM", "FLT1",
    ## MuSC
    "DLK1", "MYF5", "PAX7", "FGFR4",
    ## Skeletal Muscle
    "ACTN2", "ACTA1", "DES", "MB",
    ## RBC
    "HBB",
    ## CYCLING
    "TOP2A", "MKI67", "SUSD2",
    ## mesenchymal cell markers
    "THY1", "NT5E", "ENG",
    "PROM1", "LGR5", "ENTPD2",
    ## Others    
    "IGHG2", "CPA3"
    )

