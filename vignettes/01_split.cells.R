## consider removing plots b/c it's taken care of by the shiny app

## ---- NOTE TO SELF:
## use new p11044x -- contains new columns not in p11044s

## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
source("vignettes/main.R")

## ----output_directory_setup---------------------------------------------------
sample.name <- "01.p11044.split_cells"
tabDir <- file.path("tables", sample.name)
figDir <- file.path("figures", sample.name)
tmpDir <- file.path("tmp", sample.name)
tmpImage <- "p11044.now.rda"
## create directories
sapply(c(tabDir, figDir, tmpDir), usethis::use_directory)


## ----load_data----------------------------------------------------------------
## data(p11044s)

## 
p11044s@meta.data$sample.id <- naturalsort::naturalfactor(
                                 paste(p11044s@meta.data$patient.id,
                                       p11044s@meta.data$cycle, sep = "_"))

## set up paired and unpaired data
p11044s@meta.data <- p11044s[[]] %>%
    mutate(
        trt  = cycle,
        sample.id = paste(patient.id, cycle, sep = "_"),
        paired.id = ifelse(patient.id %in% c("P03", "P12", "P13"),
                           patient.id, "Unpaired"),
        paired.id = factor(paired.id, levels = c("P03", "P12", "P13",
                                                 "Unpaired")),
        paired.cycle = paste(paired.id, cycle, sep = "_")
    )


## quantify cell types by predicted lineages based on the
## adiposeref, and pbmcref
## assumption that we would remove any leftover doublets
p11044s@meta.data %>%
    group_by(predicted.celltype.l1, predicted.celltype.l1.pbmcref) %>%
    tally %>% print(n = Inf)

p11044s@meta.data %>%
    group_by(combined.predicted.adipose.pbmc) %>%
    tally %>% print(n = Inf)

## ----selected_dominant_groups-------------------------------------------------
dominant.cell.type.combinations <- c(
    "Adipocyte:B", "ASPC:other", "B:B",
    "Dendritic:B", "Dendritic:DC", "Dendritic:Mono",
    "Endothelial:other", "Macrophage:Mono", "Macrophage:other",
    "Mast:other", "Monocyte:Mono", "NK:CD8 T", "NK:NK",
    "Pericyte:other", "Smooth Muscle:other",
    "T:CD4 T", "T:CD8 T", "T:other T")

## use new p11044x -- contains new columns not in p11044s
## select dominant cell type combinations
p11044x <- p11044s %>%
    subset(combined.predicted.adipose.pbmc %in%
           dominant.cell.type.combinations) %>%
    DietSeurat()

## ---- analyze through seurat pipeline
p11044x <- p11044x %>%
    SCTransform(vars.to.regress = c("percent.mt"),
                ncells = floor(ncol(p11044x)*.15)) %>%
    FindVariableFeatures() %>%
    RunPCA(npcs = 35) %>% 
    FindNeighbors(dims = 1:35) %>%
    FindClusters(resolution = seq(0.4, 1.4, by = .2)) %>%
    RunUMAP(dims = 1:35)



## ----final_cell_type_curation-------------------------------------------------
p11044x$curated.celltype.l1 <- p11044x$predicted.celltype.l1 %>%
    recode("ASPC" =      "Cancer Cells",
           "Adipocyte" = "Plasma Cells",
           "Pericyte" =  "Fibroblasts",
           "Monocyte" =  "Macrophage",
           "Dendritic" = "Macrophage",
           "Smooth Muscle" =  "Fibroblasts") %>%
    factor(levels = c("Cancer Cells", "T", "B", "Macrophage",
                      "NK", "Endothelial", "Fibroblasts",
                      "Plasma Cells", "Mast"))

## 
p11044xL <- p11044x %>%
    SplitObject(split.by = "curated.celltype.l1")


## ----subset_putative_tumor_cells----------------------------------------------
## use curated cell type rather than predicted.celltype.l1 %in% ASPC and Adipocyte
putative.tumor.cells <- p11044x[[]] %>%
    filter(
        curated.celltype.l1 %in% c("Cancer Cells")) %>%
    rownames

## ---- analyze through seurat pipeline
## High MDM2, CDK4, & HMGA2
t11044x <- p11044x[, putative.tumor.cells] %>%
    DietSeurat %>%
    SCTransform(vars.to.regress = c("percent.mt"),
                ncells = 10000) %>%
    FindVariableFeatures() %>%
    RunPCA(npcs = 35) %>% 
    FindNeighbors(dims = 1:35) %>%
    FindClusters(resolution = c(seq(0.6, 1.4, by = .2), 0.4)) %>%
    RunUMAP(dims = 1:35)

## ----subset_putative_immune_cells---------------------------------------------
## immune cell subset
putative.immune.cells <- p11044x[[]] %>%
    filter(curated.celltype.l1 %in% c("Dendritic", "Macrophage",
                                        "B", "Mast", "Monocyte",
                                        "NK", "Neutrophil", "T")) %>%
    rownames

## ---- analyze through seurat pipeline
i11044x <- p11044x[, putative.immune.cells] %>%
    DietSeurat %>%
    SCTransform(vars.to.regress = c("percent.mt"),
                ncells = 10000) %>%
    FindVariableFeatures() %>%
    RunPCA(npcs = 35) %>% 
    FindNeighbors(dims = 1:35) %>%
    FindClusters(resolution = seq(0.4, 1.4, by = .2)) %>%
    RunUMAP(dims = 1:35)

## 
p11044xL <- p11044x %>%
    SplitObject(split.by = "curated.celltype.l1")


## Remov MAST from list ; only 17 cells
p11044xL <- p11044xL[-grep("Mast",names(p11044xL))]

## 
for(i in names(p11044xL)) {
    print(paste(i, "n =", ncol(p11044xL[[i]])))
    p11044xL[[i]] <- p11044xL[[i]] %>%
            SCTransform(vars.to.regress = c("percent.mt"),
                ncells = floor(ncol(p11044x)*.15)) %>%
        FindVariableFeatures() %>%
        RunPCA(npcs = 35) %>% 
        FindNeighbors(dims = 1:35) %>%
        FindClusters(resolution = seq(0.4, 1.4, by = .2)) %>%
        RunUMAP(dims = 1:35)
}


## Remov MAST from list ; only 17 cells
p11044xL <- p11044xL[-grep("Mast",names(p11044xL))]

## 
for(i in names(p11044xL)) {
    print(paste(i, "n =", ncol(p11044xL[[i]])))
    p11044xL[[i]] <- p11044xL[[i]] %>%
            SCTransform(vars.to.regress = c("percent.mt"),
                ncells = floor(ncol(p11044x)*.15)) %>%
        FindVariableFeatures() %>%
        RunPCA(npcs = 35) %>% 
        FindNeighbors(dims = 1:35) %>%
        FindClusters(resolution = seq(0.4, 1.4, by = .2)) %>%
        RunUMAP(dims = 1:35)
}



usethis::use_data(p11044s, p11044x, p11044xL, i11044x, t11044x, overwrite = TRUE)
usethis::use_data(p11044x, p11044xL, overwrite = TRUE)

## Build Comprehensive Shiny App
DefaultAssay(p11044x) <- "SCT"
shiny.dir <- "shiny/new_p11044x/"
usethis::use_directory(file.path(shiny.dir))

## data set 1 :
scConf1 = createConfig(p11044x)
makeShinyFiles(p11044x, scConf1,
               gex.assay = "SCT",
               gene.mapping = TRUE,
               shiny.prefix = "sc1",
               shiny.dir = shiny.dir,
               default.gene1 = "MDM2", default.gene2 = "CD274",
               default.dimred = c("UMAP_1", "UMAP_2"))
## data set 2 :
scConf2 = createConfig(t11044x)
makeShinyFiles(t11044x, scConf2,
               gex.assay = "SCT",
               gene.mapping = TRUE,
               shiny.prefix = "sc2",
               shiny.dir = shiny.dir,
               default.gene1 = "MDM2", default.gene2 = "CD274",
               default.dimred = c("UMAP_1", "UMAP_2"))
## data set 3 :
scConf3 = createConfig(i11044x)
makeShinyFiles(i11044x, scConf3,
               gex.assay = "SCT",
               gene.mapping = TRUE,
               shiny.prefix = "sc3",
               shiny.dir = shiny.dir,
               default.gene1 = "MDM2", default.gene2 = "CD274",
               default.dimred = c("UMAP_1", "UMAP_2"))

#' cell type specific subsets
#' @param so.list seurat object list
#' @param shinly.dir shiny output directory
#' @param skip.n number of shiny's to skip.  The starting 'sc#' on
#' the output will be (skip.n + 1) e.g. skip.n = 3 
#' 
make_shiny_from_so_list <- function(so.list, shiny.dir, skip.n = 0) {
    ## 
    nx <- 1:length(so.list)
    ## 
    lapply(nx, function(i) {
    ## 
        scConfx <- createConfig(so.list[[i]])
        shiny.n <- i + skip.n
        makeShinyFiles(so.list[[i]], scConfx,
                       gex.assay = "SCT",
                       gene.mapping = TRUE,
                       shiny.prefix = paste0("sc", shiny.n),
                       shiny.dir = shiny.dir,
                       default.gene1 = "MDM2", default.gene2 = "PTPRC",
                       default.dimred = c("UMAP_1", "UMAP_2"))
        })
    ## 
    return(NULL)
}

## run fx
make_shiny_from_so_list(p11044xL, shiny.dir = shiny.dir, skip.n = 3)

## build the app to contain multiple data sets
makeShinyCodesMulti(
    shiny.title = "Project 11044 : Palbociclib + I/O Clinical Trial, P01>P13",
    shiny.footnotes = "Singer Laboratory / Confidential",
    shiny.prefix = paste0("sc", 1:14),
    shiny.headers = c("All Cells", "Tumor Cells",
                      "Immune Cells",
                      paste(names(p11044xL), "Cells")),
    shiny.dir = shiny.dir)
 
#% 
#% ## plots
#% ## plot construction
#% celltypes.adiposeref.l1 <- p11044s$predicted.celltype.l1 %>%
#%     unique %>% sort 
#% colors.celltypes.adiposeref.l1 <- pals::brewer.paired(n = length(celltypes.adiposeref.l1))
#% names(colors.celltypes.adiposeref.l1) <- celltypes.adiposeref.l1
#% 
#% colors.celltypes.pbmc.l1 <- c(
#%     "other" = "#202020",
#%     "CD4 T" = "#005824",
#%     "CD8 T" = "#238B45",
#%     "other T" = "#66C2A4",
#%     "NK" = "#74A9CF",
#%     "B" = "#F16913",
#%     "Mono" = "#88419D",
#%     "DC" = "#6E016B")
#% 
#% 
#% ## cluster umap
#% umap.cluster <- DimPlot(p11044x, group.by = "seurat_clusters", label = TRUE) +
#%     theme_light() + NoLegend()
#% ## predicted cell types - Level 1
#% umap.celltypes.l1 <- DimPlot(p11044x, group.by = "predicted.celltype.l1") +
#%     theme_light() +
#%     scale_colour_manual(name = "Cell Type",
#%                         values = colors.celltypes.adiposeref.l1)
#% ## liposarcoma cell score
#% umap.lps.score <- FeaturePlot(p11044x, features = "wd.dd.liposarcoma.score",
#%                               cols = pals::brewer.ylgnbu(n = 5)) +
#%     theme_light()
#% ## cd4 t-cell
#% umap.cd4.t.score <- FeaturePlot(p11044x, feature = "T.CD4",
#%                                  cols = pals::brewer.ylgnbu(n = 5)) +
#%     theme_light()+ NoLegend()
#% ## cd8 t-cell
#% umap.cd8.t.score <- FeaturePlot(p11044x, feature = "T.CD8",
#%                                  cols = pals::brewer.ylgnbu(n = 5)) +
#%     theme_light()+ NoLegend()
#% ## NK cell
#% umap.nk.score <- FeaturePlot(p11044x, feature = "NK.CELL",
#%                                 cols = pals::brewer.ylgnbu(n = 5)) +
#%     theme_light() + NoLegend()
#% ## b-cell
#% umap.b.cell.score <- FeaturePlot(p11044x, feature = "B.CELL",
#%                                 cols = pals::brewer.ylgnbu(n = 5)) +
#%     theme_light() + NoLegend()
#% ## macrophages
#% umap.macrophage.score <- FeaturePlot(p11044x, feature = "MACROPHAGE",
#%                                 cols = pals::brewer.ylgnbu(n = 5)) +
#%     theme_light() + NoLegend()
#% 
#% ## pdf(file.path(tmpDir, paste0(sample.name, ".umap.cell_types_comprehensive.pdf")),
#% ##     width = 420/25.4, height = 210/25.4)
#% png(file.path(tmpDir, paste0(sample.name, ".umap.cell_types_comprehensive.png")),
#%     width = 420, height = 210, units = "mm", res = 300)
#% (( umap.cluster | umap.celltypes.l1 | umap.lps.score | umap.cd4.t.score ) /
#%     ( umap.cd8.t.score | umap.nk.score | umap.b.cell.score | umap.macrophage.score)) + plot_layout(guides = "collect")
#% dev.off()


## correct meta.data files from shiny app

( meta.files <- list.files("shiny/p11044x/", pattern = "meta.rds") %>% naturalsort::naturalsort(.))
meta.L <- lapply(meta.files, function(rds) readRDS(file.path("shiny/p11044x/", rds)))

all(sapply(meta.L[-1], function(x) all(names(meta.L[[1]]) %in% names(x))))
all(sapply(meta.L[-1], function(x) all(names(meta.L[[1]]) == names(x))))

setdiff(names(meta.L[[1]]), names(meta.data.current))
setdiff(names(meta.data.current), names(meta.L[[1]]))
intersect(names(meta.data.current), names(meta.L[[1]]))


columns.to.pass <- c(setdiff(names(meta.data.current), names(meta.L[[1]])))

meta.data.current <- meta.data.current %>%
    rownames_to_column(var = "sampleID")

meta.Lx <- lapply(meta.L, function(md) {
    out <- md %>%
        select(-(Study_ID:days.to.last.contact)) %>%
        left_join(meta.data.current[, c("sampleID", columns.to.pass)], by = "sampleID") %>%
        mutate(active.cell = factor(active.cell, c(TRUE, FALSE)),
               evaluable_for_response = factor(evaluable_for_response, c(TRUE, FALSE)),
               best_recist_response = factor(best_recist_response, c(TRUE, FALSE)),
               pr_confirmed = ifelse(is.na(pr_confirmed), 0, TRUE),
               date_first_objective = ifelse(is.na(date_first_objective), 0, date_first_objective),
               tt_first_objective_response = ifelse(is.na(tt_first_objective_response), 0, tt_first_objective_response),
               tx_ongoing = factor(tx_ongoing, c(TRUE, FALSE)),
               active_fu = factor(active_fu,  c(TRUE, FALSE)),
               last_follow_up_dt = ifelse(is.na(last_follow_up_dt), 0, last_follow_up_dt),
               dor = ifelse(is.na(dor), 0, dor),
               pt_dod = factor(pt_dod, c(TRUE, FALSE)),
               IMPACT = factor(IMPACT, c(TRUE, FALSE))) %>%
        select(sampleID:miQC.keep, active.cell,
               bioID, sample.id, paired.id, paired.cycle, trt,
               SCT_snn_res.0.4:mapping.score.pbmcref, 
               combined.predicted.adipose.pbmc,curated.celltype.l1,
               study_id_1:bor_type1, irAE.f, 
               binary.pfs.90d, binary.pfs.180d, binary.pfs.median, pfs.z,
               IMPACT, DMP.ID,
               SASP_Score:lps.down.score,
               MacMono_score:REACTOME_ZINC_TRANSPORTERS,
               JERBY_ANRON_T_CELL_EXCLUSION,
               PC_1:UMAP_2) 
})




names(meta.Lx) <- meta.files

for(i in meta.files) {
    print(i)
    saveRDS(meta.Lx[[i]], file = file.path("shiny/p11044x", i))
}

