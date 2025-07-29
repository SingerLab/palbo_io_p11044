## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----setup--------------------------------------------------------------------
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(patchwork)

## load custom functions: devtools::load_all() takes a very long time to load
sapply(list.files("R", pattern = "R$", full.names = TRUE), source)

## ----output_directory_setup---------------------------------------------------
sample.name <- "00.p11044.pre_processing"
tabDir <- file.path("tables", sample.name)
tmpDir <- file.path("tmp", sample.name)
## create directories
sapply(c(tabDir, tmpDir), usethis::use_directory)


## ----load_data----------------------------------------------------------------
data(pre.p11044)

## ----data_selection_and_cell_cleanup------------------------------------------
#' pre-filtering based on data quantity per cell
fColors <- c("nFeature_RNA" = "#7B3294", "nCount_RNA" = "#008837")

umi.gene.density <- pre.p11044@meta.data %>%
    dplyr::select(nFeature_RNA, nCount_RNA) %>%
    pivot_longer(cols = everything(),
                 names_to = "metric",
                 values_to = "n") %>%
        ggplot(aes(x = log2(n), colour = metric)) +
        geom_density(lwd = 1) +
        scale_colour_manual(values = fColors) +
        geom_vline(xintercept = c(8, 12.5, 8.5, 15),
                   col = rep(fColors, each = 2), lty = 2) +  ## log2(256)
        theme_palbo()

pdf(file.path(tmpDir, "tmp.density.umi_gene.thresholds_%03d.pdf"))
umi.gene.density
dev.off()

#' pre-filter thresholds
#' manually selected to clip the ends of reads and feature distributions
min.genes <- 2^8
max.genes <- 2^12.5
min.reads <- 2^8.5
max.reads <- 2^15


#' filter 1, focus on cells within range of genes and counts
f1 <- pre.p11044@meta.data %>% 
    filter(nCount_RNA >= min.reads, nCount_RNA <= max.reads) %>%
    filter(nFeature_RNA >= min.genes, nFeature_RNA <= max.genes) %>%
    rownames

#' scatterplot between n.genes vs percent mitochondrial, ribosomal, hemoglobin,
#' and platelet gene expression
pct.metrics <- c("percent.mt", "percent.rb", "percent.hb", "percent.plat")

#' diagnostic scatter plots
#' UMI vs number of genes
umi_gene.scatter <- pre.p11044@meta.data[f1,] %>%
    select(orig.ident, nCount_RNA, nFeature_RNA,
           contains("percent")) %>%
    ggplot(aes(nCount_RNA, nFeature_RNA)) +
    geom_point(aes(colour = percent.rb), alpha = 0.5) +
    geom_smooth(se = FALSE, method = "lm") +
    facet_wrap(. ~ orig.ident, ncol = 5) +
    theme_palbo()

png(file.path(tmpDir, "tmp.scatter.umi_gene_%03d.png"),
    width = 298, height = 596, units = 'mm', res = 200)
umi_gene.scatter +
    theme(panel.spacing = unit(0, "mm"))
dev.off()

#' percent metrics vs number of genes
pct.of.total.scatter <- pre.p11044@meta.data[f1,] %>%
    select(nFeature_RNA, contains("percent")) %>%
    pivot_longer(cols = contains("percent"),
                 names_to = "metric", values_to = "pct") %>% 
    ggplot(aes(nFeature_RNA, pct)) +
    geom_point(alpha = 0.2) +
    geom_smooth() +
    facet_wrap(~ metric, ncol = 2, scales = "free_x") +
    theme_palbo()

png(file.path(tmpDir, "tmp.scatter.pct_metrics_gene_%03d.png"),
    width = 210, height = 210, units = "mm", res = 200)
pct.of.total.scatter + coord_cartesian(expand = FALSE)
dev.off()

## select set 1 patient
hb.max <- 3
rb.max <- 50
plat.max <- 0.7
mt.max <- 20


## doublet rates
png(file.path(tmpDir, "tmp.bar.doublet_fq%03d.png"))
pre.p11044[[]] %>%
    ggplot(aes(x = orig.ident)) +
    geom_bar(aes(fill = df.class.high)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))
dev.off()

## MiQC
tx <- pre.p11044
tx@misc[["flexmix_model"]] <- NULL

tx <- RunMiQC(tx)

png(file.path(tmpDir, "tmp.scatter.miQC_%03d.png"))
PlotMiQC(tx) +
    ggplot2::scale_color_gradient(low = "grey", high = "purple")
dev.off()

rm(tx)

## Remove resected tumors 
##  most were too far out from the trial end date 
##  if needed, analyze separetly
s1patients <- c("P01", "P02", "P03", "P06", "P07",
                "P08", "P09", "P11", "P12", "P13")
filtered.cells.biopsy.use <- pre.p11044[[]] %>%
    filter(miQC.keep == "keep",
           df.class.high == "Singlet",
           percent.hb < hb.max,
           percent.rb < rb.max,
           percent.plat < plat.max,
           patient.id %in% s1patients,
           cycle != "S0")


## ----estimate_number_of_pc_based_on_cummulative_variance_explained------------
## Biopsy Only
## scale and pca 
pre.x <- pre.p11044[,rownames(filtered.cells.biopsy.use)] %>%
    SCTransform(vars.to.regress = "percent.mt",
                ncells = floor(nrow(filtered.cells.biopsy.use) * .15)) %>%
    RunPCA(npcs = 100)

prop_varexp <- pre.x@reductions$pca@stdev^2 / sum(pre.x@reductions$pca@stdev^2)
thrs <- c(80, 85, 90, 95)/100
pcs_at_thrs <- c(which(cumsum(prop_varexp) >= 0.80)[1],
                 which(cumsum(prop_varexp) >= 0.85)[1],
                 which(cumsum(prop_varexp) >= 0.90)[1],
                 which(cumsum(prop_varexp) >= 0.95)[1])
names(pcs_at_thrs) <- paste0(c(80, 85, 90, 95), "%")

pdf(file.path(tmpDir, paste0(sample.name, "_pc_variance_explained_%03d.pdf")))
ElbowPlot(pre.x, ndims = 100)
plot(cumsum(prop_varexp))
abline(h = thrs, col = "gray", lty = 3)
abline(v = pcs_at_thrs, col = "gray20", lty = 2)
text(x = pcs_at_thrs, y = rep(0.2, 4), labels = pcs_at_thrs,
     pos = 4)
dev.off()

## only a minor clean / PD affects mitochondrial expression
## not sure what to make of ribosomal
## finish pipeline
#' pcs_at_thrs =  
#' use 35
npcs <- 35 ## pcs_at_thrs["85%"]

set.seed(2021)
#' 35 PCs explains 85 % of the variance
p11044s <- pre.p11044[, rownames(filtered.cells.biopsy.use)] %>%
    SCTransform(vars.to.regress = c("percent.mt"),
                ncells = 20000) %>%
    FindVariableFeatures() %>%
    RunPCA(npcs = pcs_at_thrs["85%"]) %>% 
    FindNeighbors(dims = 1:pcs_at_thrs["85%"]) %>%
    FindClusters(resolution = seq(0.4, 1.4, by = .2)) %>%
    RunUMAP(dims = 1:pcs_at_thrs["85%"])
usethis::use_data(p11044s)


## ----plot_typed_cells---------------------------------------------------------
#' Cell typing with Azimut as suggested by Evan Seffar
#' code goes here.

## run azimuth
library(Azimuth)
library(SeuratData)
## available_data <- AvailableData()
## InstallData("adiposeref")
## InstallData("pbmcref")
#### InstallData("pbmcsca")
## Cell type annotation with PBMC reference
## provides higher resoultion of immune cells
p11044si <- RunAzimuth(p11044s, reference = "pbmcref")
## Cell type annotation with Adipose tissue reference
## provides generalized cell type for other cell types
p11044s <- RunAzimuth(p11044s, reference = "adiposeref")

pbmcref.pred <- p11044si[[]] %>% dplyr::select(nCount_refAssay:mapping.score)
names(pbmcref.pred)  <- paste(names(pbmcref.pred), "pbmcref", sep = ".")

if( all(rownames(p11044s[[]]) == rownames(pbmcref.pred)) ) {
    print(1)
    p11044s@meta.data <- cbind(p11044s@meta.data, pbmcref.pred)
    rm(p11044si)
} else {
    stop("verify row order in ref.pred")
}

p11044s@meta.data <- p11044s@meta.data %>%
    mutate(combined.predicted.adipose.pbmc = paste(
               predicted.celltype.l1,
               predicted.celltype.l1.pbmcref, sep = ":"))

## ----build_early_shiny_for_inspection-----------------------------------------
library(ShinyCell)
DefaultAssay(p11044s) <- "SCT"
shiny.dir <- "shiny/pre.p11044/"
usethis::use_directory(file.path(shiny.dir))
## setwd(shiny.dir)
scConf <- createConfig(p11044s)
makeShinyApp(p11044s, scConf = scConf, shiny.dir = shiny.dir,
             gex.assay = "SCT",
             gene.mapping = TRUE,
             shiny.title = "Palbociclib + I/O Clinical Trial (all cells)")

## ----calculate_cell_and_pathway_scores----------------------------------------
## data(cell.markers) ## data frame of markers
#' list of markers
#' data(oncoKB)
## ls.deg <- read.delim("data-raw/DEG_ANNOTATED_WDLS_DDLS_NORMAL.txt")
## ls.deg %>% filter() %>% pull(Gene.symbol)
## intersect(ls.deg %>% pull(Gene.symbol),
##           oncoKB %>% filter(Is.Tumor.Suppressor.Gene) %>%
##           pull("Hugo.Symbol"))

## data(cell.type.gene.modules, 
##      azizi.cells, jerby.anron.cell.types,
##      szabo.tcells, sasp)
data(sasp, gruel.liposarcoma.signatures)
##
##sarcoma.signatures <- list(wd.dd.liposarcoma.score = c("MDM2", "CDK4", "HMGA2"),
##                           lps.up.score = ls.deg %>%
##                               filter(STATUS == "UP REGULATED IN WDLS_DDLS") %>%
##                               pull(hgnc.symbol),
##                           wdls.up.score = ls.deg %>%
##                               filter(STATUS == "UP REGULATED IN WDLS") %>%
##                               pull(hgnc.symbol),
##                           ddls.up.score = ls.deg %>%
##                               filter(STATUS == "UP REGULATED IN DDLS") %>%
##                               pull(hgnc.symbol),
##                           lps.down.score = ls.deg %>%
##                               filter(STATUS == "UP REGULATED IN NORMAL") %>%
##                               pull(hgnc.symbol))

reactome <- msigdbr::msigdbr(category = "C2", subcategory = "REACTOME") %>%
    dplyr::select(gs_name, gene_symbol) %>%
    split(f = .$gs_name) %>%
    lapply(., function(i) as.character(i$gene_symbol))

reactome.pathways.use <- reactome[
    c("REACTOME_CELL_CYCLE",
      "REACTOME_CELL_CYCLE_CHECKPOINTS",
      "REACTOME_DNA_REPAIR",
      "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION")
]

## ----module scores
p11044s <- AddModuleScore(p11044s,
                          features = c(sasp,
                                       reactome.pathways.use,
                                       gruel.liposarcoma.signatures),
                          seed = 20210724,
                          search = TRUE)

names(p11044s@meta.data)[grep("Cluster", names(p11044s@meta.data))] <- names(
    c(sasp, reactome.pathways.use, gruel.liposarcoma.signatures))



## ----data standardized
p11044s$irAE <- factor(p11044s$irAE, c("TRUE", "FALSE"))
p11044s$irAE.f <- factor(p11044s$irAE, c("TRUE", "FALSE"))

## remove for later
## 
p11044s$binary.pfs.90d <- ifelse(p11044s$PFS >= 90, ">3 mo", "<3 mo")
p11044s$binary.pfs.180d <- ifelse(p11044s$PFS >= 180, ">6 mo", "<6 mo")


##  arrange sample.id and other metadata
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



## ----build_early_shiny_for_inspection-----------------------------------------
DefaultAssay(p11044s) <- "SCT"
shiny.dir <- "shiny/pre.p11044/"
usethis::use_directory(file.path(shiny.dir))
## setwd(shiny.dir)
scConf <- createConfig(p11044s)
makeShinyApp(p11044s, scConf = scConf, shiny.dir = shiny.dir,
             gex.assay = "SCT",
             gene.mapping = TRUE,
             shiny.title = "Palbociclib + I/O Clinical Trial (all cells)")

usethis::use_data(p11044s, p11044x, overwrite = TRUE)

message("Preprocessing Complete")

## ----session_info-------------------------------------------------------------
sessionInfo()
##> 
##> R version 4.0.5 (2021-03-31)
##> Platform: x86_64-conda-linux-gnu (64-bit)
##> Running under: CentOS Linux 7 (Core)
##> 
##> Matrix products: default
##> BLAS/LAPACK: /juno/work/singer/opt/miniconda3/envs/single-cell-rnaseq/lib/libopenblasp-r0.3.15.so
##> 
##> locale:
##>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##>  [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
##> 
##> attached base packages:
##>  [1] grid      stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
##> 
##> other attached packages:
##>  [1] ShinyCell_2.1.0             RColorBrewer_1.1-2          glue_1.4.2                 
##>  [4] gridExtra_2.3               R.utils_2.10.1              R.oo_1.24.0                
##>  [7] R.methodsS3_1.8.1           reticulate_1.20             hdf5r_1.3.3                
##> [10] Matrix_1.5-3                data.table_1.14.0           pbmcsca.SeuratData_3.0.0   
##> [13] pbmcref.SeuratData_1.0.0    adiposeref.SeuratData_1.0.0 SeuratData_0.2.2           
##> [16] Azimuth_0.4.6               shinyBS_0.61.1              patchwork_1.1.1            
##> [19] ROCR_1.0-11                 KernSmooth_2.23-20          sctransform_0.3.5          
##> [22] fields_14.1                 viridis_0.6.1               viridisLite_0.4.0          
##> [25] spam_2.9-1                  future_1.21.0               DoubletFinder_2.0.3        
##> [28] SeuratWrappers_0.3.1        SeuratObject_4.1.3          Seurat_4.3.0               
##> [31] monocle3_1.0.0              SingleCellExperiment_1.12.0 SummarizedExperiment_1.20.0
##> [34] GenomicRanges_1.42.0        GenomeInfoDb_1.26.7         IRanges_2.24.1             
##> [37] S4Vectors_0.28.1            MatrixGenerics_1.2.1        matrixStats_0.58.0         
##> [40] Biobase_2.50.0              BiocGenerics_0.36.0         forcats_0.5.1              
##> [43] stringr_1.4.0               dplyr_1.0.6                 purrr_0.3.4                
##> [46] readr_1.4.0                 tidyr_1.1.3                 tibble_3.1.2               
##> [49] ggplot2_3.5.1               tidyverse_1.3.1            
##> 
##> loaded via a namespace (and not attached):
##>   [1] shinydashboard_0.7.2   utf8_1.2.1             spatstat.explore_3.0-5 tidyselect_1.1.1      
##>   [5] htmlwidgets_1.5.3      Rtsne_0.15             munsell_0.5.0          codetools_0.2-18      
##>   [9] ica_1.0-2              DT_0.18                miniUI_0.1.1.1         withr_3.0.1           
##>  [13] spatstat.random_3.0-1  colorspace_2.0-1       progressr_0.13.0       knitr_1.33            
##>  [17] rstudioapi_0.13        tensor_1.5             listenv_0.8.0          labeling_0.4.2        
##>  [21] GenomeInfoDbData_1.2.4 polyclip_1.10-0        bit64_4.0.5            farver_2.1.0          
##>  [25] rprojroot_2.0.2        parallelly_1.25.0      vctrs_0.6.5            generics_0.1.0        
##>  [29] xfun_0.23              R6_2.5.0               rsvd_1.0.5             msigdbr_7.5.1         
##>  [33] flexmix_2.3-19         bitops_1.0-7           spatstat.utils_3.0-1   DelayedArray_0.16.3   
##>  [37] assertthat_0.2.1       promises_1.2.0.1       scales_1.3.0           googlesheets4_0.3.0   
##>  [41] nnet_7.3-16            gtable_0.3.0           globals_0.14.0         goftest_1.2-2         
##>  [45] rlang_1.1.4            splines_4.0.5          lazyeval_0.2.2         gargle_1.1.0          
##>  [49] spatstat.geom_3.0-3    broom_0.7.6            BiocManager_1.30.15    reshape2_1.4.4        
##>  [53] abind_1.4-5            modelr_0.1.8           backports_1.2.1        httpuv_1.6.1          
##>  [57] SeuratDisk_0.0.0.9020  tools_4.0.5            usethis_2.0.1          ellipsis_0.3.2        
##>  [61] ggridges_0.5.3         Rcpp_1.0.10            plyr_1.8.6             zlibbioc_1.36.0       
##>  [65] RCurl_1.98-1.3         deldir_1.0-6           pbapply_1.4-3          cowplot_1.1.1         
##>  [69] zoo_1.8-9              haven_2.4.1            ggrepel_0.9.1          cluster_2.1.2         
##>  [73] fs_1.5.0               magrittr_2.0.1         scattermore_0.7        lmtest_0.9-38         
##>  [77] reprex_2.0.0           RANN_2.6.1             googledrive_1.0.1      fitdistrplus_1.1-3    
##>  [81] shinyjs_2.0.0          hms_1.1.0              mime_0.10              xtable_1.8-4          
##>  [85] readxl_1.3.1           compiler_4.0.5         maps_3.3.0             crayon_1.4.1          
##>  [89] htmltools_0.5.1.1      mgcv_1.8-35            later_1.2.0            lubridate_1.7.10      
##>  [93] DBI_1.1.1              dbplyr_2.1.1           rappdirs_0.3.3         MASS_7.3-54           
##>  [97] babelgene_21.4         cli_3.6.0              naturalsort_0.1.3      dotCall64_1.0-2       
##> [101] igraph_1.2.6           pkgconfig_2.0.3        sp_1.6-0               plotly_4.9.3          
##> [105] spatstat.sparse_3.0-0  xml2_1.3.2             XVector_0.30.0         rvest_1.0.0           
##> [109] digest_0.6.27          RcppAnnoy_0.0.18       spatstat.data_3.0-0    cellranger_1.1.0      
##> [113] leiden_0.3.7           uwot_0.1.14            curl_4.3.1             shiny_1.6.0           
##> [117] modeltools_0.2-23      lifecycle_1.0.3        nlme_3.1-152           jsonlite_1.7.2        
##> [121] desc_1.3.0             fansi_0.4.2            pillar_1.6.1           lattice_0.20-44       
##> [125] fastmap_1.1.0          httr_1.4.2             survival_3.2-11        remotes_2.3.0         
##> [129] png_0.1-7              bit_4.0.4              presto_1.0.0           stringi_1.6.2         
##> [133] irlba_2.3.3            future.apply_1.7.0    

## ___EOF___

