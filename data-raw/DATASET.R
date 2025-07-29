## =============================================================================
## This code prepares three data sets :
##   - a native seurat with default import parameters
##   - a seurat object where each sample is processed by doublet finder
##   - a native monocle3 with default import parameters
## 
## All are prefixed with `pre.` in their object name implying minimal data
##  curation. These contain a simple annotation and are made to match as
##  close as possible.
## 
## pre.cds has an min.expr = 0.1 by default which removes some cells that Seurat
##  does not.  This is resolved in subsequent scripts.
##
## pre.seu.dbf is the most comprehensive and will be the default data set used.
##
## The trial had two arms:
##   - The first arm containing patients P01 -> P13, except P05.  This dataset 
##     is refered by it's IGO Project ID as p11044.
##
##   - The second arm contains patients P14 -> P40.  This dataset is refered
##     as p12606.
##  
## =============================================================================
Purposely placed this text unquoted to prevent an accidental .ess.source

## ----libraries----------------------------------------------------------------
library(tidyverse)
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(DoubletFinder)
sapply(list.files("R", pattern = "R$", full.names = TRUE), source)

## other parameters
## set for 500 gb of ram
library(future)
options(future.globals.maxSize = 180 * 1024^3,
        future.seed=TRUE)
#future::plan("multicore", workers=14)

## ----setup--------------------------------------------------------------------
datTmp <- "data-raw/tmp"
usethis::use_directory(datTmp)

set.seed(2024)


## ----retrieve_sample_names_and_10x_directory_paths----------------------------
( sc.files <- list.files("data-raw/p11044_h5/h5/", pattern = "^P.*",
                         include.dirs = TRUE, full.names = TRUE) )
( sample.names <- gsub(".*/(P.*).h5", "\\1", sc.files) )
( sc.files <- setNames(sc.files, sample.names) )


## ---- collect run metrics
( metrics.files <- list.files("data-raw/p11044_h5/metrics_summary/",
                            pattern = "^P.*csv$",
                            include.dirs = TRUE, full.names = TRUE) )
( metrics.files <- setNames(metrics.files, sample.names) )

run.metrics <- do.call(rbind,
                       lapply(metrics.files,
                              function(i) {
                                  mtx <- read.csv(i)
                                  if(ncol(mtx) == 19) {
                                      mtx$Q30.Bases.in.Sample.Index <- NA
                                  }
                                  mtx
                              }))

for(i in 1:ncol(run.metrics)) {
    run.metrics[, i] <- as.numeric(gsub(",|%", "", run.metrics[, i]))
}

##run.metrics[,c(5:17)] <- run.metrics[,c(5:17)]/100
run.metrics <- run.metrics %>%
    rownames_to_column(var = "sample.name") %>%
    transform(row.names = sample.name) %>%
    mutate(use = ifelse(Estimated.Number.of.Cells > 100, TRUE, FALSE))

run.metrics %>% head

## PUBLICATION DATA 18.JUL.2025
usethis::use_data(run.metrics, overwrite = FALSE)

sc.files.use <- sc.files[run.metrics %>%
                         filter(Estimated.Number.of.Cells > 95) %>%
                         pull(sample.name)]


## ----prepare_seurat_objects---------------------------------------------------
## a.data is a list of 10x gene read count matrices
a.data <- read.m.10x(sc.files.use, names(sc.files.use))

## test.data <- lapply(a.data, function(ax) ax[,1:140])
                    
## ----run_doublet_finder_and_prepare_prelim_seurat_object
pre.seu.dbf <- prepare_seu(a.data, version = "v4", merge = TRUE,
                             dbf = TRUE, nbin = 20)

## ----reduce_metadata
original.pre.seu.dbf.meta.data <- pre.seu.dbf[[]]

write.table(original.pre.seu.dbf.meta.data,
            file = gzfile(file.path("data-raw", "pre.seu.dbf.meta.data.txt.gz")),
            sep = "\t", quote = FALSE, row.names = TRUE)
usethis::use_data(original.pre.seu.dbf.meta.data, overwrite = FALSE) ## 14.Aug.2024

reduced.pre.seu.dbf.meta.data <- original.pre.seu.dbf.meta.data %>%
    select(-contains("DF.classifications"), -contains("pANN_"))

write.table(reduced.pre.seu.dbf.meta.data,
            file = gzfile(file.path("data-raw", "pre.seu.dbf.reduced.meta.data.txt.gz")),
            sep = "\t", quote = FALSE, row.names = TRUE)
usethis::use_data(reduced.pre.seu.dbf.meta.data, overwrite = FALSE) ## 14.Aug.2024

pre.seu.dbf@meta.data <- pre.seu.dbf@meta.data %>%
    select(-contains("DF.classifications"), -contains("pANN_")) %>%
    mutate(cycle = ifelse(cycle == "C0", "C1", cycle))

## QC metrics
patient.list <- unique(pre.seu.dbf$patient.id)


## usethis::use_data(pre.seu.L, overwrite = FALSE) ## 14.Aug.2024
## 14.Aug.2024 data freeze
usethis::use_data(pre.seu.dbf, overwrite = FALSE) 


## ----run_MiQC-----------------------------------------------------------------
pre.seu.dbf <- pre.seu.dbf %>% RunMiQC


## ----MiQC_plots_mitochondria--------------------------------------------------
png(file.path(datTmp, "tmp.scatter.miqc_%03d.png"),
    width = 210, height = 210, units = "mm", res = 200)
PlotMiQC(pre.seu.dbf)
dev.off()


## ----append_clinical_outcomes_data_for_project_p11044-------------------------
## merge w/ clinical data
## updated phenotypes

## use median PFS as 213 ; estimated by Li-Xuan
clinical.pheno.p11044 <- read.delim("data-raw/p11044_clinical_phenotypes_NO_PHI.txt") %>%
    mutate(cutoff.date = "2024-05-01")

## not included pre.rep, pre.seu.rep
usethis::use_data(clinical.pheno.p11044,
                  overwrite = TRUE, compress = "xz")




## merge clinical data of p11044
xA <- pre.seu.dbf@meta.data %>%
    rownames_to_column(var = "cellID") %>%
    merge(clinical.pheno.p11044, by = "patient.id",
          all.x = TRUE, sort = FALSE) %>%
    transform(row.names = cellID,
              cellID = NULL)

nrow(pre.seu.dbf@meta.data)
nrow(xA)




## confirm orders & coerce order
if(all(rownames(pre.seu.dbf@meta.data) %in% rownames(xA) &
    nrow(pre.seu.dbf@meta.data) == nrow(xA))) {
    print(1)
    xA <- xA[rownames(pre.seu.dbf@meta.data), ]
} else {
    warning("not all cells are present, double-check")
}


## don't remove data from xA, it will be used to make a copy in xB
if(all(rownames(pre.seu.dbf@meta.data) == rownames(xA))) {
    clin.vars <- setdiff(names(xA), names(pre.seu.dbf@meta.data))
    pre.seu.dbf@meta.data <- cbind(pre.seu.dbf@meta.data, xA[,clin.vars])
    print(1)
} else {
    warning("not all cells are present, double-check")
}

## save completed object as project name with all clinical data
pre.p11044 <- pre.seu.dbf

## save pre.p11044 -- frozen !! -- do not overwrite
usethis::use_data(pre.p11044, overwrite = TRUE)


## __EOF__
