#' read multiple 10x outputs
#' @param sample.dir sample directory containing  `outs/filtered_feature_bc_matrix`
#' @param sample.names names to use for the samples e.g. directory name w/o path
read.m.10x <- function(files, sample.names) {

    if(is.null(sample.names)) {
        sample.names <- names(files)
    }
    
    ## read all files
    message(date(), "  start")
    a.data <- lapply(files,
                  function(i) {
                      out <- Read10X_h5(i)
                      message(date(), "  finished reading: ", i)
                      return(out)
                  })
    
    names(a.data) <- sample.names

    return(a.data)
}


#' parse orig.ident
#'
#' 
parse.original.identity <- function(pso) {

    ## parse sample info from orig.ident
    cx <- as.character(pso@meta.data$orig.ident)
    sAnnot <- data.frame(do.call(rbind, strsplit(cx, split = "_")))
    
    ## names depending if it has 3, 4, or 5 columns
    if(ncol(sAnnot) == 3) {
        names(sAnnot) <- c("patient.id", "cycle", "day")
    } else {
        if(ncol(sAnnot) == 4) {
            names(sAnnot) <- c("patient.id", "cycle", "day", "replicate")
        } else {
            if(ncol(sAnnot) == 5) {
                names(sAnnot) <- c("patient.id", "cycle", "day",
                                   "region", "replicate")
            } else {
                if(ncol(sAnnot) == 6) {
                    names(sAnnot) <- c("patient.id", "cycle", "day",
                                       "region", "tissue", "replicate")
                }
            }
        }
    }
    columns.remove <- c("percent.mt", "nCount_SCT", "nFeature_SCT",
                        "SCT_snn_res.0.8", "seurat_clusters")
    columns.keep <- setdiff(names(pso@meta.data), columns.remove)
    
    pso@meta.data <- pso@meta.data[, columns.keep]
    
    pso@meta.data <- cbind(pso@meta.data, sAnnot)
    
    pso@meta.data <- pso@meta.data %>%
        mutate(
            bioID = patient.id,
            region = ifelse("region" %in% names(.), region, NA),
            tissue = ifelse("tissue" %in% names(.), tissue, 
                     ifelse(cycle == "SO", "resection", 
                     ifelse(region == "PBMC", "PBMC", "biopsy"))),
            normalized.days = as.numeric(gsub("D", "", day)) -7,
            ## new columns
            biopsy.id = ifelse(is.na(region),
                               paste(patient.id,
                                     cycle, day, sep = "_"),
                               paste(paste(patient.id, cycle,
                                           day, region, sep = "_")))
        )

    return(pso)

}

#' build s.data ; list of seurat objects
#' @param a.data read counts from read.m.10x list form
#' @param sample.names name of the samples if different than a.data, default is NULL
#' @param dbf run doublet finder
#' @param version version of seurat to use v5 or v4, default is v4
#' @importFrom DoubletFinder dubFinder
#' @export
build.seurat.list <- function(a.data, sample.names = NULL, dbf = TRUE,
                              version = "v4", 
                              max.mt.pct = NULL, min.genes = 200) {
                              
    if(is.null(sample.names)) {
        sample.names <- names(a.data)
    }
    
    s.data <- lapply(sample.names, function(sm) {
        print(sm)
        ## run doublet finder or jsut CreateSeuratObject
        if(dbf) {
            ## create dbf temporary output directory
            dbfPath = "tmp/00.dbf/"
            if(!dir.exists(dbfPath)) { dir.create(dbfPath, recursive = TRUE) }
            ## run dbf
            pso <- dubFinder(a.data[[sm]], sm, png.path = dbfPath,
                             version = version,
                             max.mt.pct = max.mt.pct, min.genes = min.genes)
        } else {
            ## if no dbf, checkv ersion and run either V5 or other
            if(version == "v5") {
                v5.a.data <- CreateAssay5Object(a.data[[sm]])
                pso <- CreateSeuratObject(v5.a.data, project = sm)
            } else {
                pso <- CreateSeuratObject(a.data[[sm]], project = sm)
            }
        }
        
        ## parse sample info from orig.ident
        pso <- parse.original.identity(pso = pso)
            
        
        DefaultAssay(pso) <- "RNA"
        pso <- DietSeurat(pso, assays = "RNA",
                          misc = FALSE)
        
        
        
        return(pso)
    })

    names(s.data) <- sample.names
    
    return(s.data)
    
}

#' prepare a joined seurat object from a list of files
#'
#' @export
prepare_seu <- function(a.data, sample.names = NULL,
                        version = "v4",
                        merge = TRUE,
                        project = "Palbo_aPD1",
                        dbf = FALSE,
                        nbin = 24,
                        seed = 2021, ...) {
    
    if(is.null(sample.names)) {
        sample.names <- names(a.data)
    }
    
    message(Sys.time(), "  start")
    s.data <- build.seurat.list(a.data, sample.names = sample.names,
                                version = version, dbf = dbf, ...)
    
    message(Sys.time(), "  merging seurat")
    
    if(merge) {
        ## merge data w/o doublet tagging
        pre.seu <- merge(s.data[[1]], s.data[2:length(s.data)],
                         add.cell.ids = sample.names,
                         project = project)
        message(Sys.time(), "  estimate_qc_metrics")
        if(version == "v5") {
            ## specific to version v5
            pre.seu[["RNA"]] <- JoinLayers(pre.seu.[["RNA"]])
            pre.seu <- pre.seu %>% NormalizeData
        }
        
        pre.seu <- pre.seu %>%
            estimate_qc_metrics(nbin = nbin, seed = seed)

        message(Sys.time(), "  categorized_qc")
        pre.seu <- pre.seu %>%
            categorize_qc()
                message(Sys.time(), "  new metadata completed")
        
    } else {

        pre.seu <- merge_by_patient(s.data)
        
        message(Sys.time(), "  estimate_qc_metrics")
        if(version == "v5") {
            for(i in 1:length(pre.seu)) {
                pre.seu[[i]][["RNA"]] <- JoinLayers(pre.seu[[i]][["RNA"]])
                pre.seu[[i]] <- NormalizeData(pre.seu[[i]])
            }
        }
            
        pre.seu <- lapply(pre.seu, estimate_qc_metrics, nbin = nbin, seed = seed)

        message(Sys.time(), "  categorized_qc")
        pre.seu <- lapply(pre.seu, categorize_qc)

        message(Sys.time(), "  new metadata completed")
        
    }
    
    return(pre.seu)
} ## end prepare_seu


#' ----build_seurat_object_with_doublet_detection-------------------------------
prepare_seu_dbf <- function(a.data,  ...) {
    pre.seu.out <- prepare_seu(a.data = a.data, dbf = TRUE,
                               ...)
    return(pre.seu.out)
}
                            

#' ----estimate QC metrics of seurat object-------------------------------------
#' Estimates percent of expressed mitochondrial, ribosomal, hemoglobin, and
#'  platelet genes relative to all genes.  In addition, it runs CellCycleScoring
#'  to estimate S and G2m scores, followed by G1, S, G2m phase prediction. By
#'  default use 24 bins for cell cycle scores and cell cycle phase prediction
#'
#' @param seu a seurat object
#' @param nbin number of bins for cell cycle scoring, default 24
#' @param seed seed for cell cycle scoring, default 2024
#'
#' @export
estimate_qc_metrics <- function(seu, nbin = 24, seed = 2024, ...) {
    seu <- seu %>%
        PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
        PercentageFeatureSet(pattern = "^RP[SL]", col.name = "percent.rb") %>%
        PercentageFeatureSet(pattern = "^HB[^(P)]", col.name = "percent.hb") %>%
        PercentageFeatureSet(pattern = "PECAM1|PF4",
                             col.name = "percent.plat")
        ## CellCycleScoring(s.features = cc.genes.updated.2019$s.genes,
        ##                  g2m.features = cc.genes.updated.2019$g2m.genes,
        ##                  nbin = nbin,
        ##                  seed = seed)
    ## seu$CC.Difference <- seu$S.Score - seu$G2M.Score
    ##

    seu@meta.data <- seu@meta.data %>%
        mutate(across(c("percent.mt", "percent.rb", "percent.hb", "percent.plat"), ~ round(., 4)))
        
    return(seu)
}

#' define active cells for standard, and spliced/unspliced counts
#'
#'  Because spliced/unspliced/ambigous counts are likely coming from the same set of genes,
#'   we use the  max UMI, and max Gene counts across the three categories as the `set of genes`, and
#'   assume this approximates the total number of unique umi and unique genes
define.active.cell <- function(md, min.genes = 200, max.genes = 6000, spliced.counts = FALSE) {

    if(spliced.counts) {
        nCount_RNA <- md %>%
            select(contains("nCount_")) %>%
            rowSums()
        nFeature_RNA <- md %>%
            select(contains("nFeature_")) %>%
            rowSums()

        md$nCount_RNA <- nCount_RNA
        nFeature_RNA <- nFeature_RNA
    }
    
    md = md %>%
        mutate(active.cell = ifelse(nFeature_RNA > min.genes & nFeature_RNA < max.genes,
                                    TRUE, FALSE))
    return(md$active.cell)
    
}



#' ----categorize QC metrics----------------------------------------------------
#' Sets categories for percent.mt, percent.rb, percent.hb, and percent.plat
#' 
#' Default thresholds and category labels are:
#' | metric       | threshold (%) | label |
#' |--------------+---------------+-------|
#' | percent.mt   | <24           |   <24 |
#' |              | >=24          |   >24 |
#' |              | >= 60         |   >60 |
#' | percent.rb   | < 50          |   <50 |
#' |              | >=50          |   >50 |
#' | percent.hb   | < 3           |    <3 |
#' |              | >=3           |    >3 |
#' | percent.plat | < 1           |    <1 |
#' |              | >=1           |    >1 |
#'
#' 
#' In addition categories are set as factor, and ordered from smallest to
#'  largest, and cell cycle phase as G1, S, G2M.
#' 
categorize_qc <- function(seu) {
    seu@meta.data <- seu@meta.data %>%
        mutate(
            ## convert qc to categories
            mt.cat = ifelse(as.numeric(percent.mt) < 24, "<24",
                     ifelse(as.numeric(percent.mt) >= 60, ">60",
                            ">24")),
            rb.cat = ifelse(as.numeric(percent.rb) >= 50, ">50", "<50"),
            hb.cat = ifelse(as.numeric(percent.hb) >= 3, ">3", "<3"),
            plat.cat = ifelse(as.numeric(percent.plat) >= 1, ">1", "<1"),
            active.cell = NA, ## ifelse(nFeature_RNA > 200 & nFeature_RNA < 6000, TRUE, FALSE)
            log10GenesPerUMI = NA, ## log10(nFeature_RNA) / log10(nCount_RNA),
            ## class change
            ## Phase = factor(Phase, c("G1", "S", "G2M")),
            mt.cat = factor(mt.cat, c("<24", ">24", ">60")),
            rb.cat = factor(rb.cat, c(">50", "<50")),
            hb.cat = factor(hb.cat, c(">3", "<3")),
            plat.cat = factor(plat.cat, c(">1", "<1"))
        )

    
    return(seu)
}

#' ----merging seurat data list by patient id-----------------------------------
#' Function is Palbo+IO clincial trial specific.  Created in response to not
#' being able to merge the whole object.
#'
#' List names should be as "P[0-9]+_[CS][0-9]_.*".  Patient id is captured from
#'  the first element of this pattern e.g. "P[0-9]+".  All matching list
#'  elements that begin with the patient ID will be merged as a single object.
#' 
#' @param s.data.list seurat data object list.  
#'
#' @export
merge_by_patient <- function(s.data.list) {
    
    patient.ids <- names(s.data.list) %>% gsub("(P[0-9]+)_[CS].*", "\\1", .)
    
    patient.list <- lapply(unique(patient.ids), function(i) {
        names(s.data.list)[which(patient.ids == i)]
    })
    names(patient.list) <- unique(patient.ids)
    
    pre.seu.L <- lapply(patient.list, function(i) {
        print(i)
        ##
        if(length(i) == 1) {
            out = s.data.list[[i]]
        }
        if(length(i) == 2) {
            out = merge(s.data.list[[i[1]]],
                        s.data.list[[i[2]]],
                        add.cell.ids = i)
        }
        if(length(i) > 2) {
            out = merge(s.data.list[[i[1]]],
                        s.data.list[i[-1]],
                        add.cell.ids = i)
        }
        return(out)
    })
    
}



#' ----convert_to_monocle, eval=TRUE--------------------------------------------
#'
#' @param seu seurat object
#' @param align logical, align using align_cds
#' @param alignment_group alignment variable
#' @param k  k 
#' @param min_expr minimum expresssion
#' @param round_expr round expression
#' 
#' @export
convert_to_monocle <- function(seu,
                               k = 10, min_expr = 0.01,
                               round_exprs = FALSE,  ...) {
    ## convert seurat object into cds ; cluster cells, and learn_graph
    cds <- as.cell_data_set(seu) %>%
        cluster_cells(k = k) %>%
        learn_graph()
    ## these were missing from the conversion and are required for
    ## monocle to function
    ## detect expressed genes ## by default 0.01
    cds <- detect_genes(cds, min_expr = min_expr)
    ## estimate size factors for samples
    cds <- estimate_size_factors(cds, round_exprs = round_exprs)
    ## row data must have a column named gene_short_name, use rownames
    ## we loose ensembl annotation which is more comprehensive than seurat
    rowData(cds)$gene_short_name <- rownames(rowData(cds))
    ## append cluster and partition to the metadata
    pData(cds)$cluster <- clusters(cds)
    pData(cds)$partition <- partitions(cds)

    return(cds)
}


#' combine doublet finder processed metadata
#'
#' @param meta.data.list list of metadata
#' @param reduced logical, weather to remove the "DF.classifications"
#'   and/or "pANN_" annotations, default FALSE
#' 
#' @param ... optional parameters, unused
#'
#' @importFrom 
#' 
combine_dbf_meta.data <- function(meta.data.list,
                                  reduced = FALSE, ...) {
                                  

    lx <- lapply(meta.data.list,
                 tibble::rownames_to_column,
                 var = "barcode")

    out <- plyr::rbind.fill(lx) %>%
        transform(row.names = paste(orig.ident, barcode, sep = "_"))
    
    
    if(reduced) {
        out <- out %>%
            select(-contains("DF.classifications"), -contains("pANN_"))
    }
    
    return(out)

}

## simple remove
remove_dbf_extended <- function(so, ... ) {

    out <- so@meta.data %>%
        select(-contains("DF.classifications"),
               -contains("pANN_"))
    
    so@meta.data <- out

    return(so)

}



## =============================================================================
##
## velocyto/loom functions
## 
## =============================================================================

#' read loom files
#'
#' @export
read.m.loom <- function(loom.files, sample.names = NULL) {
    message(date(), "  start")

    if(is.null(sample.names)) {
        sample.names <- names(loom.files)
    }
    
    ldat <- lapply(loom.files, ReadVelocity)

    if(is.null(names(ldat))) {
        names(ldat) <- sample.names
    }
    
    message(date(), "  end")
    
    return(ldat)
}


#'
#' 
prepare_loom <- function(loom.data, sample.names =NULL, merge = FALSE, nbin = 24) {
    message(date(), "  start")

    if(is.null(sample.names)) {
        sample.names <- names(loom.data)
    }
    
    ## coerce to seurat
    pre.seu <- lapply(loom.data, function(ldat) {
        pso <- as.Seurat(ldat)
        pso$orig.ident <- gsub("(.*):[ATCG]+x", "\\1", colnames(pso))

        pso <- parse.original.identity(pso) %>%
            estimate_qc_metrics(nbin = nbin) %>%
            categorize_qc()        
        
        return(pso)
    })
    names(pre.seu) <- sample.names
    
    
    ## assert lengths match
    nL <- length(loom.data)
    assertthat::assert_that(nL == length(loom.data))
    assertthat::assert_that(nL == length(pre.seu))
    message(date(), "  end")
    
    if(merge) {
        ## merge data
        ## estimate percentage of mitochondrial expression
        ## cell cycle scoring and prediction
        pre.seu <- merge(pre.seu[[1]], pre.seu[2:nL])
    }
    
    return(pre.seu)
} ## end prepare_loom




#' pass cds metadata to seurat object
#'
#' @export
pass_cds_metadata <- function(so, cds) {

    seu_columns <- names(so@meta.data)
    cds_columns <- names(pData(cds))

    columns_to_pass <- setdiff(cds_columns, seu_columns) %>%
        grep("ident", ., invert = TRUE)

    cds_new <- data.frame(pData(cds)[,columns_to_pass])

    so <- AddMetaData(so, metadata = cds_new)
    
    return(so)
}
