#' z value estimation for figure 1B
#' @param pfs PFS
#' @export
z.value <- function(pfs) {
    mu <- mean(pfs, na.rm = TRUE)
    s <- sd(pfs, na.rm = TRUE)

    z <- (pfs - mu) / s

    return(z)
}


#' gene effect plot grouped by cycle (or other paired that can be paired)
#'
#' @param so seurat object
#' @param gene gene of interst (only 1)
#' @param group grouping variable
#' @param comparison groups to compare
gene_cycle_effect_plot_paired <- function(so, gene, group = "cycle",
                                          comparison = list(c("C1", "C3"))) {
    dat <- FetchData(so, vars = c("paired.id", "curated.celltype.l1", group, gene))
    names(dat)[3:4] <- c("group", "gene")
    
    dat <- dat %>%
        filter(group %in% comparison[[1]],
               curated.celltype.l1 == "Cancer Cells")
    
    smry <- summarySE(dat, measurevar = "gene",
                      groupvars = c("paired.id", "group"))
    
    ## The errorbars overlapped, so use position_dodge to move them horizontally
    pd <- position_dodge(0.1) # move them .05 to the left and right

    p1 <- smry %>%
        ggplot(aes(x = group, y = gene,
                   shape = paired.id,
                   group = paired.id,
                   colour = paired.id)) + 
        geom_errorbar(aes(ymin = gene - se,
                          ymax = gene + se),
                      width = 0.1, position = pd) +
        ylab(label = gene) +
        geom_line(position = pd) + 
        geom_point(position = pd)
    
    return(p1)
}

#' gene violin plot grouped by cycle (or other paired that can be paired)
#'
#' @param so seurat object
#' @param gene gene of interst (only 1)
#' @param group grouping variable
#' @param comparison groups to compare
gene_cycle_violin_plot_paired <- function(so, gene, group = "cycle",
                                          comparison = list(c("C1", "C3")),
                                          dodge = 0.8) {

    dat <- FetchData(so, vars = c("paired.id", "curated.celltype.l1", group, gene))
    names(dat)[3:4] <- c("group", "gene")
    
    dat <- dat %>%
        filter(group %in% comparison[[1]],
               curated.celltype.l1 == "Cancer Cells")
    
    smry <- summarySE(dat, measurevar = "gene",
                      groupvars = c("paired.id", "group"))
    
    ## The errorbars overlapped, so use position_dodge to move them horizontally
    pd <- position_dodge(width = dodge) # move them .05 to the left and right

    p1 <- dat %>%
        ggplot(aes(x = group, y = gene,
                   fill = paired.id)) +
        geom_violin(
            trim = FALSE,
            position = pd,
            linewidth = 0.1) +
        geom_boxplot(width = 0.08,
                     linewidth = 0.1,
                     position = pd,
                     ##colour = "black",
                     ##fill = "white",
                     outlier.shape = NA,
                     outlier.color = "#00000088") +
        ylab(label = gene)  +
        stat_compare_means(comparisons = comparison,
                           method = "wilcox.test",
                           label = "p.signif",
                           bracket.size = .1) 
##        stat_summary(fun.data = give.n0, geom = "text",
##                     position = pd,
##                     angle = 0, 
##                     size = 1.4)
    
    return(p1)
}

gene_cycle_by_patient_violin_paired_plot <- function(so, gene,
                                                     group = "patient.id",
                                                     comparison = list(c("C1", "C3")),
                                                     dodge = 0.8) {

    dat <- FetchData(so, vars = unique(c("paired.id", "curated.celltype.l1",
                                  "cycle", group, gene)))
                                  
    names(dat)[4:5] <- c("group", "gene")
    
    dat <- dat %>%
        filter(group %in% comparison[[1]],
               curated.celltype.l1 == "Cancer Cells")
    
    smry <- summarySE(dat, measurevar = "gene",
                      groupvars = c("group", "cycle"))
    
    ## The errorbars overlapped, so use position_dodge to move them horizontally
    pd <- position_dodge(width = dodge) # move them .05 to the left and right

    p1 <- dat %>%
        ggplot(aes(x = group, y = gene,
                   fill = cycle)) +
        geom_violin(
            trim = FALSE,
            position = pd,
            linewidth = 0.1) +
        geom_boxplot(width = 0.08,
                     linewidth = 0.1,
                     position = pd,
                     ##colour = "black",
                     ##fill = "white",
                     outlier.shape = NA,
                     outlier.color = "#00000088") +
        ylab(label = gene)  +
##        stat_summary(fun.data = give.n0, geom = "text",
##                     position = pd,
##                     angle = 0, 
##                     size = 1.4)
    
    return(p1)

}
    

## Fetch Best Overall Response data
bor.data.fetch <- function(so, list.x) {
    bor.data.df <- FetchData(p11044x,
                             vars = c("best.overall.response", "curated.celltype.l1",
                                      "cycle", list.x)) %>%
        pivot_longer(cols = dx.list, names_to = "gene", values_to = "Expression") %>%
        filter(cycle == "C1",
               best.overall.response != "NE") %>%           
        mutate(gene = factor(gene, list.x),
               best.overall.response = factor(best.overall.response,
                                              c("SD", "PD")))
    bor.data.df <- bor.data.df %>%
        mutate(celltype_x_gene = paste(curated.celltype.l1, gene, sep = ": "))

    bor.data.L <- split(bor.data.df, f = bor.data.df$celltype_x_gene)

    return(bor.data.L)
}


