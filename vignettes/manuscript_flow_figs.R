## =============================================================================
#'
#' This code analyzes our patient PBMC flow cytometry data, performes the
#'  differential abundance of PBMC cell types
#'
#' Code contributors: Jasme Lee, MS
#' Department of Epidemiology & Biostatistics, Memorial Sloan Kettering Cancer Center
#' 2025-07-31
#' Manuscript Figs 4B-D (Flow analysis)
#'
#' 
## =============================================================================


## Set up ---------------------------------------
library(tidyverse)
library(patchwork)

## set ggplot theme
theme_set(theme_classic() + 
            theme(text = element_text(size = 12), 
                  axis.text = element_text(color = "black"), 
                  legend.position = 'bottom') )

col_bor = c("PR" = "#00BFC4", "SD" = "#4C86C6", "PD" = "#A3D276")

## Read in data of perc change for selected subpops
## Assumes same directory as script- change file path here.
flow_dat = 
  readRDS("./flow_perc_change_dat.Rds")

#' nrows = n pts*n subpops = 10*63 = 630
#' ncols = 5
#'  Columns are: 
#'  (1) ID; (2) BOR; (3) full subpopulation name, 
#'  (4) shortened/prettified subpop label for plots; 
#'  (5) the percent diff from baseline to on-trt
dim(flow_dat)


###########################################################
## Fig 4B ---------------------
# ordered by the median of each subpop distribution
flow_dat_w_summary = 
  flow_dat %>%
  # first calculate the medians per subpop
  dplyr::group_by(subpopulation) %>% 
  dplyr::mutate(
    perc_diff_median = 
      # There should not be any NAs
      median(perc_diff, na.rm = TRUE)
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    # Order the subpops by median
    # need to make factor so that the subpops are plotted in order 
    subpop_label_ordered = 
      forcats::fct_reorder(subpop_label, perc_diff_median),
    subpop_ordered = 
      forcats::fct_reorder(subpopulation, perc_diff_median)
  )

# Calculate the FDR adjusted p-values to add the * to fig
p_fdr_perc_change =
  flow_dat_w_summary %>% 
  dplyr::group_by(subpopulation) %>% 
  dplyr::summarise(
    p = stats::wilcox.test(perc_diff, 
                           mu = 0, 
                           paired = FALSE)$p.value
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    # FDR correction
    q = stats::p.adjust(p, "fdr")
  )

# Boxplots
fig_4b_boxplots =
  flow_dat_w_summary %>% 
  ggplot() + 
  # 0 line for no change
  geom_hline(yintercept = 0, linetype = "dotted") + 
  geom_boxplot(aes(x = forcats::fct_rev(subpop_label_ordered),
                   y = perc_diff),
               outlier.size = 0.7) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = NULL,
       y = "Percentage change from baseline to on-treatment")

# The * for q values sig
fig_4b_q = 
  p_fdr_perc_change %>% 
  dplyr::left_join(
    .,
    flow_dat_w_summary,
    by = c("subpopulation")
  ) %>%
  dplyr::mutate(
    sig_display = 
      dplyr::case_when(
        q < 0.001 ~ "***",
        q < 0.01 ~ "**",
        q < 0.05 ~ "*",
        .default = ""
      )
  ) %>% 
  ggplot() + 
  geom_text(aes(x = forcats::fct_rev(subpop_label_ordered), 
                y = "P", 
                label = sig_display), 
            size = 3) + 
  theme(axis.text.x = element_text(angle = 90),
        axis.ticks.y = element_blank(), 
        axis.text.y =  element_blank(),
        axis.line.y = element_blank()) + 
  labs(x = NULL, y = NULL)

# put together
fig_4b = 
  fig_4b_boxplots + fig_4b_q + 
  plot_layout(ncol = 1, 
              heights = c(6, 0.1), 
              guides = "collect",
              axes = "collect")

fig_4b

###########################################################
## 4C and 4D ---------------------
# Only the CD8 and CD4 cell states are included 
# Subset to the 4 combos of CD45RA, CCR7 for CD8+s, then CD4+s
diff_dat_cd8 = 
  flow_dat_w_summary %>%  
  dplyr::filter(subpopulation %in%
                  c("CD45RA_NEG_CCR7_NEG_OF_CD4_NEG_CD8_POS",
                    "CD45RA_NEG_CCR7_POS_OF_CD4_NEG_CD8_POS", 
                    "CD45RA_POS_CCR7_NEG_OF_CD4_NEG_CD8_POS",
                    "CD45RA_POS_CCR7_POS_OF_CD4_NEG_CD8_POS"
                  )) %>% 
  # To make the color labels a bit nicer
  dplyr::mutate(
    subpop_label = forcats::fct_relevel(
      stringr::str_remove(subpop_label, "CD8\\+ "),
      "Tn", "Tcm", "Tem", "Temra"
    ),
    # scale the perc_diff for viz
    perc_diff_scale = as.vector(scale(perc_diff, center = FALSE, scale = TRUE))
  ) %>% 
  dplyr::group_by(
    ctms_id
  ) %>% 
  dplyr::mutate(
    # Use change of CD8+ Tn for sorting pts
    sort_perc = 
      perc_diff[subpop_label == "Tn"]
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    # Factor to make x axis show the order 
    ctms_id_ordered = forcats::fct_reorder(ctms_id, sort_perc)
  ) 

# CD4+s
diff_dat_cd4 = 
  flow_dat_w_summary %>%
  dplyr::filter(subpopulation %in%
                  c("CD45RA_NEG_CCR7_NEG_OF_CD4_POS_CD8_NEG",
                    "CD45RA_NEG_CCR7_POS_OF_CD4_POS_CD8_NEG", 
                    "CD45RA_POS_CCR7_NEG_OF_CD4_POS_CD8_NEG",
                    "CD45RA_POS_CCR7_POS_OF_CD4_POS_CD8_NEG"
                  )) %>% 
  # To make the color labels a bit nicer
  dplyr::mutate(
    subpop_label = forcats::fct_relevel(
      stringr::str_remove(subpop_label, "CD4\\+ "),
      "Tn", "Tcm", "Tem", "Temra"
    ),
    # scale the perc_diff for viz
    perc_diff_scale = as.vector(scale(perc_diff, center = FALSE, scale = TRUE)),
    # Retain the order in CD8+s for consistency/easier viz comparison
    ctms_id_ordered = forcats::fct_relevel(ctms_id, 
                                           levels(diff_dat_cd8$ctms_id_ordered))
  ) 

# Plot for CD4+ diff bars
fig_4c_bar = 
  diff_dat_cd4 %>% 
  ggplot() + 
  geom_bar(aes(x = ctms_id_ordered, fill = subpop_label, y = perc_diff_scale), 
           stat = "identity") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 5, "Paired")[2:5]) + 
  scale_y_continuous(breaks = seq(-3, 3, 1), limits = c(-3.2, 3.2)) + 
  scale_color_manual(values = col_bor) +
  labs(subtitle = "CD4+", 
       fill = "T-cell state",
       y = "Scaled percentage change",,
       x = NULL)

# BORs annotations
fig_4c_bor = 
  diff_dat_cd4 %>% 
  ggplot() + 
  geom_point(aes(x = ctms_id_ordered, y = "BOR", color = BOR_TYPE), size = 3) + 
  scale_color_manual(values = col_bor) + 
  theme_minimal() + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(angle = 0),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) + 
  labs(x = NULL, 
       y = "BOR", 
       color = "Best overall response") + 
  theme(legend.direction = "horizontal")

# Put together Fig 4C
fig_4c = 
  fig_4c_bar + fig_4c_bor +
  plot_layout(ncol = 1,
              heights = c(6, 0.5),
              guides = "collect",
              axes = "collect") + 
  theme(legend.direction = "horizontal") 


fig_4c

# Plot for CD8+ bars
fig_4d_bar = 
  diff_dat_cd8 %>% 
  ggplot() + 
  geom_bar(aes(x = ctms_id_ordered, fill = subpop_label, y = perc_diff_scale), 
           stat = "identity") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 5, "Paired")[2:5]) + 
  scale_y_continuous(breaks = seq(-3, 3, 1), limits = c(-3.2, 3.2)) + 
  scale_color_manual(values = col_bor) +
  labs(subtitle = "CD8+", 
       fill = "T-cell state",
       y = "Scaled percentage change",,
       x = NULL)

# BOR labels
fig_4d_bor = 
  diff_dat_cd8 %>% 
  ggplot() + 
  geom_point(aes(x = ctms_id_ordered, y = "BOR", color = BOR_TYPE), size = 3) + 
  scale_color_manual(values = col_bor) + 
  theme_minimal() + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(angle = 0),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) + 
  labs(x = NULL, 
       y = "BOR", 
       color = "Best overall response") + 
  theme(legend.direction = "horizontal")


fig_4d = 
  fig_4d_bar + fig_4d_bor +
  plot_layout(ncol = 1,
              heights = c(6, 0.5),
              guides = "collect",
              axes = "collect") + 
  theme(legend.direction = "horizontal") 


fig_4d
