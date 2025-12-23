#' Plot Median Values Across Treatment Cycles
#'
#' This function generates a plot that visualizes the median values of a specified variable 
#' across treatment cycles, grouped by patient and pathway clusters.
#'
#' @param data data frame containing the input data. It must have columns `SCT_snn_res.0.4_swapped`,
#'             `bioID`, `cycle`, and `Pathway`.
#' @param variable_name A string specifying the column name in `data` for which the median values are calculated.
#' @param color_scheme A named vector color map
#' @param title plot title
#' @param ylab y-axis title
#'
#' @return A ggplot object representing the plot of median values across cycles for different patient and pathway clusters.
#' @import dplyr
#' @import ggplot2
#' @importFrom tidyr pivot_wider
#' 
#' @export
plot_median_values <- function(data,
                               variable_name,
                               color_scheme,
                               title,
                               ylab) {
  # Aggregating median values
  Median_per_Cycle_Patient_Cluster <- data %>%
      dplyr::group_by(SCT_snn_res.0.4_swapped, bioID, cycle, Pathway) %>%
      dplyr::summarise(Value = median(.data[[variable_name]], na.rm = TRUE)) %>%
      pivot_wider(names_from = cycle, values_from = Value)
    
  # Plotting
  plot <- ggplot(Median_per_Cycle_Patient_Cluster) +
    geom_segment(aes(
      x = SCT_snn_res.0.4_swapped,
      xend = SCT_snn_res.0.4_swapped,
      y = C1,
      yend = C3
    ),
    color = "grey") +
    geom_point(
      aes(x = SCT_snn_res.0.4_swapped, y = C1),
      color = color_scheme["C1"],
      fill = color_scheme["C1"],
      size = 3,
      shape = 16
    ) +
    geom_point(
      aes(x = SCT_snn_res.0.4_swapped, y = C3),
      color = color_scheme["C3"],
      fill = color_scheme["C3"],
      size = 3,
      shape = 24
    ) +
    coord_flip() +
    theme_linedraw(base_size = 6) +
    theme(
      strip.background = element_blank(),
      axis.text.y = element_text(face = "bold"),
      strip.text = element_text(colour = "black", face = "bold"),
      legend.position = "none"
    ) +
    xlab("") +
    ggtitle(title) +
    facet_grid(bioID ~ Pathway, scales = "free_y") +
    ylab(ylab)
  
  return(plot)
}


#' Function to assign colors based on fold change
#' @param fold_change  DEG fold cha
assign_color <- function(fold_change) {
    if (is.na(fold_change)) {
        return(NA)
    } else if (fold_change >= 2) {
        return("#9E0101")
    } else if (fold_change >= 1) {
        return("#E63535")
    } else if (fold_change >= 0) {
        return("#F78A8A")
    } else if (fold_change >= -1) {
        return("#86D3FA")
    } else if (fold_change >= -2) {
        return("#07A9FA")
    } else {
        return("#1E6284")
    }
}
