library(ggplot2)
library(ggpubr)
library(cowplot)
library(latex2exp)
library(xtable)
library(rlang)

### GLOBAL PLOTTING SETTINGS

ideal_order <- c("BDML_r2d2", "BDML_b2", "FDML_full", "FDML_split", "hahn", "linero", "naive")


### AUXILIARY FUNCTIONS

make_results_table <- function(results){
  results_table <- results %>% 
  group_by(Method, setting, sigma, N, P) %>% 
  summarise(coverage = mean(catch), 
            rmse = sqrt(mean(squared_error)), 
            width = mean(interval_width))
}

sort_methods <- function(results){
  # get the actual order of the methods and sort them according to ideal order
  results %>% mutate(Method = factor(Method, levels = ideal_order))
}

### PLOTTING CODE

# Baseline function to generate individual plots
get_individual_plot <- function(results, y_var, y_label, scale_y_log = FALSE, custom_colors = NULL, custom_shapes = NULL) {

  data <- results %>% make_results_table() %>% sort_methods()
  
  plot <- data %>%
    filter(setting != "random") %>%
    ggplot(aes(x = sigma, y = .data[[y_var]], color = Method, shape = Method)) +
    geom_point() + geom_line() +
    theme_bw() +
    xlab(TeX("$\\sigma$ (log scale)")) + ylab(y_label) +
    theme(legend.position = "none") +
    labs(color = "", shape = "") +
    guides(color = guide_legend(nrow = 1, byrow = TRUE),
           shape = guide_legend(nrow = 1, byrow = TRUE)) +
    scale_x_log10(breaks = unique(data$sigma))

  if (scale_y_log) {
    plot <- plot + scale_y_log10()
  }

  if (!is.null(custom_colors)) {
    plot <- plot + scale_color_manual(values = custom_colors)
  }

  if (!is.null(custom_shapes)) {
    plot <- plot + scale_shape_manual(values = custom_shapes)
  }

  return(plot)
}


# Wrapper function to combine the plots for Coverage, Interval Width, and RMSE
get_combined_plots <- function(results, save=TRUE){
    p_1 <- get_individual_plot(results, "coverage", "Coverage")
    p_2 <- get_individual_plot(results, "width", "Interval Width (log scale)", scale_y_log = TRUE)
    p_3 <- get_individual_plot(results, "rmse", "RMSE (log scale)", scale_y_log = TRUE)
    
    combined_plots <- plot_grid(p_1, p_2, p_3, ncol = 3)
    
    # Create a common legend from one of the plots
    legend <- ggpubr::get_legend(p_1 + theme(legend.position = "bottom"))

    # Combine the plots and the legend in a vertical layout
    final_plot <- plot_grid(combined_plots, legend, ncol = 1, rel_heights = c(1, 0.1))
    
    if (save) {
        ggsave(paste0("results/plot_", format(Sys.time(), "%Y%m%d%H"), ".pdf"), final_plot, width = 9, height = 3.5)
    }

    return(final_plot)
}

get_combined_plots_zoom <- function(results, save=TRUE, zoom_in = c("naive", "hahn", "FDML_split")){

    p_1 <- get_individual_plot(results, "coverage", "Coverage")
    # get the factor levels for the methods in zoom_in to extract the colors and shapes
    zoom_in_index <- match(ideal_order, zoom_in) %>% order(na.last = NA)
    extracted_colors <- unique(ggplot_build(p_1)$data[[1]]$colour)[zoom_in_index]
    extracted_shapes <- unique(ggplot_build(p_1)$data[[1]]$shape)[zoom_in_index]

    results_filtered <- results %>% filter(Method %in% zoom_in)

    p_1_zoom <- get_individual_plot(results_filtered, "coverage", "Coverage", custom_colors = extracted_colors, custom_shapes = extracted_shapes)
    p_2_zoom <- get_individual_plot(results_filtered, "width", "Interval Width (log scale)", scale_y_log = TRUE, custom_colors = extracted_colors, custom_shapes = extracted_shapes)
    p_3_zoom <- get_individual_plot(results_filtered, "rmse", "RMSE (log scale)", scale_y_log = TRUE, custom_colors = extracted_colors, custom_shapes = extracted_shapes)
    
    combined_plots <- plot_grid(p_1_zoom, p_2_zoom, p_3_zoom, ncol = 3)
    
    # Create a common legend from one of the plots
    legend <- ggpubr::get_legend(p_1_zoom + theme(legend.position = "bottom"))

    # Combine the plots and the legend in a vertical layout
    final_plot <- plot_grid(combined_plots, legend, ncol = 1, rel_heights = c(1, 0.1))
    
    if (save) {
        ggsave(paste0("results/plot_zoom_", format(Sys.time(), "%Y%m%d%H"), ".pdf"), final_plot, width = 9, height = 3.5)
    }

    return(final_plot)
}
