library(ggplot2)
library(ggpubr)
library(cowplot)
library(latex2exp)
library(xtable)
library(rlang)
library(tidyverse)

### GLOBAL PLOTTING SETTINGS
ideal_order <- c(#"BDML-R2D2", 
                 "BDML-IW-HP",
                 "BDML-LKJ-HP", 
                 "BDML-IW",
                 "BDML-LKJ", 
                 "Linero", 
                 "HCPH", 
                 "Naive", 
                 "FDML-Full", 
                 # "FDML-Split",
                 "FDML-XFit",
                 "FDML-Alt",
                 "OLS",
                 "Oracle")
shape_values <- c(19,
                  18, 
                  17, 
                  15, 
                  1, 
                  2, 
                  0, 
                  3, 
                  8,
                  7,
                  NA,
                  NA)
color_values <- c("firebrick4", "#F8766D", "darkorange2", "orange", "#00BA38", "green", 
                  "#619CFF", "steelblue4", "purple", "#F564E3", grey(0.5),"black")

style_mapping <- tibble(Method = ideal_order, shape = shape_values, color = color_values)

### AUXILIARY FUNCTIONS

make_results_table <- function(results){
  results_table <- results %>% 
  group_by(Method, R_Y2, R_D2, rho, alpha, n, p) %>% 
  summarise(coverage = mean(catch), 
            rmse = sqrt(mean(squared_error)), 
            width = mean(interval_width))
}

sort_methods <- function(results_table){
  # get the actual order of the methods and sort them according to ideal order
  results_table %>% mutate(Method = factor(Method, levels = ideal_order))
}

### PLOTTING CODE
# Baseline function to generate individual plots
get_individual_plot <- function(results, y_var, y_label, scale_y_log = FALSE, custom_colors = NULL, custom_shapes = NULL) {

  data <- results %>% make_results_table() %>% sort_methods()
  methods_present <- intersect(unique(data$Method), style_mapping$Method)
  
  plot <- data %>%
    ggplot(aes(x = R_Y2, y = .data[[y_var]], color = Method, shape = Method)) +
    geom_point() + geom_line() +
    theme_bw() +
    xlab(TeX("$Partial R^2_Y$")) + ylab(y_label) +
    theme(legend.position = "none") +
    labs(color = "", shape = "") +
    guides(color = guide_legend(nrow = 1, byrow = TRUE),
           shape = guide_legend(nrow = 1, byrow = TRUE)) +
    scale_x_continuous(breaks = unique(data$sigma))

  if (scale_y_log) {
    plot <- plot + scale_y_log10()
  }

  if (!is.null(custom_colors)) {
    plot <- plot + scale_color_manual(values = custom_colors)
  } else {
    plot <- plot + scale_color_manual(values = setNames(style_mapping$color[style_mapping$Method %in% methods_present],
                                                        style_mapping$Method[style_mapping$Method %in% methods_present]))
  }

  if (!is.null(custom_shapes)) {
    plot <- plot + scale_shape_manual(values = custom_shapes)
  } else {
    plot <- plot + scale_shape_manual(values = setNames(style_mapping$shape[style_mapping$Method %in% methods_present],
                                                        style_mapping$Method[style_mapping$Method %in% methods_present]))
  }


  # return also the mapping of the colors and shapes to the methods
  mapping <- data.frame(
    unique_Method = unique(data$Method),
    mapped_colors = unique(ggplot_build(plot)$data[[1]]$colour),
    mapped_shapes = unique(ggplot_build(plot)$data[[1]]$shape)
  )

  return(list(plot = plot, mapping = mapping))
}


# Wrapper function to combine the plots for Coverage, Interval Width, and RMSE
get_combined_plots <- function(results, 
                               save="results/", 
                               datetimespec_tag){
    p_1 <- get_individual_plot(results, "coverage", "Coverage")$plot
    p_2 <- get_individual_plot(results, "width", "Interval Width (log scale)", scale_y_log = TRUE)$plot
    p_3 <- get_individual_plot(results, "rmse", "RMSE (log scale)", scale_y_log = TRUE)$plot

    # add a horizontal line at 0.95 to the coverage plot
    p_1 <- p_1 + geom_hline(yintercept = 0.95, linetype = "dashed", color = "black")
    
    combined_plots <- plot_grid(p_1, p_2, p_3, ncol = 3)
    
    # Create a common legend from one of the plots
    legend <- ggpubr::get_legend(p_1 + theme(legend.position = "bottom",
                                             legend.text = element_text(size = 12)))

    # Combine the plots and the legend in a vertical layout
    final_plot <- plot_grid(combined_plots, legend, ncol = 1, rel_heights = c(1, 0.1))
    
    if (!is.null(save)) {
      # create the directory if it doesn't exist
      if (!dir.exists(save)) {dir.create(save)}
        ggsave(
          filename = file.path(save, paste0("plot_", datetimespec_tag, ".pdf")), final_plot, width = 9, height = 3.5)
    }

    return(final_plot)
}

get_combined_plots_zoom <- function(results, 
                                    save = "results/",
                                    zoom_in = c("BDML-IW-HP", "BDML-LKJ-HP", "Linero"),
                                    datetimespec_tag) {

    # extract the colors and shapes for the methods in zoom_in
    plot_mappings <- get_individual_plot(results, "coverage", "Coverage")$mapping
    zoom_mapping <- plot_mappings %>% 
      filter(unique_Method %in% zoom_in) %>% 
      arrange(factor(unique_Method, levels = zoom_in))
    # get shapes and colors for the methods in zoom_in
    extracted_colors <- zoom_mapping$mapped_colors
    extracted_shapes <- zoom_mapping$mapped_shapes
    
    results_filtered <- results %>% 
      filter(Method %in% zoom_in)
    
    p_1_zoom <- get_individual_plot(results_filtered, "coverage", "Coverage", custom_colors = extracted_colors, custom_shapes = extracted_shapes)$plot
    p_2_zoom <- get_individual_plot(results_filtered, "width", "Interval Width (log scale)", scale_y_log = TRUE, custom_colors = extracted_colors, custom_shapes = extracted_shapes)$plot
    p_3_zoom <- get_individual_plot(results_filtered, "rmse", "RMSE (log scale)", scale_y_log = TRUE, custom_colors = extracted_colors, custom_shapes = extracted_shapes)$plot

    # add a horizontal line at 0.95 to the coverage plot
    p_1_zoom <- p_1_zoom + geom_hline(yintercept = 0.95, linetype = "dashed", color = "black")
    
    combined_plots <- plot_grid(p_1_zoom, p_2_zoom, p_3_zoom, ncol = 3)
    
    # Create a common legend from one of the plots
    legend <- ggpubr::get_legend(p_1_zoom + theme(legend.position = "bottom",
                                                  legend.text = element_text(size = 12)))

    # Combine the plots and the legend in a vertical layout
    final_plot <- plot_grid(combined_plots, legend, ncol = 1, rel_heights = c(1, 0.1))
    
    if (!is.null(save)) {
      # create the directory if it doesn't exist
      if (!dir.exists(save)) {dir.create(save)}
        ggsave(
          filename = file.path(save, paste0("plot_zoom_", datetimespec_tag, ".pdf")), final_plot, width = 9, height = 3.5)
    }

    return(final_plot)
}
