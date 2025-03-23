library(ggplot2)
library(ggpubr)
library(cowplot)
library(latex2exp)
library(xtable)
library(rlang)
library(tidyverse)

### GLOBAL PLOTTING SETTINGS
ideal_order <- c(#"BDML-R2D2", 
                 "BDML-HP-IW",
                 "BDML-HP-LKJ", 
                 "BDML-IW",
                 "BDML-LKJ", 
                 "Linero", 
                 "HCPH", 
                 "Naive", 
                 "FDML-Full", 
                 "FDML-Split")
shape_values <- c(19,
                  18, 
                  17, 
                  15, 
                  1, 
                  2, 
                  0, 
                  3, 
                  8)
color_values <- c("firebrick4", "#F8766D", "darkorange2", "orange", "#00BA38", "green", "#619CFF", "purple", "#F564E3")

style_mapping <- tibble(Method = ideal_order, shape = shape_values, color = color_values)

### AUXILIARY FUNCTIONS

make_results_table <- function(results){
  results_table <- results %>% 
  group_by(Method, setting, sigma, N, P) %>% 
  summarise(coverage = mean(catch), 
            rmse = sqrt(mean(squared_error)), 
            width = mean(interval_width)) %>%
  mutate(Method = case_when(
    Method == "BDML_r2d2" ~ "BDML-R2D2",
    Method == "BDML_b2_iw" ~ "BDML-HP-IW",
    Method == "BDML_b" ~ "BDML-LKJ",
    Method == "BDML_b2" ~ "BDML-HP-LKJ",
    Method == "BDML_iw" ~ "BDML-IW",
    Method == "FDML_full" ~ "FDML-Full",
    Method == "FDML_split" ~ "FDML-Split",
    Method == "hahn" ~ "HCPH",
    Method == "linero" ~ "Linero",
    Method == "naive" ~ "Naive",
    TRUE ~ Method  # keep other values unchanged
  ))
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
get_combined_plots <- function(results, save=TRUE){
    p_1 <- get_individual_plot(results, "coverage", "Coverage")$plot
    p_2 <- get_individual_plot(results, "width", "Interval Width (log scale)", scale_y_log = TRUE)$plot
    p_3 <- get_individual_plot(results, "rmse", "RMSE (log scale)", scale_y_log = TRUE)$plot

    # add a horizontal line at 0.95 to the coverage plot
    p_1 <- p_1 + geom_hline(yintercept = 0.95, linetype = "dashed", color = "black")
    
    combined_plots <- plot_grid(p_1, p_2, p_3, ncol = 3)
    
    # Create a common legend from one of the plots
    legend <- ggpubr::get_legend(p_1 + theme(legend.position = "bottom"))

    # Combine the plots and the legend in a vertical layout
    final_plot <- plot_grid(combined_plots, legend, ncol = 1, rel_heights = c(1, 0.1))
    
    if (save) {
        ggsave(paste0("results/plot_", format(Sys.time(), "%Y%m%d-%H%M"), ".pdf"), final_plot, width = 9, height = 3.5)
    }

    return(final_plot)
}

# Wrapper function to combine the plots for Coverage, Interval Width, and RMSE for a subset of methods
get_combined_plots_zoom <- function(results, save=TRUE, zoom_in = c("BDML-R2D2", "BDML-Hier", "BDML-Basic", "Linero"), suffix = ""){

    # extract the colors and shapes for the methods in zoom_in
    zoom_methods <- intersect(zoom_in, style_mapping$Method)
    extracted_colors <- setNames(style_mapping$color, style_mapping$Method)[zoom_methods]
    extracted_shapes <- setNames(style_mapping$shape, style_mapping$Method)[zoom_methods]


    results_filtered <- results %>% 
        # cheeky double rename to be able to filter using zoom_in
        mutate(
          Method = case_when(
          Method == "BDML_r2d2" ~ "BDML-R2D2",
          Method == "BDML_b2_iw" ~ "BDML-IW-Hier",
          Method == "BDML_b" ~ "BDML-Basic",
          Method == "BDML_b2" ~ "BDML-Hier",
          Method == "BDML_iw" ~ "BDML-IW",
          Method == "FDML_full" ~ "FDML-Full",
          Method == "FDML_split" ~ "FDML-Split",
          Method == "hahn" ~ "HCPH",
          Method == "linero" ~ "Linero",
          Method == "naive" ~ "Naive",
          TRUE ~ Method  # keep other values unchanged
          )) %>%
        filter(Method %in% zoom_in)

    p_1_zoom <- get_individual_plot(results_filtered, "coverage", "Coverage", custom_colors = extracted_colors, custom_shapes = extracted_shapes)$plot
    p_2_zoom <- get_individual_plot(results_filtered, "width", "Interval Width (log scale)", scale_y_log = TRUE, custom_colors = extracted_colors, custom_shapes = extracted_shapes)$plot
    p_3_zoom <- get_individual_plot(results_filtered, "rmse", "RMSE (log scale)", scale_y_log = TRUE, custom_colors = extracted_colors, custom_shapes = extracted_shapes)$plot

    # add a horizontal line at 0.95 to the coverage plot
    p_1_zoom <- p_1_zoom + geom_hline(yintercept = 0.95, linetype = "dashed", color = "black")
    
    combined_plots <- plot_grid(p_1_zoom, p_2_zoom, p_3_zoom, ncol = 3)
    
    # Create a common legend from one of the plots
    legend <- ggpubr::get_legend(p_1_zoom + theme(legend.position = "bottom"))

    # Combine the plots and the legend in a vertical layout
    final_plot <- plot_grid(combined_plots, legend, ncol = 1, rel_heights = c(1, 0.1))
    
    if (save) {
        ggsave(paste0("results/plot_zoom_", format(Sys.time(), "%Y%m%d-%H%M"), suffix, ".pdf"), final_plot, width = 9, height = 3.5)
    }

    return(final_plot)
}
