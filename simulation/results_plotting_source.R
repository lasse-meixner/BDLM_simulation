library(ggplot2)
library(ggpubr)
library(cowplot)
library(purrr)
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
                 "BDML-IW-JS-MAT",
                 "BDML-IW-JS-I",
                 "Linero", 
                 "HCPH", 
                 "Naive", 
                 "FDML-Full",
                 # Doesnt include "FDML-Split" currently 
                 "FDML-XFit",
                 "FDML-Alt",
                 "OLS",
                 "Oracle")
shape_values <- c(19,
                  18, 
                  17, 
                  15, 
                  #16,
                  10,
                  12,
                  1, 
                  2, 
                  0, 
                  3, 
                  8,
                  7,
                  4,
                  5)
color_values <- c("firebrick4", 
                  "#F8766D",
                  "darkorange2",
                  "orange",
                  "black",
                  "grey", 
                  "#00BA38",
                  "green",
                  "#619CFF",
                  "steelblue4",
                  "purple",
                  "#F564E3",
                  "lightgrey",
                  "black")

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
    xlab(TeX("$Partial R^2_Y$")) + 
    ylab(y_label) +
    theme(legend.position = "none") +
    labs(color = "", shape = "") +
    guides(color = guide_legend(nrow = 1, byrow = TRUE),
           shape = guide_legend(nrow = 1, byrow = TRUE)) +
    scale_x_continuous(breaks = unique(data$R_Y2))

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
                               suffix = ""){
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
          filename = file.path(save, paste0("plot_", suffix, ".pdf")), final_plot, width = 9, height = 3.5)
    }

    return(final_plot)
}

get_combined_plots_zoom <- function(results, 
                                    save = "results/",
                                    zoom_in = c("BDML-IW-HP", "BDML-LKJ-HP", "Linero"),
                                    suffix = "") {

    # extract the colors and shapes for the methods in zoom_in
    plot_mappings <- get_individual_plot(results, "coverage", "Coverage")$mapping
    zoom_mapping <- plot_mappings %>% 
      filter(unique_Method %in% zoom_in) %>% 
      arrange(factor(unique_Method, levels = zoom_in))
  
    # name the colors & shapes by the method so ggplot matches by name
    extracted_colors <- setNames(
      zoom_mapping$mapped_colors,
      zoom_mapping$unique_Method
    )
    extracted_shapes <- setNames(
      zoom_mapping$mapped_shapes,
      zoom_mapping$unique_Method
    )


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
          filename = file.path(save, paste0("plot_zoom_", suffix, ".pdf")), final_plot, width = 9, height = 3.5)
    }

    return(final_plot)
}

# ---------------------------------------------------------------------
#' Generate and save “stacked” R_D2 panels for each (rho, alpha),
#' using only a specified list of zoom_in methods.
#'
#' @param results        data.frame of your combined results
#' @param datetime_tag   character; subdirectory (and filename) tag,
#'                       e.g. "20250613-1430" (defaults to current time)
#' @param save_dir       character; top‐level results directory
#'                       (defaults to "results")
#' @param zoom_in        character vector of Methods to include in each panel
#'                       (defaults to c("BDML-LKJ-HP","BDML-IW-HP","BDML-IW","Linero","Oracle"))
#' @export
generate_stacked_plots <- function(results,
                                   datetime_tag = format(Sys.time(), "%Y%m%d-%H%M"),
                                   save_dir     = "results",
                                   zoom_in      = c("BDML-LKJ",
                                                    "BDML-LKJ-HP",
                                                    "BDML-IW",
                                                    "BDML-IW-HP",
                                                    "BDML-IW-JS-MAT",
                                                    "BDML-IW-JS-I",
                                                    "BLRs-baseline",
                                                    "BLRs-FDML",
                                                    "BLRs-OLS-oracle")) {
  # create output folder
  out_dir <- file.path(save_dir, datetime_tag)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # grid over (rho, alpha)
  grid_ra <- expand.grid(
    rho   = unique(results$rho),
    alpha = unique(results$alpha)
  )
  rd2_vals <- unique(results$R_D2)

  for (i in seq_len(nrow(grid_ra))) {
    rho_i   <- grid_ra$rho[i]
    alpha_i <- grid_ra$alpha[i]

    # build one “row” per R_D2
    make_row <- function(r_d2) {
      subres <- results %>%
        filter(R_D2   == r_d2,
               rho    == rho_i,
               alpha  == alpha_i,
               Method %in% zoom_in)

      # 3-panel without legend
      row_plot <- get_combined_plots(subres, save = NULL) +
        theme(legend.position = "none")

      # label strip “R_D2 = …”
      label_strip <- ggdraw() +
        draw_label(
          paste0("R_D2 = ", r_d2),
          fontface = "bold", x = 0, hjust = 0, size = 12
        )

      plot_grid(label_strip, row_plot, ncol = 1,
                rel_heights = c(0.05, 0.95))
    }

    # stack them all
    all_rows    <- purrr::map(rd2_vals, make_row)
    stacked_rows <- plot_grid(plotlist = all_rows, ncol = 1)

    # extract a shared legend from a dummy coverage plot
    dummy_sub <- results %>%
      filter(R_D2   == rd2_vals[1],
             rho    == rho_i,
             alpha  == alpha_i,
             Method %in% zoom_in)

    dummy_cov <- get_individual_plot(dummy_sub, "coverage", "Coverage")$plot +
      geom_hline(yintercept = 0.95, linetype = "dashed") +
      theme(legend.position = "bottom",
            legend.text     = element_text(size = 10))

    shared_legend <- ggpubr::get_legend(dummy_cov)

    # combine panels + legend
    final_plot <- plot_grid(
      stacked_rows, shared_legend,
      ncol        = 1,
      rel_heights = c(1, 0.08)
    )

    # save
    out_name <- paste0("stacked_R_D2_",
                       datetime_tag,
                       "_rho_",   rho_i,
                       "_alpha_", alpha_i,
                       ".pdf")

    ggsave(
      filename = file.path(out_dir, out_name),
      plot     = final_plot,
      width    = 9,
      height   = 3.5 * length(rd2_vals)
    )
  }
}
