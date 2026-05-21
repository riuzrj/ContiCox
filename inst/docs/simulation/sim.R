source("/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls/simulation/simulation_gaussian_latent_factor.R")


res_I <- run_gaussian_scenario(
  scenario = "I",
  R = 10,
  n_cores = 10,
  n_samples = 100,
  n_total_genes = 2000,
  signal_strength = 1.0,
  censoring = 0.1,
  auc_method = "NNE"
)

res_I$paired_difference_plots <- plot_paired_difference(res_I$mc_res_par)

res_I$mc_res_par$draws
res_I_stats <- summarize_mc_iAUC(res_I$mc_res_par)
res_I$selected_alpha_stats$summary
res_I$selected_alpha_stats$frequency
res_I$selected_alpha_plot

# supplementary
res_I$mc_res_par_plots$violin
res_I$mc_res_par_plots$boxplot

# main
res_I$mc_res_main_plots <- plot_main_iAUC(
  res_I$mc_res_par,
  y_limits = c(0.5, 1.0),
  y_breaks = seq(0.5, 1.0, by = 0.1)
)
res_I$mc_res_main_plots$violin
res_I$mc_res_main_plots$boxplot
res_I$mc_res_par_stats$paired_difference
res_I$paired_difference_plots$forest
res_I$paired_difference_plots$violin



saveRDS(res_I, file = "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls/simulation/sim_res_I")


res_II <- run_gaussian_scenario(
  scenario = "II",
  R = 10,
  n_cores = 10,
  n_samples = 150,
  n_total_genes = 1000,
  n_factors = 12,
  lowvar_index = 10,
  signal_strength = 1.0,
  censoring = 0.1,
  auc_method = "NNE"
)

res_II$paired_difference_plots <- plot_paired_difference(res_II$mc_res_par)

res_II$mc_res_par$draws
res_II_stats <- summarize_mc_iAUC(res_II$mc_res_par)
res_II$selected_alpha_stats$summary
res_II$selected_alpha_stats$frequency
res_II$selected_alpha_plot

ablation_II <- res_II$mc_res_par$draws[, c("CPPC_alpha0", "CPPC_alpha1", "DRPCA_PLS")]

colMeans(ablation_II, na.rm = TRUE)


# supplementary
res_II$mc_res_par_plots$violin
res_II$mc_res_par_plots$boxplot

# main
res_II$mc_res_main_plots <- plot_main_iAUC(
  res_II$mc_res_par,
  y_limits = c(0.3, 1.0),
  y_breaks = seq(0.5, 1.0, by = 0.1)
)
res_II$mc_res_main_plots$violin
res_II$mc_res_main_plots$boxplot
res_II$mc_res_par_stats$paired_difference
res_II$paired_difference_plots$forest
res_II$paired_difference_plots$violin

saveRDS(res_II, file = "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls/simulation/sim_res_II")


res_III <- run_gaussian_scenario(
  scenario = "III",
  R = 10,
  n_cores = 10,
  n_samples = 100,
  n_total_genes = 1000,
  signal_strength = 1.0,
  censoring = 0.4,
  auc_method = "NNE"
)

res_III$paired_difference_plots <- plot_paired_difference(res_III$mc_res_par)

res_III$mc_res_par$draws
res_III_stats <- summarize_mc_iAUC(res_III$mc_res_par)
res_III$selected_alpha_stats$summary
res_III$selected_alpha_stats$frequency
res_III$selected_alpha_plot

# supplementary
res_III$mc_res_par_plots$violin
res_III$mc_res_par_plots$boxplot

# main
res_III$mc_res_main_plots <- plot_main_iAUC(
  res_III$mc_res_par,
  y_limits = c(0.5, 1.0),
  y_breaks = seq(0.5, 1.0, by = 0.1)
)
res_III$mc_res_main_plots$violin
res_III$mc_res_main_plots$boxplot
res_III$mc_res_par_stats$paired_difference
res_III$paired_difference_plots$forest
res_III$paired_difference_plots$violin


saveRDS(res_III, file = "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls/simulation/sim_res_III")

res_I <- readRDS("/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls/simulation/sim_res_I")
res_II <- readRDS("/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls/simulation/sim_res_II")
res_III <- readRDS("/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls/simulation/sim_res_III")


export_main_iAUC_boxplots <- function(
    res_I,
    res_II,
    res_III,
    eigengene_plot,
    output_dir = "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls/simulation/figures",
    width = 4.8,
    height = 4.2,
    dpi = 300
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required to export plots.")
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  plots <- list(
    scenario_I = res_I$mc_res_main_plots$boxplot,
    scenario_II = res_II$mc_res_main_plots$boxplot,
    scenario_III = res_III$mc_res_main_plots$boxplot,
    scenario_IV_eigengene = eigengene_plot
  )
  
  for (nm in names(plots)) {
    ggplot2::ggsave(
      filename = file.path(output_dir, paste0(nm, "_iAUC_boxplot.pdf")),
      plot = plots[[nm]],
      width = width,
      height = height,
      units = "in",
      device = "pdf",
      limitsize = FALSE
    )
    ggplot2::ggsave(
      filename = file.path(output_dir, paste0(nm, "_iAUC_boxplot.png")),
      plot = plots[[nm]],
      width = width,
      height = height,
      units = "in",
      device = "png",
      dpi = dpi,
      bg = "white",
      limitsize = FALSE
    )
  }
  
  invisible(file.path(output_dir, paste0(names(plots), "_iAUC_boxplot.pdf")))
}

# Run after the eigengene main plot has been generated in simulation_eigengene.R:
exported_main_figures <- export_main_iAUC_boxplots(
  res_I = res_I,
  res_II = res_II,
  res_III = res_III,
  eigengene_plot = mc_res_main_plots$boxplot
)
