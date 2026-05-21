project_dir <- "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls"

plot_saved_tcga_auc <- function(rdata_file, file_base, ylim_use, legend_cex = 0.55) {
  env <- new.env(parent = emptyenv())
  load(rdata_file, envir = env)

  time_points <- env$tcga_auc_results$auc_time_grid
  auc_curve_values <- list(
    PLSDR = env$auc_curves_list_drpls$best_auc_curve,
    ContiCox = env$auc_curves_list_drpcapls$best_auc_curve,
    pCox = env$auc_curves_list_pcox$best_auc_curve,
    `Ridge-Cox` = env$auc_curves_list_ridgecox$best_auc_curve,
    `EN-Cox` = env$auc_curves_list_enetcox$best_auc_curve
  )

  cols <- c(
    "PLSDR" = "#4C78A8",
    "ContiCox" = "#E45756",
    "pCox" = "#72B7B2",
    "Ridge-Cox" = "#B279A2",
    "EN-Cox" = "#79706E"
  )
  ltys <- c(
    "PLSDR" = 1,
    "ContiCox" = 1,
    "pCox" = 1,
    "Ridge-Cox" = 2,
    "EN-Cox" = 3
  )
  lwds <- c(
    "PLSDR" = 2.2,
    "ContiCox" = 2.8,
    "pCox" = 2.2,
    "Ridge-Cox" = 1.8,
    "EN-Cox" = 1.8
  )
  method_names <- names(auc_curve_values)

  draw_plot <- function() {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    par(mar = c(5.1, 5.1, 4.1, 2.1), xpd = FALSE)
    plot(
      time_points,
      auc_curve_values[[1]],
      type = "n",
      ylim = ylim_use,
      xlab = "Time (days)",
      ylab = "Predictive accuracy (AUC)",
      main = "Predictive accuracy (AUC)"
    )
    abline(h = 0.5, lty = 3, col = "grey80")

    for (method in method_names) {
      lines(
        time_points,
        auc_curve_values[[method]],
        col = cols[method],
        lty = ltys[method],
        lwd = lwds[method]
      )
    }

    legend(
      "topright",
      legend = method_names,
      col = cols[method_names],
      lty = ltys[method_names],
      lwd = lwds[method_names],
      bty = "n",
      cex = legend_cex,
      seg.len = 2.4,
      x.intersp = 0.85,
      y.intersp = 1.05,
      inset = c(0.015, 0.02)
    )
  }

  dir.create(dirname(file_base), recursive = TRUE, showWarnings = FALSE)

  grDevices::pdf(paste0(file_base, ".pdf"), width = 6.5, height = 4.2,
                 pointsize = 10, useDingbats = FALSE)
  draw_plot()
  grDevices::dev.off()

  grDevices::png(paste0(file_base, ".png"), width = 6.5, height = 4.2,
                 units = "in", res = 600, pointsize = 10)
  draw_plot()
  grDevices::dev.off()
}

plot_saved_tcga_auc(
  rdata_file = file.path(project_dir, "TCGA_GBM_validation_results_HNSC.RData"),
  file_base = file.path(project_dir, "results", "TCGA_GBM_auc_curve_0517"),
  ylim_use = c(0.48, 0.95)
)

plot_saved_tcga_auc(
  rdata_file = file.path(project_dir, "TCGA_KIRC_validation_results_BRCA.RData"),
  file_base = file.path(project_dir, "results", "TCGA_KIRC_auc_curve_0519"),
  ylim_use = c(0.58, 0.92)
)
