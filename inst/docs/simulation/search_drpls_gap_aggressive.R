project_dir <- "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls"
setwd(project_dir)

script_lines <- readLines(file.path(project_dir, "simulation", "search_drpls_gap_settings.R"), warn = FALSE)
settings_start <- grep("^settings <- list", script_lines)[1L]
eval(parse(text = script_lines[seq_len(settings_start - 1L)]), envir = .GlobalEnv)

settings <- list(
  make_setting(
    "var_signal_sparse5_scale8_n80_p3000_g07_c45",
    "simulate_variance_signal",
    list(n_samples = 80, n_total_genes = 3000,
         signal_genes = 5, signal_scale = 8.0,
         signal_rho = 0.95, signal_strength = 0.7,
         decoy_blocks = 12, decoy_genes = 50,
         censoring = 0.45),
    401L
  ),
  make_setting(
    "var_signal_sparse8_scale8_n80_p3000_g07_c45",
    "simulate_variance_signal",
    list(n_samples = 80, n_total_genes = 3000,
         signal_genes = 8, signal_scale = 8.0,
         signal_rho = 0.95, signal_strength = 0.7,
         decoy_blocks = 12, decoy_genes = 50,
         censoring = 0.45),
    402L
  ),
  make_setting(
    "var_signal_sparse10_scale6_n80_p3000_g075_c45",
    "simulate_variance_signal",
    list(n_samples = 80, n_total_genes = 3000,
         signal_genes = 10, signal_scale = 6.0,
         signal_rho = 0.95, signal_strength = 0.75,
         decoy_blocks = 12, decoy_genes = 50,
         censoring = 0.45),
    403L
  ),
  make_setting(
    "var_signal_sparse8_scale10_n70_p3000_g07_c50",
    "simulate_variance_signal",
    list(n_samples = 70, n_total_genes = 3000,
         signal_genes = 8, signal_scale = 10.0,
         signal_rho = 0.95, signal_strength = 0.7,
         decoy_blocks = 12, decoy_genes = 50,
         censoring = 0.50),
    404L
  ),
  make_setting(
    "var_signal_sparse12_scale8_n70_p4000_g07_c50",
    "simulate_variance_signal",
    list(n_samples = 70, n_total_genes = 4000,
         signal_genes = 12, signal_scale = 8.0,
         signal_rho = 0.95, signal_strength = 0.7,
         decoy_blocks = 14, decoy_genes = 60,
         censoring = 0.50),
    405L
  ),
  make_setting(
    "var_signal_sparse10_scale10_n60_p3000_g065_c50",
    "simulate_variance_signal",
    list(n_samples = 60, n_total_genes = 3000,
         signal_genes = 10, signal_scale = 10.0,
         signal_rho = 0.95, signal_strength = 0.65,
         decoy_blocks = 12, decoy_genes = 50,
         censoring = 0.50),
    406L
  )
)

summaries <- data.frame()
for (setting in settings) {
  cat("\n--- Aggressive DRPLS-gap coarse screening:", setting$name, "---\n")
  flush.console()
  ans <- run_core_setting(setting, R = 5L, n_cores = 5L, seed = 17100)
  summaries <- rbind(summaries, ans)
  print(ans, digits = 4)
  flush.console()
}

cat("\n=== Aggressive DRPLS-gap coarse screening summary ===\n")
summaries <- summaries[order(-summaries$diff_vs_DRPLS), ]
print(summaries, digits = 4)

best <- settings[[match(summaries$setting[1L], vapply(settings, `[[`, "", "name"))]]
cat("\n--- Confirming best aggressive DRPLS-gap setting with R = 10:", best$name, "---\n")
confirm <- run_core_setting(best, R = 10L, n_cores = 10L, seed = 17900)
print(confirm, digits = 4)

