project_dir <- "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls"
setwd(project_dir)

script_lines <- readLines(file.path(project_dir, "simulation", "search_drpls_gap_aggressive.R"), warn = FALSE)
settings_start <- grep("^settings <- list", script_lines)[1L]
eval(parse(text = script_lines[seq_len(settings_start - 1L)]), envir = .GlobalEnv)

setting <- make_setting(
  "var_signal_sparse10_scale10_n60_p3000_g065_c50",
  "simulate_variance_signal",
  list(n_samples = 60, n_total_genes = 3000,
       signal_genes = 10, signal_scale = 10.0,
       signal_rho = 0.95, signal_strength = 0.65,
       decoy_blocks = 12, decoy_genes = 50,
       censoring = 0.50),
  406L
)

cat("\n--- Confirming balanced DRPLS-gap candidate with R = 10:", setting$name, "---\n")
flush.console()
ans <- run_core_setting(setting, R = 10L, n_cores = 10L, seed = 18100)
print(ans, digits = 4)

