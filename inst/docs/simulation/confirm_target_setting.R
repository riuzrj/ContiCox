project_dir <- "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls"
setwd(project_dir)

script_lines <- readLines(file.path(project_dir, "simulation", "search_target_gap_settings.R"), warn = FALSE)
settings_start <- grep("^settings <- list", script_lines)[1L]
eval(parse(text = script_lines[seq_len(settings_start - 1L)]), envir = .GlobalEnv)

settings <- list(
  make_setting(
    "two_source_n80_p2000_balanced_hard",
    "simulate_two_source_factor",
    list(n_samples = 80, n_total_genes = 2000,
         stable_genes = 60, supervised_genes = 12,
         stable_rho = 0.92, supervised_rho = 0.30,
         signal_strength = 0.75, stable_weight = 0.70,
         supervised_weight = 0.70, n_decoy_blocks = 12,
         decoy_genes = 40, decoy_rho = 0.90,
         censoring = 0.50),
    201L
  ),
  make_setting(
    "block_n80_p1500_r09_g08_c45",
    "simulate_block_gaussian",
    list(n_samples = 80, n_total_genes = 1500, rho = 0.9,
         signal_strength = 0.8, censoring = 0.45,
         signal_signs = c(1, -1)),
    202L
  )
)

for (setting in settings) {
  cat("\n--- Confirming candidate with R = 10:", setting$name, "---\n")
  flush.console()
  ans <- run_core_setting(setting, R = 10L, n_cores = 10L, seed = 13100)
  print(ans, digits = 4)
}

