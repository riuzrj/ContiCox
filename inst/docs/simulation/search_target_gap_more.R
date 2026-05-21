project_dir <- "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls"
setwd(project_dir)

script_lines <- readLines(file.path(project_dir, "simulation", "search_target_gap_settings.R"), warn = FALSE)
settings_start <- grep("^settings <- list", script_lines)[1L]
eval(parse(text = script_lines[seq_len(settings_start - 1L)]), envir = .GlobalEnv)

settings <- list(
  make_setting(
    "block_n70_p2000_r09_g07_c50",
    "simulate_block_gaussian",
    list(n_samples = 70, n_total_genes = 2000, rho = 0.9,
         signal_strength = 0.7, censoring = 0.50,
         signal_signs = c(1, -1)),
    101L
  ),
  make_setting(
    "block_n80_p3000_r09_g07_c45",
    "simulate_block_gaussian",
    list(n_samples = 80, n_total_genes = 3000, rho = 0.9,
         signal_strength = 0.7, censoring = 0.45,
         signal_signs = c(1, -1)),
    102L
  ),
  make_setting(
    "block_n80_p1500_r09_g08_bigblock",
    "simulate_block_gaussian",
    list(n_samples = 80, n_total_genes = 1500, n_blocks = 6,
         genes_per_block = 50, rho = 0.9, signal_strength = 0.8,
         censoring = 0.45, signal_signs = c(1, -1)),
    103L
  ),
  make_setting(
    "two_source_n70_p2000_stable_dominant",
    "simulate_two_source_factor",
    list(n_samples = 70, n_total_genes = 2000,
         stable_genes = 80, supervised_genes = 8,
         stable_rho = 0.95, supervised_rho = 0.25,
         signal_strength = 0.8, stable_weight = 0.90,
         supervised_weight = 0.35, n_decoy_blocks = 10,
         decoy_genes = 40, decoy_rho = 0.90,
         censoring = 0.50),
    104L
  ),
  make_setting(
    "two_source_n80_p2000_stable_dominant",
    "simulate_two_source_factor",
    list(n_samples = 80, n_total_genes = 2000,
         stable_genes = 80, supervised_genes = 10,
         stable_rho = 0.95, supervised_rho = 0.30,
         signal_strength = 0.8, stable_weight = 0.85,
         supervised_weight = 0.45, n_decoy_blocks = 10,
         decoy_genes = 40, decoy_rho = 0.90,
         censoring = 0.45),
    105L
  ),
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
    106L
  )
)

summaries <- data.frame()
for (setting in settings) {
  cat("\n--- Target-gap extended coarse screening:", setting$name, "---\n")
  flush.console()
  ans <- run_core_setting(setting, R = 5L, n_cores = 5L, seed = 11200)
  summaries <- rbind(summaries, ans)
  print(ans, digits = 4)
  flush.console()
}

cat("\n=== Target-gap extended coarse screening summary ===\n")
summaries <- summaries[order(
  -(2 * summaries$diff_vs_DRPLS + summaries$diff_vs_partialCox)
), ]
print(summaries, digits = 4)

best <- settings[[match(summaries$setting[1L], vapply(settings, `[[`, "", "name"))]]
cat("\n--- Confirming best extended setting with R = 10:", best$name, "---\n")
confirm <- run_core_setting(best, R = 10L, n_cores = 10L, seed = 11900)
print(confirm, digits = 4)

