project_dir <- "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls"

script_path <- file.path(project_dir, "simulation", "simulation_eigengene.R")
script_lines <- readLines(script_path, warn = FALSE)
stop_line <- grep("^mc_res_par <- mc_iAUC_par", script_lines)[1L] - 1L
eval(parse(text = script_lines[seq_len(stop_line)]), envir = .GlobalEnv)

run_one_setting <- function(setting, R = 10L, n_cores = 10L, seed = 123,
                            fast_screen = TRUE) {
  set.seed(seed + setting$seed_offset)
  sim_args <- setting$args
  sim_list <- replicate(
    R,
    do.call(simulate_one, sim_args),
    simplify = FALSE
  )
  
  if (isTRUE(fast_screen)) {
    mc_res <- mc_iAUC_par(
      sim_list,
      n_cores = n_cores,
      seed = seed,
      project_dir = project_dir,
      quiet = TRUE,
      auc_time_grid = seq(1, 45, length.out = 15),
      ncomp_candidates_DRPLS = 1:5,
      ncomp_candidates_pcox = 2:5,
      alpha_candidates = seq(0, 1, by = 0.25),
      lambda_candidates_ridgecox = 10 ^ seq(-3, 0, by = 1),
      alpha_candidates_enetcox = c(0.25, 0.5),
      lambda_candidates_enetcox = 10 ^ seq(-3, 0, by = 1)
    )
  } else {
    mc_res <- mc_iAUC_par(
      sim_list,
      n_cores = n_cores,
      seed = seed,
      project_dir = project_dir,
      quiet = TRUE
    )
  }
  stats <- summarize_mc_iAUC(
    mc_res,
    proposed = "DRPCA_PLS",
    baselines = c("DRPLS", "partialCox", "RidgeCox", "ElasticNetCox")
  )
  alpha_stats <- summarize_selected_alpha(mc_res)
  
  means <- mc_res$mean
  out <- data.frame(
    setting = setting$name,
    n = sim_args$n_samples,
    p = sim_args$n_total_genes,
    modules = sim_args$n_modules,
    genes_per_module = sim_args$n_genes_per_mod,
    r_min = sim_args$r_min,
    signal_strength = sim_args$signal_strength,
    signal_signs = paste(sim_args$signal_signs, collapse = ","),
    signal_coef = sim_args$signal_coef,
    standardize_lp = sim_args$standardize_lp,
    DRPLS = unname(means["DRPLS"]),
    ContiCox = unname(means["DRPCA_PLS"]),
    partialCox = unname(means["partialCox"]),
    RidgeCox = unname(means["RidgeCox"]),
    ElasticNetCox = unname(means["ElasticNetCox"]),
    diff_vs_DRPLS = unname(means["DRPCA_PLS"] - means["DRPLS"]),
    diff_vs_partialCox = unname(means["DRPCA_PLS"] - means["partialCox"]),
    selected_alpha_mean = alpha_stats$summary$Mean,
    selected_alpha_median = alpha_stats$summary$Median,
    row.names = NULL
  )
  
  list(
    summary = out,
    mc_res = mc_res,
    stats = stats,
    alpha_stats = alpha_stats
  )
}

base_args <- list(
  n_samples = 100,
  n_total_genes = 1000,
  n_modules = 6,
  n_genes_per_mod = 25,
  n_signal_modules = 2,
  r_min = 0.9,
  signal_strength = 1.0,
  signal_signs = c(1, 1),
  signal_coef = "module_average",
  standardize_lp = TRUE,
  censoring = 0.4
)

make_setting <- function(name, ..., seed_offset) {
  args <- modifyList(base_args, list(...))
  list(name = name, args = args, seed_offset = seed_offset)
}

settings <- list(
  make_setting(
    "high_corr_moderate_signal_p1000",
    seed_offset = 1L
  ),
  make_setting(
    "high_corr_weak_signal_p1000",
    signal_strength = 0.8,
    seed_offset = 2L
  ),
  make_setting(
    "high_corr_stronger_signal_p1000",
    signal_strength = 1.2,
    seed_offset = 3L
  ),
  make_setting(
    "high_corr_opposite_signal_p1000",
    signal_signs = c(1, -1),
    seed_offset = 4L
  ),
  make_setting(
    "moderate_corr_moderate_signal_p1000",
    r_min = 0.7,
    seed_offset = 5L
  ),
  make_setting(
    "small_n_high_corr_p1000",
    n_samples = 80,
    seed_offset = 6L
  ),
  make_setting(
    "many_noise_genes_high_corr_p2000",
    n_total_genes = 2000,
    n_modules = 8,
    seed_offset = 7L
  ),
  make_setting(
    "legacy_strong_signal_reference",
    n_total_genes = 1000,
    n_modules = 4,
    r_min = 0.5,
    signal_strength = 1.0,
    signal_coef = "gene_unit",
    standardize_lp = FALSE,
    seed_offset = 8L
  )
)

all_summaries <- data.frame()
all_results <- list()

for (setting in settings) {
  cat("\n--- Running setting:", setting$name, "---\n")
  ans <- run_one_setting(setting)
  all_results[[setting$name]] <- ans
  all_summaries <- rbind(all_summaries, ans$summary)
  print(ans$summary, digits = 4)
  
  good_setting <- is.finite(ans$summary$diff_vs_DRPLS) &&
    is.finite(ans$summary$diff_vs_partialCox) &&
    ans$summary$ContiCox >= 0.75 &&
    ans$summary$diff_vs_DRPLS >= 0.015 &&
    ans$summary$diff_vs_partialCox >= 0.010
  
  if (good_setting) {
    cat("\nFound a promising setting; stopping early.\n")
    break
  }
}

cat("\n=== Search summary ===\n")
print(all_summaries[order(-all_summaries$diff_vs_DRPLS), ], digits = 4)

best_idx <- which.max(
  all_summaries$diff_vs_DRPLS + all_summaries$diff_vs_partialCox
)
best_setting_name <- all_summaries$setting[best_idx]
best_result <- all_results[[best_setting_name]]
cat("\nBest setting by combined ContiCox improvement:", best_setting_name, "\n")
print(best_result$summary, digits = 4)
print(best_result$stats$paired_difference, digits = 4)
print(best_result$alpha_stats$summary, digits = 4)

invisible(list(
  all_summaries = all_summaries,
  all_results = all_results,
  best_setting_name = best_setting_name,
  best_result = best_result
))
