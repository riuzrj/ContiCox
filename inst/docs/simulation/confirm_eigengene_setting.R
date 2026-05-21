project_dir <- "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls"
script_path <- file.path(project_dir, "simulation", "simulation_eigengene.R")
script_lines <- readLines(script_path, warn = FALSE)
stop_line <- grep("^mc_res_par <- mc_iAUC_par", script_lines)[1L] - 1L
eval(parse(text = script_lines[seq_len(stop_line)]), envir = .GlobalEnv)

set.seed(123)
sim_list_confirm <- replicate(
  10,
  simulate_one(
    n_samples = 100,
    n_total_genes = 1000,
    n_modules = 6,
    n_genes_per_mod = 25,
    n_signal_modules = 2,
    r_min = 0.9,
    signal_strength = 1.0,
    signal_signs = c(1, -1),
    signal_coef = "module_average",
    standardize_lp = TRUE,
    censoring = 0.4
  ),
  simplify = FALSE
)

mc_res_confirm <- mc_iAUC_par(
  sim_list_confirm,
  n_cores = 10,
  seed = 123,
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

mc_res_confirm$draws
mc_res_confirm$mean
mc_res_confirm_stats <- summarize_mc_iAUC(mc_res_confirm)
mc_res_confirm_stats$summary
mc_res_confirm_stats$paired_difference
summarize_selected_alpha(mc_res_confirm)
