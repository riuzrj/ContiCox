project_dir <- "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls"
setwd(project_dir)

suppressPackageStartupMessages({
  source(file.path(project_dir, "simulation", "simulation_gaussian_latent_factor.R"))
  source(file.path(project_dir, "R", "pcapls_CR.R"))
  source(file.path(project_dir, "R", "DRPLS.R"))
  source(file.path(project_dir, "R", "DRPCAPLS.R"))
  source(file.path(project_dir, "R", "partial_cox.R"))
})

quiet_value <- function(expr) {
  value <- NULL
  invisible(utils::capture.output(
    value <- suppressWarnings(suppressMessages(eval.parent(substitute(expr)))),
    type = "output"
  ))
  value
}

get_iAUC <- function(res) {
  if (is.null(res) || is.null(res$best_iAUC)) return(NA_real_)
  val <- res$best_iAUC
  if (!is.finite(val)) return(NA_real_)
  val
}

eval_core_one <- function(dat, seed = 123, k = 5) {
  X <- dat$X
  time <- dat$surv$Time
  status <- dat$surv$Status
  set.seed(seed)
  folds <- sample(rep(seq_len(k), length.out = nrow(X)))
  
  common_grid <- list(
    k = k,
    auc_time_grid = seq(1, 45, length.out = 15),
    use_preselection = TRUE,
    folds = folds,
    seed = NULL
  )
  
  res_drpls <- tryCatch(
    quiet_value(do.call(
      val_DRPLS_cv,
      c(list(X = X, time = time, status = status,
             ncomp_candidates = 1:4), common_grid)
    )),
    error = function(e) NULL
  )
  
  res_conticox <- tryCatch(
    quiet_value(do.call(
      val_DRPCAPLS_cv,
      c(list(X = X, time = time, status = status,
             ncomp_candidates = 2:5,
             alpha_candidates = seq(0, 1, by = 0.25)), common_grid)
    )),
    error = function(e) NULL
  )
  
  res_pcox <- tryCatch(
    quiet_value(do.call(
      val_partial_cox_cv,
      c(list(X = X, time = time, status = status,
             ncomp_candidates = 2:5), common_grid)
    )),
    error = function(e) NULL
  )
  
  c(
    DRPLS = get_iAUC(res_drpls),
    ContiCox = get_iAUC(res_conticox),
    partialCox = get_iAUC(res_pcox),
    selected_alpha = if (!is.null(res_conticox) && is.finite(res_conticox$best_alpha)) {
      res_conticox$best_alpha
    } else {
      NA_real_
    }
  )
}

run_core_setting <- function(setting, R = 10L, n_cores = 10L, seed = 123) {
  set.seed(seed + setting$seed_offset)
  sim_list <- replicate(
    R,
    do.call(simulate_gaussian_latent_factor, setting$args),
    simplify = FALSE
  )
  
  res <- parallel::mclapply(
    seq_along(sim_list),
    function(i) eval_core_one(sim_list[[i]], seed = seed + i),
    mc.cores = min(n_cores, R)
  )
  mat <- do.call(rbind, res)
  means <- colMeans(mat, na.rm = TRUE)
  sds <- apply(mat, 2, stats::sd, na.rm = TRUE)
  
  data.frame(
    setting = setting$name,
    scenario = setting$args$scenario,
    n = setting$args$n_samples,
    p = setting$args$n_total_genes,
    n_factors = setting$args$n_factors,
    lowvar_index = setting$args$lowvar_index,
    eigen_decay = setting$args$eigen_decay,
    noise_eigenvalue = setting$args$noise_eigenvalue,
    signal_strength = setting$args$signal_strength,
    mix_a = setting$args$mix_a,
    mix_b = setting$args$mix_b,
    censoring = setting$args$censoring,
    DRPLS = means["DRPLS"],
    ContiCox = means["ContiCox"],
    partialCox = means["partialCox"],
    diff_vs_DRPLS = means["ContiCox"] - means["DRPLS"],
    diff_vs_partialCox = means["ContiCox"] - means["partialCox"],
    sd_ContiCox = sds["ContiCox"],
    selected_alpha_mean = means["selected_alpha"],
    selected_alpha_sd = sds["selected_alpha"],
    row.names = NULL
  )
}

base_args <- list(
  scenario = "III",
  n_samples = 100,
  n_total_genes = 1000,
  n_factors = 12,
  lowvar_index = 8,
  signal_strength = 0.9,
  mix_a = 0.6,
  mix_b = 0.8,
  leading_eigenvalue = 8,
  eigen_decay = 0.55,
  eigen_floor = 0.10,
  noise_eigenvalue = 0.10,
  base_rate = 0.1,
  censoring = 0.4
)

make_setting <- function(name, ..., seed_offset) {
  list(name = name, args = modifyList(base_args, list(...)), seed_offset = seed_offset)
}

settings <- list(
  make_setting("S3_mix_06_08_L8_gamma09_decay055", seed_offset = 1L),
  make_setting("S3_mix_04_09_L8_gamma09_decay055", mix_a = 0.4, mix_b = 0.9, seed_offset = 2L),
  make_setting("S3_mix_08_06_L8_gamma09_decay055", mix_a = 0.8, mix_b = 0.6, seed_offset = 3L),
  make_setting("S3_mix_06_08_L10_gamma09_decay055", lowvar_index = 10, seed_offset = 4L),
  make_setting("S3_mix_06_08_L10_gamma08_decay05", lowvar_index = 10, signal_strength = 0.8, eigen_decay = 0.50, seed_offset = 5L),
  make_setting("S3_mix_06_08_L10_gamma10_decay05", lowvar_index = 10, signal_strength = 1.0, eigen_decay = 0.50, seed_offset = 6L),
  make_setting("S2_L10_gamma10_decay05", scenario = "II", lowvar_index = 10, signal_strength = 1.0, eigen_decay = 0.50, seed_offset = 7L),
  make_setting("S2_L10_gamma12_decay05_censor03", scenario = "II", lowvar_index = 10, signal_strength = 1.2, eigen_decay = 0.50, censoring = 0.3, seed_offset = 8L),
  make_setting("S1_gamma08_decay055", scenario = "I", signal_strength = 0.8, eigen_decay = 0.55, seed_offset = 9L),
  make_setting("S3_p2000_mix_06_08_L10_gamma09", n_total_genes = 2000, lowvar_index = 10, signal_strength = 0.9, eigen_decay = 0.50, seed_offset = 10L),
  make_setting("S3_n80_p2000_mix_06_08_L10_gamma09", n_samples = 80, n_total_genes = 2000, lowvar_index = 10, signal_strength = 0.9, eigen_decay = 0.50, seed_offset = 11L)
)

summaries <- data.frame()
stop_early <- FALSE
for (setting in settings) {
  cat("\n--- Gaussian core screening:", setting$name, "---\n")
  flush.console()
  ans <- run_core_setting(setting)
  summaries <- rbind(summaries, ans)
  print(ans, digits = 4)
  flush.console()
  
  if (isTRUE(stop_early) &&
      is.finite(ans$diff_vs_DRPLS) &&
      is.finite(ans$diff_vs_partialCox) &&
      ans$ContiCox >= 0.72 &&
      ans$diff_vs_DRPLS >= 0.015 &&
      ans$diff_vs_partialCox >= 0.010) {
    cat("\nFound a promising Gaussian setting; stopping early.\n")
    break
  }
}

cat("\n=== Gaussian core screening summary ===\n")
summaries <- summaries[order(
  -(summaries$diff_vs_DRPLS + summaries$diff_vs_partialCox)
), ]
print(summaries, digits = 4)

best_setting <- summaries$setting[1L]
cat("\nBest Gaussian setting:", best_setting, "\n")
