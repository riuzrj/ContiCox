project_dir <- "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls"
setwd(project_dir)

script_lines <- readLines(file.path(project_dir, "simulation", "simulation_eigengene.R"), warn = FALSE)
simulate_one_start <- grep("^simulate_one <-", script_lines)[1L]
simulate_one_end <- grep("^estimate_exponential_censor_rate <-", script_lines)[1L] - 1L
eval(parse(text = script_lines[simulate_one_start:simulate_one_end]), envir = .GlobalEnv)

suppressPackageStartupMessages({
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
    do.call(simulate_one, setting$args),
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
    n = setting$args$n_samples,
    p = setting$args$n_total_genes,
    modules = setting$args$n_modules,
    r_min = setting$args$r_min,
    signal_strength = setting$args$signal_strength,
    signal_signs = paste(setting$args$signal_signs, collapse = ","),
    signal_coef = setting$args$signal_coef,
    standardize_lp = setting$args$standardize_lp,
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
  list(name = name, args = modifyList(base_args, list(...)), seed_offset = seed_offset)
}

settings <- list(
  make_setting("p1000_r09_gamma08_same", signal_strength = 0.8, seed_offset = 1L),
  make_setting("p1000_r09_gamma10_same", signal_strength = 1.0, seed_offset = 2L),
  make_setting("p1000_r09_gamma12_same", signal_strength = 1.2, seed_offset = 3L),
  make_setting("p1000_r07_gamma10_same", r_min = 0.7, signal_strength = 1.0, seed_offset = 4L),
  make_setting("p1000_r09_gamma10_opposite", signal_signs = c(1, -1), seed_offset = 5L),
  make_setting("p2000_r09_gamma08_same", n_total_genes = 2000, signal_strength = 0.8, seed_offset = 6L),
  make_setting("p2000_r09_gamma10_same", n_total_genes = 2000, signal_strength = 1.0, seed_offset = 7L),
  make_setting("n80_p2000_r09_gamma10_same", n_samples = 80, n_total_genes = 2000, seed_offset = 8L),
  make_setting("p2000_modules8_r09_gamma10_same", n_total_genes = 2000, n_modules = 8, seed_offset = 9L),
  make_setting("legacy_geneunit_unscaled", n_total_genes = 1000, n_modules = 4,
               r_min = 0.5, signal_strength = 1.0,
               signal_coef = "gene_unit", standardize_lp = FALSE,
               seed_offset = 10L)
)

summaries <- data.frame()
for (setting in settings) {
  cat("\n--- Core screening:", setting$name, "---\n")
  flush.console()
  ans <- run_core_setting(setting)
  summaries <- rbind(summaries, ans)
  print(ans, digits = 4)
  flush.console()
  
  if (is.finite(ans$diff_vs_DRPLS) &&
      is.finite(ans$diff_vs_partialCox) &&
      ans$ContiCox >= 0.72 &&
      ans$diff_vs_DRPLS >= 0.015 &&
      ans$diff_vs_partialCox >= 0.010) {
    cat("\nFound a promising core setting; stopping early.\n")
    break
  }
}

cat("\n=== Core screening summary ===\n")
summaries <- summaries[order(
  -(summaries$diff_vs_DRPLS + summaries$diff_vs_partialCox)
), ]
print(summaries, digits = 4)

best_setting <- summaries$setting[1L]
cat("\nBest core setting:", best_setting, "\n")
