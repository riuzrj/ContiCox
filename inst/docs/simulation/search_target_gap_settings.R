project_dir <- "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls"
setwd(project_dir)

suppressPackageStartupMessages({
  source(file.path(project_dir, "R", "pcapls_CR.R"))
  source(file.path(project_dir, "R", "DRPLS.R"))
  source(file.path(project_dir, "R", "DRPCAPLS.R"))
  source(file.path(project_dir, "R", "partial_cox.R"))
})

estimate_exponential_censor_rate <- function(event_time, target_censoring) {
  event_time <- event_time[is.finite(event_time) & event_time > 0]
  objective <- function(rate) mean(1 - exp(-rate * event_time)) - target_censoring
  upper <- 1 / stats::median(event_time)
  if (!is.finite(upper) || upper <= 0) upper <- 1
  while (objective(upper) < 0) upper <- upper * 2
  stats::uniroot(objective, lower = 0, upper = upper)$root
}

simulate_block_gaussian <- function(
    n_samples = 100,
    n_total_genes = 1000,
    n_blocks = 6,
    genes_per_block = 25,
    n_signal_blocks = 2,
    rho = 0.9,
    signal_strength = 1.0,
    signal_signs = c(1, -1),
    base_rate = 0.1,
    censoring = 0.4
) {
  stopifnot(n_blocks * genes_per_block <= n_total_genes)
  stopifnot(n_signal_blocks <= n_blocks)
  stopifnot(length(signal_signs) == n_signal_blocks)

  X <- matrix(rnorm(n_samples * n_total_genes), nrow = n_samples)
  colnames(X) <- paste0("Gene_", seq_len(n_total_genes))
  block_scores <- matrix(NA_real_, nrow = n_samples, ncol = n_blocks)

  current_col <- 1L
  for (block in seq_len(n_blocks)) {
    z <- rnorm(n_samples)
    block_scores[, block] <- z
    for (k in seq_len(genes_per_block)) {
      X[, current_col] <- sqrt(rho) * z + sqrt(1 - rho) * rnorm(n_samples)
      current_col <- current_col + 1L
    }
  }
  X <- as.matrix(scale(X, center = TRUE, scale = FALSE))

  raw_score <- rowSums(sweep(
    block_scores[, seq_len(n_signal_blocks), drop = FALSE],
    2,
    signal_signs,
    `*`
  ))
  raw_score <- as.numeric(scale(raw_score))
  lp <- signal_strength * raw_score

  event_time <- -log(runif(n_samples)) / (base_rate * exp(lp))
  censor_rate <- estimate_exponential_censor_rate(event_time, censoring)
  censor_time <- rexp(n_samples, rate = censor_rate)

  list(
    X = X,
    surv = data.frame(
      Time = pmin(event_time, censor_time),
      Status = as.integer(event_time <= censor_time)
    ),
    latent = block_scores,
    linear_predictor = lp,
    params = list(
      scenario = "block_gaussian_target_gap",
      n_samples = n_samples,
      n_total_genes = n_total_genes,
      n_blocks = n_blocks,
      genes_per_block = genes_per_block,
      n_signal_blocks = n_signal_blocks,
      rho = rho,
      signal_strength = signal_strength,
      signal_signs = signal_signs,
      censoring = censoring
    )
  )
}

simulate_two_source_factor <- function(
    n_samples = 100,
    n_total_genes = 1000,
    stable_genes = 40,
    supervised_genes = 20,
    stable_rho = 0.9,
    supervised_rho = 0.45,
    signal_strength = 1.0,
    stable_weight = 0.55,
    supervised_weight = 0.85,
    n_decoy_blocks = 6,
    decoy_genes = 30,
    decoy_rho = 0.85,
    base_rate = 0.1,
    censoring = 0.4
) {
  stopifnot(stable_genes + supervised_genes + n_decoy_blocks * decoy_genes <= n_total_genes)

  X <- matrix(rnorm(n_samples * n_total_genes), nrow = n_samples)
  colnames(X) <- paste0("Gene_", seq_len(n_total_genes))

  z_stable <- rnorm(n_samples)
  z_supervised <- rnorm(n_samples)
  current_col <- 1L

  for (j in seq_len(stable_genes)) {
    X[, current_col] <- sqrt(stable_rho) * z_stable +
      sqrt(1 - stable_rho) * rnorm(n_samples)
    current_col <- current_col + 1L
  }
  for (j in seq_len(supervised_genes)) {
    X[, current_col] <- sqrt(supervised_rho) * z_supervised +
      sqrt(1 - supervised_rho) * rnorm(n_samples)
    current_col <- current_col + 1L
  }
  for (block in seq_len(n_decoy_blocks)) {
    z_decoy <- rnorm(n_samples)
    for (j in seq_len(decoy_genes)) {
      X[, current_col] <- sqrt(decoy_rho) * z_decoy +
        sqrt(1 - decoy_rho) * rnorm(n_samples)
      current_col <- current_col + 1L
    }
  }

  X <- as.matrix(scale(X, center = TRUE, scale = FALSE))
  raw_score <- stable_weight * as.numeric(scale(z_stable)) +
    supervised_weight * as.numeric(scale(z_supervised))
  raw_score <- as.numeric(scale(raw_score))
  lp <- signal_strength * raw_score

  event_time <- -log(runif(n_samples)) / (base_rate * exp(lp))
  censor_rate <- estimate_exponential_censor_rate(event_time, censoring)
  censor_time <- rexp(n_samples, rate = censor_rate)

  list(
    X = X,
    surv = data.frame(
      Time = pmin(event_time, censor_time),
      Status = as.integer(event_time <= censor_time)
    ),
    latent = data.frame(stable = z_stable, supervised = z_supervised),
    linear_predictor = lp,
    params = list(
      scenario = "two_source_factor_target_gap",
      n_samples = n_samples,
      n_total_genes = n_total_genes,
      stable_genes = stable_genes,
      supervised_genes = supervised_genes,
      stable_rho = stable_rho,
      supervised_rho = supervised_rho,
      signal_strength = signal_strength,
      stable_weight = stable_weight,
      supervised_weight = supervised_weight,
      n_decoy_blocks = n_decoy_blocks,
      decoy_genes = decoy_genes,
      decoy_rho = decoy_rho,
      censoring = censoring
    )
  )
}

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

make_setting <- function(name, generator, args, seed_offset) {
  list(name = name, generator = generator, args = args, seed_offset = seed_offset)
}

run_core_setting <- function(setting, R = 5L, n_cores = 5L, seed = 9100) {
  set.seed(seed + setting$seed_offset)
  sim_list <- replicate(
    R,
    do.call(setting$generator, setting$args),
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
    generator = setting$generator,
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

settings <- list(
  make_setting(
    "block_n80_p1500_r09_g08_c45",
    "simulate_block_gaussian",
    list(n_samples = 80, n_total_genes = 1500, rho = 0.9,
         signal_strength = 0.8, censoring = 0.45,
         signal_signs = c(1, -1)),
    1L
  ),
  make_setting(
    "block_n80_p2000_r095_g08_c45",
    "simulate_block_gaussian",
    list(n_samples = 80, n_total_genes = 2000, rho = 0.95,
         signal_strength = 0.8, censoring = 0.45,
         signal_signs = c(1, -1)),
    2L
  ),
  make_setting(
    "block_n100_p2000_r095_g10_c40",
    "simulate_block_gaussian",
    list(n_samples = 100, n_total_genes = 2000, rho = 0.95,
         signal_strength = 1.0, censoring = 0.40,
         signal_signs = c(1, -1)),
    3L
  ),
  make_setting(
    "two_source_n80_p1500_decoy",
    "simulate_two_source_factor",
    list(n_samples = 80, n_total_genes = 1500,
         stable_genes = 45, supervised_genes = 12,
         stable_rho = 0.92, supervised_rho = 0.35,
         signal_strength = 0.9, stable_weight = 0.60,
         supervised_weight = 0.80, n_decoy_blocks = 8,
         decoy_genes = 35, decoy_rho = 0.88,
         censoring = 0.45),
    4L
  ),
  make_setting(
    "two_source_n100_p2000_decoy",
    "simulate_two_source_factor",
    list(n_samples = 100, n_total_genes = 2000,
         stable_genes = 50, supervised_genes = 15,
         stable_rho = 0.95, supervised_rho = 0.35,
         signal_strength = 1.0, stable_weight = 0.55,
         supervised_weight = 0.85, n_decoy_blocks = 10,
         decoy_genes = 35, decoy_rho = 0.90,
         censoring = 0.40),
    5L
  ),
  make_setting(
    "two_source_n80_p2000_sparse_supervised",
    "simulate_two_source_factor",
    list(n_samples = 80, n_total_genes = 2000,
         stable_genes = 60, supervised_genes = 8,
         stable_rho = 0.95, supervised_rho = 0.25,
         signal_strength = 0.9, stable_weight = 0.55,
         supervised_weight = 0.90, n_decoy_blocks = 10,
         decoy_genes = 40, decoy_rho = 0.90,
         censoring = 0.45),
    6L
  )
)

summaries <- data.frame()
for (setting in settings) {
  cat("\n--- Target-gap coarse screening:", setting$name, "---\n")
  flush.console()
  ans <- run_core_setting(setting, R = 5L, n_cores = 5L)
  summaries <- rbind(summaries, ans)
  print(ans, digits = 4)
  flush.console()
}

cat("\n=== Target-gap coarse screening summary ===\n")
summaries <- summaries[order(
  -(summaries$diff_vs_DRPLS + summaries$diff_vs_partialCox)
), ]
print(summaries, digits = 4)

best <- settings[[match(summaries$setting[1L], vapply(settings, `[[`, "", "name"))]]
cat("\n--- Confirming best setting with R = 10:", best$name, "---\n")
confirm <- run_core_setting(best, R = 10L, n_cores = 10L, seed = 9900)
print(confirm, digits = 4)

