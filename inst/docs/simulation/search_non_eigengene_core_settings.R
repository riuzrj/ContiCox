project_dir <- "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls"
setwd(project_dir)

script_lines <- readLines(file.path(project_dir, "simulation", "simulation_eigengene.R"), warn = FALSE)
simulate_one_start <- grep("^simulate_one <-", script_lines)[1L]
simulate_one_end <- grep("^estimate_exponential_censor_rate <-", script_lines)[1L] - 1L
cluster_start <- grep("^simulate_cluster_complex <-", script_lines)[1L]
cluster_end <- grep("^R <- 100", script_lines[cluster_start:length(script_lines)])[1L] +
  cluster_start - 2L
factorial_start <- grep("^simulate_factorial <-", script_lines)[1L]
factorial_end <- grep("^R <- 100", script_lines[factorial_start:length(script_lines)])[1L] +
  factorial_start - 2L
eval(parse(text = script_lines[simulate_one_start:simulate_one_end]), envir = .GlobalEnv)
eval(parse(text = script_lines[cluster_start:cluster_end]), envir = .GlobalEnv)
eval(parse(text = script_lines[factorial_start:factorial_end]), envir = .GlobalEnv)

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
    rho = 0.8,
    signal_strength = 1.0,
    signal_signs = NULL,
    decoy_strength = 0,
    base_rate = 0.1,
    censoring = 0.4
) {
  stopifnot(n_blocks * genes_per_block <= n_total_genes)
  stopifnot(n_signal_blocks <= n_blocks)
  if (is.null(signal_signs)) {
    signal_signs <- rep(1, n_signal_blocks)
  }
  stopifnot(length(signal_signs) == n_signal_blocks)
  
  X <- matrix(rnorm(n_samples * n_total_genes), nrow = n_samples)
  colnames(X) <- paste0("Gene_", seq_len(n_total_genes))
  block_scores <- matrix(NA_real_, nrow = n_samples, ncol = n_blocks)
  
  current_col <- 1L
  for (block in seq_len(n_blocks)) {
    z <- rnorm(n_samples)
    if (block > n_signal_blocks && decoy_strength > 0) {
      z <- z * decoy_strength
    }
    block_scores[, block] <- z
    for (k in seq_len(genes_per_block)) {
      X[, current_col] <- sqrt(rho) * z + sqrt(1 - rho) * rnorm(n_samples)
      current_col <- current_col + 1L
    }
  }
  X <- scale(X, center = TRUE, scale = FALSE)
  X <- as.matrix(X)
  
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
    params = list(
      scenario = "block_gaussian",
      n_samples = n_samples,
      n_total_genes = n_total_genes,
      n_blocks = n_blocks,
      genes_per_block = genes_per_block,
      n_signal_blocks = n_signal_blocks,
      rho = rho,
      signal_strength = signal_strength,
      signal_signs = signal_signs,
      decoy_strength = decoy_strength,
      censoring = censoring
    )
  )
}

simulate_factorial_signed <- function(
    n_samples = 100,
    n_total_genes = 1000,
    n_groups = 6,
    genes_per_group = 25,
    rho_within = 0.8,
    n_signal_groups = 2,
    signal_strength = 1.0,
    signal_signs = NULL,
    base_rate = 0.1,
    censoring = 0.4
) {
  stopifnot(n_groups * genes_per_group <= n_total_genes)
  if (is.null(signal_signs)) {
    signal_signs <- rep(1, n_signal_groups)
  }
  stopifnot(length(signal_signs) == n_signal_groups)
  
  dat <- simulate_factorial(
    n_samples = n_samples,
    n_total_genes = n_total_genes,
    n_groups = n_groups,
    genes_per_group = genes_per_group,
    rho_within = rho_within,
    n_signal_groups = n_signal_groups,
    beta_signal = 0.0,
    base_rate = base_rate,
    censor_rate = censoring
  )
  X <- dat$X
  group_score <- matrix(NA_real_, nrow = n_samples, ncol = n_signal_groups)
  for (g in seq_len(n_signal_groups)) {
    idx <- ((g - 1L) * genes_per_group + 1L):(g * genes_per_group)
    group_score[, g] <- rowMeans(scale(X[, idx, drop = FALSE]))
  }
  raw_score <- rowSums(sweep(group_score, 2, signal_signs, `*`))
  raw_score <- as.numeric(scale(raw_score))
  lp <- signal_strength * raw_score
  
  event_time <- -log(runif(n_samples)) / (base_rate * exp(lp))
  censor_rate <- estimate_exponential_censor_rate(event_time, censoring)
  censor_time <- rexp(n_samples, rate = censor_rate)
  dat$surv <- data.frame(
    Time = pmin(event_time, censor_time),
    Status = as.integer(event_time <= censor_time)
  )
  dat$linear_predictor <- lp
  dat$params$scenario <- "factorial_signed"
  dat$params$signal_strength <- signal_strength
  dat$params$signal_signs <- signal_signs
  dat$params$censoring <- censoring
  dat
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

run_core_setting <- function(setting, R = 10L, n_cores = 10L, seed = 123) {
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
    "block_r08_gamma08_opposite",
    "simulate_block_gaussian",
    list(rho = 0.8, signal_strength = 0.8, signal_signs = c(1, -1)),
    1L
  ),
  make_setting(
    "block_r08_gamma10_opposite",
    "simulate_block_gaussian",
    list(rho = 0.8, signal_strength = 1.0, signal_signs = c(1, -1)),
    2L
  ),
  make_setting(
    "block_r09_gamma08_opposite",
    "simulate_block_gaussian",
    list(rho = 0.9, signal_strength = 0.8, signal_signs = c(1, -1)),
    3L
  ),
  make_setting(
    "block_r09_gamma10_opposite",
    "simulate_block_gaussian",
    list(rho = 0.9, signal_strength = 1.0, signal_signs = c(1, -1)),
    4L
  ),
  make_setting(
    "block_r09_gamma10_same_decoy",
    "simulate_block_gaussian",
    list(rho = 0.9, signal_strength = 1.0, signal_signs = c(1, 1),
         n_blocks = 8, decoy_strength = 1.5),
    5L
  ),
  make_setting(
    "factorial_r08_gamma10_opposite",
    "simulate_factorial_signed",
    list(rho_within = 0.8, signal_strength = 1.0, signal_signs = c(1, -1)),
    6L
  ),
  make_setting(
    "factorial_r09_gamma10_opposite",
    "simulate_factorial_signed",
    list(rho_within = 0.9, signal_strength = 1.0, signal_signs = c(1, -1)),
    7L
  ),
  make_setting(
    "cluster_default",
    "simulate_cluster_complex",
    list(),
    8L
  ),
  make_setting(
    "cluster_highcorr",
    "simulate_cluster_complex",
    list(r_min = 0.9, beta_cont = 0.8, delta1 = 1.0, delta2 = 1.0),
    9L
  )
)

summaries <- data.frame()
for (setting in settings) {
  cat("\n--- Non-eigengene core screening:", setting$name, "---\n")
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
    cat("\nFound a promising non-eigengene setting; stopping early.\n")
    break
  }
}

cat("\n=== Non-eigengene core screening summary ===\n")
summaries <- summaries[order(
  -(summaries$diff_vs_DRPLS + summaries$diff_vs_partialCox)
), ]
print(summaries, digits = 4)

best_setting <- summaries$setting[1L]
cat("\nBest non-eigengene core setting:", best_setting, "\n")
