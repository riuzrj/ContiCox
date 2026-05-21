method_files <- c(
  file.path("R", "auc_utils.R"),
  file.path("R", "pcapls_CR.R"),
  file.path("R", "DRPLS.R"),
  file.path("R", "DRPCAPLS.R"),
  file.path("R", "partial_cox.R"),
  file.path("R", "penalized_cox_baselines.R"),
  file.path("R", "train_test_validation.R")
)

script_candidate_dirs <- function() {
  ofiles <- vapply(
    sys.frames(),
    function(frame) {
      val <- frame$ofile
      if (is.null(val)) "" else as.character(val)[1L]
    },
    character(1)
  )
  ofiles <- ofiles[nzchar(ofiles)]
  unique(dirname(normalizePath(ofiles, mustWork = FALSE)))
}

locate_project_dir <- function(project_dir = NULL, required_files = method_files) {
  starts <- unique(c(
    project_dir,
    getwd(),
    script_candidate_dirs(),
    "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls"
  ))
  starts <- starts[!is.na(starts) & nzchar(starts)]
  
  candidates <- unique(unlist(lapply(starts, function(path) {
    path <- normalizePath(path, mustWork = FALSE)
    c(path, dirname(path), dirname(dirname(path)))
  })))
  
  has_all_files <- vapply(
    candidates,
    function(path) all(file.exists(file.path(path, required_files))),
    logical(1)
  )
  
  if (!any(has_all_files)) {
    stop(
      "Cannot locate project root containing method files: ",
      paste(required_files, collapse = ", "),
      ". Please call run_gaussian_scenario(..., project_dir = ",
      "\"/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls\")."
    )
  }
  
  normalizePath(candidates[which(has_all_files)[1L]], mustWork = TRUE)
}

source_method_files <- function(project_dir = NULL, envir = parent.frame()) {
  project_dir <- locate_project_dir(project_dir)
  required_files <- file.path(project_dir, method_files)
  missing_files <- method_files[!file.exists(required_files)]
  if (length(missing_files) > 0L) {
    stop("Cannot find method files: ", paste(missing_files, collapse = ", "))
  }
  for (method_file in required_files) {
    source(method_file, local = envir)
  }
  invisible(TRUE)
}

ensure_method_functions <- function(project_dir = NULL, envir = parent.frame()) {
  required_functions <- c(
    "val_drpls_cv_test",
    "val_conticox_cv_test",
    "val_partial_cox_cv_test",
    "val_ridge_cox_cv_test",
    "val_elastic_net_cox_cv_test"
  )
  available <- vapply(
    required_functions,
    exists,
    logical(1),
    envir = envir,
    mode = "function",
    inherits = TRUE
  )
  if (!all(available)) {
    source_method_files(project_dir = project_dir, envir = envir)
  }
  invisible(TRUE)
}

.quiet_eval <- function(expr, quiet = TRUE) {
  if (!isTRUE(quiet)) {
    return(eval.parent(substitute(expr)))
  }
  value <- NULL
  invisible(utils::capture.output(
    value <- suppressWarnings(suppressMessages(eval.parent(substitute(expr)))),
    type = "output"
  ))
  value
}

estimate_exponential_censor_rate <- function(event_time, target_censoring) {
  stopifnot(target_censoring > 0, target_censoring < 1)
  event_time <- event_time[is.finite(event_time) & event_time > 0]
  if (length(event_time) < 1L) {
    stop("event_time must contain positive finite values.")
  }
  
  objective <- function(rate) {
    mean(1 - exp(-rate * event_time)) - target_censoring
  }
  
  upper <- 1 / stats::median(event_time)
  if (!is.finite(upper) || upper <= 0) {
    upper <- 1
  }
  while (objective(upper) < 0) {
    upper <- upper * 2
  }
  
  stats::uniroot(objective, lower = 0, upper = upper)$root
}

make_orthonormal_directions <- function(n_total_genes, n_factors) {
  raw_basis <- matrix(
    stats::rnorm(n_total_genes * n_factors),
    nrow = n_total_genes,
    ncol = n_factors
  )
  qr.Q(qr(raw_basis))
}

simulate_gaussian_latent_factor <- function(
    n_samples = 100,
    n_total_genes = 1000,
    scenario = c("I", "II", "III"),
    n_factors = 12,
    lowvar_index = 6,
    signal_strength = 1.0,
    mix_a = 0.6,
    mix_b = 0.8,
    leading_eigenvalue = 8,
    eigen_decay = 0.60,
    eigen_floor = 0.10,
    noise_eigenvalue = 0.10,
    base_rate = 0.1,
    censoring = 0.4
) {
  scenario <- match.arg(scenario)
  if (!is.numeric(n_samples) || length(n_samples) != 1L ||
      !is.finite(n_samples) || n_samples <= 2) {
    stop("n_samples must be a single numeric value greater than 2.")
  }
  if (!is.numeric(n_total_genes) || length(n_total_genes) != 1L ||
      !is.finite(n_total_genes) || n_total_genes < 2) {
    stop("n_total_genes must be a single numeric value at least 2.")
  }
  stopifnot(n_factors >= lowvar_index)
  stopifnot(n_total_genes >= n_factors)
  stopifnot(lowvar_index >= 2)
  stopifnot(signal_strength > 0)
  stopifnot(eigen_decay > 0, eigen_decay < 1)
  stopifnot(noise_eigenvalue > 0)
  stopifnot(base_rate > 0)
  stopifnot(censoring > 0, censoring < 1)
  n_samples <- as.integer(n_samples)
  n_total_genes <- as.integer(n_total_genes)
  
  Q <- make_orthonormal_directions(n_total_genes, n_factors)
  lambda_factors <- leading_eigenvalue * eigen_decay ^ seq(0, n_factors - 1L) +
    eigen_floor
  
  factor_scores <- matrix(
    stats::rnorm(n_samples * n_factors),
    nrow = n_samples,
    ncol = n_factors
  )
  factor_scores <- sweep(factor_scores, 2, sqrt(lambda_factors), `*`)
  residual_noise <- matrix(
    stats::rnorm(n_samples * n_total_genes, sd = sqrt(noise_eigenvalue)),
    nrow = n_samples,
    ncol = n_total_genes
  )
  X <- factor_scores %*% t(Q) + residual_noise
  X <- scale(X, center = TRUE, scale = FALSE)
  X <- as.matrix(X)
  colnames(X) <- paste0("Gene_", seq_len(n_total_genes))
  
  beta_direction <- switch(
    scenario,
    I = Q[, 1L],
    II = Q[, lowvar_index],
    III = mix_a * Q[, 1L] + mix_b * Q[, lowvar_index]
  )
  beta_direction <- as.numeric(beta_direction / sqrt(sum(beta_direction ^ 2)))
  names(beta_direction) <- colnames(X)
  
  raw_score <- as.numeric(X %*% beta_direction)
  raw_sd <- stats::sd(raw_score)
  if (is.finite(raw_sd) && raw_sd > 0) {
    raw_score <- as.numeric(scale(raw_score))
    beta_true <- (signal_strength / raw_sd) * beta_direction
  } else {
    beta_true <- signal_strength * beta_direction
  }
  names(beta_true) <- colnames(X)
  eta <- signal_strength * raw_score
  
  event_time <- -log(stats::runif(n_samples)) / (base_rate * exp(eta))
  censor_rate <- estimate_exponential_censor_rate(event_time, censoring)
  censor_time <- stats::rexp(n_samples, rate = censor_rate)
  
  time <- pmin(event_time, censor_time)
  status <- as.integer(event_time <= censor_time)
  
  list(
    X = X,
    surv = data.frame(Time = time, Status = status),
    beta_direction = beta_direction,
    beta_true = beta_true,
    eigenvectors = Q,
    eigenvalues = lambda_factors + noise_eigenvalue,
    noise_eigenvalue = noise_eigenvalue,
    linear_predictor = eta,
    params = list(
      scenario = scenario,
      n_samples = n_samples,
      n_total_genes = n_total_genes,
      n_factors = n_factors,
      lowvar_index = lowvar_index,
      signal_strength = signal_strength,
      mix_a = mix_a,
      mix_b = mix_b,
      leading_eigenvalue = leading_eigenvalue,
      eigen_decay = eigen_decay,
      eigen_floor = eigen_floor,
      noise_eigenvalue = noise_eigenvalue,
      base_rate = base_rate,
      censoring = censoring,
      censor_rate = censor_rate
    )
  )
}

make_gaussian_scenario_list <- function(
    scenario = c("I", "II", "III"),
    R = 100,
    seed = 123,
    n_samples = 100,
    n_total_genes = 1000,
    n_factors = 12,
    lowvar_index = 6,
    signal_strength = 1.0,
    mix_a = 0.6,
    mix_b = 0.8,
    leading_eigenvalue = 8,
    eigen_decay = 0.60,
    eigen_floor = 0.10,
    noise_eigenvalue = 0.10,
    base_rate = 0.1,
    censoring = 0.4
) {
  scenario <- match.arg(scenario)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  replicate(
    R,
    simulate_gaussian_latent_factor(
      scenario = scenario,
      n_samples = n_samples,
      n_total_genes = n_total_genes,
      n_factors = n_factors,
      lowvar_index = lowvar_index,
      signal_strength = signal_strength,
      mix_a = mix_a,
      mix_b = mix_b,
      leading_eigenvalue = leading_eigenvalue,
      eigen_decay = eigen_decay,
      eigen_floor = eigen_floor,
      noise_eigenvalue = noise_eigenvalue,
      base_rate = base_rate,
      censoring = censoring
    ),
    simplify = FALSE
  )
}

gaussian_latent_factor_grid <- expand.grid(
  scenario = c("I", "II", "III"),
  n_total_genes = c(1000L, 2000L),
  signal_strength = c(0.5, 1.0, 1.5),
  lowvar_index = c(6L),
  KEEP.OUT.ATTRS = FALSE
)

eval_iAUC_one_dataset <- function(
    X, time, status,
    k = 5,
    test_prop = 0.3,
    auc_time_grid = seq(1, 45, length.out = 30),
    ncomp_candidates_DRPLS = 1:10,
    ncomp_candidates_pcox  = 2:10,
    alpha_candidates       = seq(0, 1, by = 0.1),
    lambda_candidates_ridgecox = 10 ^ seq(-4, 1, by = 1),
    alpha_candidates_enetcox = c(0.25, 0.5, 0.75),
    lambda_candidates_enetcox = 10 ^ seq(-4, 1, by = 1),
    use_preselection       = TRUE,
    auc_method = "IPCW",
    train_idx = NULL,
    test_idx = NULL,
    seed = NULL
) {
  ensure_method_functions(envir = parent.env(environment()))
  auc_method <- match.arg(auc_method, choices = c("IPCW", "NNE", "KM"))
  stopifnot(nrow(X) == length(time), length(time) == length(status))

  if (is.null(train_idx) || is.null(test_idx)) {
    split <- .cvtest_make_split(status, test_prop = test_prop, seed = seed)
    train_idx <- split$train_idx
    test_idx <- split$test_idx
  }
  
  get_iAUC <- function(res) {
    if (is.null(res) || is.null(res$test_iAUC)) return(NA_real_)
    val <- res$test_iAUC
    if (!is.finite(val)) return(NA_real_)
    val
  }
  get_cv_iAUC <- function(res) {
    if (is.null(res) || is.null(res$best_cv_iAUC)) return(NA_real_)
    val <- res$best_cv_iAUC
    if (!is.finite(val)) return(NA_real_)
    val
  }
  get_selected_alpha <- function(res) {
    if (is.null(res) || is.null(res$best_alpha)) return(NA_real_)
    val <- res$best_alpha
    if (!is.finite(val)) return(NA_real_)
    val
  }
  get_selected_ncomp <- function(res) {
    if (is.null(res) || is.null(res$best_ncomp)) return(NA_real_)
    val <- res$best_ncomp
    if (!is.finite(val)) return(NA_real_)
    val
  }
  
  res_drpls <- tryCatch(
    val_drpls_cv_test(
      X, time, status,
      test_prop = test_prop,
      train_idx = train_idx,
      test_idx = test_idx,
      k = k,
      ncomp_candidates = ncomp_candidates_DRPLS,
      auc_time_grid    = auc_time_grid,
      use_preselection = use_preselection,
      auc_method       = auc_method,
      seed             = seed
    ),
    error = function(e) NULL
  )

  res_cppc_alpha0 <- tryCatch(
    val_conticox_cv_test(
      X, time, status,
      test_prop = test_prop,
      train_idx = train_idx,
      test_idx = test_idx,
      k = k,
      ncomp_candidates = ncomp_candidates_DRPLS,
      alpha_candidates = 0,
      auc_time_grid    = auc_time_grid,
      use_preselection = use_preselection,
      auc_method       = auc_method,
      seed             = seed
    ),
    error = function(e) NULL
  )

  res_cppc_alpha1 <- tryCatch(
    val_conticox_cv_test(
      X, time, status,
      test_prop = test_prop,
      train_idx = train_idx,
      test_idx = test_idx,
      k = k,
      ncomp_candidates = ncomp_candidates_DRPLS,
      alpha_candidates = 1,
      auc_time_grid    = auc_time_grid,
      use_preselection = use_preselection,
      auc_method       = auc_method,
      seed             = seed
    ),
    error = function(e) NULL
  )
  
  res_drpcapls <- tryCatch(
    val_conticox_cv_test(
      X, time, status,
      test_prop = test_prop,
      train_idx = train_idx,
      test_idx = test_idx,
      k = k,
      ncomp_candidates = ncomp_candidates_DRPLS,
      alpha_candidates = sort(unique(c(0, 1, alpha_candidates))),
      auc_time_grid    = auc_time_grid,
      use_preselection = use_preselection,
      auc_method       = auc_method,
      seed             = seed
    ),
    error = function(e) NULL
  )
  
  res_pcox <- tryCatch(
    val_partial_cox_cv_test(
      X, time, status,
      test_prop = test_prop,
      train_idx = train_idx,
      test_idx = test_idx,
      k = k,
      ncomp_candidates = ncomp_candidates_pcox,
      auc_time_grid    = auc_time_grid,
      use_preselection = use_preselection,
      auc_method       = auc_method,
      seed             = seed
    ),
    error = function(e) NULL
  )
  
  res_ridgecox <- tryCatch(
    val_ridge_cox_cv_test(
      X, time, status,
      test_prop = test_prop,
      train_idx = train_idx,
      test_idx = test_idx,
      k = k,
      lambda_candidates = lambda_candidates_ridgecox,
      auc_time_grid    = auc_time_grid,
      use_preselection = use_preselection,
      auc_method       = auc_method,
      seed             = seed
    ),
    error = function(e) {
      warning("RidgeCox failed: ", e$message)
      NULL
    }
  )
  
  res_enetcox <- tryCatch(
    val_elastic_net_cox_cv_test(
      X, time, status,
      test_prop = test_prop,
      train_idx = train_idx,
      test_idx = test_idx,
      k = k,
      alpha_candidates = alpha_candidates_enetcox,
      lambda_candidates = lambda_candidates_enetcox,
      auc_time_grid    = auc_time_grid,
      use_preselection = use_preselection,
      auc_method       = auc_method,
      seed             = seed
    ),
    error = function(e) {
      warning("ElasticNetCox failed: ", e$message)
      NULL
    }
  )
  
  list(
    iAUC = c(
      DRPLS          = get_iAUC(res_drpls),
      CPPC_alpha0    = get_iAUC(res_cppc_alpha0),
      CPPC_alpha1    = get_iAUC(res_cppc_alpha1),
      DRPCA_PLS      = get_iAUC(res_drpcapls),
      partialCox     = get_iAUC(res_pcox),
      RidgeCox       = get_iAUC(res_ridgecox),
      ElasticNetCox  = get_iAUC(res_enetcox)
    ),
    cv_iAUC = c(
      DRPLS          = get_cv_iAUC(res_drpls),
      CPPC_alpha0    = get_cv_iAUC(res_cppc_alpha0),
      CPPC_alpha1    = get_cv_iAUC(res_cppc_alpha1),
      DRPCA_PLS      = get_cv_iAUC(res_drpcapls),
      partialCox     = get_cv_iAUC(res_pcox),
      RidgeCox       = get_cv_iAUC(res_ridgecox),
      ElasticNetCox  = get_cv_iAUC(res_enetcox)
    ),
    selected_alpha = get_selected_alpha(res_drpcapls),
    selected_ncomp = get_selected_ncomp(res_drpcapls),
    train_idx = train_idx,
    test_idx = test_idx,
    auc_method = auc_method
  )
}

mc_iAUC_seq <- function(
    sim_list,
    k = 5,
    test_prop = 0.3,
    auc_time_grid = seq(1, 45, length.out = 30),
    ncomp_candidates_DRPLS = 1:10,
    ncomp_candidates_pcox  = 2:10,
    alpha_candidates       = seq(0, 1, by = 0.1),
    lambda_candidates_ridgecox = 10 ^ seq(-4, 1, by = 1),
    alpha_candidates_enetcox = c(0.25, 0.5, 0.75),
    lambda_candidates_enetcox = 10 ^ seq(-4, 1, by = 1),
    use_preselection       = TRUE,
    auc_method = "IPCW",
    seed = 123,
    quiet = TRUE
) {
  auc_method <- match.arg(auc_method, choices = c("IPCW", "NNE", "KM"))
  R <- length(sim_list)
  method_names <- c(
    "DRPLS", "CPPC_alpha0", "CPPC_alpha1", "DRPCA_PLS",
    "partialCox", "RidgeCox", "ElasticNetCox"
  )
  auc_mat <- matrix(NA_real_, nrow = R, ncol = length(method_names))
  cv_auc_mat <- matrix(NA_real_, nrow = R, ncol = length(method_names))
  selected_alpha <- rep(NA_real_, R)
  selected_ncomp <- rep(NA_real_, R)
  colnames(auc_mat) <- method_names
  colnames(cv_auc_mat) <- method_names
  
  for (r in seq_len(R)) {
    dat <- sim_list[[r]]
    X <- dat$X
    time <- dat$surv$Time
    status <- dat$surv$Status
    split_seed <- if (is.null(seed)) NULL else seed + r - 1L
    
    eval_res <- .quiet_eval(eval_iAUC_one_dataset(
      X, time, status,
      k = k,
      test_prop = test_prop,
      auc_time_grid = auc_time_grid,
      ncomp_candidates_DRPLS = ncomp_candidates_DRPLS,
      ncomp_candidates_pcox  = ncomp_candidates_pcox,
      alpha_candidates       = alpha_candidates,
      lambda_candidates_ridgecox = lambda_candidates_ridgecox,
      alpha_candidates_enetcox = alpha_candidates_enetcox,
      lambda_candidates_enetcox = lambda_candidates_enetcox,
      use_preselection       = use_preselection,
      auc_method             = auc_method,
      seed                   = split_seed
    ), quiet = quiet)
    auc_mat[r, ] <- eval_res$iAUC[method_names]
    cv_auc_mat[r, ] <- eval_res$cv_iAUC[method_names]
    selected_alpha[r] <- eval_res$selected_alpha
    selected_ncomp[r] <- eval_res$selected_ncomp
  }
  
  mean_iAUC <- colMeans(auc_mat, na.rm = TRUE)
  sd_iAUC   <- apply(auc_mat, 2, sd, na.rm = TRUE)
  n_eff     <- colSums(is.finite(auc_mat))
  mcse_iAUC <- sd_iAUC / sqrt(n_eff)
  
  ci_mean <- cbind(
    CI_lower = mean_iAUC - 1.96 * mcse_iAUC,
    CI_upper = mean_iAUC + 1.96 * mcse_iAUC
  )
  
  list(draws = auc_mat, test_iAUC = auc_mat, cv_iAUC = cv_auc_mat,
       auc_method = auc_method, selected_alpha = selected_alpha,
       selected_ncomp = selected_ncomp, mean = mean_iAUC, sd = sd_iAUC,
       mcse = mcse_iAUC, ci_mean = ci_mean, n_eff = n_eff, R = R)
}

.format_mean_sd <- function(mean, sd, digits = 3) {
  ifelse(
    is.finite(mean) & is.finite(sd),
    paste0(formatC(mean, format = "f", digits = digits),
           " \u00b1 ",
           formatC(sd, format = "f", digits = digits)),
    NA_character_
  )
}

.format_p_value <- function(p, digits = 3) {
  ifelse(
    is.finite(p),
    ifelse(p < 10 ^ (-digits),
           paste0("<", formatC(10 ^ (-digits), format = "f", digits = digits)),
           formatC(p, format = "f", digits = digits)),
    NA_character_
  )
}

summarize_mc_iAUC <- function(mc_res,
                              proposed = "DRPCA_PLS",
                              baselines = NULL,
                              digits = 3) {
  draws <- as.data.frame(mc_res$draws)
  methods <- colnames(draws)
  if (is.null(baselines)) {
    baselines <- setdiff(methods, proposed)
  }
  
  mean_vals <- colMeans(mc_res$draws, na.rm = TRUE)
  sd_vals <- apply(mc_res$draws, 2, stats::sd, na.rm = TRUE)
  n_eff <- colSums(is.finite(mc_res$draws))
  mcse_vals <- sd_vals / sqrt(n_eff)
  ci_lower <- mean_vals - 1.96 * mcse_vals
  ci_upper <- mean_vals + 1.96 * mcse_vals
  
  summary_table <- data.frame(
    Method = methods,
    Mean = as.numeric(mean_vals[methods]),
    SD = as.numeric(sd_vals[methods]),
    Mean_SD = .format_mean_sd(mean_vals[methods], sd_vals[methods], digits),
    MCSE = as.numeric(mcse_vals[methods]),
    CI_lower = as.numeric(ci_lower[methods]),
    CI_upper = as.numeric(ci_upper[methods]),
    N_eff = as.integer(n_eff[methods]),
    row.names = NULL
  )
  
  paired_table <- do.call(rbind, lapply(baselines, function(base) {
    if (!(proposed %in% methods) || !(base %in% methods)) {
      return(NULL)
    }
    keep <- is.finite(draws[[proposed]]) & is.finite(draws[[base]])
    diff_vals <- draws[[proposed]][keep] - draws[[base]][keep]
    n_pair <- length(diff_vals)
    diff_mean <- if (n_pair > 0L) mean(diff_vals) else NA_real_
    diff_sd <- if (n_pair > 1L) stats::sd(diff_vals) else NA_real_
    diff_mcse <- if (n_pair > 1L) diff_sd / sqrt(n_pair) else NA_real_
    p_val <- if (n_pair > 1L && any(diff_vals != 0)) {
      tryCatch(
        stats::wilcox.test(draws[[proposed]][keep],
                           draws[[base]][keep],
                           paired = TRUE,
                           exact = FALSE)$p.value,
        error = function(e) NA_real_
      )
    } else {
      NA_real_
    }
    data.frame(
      Proposed = proposed,
      Baseline = base,
      N_pair = n_pair,
      Mean_diff = diff_mean,
      SD_diff = diff_sd,
      MCSE_diff = diff_mcse,
      CI_lower = diff_mean - 1.96 * diff_mcse,
      CI_upper = diff_mean + 1.96 * diff_mcse,
      Wilcoxon_p = p_val,
      Wilcoxon_p_fmt = .format_p_value(p_val, digits),
      row.names = NULL
    )
  }))
  
  if (is.null(paired_table)) {
    paired_table <- data.frame(
      Proposed = character(),
      Baseline = character(),
      N_pair = integer(),
      Mean_diff = numeric(),
      SD_diff = numeric(),
      MCSE_diff = numeric(),
      CI_lower = numeric(),
      CI_upper = numeric(),
      Wilcoxon_p = numeric(),
      Wilcoxon_p_fmt = character()
    )
  }
  
  list(
    summary = summary_table,
    paired_difference = paired_table
  )
}

plot_paired_difference <- function(mc_res,
                                   proposed = "DRPCA_PLS",
                                   baselines = c(
                                     "DRPLS", "partialCox",
                                     "RidgeCox", "ElasticNetCox"
                                   ),
                                   method_labels = c(
                                     DRPLS = "DRPLS",
                                     CPPC_alpha0 = "CPPC-Cox alpha=0",
                                     CPPC_alpha1 = "CPPC-Cox alpha=1",
                                     DRPCA_PLS = "ContiCox",
                                     partialCox = "partialCox",
                                     RidgeCox = "Ridge Cox",
                                     ElasticNetCox = "Elastic-net Cox"
                                   ),
                                   diff_limits = NULL,
                                   diff_breaks = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plot_paired_difference().")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package `tidyr` is required for plot_paired_difference().")
  }
  
  draws <- as.data.frame(mc_res$draws)
  if (!(proposed %in% colnames(draws))) {
    stop("proposed method not found in mc_res$draws: ", proposed)
  }
  baselines <- baselines[baselines %in% colnames(draws)]
  if (length(baselines) < 1L) {
    stop("No requested baselines were found in mc_res$draws.")
  }
  
  paired_table <- summarize_mc_iAUC(
    mc_res,
    proposed = proposed,
    baselines = baselines
  )$paired_difference
  paired_table$Baseline_label <- ifelse(
    paired_table$Baseline %in% names(method_labels),
    unname(method_labels[paired_table$Baseline]),
    paired_table$Baseline
  )
  baseline_order <- ifelse(
    baselines %in% names(method_labels),
    unname(method_labels[baselines]),
    baselines
  )
  paired_table$Baseline_label <- factor(
    paired_table$Baseline_label,
    levels = rev(baseline_order)
  )
  
  diff_wide <- as.data.frame(lapply(baselines, function(base) {
    draws[[proposed]] - draws[[base]]
  }))
  colnames(diff_wide) <- baselines
  diff_wide$rep <- seq_len(nrow(diff_wide))
  diff_long <- tidyr::pivot_longer(
    diff_wide,
    cols = -rep,
    names_to = "Baseline",
    values_to = "Difference"
  )
  diff_long <- diff_long[is.finite(diff_long$Difference), , drop = FALSE]
  diff_long$Baseline_label <- ifelse(
    diff_long$Baseline %in% names(method_labels),
    unname(method_labels[diff_long$Baseline]),
    diff_long$Baseline
  )
  diff_long$Baseline_label <- factor(
    diff_long$Baseline_label,
    levels = baseline_order
  )
  
  x_scale <- list()
  y_scale <- list()
  if (!is.null(diff_breaks)) {
    x_scale <- c(x_scale, list(ggplot2::scale_x_continuous(breaks = diff_breaks)))
    y_scale <- c(y_scale, list(ggplot2::scale_y_continuous(breaks = diff_breaks)))
  }
  if (!is.null(diff_limits)) {
    x_scale <- c(x_scale, list(ggplot2::coord_cartesian(xlim = diff_limits)))
    y_scale <- c(y_scale, list(ggplot2::coord_cartesian(ylim = diff_limits)))
  }
  
  forest <- ggplot2::ggplot(paired_table, ggplot2::aes(y = Baseline_label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_segment(
      ggplot2::aes(x = CI_lower, xend = CI_upper,
                   yend = Baseline_label),
      linewidth = 0.55
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = Mean_diff),
      shape = 21, size = 2.4, fill = "red", color = "black"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    x_scale +
    ggplot2::labs(
      x = "Paired iAUC difference: ContiCox - baseline",
      y = "Baseline"
    )
  
  violin <- ggplot2::ggplot(
    diff_long,
    ggplot2::aes(x = Baseline_label, y = Difference, fill = Baseline_label)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.35, color = NA, width = 0.55) +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA,
                          fill = "white", color = "black", linewidth = 0.35) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 21,
                          size = 2.2, fill = "red", color = "black") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 15, hjust = 1)) +
    y_scale +
    ggplot2::labs(
      x = "Baseline",
      y = "Paired iAUC difference: ContiCox - baseline"
    )
  
  list(
    summary = paired_table,
    differences = diff_long,
    forest = forest,
    violin = violin
  )
}

summarize_selected_alpha <- function(mc_res, digits = 3) {
  if (is.null(mc_res$selected_alpha)) {
    stop("mc_res does not contain selected_alpha.")
  }
  
  alpha <- mc_res$selected_alpha
  alpha <- alpha[is.finite(alpha)]
  if (length(alpha) < 1L) {
    empty_summary <- data.frame(
      N = 0L,
      Mean = NA_real_,
      SD = NA_real_,
      Mean_SD = NA_character_,
      Median = NA_real_,
      Q25 = NA_real_,
      Q75 = NA_real_,
      Prop_low = NA_real_,
      Prop_mid = NA_real_,
      Prop_high = NA_real_
    )
    return(list(summary = empty_summary, frequency = data.frame()))
  }
  
  summary_table <- data.frame(
    N = length(alpha),
    Mean = mean(alpha),
    SD = stats::sd(alpha),
    Mean_SD = .format_mean_sd(mean(alpha), stats::sd(alpha), digits),
    Median = stats::median(alpha),
    Q25 = as.numeric(stats::quantile(alpha, 0.25, names = FALSE)),
    Q75 = as.numeric(stats::quantile(alpha, 0.75, names = FALSE)),
    Prop_low = mean(alpha <= 0.3),
    Prop_mid = mean(alpha > 0.3 & alpha < 0.7),
    Prop_high = mean(alpha >= 0.7),
    row.names = NULL
  )
  
  alpha_levels <- sort(unique(c(seq(0, 1, by = 0.1), alpha)))
  frequency <- as.data.frame(table(
    factor(alpha, levels = alpha_levels),
    useNA = "no"
  ))
  colnames(frequency) <- c("selected_alpha", "Freq")
  frequency$selected_alpha <- as.numeric(as.character(frequency$selected_alpha))
  frequency$Prop <- frequency$Freq / sum(frequency$Freq)
  frequency <- frequency[frequency$Freq > 0, , drop = FALSE]
  rownames(frequency) <- NULL
  
  list(summary = summary_table, frequency = frequency)
}

scientific_iAUC_palette <- function(method_levels) {
  base_palette <- c(
    "DRPLS" = "#0072B2",
    "CPPC-Cox alpha=0" = "#56B4E9",
    "CPPC-Cox alpha=1" = "#E69F00",
    "CPPC-Cox CV" = "#009E73",
    "ContiCox" = "#009E73",
    "partialCox" = "#D55E00",
    "Ridge Cox" = "#CC79A7",
    "Elastic-net Cox" = "#999999"
  )
  missing_levels <- setdiff(method_levels, names(base_palette))
  if (length(missing_levels) > 0L) {
    base_palette <- c(
      base_palette,
      stats::setNames(
        grDevices::hcl.colors(length(missing_levels), palette = "Dark 3"),
        missing_levels
      )
    )
  }
  base_palette[method_levels]
}

plot_mc_iAUC <- function(mc_res,
                         method_labels = c(
                           DRPLS = "DRPLS",
                           CPPC_alpha0 = "CPPC-Cox alpha=0",
                           CPPC_alpha1 = "CPPC-Cox alpha=1",
                           DRPCA_PLS = "CPPC-Cox CV",
                           partialCox = "partialCox",
                           RidgeCox = "Ridge Cox",
                           ElasticNetCox = "Elastic-net Cox"
                         ),
                         method_order = names(method_labels),
                         y_limits = c(0.5, 1.0),
                         y_breaks = seq(0.5, 1.0, by = 0.1),
                         method_palette = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plot_mc_iAUC().")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package `tidyr` is required for plot_mc_iAUC().")
  }
  
  draws <- as.data.frame(mc_res$draws)
  draws$rep <- seq_len(nrow(draws))
  df_long <- tidyr::pivot_longer(
    draws,
    cols = -rep,
    names_to = "Method",
    values_to = "iAUC"
  )
  df_long <- df_long[is.finite(df_long$iAUC), , drop = FALSE]
  df_long$Method <- ifelse(
    df_long$Method %in% names(method_labels),
    unname(method_labels[df_long$Method]),
    df_long$Method
  )
  order_labels <- unname(method_labels[method_order[method_order %in% names(method_labels)]])
  df_long$Method <- factor(df_long$Method, levels = order_labels)
  y_scale <- list()
  if (!is.null(y_breaks)) {
    y_scale <- c(y_scale, list(ggplot2::scale_y_continuous(breaks = y_breaks)))
  }
  if (!is.null(y_limits)) {
    y_scale <- c(y_scale, list(ggplot2::coord_cartesian(ylim = y_limits)))
  }
  if (is.null(method_palette)) {
    method_palette <- scientific_iAUC_palette(levels(df_long$Method))
  } else {
    method_palette <- method_palette[levels(df_long$Method)]
  }
  fill_scale <- ggplot2::scale_fill_manual(values = method_palette, drop = FALSE)
  
  boxplot <- ggplot2::ggplot(df_long, ggplot2::aes(x = Method, y = iAUC, fill = Method)) +
    ggplot2::geom_boxplot(outlier.size = 0.8, width = 0.55) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 21,
                          size = 2.2, fill = "white", color = "black") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 15, hjust = 1)) +
    fill_scale +
    y_scale +
    ggplot2::labs(x = "Method", y = "iAUC")
  
  violin <- ggplot2::ggplot(df_long, ggplot2::aes(x = Method, y = iAUC, fill = Method)) +
    ggplot2::geom_violin(trim = FALSE, scale = "width", alpha = 0.35,
                         color = NA, width = 0.55) +
    ggplot2::geom_boxplot(width = 0.10, outlier.shape = NA,
                          fill = "white", color = "black", linewidth = 0.35) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 21,
                          size = 2.2, fill = "white", color = "black") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 15, hjust = 1)) +
    fill_scale +
    y_scale +
    ggplot2::labs(x = "Method", y = "iAUC")
  
  list(data = df_long, boxplot = boxplot, violin = violin)
}

plot_main_iAUC <- function(mc_res,
                           method_labels = c(
                             DRPLS = "PLSDR",
                             DRPCA_PLS = "ContiCox",
                             partialCox = "pCox",
                             RidgeCox = "Ridge-Cox",
                             ElasticNetCox = "EN-Cox"
                           ),
                           method_order = names(method_labels),
                           y_limits = c(0.5, 1.0),
                           y_breaks = seq(0.5, 1.0, by = 0.1),
                           method_palette = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plot_main_iAUC().")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package `tidyr` is required for plot_main_iAUC().")
  }
  
  draws <- as.data.frame(mc_res$draws)
  keep_methods <- method_order[method_order %in% colnames(draws)]
  draws <- draws[, keep_methods, drop = FALSE]
  draws$rep <- seq_len(nrow(draws))
  df_long <- tidyr::pivot_longer(
    draws,
    cols = -rep,
    names_to = "Method",
    values_to = "iAUC"
  )
  df_long <- df_long[is.finite(df_long$iAUC), , drop = FALSE]
  df_long$Method <- unname(method_labels[df_long$Method])
  order_labels <- unname(method_labels[keep_methods])
  df_long$Method <- factor(df_long$Method, levels = order_labels)
  y_scale <- list()
  if (!is.null(y_breaks)) {
    y_scale <- c(y_scale, list(ggplot2::scale_y_continuous(breaks = y_breaks)))
  }
  if (!is.null(y_limits)) {
    y_scale <- c(y_scale, list(ggplot2::coord_cartesian(ylim = y_limits)))
  }
  if (is.null(method_palette)) {
    method_palette <- c(
      "PLSDR" = "#4C78A8",
      "ContiCox" = "#E45756",
      "pCox" = "#72B7B2",
      "Ridge-Cox" = "#B279A2",
      "EN-Cox" = "#79706E"
    )[levels(df_long$Method)]
  } else {
    method_palette <- method_palette[levels(df_long$Method)]
  }
  fill_scale <- ggplot2::scale_fill_manual(values = method_palette, drop = FALSE)
  
  boxplot <- ggplot2::ggplot(df_long, ggplot2::aes(x = Method, y = iAUC, fill = Method)) +
    ggplot2::geom_boxplot(outlier.size = 0.8, width = 0.55) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 21,
                          size = 2.2, fill = "white", color = "black") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 15, hjust = 1)) +
    fill_scale +
    y_scale +
    ggplot2::labs(x = "Method", y = "iAUC")
  
  violin <- ggplot2::ggplot(df_long, ggplot2::aes(x = Method, y = iAUC, fill = Method)) +
    ggplot2::geom_violin(trim = FALSE, scale = "width", alpha = 0.35,
                         color = NA, width = 0.55) +
    ggplot2::geom_boxplot(width = 0.10, outlier.shape = NA,
                          fill = "white", color = "black", linewidth = 0.35) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 21,
                          size = 2.2, fill = "white", color = "black") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 15, hjust = 1)) +
    fill_scale +
    y_scale +
    ggplot2::labs(x = "Method", y = "iAUC")
  
  list(data = df_long, boxplot = boxplot, violin = violin)
}

plot_selected_alpha <- function(mc_res) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plot_selected_alpha().")
  }
  if (is.null(mc_res$selected_alpha)) {
    stop("mc_res does not contain selected_alpha.")
  }
  
  df_alpha <- data.frame(
    rep = seq_along(mc_res$selected_alpha),
    selected_alpha = mc_res$selected_alpha
  )
  df_alpha <- df_alpha[is.finite(df_alpha$selected_alpha), , drop = FALSE]
  
  ggplot2::ggplot(df_alpha, ggplot2::aes(x = selected_alpha)) +
    ggplot2::geom_bar(width = 0.08, fill = "#4C78A8", color = "white") +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, by = 0.1),
                                limits = c(-0.05, 1.05)) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(x = "Selected alpha", y = "Frequency")
}

mc_iAUC_par <- function(
    sim_list,
    k = 5,
    test_prop = 0.3,
    auc_time_grid = seq(1, 45, length.out = 30),
    ncomp_candidates_DRPLS = 1:10,
    ncomp_candidates_pcox  = 2:10,
    alpha_candidates       = seq(0, 1, by = 0.1),
    lambda_candidates_ridgecox = 10 ^ seq(-4, 1, by = 1),
    alpha_candidates_enetcox = c(0.25, 0.5, 0.75),
    lambda_candidates_enetcox = 10 ^ seq(-4, 1, by = 1),
    use_preselection       = TRUE,
    auc_method = "IPCW",
    n_cores = max(1L, parallel::detectCores() - 2L),
    seed = 123,
    project_dir = NULL,
    quiet = TRUE
) {
  auc_method <- match.arg(auc_method, choices = c("IPCW", "NNE", "KM"))
  if (!requireNamespace("foreach", quietly = TRUE)) {
    stop("Package `foreach` is required for mc_iAUC_par().")
  }
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("Package `doParallel` is required for mc_iAUC_par().")
  }
  
  `%dopar%` <- foreach::`%dopar%`
  R <- length(sim_list)
  if (R < 1L) {
    stop("sim_list must contain at least one simulated dataset.")
  }
  
  n_cores <- max(1L, min(as.integer(n_cores), R))
  required_files <- method_files
  project_dir <- locate_project_dir(project_dir, required_files = required_files)
  missing_files <- required_files[!file.exists(file.path(project_dir, required_files))]
  if (length(missing_files) > 0L) {
    stop(
      "project_dir does not contain required files: ",
      paste(missing_files, collapse = ", ")
    )
  }
  
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)
  
  parallel::clusterCall(cl, function(project_dir, required_files, quiet) {
    load_worker <- function() {
      setwd(project_dir)
      missing_files <- required_files[!file.exists(required_files)]
      if (length(missing_files) > 0L) {
        stop("Worker cannot find files: ", paste(missing_files, collapse = ", "))
      }
      for (method_file in required_files) {
        source(method_file)
      }
      NULL
    }
    if (isTRUE(quiet)) {
      invisible(utils::capture.output(
        suppressWarnings(suppressMessages(load_worker())),
        type = "output"
      ))
    } else {
      load_worker()
    }
    NULL
  }, project_dir, required_files, quiet)
  
  auc_mat <- foreach::foreach(
    r = seq_len(R),
    .combine = rbind,
    .export = c(
      "eval_iAUC_one_dataset",
      "ensure_method_functions",
      "source_method_files",
      "locate_project_dir",
      "script_candidate_dirs",
      ".quiet_eval",
      "method_files"
    ),
    .packages = c(
      "survival", "pls", "timeROC", "survivalROC", "prodlim",
      "survcomp", "survAUC", "RSpectra", "glmnet"
    )
  ) %dopar% {
    dat <- sim_list[[r]]
    X <- dat$X
    time <- dat$surv$Time
    status <- dat$surv$Status
    split_seed <- if (is.null(seed)) NULL else seed + r - 1L
    
    eval_res <- .quiet_eval(eval_iAUC_one_dataset(
      X, time, status,
      k = k,
      test_prop = test_prop,
      auc_time_grid = auc_time_grid,
      ncomp_candidates_DRPLS = ncomp_candidates_DRPLS,
      ncomp_candidates_pcox  = ncomp_candidates_pcox,
      alpha_candidates       = alpha_candidates,
      lambda_candidates_ridgecox = lambda_candidates_ridgecox,
      alpha_candidates_enetcox = alpha_candidates_enetcox,
      lambda_candidates_enetcox = lambda_candidates_enetcox,
      use_preselection       = use_preselection,
      auc_method             = auc_method,
      seed                   = split_seed
    ), quiet = quiet)
    c(eval_res$iAUC,
      stats::setNames(eval_res$cv_iAUC, paste0("cv_", names(eval_res$cv_iAUC))),
      selected_alpha = eval_res$selected_alpha,
      selected_ncomp = eval_res$selected_ncomp)
  }
  
  method_names <- c(
    "DRPLS", "CPPC_alpha0", "CPPC_alpha1", "DRPCA_PLS",
    "partialCox", "RidgeCox", "ElasticNetCox"
  )
  cv_method_names <- paste0("cv_", method_names)
  worker_names <- c(method_names, cv_method_names, "selected_alpha", "selected_ncomp")
  auc_mat <- if (is.null(dim(auc_mat))) {
    matrix(auc_mat, nrow = 1L, byrow = TRUE)
  } else {
    as.matrix(auc_mat)
  }
  if (is.null(colnames(auc_mat))) {
    colnames(auc_mat) <- worker_names
  }
  selected_alpha <- as.numeric(auc_mat[, "selected_alpha"])
  selected_ncomp <- as.numeric(auc_mat[, "selected_ncomp"])
  cv_auc_mat <- auc_mat[, cv_method_names, drop = FALSE]
  colnames(cv_auc_mat) <- method_names
  auc_mat <- auc_mat[, method_names, drop = FALSE]
  colnames(auc_mat) <- method_names
  
  mean_iAUC <- colMeans(auc_mat, na.rm = TRUE)
  sd_iAUC   <- apply(auc_mat, 2, sd, na.rm = TRUE)
  n_eff     <- colSums(is.finite(auc_mat))
  mcse_iAUC <- sd_iAUC / sqrt(n_eff)
  
  ci_mean <- cbind(
    CI_lower = mean_iAUC - 1.96 * mcse_iAUC,
    CI_upper = mean_iAUC + 1.96 * mcse_iAUC
  )
  
  list(draws = auc_mat, test_iAUC = auc_mat, cv_iAUC = cv_auc_mat,
       auc_method = auc_method, selected_alpha = selected_alpha,
       selected_ncomp = selected_ncomp, mean = mean_iAUC, sd = sd_iAUC,
       mcse = mcse_iAUC, ci_mean = ci_mean, n_eff = n_eff, R = R)
}

run_gaussian_scenario <- function(
    scenario = c("I", "II", "III"),
    R = 100,
    n_cores = 10,
    n_samples = 100,
    n_total_genes = 1000,
    n_factors = 12,
    lowvar_index = 6,
    signal_strength = 1.0,
    mix_a = 0.6,
    mix_b = 0.8,
    leading_eigenvalue = 8,
    eigen_decay = 0.60,
    eigen_floor = 0.10,
    noise_eigenvalue = 0.10,
    base_rate = 0.1,
    censoring = 0.4,
    k = 5,
    test_prop = 0.3,
    auc_time_grid = seq(1, 45, length.out = 30),
    ncomp_candidates_DRPLS = 1:10,
    ncomp_candidates_pcox  = 2:10,
    alpha_candidates       = seq(0, 1, by = 0.1),
    lambda_candidates_ridgecox = 10 ^ seq(-4, 1, by = 1),
    alpha_candidates_enetcox = c(0.25, 0.5, 0.75),
    lambda_candidates_enetcox = 10 ^ seq(-4, 1, by = 1),
    use_preselection       = TRUE,
    auc_method = "IPCW",
    run_seq = FALSE,
    quiet = TRUE,
    seed = 123,
    project_dir = NULL
) {
  scenario <- match.arg(scenario)
  auc_method <- match.arg(auc_method, choices = c("IPCW", "NNE", "KM"))
  sim_list <- make_gaussian_scenario_list(
    scenario = scenario,
    R = R,
    seed = seed,
    n_samples = n_samples,
    n_total_genes = n_total_genes,
    n_factors = n_factors,
    lowvar_index = lowvar_index,
    signal_strength = signal_strength,
    mix_a = mix_a,
    mix_b = mix_b,
    leading_eigenvalue = leading_eigenvalue,
    eigen_decay = eigen_decay,
    eigen_floor = eigen_floor,
    noise_eigenvalue = noise_eigenvalue,
    base_rate = base_rate,
    censoring = censoring
  )
  
  mc_res_seq <- NULL
  mc_res_seq_stats <- NULL
  if (isTRUE(run_seq)) {
    mc_res_seq <- mc_iAUC_seq(
      sim_list[1],
      k = k,
      test_prop = test_prop,
      auc_time_grid = auc_time_grid,
      ncomp_candidates_DRPLS = ncomp_candidates_DRPLS,
      ncomp_candidates_pcox = ncomp_candidates_pcox,
      alpha_candidates = alpha_candidates,
      lambda_candidates_ridgecox = lambda_candidates_ridgecox,
      alpha_candidates_enetcox = alpha_candidates_enetcox,
      lambda_candidates_enetcox = lambda_candidates_enetcox,
      use_preselection = use_preselection,
      auc_method = auc_method,
      seed = seed,
      quiet = quiet
    )
    mc_res_seq_stats <- summarize_mc_iAUC(mc_res_seq)
  }
  
  mc_res_par <- mc_iAUC_par(
    sim_list,
    k = k,
    test_prop = test_prop,
    auc_time_grid = auc_time_grid,
    ncomp_candidates_DRPLS = ncomp_candidates_DRPLS,
    ncomp_candidates_pcox = ncomp_candidates_pcox,
    alpha_candidates = alpha_candidates,
    lambda_candidates_ridgecox = lambda_candidates_ridgecox,
    alpha_candidates_enetcox = alpha_candidates_enetcox,
    lambda_candidates_enetcox = lambda_candidates_enetcox,
    use_preselection = use_preselection,
    auc_method = auc_method,
    n_cores = n_cores,
    seed = seed,
    project_dir = project_dir,
    quiet = quiet
  )
  mc_res_par_stats <- summarize_mc_iAUC(mc_res_par)
  mc_res_par_plots <- plot_mc_iAUC(mc_res_par)
  mc_res_main_plots <- plot_main_iAUC(
    mc_res_par,
    y_limits = c(0.5, 1.0),
    y_breaks = seq(0.5, 1.0, by = 0.1)
  )
  paired_difference_plots <- plot_paired_difference(mc_res_par)
  selected_alpha_stats <- summarize_selected_alpha(mc_res_par)
  selected_alpha_plot <- plot_selected_alpha(mc_res_par)
  
  list(
    sim_list = sim_list,
    mc_res_seq = mc_res_seq,
    mc_res_seq_stats = mc_res_seq_stats,
    mc_res_par = mc_res_par,
    mc_res_par_stats = mc_res_par_stats,
    mc_res_par_plots = mc_res_par_plots,
    mc_res_main_plots = mc_res_main_plots,
    paired_difference_plots = paired_difference_plots,
    selected_alpha_stats = selected_alpha_stats,
    selected_alpha_plot = selected_alpha_plot
  )
}

if (interactive() && identical(Sys.getenv("RUN_GAUSSIAN_SCENARIOS"), "1")) {
  set.seed(123)
  
  sim_list_scenario_I <- make_gaussian_scenario_list(
    scenario = "I",
    R = 100,
    seed = 123,
    n_samples = 100,
    n_total_genes = 1000,
    lowvar_index = 6,
    signal_strength = 1.0,
    censoring = 0.4
  )
  sim_list <- sim_list_scenario_I
  
  mc_res_seq <- mc_iAUC_seq(sim_list[1], seed = 123)
  mc_res_seq$draws
  mc_res_seq$mean
  mc_res_seq_stats <- summarize_mc_iAUC(mc_res_seq)
  mc_res_seq_stats$summary
  mc_res_seq_stats$paired_difference
  
  mc_res_par <- mc_iAUC_par(
    sim_list[1:100],
    n_cores = 10,
    seed = 123,
    project_dir = NULL
  )
  mc_res_par$draws
  mc_res_par$mean
  mc_res_par_stats <- summarize_mc_iAUC(mc_res_par)
  mc_res_par_stats$summary
  mc_res_par_stats$paired_difference
  selected_alpha_stats <- summarize_selected_alpha(mc_res_par)
  selected_alpha_stats$summary
  selected_alpha_stats$frequency
  mc_res_par_plots <- plot_mc_iAUC(mc_res_par)
  mc_res_par_plots$boxplot
  mc_res_par_plots$violin
  mc_res_main_plots <- plot_main_iAUC(
    mc_res_par,
    y_limits = c(0.5, 1.0),
    y_breaks = seq(0.5, 1.0, by = 0.1)
  )
  mc_res_main_plots$boxplot
  mc_res_main_plots$violin
  paired_difference_plots <- plot_paired_difference(mc_res_par)
  paired_difference_plots$forest
  paired_difference_plots$violin
  selected_alpha_plot <- plot_selected_alpha(mc_res_par)
  selected_alpha_plot
}
