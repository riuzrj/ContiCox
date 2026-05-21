# -----------------------------------------------------------------------------
# Train/test evaluation wrappers
# -----------------------------------------------------------------------------
# These helpers separate hyperparameter tuning from final performance reporting:
#   1. split data into training and independent test sets;
#   2. tune hyperparameters by cross-validation within the training set;
#   3. refit the selected model on the full training set;
#   4. evaluate time-dependent AUC/iAUC on the held-out test set only.
# -----------------------------------------------------------------------------

if (!exists(".coxpcapls_auc_curve", mode = "function")) {
  for (.auc_utils_path in c("R/auc_utils.R", "auc_utils.R")) {
    if (file.exists(.auc_utils_path)) {
      source(.auc_utils_path)
      break
    }
  }
  rm(.auc_utils_path)
}
if (!exists(".coxpcapls_auc_curve", mode = "function")) {
  stop("Cannot locate R/auc_utils.R; source it before train_test_validation.R.")
}

.cvtest_check_data <- function(X, time, status) {
  X <- as.matrix(X)
  time <- as.numeric(time)
  status <- as.integer(status != 0L)

  if (nrow(X) != length(time) || length(time) != length(status)) {
    stop("nrow(X), length(time), and length(status) must be identical.")
  }

  ok <- is.finite(time) & is.finite(status) & time > 0 &
    rowSums(is.finite(X)) == ncol(X)

  list(
    X = X[ok, , drop = FALSE],
    time = time[ok],
    status = status[ok],
    kept = which(ok)
  )
}

.cvtest_make_split <- function(status, test_prop = 0.3, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  idx_event <- which(status == 1L)
  idx_cens <- which(status == 0L)

  sample_group <- function(idx) {
    if (length(idx) == 0L) {
      return(integer(0))
    }
    n_test <- max(1L, floor(length(idx) * test_prop))
    n_test <- min(n_test, length(idx))
    sample(idx, n_test)
  }

  test_idx <- sort(c(sample_group(idx_event), sample_group(idx_cens)))
  train_idx <- setdiff(seq_along(status), test_idx)

  if (length(train_idx) == 0L || length(test_idx) == 0L) {
    stop("Train/test split failed; adjust test_prop.")
  }
  if (sum(status[train_idx] == 1L) < 3L || sum(status[test_idx] == 1L) < 1L) {
    stop("Not enough observed events in training or test set.")
  }

  list(train_idx = train_idx, test_idx = test_idx)
}

.cvtest_select_features <- function(X, time, status,
                                    use_preselection = FALSE,
                                    p_thresh = 0.05) {
  safe_var <- function(x) {
    v <- stats::var(x, na.rm = TRUE)
    if (is.na(v) || !is.finite(v)) 0 else v
  }

  raw_var <- apply(X, 2, safe_var)
  valid_idx <- which(is.finite(raw_var) & raw_var > 0)
  if (length(valid_idx) == 0L) {
    stop("No non-constant predictors remain.")
  }

  if (!use_preselection) {
    return(valid_idx)
  }

  pvals <- rep(Inf, length(valid_idx))
  for (ii in seq_along(valid_idx)) {
    j <- valid_idx[ii]
    fit <- tryCatch(
      suppressWarnings(survival::coxph(survival::Surv(time, status) ~ X[, j])),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      coef_tab <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
      if (!is.null(coef_tab)) {
        pval <- as.numeric(coef_tab[1, "Pr(>|z|)"])
        if (is.finite(pval)) {
          pvals[ii] <- pval
        }
      }
    }
  }

  keep_local <- which(is.finite(pvals) & pvals < p_thresh)
  if (length(keep_local) == 0L) {
    stop("No variables pass univariate Cox screening at p < ", p_thresh, ".")
  }

  valid_idx[keep_local]
}

.cvtest_auc_curve <- function(time, status, marker, auc_time_grid,
                              auc_method = c("IPCW", "NNE", "KM")) {
  .coxpcapls_auc_curve(
    time = time,
    status = status,
    marker = marker,
    auc_time_grid = auc_time_grid,
    auc_method = auc_method
  )
}

.cvtest_predict_conticox <- function(model, X_new) {
  T_new <- as.data.frame(as.matrix(X_new) %*% model$A)
  colnames(T_new) <- paste0("comp_", seq_len(ncol(T_new)))

  if (is.null(model$final_cox_model)) {
    stop("The fitted ContiCox Cox model is NULL.")
  }

  as.numeric(stats::predict(model$final_cox_model, newdata = T_new, type = "lp"))
}

val_conticox_cv_test <- function(
    X, time, status,
    test_prop = 0.3,
    train_idx = NULL,
    test_idx = NULL,
    k = 5,
    ncomp_candidates = 1:5,
    alpha_candidates = seq(0, 1, by = 0.1),
    auc_time_grid = NULL,
    use_preselection = FALSE,
    p_thresh = 0.05,
    auc_method = c("IPCW", "NNE", "KM"),
    seed = 123,
    residual_type = "deviance",
    normalize_alpha_terms = TRUE
) {
  auc_method <- match.arg(auc_method)

  dat <- .cvtest_check_data(X, time, status)
  X <- dat$X
  time <- dat$time
  status <- dat$status

  if (is.null(train_idx) || is.null(test_idx)) {
    split <- .cvtest_make_split(status, test_prop = test_prop, seed = seed)
    train_idx <- split$train_idx
    test_idx <- split$test_idx
  }

  X_train <- X[train_idx, , drop = FALSE]
  time_train <- time[train_idx]
  status_train <- status[train_idx]
  X_test <- X[test_idx, , drop = FALSE]
  time_test <- time[test_idx]
  status_test <- status[test_idx]

  if (is.null(auc_time_grid)) {
    auc_time_grid <- as.numeric(stats::quantile(
      time_train,
      probs = seq(0.1, 0.8, length.out = 15),
      na.rm = TRUE
    ))
  }

  feature_idx <- .cvtest_select_features(
    X = X_train,
    time = time_train,
    status = status_train,
    use_preselection = use_preselection,
    p_thresh = p_thresh
  )
  X_train_work <- X_train[, feature_idx, drop = FALSE]
  X_test_work <- X_test[, feature_idx, drop = FALSE]

  if (!is.null(seed)) {
    set.seed(seed + 1L)
  }
  inner_folds <- sample(rep(seq_len(k), length.out = nrow(X_train_work)))

  cv_res <- val_DRPCAPLS_cv(
    X = X_train_work,
    time = time_train,
    status = status_train,
    k = k,
    ncomp_candidates = ncomp_candidates,
    alpha_candidates = alpha_candidates,
    auc_time_grid = auc_time_grid,
    use_preselection = FALSE,
    auc_method = auc_method,
    folds = inner_folds,
    seed = NULL,
    residual_type = residual_type,
    normalize_alpha_terms = normalize_alpha_terms
  )

  final_model <- cox_pca_pls_dr(
    X = X_train_work,
    time = time_train,
    status = status_train,
    residual_type = residual_type,
    n_components = cv_res$best_ncomp,
    alpha = cv_res$best_alpha,
    normalize_alpha_terms = normalize_alpha_terms,
    return_all = TRUE
  )

  test_marker <- .cvtest_predict_conticox(final_model, X_test_work)
  test_auc_curve <- .cvtest_auc_curve(
    time = time_test,
    status = status_test,
    marker = test_marker,
    auc_time_grid = auc_time_grid,
    auc_method = auc_method
  )
  test_iAUC <- mean(test_auc_curve, na.rm = TRUE)

  list(
    best_model = final_model$final_cox_model,
    final_fit = final_model,
    cv_result = cv_res,
    best_ncomp = cv_res$best_ncomp,
    best_alpha = cv_res$best_alpha,
    best_cv_iAUC = cv_res$best_iAUC,
    best_iAUC = test_iAUC,
    test_iAUC = test_iAUC,
    best_auc_curve = test_auc_curve,
    test_auc_curve = test_auc_curve,
    auc_time_grid = auc_time_grid,
    test_marker = test_marker,
    metrics = list(iAUC = test_iAUC, CV_iAUC = cv_res$best_iAUC),
    feature_idx = feature_idx,
    train_idx = train_idx,
    test_idx = test_idx
  )
}

.cvtest_predict_drpls <- function(model, X_new) {
  if (is.null(model$final_cox_model)) {
    stop("The fitted DRPLS Cox model is NULL.")
  }

  X_new <- as.matrix(X_new)
  center <- model$X_center
  scale <- model$X_scale
  scale[!is.finite(scale) | scale == 0] <- 1
  center[!is.finite(center)] <- 0
  X_new <- sweep(sweep(X_new, 2, center, "-"), 2, scale, "/")

  T_new <- as.data.frame(X_new %*% model$A)
  colnames(T_new) <- paste0("comp_", seq_len(ncol(T_new)))
  as.numeric(stats::predict(model$final_cox_model, newdata = T_new, type = "lp"))
}

val_drpls_cv_test <- function(
    X, time, status,
    test_prop = 0.3,
    train_idx = NULL,
    test_idx = NULL,
    k = 5,
    ncomp_candidates = 1:5,
    auc_time_grid = NULL,
    use_preselection = FALSE,
    p_thresh = 0.05,
    auc_method = c("IPCW", "NNE", "KM"),
    seed = 123
) {
  auc_method <- match.arg(auc_method)

  dat <- .cvtest_check_data(X, time, status)
  X <- dat$X
  time <- dat$time
  status <- dat$status

  if (is.null(train_idx) || is.null(test_idx)) {
    split <- .cvtest_make_split(status, test_prop = test_prop, seed = seed)
    train_idx <- split$train_idx
    test_idx <- split$test_idx
  }

  X_train <- X[train_idx, , drop = FALSE]
  time_train <- time[train_idx]
  status_train <- status[train_idx]
  X_test <- X[test_idx, , drop = FALSE]
  time_test <- time[test_idx]
  status_test <- status[test_idx]

  if (is.null(auc_time_grid)) {
    auc_time_grid <- as.numeric(stats::quantile(time_train, probs = seq(0.1, 0.8, length.out = 15), na.rm = TRUE))
  }

  feature_idx <- .cvtest_select_features(X_train, time_train, status_train, use_preselection, p_thresh)
  X_train_work <- X_train[, feature_idx, drop = FALSE]
  X_test_work <- X_test[, feature_idx, drop = FALSE]

  if (!is.null(seed)) set.seed(seed + 1L)
  inner_folds <- sample(rep(seq_len(k), length.out = nrow(X_train_work)))

  cv_res <- val_DRPLS_cv(
    X = X_train_work,
    time = time_train,
    status = status_train,
    k = k,
    ncomp_candidates = ncomp_candidates,
    auc_time_grid = auc_time_grid,
    use_preselection = FALSE,
    auc_method = auc_method,
    folds = inner_folds,
    seed = NULL
  )

  final_model <- cox_pls_dr(
    X = X_train_work,
    time = time_train,
    status = status_train,
    n_components = cv_res$best_ncomp,
    return_all = TRUE
  )
  test_marker <- .cvtest_predict_drpls(final_model, X_test_work)
  test_auc_curve <- .cvtest_auc_curve(time_test, status_test, test_marker, auc_time_grid, auc_method)

  list(
    best_model = final_model$final_cox_model,
    final_fit = final_model,
    cv_result = cv_res,
    best_ncomp = cv_res$best_ncomp,
    best_cv_iAUC = cv_res$best_iAUC,
    best_iAUC = mean(test_auc_curve, na.rm = TRUE),
    test_iAUC = mean(test_auc_curve, na.rm = TRUE),
    best_auc_curve = test_auc_curve,
    test_auc_curve = test_auc_curve,
    auc_time_grid = auc_time_grid,
    test_marker = test_marker,
    metrics = list(iAUC = mean(test_auc_curve, na.rm = TRUE), CV_iAUC = cv_res$best_iAUC),
    feature_idx = feature_idx,
    train_idx = train_idx,
    test_idx = test_idx
  )
}

.cvtest_predict_partial_cox <- function(model, X_new) {
  T_new <- transform_partial_cox(
    X_new = X_new,
    betas_list = model$betas_list,
    w_list = model$w_list
  )
  T_new <- as.data.frame(T_new)
  colnames(T_new) <- paste0("T", seq_len(ncol(T_new)))
  as.numeric(stats::predict(model$final_model, newdata = T_new, type = "lp"))
}

val_partial_cox_cv_test <- function(
    X, time, status,
    test_prop = 0.3,
    train_idx = NULL,
    test_idx = NULL,
    k = 5,
    ncomp_candidates = 1:5,
    auc_time_grid = NULL,
    use_preselection = FALSE,
    p_thresh = 0.05,
    auc_method = c("IPCW", "NNE", "KM"),
    seed = 123
) {
  auc_method <- match.arg(auc_method)

  dat <- .cvtest_check_data(X, time, status)
  X <- dat$X
  time <- dat$time
  status <- dat$status

  if (is.null(train_idx) || is.null(test_idx)) {
    split <- .cvtest_make_split(status, test_prop = test_prop, seed = seed)
    train_idx <- split$train_idx
    test_idx <- split$test_idx
  }

  X_train <- X[train_idx, , drop = FALSE]
  time_train <- time[train_idx]
  status_train <- status[train_idx]
  X_test <- X[test_idx, , drop = FALSE]
  time_test <- time[test_idx]
  status_test <- status[test_idx]

  if (is.null(auc_time_grid)) {
    auc_time_grid <- as.numeric(stats::quantile(time_train, probs = seq(0.1, 0.8, length.out = 15), na.rm = TRUE))
  }

  feature_idx <- .cvtest_select_features(X_train, time_train, status_train, use_preselection, p_thresh)
  X_train_work <- X_train[, feature_idx, drop = FALSE]
  X_test_work <- X_test[, feature_idx, drop = FALSE]

  if (!is.null(seed)) set.seed(seed + 1L)
  inner_folds <- sample(rep(seq_len(k), length.out = nrow(X_train_work)))

  cv_res <- val_partial_cox_cv(
    X = X_train_work,
    time = time_train,
    status = status_train,
    k = k,
    ncomp_candidates = ncomp_candidates,
    auc_time_grid = auc_time_grid,
    use_preselection = FALSE,
    auc_method = auc_method,
    folds = inner_folds,
    seed = NULL
  )

  final_model <- partial_cox(
    X = X_train_work,
    time = time_train,
    status = status_train,
    n_components = cv_res$best_ncomp
  )
  test_marker <- .cvtest_predict_partial_cox(final_model, X_test_work)
  test_auc_curve <- .cvtest_auc_curve(time_test, status_test, test_marker, auc_time_grid, auc_method)

  list(
    best_model = final_model$final_model,
    final_fit = final_model,
    cv_result = cv_res,
    best_ncomp = cv_res$best_ncomp,
    best_cv_iAUC = cv_res$best_iAUC,
    best_iAUC = mean(test_auc_curve, na.rm = TRUE),
    test_iAUC = mean(test_auc_curve, na.rm = TRUE),
    best_auc_curve = test_auc_curve,
    test_auc_curve = test_auc_curve,
    auc_time_grid = auc_time_grid,
    test_marker = test_marker,
    metrics = list(iAUC = mean(test_auc_curve, na.rm = TRUE), CV_iAUC = cv_res$best_iAUC),
    feature_idx = feature_idx,
    train_idx = train_idx,
    test_idx = test_idx
  )
}

val_penalized_cox_cv_test <- function(
    X, time, status,
    test_prop = 0.3,
    train_idx = NULL,
    test_idx = NULL,
    k = 5,
    alpha_candidates = 0,
    lambda_candidates = 10 ^ seq(-4, 1, by = 1),
    auc_time_grid = NULL,
    use_preselection = FALSE,
    p_thresh = 0.05,
    max_features = NULL,
    auc_method = c("IPCW", "NNE", "KM"),
    seed = 123,
    standardize = TRUE
) {
  auc_method <- match.arg(auc_method)

  dat <- .cvtest_check_data(X, time, status)
  X <- dat$X
  time <- dat$time
  status <- dat$status

  if (is.null(train_idx) || is.null(test_idx)) {
    split <- .cvtest_make_split(status, test_prop = test_prop, seed = seed)
    train_idx <- split$train_idx
    test_idx <- split$test_idx
  }

  X_train <- X[train_idx, , drop = FALSE]
  time_train <- time[train_idx]
  status_train <- status[train_idx]
  X_test <- X[test_idx, , drop = FALSE]
  time_test <- time[test_idx]
  status_test <- status[test_idx]

  if (is.null(auc_time_grid)) {
    auc_time_grid <- as.numeric(stats::quantile(time_train, probs = seq(0.1, 0.8, length.out = 15), na.rm = TRUE))
  }

  if (!is.null(seed)) set.seed(seed + 1L)
  inner_folds <- sample(rep(seq_len(k), length.out = nrow(X_train)))

  cv_res <- val_penalized_cox_cv(
    X = X_train,
    time = time_train,
    status = status_train,
    k = k,
    alpha_candidates = alpha_candidates,
    lambda_candidates = lambda_candidates,
    auc_time_grid = auc_time_grid,
    use_preselection = use_preselection,
    p_thresh = p_thresh,
    max_features = max_features,
    auc_method = auc_method,
    folds = inner_folds,
    seed = NULL,
    standardize = standardize
  )

  final_fit <- fit_penalized_cox(
    X = X_train,
    time = time_train,
    status = status_train,
    alpha = cv_res$best_alpha,
    lambda = cv_res$best_lambda,
    use_preselection = use_preselection,
    p_thresh = p_thresh,
    max_features = max_features,
    standardize = standardize
  )
  test_marker <- predict_penalized_cox(final_fit, X_test, type = "lp")
  test_auc_curve <- .cvtest_auc_curve(time_test, status_test, test_marker, auc_time_grid, auc_method)

  list(
    best_model = final_fit$glmnet_fit,
    final_fit = final_fit,
    cv_result = cv_res,
    best_alpha = cv_res$best_alpha,
    best_lambda = cv_res$best_lambda,
    best_cv_iAUC = cv_res$best_iAUC,
    best_iAUC = mean(test_auc_curve, na.rm = TRUE),
    test_iAUC = mean(test_auc_curve, na.rm = TRUE),
    best_auc_curve = test_auc_curve,
    test_auc_curve = test_auc_curve,
    auc_time_grid = auc_time_grid,
    test_marker = test_marker,
    metrics = list(iAUC = mean(test_auc_curve, na.rm = TRUE), CV_iAUC = cv_res$best_iAUC),
    train_idx = train_idx,
    test_idx = test_idx
  )
}

val_ridge_cox_cv_test <- function(..., alpha_candidates = 0) {
  val_penalized_cox_cv_test(..., alpha_candidates = 0)
}

val_elastic_net_cox_cv_test <- function(...,
                                        alpha_candidates = c(0.25, 0.5, 0.75)) {
  val_penalized_cox_cv_test(..., alpha_candidates = alpha_candidates)
}
