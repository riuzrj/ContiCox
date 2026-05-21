# User-facing wrappers for the paper methods.

auc_curve <- function(time, status, marker, auc_time_grid,
                      auc_method = c("IPCW", "NNE", "KM"),
                      ipcw_weighting = "marginal",
                      nne_span = 0.30) {
  .coxpcapls_auc_curve(
    time = time,
    status = status,
    marker = marker,
    auc_time_grid = auc_time_grid,
    auc_method = match.arg(auc_method),
    ipcw_weighting = ipcw_weighting,
    nne_span = nne_span
  )
}

integrated_auc <- function(auc_values) {
  .coxpcapls_iAUC(auc_values)
}

conticox_fit <- function(X, time, status,
                         n_components = min(7, ncol(as.matrix(X))),
                         alpha = 0.5,
                         residual_type = "deviance",
                         normalize_alpha_terms = TRUE,
                         ...) {
  cox_pca_pls_dr(
    X = X,
    time = time,
    status = status,
    n_components = n_components,
    alpha = alpha,
    residual_type = residual_type,
    normalize_alpha_terms = normalize_alpha_terms,
    return_all = TRUE,
    ...
  )
}

conticox_cv <- function(..., auc_method = c("IPCW", "NNE", "KM")) {
  val_DRPCAPLS_cv(..., auc_method = match.arg(auc_method))
}

conticox_cv_test <- function(..., auc_method = c("IPCW", "NNE", "KM")) {
  val_conticox_cv_test(..., auc_method = match.arg(auc_method))
}

plsdr_fit <- function(X, time, status,
                      n_components = min(7, ncol(as.matrix(X))),
                      ...) {
  cox_pls_dr(
    X = X,
    time = time,
    status = status,
    n_components = n_components,
    return_all = TRUE,
    ...
  )
}

plsdr_cv <- function(..., auc_method = c("IPCW", "NNE", "KM")) {
  val_DRPLS_cv(..., auc_method = match.arg(auc_method))
}

plsdr_cv_test <- function(..., auc_method = c("IPCW", "NNE", "KM")) {
  val_drpls_cv_test(..., auc_method = match.arg(auc_method))
}

partialcox_fit <- function(X, time, status,
                           n_components = min(7, ncol(as.matrix(X))),
                           ...) {
  partial_cox(
    X = X,
    time = time,
    status = status,
    n_components = n_components,
    ...
  )
}

partialcox_cv <- function(..., auc_method = c("IPCW", "NNE", "KM")) {
  val_partial_cox_cv(..., auc_method = match.arg(auc_method))
}

partialcox_cv_test <- function(..., auc_method = c("IPCW", "NNE", "KM")) {
  val_partial_cox_cv_test(..., auc_method = match.arg(auc_method))
}

penalizedcox_cv <- function(..., auc_method = c("IPCW", "NNE", "KM")) {
  val_penalized_cox_cv(..., auc_method = match.arg(auc_method))
}

penalizedcox_cv_test <- function(..., auc_method = c("IPCW", "NNE", "KM")) {
  val_penalized_cox_cv_test(..., auc_method = match.arg(auc_method))
}

ridgecox_cv <- function(...) {
  val_ridge_cox_cv(...)
}

ridgecox_cv_test <- function(...) {
  val_ridge_cox_cv_test(...)
}

enetcox_cv <- function(...) {
  val_elastic_net_cox_cv(...)
}

enetcox_cv_test <- function(...) {
  val_elastic_net_cox_cv_test(...)
}

.conticox_match_methods <- function(methods) {
  aliases <- c(
    conticox = "ContiCox",
    drpcapls = "ContiCox",
    plsdr = "DRPLS",
    drpls = "DRPLS",
    partialcox = "partialCox",
    partial_cox = "partialCox",
    pcox = "partialCox",
    ridgecox = "RidgeCox",
    ridge = "RidgeCox",
    elasticnetcox = "ElasticNetCox",
    enetcox = "ElasticNetCox",
    enet = "ElasticNetCox"
  )

  key <- gsub("[^[:alnum:]_]", "", tolower(methods))
  out <- unname(aliases[key])
  bad <- is.na(out)
  if (any(bad)) {
    stop("Unknown method(s): ", paste(methods[bad], collapse = ", "))
  }
  unique(out)
}

.conticox_extract_summary_row <- function(method, result) {
  if (!is.null(result$error)) {
    return(data.frame(
      method = method,
      best_cv_iAUC = NA_real_,
      test_iAUC = NA_real_,
      best_ncomp = NA_integer_,
      best_alpha = NA_real_,
      best_lambda = NA_real_,
      error = result$error,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    method = method,
    best_cv_iAUC = if (!is.null(result$best_cv_iAUC)) result$best_cv_iAUC else NA_real_,
    test_iAUC = if (!is.null(result$test_iAUC)) result$test_iAUC else NA_real_,
    best_ncomp = if (!is.null(result$best_ncomp)) result$best_ncomp else NA_integer_,
    best_alpha = if (!is.null(result$best_alpha)) result$best_alpha else NA_real_,
    best_lambda = if (!is.null(result$best_lambda)) result$best_lambda else NA_real_,
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

compare_survival_methods <- function(
    X, time, status,
    methods = c("ContiCox", "DRPLS", "partialCox", "RidgeCox", "ElasticNetCox"),
    test_prop = 0.3,
    train_idx = NULL,
    test_idx = NULL,
    k = 5,
    ncomp_candidates = 1:5,
    conticox_alpha_candidates = seq(0, 1, by = 0.1),
    ridge_lambda_candidates = 10 ^ seq(-4, 1, by = 1),
    enet_alpha_candidates = c(0.25, 0.5, 0.75),
    enet_lambda_candidates = 10 ^ seq(-4, 1, by = 1),
    auc_time_grid = NULL,
    use_preselection = FALSE,
    p_thresh = 0.05,
    max_features = NULL,
    auc_method = c("IPCW", "NNE", "KM"),
    seed = 123,
    residual_type = "deviance",
    normalize_alpha_terms = TRUE,
    standardize = TRUE,
    stop_on_error = FALSE,
    verbose = TRUE) {
  auc_method <- match.arg(auc_method)
  methods <- .conticox_match_methods(methods)

  dat <- .cvtest_check_data(X, time, status)
  X <- dat$X
  time <- dat$time
  status <- dat$status

  if (is.null(train_idx) || is.null(test_idx)) {
    split <- .cvtest_make_split(status, test_prop = test_prop, seed = seed)
    train_idx <- split$train_idx
    test_idx <- split$test_idx
  }

  if (is.null(auc_time_grid)) {
    auc_time_grid <- as.numeric(stats::quantile(
      time[train_idx],
      probs = seq(0.1, 0.8, length.out = 15),
      na.rm = TRUE
    ))
  }

  common_args <- list(
    X = X,
    time = time,
    status = status,
    train_idx = train_idx,
    test_idx = test_idx,
    k = k,
    auc_time_grid = auc_time_grid,
    use_preselection = use_preselection,
    p_thresh = p_thresh,
    auc_method = auc_method,
    seed = seed
  )

  run_one <- function(method) {
    if (isTRUE(verbose)) {
      message("Running ", method, " with ", auc_method, " iAUC...")
    }

    expr <- switch(
      method,
      ContiCox = do.call(val_conticox_cv_test, c(common_args, list(
        ncomp_candidates = ncomp_candidates,
        alpha_candidates = conticox_alpha_candidates,
        residual_type = residual_type,
        normalize_alpha_terms = normalize_alpha_terms
      ))),
      DRPLS = do.call(val_drpls_cv_test, c(common_args, list(
        ncomp_candidates = ncomp_candidates
      ))),
      partialCox = do.call(val_partial_cox_cv_test, c(common_args, list(
        ncomp_candidates = ncomp_candidates
      ))),
      RidgeCox = do.call(val_ridge_cox_cv_test, c(common_args, list(
        lambda_candidates = ridge_lambda_candidates,
        max_features = max_features,
        standardize = standardize
      ))),
      ElasticNetCox = do.call(val_elastic_net_cox_cv_test, c(common_args, list(
        alpha_candidates = enet_alpha_candidates,
        lambda_candidates = enet_lambda_candidates,
        max_features = max_features,
        standardize = standardize
      )))
    )

    expr
  }

  results <- stats::setNames(vector("list", length(methods)), methods)
  for (method in methods) {
    results[[method]] <- tryCatch(
      run_one(method),
      error = function(e) {
        if (isTRUE(stop_on_error)) {
          stop(e)
        }
        list(error = conditionMessage(e))
      }
    )
  }

  summary <- do.call(
    rbind,
    Map(.conticox_extract_summary_row, names(results), results)
  )
  rownames(summary) <- NULL

  list(
    summary = summary,
    results = results,
    train_idx = train_idx,
    test_idx = test_idx,
    auc_time_grid = auc_time_grid,
    auc_method = auc_method
  )
}
