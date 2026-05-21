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
  stop("Cannot locate R/auc_utils.R; source it before penalized_cox_baselines.R.")
}


# -----------------------------------------------------------------------------
# Penalized Cox baselines: ridge Cox and elastic-net Cox
# -----------------------------------------------------------------------------
# These are standard Cox partial likelihood baselines fitted with glmnet:
#   alpha = 0       ridge Cox
#   0 < alpha < 1   elastic-net Cox
#   alpha = 1       lasso Cox, available through the shared helper if needed
#
# The CV wrappers mirror the repository's val_*_cv return format so simulations
# can compare them directly with DRPLS, DRPCAPLS, partial Cox, and PPLS-Cox.
# -----------------------------------------------------------------------------


.pencox_clean_x <- function(X) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  X[!is.finite(X)] <- NA_real_

  center <- colMeans(X, na.rm = TRUE)
  center[!is.finite(center)] <- 0
  for (j in seq_len(ncol(X))) {
    miss <- !is.finite(X[, j])
    if (any(miss)) {
      X[miss, j] <- center[j]
    }
  }

  X
}


.pencox_safe_var <- function(x) {
  v <- stats::var(x, na.rm = TRUE)
  if (is.na(v) || !is.finite(v)) 0 else v
}


.pencox_positive_time <- function(time, eps = 1e-8) {
  time <- as.numeric(time)
  if (any(!is.finite(time))) {
    stop("time contains non-finite values.")
  }
  positive <- time > 0
  if (!any(positive)) {
    stop("At least one positive survival time is required.")
  }
  if (any(!positive)) {
    replacement <- min(min(time[positive]) / 2, eps)
    time[!positive] <- replacement
  }
  time
}


.pencox_univariate_pvalues <- function(X, time, status) {
  pvals <- rep(Inf, ncol(X))

  for (j in seq_len(ncol(X))) {
    if (.pencox_safe_var(X[, j]) <= 0) {
      next
    }
    fit <- tryCatch(
      suppressWarnings(survival::coxph(survival::Surv(time, status) ~ X[, j])),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      next
    }
    coef_mat <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
    if (!is.null(coef_mat)) {
      pval <- as.numeric(coef_mat[1, "Pr(>|z|)"])
      if (is.finite(pval)) {
        pvals[j] <- pval
      }
    }
  }

  pvals
}


.pencox_select_features <- function(X,
                                    time,
                                    status,
                                    use_preselection = FALSE,
                                    p_thresh = 0.05,
                                    max_features = NULL,
                                    min_features = 2L) {
  raw_var <- apply(X, 2, .pencox_safe_var)
  valid_idx <- which(is.finite(raw_var) & raw_var > 0)
  if (length(valid_idx) == 0L) {
    stop("No non-constant predictors remain after preprocessing.")
  }

  if (use_preselection) {
    pvals <- .pencox_univariate_pvalues(X[, valid_idx, drop = FALSE], time, status)
    keep_local <- which(is.finite(pvals) & pvals < p_thresh)
    if (length(keep_local) == 0L) {
      keep_local <- order(pvals, decreasing = FALSE, na.last = NA)
    }
    if (length(keep_local) == 0L) {
      keep_local <- order(raw_var[valid_idx], decreasing = TRUE)
    }
    keep_local <- keep_local[order(pvals[keep_local], decreasing = FALSE, na.last = NA)]
    keep_idx <- valid_idx[keep_local]
  } else {
    keep_idx <- valid_idx
  }

  if (!is.null(max_features) && length(keep_idx) > max_features) {
    ord <- order(raw_var[keep_idx], decreasing = TRUE)
    keep_n <- max(min_features, min(max_features, length(keep_idx)))
    keep_idx <- keep_idx[ord[seq_len(keep_n)]]
  }

  keep_idx
}


fit_penalized_cox <- function(X,
                              time,
                              status,
                              alpha = 0,
                              lambda = NULL,
                              use_preselection = FALSE,
                              p_thresh = 0.05,
                              max_features = NULL,
                              standardize = TRUE,
                              selected_idx = NULL,
                              ...) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package `glmnet` is required for penalized Cox baselines.")
  }

  X <- .pencox_clean_x(X)
  n <- nrow(X)
  if (length(time) != n || length(status) != n) {
    stop("length(time) and length(status) must match nrow(X).")
  }
  time <- .pencox_positive_time(time)
  status <- as.integer(status != 0L)
  if (sum(status == 1L, na.rm = TRUE) < 2L) {
    stop("At least two observed events are required.")
  }
  if (!is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop("alpha must lie in [0, 1].")
  }

  if (is.null(selected_idx)) {
    selected_idx <- .pencox_select_features(
      X = X,
      time = time,
      status = status,
      use_preselection = use_preselection,
      p_thresh = p_thresh,
      max_features = max_features
    )
  }
  X_work <- X[, selected_idx, drop = FALSE]

  if (!is.null(lambda)) {
    lambda <- sort(unique(as.numeric(lambda[is.finite(lambda) & lambda > 0])), decreasing = TRUE)
    if (length(lambda) == 0L) {
      stop("lambda must contain at least one positive finite value.")
    }
  }

  y <- survival::Surv(time, status)
  fit <- glmnet::glmnet(
    x = X_work,
    y = y,
    family = "cox",
    alpha = alpha,
    lambda = lambda,
    standardize = standardize,
    ...
  )

  lambda_use <- if (is.null(lambda)) {
    fit$lambda[length(fit$lambda)]
  } else {
    min(lambda)
  }

  beta_work <- as.numeric(stats::coef(fit, s = lambda_use))
  beta_full <- numeric(ncol(X))
  beta_full[selected_idx] <- beta_work
  names(beta_full) <- colnames(X)
  score <- as.numeric(stats::predict(fit, newx = X_work, s = lambda_use, type = "link"))

  list(
    glmnet_fit = fit,
    final_model = fit,
    beta = beta_full,
    beta_std = beta_full,
    score = score,
    lambda = lambda_use,
    alpha = alpha,
    selected_idx = selected_idx,
    feature_names = colnames(X),
    standardize = standardize,
    use_preselection = use_preselection,
    p_thresh = p_thresh,
    max_features = max_features,
    method = if (alpha == 0) "ridge_cox" else "elastic_net_cox"
  )
}


predict_penalized_cox <- function(object,
                                  X_new,
                                  type = c("lp", "risk", "score")) {
  type <- match.arg(type)
  X_new <- .pencox_clean_x(X_new)

  if (!is.null(object$feature_names) && !is.null(colnames(X_new))) {
    missing_names <- setdiff(object$feature_names, colnames(X_new))
    if (length(missing_names) > 0L) {
      stop("X_new is missing ", length(missing_names), " training features.")
    }
    X_new <- X_new[, object$feature_names, drop = FALSE]
  }

  X_work <- X_new[, object$selected_idx, drop = FALSE]
  lp <- as.numeric(stats::predict(
    object$glmnet_fit,
    newx = X_work,
    s = object$lambda,
    type = "link"
  ))

  if (type %in% c("lp", "score")) {
    return(lp)
  }
  exp(lp)
}


.pencox_metrics <- function(time, status, score) {
  harrell_c <- NA_real_
  uno_c <- NA_real_

  if (requireNamespace("survcomp", quietly = TRUE)) {
    harrell_c <- tryCatch(
      survcomp::concordance.index(
        x = score,
        surv.time = time,
        surv.event = status,
        method = "noether"
      )$c.index,
      error = function(e) NA_real_
    )
  }

  if (requireNamespace("survAUC", quietly = TRUE)) {
    uno_c <- tryCatch(
      survAUC::UnoC(
        survival::Surv(time, status),
        survival::Surv(time, status),
        score
      )$C,
      error = function(e) NA_real_
    )
  }

  list(Harrell_Cindex = harrell_c, Uno_Cindex = uno_c)
}


val_penalized_cox_cv <- function(X,
                                 time,
                                 status,
                                 k = 5,
                                 alpha_candidates = 0,
                                 lambda_candidates = 10 ^ seq(-4, 1, by = 1),
                                 auc_time_grid = NULL,
                                 use_preselection = FALSE,
                                 p_thresh = 0.05,
                                 max_features = NULL,
                                 auc_method = c("IPCW", "NNE", "KM"),
                                 folds = NULL,
                                 seed = NULL,
                                 standardize = TRUE) {
  auc_method <- match.arg(auc_method)

  X <- .pencox_clean_x(X)
  n <- nrow(X)
  time <- .pencox_positive_time(time)
  status <- as.integer(status != 0L)

  alpha_candidates <- sort(unique(as.numeric(
    alpha_candidates[is.finite(alpha_candidates) & alpha_candidates >= 0 & alpha_candidates <= 1]
  )))
  if (length(alpha_candidates) == 0L) {
    stop("At least one alpha candidate in [0, 1] is required.")
  }

  lambda_candidates <- sort(unique(as.numeric(
    lambda_candidates[is.finite(lambda_candidates) & lambda_candidates > 0]
  )))
  if (length(lambda_candidates) == 0L) {
    stop("At least one positive lambda candidate is required.")
  }

  if (is.null(folds)) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    folds <- sample(rep(seq_len(k), length.out = n))
  } else if (length(folds) != n) {
    stop("length(folds) must equal nrow(X).")
  }

  if (is.null(auc_time_grid)) {
    auc_time_grid <- as.numeric(stats::quantile(
      time,
      probs = seq(0.1, 0.8, length.out = 15),
      na.rm = TRUE
    ))
  }
  auc_time_grid <- sort(unique(auc_time_grid[is.finite(auc_time_grid)]))

  results <- data.frame(
    alpha = numeric(),
    lambda = numeric(),
    iAUC = numeric(),
    auc_curve = I(list())
  )
  candidate_auc_curves <- list()

  best_iAUC <- -Inf
  best_alpha <- NA_real_
  best_lambda <- NA_real_
  best_auc_curve <- NULL
  fit_errors <- character()

  for (a in alpha_candidates) {
    for (lam in lambda_candidates) {
      cat("Evaluating penalized Cox alpha =", a, "lambda =", lam, "with", auc_method, "...\n")
      fold_auc_curves <- matrix(NA_real_, nrow = k, ncol = length(auc_time_grid))
      fold_iAUCs <- rep(NA_real_, k)

      for (fold in seq_len(k)) {
        test_idx <- which(folds == fold)
        train_idx <- setdiff(seq_len(n), test_idx)

        model <- tryCatch(
          fit_penalized_cox(
            X = X[train_idx, , drop = FALSE],
            time = time[train_idx],
            status = status[train_idx],
            alpha = a,
            lambda = lam,
            use_preselection = use_preselection,
            p_thresh = p_thresh,
            max_features = max_features,
            standardize = standardize
          ),
          error = function(e) {
            fit_errors <<- c(fit_errors, conditionMessage(e))
            NULL
          }
        )
        if (is.null(model)) {
          next
        }

        marker_lp <- tryCatch(
          predict_penalized_cox(model, X[test_idx, , drop = FALSE], type = "lp"),
          error = function(e) rep(NA_real_, length(test_idx))
        )
        if (all(!is.finite(marker_lp))) {
          next
        }

        auc_vals <- .coxpcapls_auc_curve(
          time = time[test_idx],
          status = status[test_idx],
          marker = marker_lp,
          auc_time_grid = auc_time_grid,
          auc_method = auc_method
        )

        fold_auc_curves[fold, ] <- auc_vals
        fold_iAUCs[fold] <- mean(auc_vals, na.rm = TRUE)
      }

      mean_curve <- colMeans(fold_auc_curves, na.rm = TRUE)
      mean_iAUC <- mean(fold_iAUCs, na.rm = TRUE)
      if (!is.finite(mean_iAUC)) {
        mean_iAUC <- NA_real_
      }

      results <- rbind(
        results,
        data.frame(alpha = a, lambda = lam, iAUC = mean_iAUC, auc_curve = I(list(mean_curve)))
      )
      key <- paste("alpha", a, "lambda", lam, sep = "_")
      candidate_auc_curves[[key]] <- mean_curve

      if (is.finite(mean_iAUC) && mean_iAUC > best_iAUC) {
        best_iAUC <- mean_iAUC
        best_alpha <- a
        best_lambda <- lam
        best_auc_curve <- mean_curve
      }
    }
  }

  if (!is.finite(best_iAUC)) {
    detail <- if (length(fit_errors) > 0L) {
      paste(unique(fit_errors)[seq_len(min(3L, length(unique(fit_errors))))], collapse = "; ")
    } else {
      "all fold-level AUC values were NA; check whether auc_time_grid has enough events and controls within each test fold."
    }
    stop("Penalized Cox CV failed for all parameter combinations. ", detail)
  }

  final_fit <- fit_penalized_cox(
    X = X,
    time = time,
    status = status,
    alpha = best_alpha,
    lambda = best_lambda,
    use_preselection = use_preselection,
    p_thresh = p_thresh,
    max_features = max_features,
    standardize = standardize
  )
  marker_all <- predict_penalized_cox(final_fit, X, type = "lp")
  metrics <- .pencox_metrics(time, status, marker_all)
  metrics$iAUC <- best_iAUC

  list(
    best_model = final_fit$glmnet_fit,
    best_fit = final_fit,
    best_alpha = best_alpha,
    best_lambda = best_lambda,
    best_iAUC = best_iAUC,
    results = results,
    candidate_auc_curves = candidate_auc_curves,
    auc_time_grid = auc_time_grid,
    best_auc_curve = best_auc_curve,
    metrics = metrics
  )
}


val_ridge_cox_cv <- function(X,
                             time,
                             status,
                             k = 5,
                             lambda_candidates = 10 ^ seq(-4, 1, by = 1),
                             auc_time_grid = NULL,
                             use_preselection = FALSE,
                             p_thresh = 0.05,
                             max_features = NULL,
                             auc_method = c("IPCW", "NNE", "KM"),
                             folds = NULL,
                             seed = NULL,
                             standardize = TRUE) {
  val_penalized_cox_cv(
    X = X,
    time = time,
    status = status,
    k = k,
    alpha_candidates = 0,
    lambda_candidates = lambda_candidates,
    auc_time_grid = auc_time_grid,
    use_preselection = use_preselection,
    p_thresh = p_thresh,
    max_features = max_features,
    auc_method = auc_method,
    folds = folds,
    seed = seed,
    standardize = standardize
  )
}


val_elastic_net_cox_cv <- function(X,
                                   time,
                                   status,
                                   k = 5,
                                   alpha_candidates = c(0.25, 0.5, 0.75),
                                   lambda_candidates = 10 ^ seq(-4, 1, by = 1),
                                   auc_time_grid = NULL,
                                   use_preselection = FALSE,
                                   p_thresh = 0.05,
                                   max_features = NULL,
                                   auc_method = c("IPCW", "NNE", "KM"),
                                   folds = NULL,
                                   seed = NULL,
                                   standardize = TRUE) {
  alpha_candidates <- alpha_candidates[alpha_candidates > 0 & alpha_candidates < 1]
  if (length(alpha_candidates) == 0L) {
    stop("Elastic-net Cox requires alpha_candidates strictly between 0 and 1.")
  }

  val_penalized_cox_cv(
    X = X,
    time = time,
    status = status,
    k = k,
    alpha_candidates = alpha_candidates,
    lambda_candidates = lambda_candidates,
    auc_time_grid = auc_time_grid,
    use_preselection = use_preselection,
    p_thresh = p_thresh,
    max_features = max_features,
    auc_method = auc_method,
    folds = folds,
    seed = seed,
    standardize = standardize
  )
}
