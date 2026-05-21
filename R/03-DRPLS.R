

##### 2008 PLSDR referenced Philippe Bastien

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
  stop("Cannot locate R/auc_utils.R; source it before DRPLS.R.")
}

# cox_pls_dr <- function(X, time, status, time2 = NULL,
#                        residual_type = "deviance",
#                        scale_x = TRUE, scale_time = TRUE,
#                        n_components = min(7, ncol(X)),
#                        return_all = FALSE) {
#   
#   #  1. 
#   if (scale_x) {
#     X_scaled <- scale(X)
#     X_scale <- attr(X_scaled, "scaled:scale")
#     X_center <- attr(X_scaled, "scaled:center")
#     X <- as.data.frame(X_scaled)
#   } else {
#     X <- as.data.frame(X)
#     X_scale <- rep(1, ncol(X))
#     X_center <- rep(0, ncol(X))
#   }
#   
#   if (scale_time && is.null(time2)) {
#     time <- scale(time)
#   }
#   
#   #  2. 
#   surv_obj <- if (is.null(time2)) {
#     Surv(time, status)
#   } else {
#     Surv(time, time2, status)
#   }
#   #print('surv_obj')
#   #  3.  Cox 
#   null_model <- coxph(surv_obj ~ 1)
#   deviance_residuals <- residuals(null_model, type = residual_type)
#   #print('deviance_residuals')
#   
#   
#   #  4. PLS
#   pls_fit <- plsr(deviance_residuals ~ as.matrix(X),
#                   ncomp = n_components,
#                   scale = FALSE)
#   #print('pls_fit')
#   A <- pls_fit$projection
#   pls_scores <- as.data.frame(as.matrix(X) %*% A)
#   #pls_scores <- as.data.frame(scores(pls_fit, type="scores")[, 1:n_components, drop=FALSE])
#   colnames(pls_scores) <- paste0("comp_", 1:ncol(pls_scores))
#   #print('pls_scores')
#   #  5. Cox
#   final_cox_model <- coxph(surv_obj ~ ., data = pls_scores, 
#                            control = coxph.control(iter.max = 100), x = TRUE,
#                            y = TRUE)
#   #print('final_cox_model')
#   #  6.  - 
#   if (return_all) {
#     coef_matrix <- matrix(NA, nrow = n_components, ncol = n_components)
#     actual_ncomp <- ncol(pls_scores) 
#     
#     for (i in 1:actual_ncomp) {
#       cols_to_use <- paste0("comp_", 1:i)
#       cox_i <- coxph(surv_obj ~ ., data = pls_scores[, cols_to_use, drop = FALSE], 
#                      control = coxph.control(iter.max = 100), x = TRUE, y = TRUE)
#       coefs <- coef(cox_i)
#       coef_matrix[1:length(coefs), i] <- coefs
#     }
#     
#     return(list(pls_scores = pls_scores,
#                 final_cox_model = final_cox_model,
#                 pls_model = pls_fit,
#                 X_scale = X_scale,
#                 X_center = X_center,
#                 coef_matrix = coef_matrix,
#                 A = A))
#   }
# }


cox_pls_dr <- function(X, time, status, time2 = NULL,
                       residual_type = "deviance",
                       scale_x = TRUE, scale_time = TRUE,
                       n_components = min(7, ncol(X)),
                       return_all = FALSE) {
  
  #  1. 
  if (scale_x) {
    X_scaled  <- scale(X)
    X_scale   <- attr(X_scaled, "scaled:scale")
    X_center  <- attr(X_scaled, "scaled:center")
    X         <- as.data.frame(X_scaled)
  } else {
    X         <- as.data.frame(X)
    X_scale   <- rep(1, ncol(X))
    X_center  <- rep(0, ncol(X))
  }
  
  if (scale_time && is.null(time2)) {
    time <- scale(time)
  }
  
  #  2. 
  surv_obj <- if (is.null(time2)) {
    Surv(time, status)
  } else {
    Surv(time, time2, status)
  }
  
  #  3.  Cox 
  null_model <- coxph(surv_obj ~ 1)
  deviance_residuals <- residuals(null_model, type = residual_type)
  
  #  4.  PLS
  pls_fit <- plsr(deviance_residuals ~ as.matrix(X),
                  ncomp = n_components,
                  scale = FALSE)
  A <- pls_fit$projection
  pls_scores <- as.data.frame(as.matrix(X) %*% A)
  colnames(pls_scores) <- paste0("comp_", 1:ncol(pls_scores))
  
  #  5.  Cox   ----  tryCatch
  final_cox_model <- tryCatch({
    coxph(surv_obj ~ .,
          data    = pls_scores,
          control = coxph.control(iter.max = 100),
          x = TRUE, y = TRUE)
  }, error = function(e) {
    warning("final Cox model (", n_components, " components) failed: ",
            e$message)
    NULL  #  NULL predict  tryCatch 
  })
  
  # 
  if (!return_all) {
    return(list(
      pls_scores = pls_scores,
      final_cox_model = final_cox_model,
      pls_model = pls_fit,
      X_scale = X_scale,
      X_center = X_center,
      A = A
    ))
  }
  
  #  6.  Cox coef_matrix  ----  i  tryCatch
  coef_matrix   <- matrix(NA_real_, nrow = n_components, ncol = n_components)
  actual_ncomp  <- ncol(pls_scores)
  
  for (i in 1:actual_ncomp) {
    cols_to_use <- paste0("comp_", 1:i)
    
    cox_i <- tryCatch({
      coxph(surv_obj ~ .,
            data    = pls_scores[, cols_to_use, drop = FALSE],
            control = coxph.control(iter.max = 100),
            x = TRUE, y = TRUE)
    }, error = function(e) {
      warning("Cox model with ", i, " component(s) failed: ", e$message)
      NULL
    })
    
    if (!is.null(cox_i)) {
      coefs <- coef(cox_i)
      if (length(coefs) > 0) {
        coef_matrix[1:length(coefs), i] <- coefs
      }
    }
    #  cox_i  NULL NA
  }
  
  return(list(
    pls_scores     = pls_scores,
    final_cox_model = final_cox_model,  #  NULL
    pls_model      = pls_fit,
    X_scale        = X_scale,
    X_center       = X_center,
    coef_matrix    = coef_matrix,
    A              = A
  ))
}


# validate the model using cv "NNE"
# val_DRPLS_cv <- function(X, time, status,
#                          k = 5,
#                          ncomp_candidates = 1:5,
#                          auc_time_grid = NULL) {
#   n <- nrow(X)
#   set.seed(123)
#   folds <- sample(rep(1:k, length.out = n))
#   
#   # if (is.null(auc_time_grid)) {
#   #   auc_time_grid <- seq(min(time), max(time), length.out = 10)
#   # }
#   
#   if (is.null(auc_time_grid)) {
#     #  time  10%  90% 10
#     auc_time_grid <- quantile(time, probs = seq(0.1, 0.8, length.out = 15))
#   }
#   
#   results <- data.frame(ncomp = integer(0), iAUC = numeric(0), auc_curve = I(list()))
#   candidate_auc_curves <- list()
#   
#   best_iAUC <- -Inf
#   best_ncomp <- NA
#   best_model <- NULL
#   best_A <- NULL
#   best_auc_curve <- NULL
#   
#   for (nc in ncomp_candidates) {
#     cat("Evaluating ncomp =", nc, "\n")
#     
#     #  AUC  iAUC
#     fold_auc_curves <- matrix(NA, nrow = k, ncol = length(auc_time_grid))
#     fold_iAUCs <- numeric(k)
#     
#     for (f in 1:k) {
#       test_idx <- which(folds == f)
#       train_idx <- setdiff(1:n, test_idx)
#       
#       model <- cox_pls_dr(X = X[train_idx, ],
#                           time = time[train_idx],
#                           status = status[train_idx],
#                           n_components = nc,
#                           return_all = TRUE)
#       
#       A <- model$A
#       T_test <- as.data.frame(as.matrix(X[test_idx, ]) %*% A)
#       colnames(T_test) <- paste0("comp_", 1:nc)
#       marker_lp <- predict(model$final_cox_model, newdata = T_test, type = "lp")
#       
#       #  AUC 
#       # auc_vals <- sapply(auc_time_grid, function(tp) {
#       #   roc_obj <- survivalROC(Stime = time[test_idx], status = status[test_idx],
#       #                          marker = marker_lp, predict.time = tp, method = "NNE",
#       #                          span = 0.3)
#       #   return(roc_obj$AUC)
#       # })
#       
#       auc_vals <- sapply(auc_time_grid, function(tp) {
#         res <- tryCatch({
#           roc_obj <- survivalROC(
#             Stime        = time[test_idx],
#             status       = status[test_idx],
#             marker       = marker_lp,
#             predict.time = tp,
#             method       = "NNE",
#             span         = 0.30
#           )
#           roc_obj$AUC
#         }, error = function(e) {
#           warning(paste("tp =", tp, ":", e$message))
#           return(NA)
#         })
#         return(res)
#       })
#       
#       
#       fold_auc_curves[f, ] <- auc_vals
#       fold_iAUCs[f] <- mean(auc_vals, na.rm = TRUE)
#     }
#     
#     #  ncomp  AUC  iAUC
#     mean_auc_curve <- colMeans(fold_auc_curves, na.rm = TRUE)
#     mean_iAUC <- mean(fold_iAUCs, na.rm = TRUE)
#     
#     results <- rbind(results, data.frame(ncomp = nc, iAUC = mean_iAUC,
#                                          auc_curve = I(list(mean_auc_curve))))
#     candidate_auc_curves[[as.character(nc)]] <- mean_auc_curve
#     
#     if (mean_iAUC > best_iAUC) {
#       best_iAUC <- mean_iAUC
#       best_ncomp <- nc
#       best_model <- model$final_cox_model
#       best_A <- A
#       best_auc_curve <- mean_auc_curve
#     }
#   }
#   
#   return(list(best_model = best_model,
#               best_ncomp = best_ncomp,
#               best_iAUC = best_iAUC,
#               results = results,
#               best_A = best_A,
#               best_auc_curve = best_auc_curve,
#               candidate_auc_curves = candidate_auc_curves,
#               auc_time_grid = auc_time_grid))
# }
# 
# 
# # "KM"
# val_DRPLS_cv <- function(X, time, status,
#                          k = 5,
#                          ncomp_candidates = 1:5,
#                          auc_time_grid = NULL) {
#   n <- nrow(X)
#   set.seed(123)
#   folds <- sample(rep(1:k, length.out = n))
#   
#   # if (is.null(auc_time_grid)) {
#   #   auc_time_grid <- seq(min(time), max(time), length.out = 10)
#   # }
#   
#   if (is.null(auc_time_grid)) {
#     auc_time_grid <- quantile(time, probs = seq(0.1, 0.8, length.out = 10))
#   }
#   
#   results <- data.frame(ncomp = integer(0), iAUC = numeric(0), auc_curve = I(list()))
#   candidate_auc_curves <- list()
#   
#   best_iAUC <- -Inf
#   best_ncomp <- NA
#   best_model <- NULL
#   best_A <- NULL
#   best_auc_curve <- NULL
#   
#   for (nc in ncomp_candidates) {
#     cat("Evaluating ncomp =", nc, "\n")
#     
#     #  AUC  iAUC
#     fold_auc_curves <- matrix(NA, nrow = k, ncol = length(auc_time_grid))
#     fold_iAUCs <- numeric(k)
#     
#     for (f in 1:k) {
#       test_idx <- which(folds == f)
#       train_idx <- setdiff(1:n, test_idx)
#       
#       model <- cox_pls_dr(X = X[train_idx, ],
#                           time = time[train_idx],
#                           status = status[train_idx],
#                           n_components = nc,
#                           return_all = TRUE)
#       
#       A <- model$A
#       T_test <- as.data.frame(as.matrix(X[test_idx, ]) %*% A)
#       colnames(T_test) <- paste0("comp_", 1:nc)
#       marker_lp <- predict(model$final_cox_model, newdata = T_test, type = "lp")
# 
#       auc_vals <- sapply(auc_time_grid, function(tp) {
#         roc_obj <- survivalROC(Stime = time[test_idx], status = status[test_idx],
#                                marker = marker_lp, predict.time = tp, method = "KM")
#         return(roc_obj$AUC)
#       })
#       
#       
#       fold_auc_curves[f, ] <- auc_vals
#       fold_iAUCs[f] <- mean(auc_vals, na.rm = TRUE)
#     }
#     
#     #  ncomp  AUC  iAUC
#     mean_auc_curve <- colMeans(fold_auc_curves, na.rm = TRUE)
#     mean_iAUC <- mean(fold_iAUCs, na.rm = TRUE)
#     
#     results <- rbind(results, data.frame(ncomp = nc, iAUC = mean_iAUC,
#                                          auc_curve = I(list(mean_auc_curve))))
#     candidate_auc_curves[[as.character(nc)]] <- mean_auc_curve
#     
#     if (mean_iAUC > best_iAUC) {
#       best_iAUC <- mean_iAUC
#       best_ncomp <- nc
#       best_model <- model$final_cox_model
#       best_A <- A
#       best_auc_curve <- mean_auc_curve
#     }
#   }
#   
#   return(list(best_model = best_model,
#               best_ncomp = best_ncomp,
#               best_iAUC = best_iAUC,
#               results = results,
#               best_A = best_A,
#               best_auc_curve = best_auc_curve,
#               candidate_auc_curves = candidate_auc_curves,
#               auc_time_grid = auc_time_grid))
# }


# pre-selesct features for NNE
# val_DRPLS_cv <- function(
#     X, time, status,
#     k = 5,
#     ncomp_candidates = 1:5,
#     auc_time_grid = NULL,
#     use_preselection = FALSE,
#     p_thresh = 0.05,
#     auc_method = c("NNE", "KM")
# ) {
#   #  auc_method
#   auc_method <- match.arg(auc_method)
#   
#   # --- Step 0 ---
#   if (use_preselection) {
#     p <- ncol(X)
#     significant_idx <- integer(0)
#     for (j in seq_len(p)) {
#       fit <- tryCatch({
#         survival::coxph(survival::Surv(time, status) ~ X[, j])
#       }, error = function(e) NULL)
#       if (!is.null(fit)) {
#         pval <- summary(fit)$coefficients[ , "Pr(>|z|)"]
#         if (!is.na(pval) && pval < p_thresh) {
#           significant_idx <- c(significant_idx, j)
#         }
#       }
#     }
#     if (length(significant_idx) == 0) {
#       stop("No variables pass univariate Cox at p <", p_thresh)
#     }
#     X <- X[, significant_idx, drop = FALSE]
#     message(length(significant_idx),
#             " predictors retained after preselection (p < ", p_thresh, ").")
#   }
#   
#   # ---  ---
#   n <- nrow(X)
#   set.seed(123)
#   folds <- sample(rep(1:k, length.out = n))
#   
#   if (is.null(auc_time_grid)) {
#     auc_time_grid <- quantile(time,
#                               probs = seq(0.1, 0.8, length.out = 15))
#   }
#   
#   results <- data.frame(
#     ncomp     = integer(0),
#     iAUC      = numeric(0),
#     auc_curve = I(list())
#   )
#   candidate_auc_curves <- list()
#   best_iAUC   <- -Inf
#   best_ncomp  <- NA
#   best_model  <- NULL
#   best_A      <- NULL
#   best_curve  <- NULL
#   
#   for (nc in ncomp_candidates) {
#     cat("Evaluating ncomp =", nc, "with", auc_method, "...\n")
#     
#     fold_auc_curves <- matrix(NA, nrow = k, ncol = length(auc_time_grid))
#     fold_iAUCs     <- numeric(k)
#     
#     for (f in seq_len(k)) {
#       test_idx  <- which(folds == f)
#       train_idx <- setdiff(seq_len(n), test_idx)
#       
#       model <- cox_pls_dr(
#         X            = X[train_idx, , drop = FALSE],
#         time         = time[train_idx],
#         status       = status[train_idx],
#         n_components = nc,
#         return_all   = TRUE
#       )
#       
#       A        <- model$A
#       T_test   <- as.data.frame(as.matrix(X[test_idx, , drop = FALSE]) %*% A)
#       colnames(T_test) <- paste0("comp_", seq_len(nc))
#       
#       marker_lp <- predict(model$final_cox_model,
#                            newdata = T_test,
#                            type   = "lp")
#       
#       #   auc_method   
#       if (auc_method == "NNE") {
#         auc_vals <- sapply(auc_time_grid, function(tp) {
#           tryCatch({
#             roc_obj <- survivalROC::survivalROC(
#               Stime        = time[test_idx],
#               status       = status[test_idx],
#               marker       = marker_lp,
#               predict.time = tp,
#               method       = "NNE",
#               span         = 0.30
#             )
#             roc_obj$AUC
#           }, error = function(e) {
#             warning("tp =", tp, " error: ", e$message)
#             NA
#           })
#         })
#       } else {  # auc_method == "KM"
#         auc_vals <- sapply(auc_time_grid, function(tp) {
#           roc_obj <- survivalROC::survivalROC(
#             Stime        = time[test_idx],
#             status       = status[test_idx],
#             marker       = marker_lp,
#             predict.time = tp,
#             method       = "KM"
#           )
#           roc_obj$AUC
#         })
#       }
#       
#       fold_auc_curves[f, ] <- auc_vals
#       fold_iAUCs[f]       <- mean(auc_vals, na.rm = TRUE)
#     }
#     
#     mean_curve <- colMeans(fold_auc_curves, na.rm = TRUE)
#     mean_iAUC  <- mean(fold_iAUCs, na.rm = TRUE)
#     
#     results <- rbind(
#       results,
#       data.frame(ncomp = nc, iAUC = mean_iAUC, auc_curve = I(list(mean_curve)))
#     )
#     candidate_auc_curves[[as.character(nc)]] <- mean_curve
#     
#     if (mean_iAUC > best_iAUC) {
#       best_iAUC  <- mean_iAUC
#       best_ncomp <- nc
#       best_model <- model$final_cox_model
#       best_A     <- A
#       best_curve <- mean_curve
#     }
#   }
#   
#   list(
#     best_model         = best_model,
#     best_ncomp         = best_ncomp,
#     best_iAUC          = best_iAUC,
#     results            = results,
#     best_A             = best_A,
#     best_auc_curve     = best_curve,
#     candidate_auc_curves = candidate_auc_curves,
#     auc_time_grid      = auc_time_grid
#   )
# }


# multi-index
val_DRPLS_cv <- function(
    X, time, status,
    k = 5,
    ncomp_candidates = 1:5,
    auc_time_grid = NULL,
    use_preselection = FALSE,
    p_thresh = 0.05,
    auc_method = c("IPCW", "NNE", "KM"),
    folds = NULL,
    seed = 123
) {
  
  auc_method <- match.arg(auc_method)
  
  n <- nrow(X)
  
  if (is.null(folds)) {
    if (!is.null(seed)) set.seed(seed)
    folds <- sample(rep(1:k, length.out = n))
  } else {
    # 
    if (length(folds) != n)
      stop("length(folds) (", length(folds), ") != nrow(X) (", n, ")")
    if (!all(folds %in% 1:k))
      stop("folds  1:k ")
  }
  
  # --- Step 0 ---
  if (use_preselection) {
    p <- ncol(X)
    significant_idx <- integer(0)
    for (j in seq_len(p)) {
      fit <- tryCatch({
        survival::coxph(survival::Surv(time, status) ~ X[, j])
      }, error = function(e) NULL)
      if (!is.null(fit)) {
        pval <- summary(fit)$coefficients[, "Pr(>|z|)"]
        if (!is.na(pval) && pval < p_thresh) {
          significant_idx <- c(significant_idx, j)
        }
      }
    }
    if (length(significant_idx) == 0) {
      stop("No variables pass univariate Cox at p <", p_thresh)
    }
    X <- X[, significant_idx, drop = FALSE]
    message(length(significant_idx),
            " predictors retained after preselection (p < ", p_thresh, ").")
  }
  
  if (is.null(auc_time_grid)) {
    auc_time_grid <- quantile(time, probs = seq(0.1, 0.8, length.out = 15))
  }
  
  results <- data.frame(
    ncomp     = integer(0),
    iAUC      = numeric(0),
    auc_curve = I(list())
  )
  candidate_auc_curves <- list()
  
  best_iAUC  <- -Inf
  best_ncomp <- NA
  best_model <- NULL
  best_A     <- NULL
  best_curve <- NULL
  
  # ---  ---
  for (nc in ncomp_candidates) {
    cat("Evaluating ncomp =", nc, "with", auc_method, "...\n")
    
    fold_auc_curves <- matrix(NA, nrow = k, ncol = length(auc_time_grid))
    fold_iAUCs     <- numeric(k)
    
    for (f in seq_len(k)) {
      test_idx  <- which(folds == f)
      train_idx <- setdiff(seq_len(n), test_idx)
      
      model <- cox_pls_dr(
        X            = X[train_idx, , drop = FALSE],
        time         = time[train_idx],
        status       = status[train_idx],
        n_components = nc,
        return_all   = TRUE
      )
      
      A        <- model$A
      T_test   <- as.data.frame(as.matrix(X[test_idx, , drop = FALSE]) %*% A)
      
      # ---  ---
      if (!is.null(model$final_cox_model$coefficients)) {
        colnames(T_test) <- names(model$final_cox_model$coefficients)
      } else {
        colnames(T_test) <- paste0("comp_", seq_len(nc))
      }
      
      marker_lp <- tryCatch({
        predict(model$final_cox_model,
                newdata = T_test,
                type = "lp")
      }, error = function(e) {
        warning("Fold ", f, ": prediction failed -> ", e$message)
        rep(NA_real_, nrow(T_test))
      })
      
      auc_vals <- .coxpcapls_auc_curve(
        time = time[test_idx],
        status = status[test_idx],
        marker = marker_lp,
        auc_time_grid = auc_time_grid,
        auc_method = auc_method
      )
      
      fold_auc_curves[f, ] <- auc_vals
      fold_iAUCs[f]       <- mean(auc_vals, na.rm = TRUE)
    }
    
    mean_curve <- colMeans(fold_auc_curves, na.rm = TRUE)
    mean_iAUC  <- mean(fold_iAUCs, na.rm = TRUE)
    
    if (!any(is.finite(fold_iAUCs))) {
      mean_iAUC <- NA_real_
    }
    
    results <- rbind(
      results,
      data.frame(ncomp = nc, iAUC = mean_iAUC, auc_curve = I(list(mean_curve)))
    )
    candidate_auc_curves[[as.character(nc)]] <- mean_curve
    
    if (is.finite(mean_iAUC) && mean_iAUC > best_iAUC) {
      best_iAUC  <- mean_iAUC
      best_ncomp <- nc
      best_model <- model$final_cox_model
      best_A     <- A
      best_curve <- mean_curve
    }
  }
  
  # --- Step 2:  ---
  if (!is.finite(best_iAUC) || is.na(best_ncomp)) {
    stop("val_DRPLS_cv failed for all component candidates.")
  }
  
  final_model <- cox_pls_dr(
    X            = X,
    time         = time,
    status       = status,
    n_components = best_ncomp,
    return_all   = TRUE
  )
  
  T_all <- as.data.frame(as.matrix(X) %*% final_model$A)
  
  # ---  comp_1 ---
  if (!is.null(final_model$final_cox_model$coefficients)) {
    colnames(T_all) <- names(final_model$final_cox_model$coefficients)
  } else {
    colnames(T_all) <- paste0("comp_", seq_len(best_ncomp))
  }
  
  marker_all <- tryCatch({
    predict(final_model$final_cox_model,
            newdata = T_all,
            type = "lp")
  }, error = function(e) {
    warning("Final model prediction failed: ", e$message)
    rep(NA_real_, nrow(T_all))
  })
  
  # --- Step 3:  ---
  cindex_harrell <- survcomp::concordance.index(
    x = marker_all, surv.time = time, surv.event = status, method = "noether"
  )$c.index
  
  cindex_uno <- tryCatch({
    survAUC::UnoC(
      Surv.rsp = Surv(time, status),
      Surv.rsp.new = Surv(time, status),
      lpnew = marker_all
    )$C
  }, error = function(e) {
    warning("UnoC failed: ", e$message)
    NA
  })
  
  # =====================  IPCW-Brier pec =====================

  cat(" Brier ...\n")
  
  # 0) 
  df_brier <- cbind(T_all, time = time, status = status)
  comp_names <- paste0("comp_", seq_len(ncol(T_all)))
  colnames(df_brier)[1:ncol(T_all)] <- comp_names
  df_brier <- as.data.frame(df_brier, stringsAsFactors = FALSE, check.names = FALSE)
  df_brier$time   <- as.numeric(df_brier$time)
  df_brier$status <- as.integer(df_brier$status != 0L)
  stopifnot(all(df_brier$status %in% c(0L,1L)))
  
  # 1) 
  event_times <- df_brier$time[df_brier$status == 1L]
  if (length(event_times) >= 5L) {
    times_use <- as.numeric(stats::quantile(event_times, probs = seq(0.1, 0.9, length.out = 15), na.rm = TRUE))
  } else {
    times_use <- as.numeric(stats::quantile(df_brier$time, probs = seq(0.1, 0.9, length.out = 15), na.rm = TRUE))
  }
  max_event_time <- suppressWarnings(max(event_times, na.rm = TRUE))
  times_use <- unique(times_use[!is.na(times_use) & times_use > 0])
  if (is.finite(max_event_time)) times_use <- sort(times_use[times_use <= max_event_time])
  if (length(times_use) == 0L) {
    alt_t <- min(stats::median(df_brier$time, na.rm = TRUE), max(df_brier$time, na.rm = TRUE))
    times_use <- alt_t[is.finite(alt_t)]
  }
  
  cat(" :", nrow(df_brier), ",", ncol(df_brier), "\n")
  cat(" :", min(df_brier$time), "-", max(df_brier$time), "\n")
  cat(" :", sum(df_brier$status), "\n")
  cat(" time grid(<=) :", length(times_use), "\n")
  
  # 2)  df_brier  Cox x/y
  cox_for_brier <- tryCatch({
    survival::coxph(
      formula = as.formula(paste0("survival::Surv(time, status) ~ ", paste(comp_names, collapse = " + "))),
      data = df_brier, x = TRUE, y = TRUE, ties = "breslow", na.action = na.exclude
    )
  }, error = function(e) { stop(" Cox : ", e$message)
    NA})
  
  # 3)  times_use  S_i(t)
  n  <- nrow(df_brier)
  Tn <- length(times_use)
  S_mat <- matrix(NA_real_, nrow = n, ncol = Tn)
  for (i in seq_len(n)) {
    sfi <- survival::survfit(cox_for_brier, newdata = df_brier[i, comp_names, drop = FALSE])
    si  <- summary(sfi, times = times_use, extend = TRUE)
    if (length(si$surv) == Tn) {
      S_mat[i, ] <- si$surv
    } else {
      # 
      S_mat[i, ] <- approx(x = si$time, y = si$surv, xout = times_use, rule = 2, ties = "ordered")$y
    }
  }
  
  # 4)  prodlim  G(t) G(t)G(T_i-)
  fitG <- prodlim::prodlim(prodlim::Hist(time, 1L - status) ~ 1, data = df_brier)
  predict_prodlim <- getS3method("predict", "prodlim")   #  predict.prodlim  S3 
  G_t      <- as.numeric(predict_prodlim(fitG, times = times_use,                  type = "surv"))
  eps <- 1e-8
  G_Tminus <- as.numeric(predict_prodlim(fitG, times = pmax(df_brier$time - eps, 0), type = "surv"))
  
  # /
  G_t[G_t <= 1e-8 | !is.finite(G_t)]         <- NA_real_
  G_Tminus[G_Tminus <= 1e-8 | !is.finite(G_Tminus)] <- NA_real_
  
  # 5) IPCW-Brier 
  TT <- df_brier$time
  DD <- df_brier$status
  briers <- numeric(Tn)
  for (j in seq_len(Tn)) {
    t0 <- times_use[j]
    # I(T>t)/G(t) * (0 - S_i(t))^2
    term1 <- (TT >  t0) * ((S_mat[, j])^2) / G_t[j]
    # I(T<=t, =1)/G(T_i-) * (1 - S_i(t))^2
    term2 <- (TT <= t0 & DD == 1L) * ((1 - S_mat[, j])^2) / G_Tminus
    briers[j] <- mean(c(term1, term2), na.rm = TRUE)
  }
  brier_mean <- mean(briers, na.rm = TRUE)
  
  cat(" Brier :", brier_mean, "\n")
  #  `briers`  `times_use` 
  #  brier_curve = briers, brier_times = times_use
  # =====================  IPCW-Brier  =====================
  
    
  metrics <- list(
    Harrell_Cindex = cindex_harrell,
    Uno_Cindex     = if (!is.null(cindex_uno)) cindex_uno else NA,
    iAUC           = best_iAUC,
    Brier          = if (!is.na(brier_mean)) brier_mean else NA
  )
  
  
  # --- Step 4:  ---
  list(
    best_model          = final_model$final_cox_model,
    best_ncomp          = best_ncomp,
    best_iAUC           = best_iAUC,
    results             = results,
    best_A              = final_model$A,
    best_auc_curve      = best_curve,
    candidate_auc_curves = candidate_auc_curves,
    auc_time_grid       = auc_time_grid,
    metrics             = metrics
  )
}









# val_DRPLS_cv <- function(X, time, status,
#                          k = 5,
#                          ncomp_candidates = 1:5,
#                          auc_time_grid = NULL) {
#   # X: 
#   # time: 
#   # status: 1 0 
#   # k: 
#   # ncomp_candidates: 
#   # auc_time_grid:  ROC 
#   #
#   #  ncomp iAUC 
# 
#   n <- nrow(X)
#   #  k 
#   set.seed(123)  # 
#   folds <- sample(rep(1:k, length.out = n))
# 
#   #  auc_time_grid
#   #  min  max 
#   # 
#   if (is.null(auc_time_grid)) {
#     auc_time_grid <- seq(min(time), max(time), length.out = 50)
#     # 
#     # auc_time_grid <- seq(quantile(time, 0.1), quantile(time, 0.9), length.out = 50)
#   }
# 
#   #  AUC  AUC 
#   results <- data.frame(ncomp = integer(0), iAUC = numeric(0), auc_curve = I(list()))
#   candidate_auc_curves <- list()
# 
#   best_iAUC <- -Inf
#   best_ncomp <- NA
#   best_model <- NULL
#   best_A <- NULL
#   best_auc_curve <- NULL
# 
#   # 
#   for (nc in ncomp_candidates) {
#     # 
#     pooled_time <- c()
#     pooled_status <- c()
#     pooled_marker <- c()
# 
#     #  k 
#     for(f in 1:k) {
#       test_idx <- which(folds == f)
#       train_idx <- setdiff(1:n, test_idx)
# 
#       #  cox_pls_dr 
#       #  final_cox_model Cox  A
#       model <- cox_pls_dr(X = X[train_idx, ],
#                           time = time[train_idx],
#                           status = status[train_idx],
#                           n_components = nc,
#                           return_all = TRUE)
# 
#       A <- model$A  # 
# 
#       # 
#       T_test <- as.data.frame(as.matrix(X[test_idx, ]) %*% A)
#       colnames(T_test) <- paste0("comp_", 1:nc)
# 
#       #  Cox  marker_lp
#       marker_lp <- predict(model$final_cox_model, newdata = T_test, type = "lp")
# 
#       # 
#       pooled_time <- c(pooled_time, time[test_idx])
#       pooled_status <- c(pooled_status, status[test_idx])
#       pooled_marker <- c(pooled_marker, marker_lp)
#     }  # end fold loop
# 
#     #  AUC 
#     auc_vals <- sapply(auc_time_grid, function(tp) {
#       roc_obj <- survivalROC(Stime = pooled_time, status = pooled_status,
#                              marker = pooled_marker, predict.time = tp, method = "KM")
#       return(roc_obj$AUC)
#     })
# 
#     #  AUC AUC 
#     iAUC <- mean(auc_vals, na.rm = TRUE)
# 
#     #  results 
#     results <- rbind(results, data.frame(ncomp = nc, iAUC = iAUC, auc_curve = I(list(auc_vals))))
#     candidate_auc_curves[[as.character(nc)]] <- auc_vals
# 
#     #  iAUC 
#     if (iAUC > best_iAUC) {
#       best_iAUC <- iAUC
#       best_ncomp <- nc
#       #  best_model  best_A 
#       # 
#       # ncomp
#       best_model <- model$final_cox_model
#       best_A <- A
#       best_auc_curve <- auc_vals
#     }
#   }  # end candidate loop
# 
#   return(list(best_model = best_model,
#               best_ncomp = best_ncomp,
#               best_iAUC = best_iAUC,
#               results = results,
#               best_A = best_A,
#               best_auc_curve = best_auc_curve,
#               candidate_auc_curves = candidate_auc_curves,
#               auc_time_grid = auc_time_grid))
# }



# val_DRPLS <- function(X_train, time_train, status_train,
#                       X_test, time_test, status_test,
#                       ncomp_candidates = 1:5,
#                       auc_time_grid = NULL) {
# 
#   if (is.null(auc_time_grid)) {
#     auc_time_grid <- seq(min(c(time_train, time_test)), max(c(time_train, time_test)), length.out = 50)
#   }
#   
# 
#   #  results  ncompiAUC  auc_curve
#   results <- data.frame(ncomp = integer(0), iAUC = numeric(0), auc_curve = I(list()))
#   candidate_auc_curves <- list()  #  auc_time_grid  AUC 
#   best_iAUC <- -Inf
#   best_ncomp <- NA
#   best_model <- NULL
#   best_A <- NULL
#   best_auc_curve <- NULL
# 
#   # 
#   for (nc in ncomp_candidates) {
#     #  cox_pls_dr  final_cox_model  A
#     model <- cox_pls_dr(X = X_train, time = time_train, status = status_train,
#                         n_components = nc, return_all = TRUE)
#     A <- model$A  # 
# 
#     # 
#     T_test <- as.data.frame(as.matrix(X_test) %*% A)
#     colnames(T_test) <- paste0("comp_", 1:nc)
# 
#     #  Cox 
#     marker_lp <- predict(model$final_cox_model, newdata = T_test, type = "lp")
# 
#     #  AUC
#     auc_vals <- numeric(length(auc_time_grid))
#     for (i in seq_along(auc_time_grid)) {
#       tp <- auc_time_grid[i]
#       roc_obj <- survivalROC(Stime = time_test, status = status_test,
#                              marker = marker_lp, predict.time = tp, method = "KM")
#       auc_vals[i] <- roc_obj$AUC
#       print(paste0("ncomp: ", nc, ", time: ", tp, ", AUC: ", auc_vals[i]))
#     }
#     #  AUC AUC 
#     iAUC <- mean(auc_vals, na.rm = TRUE)
#     #print(paste0("ncomp: ", nc, ", iAUC: ", iAUC))
#     #  results  auc_curve
#     results <- rbind(results, data.frame(ncomp = nc, iAUC = iAUC, auc_curve = I(list(auc_vals))))
#     candidate_auc_curves[[as.character(nc)]] <- auc_vals
# 
#     #  AUC 
#     if (iAUC > best_iAUC) {
#       best_iAUC <- iAUC
#       best_ncomp <- nc
#       best_model <- model$final_cox_model
#       best_A <- A
#       best_auc_curve <- auc_vals
#     }
#   }
# 
#   return(list(best_model = best_model,
#               best_ncomp = best_ncomp,
#               best_iAUC = best_iAUC,
#               results = results,
#               best_A = best_A,
#               best_auc_curve = best_auc_curve,
#               candidate_auc_curves = candidate_auc_curves,
#               auc_time_grid = auc_time_grid))
# }



# :
if (FALSE) {
  # 1. 
  set.seed(123)
  n <- 100
  p <- 10
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  time <- rexp(n, 0.1)
  status <- sample(0:1, n, replace = TRUE)
  result <- cox_pls_dr(X = X, time = time, status = status, n_components = 5, return_all = TRUE)
  
  summary(result$final_cox_model)
  result$c_index
  
  print(dim(X))
  
  set.seed(123)
  
  # 
  n <- 200      # 
  p <- 15       # 
  
  # X ()
  Sigma <- diag(p)
  Sigma[1:5, 1:5] <- 0.5  # 5
  diag(Sigma) <- 1
  
  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  
  #  ()
  beta_true <- c(1.5, -2, 0.8, 0, 0, rep(0, p - 5))
  
  # Weibull
  baseline_hazard <- 0.01
  shape_param <- 1.5  # Weibull
  
  # 
  linear_pred <- X %*% beta_true
  
  # Weibull 
  U <- runif(n)
  time <- (-log(U) / (baseline_hazard * exp(linear_pred)))^(1 / shape_param)
  
  # 
  censor_time <- rexp(n, rate = 0.01) # 
  
  # 
  observed_time <- pmin(time, censor_time)
  status <- as.numeric(time <= censor_time)
  
  # 
  cat(":", mean(status == 0), "\n")
  
  # 
  simulated_data <- data.frame(time = observed_time, status = status, X)
  
  result <- cox_pls_dr(
    simulated_data[, -(1:2)],
    time = simulated_data$time,
    status = simulated_data$status,
    n_components = 5,
    return_all = TRUE
  )
  
  # 
  summary(result$final_cox_model)
}
