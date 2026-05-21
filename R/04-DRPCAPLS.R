

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
  stop("Cannot locate R/auc_utils.R; source it before DRPCAPLS.R.")
}

if (!exists("acc_pca_pls_with_alpha", mode = "function")) {
  helper_path <- file.path(getwd(), "pcapls_CR.R")
  if (file.exists(helper_path)) {
    source(helper_path)
  }
}

if (!exists("acc_pca_pls_with_alpha", mode = "function")) {
  stop("`acc_pca_pls_with_alpha()` is unavailable. Please load `pcapls_CR.R` first.")
}

##### add PCA information to the PLSDR
# cox_pca_pls_dr <- function(X, time, status, time2 = NULL,
#                        residual_type = "deviance",
#                        scale_time = TRUE,
#                        n_components = min(7, ncol(X)),
#                        alpha = 0.5,  # alpha
#                        return_all = FALSE) {
# 
#   
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
#   
#   #  3.  Cox 
#   null_model <- coxph(surv_obj ~ 1)
#   deviance_residuals <- residuals(null_model, type = residual_type)
# 
# 
#   #  4. 
#   #  pca_pls_with_alpha
#   pls_fit <- acc_pca_pls_with_alpha(X = as.matrix(X), 
#                                 Y = deviance_residuals,
#                                 ncomp = n_components, 
#                                 alpha = alpha)
# 
#   #  scores
#   A <- pls_fit$Projection
#   pls_scores <- as.data.frame(as.matrix(X) %*% A)
#   #pls_scores <- as.data.frame(pls_fit$scores[, 1:n_components, drop=FALSE])
#   colnames(pls_scores) <- paste0("comp_", 1:ncol(pls_scores))
# 
#   #  5.  Cox 
#   # final_cox_model <- coxph(surv_obj ~ ., data = pls_scores, x = TRUE, y = TRUE)
# 
#   final_cox_model <- tryCatch({
#     coxph(surv_obj ~ .,
#           data    = pls_scores,
#           control = coxph.control(iter.max = 100),
#           x = TRUE, y = TRUE)
#   }, error = function(e) {
#     warning("final Cox model (", n_components, " components) failed: ",
#             e$message)
#     NULL  #  NULL predict  tryCatch 
#   })
#   
#   #  6.  - 
#   if (return_all) {
#     coef_matrix <- matrix(NA, nrow = n_components, ncol = n_components)
#     for (i in 1:n_components) {
#       cox_i <- coxph(surv_obj ~ ., 
#                      data = pls_scores[, 1:i, drop = FALSE], x = TRUE, y = TRUE)
#       coefs <- coef(cox_i)
#       coef_matrix[1:length(coefs), i] <- coefs
#     }
#     
#     return(list(pls_scores = pls_scores,
#                 final_cox_model = final_cox_model,
#                 pls_model = pls_fit,
#                 X = X,
#                 coef_matrix = coef_matrix,
#                 A = A))
#   } else {
#     return(final_cox_model)
#   }
# }


cox_pca_pls_dr <- function(X, time, status, time2 = NULL,
                           residual_type = "deviance",
                           scale_time = TRUE,
                           n_components = min(7, ncol(X)),
                           alpha = 0.5,
                           normalize_alpha_terms = TRUE,
                           return_all = FALSE,
                           var_eps = 1e-10) {
  
  # ---- 0)  ----
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  stopifnot(length(time) == n, length(status) == n)
  
  #  scale 
  if (scale_time && is.null(time2)) {
    time <- as.numeric(scale(time))
  } else {
    time <- as.numeric(time)
  }
  status <- as.integer(status != 0L)
  
  #  Surv
  surv_obj <- if (is.null(time2)) {
    survival::Surv(time, status)
  } else {
    survival::Surv(time, time2, status)
  }
  
  # ---- 1) null Cox +  ----
  null_model <- tryCatch(
    survival::coxph(surv_obj ~ 1),
    error = function(e) stop("null Cox failed: ", e$message)
  )
  deviance_residuals <- residuals(null_model, type = residual_type)
  
  # ---- 2)  alpha- ----
  pls_fit <- acc_pca_pls_with_alpha(
    X = X,
    Y = deviance_residuals,
    ncomp = n_components,
    alpha = alpha,
    normalize_terms = normalize_alpha_terms
  )
  
  A <- pls_fit$Projection
  if (is.null(A)) stop("acc_pca_pls_with_alpha did not return Projection.")
  
  pls_scores <- as.data.frame(X %*% A)
  colnames(pls_scores) <- paste0("comp_", seq_len(ncol(pls_scores)))
  
  # ---- 3)  Cox  ----
  safe_cox_on_scores <- function(scores_df) {
    # 
    n_event <- sum(status == 1L, na.rm = TRUE)
    if (!is.finite(n_event) || n_event < 3L) return(NULL)
    
    df <- as.data.frame(scores_df)
    
    # 
    df[] <- lapply(df, function(x) {
      x <- as.numeric(x)
      x[!is.finite(x)] <- NA_real_
      x
    })
    
    # 
    v <- sapply(df, var, na.rm = TRUE)
    keep <- which(is.finite(v) & v > var_eps)
    if (length(keep) == 0) return(NULL)
    df <- df[, keep, drop = FALSE]
    
    #  < 
    if (ncol(df) >= n_event) {
      df <- df[, seq_len(max(1, n_event - 1)), drop = FALSE]
    }
    
    fit <- tryCatch(
      survival::coxph(
        surv_obj ~ .,
        data = df,
        x = TRUE, y = TRUE,
        ties = "breslow",
        na.action = na.exclude,
        control = survival::coxph.control(
          iter.max = 100,
          toler.chol = 1e-10
        )
      ),
      error = function(e) NULL
    )
    fit
  }
  
  # ---- 4)  Cox ----
  final_cox_model <- safe_cox_on_scores(pls_scores)
  
  # ---- 5) return_all ----
  if (return_all) {
    K <- n_components
    coef_matrix <- matrix(NA_real_, nrow = K, ncol = K)
    
    #  i 
    for (i in seq_len(K)) {
      scores_i <- pls_scores[, 1:i, drop = FALSE]
      cox_i <- safe_cox_on_scores(scores_i)
      
      if (is.null(cox_i)) {
        #  i
        next
      }
      coefs <- coef(cox_i)
      coef_matrix[seq_along(coefs), i] <- coefs
    }
    
    return(list(
      pls_scores = pls_scores,
      final_cox_model = final_cox_model,
      pls_model = pls_fit,
      X = X,
      coef_matrix = coef_matrix,
      A = A,
      normalize_alpha_terms = normalize_alpha_terms
    ))
  }
  
  final_cox_model
}

# # "NNE"
# val_DRPCAPLS_cv <- function(X, time, status,
#                             k = 5,
#                             ncomp_candidates = 1:5,
#                             alpha_candidates = seq(0, 1, by = 0.1),
#                             auc_time_grid = NULL) {
#   
#   n <- nrow(X)
#   set.seed(123)
#   folds <- sample(rep(1:k, length.out = n))
#   
#   # if (is.null(auc_time_grid)) {
#   #   auc_time_grid <- seq(min(time), max(time), length.out = 50)
#   # }
#   
#   if (is.null(auc_time_grid)) {
#     #  time  10%  90% 10
#     auc_time_grid <- quantile(time, probs = seq(0.1, 0.8, length.out = 15))
#   }
#   
#   results <- data.frame(ncomp = integer(), alpha = numeric(), iAUC = numeric(),
#                         auc_curve = I(list()))
#   candidate_auc_curves <- list()
#   
#   best_iAUC <- -Inf
#   best_ncomp <- NA
#   best_alpha <- NA
#   best_model <- NULL
#   best_A <- NULL
#   best_auc_curve <- NULL
#   
#   for (nc in ncomp_candidates) {
#     print(paste("Evaluating ncomp =", nc))
#     for (a in alpha_candidates) {
#       print(paste("Evaluating alpha =", a))
#       
#       fold_auc_curves <- matrix(NA, nrow = k, ncol = length(auc_time_grid))
#       fold_iAUCs <- numeric(k)
#       
#       for (f in 1:k) {
#         test_idx <- which(folds == f)
#         train_idx <- setdiff(1:n, test_idx)
#         
#         model <- cox_pca_pls_dr(X = X[train_idx, ],
#                                 time = time[train_idx],
#                                 status = status[train_idx],
#                                 n_components = nc,
#                                 alpha = a,
#                                 return_all = TRUE)
#         
#         A <- model$A
#         T_test <- as.data.frame(as.matrix(X[test_idx, ]) %*% A)
#         colnames(T_test) <- paste0("comp_", 1:nc)
#         
#         marker_lp <- predict(model$final_cox_model, newdata = T_test, type = "lp")
#         
#         #  AUC  iAUC
#         # auc_vals <- numeric(length(auc_time_grid))
#         # for (i in seq_along(auc_time_grid)) {
#         #   tp <- auc_time_grid[i]
#         #   roc_obj <- survivalROC(Stime = time[test_idx], status = status[test_idx],
#         #                          marker = marker_lp, predict.time = tp, method = "NNE",
#         #                          span = 0.25)
#         #   auc_vals[i] <- roc_obj$AUC
#         # }
#         
#         auc_vals <- sapply(auc_time_grid, function(tp) {
#           res <- tryCatch({
#             roc_obj <- survivalROC(
#               Stime        = time[test_idx],
#               status       = status[test_idx],
#               marker       = marker_lp,
#               predict.time = tp,
#               method       = "NNE",
#               span         = 0.30
#             )
#             roc_obj$AUC
#           }, error = function(e) {
#             warning(paste("tp =", tp, ":", e$message))
#             return(NA)
#           })
#           return(res)
#         })
#         
#         fold_auc_curves[f, ] <- auc_vals
#         fold_iAUCs[f] <- mean(auc_vals, na.rm = TRUE)
#       }
#       
#       #  AUC  iAUC
#       mean_auc_curve <- colMeans(fold_auc_curves, na.rm = TRUE)
#       mean_iAUC <- mean(fold_iAUCs, na.rm = TRUE)
#       
#       # 
#       results <- rbind(results, data.frame(ncomp = nc, alpha = a, iAUC = mean_iAUC,
#                                            auc_curve = I(list(mean_auc_curve))))
#       key <- paste("ncomp", nc, "alpha", a, sep = "_")
#       candidate_auc_curves[[key]] <- mean_auc_curve
#       
#       if (mean_iAUC > best_iAUC) {
#         best_iAUC <- mean_iAUC
#         best_ncomp <- nc
#         best_alpha <- a
#         best_model <- model$final_cox_model
#         best_A <- A
#         best_auc_curve <- mean_auc_curve
#       }
#     }
#   }
#   
#   return(list(best_model = best_model,
#               best_ncomp = best_ncomp,
#               best_alpha = best_alpha,
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
# val_DRPCAPLS_cv <- function(X, time, status,
#                             k = 5,
#                             ncomp_candidates = 1:5,
#                             alpha_candidates = seq(0, 1, by = 0.1),
#                             auc_time_grid = NULL) {
#   
#   n <- nrow(X)
#   set.seed(123)
#   folds <- sample(rep(1:k, length.out = n))
# 
#   if (is.null(auc_time_grid)) {
#     auc_time_grid <- seq(min(time), max(time), length.out = 50)
#   }
#   
#   results <- data.frame(ncomp = integer(), alpha = numeric(), iAUC = numeric(),
#                         auc_curve = I(list()))
#   candidate_auc_curves <- list()
#   
#   best_iAUC <- -Inf
#   best_ncomp <- NA
#   best_alpha <- NA
#   best_model <- NULL
#   best_A <- NULL
#   best_auc_curve <- NULL
#   
#   for (nc in ncomp_candidates) {
#     print(paste("Evaluating ncomp =", nc))
#     for (a in alpha_candidates) {
#       print(paste("Evaluating alpha =", a))
#       
#       fold_auc_curves <- matrix(NA, nrow = k, ncol = length(auc_time_grid))
#       fold_iAUCs <- numeric(k)
#       
#       for (f in 1:k) {
#         test_idx <- which(folds == f)
#         train_idx <- setdiff(1:n, test_idx)
#         
#         model <- cox_pca_pls_dr(X = X[train_idx, ],
#                                 time = time[train_idx],
#                                 status = status[train_idx],
#                                 n_components = nc,
#                                 alpha = a,
#                                 return_all = TRUE)
#         
#         A <- model$A
#         T_test <- as.data.frame(as.matrix(X[test_idx, ]) %*% A)
#         colnames(T_test) <- paste0("comp_", 1:nc)
#         
#         marker_lp <- predict(model$final_cox_model, newdata = T_test, type = "lp")
#         
#         #  AUC  iAUC
#         auc_vals <- numeric(length(auc_time_grid))
#         for (i in seq_along(auc_time_grid)) {
#           tp <- auc_time_grid[i]
#           roc_obj <- survivalROC(Stime = time[test_idx], status = status[test_idx],
#                                  marker = marker_lp, predict.time = tp, method = "KM")
#           auc_vals[i] <- roc_obj$AUC
#         }
#         
#         fold_auc_curves[f, ] <- auc_vals
#         fold_iAUCs[f] <- mean(auc_vals, na.rm = TRUE)
#       }
#       
#       #  AUC  iAUC
#       mean_auc_curve <- colMeans(fold_auc_curves, na.rm = TRUE)
#       mean_iAUC <- mean(fold_iAUCs, na.rm = TRUE)
#       
#       # 
#       results <- rbind(results, data.frame(ncomp = nc, alpha = a, iAUC = mean_iAUC,
#                                            auc_curve = I(list(mean_auc_curve))))
#       key <- paste("ncomp", nc, "alpha", a, sep = "_")
#       candidate_auc_curves[[key]] <- mean_auc_curve
#       
#       if (mean_iAUC > best_iAUC) {
#         best_iAUC <- mean_iAUC
#         best_ncomp <- nc
#         best_alpha <- a
#         best_model <- model$final_cox_model
#         best_A <- A
#         best_auc_curve <- mean_auc_curve
#       }
#     }
#   }
#   
#   return(list(best_model = best_model,
#               best_ncomp = best_ncomp,
#               best_alpha = best_alpha,
#               best_iAUC = best_iAUC,
#               results = results,
#               best_A = best_A,
#               best_auc_curve = best_auc_curve,
#               candidate_auc_curves = candidate_auc_curves,
#               auc_time_grid = auc_time_grid))
# }



# pre select features 
# val_DRPCAPLS_cv <- function(
#     X, time, status,
#     k = 5,
#     ncomp_candidates = 1:5,
#     alpha_candidates = seq(0, 1, by = 0.1),
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
#         pval <- summary(fit)$coefficients[, "Pr(>|z|)"]
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
#   # ---  CV ---
#   n <- nrow(X)
#   set.seed(123)
#   folds <- sample(rep(1:k, length.out = n))
#   
#   if (is.null(auc_time_grid)) {
#     auc_time_grid <- quantile(time, probs = seq(0.1, 0.8, length.out = 15))
#   }
#   
#   results <- data.frame(
#     ncomp     = integer(),
#     alpha     = numeric(),
#     iAUC      = numeric(),
#     auc_curve = I(list())
#   )
#   candidate_auc_curves <- list()
#   
#   best_iAUC     <- -Inf
#   best_ncomp    <- NA
#   best_alpha    <- NA
#   best_model    <- NULL
#   best_A        <- NULL
#   best_curve    <- NULL
#   
#   for (nc in ncomp_candidates) {
#     cat("Evaluating ncomp =", nc, "...\n")
#     for (a in alpha_candidates) {
#       cat("  Evaluating alpha =", a, "with", auc_method, "...\n")
#       
#       fold_auc_curves <- matrix(NA, nrow = k, ncol = length(auc_time_grid))
#       fold_iAUCs     <- numeric(k)
#       
#       for (f in seq_len(k)) {
#         test_idx  <- which(folds == f)
#         train_idx <- setdiff(seq_len(n), test_idx)
#         
#         model <- cox_pca_pls_dr(
#           X            = X[train_idx, , drop = FALSE],
#           time         = time[train_idx],
#           status       = status[train_idx],
#           n_components = nc,
#           alpha        = a,
#           return_all   = TRUE
#         )
#         
#         A        <- model$A
#         T_test   <- as.data.frame(as.matrix(X[test_idx, , drop = FALSE]) %*% A)
#         colnames(T_test) <- paste0("comp_", seq_len(nc))
#         
#         marker_lp <- predict(model$final_cox_model,
#                              newdata = T_test,
#                              type   = "lp")
#         
#         #   auc_method  AUC   
#         if (auc_method == "NNE") {
#           auc_vals <- sapply(auc_time_grid, function(tp) {
#             tryCatch({
#               roc_obj <- survivalROC::survivalROC(
#                 Stime        = time[test_idx],
#                 status       = status[test_idx],
#                 marker       = marker_lp,
#                 predict.time = tp,
#                 method       = "NNE",
#                 span         = 0.30
#               )
#               roc_obj$AUC
#             }, error = function(e) {
#               warning("tp =", tp, " error: ", e$message)
#               NA
#             })
#           })
#         } else {
#           auc_vals <- sapply(auc_time_grid, function(tp) {
#             roc_obj <- survivalROC::survivalROC(
#               Stime        = time[test_idx],
#               status       = status[test_idx],
#               marker       = marker_lp,
#               predict.time = tp,
#               method       = "KM"
#             )
#             roc_obj$AUC
#           })
#         }
#         
#         fold_auc_curves[f, ] <- auc_vals
#         fold_iAUCs[f]       <- mean(auc_vals, na.rm = TRUE)
#       }
#       
#       mean_curve <- colMeans(fold_auc_curves, na.rm = TRUE)
#       mean_iAUC  <- mean(fold_iAUCs, na.rm = TRUE)
#       
#       results <- rbind(
#         results,
#         data.frame(ncomp = nc, alpha = a, iAUC = mean_iAUC, auc_curve = I(list(mean_curve)))
#       )
#       key <- paste("ncomp", nc, "alpha", a, sep = "_")
#       candidate_auc_curves[[key]] <- mean_curve
#       
#       if (mean_iAUC > best_iAUC) {
#         best_iAUC  <- mean_iAUC
#         best_ncomp <- nc
#         best_alpha <- a
#         best_model <- model$final_cox_model
#         best_A     <- A
#         best_curve <- mean_curve
#       }
#     }
#   }
#   
#   list(
#     best_model         = best_model,
#     best_ncomp         = best_ncomp,
#     best_alpha         = best_alpha,
#     best_iAUC          = best_iAUC,
#     results            = results,
#     best_A             = best_A,
#     best_auc_curve     = best_curve,
#     candidate_auc_curves = candidate_auc_curves,
#     auc_time_grid      = auc_time_grid
#   )
# }


# multi-index
val_DRPCAPLS_cv <- function(
    X, time, status,
    k = 5,
    ncomp_candidates = 1:5,
    alpha_candidates = seq(0, 1, by = 0.1),
    auc_time_grid = NULL,
    use_preselection = FALSE,
    p_thresh = 0.05,
    auc_method = c("IPCW", "NNE", "KM"),
    folds = NULL,
    seed = 123,
    residual_type = "deviance",
    normalize_alpha_terms = TRUE
) {
  auc_method <- match.arg(auc_method)
  n <- nrow(X)
  
  if (is.null(folds)) {
    if (!is.null(seed)) set.seed(seed)
    folds <- sample(rep(1:k, length.out = n))
  } else {
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
  
  # --- Step 1 CV ---
  if (is.null(auc_time_grid)) {
    auc_time_grid <- quantile(time, probs = seq(0.1, 0.8, length.out = 15))
  }
  
  results <- data.frame(
    ncomp     = integer(),
    alpha     = numeric(),
    iAUC      = numeric(),
    auc_curve = I(list())
  )
  candidate_auc_curves <- list()
  
  best_iAUC  <- -Inf
  best_ncomp <- NA
  best_alpha <- NA
  best_model <- NULL
  best_A     <- NULL
  best_curve <- NULL
  
  # --- Step 2 ---
  for (nc in ncomp_candidates) {
    cat("Evaluating ncomp =", nc, "...\n")
    for (a in alpha_candidates) {
      cat("  Evaluating alpha =", a, "with", auc_method, "...\n")
      
      fold_auc_curves <- matrix(NA, nrow = k, ncol = length(auc_time_grid))
      fold_iAUCs     <- numeric(k)
      
      for (f in seq_len(k)) {
        test_idx  <- which(folds == f)
        train_idx <- setdiff(seq_len(n), test_idx)
        
        model <- cox_pca_pls_dr(
          X            = X[train_idx, , drop = FALSE],
          time         = time[train_idx],
          status       = status[train_idx],
          residual_type = residual_type,
          n_components = nc,
          alpha        = a,
          normalize_alpha_terms = normalize_alpha_terms,
          return_all   = TRUE
        )
        
        A        <- model$A
        T_test   <- as.data.frame(as.matrix(X[test_idx, , drop = FALSE]) %*% A)
        
        # ----  ----
        if (!is.null(model$final_cox_model$coefficients)) {
          colnames(T_test) <- names(model$final_cox_model$coefficients)
        } else {
          colnames(T_test) <- paste0("comp_", seq_len(nc))
        }
        
        marker_lp <- tryCatch({
          predict(model$final_cox_model,
                  newdata = T_test,
                  type   = "lp")
        }, error = function(e) {
          warning("Fold ", f, ": prediction failed -> ", e$message)
          rep(NA_real_, length(test_idx))
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
        data.frame(ncomp = nc, alpha = a, iAUC = mean_iAUC, auc_curve = I(list(mean_curve)))
      )
      key <- paste("ncomp", nc, "alpha", a, sep = "_")
      candidate_auc_curves[[key]] <- mean_curve
      
      if (is.finite(mean_iAUC) && mean_iAUC > best_iAUC) {
        best_iAUC  <- mean_iAUC
        best_ncomp <- nc
        best_alpha <- a
        best_model <- model$final_cox_model
        best_A     <- A
        best_curve <- mean_curve
      }
    }
  }
  
  # --- Step 3 ---
  if (!is.finite(best_iAUC) || is.na(best_ncomp) || is.na(best_alpha)) {
    stop("val_DRPCAPLS_cv failed for all parameter combinations.")
  }
  
  final_model <- cox_pca_pls_dr(
    X            = X,
    time         = time,
    status       = status,
    residual_type = residual_type,
    n_components = best_ncomp,
    alpha        = best_alpha,
    normalize_alpha_terms = normalize_alpha_terms,
    return_all   = TRUE
  )
  
  T_all <- as.data.frame(as.matrix(X) %*% final_model$A)
  
  # ----  ----
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
    stop("Final model prediction failed: ", e$message)
  })
  
  # --- Step 4 ---
  cindex_harrell <- survcomp::concordance.index(
    x = marker_all, surv.time = time, surv.event = status, method = "noether"
  )$c.index
  
  cindex_uno <- tryCatch({
    survAUC::UnoC(Surv(time, status),
                  Surv(time, status),
                  marker_all)$C
  }, error = function(e) NA)
  
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
    NA })
  
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
  
  
  # --- Step 5 ---
  list(
    best_model          = final_model$final_cox_model,
    best_ncomp          = best_ncomp,
    best_alpha          = best_alpha,
    best_iAUC           = best_iAUC,
    results             = results,
    best_A              = final_model$A,
    best_auc_curve      = best_curve,
    candidate_auc_curves = candidate_auc_curves,
    auc_time_grid       = auc_time_grid,
    residual_type       = residual_type,
    normalize_alpha_terms = normalize_alpha_terms,
    metrics             = metrics
  )
}






# val_DRPCAPLS_cv <- function(X, time, status,
#                             k = 5,
#                             ncomp_candidates = 1:5,
#                             alpha_candidates = seq(0, 1, by = 0.1),
#                             auc_time_grid = NULL) {
# 
#   n <- nrow(X)
#   set.seed(123)
#   folds <- sample(rep(1:k, length.out = n))
# 
#   # 
#   if (is.null(auc_time_grid)) {
#     auc_time_grid <- seq(min(time), max(time), length.out = 50)
#   }
# 
#   results <- data.frame(ncomp = integer(), alpha = numeric(), iAUC = numeric(),
#                         auc_curve = I(list()))
#   candidate_auc_curves <- list()
# 
#   best_iAUC <- -Inf
#   best_ncomp <- NA
#   best_alpha <- NA
#   best_model <- NULL
#   best_A <- NULL
#   best_auc_curve <- NULL
# 
#   for (nc in ncomp_candidates) {
#     print(paste("Evaluating ncomp =", nc))
#     for (a in alpha_candidates) {
#       print(paste("Evaluating alpha =", a))
#       pooled_time <- c()
#       pooled_status <- c()
#       pooled_marker <- c()
# 
#       # k-fold cross-validation
#       for (f in 1:k) {
# 
#         test_idx <- which(folds == f)
#         train_idx <- setdiff(1:n, test_idx)
# 
#         model <- cox_pca_pls_dr(X = X[train_idx, ],
#                                 time = time[train_idx],
#                                 status = status[train_idx],
#                                 n_components = nc,
#                                 alpha = a,
#                                 return_all = TRUE)
# 
#         A <- model$A
# 
#         T_test <- as.data.frame(as.matrix(X[test_idx, ]) %*% A)
#         colnames(T_test) <- paste0("comp_", 1:nc)
# 
#         marker_lp <- predict(model$final_cox_model, newdata = T_test, type = "lp")
# 
#         pooled_time <- c(pooled_time, time[test_idx])
#         pooled_status <- c(pooled_status, status[test_idx])
#         pooled_marker <- c(pooled_marker, marker_lp)
#       }
# 
#       # AUC
#       auc_vals <- sapply(auc_time_grid, function(tp) {
#         roc_obj <- survivalROC(Stime = pooled_time, status = pooled_status,
#                                marker = pooled_marker, predict.time = tp, method = "KM")
#         return(roc_obj$AUC)
#       })
# 
#       iAUC <- mean(auc_vals, na.rm = TRUE)
# 
#       results <- rbind(results, data.frame(ncomp = nc, alpha = a, iAUC = iAUC,
#                                            auc_curve = I(list(auc_vals))))
#       key <- paste("ncomp", nc, "alpha", a, sep = "_")
#       candidate_auc_curves[[key]] <- auc_vals
# 
#       if (iAUC > best_iAUC) {
#         best_iAUC <- iAUC
#         best_ncomp <- nc
#         best_alpha <- a
#         best_model <- model$final_cox_model
#         best_A <- A
#         best_auc_curve <- auc_vals
#       }
#     }
#   }
# 
#   return(list(best_model = best_model,
#               best_ncomp = best_ncomp,
#               best_alpha = best_alpha,
#               best_iAUC = best_iAUC,
#               results = results,
#               best_A = best_A,
#               best_auc_curve = best_auc_curve,
#               candidate_auc_curves = candidate_auc_curves,
#               auc_time_grid = auc_time_grid))
# }





# val_DRPCAPLS <- function(X_train, time_train, status_train,
#                          X_test, time_test, status_test,
#                          ncomp_candidates = 1:5,
#                          alpha_candidates = seq(0, 1, by = 0.1),
#                          auc_time_grid = NULL) {
# 
#   if (is.null(auc_time_grid)) {
#     auc_time_grid <- seq(min(c(time_train, time_test)), max(c(time_train, time_test)), length.out = 50)
#   }
#   
# 
#   #  auc_curve  AUC  list 
#   results <- data.frame(ncomp = integer(0), alpha = numeric(0), iAUC = numeric(0),
#                         auc_curve = I(list()))
#   candidate_auc_curves <- list()  #  auc_time_grid  AUC 
# 
#   best_iAUC <- -Inf
#   best_ncomp <- NA
#   best_alpha <- NA
#   best_model <- NULL
#   best_A <- NULL
#   best_auc_curve <- NULL
# 
#   #  alpha 
#   for (nc in ncomp_candidates) {
#     for (a in alpha_candidates) {
#       #  cox_pls_dr 
#       model <- cox_pca_pls_dr(X = X_train, time = time_train, status = status_train,
#                           n_components = nc, return_all = TRUE, alpha = a)
#       A <- model$A   # 
# 
#       #  X_test  T_test
#       T_test <- as.data.frame(as.matrix(X_test) %*% A)
#       colnames(T_test) <- paste0("comp_", 1:nc)
# 
#       #  Cox 
#       marker_lp <- predict(model$final_cox_model, newdata = T_test, type = "lp")
# 
#       #  AUC
#       auc_vals <- numeric(length(auc_time_grid))
#       for (i in seq_along(auc_time_grid)) {
#         tp <- auc_time_grid[i]
#         roc_obj <- survivalROC(Stime = time_test, status = status_test,
#                                marker = marker_lp, predict.time = tp, method = "KM")
#         auc_vals[i] <- roc_obj$AUC
#       }
#       #  AUC AUC 
#       iAUC <- mean(auc_vals, na.rm = TRUE)
# 
#       #  auc_vals  auc_curve  list 
#       results <- rbind(results, data.frame(ncomp = nc, alpha = a, iAUC = iAUC,
#                                            auc_curve = I(list(auc_vals))))
#       key <- paste("ncomp", nc, "alpha", a, sep = "_")
#       candidate_auc_curves[[key]] <- auc_vals
# 
#       #  AUC 
#       if (iAUC > best_iAUC) {
#         best_iAUC <- iAUC
#         best_ncomp <- nc
#         best_alpha <- a
#         best_model <- model$final_cox_model
#         best_A <- A
#         best_auc_curve <- auc_vals
#       }
#     }
#   }
# 
#   return(list(best_model = best_model,
#               best_ncomp = best_ncomp,
#               best_alpha = best_alpha,
#               best_iAUC = best_iAUC,
#               results = results,
#               best_A = best_A,
#               best_auc_curve = best_auc_curve,
#               candidate_auc_curves = candidate_auc_curves,
#               auc_time_grid = auc_time_grid))
# }



# select the best alpha
# cox_pca_pls_dr_select_alpha <- function(X, time, status, time2 = NULL,
#                                         residual_type = "deviance",
#                                         scale_time = TRUE,
#                                         n_components = min(7, ncol(X)),
#                                         alpha_vals = seq(0, 1, by = 0.1),
#                                         return_all = TRUE) {
#   best_cindex <- -Inf    #  c-index
#   best_alpha <- NA       #  alpha
#   best_model <- NULL     # 
#   results <- data.frame(alpha = numeric(0), cindex = numeric(0))  #  alpha  c-index
#   
#   #  alpha 
#   for (a in alpha_vals) {
#     # 
#     model <- cox_pca_pls_dr(X = X, time = time, status = status, time2 = time2,
#                             residual_type = residual_type,
#                             scale_time = scale_time,
#                             n_components = n_components,
#                             alpha = a,
#                             return_all = return_all)
#     
#     #  return_all 
#     final_cox_model <- if (return_all) model$final_cox_model else model
#     
#     #  summary()  concordance indexc-index
#     cindex <- summary(final_cox_model)$concordance[1]
#     
#     #  alpha  c-index 
#     results <- rbind(results, data.frame(alpha = a, cindex = cindex))
#     
#     #  c-index 
#     if (cindex > best_cindex) {
#       best_cindex <- cindex
#       best_alpha <- a
#       best_model <- model
#     }
#   }
#   
#   #  alpha c-index  alpha 
#   return(list(best_model = if (return_all) best_model$final_cox_model else best_model,
#               A = best_model$A,
#               best_alpha = best_alpha,
#               best_cindex = best_cindex,
#               results = results))
# }
