

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
  stop("Cannot locate R/auc_utils.R; source it before partial_cox.R.")
}

safe_cox <- function(formula, data) {
  ctrl <- coxph.control(toler.chol = 1e-9)   #  
  tryCatch(
    coxph(formula,
          data    = data,
          control = ctrl,
          singular.ok = TRUE),               #  
    error   = function(e) NULL,
    warning = function(w) {
      if (grepl("coxph.wtest", w$message)) invokeRestart("muffleWarning")
      NULL
    }
  )
}


partial_cox <- function(X, time, status,
                        n_components,
                        var_eps = 1e-8) {
  
  ## -------- Step 0:  --------
  X_means <- colMeans(X)
  V <- sweep(as.matrix(X), 2, X_means, "-")   # V_1j
  n <- nrow(V); p <- ncol(V)
  
  ## --------  --------
  T_matrix   <- matrix(0, n, n_components)
  colnames(T_matrix) <- paste0("T", seq_len(n_components))
  betas_list <- vector("list", n_components)
  w_list     <- vector("list", n_components)
  
  ## -------- Component 1 --------
  beta_1j <- numeric(p)
  for (j in seq_len(p)) {
    beta_1j[j] <- coxph(Surv(time, status) ~ V[, j])$coef
  }
  
  w1 <- apply(V, 2, var)
  w1 <- w1 / sum(w1)
  T_matrix[, 1] <- V %*% (w1 * beta_1j)
  
  betas_list[[1]] <- beta_1j
  w_list[[1]]     <- w1
  
  safe_var <- function(x) {                        
    v <- var(x, na.rm = TRUE)
    if (is.na(v)) 0 else v
  }
  
  ## -------- Components 2  K --------
  if (n_components >= 2L) for (k in 2:n_components) {
    
    ## (a)  ( T_{k-1}  )
    T_prev     <- T_matrix[, k - 1]
    denom      <- sum(T_prev^2)
    V_residual <- V                                     # 
    for (j in seq_len(p)) {
      if (denom == 0) {                             ## <<<  0/0
        V_residual[, j] <- V[, j]
      } else {
        proj_coef <- sum(V[, j] * T_prev) / denom
        if (is.nan(proj_coef)) proj_coef <- 0       ## <<< 
        V_residual[, j] <- V[, j] - proj_coef * T_prev
      }
    }
    
    ## (b)  CoxT1T_{k-1} + Vj
    beta_kj <- numeric(p)
    prev_T  <- as.data.frame(T_matrix[, 1:(k - 1), drop = FALSE])
    names(prev_T) <- paste0("T", seq_len(k - 1))
    
    for (j in seq_len(p)) {
      vj_var <- safe_var(V_residual[, j])           ## <<<  safe_var
      if (vj_var < var_eps) {                       ## <<<  0 
        beta_kj[j] <- 0
        next
      }
      
      df  <- cbind(prev_T,
                   Vj     = V_residual[, j],
                   time   = time,
                   status = status)
      fit <- safe_cox(Surv(time, status) ~ ., data = df)
      beta_kj[j] <- if (is.null(fit)) 0 else coef(fit)["Vj"]
    }
    
    ## (c) 
    wk <- apply(V_residual, 2, safe_var)
    wk <- wk / sum(wk)
    
    ## (d)  T_k
    T_matrix[, k] <- V_residual %*% (wk * beta_kj)
    
    betas_list[[k]] <- beta_kj
    w_list[[k]]     <- wk
    
    ## (e)  V
    V <- V_residual
  }
  
  ## --------  Cox --------
  final_df    <- data.frame(time = time, status = status, T_matrix)
  final_model <- coxph(Surv(time, status) ~ ., data = final_df, x = TRUE,
                       y = TRUE)
  
  ## --------  --------
  list(
    components  = T_matrix,
    final_model = final_model,
    betas_list  = betas_list,
    w_list      = w_list,
    X_means     = X_means        # 
  )
}


partial_cox <- function(X, time, status, n_components) {
  # X:          n  p 
  # time:       
  # status:      (0=,1=)
  # n_components:  K
  
  n <- nrow(X)
  p <- ncol(X)
  
  V <- as.matrix(X)
  
  ##     
  T_matrix   <- matrix(0, n, n_components)
  colnames(T_matrix) <- paste0("T", seq_len(n_components))
  betas_list <- vector("list", n_components)  #  _{kj}
  w_list     <- vector("list", n_components)  #  w_{kj}
  
  ##  Component 1   
  # (1)  Cox 
  beta_1j <- numeric(p)
  for (j in seq_len(p)) {
    beta_1j[j] <- coxph(Surv(time, status) ~ V[, j])$coef
  }
  # (2) 
  w1 <- apply(V, 2, var)
  w1 <- w1 / sum(w1)
  # (3)  T1
  T_matrix[, 1] <- V %*% (w1 * beta_1j)
  betas_list[[1]] <- beta_1j
  w_list[[1]]     <- w1
  
  ##  Components 2  K   
  if (n_components >= 2L) for (k in 2:n_components) {
    # (a)  V_residual
    T_prev <- T_matrix[, k - 1]
    V_residual <- V  # 
    for (j in seq_len(p)) {
      proj_coef <- sum(V[, j] * T_prev) / sum(T_prev^2)
      V_residual[, j] <- V[, j] - proj_coef * T_prev
    }
    
    # (b)  Cox T1T_{k-1} +  Vj
    beta_kj <- numeric(p)
    prev_T  <- as.data.frame(T_matrix[, 1:(k - 1), drop = FALSE])
    names(prev_T) <- paste0("T", 1:(k - 1))
    for (j in seq_len(p)) {
      df <- cbind(prev_T,
                  Vj    = V_residual[, j],
                  time  = time,
                  status = status)
      fit_joint <- coxph(Surv(time, status) ~ ., data = df)
      beta_kj[j] <- coef(fit_joint)["Vj"]
    }
    
    # (c) 
    wk <- apply(V_residual, 2, var)
    wk <- wk / sum(wk)
    # (d)  T_k
    T_matrix[, k] <- V_residual %*% (wk * beta_kj)
    
    # 
    betas_list[[k]] <- beta_kj
    w_list[[k]]     <- wk
    
    # (e)  V
    V <- V_residual
  }
  
  ##   Cox    
  final_df   <- data.frame(time = time,
                           status = status,
                           T_matrix)
  # final_model <- coxph(Surv(time, status) ~ ., data = final_df)
  
  final_model <- tryCatch(
    {
      survival::coxph(
        survival::Surv(time, status) ~ .,
        data    = final_df,
        x       = TRUE,
        y       = TRUE
      )
    },
    error = function(e) {
      warning("Final Cox model fitting failed: ", e$message)
      NULL
    }
  )
  
  
  ##     
  list(
    components  = T_matrix,
    final_model = final_model,
    betas_list  = betas_list,
    w_list      = w_list
  )
}



partial_cox <- function(X, time, status, n_components) {
  n <- nrow(X)
  p <- ncol(X)

  #  ()
  V <- as.matrix(X)

  # 
  T_matrix <- matrix(0, n, n_components)
  colnames(T_matrix) <- paste0("T", 1:n_components)

  # 
  betas_list <- vector("list", n_components)  # pbeta
  w_list     <- vector("list", n_components)  # w

  # 1
  beta_1j <- numeric(p)
  for (j in 1:p) {
    fit_j <- coxph(Surv(time, status) ~ V[, j])
    beta_1j[j] <- fit_j$coef
  }

  w1 <- abs(beta_1j) / sum(abs(beta_1j))   # 
  T1 <- V %*% (w1 * beta_1j)
  T_matrix[, 1] <- T1

  betas_list[[1]] <- beta_1j
  w_list[[1]]     <- w1

  # 
  if (n_components >= 2L) for (k in 2:n_components) {
    beta_kj <- numeric(p)
    V_residual <- matrix(0, n, p)

    # 
    for (j in 1:p) {
      proj_coef <- sum(V[, j] * T_matrix[, k - 1]) / sum(T_matrix[, k - 1]^2)
      V_residual[, j] <- V[, j] - proj_coef * T_matrix[, k - 1]
    }

    # Cox
    for (j in 1:p) {
      df <- data.frame(T_prev = T_matrix[, k - 1], Vj = V_residual[, j], time, status)
      fit_joint <- coxph(Surv(time, status) ~ T_prev + Vj, data = df)
      beta_kj[j] <- coef(fit_joint)["Vj"]
    }

    wk <- abs(beta_kj) / sum(abs(beta_kj))
    T_k <- V_residual %*% (wk * beta_kj)
    T_matrix[, k] <- T_k

    # 
    betas_list[[k]] <- beta_kj
    w_list[[k]]     <- wk

    #  V
    V <- V_residual
  }

  final_df <- data.frame(
    time   = time,
    status = status,
    T_matrix
  )


  final_model <- coxph(Surv(time, status) ~ ., data = final_df)

  # 
  list(
    components = T_matrix,
    final_model = final_model,
    betas_list = betas_list,
    w_list = w_list,
    X_means = colMeans(X)
  )
}





transform_partial_cox <- function(X_new, betas_list, w_list) {
  #  X_new 
  V_test <- as.matrix(X_new)
  
  n <- nrow(X_new)
  p <- ncol(X_new)
  ncomp <- length(betas_list)
  
  T_test <- matrix(0, n, ncomp)
  colnames(T_test) <- paste0("T", 1:ncomp)
  
  # 
  beta_1j <- betas_list[[1]]
  w1      <- w_list[[1]]
  T1_test <- V_test %*% (w1 * beta_1j)
  T_test[, 1] <- T1_test
  
  # 
  if (ncomp >= 2L) for (k in 2:ncomp) {
    beta_kj <- betas_list[[k]]
    wk      <- w_list[[k]]
    
    #  T_{k-1} 
    V_residual <- matrix(0, n, p)
    for (j in 1:p) {
      proj_coef <- sum(V_test[, j] * T_test[, k - 1]) / sum(T_test[, k - 1]^2)
      V_residual[, j] <- V_test[, j] - proj_coef * T_test[, k - 1]
    }
    
    T_k_test <- V_residual %*% (wk * beta_kj)
    T_test[, k] <- T_k_test
    
    #  V_test 
    V_test <- V_residual
  }
  
  T_test
}



# pre select features
# val_partial_cox_cv <- function(
#     X, time, status,
#     k = 5,
#     ncomp_candidates = 1:5,
#     auc_time_grid = NULL,
#     use_preselection = FALSE,
#     p_thresh = 0.05,
#     auc_method = c("NNE", "KM"),
# ) {
#   #  auc_method
#   auc_method <- match.arg(auc_method)
#   
#   
#   # ---------------------------
#   # Step 0 Coxp < p_thresh
#   # ---------------------------
#   if (use_preselection) {
#     p_total <- ncol(X)
#     sig_idx <- integer(0)
#     for (j in seq_len(p_total)) {
#       fit <- tryCatch({
#         survival::coxph(survival::Surv(time, status) ~ X[, j])
#       }, error = function(e) NULL)
#       if (!is.null(fit)) {
#         pval <- summary(fit)$coefficients[,"Pr(>|z|)"]
#         if (!is.na(pval) && pval < p_thresh) {
#           sig_idx <- c(sig_idx, j)
#         }
#       }
#     }
#     if (length(sig_idx) == 0) {
#       stop("No variables pass univariate Cox at p <", p_thresh)
#     }
#     X <- X[, sig_idx, drop = FALSE]
#     message(length(sig_idx),
#             " predictors retained after preselection (p < ", p_thresh, ").")
#   }
#   
#   # ---------------------------
#   # Step 1 CV 
#   # ---------------------------
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
#     iAUC      = numeric(),
#     auc_curve = I(list())
#   )
#   candidate_auc_curves <- list()
#   
#   best_iAUC       <- -Inf
#   best_ncomp      <- NA
#   best_model      <- NULL
#   best_betas_list <- NULL
#   best_w_list     <- NULL
#   best_curve      <- NULL
#   
#   # ---------------------------
#   # Step 2
#   # ---------------------------
#   for (nc in ncomp_candidates) {
#     cat("Evaluating ncomp =", nc, "with", auc_method, "...\n")
#     fold_auc_curves <- matrix(NA, nrow = k, ncol = length(auc_time_grid))
#     fold_iAUCs     <- numeric(k)
#     
#     for (fold in seq_len(k)) {
#       test_idx  <- which(folds == fold)
#       train_idx <- setdiff(seq_len(n), test_idx)
#       
#       model <- partial_cox(
#         X           = X[train_idx, , drop = FALSE],
#         time        = time[train_idx],
#         status      = status[train_idx],
#         n_components= nc
#       )
#       
#       betas_list <- model$betas_list
#       w_list     <- model$w_list
#       
#       T_test <- transform_partial_cox(
#         X_new      = X[test_idx, , drop = FALSE],
#         betas_list = betas_list,
#         w_list     = w_list
#       )
#       T_test <- as.data.frame(T_test)
#       colnames(T_test) <- paste0("T", seq_len(ncol(T_test)))
#       
#       marker_lp <- predict(model$final_model,
#                            newdata = T_test,
#                            type   = "lp")
#       
#       #   auc_method   
#       if (auc_method == "NNE") {
#         auc_vals <- sapply(auc_time_grid, function(tp) {
#           tryCatch({
#             roc <- survivalROC::survivalROC(
#               Stime        = time[test_idx],
#               status       = status[test_idx],
#               marker       = marker_lp,
#               predict.time = tp,
#               method       = "NNE",
#               span         = 0.30
#             )
#             roc$AUC
#           }, error = function(e) {
#             warning("tp =", tp, " error: ", e$message)
#             NA
#           })
#         })
#       } else {
#         auc_vals <- sapply(auc_time_grid, function(tp) {
#           roc <- survivalROC::survivalROC(
#             Stime        = time[test_idx],
#             status       = status[test_idx],
#             marker       = marker_lp,
#             predict.time = tp,
#             method       = "KM"
#           )
#           roc$AUC
#         })
#       }
#       
#       fold_auc_curves[fold, ] <- auc_vals
#       fold_iAUCs[fold]       <- mean(auc_vals, na.rm = TRUE)
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
#       best_iAUC       <- mean_iAUC
#       best_ncomp      <- nc
#       best_model      <- model$final_model
#       best_betas_list <- betas_list
#       best_w_list     <- w_list
#       best_curve      <- mean_curve
#     }
#   }
#   
#   # ---------------------------
#   # 
#   # ---------------------------
#   list(
#     best_model         = best_model,
#     best_ncomp         = best_ncomp,
#     best_iAUC          = best_iAUC,
#     results            = results,
#     candidate_auc_curves = candidate_auc_curves,
#     auc_time_grid      = auc_time_grid,
#     best_betas_list    = best_betas_list,
#     best_w_list        = best_w_list,
#     best_auc_curve     = best_curve
#   )
# }


# multi-index
val_partial_cox_cv <- function(
    X, time, status,
    k = 5,
    ncomp_candidates = 1:5,
    auc_time_grid = NULL,
    use_preselection = FALSE,
    p_thresh = 0.05,
    auc_method = c("IPCW", "NNE", "KM"),
    folds = NULL,
    seed = NULL
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
  
  # ---------------------------
  # Step 0 Cox
  # ---------------------------
  if (use_preselection) {
    p_total <- ncol(X)
    sig_idx <- integer(0)
    for (j in seq_len(p_total)) {
      fit <- tryCatch({
        survival::coxph(survival::Surv(time, status) ~ X[, j])
      }, error = function(e) NULL)
      if (!is.null(fit)) {
        pval <- summary(fit)$coefficients[, "Pr(>|z|)"]
        if (!is.na(pval) && pval < p_thresh) {
          sig_idx <- c(sig_idx, j)
        }
      }
    }
    if (length(sig_idx) == 0) {
      stop("No variables pass univariate Cox at p <", p_thresh)
    }
    X <- X[, sig_idx, drop = FALSE]
    message(length(sig_idx),
            " predictors retained after preselection (p < ", p_thresh, ").")
  }
  
  # ---------------------------
  # Step 1CV 
  # ---------------------------
  n <- nrow(X)
  
  if (is.null(auc_time_grid)) {
    auc_time_grid <- quantile(time, probs = seq(0.1, 0.8, length.out = 15))
  }
  
  results <- data.frame(
    ncomp     = integer(),
    iAUC      = numeric(),
    auc_curve = I(list())
  )
  candidate_auc_curves <- list()
  
  best_iAUC       <- -Inf
  best_ncomp      <- NA
  best_model      <- NULL
  best_betas_list <- NULL
  best_w_list     <- NULL
  best_curve      <- NULL
  
  # ---------------------------
  # Step 2
  # ---------------------------
  for (nc in ncomp_candidates) {
    cat("Evaluating ncomp =", nc, "with", auc_method, "...\n")
    fold_auc_curves <- matrix(NA, nrow = k, ncol = length(auc_time_grid))
    fold_iAUCs     <- numeric(k)
    
    for (fold in seq_len(k)) {
      test_idx  <- which(folds == fold)
      train_idx <- setdiff(seq_len(n), test_idx)
      
      model <- partial_cox(
        X            = X[train_idx, , drop = FALSE],
        time         = time[train_idx],
        status       = status[train_idx],
        n_components = nc
      )
      
      betas_list <- model$betas_list
      w_list     <- model$w_list
      
      T_test <- transform_partial_cox(
        X_new      = X[test_idx, , drop = FALSE],
        betas_list = betas_list,
        w_list     = w_list
      )
      T_test <- as.data.frame(T_test)
      colnames(T_test) <- paste0("T", seq_len(ncol(T_test)))
      
      marker_lp <- predict(model$final_model,
                           newdata = T_test,
                           type   = "lp")
      
      auc_vals <- .coxpcapls_auc_curve(
        time = time[test_idx],
        status = status[test_idx],
        marker = marker_lp,
        auc_time_grid = auc_time_grid,
        auc_method = auc_method
      )
      
      fold_auc_curves[fold, ] <- auc_vals
      fold_iAUCs[fold]       <- mean(auc_vals, na.rm = TRUE)
    }
    
    mean_curve <- colMeans(fold_auc_curves, na.rm = TRUE)
    mean_iAUC  <- mean(fold_iAUCs, na.rm = TRUE)
    
    results <- rbind(
      results,
      data.frame(ncomp = nc, iAUC = mean_iAUC, auc_curve = I(list(mean_curve)))
    )
    candidate_auc_curves[[as.character(nc)]] <- mean_curve
    
    if (!is.na(mean_iAUC) && mean_iAUC > best_iAUC) {
      best_iAUC       <- mean_iAUC
      best_ncomp      <- nc
      best_model      <- model$final_model
      best_betas_list <- betas_list
      best_w_list     <- w_list
      best_curve      <- mean_curve
    }
  }
  
  # ---------------------------
  # Step 3
  # ---------------------------
  final_model <- partial_cox(
    X            = X,
    time         = time,
    status       = status,
    n_components = best_ncomp
  )
  
  T_all <- transform_partial_cox(
    X_new      = X,
    betas_list = final_model$betas_list,
    w_list     = final_model$w_list
  )
  T_all <- as.data.frame(T_all)
  colnames(T_all) <- paste0("T", seq_len(ncol(T_all)))
  
  marker_all <- predict(final_model$final_model,
                        newdata = T_all,
                        type = "lp")
  
  # ---------------------------
  # Step 4C-index / Brier
  # ---------------------------
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
  
  
  # ---------------------------
  # Step 5
  # ---------------------------
  list(
    best_model          = final_model$final_model,
    best_ncomp          = best_ncomp,
    best_iAUC           = best_iAUC,
    results             = results,
    candidate_auc_curves = candidate_auc_curves,
    auc_time_grid       = auc_time_grid,
    best_betas_list     = final_model$betas_list,
    best_w_list         = final_model$w_list,
    best_auc_curve      = best_curve,
    metrics             = metrics
  )
}






# PC-PCR
val_prepca_partial_cox_cv <- function(
    X, time, status,
    k                  = 5,
    ncomp_candidates   = 1:5,
    auc_time_grid      = NULL,
    auc_method         = c("IPCW", "NNE", "KM"),
    use_pcr            = FALSE,
    use_preselection   = FALSE,   #  Cox 
    p_thresh           = 0.05     #  P 
) {
  auc_method <- match.arg(auc_method)
  
  #  0.  Cox  
  if (use_preselection) {
    p_total <- ncol(X)
    sig_idx <- integer(0)
    for (j in seq_len(p_total)) {
      fit <- tryCatch({
        survival::coxph(survival::Surv(time, status) ~ X[, j])
      }, error = function(e) NULL)
      if (!is.null(fit)) {
        pval <- summary(fit)$coefficients[, "Pr(>|z|)"]
        if (!is.na(pval) && pval < p_thresh) {
          sig_idx <- c(sig_idx, j)
        }
      }
    }
    if (length(sig_idx) == 0) {
      stop("No variables pass univariate Cox at p <", p_thresh)
    }
    X <- X[, sig_idx, drop = FALSE]
    message(length(sig_idx),
            " predictors retained after preselection (p < ", p_thresh, ").")
  }
  
  #  1.    
  n     <- nrow(X)
  set.seed(123)
  folds <- sample(rep(1:k, length.out = n))
  
  #  2.    
  if (is.null(auc_time_grid)) {
    auc_time_grid <- quantile(time, probs = seq(0.1, 0.8, length.out = 15))
  }
  
  #  3.    
  results               <- data.frame(ncomp=integer(), iAUC=numeric(), auc_curve=I(list()))
  candidate_auc_curves  <- list()
  best_iAUC             <- -Inf
  best_ncomp            <- NA
  best_model            <- NULL
  best_betas_list       <- NULL
  best_w_list           <- NULL
  best_auc_curve        <- NULL
  best_pca_obj          <- NULL
  
  #  4.    
  for (nc in ncomp_candidates) {
    cat(">>> ncomp =", nc, "\n")
    fold_iAUCs      <- numeric(k)
    fold_auc_curves <- matrix(NA, k, length(auc_time_grid))
    
    for (fold in seq_len(k)) {
      tr <- which(folds != fold)
      te <- which(folds == fold)
      
      # 4.1 PCA 
      if (use_pcr) {
        pca_obj     <- prcomp(X[tr,,drop=FALSE], center=TRUE, scale.=FALSE)
        scores_tr   <- pca_obj$x
        X_te_center <- sweep(X[te,,drop=FALSE], 2, pca_obj$center, "-")
        scores_te   <- as.matrix(X_te_center) %*% pca_obj$rotation
      } else {
        scores_tr <- X[tr,,drop=FALSE]
        scores_te <- X[te,,drop=FALSE]
      }
      
      # 4.2  Partial Cox
      model <- partial_cox(
        X            = scores_tr,
        time         = time[tr],
        status       = status[tr],
        n_components = nc
      )
      
      # 4.3  T_test
      T_test <- transform_partial_cox(
        X_new      = scores_te,
        betas_list = model$betas_list,
        w_list     = model$w_list
      )
      
      # 4.4  lp_te
      lp_te <- predict(
        model$final_model,
        newdata = as.data.frame(T_test),
        type    = "lp"
      )
      
      # 4.5  AUC
      time_te   <- time[te]
      status_te <- status[te]
      aucs <- .coxpcapls_auc_curve(
        time = time_te,
        status = status_te,
        marker = lp_te,
        auc_time_grid = auc_time_grid,
        auc_method = auc_method
      )
      
      fold_auc_curves[fold, ] <- aucs
      fold_iAUCs[fold]       <- mean(aucs, na.rm=TRUE)
    }
    
    #  5.  ncomp   
    mean_curve <- colMeans(fold_auc_curves, na.rm=TRUE)
    mean_iAUC  <- mean(fold_iAUCs, na.rm=TRUE)
    
    results <- rbind(
      results,
      data.frame(ncomp=nc, iAUC=mean_iAUC, auc_curve=I(list(mean_curve)))
    )
    candidate_auc_curves[[as.character(nc)]] <- mean_curve
    
    if (mean_iAUC > best_iAUC) {
      best_iAUC       <- mean_iAUC
      best_ncomp      <- nc
      best_model      <- model$final_model
      best_auc_curve  <- mean_curve
      best_pca_obj    <- if (use_pcr) pca_obj else NULL
      best_betas_list <- model$betas_list
      best_w_list     <- model$w_list
    }
  }
  
  #     
  list(
    best_model           = best_model,
    best_ncomp           = best_ncomp,
    best_iAUC            = best_iAUC,
    results              = results,
    candidate_auc_curves = candidate_auc_curves,
    auc_time_grid        = auc_time_grid,
    best_betas_list      = best_betas_list,
    best_w_list          = best_w_list,
    best_auc_curve       = best_auc_curve,
    best_pca_obj         = best_pca_obj
  )
}


# not include preselection PC - PCR
val_prepca_partial_cox_cv <- function(
    X, time, status,
    k                 = 5,
    ncomp_candidates  = 1:5,
    auc_time_grid     = NULL,
    auc_method        = c("IPCW", "NNE", "KM"),
    use_pcr           = FALSE
) {
  auc_method <- match.arg(auc_method)
  n          <- nrow(X)
  set.seed(123)
  folds      <- sample(rep(1:k, length.out = n))
  
  if (is.null(auc_time_grid)) {
    auc_time_grid <- quantile(time, probs = seq(0.1, 0.8, length.out = 15))
  }
  
  results            <- data.frame(ncomp=integer(), iAUC=numeric(), auc_curve=I(list()))
  candidate_auc_curves <- list()
  best_iAUC          <- -Inf
  best_ncomp         <- NA
  best_model         <- NULL
  best_betas_list    <- NULL
  best_w_list        <- NULL
  best_auc_curve     <- NULL
  best_pca_obj       <- NULL
  
  for (nc in ncomp_candidates) {
    cat(">>> ncomp =", nc, "\n")
    fold_iAUCs      <- numeric(k)
    fold_auc_curves <- matrix(NA, k, length(auc_time_grid))
    
    for (fold in seq_len(k)) {
      tr <- which(folds != fold)
      te <- which(folds == fold)
      
      # 1.PCA 
      if (use_pcr) {
        pca_obj     <- prcomp(X[tr,,drop=FALSE], center=TRUE, scale.=FALSE)
        scores_tr   <- pca_obj$x
        X_te_center <- sweep(X[te,,drop=FALSE], 2, pca_obj$center, "-")
        scores_te   <- as.matrix(X_te_center) %*% pca_obj$rotation
      } else {
        scores_tr <- X[tr,,drop=FALSE]
        scores_te <- X[te,,drop=FALSE]
      }
      
      # 2.  Partial Cox
      model <- partial_cox(
        X            = scores_tr,
        time         = time[tr],
        status       = status[tr],
        n_components = nc
      )
      
      # 3.  T_test
      T_test <- transform_partial_cox(
        X_new      = scores_te,
        betas_list = model$betas_list,
        w_list     = model$w_list
      )
      
      # 4. 
      lp_te <- predict(
        model$final_model,
        newdata = as.data.frame(T_test),
        type    = "lp"
      )
     
      time_te   <- time[te]
      status_te <- status[te]
      #valid     <- !is.na(time_te) & !is.na(status_te) & !is.na(lp_te)
      
      aucs <- .coxpcapls_auc_curve(
        time = time_te,
        status = status_te,
        marker = lp_te,
        auc_time_grid = auc_time_grid,
        auc_method = auc_method
      )
      
      fold_auc_curves[fold, ] <- aucs
      fold_iAUCs[fold]       <- mean(aucs, na.rm=TRUE)

    }
    
    #  ncomp
    mean_curve <- colMeans(fold_auc_curves, na.rm=TRUE)
    mean_iAUC  <- mean(fold_iAUCs, na.rm=TRUE)
    
    results <- rbind(
      results,
      data.frame(ncomp=nc, iAUC=mean_iAUC, auc_curve=I(list(mean_curve)))
    )
    candidate_auc_curves[[as.character(nc)]] <- mean_curve
    
    if (mean_iAUC > best_iAUC) {
      best_iAUC       <- mean_iAUC
      best_ncomp      <- nc
      best_model      <- model$final_model
      best_auc_curve  <- mean_curve
      best_pca_obj    <- if (use_pcr) pca_obj else NULL
      best_betas_list <- model$betas_list
      best_w_list     <- model$w_list
    }
  }
  
  list(
    best_model           = best_model,
    best_ncomp           = best_ncomp,
    best_iAUC            = best_iAUC,
    results              = results,
    candidate_auc_curves = candidate_auc_curves,
    auc_time_grid        = auc_time_grid,
    best_betas_list      = best_betas_list,
    best_w_list          = best_w_list,
    best_auc_curve       = best_auc_curve,
    best_pca_obj         = best_pca_obj
  )
}







# PCAPLS-PCR
val_pcapls_partial_cox_cv <- function(
    X, time, status,
    k                 = 5,
    ncomp_candidates  = 1:5,
    alpha_candidates  = seq(0, 1, by = 0.1),
    auc_time_grid     = NULL,
    use_preselection   = FALSE,   #  Cox 
    p_thresh           = 0.05,     #  Cox  P 
    auc_method        = c("IPCW", "NNE", "KM")
) {
  auc_method <- match.arg(auc_method)
  
  #  0.  Cox  
  if (use_preselection) {
    p_total <- ncol(X)
    sig_idx <- integer(0)
    for (j in seq_len(p_total)) {
      fit <- safe_cox(Surv(time, status) ~ X[, j],
                      data = data.frame(time, status, x = X[, j]))
      if (!is.null(fit)) {
        pval <- summary(fit)$coefficients[, "Pr(>|z|)"]
        if (!is.na(pval) && pval < p_thresh) {
          sig_idx <- c(sig_idx, j)
        }
      }
    }
    if (length(sig_idx) == 0) {
      stop("No variables pass univariate Cox at p <", p_thresh)
    }
    X <- X[, sig_idx, drop = FALSE]
    message(length(sig_idx),
            " predictors retained after preselection (p < ", p_thresh, ").")
  }
  
  n          <- nrow(X)
  set.seed(123)
  folds      <- sample(rep(1:k, length.out = n))
  
  if (is.null(auc_time_grid)) {
    auc_time_grid <- quantile(time, probs = seq(0.1, 0.8, length.out = 15))
  }
  
  results             <- data.frame(ncomp=integer(), alpha=numeric(),
                                    iAUC=numeric(), auc_curve=I(list()))
  candidate_auc_curves <- list()
  best_iAUC           <- -Inf
  best_ncomp          <- NA
  best_alpha          <- NA
  best_model          <- NULL
  best_betas_list     <- NULL
  best_w_list         <- NULL
  best_auc_curve      <- NULL
  best_pls_res        <- NULL
  
  for (ni in seq_along(ncomp_candidates)) {
    nc <- ncomp_candidates[ni]
    
    for (ai in seq_along(alpha_candidates)) {
      a <- alpha_candidates[ai]
      cat(sprintf("ncomp=%d, alpha=%.1f\n", nc, a))
      
      fold_iAUCs  <- numeric(k)
      fold_curves <- matrix(NA, k, length(auc_time_grid))
      
      for (fold in seq_len(k)) {
        tr <- which(folds != fold)
        te <- which(folds == fold)
        
        # 1.  deviance residuals + 
        null_mod <- coxph(Surv(time[tr], status[tr]) ~ 1)
        Y_tr     <- residuals(null_mod, type = "deviance")
        Xc_tr    <- scale(X[tr,,drop=FALSE], center=TRUE, scale=FALSE)
        cen      <- attr(Xc_tr, "scaled:center")
        
        # 2. PCAPLS  nc 
        pls_res <- acc_pca_pls_with_alpha(
          X     = Xc_tr,
          Y     = Y_tr,
          alpha = a,
          ncomp = min(nrow(Xc_tr), ncol(Xc_tr))
        )
        #  pls_res
        if (is.null(best_pls_res)) best_pls_res <- pls_res
        
        X_tr_mod <- as.matrix(Xc_tr) %*% pls_res$Projection
        Xc_te    <- sweep(X[te,,drop=FALSE], 2, cen, "-")
        X_te_mod <- as.matrix(Xc_te) %*% pls_res$Projection
        
        # 3.  Partial Cox
        model <- partial_cox(
          X            = X_tr_mod,
          time         = time[tr],
          status       = status[tr],
          n_components = nc
        )
        
        # 4.  &  AUC 
        T_test <- transform_partial_cox(
          X_new      = X_te_mod,
          betas_list = model$betas_list,
          w_list     = model$w_list
        )
        lp_te <- predict(model$final_model, newdata=as.data.frame(T_test), type="lp")
        
        time_te   <- time[te]; status_te <- status[te]
        keep      <- !is.na(time_te)&!is.na(status_te)&!is.na(lp_te)
        time_te   <- time_te[keep]; status_te <- status_te[keep]
        marker    <- lp_te[keep]
        
        aucs <- .coxpcapls_auc_curve(
          time = time_te,
          status = status_te,
          marker = marker,
          auc_time_grid = auc_time_grid,
          auc_method = auc_method
        )
        
        fold_curves[fold, ] <- aucs
        fold_iAUCs[fold]    <- mean(aucs, na.rm=TRUE)
      }
      
      mean_curve <- colMeans(fold_curves, na.rm=TRUE)
      mean_iAUC  <- mean(fold_iAUCs, na.rm=TRUE)
      
      #  results
      results <- rbind(
        results,
        data.frame(
          ncomp     = nc,
          alpha     = a,
          iAUC      = mean_iAUC,
          auc_curve = I(list(mean_curve))
        )
      )
      key <- sprintf("%d_%.1f", nc, a)
      candidate_auc_curves[[key]] <- mean_curve
      
      # 
      if (mean_iAUC > best_iAUC) {
        best_iAUC       <- mean_iAUC
        best_ncomp      <- nc
        best_alpha      <- a
        best_model      <- model$final_model
        best_auc_curve  <- mean_curve
        best_betas_list <- model$betas_list
        best_w_list     <- model$w_list
        best_pls_res    <- pls_res
      }
    }
  }
  
  list(
    best_model           = best_model,
    best_ncomp           = best_ncomp,
    best_alpha           = best_alpha,
    best_iAUC            = best_iAUC,
    results              = results,
    candidate_auc_curves = candidate_auc_curves,
    auc_time_grid        = auc_time_grid,
    best_betas_list      = best_betas_list,
    best_w_list          = best_w_list,
    best_auc_curve       = best_auc_curve,
    best_pls_res         = best_pls_res
  )
}




# not include preselection PCAPLS-PCR
val_pcapls_partial_cox_cv <- function(
    X, time, status,
    k                 = 5,
    ncomp_candidates  = 1:5,
    alpha_candidates  = seq(0, 1, by = 0.1),
    auc_time_grid     = NULL,
    auc_method        = c("IPCW", "NNE", "KM")
) {
  auc_method <- match.arg(auc_method)
  n          <- nrow(X)
  set.seed(123)
  folds      <- sample(rep(1:k, length.out = n))
  
  if (is.null(auc_time_grid)) {
    auc_time_grid <- quantile(time, probs = seq(0.1, 0.8, length.out = 15))
  }
  
  results             <- data.frame(ncomp=integer(), alpha=numeric(),
                                    iAUC=numeric(), auc_curve=I(list()))
  candidate_auc_curves <- list()
  best_iAUC           <- -Inf
  best_ncomp          <- NA
  best_alpha          <- NA
  best_model          <- NULL
  best_betas_list     <- NULL
  best_w_list         <- NULL
  best_auc_curve      <- NULL
  best_pls_res        <- NULL
  
  for (ni in seq_along(ncomp_candidates)) {
    nc <- ncomp_candidates[ni]
    
    for (ai in seq_along(alpha_candidates)) {
      a <- alpha_candidates[ai]
      cat(sprintf("ncomp=%d, alpha=%.1f\n", nc, a))
      
      fold_iAUCs  <- numeric(k)
      fold_curves <- matrix(NA, k, length(auc_time_grid))
      
      for (fold in seq_len(k)) {
        tr <- which(folds != fold)
        te <- which(folds == fold)
        
        # 1.  deviance residuals + 
        null_mod <- coxph(Surv(time[tr], status[tr]) ~ 1)
        Y_tr     <- residuals(null_mod, type = "deviance")
        Xc_tr    <- scale(X[tr,,drop=FALSE], center=TRUE, scale=FALSE)
        cen      <- attr(Xc_tr, "scaled:center")
        
        # 2. PCAPLS  nc 
        pls_res <- acc_pca_pls_with_alpha(
          X     = Xc_tr,
          Y     = Y_tr,
          alpha = a,
          ncomp = min(nrow(Xc_tr), ncol(Xc_tr))
        )
        #  pls_res
        if (is.null(best_pls_res)) best_pls_res <- pls_res
        
        X_tr_mod <- as.matrix(Xc_tr) %*% pls_res$Projection
        Xc_te    <- sweep(X[te,,drop=FALSE], 2, cen, "-")
        X_te_mod <- as.matrix(Xc_te) %*% pls_res$Projection
        
        # 3.  Partial Cox
        model <- partial_cox(
          X            = X_tr_mod,
          time         = time[tr],
          status       = status[tr],
          n_components = nc
        )
        
        # 4.  &  AUC 
        T_test <- transform_partial_cox(
          X_new      = X_te_mod,
          betas_list = model$betas_list,
          w_list     = model$w_list
        )
        lp_te <- predict(model$final_model, newdata=as.data.frame(T_test), type="lp")
        
        time_te   <- time[te]; status_te <- status[te]
        keep      <- !is.na(time_te)&!is.na(status_te)&!is.na(lp_te)
        time_te   <- time_te[keep]; status_te <- status_te[keep]
        marker    <- lp_te[keep]
        
        aucs <- .coxpcapls_auc_curve(
          time = time_te,
          status = status_te,
          marker = marker,
          auc_time_grid = auc_time_grid,
          auc_method = auc_method
        )
        
        fold_curves[fold, ] <- aucs
        fold_iAUCs[fold]    <- mean(aucs, na.rm=TRUE)
      }
      
      mean_curve <- colMeans(fold_curves, na.rm=TRUE)
      mean_iAUC  <- mean(fold_iAUCs, na.rm=TRUE)
      
      #  results
      results <- rbind(
        results,
        data.frame(
          ncomp     = nc,
          alpha     = a,
          iAUC      = mean_iAUC,
          auc_curve = I(list(mean_curve))
        )
      )
      key <- sprintf("%d_%.1f", nc, a)
      candidate_auc_curves[[key]] <- mean_curve
      
      # 
      if (mean_iAUC > best_iAUC) {
        best_iAUC       <- mean_iAUC
        best_ncomp      <- nc
        best_alpha      <- a
        best_model      <- model$final_model
        best_auc_curve  <- mean_curve
        best_betas_list <- model$betas_list
        best_w_list     <- model$w_list
        best_pls_res    <- pls_res
      }
    }
  }
  
  list(
    best_model           = best_model,
    best_ncomp           = best_ncomp,
    best_alpha           = best_alpha,
    best_iAUC            = best_iAUC,
    results              = results,
    candidate_auc_curves = candidate_auc_curves,
    auc_time_grid        = auc_time_grid,
    best_betas_list      = best_betas_list,
    best_w_list          = best_w_list,
    best_auc_curve       = best_auc_curve,
    best_pls_res         = best_pls_res
  )
}













# # "NNE"
# val_partial_cox_cv <- function(X, time, status,
#                                k = 5,
#                                ncomp_candidates = 1:5,
#                                auc_time_grid = NULL,
#                                p_thresh = 0.05) {
#   # ---------------------------
#   # Step 0 Cox p < p_thresh
#   # ---------------------------
#   p <- ncol(X)
#   significant_idx <- c()
#   for (j in 1:p) {
#     fit <- tryCatch({
#       coxph(Surv(time, status) ~ X[, j])
#     }, error = function(e) NULL)
#     
#     if (!is.null(fit)) {
#       pval <- summary(fit)$coefficients[,"Pr(>|z|)"]
#       if (!is.na(pval) && pval < p_thresh) {
#         significant_idx <- c(significant_idx, j)
#       }
#     }
#   }
#   
#   if (length(significant_idx) == 0) {
#     stop("No variables are significantly associated with survival at p < ", p_thresh)
#   }
#   
#   # 
#   X <- X[, significant_idx, drop = FALSE]
#   
#   # ---------------------------
#   # Step 1 auc 
#   # ---------------------------
#   n <- nrow(X)
#   set.seed(123)  # 
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
#   results <- data.frame(ncomp = integer(), iAUC = numeric(), auc_curve = I(list()))
#   candidate_auc_curves <- list()
#   
#   best_iAUC <- -Inf
#   best_ncomp <- NA
#   best_model <- NULL
#   best_betas_list <- NULL
#   best_w_list <- NULL
#   best_auc_curve <- NULL
#   
#   # ---------------------------
#   # Step 2 AUC  iAUC
#   # ---------------------------
#   for (nc in ncomp_candidates) {
#     cat("Evaluating ncomp =", nc, "\n")
#     
#     fold_auc_curves <- matrix(NA, nrow = k, ncol = length(auc_time_grid))
#     fold_iAUCs <- numeric(k)
#     
#     for (fold in 1:k) {
#       test_idx <- which(folds == fold)
#       train_idx <- setdiff(1:n, test_idx)
#       
#       #  X 
#       model <- partial_cox(X = X[train_idx, , drop = FALSE],
#                            time = time[train_idx],
#                            status = status[train_idx],
#                            n_components = nc)
#       
#       betas_list <- model$betas_list
#       w_list <- model$w_list
#       
#       T_test <- transform_partial_cox(X_new = X[test_idx, , drop = FALSE],
#                                       betas_list = betas_list,
#                                       w_list = w_list)
#       T_test <- as.data.frame(T_test)
#       colnames(T_test) <- paste0("T", seq_len(ncol(T_test)))
#       
#       marker_lp <- predict(model$final_model, newdata = T_test, type = "lp")
#       
#       #  auc 
#       # auc_vals <- sapply(auc_time_grid, function(tp) {
#       #   roc_obj <- survivalROC(Stime = time[test_idx], status = status[test_idx],
#       #                          marker = marker_lp, predict.time = tp, method = "NNE",
#       #                          span = 0.25)
#       #   roc_obj$AUC  })
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
#       fold_auc_curves[fold, ] <- auc_vals
#       fold_iAUCs[fold] <- mean(auc_vals, na.rm = TRUE)
#     }  # end of fold loop
#     
#     #  AUC  iAUC
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
#       best_model <- model$final_model  # 
#       best_betas_list <- betas_list
#       best_w_list <- w_list
#       best_auc_curve <- mean_auc_curve
#     }
#   }
#   
#   return(list(best_model = best_model,
#               best_ncomp = best_ncomp,
#               best_iAUC = best_iAUC,
#               results = results,
#               candidate_auc_curves = candidate_auc_curves,
#               auc_time_grid = auc_time_grid,
#               best_betas_list = best_betas_list,
#               best_w_list = best_w_list,
#               best_auc_curve = best_auc_curve))
# }
# 
# 
# 
# # "KM"
# val_partial_cox_cv <- function(X, time, status,
#                                k = 5,
#                                ncomp_candidates = 1:5,
#                                auc_time_grid = NULL,
#                                p_thresh = 0.05) {
#   # ---------------------------
#   # Step 0 Cox p < p_thresh
#   # ---------------------------
#   p <- ncol(X)
#   significant_idx <- c()
#   for (j in 1:p) {
#     fit <- tryCatch({
#       coxph(Surv(time, status) ~ X[, j])
#     }, error = function(e) NULL)
#     
#     if (!is.null(fit)) {
#       pval <- summary(fit)$coefficients[,"Pr(>|z|)"]
#       if (!is.na(pval) && pval < p_thresh) {
#         significant_idx <- c(significant_idx, j)
#       }
#     }
#   }
#   
#   if (length(significant_idx) == 0) {
#     stop("No variables are significantly associated with survival at p < ", p_thresh)
#   }
#   
#   # 
#   X <- X[, significant_idx, drop = FALSE]
#   
#   # ---------------------------
#   # Step 1 auc 
#   # ---------------------------
#   n <- nrow(X)
#   set.seed(123)  # 
#   folds <- sample(rep(1:k, length.out = n))
#   
#   if (is.null(auc_time_grid)) {
#     auc_time_grid <- seq(min(time), max(time), length.out = 50)
#   }
#   
#   results <- data.frame(ncomp = integer(), iAUC = numeric(), auc_curve = I(list()))
#   candidate_auc_curves <- list()
#   
#   best_iAUC <- -Inf
#   best_ncomp <- NA
#   best_model <- NULL
#   best_betas_list <- NULL
#   best_w_list <- NULL
#   best_auc_curve <- NULL
#   
#   # ---------------------------
#   # Step 2 AUC  iAUC
#   # ---------------------------
#   for (nc in ncomp_candidates) {
#     cat("Evaluating ncomp =", nc, "\n")
#     
#     fold_auc_curves <- matrix(NA, nrow = k, ncol = length(auc_time_grid))
#     fold_iAUCs <- numeric(k)
#     
#     for (fold in 1:k) {
#       test_idx <- which(folds == fold)
#       train_idx <- setdiff(1:n, test_idx)
#       
#       #  X 
#       model <- partial_cox(X = X[train_idx, , drop = FALSE],
#                            time = time[train_idx],
#                            status = status[train_idx],
#                            n_components = nc)
#       
#       betas_list <- model$betas_list
#       w_list <- model$w_list
#       
#       T_test <- transform_partial_cox(X_new = X[test_idx, , drop = FALSE],
#                                       betas_list = betas_list,
#                                       w_list = w_list)
#       T_test <- as.data.frame(T_test)
#       colnames(T_test) <- paste0("T", seq_len(ncol(T_test)))
#       
#       marker_lp <- predict(model$final_model, newdata = T_test, type = "lp")
#       
#       #  auc 
#       auc_vals <- sapply(auc_time_grid, function(tp) {
#         roc_obj <- survivalROC(Stime = time[test_idx], status = status[test_idx],
#                                marker = marker_lp, predict.time = tp, method = "KM")
#         roc_obj$AUC
#       })
#       
#       fold_auc_curves[fold, ] <- auc_vals
#       fold_iAUCs[fold] <- mean(auc_vals, na.rm = TRUE)
#     }  # end of fold loop
#     
#     #  AUC  iAUC
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
#       best_model <- model$final_model  # 
#       best_betas_list <- betas_list
#       best_w_list <- w_list
#       best_auc_curve <- mean_auc_curve
#     }
#   }
#   
#   return(list(best_model = best_model,
#               best_ncomp = best_ncomp,
#               best_iAUC = best_iAUC,
#               results = results,
#               candidate_auc_curves = candidate_auc_curves,
#               auc_time_grid = auc_time_grid,
#               best_betas_list = best_betas_list,
#               best_w_list = best_w_list,
#               best_auc_curve = best_auc_curve))
# }








# val_partial_cox_cv <- function(X, time, status,
#                                k = 5,
#                                ncomp_candidates = 1:5,
#                                auc_time_grid = NULL) {
# 
#   n <- nrow(X)
#   set.seed(123)  # 
#   folds <- sample(rep(1:k, length.out = n))
# 
#   # auc_time_grid50
#   if (is.null(auc_time_grid)) {
#     # 
#     auc_time_grid <- seq(min(time), max(time), length.out = 50)
#     #  auc_time_grid <- seq(quantile(time, 0.1), quantile(time, 0.9), length.out = 50)
#   }
# 
#   # 
#   results <- data.frame(ncomp = integer(), iAUC = numeric(), auc_curve = I(list()))
#   candidate_auc_curves <- list()
# 
#   best_iAUC <- -Inf
#   best_ncomp <- NA
#   best_model <- NULL
#   best_betas_list <- NULL
#   best_w_list <- NULL
#   best_auc_curve <- NULL
# 
#   # 
#   for (nc in ncomp_candidates) {
#     # 
#     pooled_time <- c()
#     pooled_status <- c()
#     pooled_marker <- c()
# 
#     # k
#     for (fold in 1:k) {
#       test_idx <- which(folds == fold)
#       train_idx <- setdiff(1:n, test_idx)
# 
#       # partial cox
#       model <- partial_cox(X = X[train_idx, ],
#                            time = time[train_idx],
#                            status = status[train_idx],
#                            n_components = nc)
# 
#       betas_list <- model$betas_list
#       w_list <- model$w_list
# 
#       # Cox
#       T_test <- transform_partial_cox(X_new = X[test_idx, ],
#                                       betas_list = betas_list,
#                                       w_list = w_list)
#       T_test <- as.data.frame(T_test)
#       colnames(T_test) <- paste0("T", seq_len(ncol(T_test)))
# 
#       # Cox
#       marker_lp <- predict(model$final_model, newdata = T_test, type = "lp")
# 
#       # 
#       pooled_time <- c(pooled_time, time[test_idx])
#       pooled_status <- c(pooled_status, status[test_idx])
#       pooled_marker <- c(pooled_marker, marker_lp)
#     }  # end of fold loop
# 
#     #  auc_time_grid  AUC 
#     auc_vals <- sapply(auc_time_grid, function(tp) {
#       roc_obj <- survivalROC(Stime = pooled_time, status = pooled_status,
#                              marker = pooled_marker, predict.time = tp, method = "KM")
#       return(roc_obj$AUC)
#     })
# 
#     # AUCAUC
#     iAUC <- mean(auc_vals, na.rm = TRUE)
# 
#     # 
#     results <- rbind(results, data.frame(ncomp = nc, iAUC = iAUC, auc_curve = I(list(auc_vals))))
#     candidate_auc_curves[[as.character(nc)]] <- auc_vals
# 
#     #  AUC 
#     if (iAUC > best_iAUC) {
#       best_iAUC <- iAUC
#       best_ncomp <- nc
#       best_model <- model$final_model
#       best_betas_list <- betas_list
#       best_w_list <- w_list
#       best_auc_curve <- auc_vals
#     }
#   }  # end candidate loop
# 
#   return(list(best_model = best_model,
#               best_ncomp = best_ncomp,
#               best_iAUC = best_iAUC,
#               results = results,
#               candidate_auc_curves = candidate_auc_curves,
#               auc_time_grid = auc_time_grid,
#               best_betas_list = best_betas_list,
#               best_w_list = best_w_list,
#               best_auc_curve = best_auc_curve))
# }




# val_partial_cox <- function(X_train, time_train, status_train,
#                             X_test, time_test, status_test,
#                             ncomp_candidates = 1:5,
#                             auc_time_grid = NULL) {
# 
#   if (is.null(auc_time_grid)) {
#     auc_time_grid <- seq(min(c(time_train, time_test)), max(c(time_train, time_test)), length.out = 50)
#   }
#   
# 
#   #  results  "auc_curve"  AUC 
#   results <- data.frame(ncomp = integer(0), iAUC = numeric(0), auc_curve = I(list()))
#   candidate_auc_curves <- list()
#   best_iAUC <- -Inf
#   best_ncomp <- NA
#   best_model <- NULL
#   best_betas_list <- NULL
#   best_w_list <- NULL
#   best_auc_curve <- NULL
# 
#   # 
#   for (nc in ncomp_candidates) {
#     
#     #  partial_cox 
#     model <- partial_cox(X = X_train, time = time_train, status = status_train,
#                          n_components = nc)
#     betas_list <- model$betas_list
#     w_list <- model$w_list
# 
#     # 
#     T_test <- transform_partial_cox(X_test, betas_list, w_list)
#     colnames(T_test) <- paste0("T", 1:ncol(T_test))
#     T_test <- as.data.frame(T_test)
# 
#     #  Cox 
#     marker_lp <- predict(model$final_model, newdata = T_test, type = "lp")
# 
#     #  AUC
#     auc_vals <- numeric(length(auc_time_grid))
#     for (i in seq_along(auc_time_grid)) {
#       tp <- auc_time_grid[i]
#       roc_obj <- survivalROC(Stime = time_test, status = status_test, marker = marker_lp,
#                              predict.time = tp, method = "KM")
#       auc_vals[i] <- roc_obj$AUC
#     }
#     #  AUC AUC 
#     iAUC <- mean(auc_vals, na.rm = TRUE)
# 
#     # auc_curve  list  auc_vals
#     results <- rbind(results, data.frame(ncomp = nc, iAUC = iAUC, auc_curve = I(list(auc_vals))))
#     candidate_auc_curves[[as.character(nc)]] <- auc_vals
# 
#     #  iAUC 
#     if (iAUC > best_iAUC) {
#       best_iAUC <- iAUC
#       best_ncomp <- nc
#       best_model <- model$final_model
#       best_betas_list <- betas_list
#       best_w_list <- w_list
#       best_auc_curve <- auc_vals
#     }
#   }
# 
#   return(list(best_model = best_model,
#               best_ncomp = best_ncomp,
#               best_iAUC = best_iAUC,
#               results = results,
#               candidate_auc_curves = candidate_auc_curves,
#               auc_time_grid = auc_time_grid,
#               best_betas_list = best_betas_list,
#               best_w_list = best_w_list,
#               best_auc_curve = best_auc_curve))
# }


# val_partial_cox <- function(X_train, time_train, status_train,
#                             X_test, time_test, status_test,
#                             ncomp_candidates = 1:5,
#                             auc_time_grid = NULL) {
#   
#   # Step 0 Cox p < 0.05
#   p <- ncol(X_train)
#   significant_idx <- c()
#   for (j in 1:p) {
#     fit <- tryCatch({
#       coxph(Surv(time_train, status_train) ~ X_train[, j])
#     }, error = function(e) NULL)
#     
#     if (!is.null(fit)) {
#       pval <- summary(fit)$coefficients[,"Pr(>|z|)"]
#       if (!is.na(pval) && pval < 0.05) {
#         significant_idx <- c(significant_idx, j)
#       }
#     }
#   }
#   
#   if (length(significant_idx) == 0) {
#     stop("No variables are significantly associated with survival at p < 0.05.")
#   }
#   
#   # 
#   X_train <- X_train[, significant_idx, drop = FALSE]
#   X_test  <- X_test[, significant_idx, drop = FALSE]
#   
#   # Step 1 AUC 
#   if (is.null(auc_time_grid)) {
#     auc_time_grid <- seq(min(c(time_train, time_test)),
#                          max(c(time_train, time_test)),
#                          length.out = 50)
#   }
#   
#   # 
#   results <- data.frame(ncomp = integer(0), iAUC = numeric(0), auc_curve = I(list()))
#   candidate_auc_curves <- list()
#   best_iAUC <- -Inf
#   best_ncomp <- NA
#   best_model <- NULL
#   best_betas_list <- NULL
#   best_w_list <- NULL
#   best_auc_curve <- NULL
#   
#   # Step 2 iAUC 
#   for (nc in ncomp_candidates) {
#     model <- partial_cox(X = X_train, time = time_train, status = status_train,
#                          n_components = nc)
#     betas_list <- model$betas_list
#     w_list <- model$w_list
#     
#     T_test <- transform_partial_cox(X_test, betas_list, w_list)
#     colnames(T_test) <- paste0("T", 1:ncol(T_test))
#     T_test <- as.data.frame(T_test)
#     
#     marker_lp <- predict(model$final_model, newdata = T_test, type = "lp")
#     
#     auc_vals <- numeric(length(auc_time_grid))
#     for (i in seq_along(auc_time_grid)) {
#       tp <- auc_time_grid[i]
#       roc_obj <- survivalROC(Stime = time_test, status = status_test,
#                              marker = marker_lp,
#                              predict.time = tp, method = "KM")
#       auc_vals[i] <- roc_obj$AUC
#     }
#     
#     iAUC <- mean(auc_vals, na.rm = TRUE)
#     
#     results <- rbind(results,
#                      data.frame(ncomp = nc, iAUC = iAUC, auc_curve = I(list(auc_vals))))
#     candidate_auc_curves[[as.character(nc)]] <- auc_vals
#     
#     if (iAUC > best_iAUC) {
#       best_iAUC <- iAUC
#       best_ncomp <- nc
#       best_model <- model$final_model
#       best_betas_list <- betas_list
#       best_w_list <- w_list
#       best_auc_curve <- auc_vals
#     }
#   }
#   
#   return(list(
#     best_model = best_model,
#     best_ncomp = best_ncomp,
#     best_iAUC = best_iAUC,
#     results = results,
#     candidate_auc_curves = candidate_auc_curves,
#     auc_time_grid = auc_time_grid,
#     best_betas_list = best_betas_list,
#     best_w_list = best_w_list,
#     best_auc_curve = best_auc_curve,
#     selected_variables = colnames(X_train)
#   ))
# }
