# Shared time-dependent AUC utilities.

.coxpcapls_auc_curve <- function(time, status, marker, auc_time_grid,
                                 auc_method = c("IPCW", "NNE", "KM"),
                                 ipcw_weighting = "marginal",
                                 nne_span = 0.30) {
  auc_method <- match.arg(auc_method)

  time <- as.numeric(time)
  status <- as.integer(status != 0L)
  marker <- as.numeric(marker)
  auc_time_grid <- sort(unique(as.numeric(auc_time_grid[is.finite(auc_time_grid)])))

  if (length(auc_time_grid) == 0L) {
    stop("auc_time_grid must contain at least one finite value.")
  }

  ok <- is.finite(time) & is.finite(status) & is.finite(marker) & time > 0
  if (sum(ok) < 3L || length(unique(status[ok])) < 2L) {
    return(rep(NA_real_, length(auc_time_grid)))
  }

  time <- time[ok]
  status <- status[ok]
  marker <- marker[ok]

  if (auc_method == "IPCW") {
    if (!requireNamespace("timeROC", quietly = TRUE)) {
      stop("Package `timeROC` is required for IPCW time-dependent AUC.")
    }
    if (!requireNamespace("survival", quietly = TRUE)) {
      stop("Package `survival` is required for IPCW time-dependent AUC.")
    }

    roc <- tryCatch(
      timeROC::timeROC(
        T = time,
        delta = status,
        marker = marker,
        cause = 1,
        weighting = ipcw_weighting,
        times = auc_time_grid,
        iid = FALSE
      ),
      error = function(e) {
        warning("IPCW AUC failed: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(roc)) {
      return(rep(NA_real_, length(auc_time_grid)))
    }

    auc <- as.numeric(roc$AUC)
    if (length(auc) == length(auc_time_grid)) {
      return(auc)
    }

    out <- rep(NA_real_, length(auc_time_grid))
    idx <- match(auc_time_grid, as.numeric(roc$times))
    keep <- !is.na(idx)
    out[keep] <- auc[idx[keep]]
    return(out)
  }

  if (!requireNamespace("survivalROC", quietly = TRUE)) {
    stop("Package `survivalROC` is required for NNE/KM time-dependent AUC.")
  }

  vapply(auc_time_grid, function(tp) {
    tryCatch({
      args <- list(
        Stime = time,
        status = status,
        marker = marker,
        predict.time = tp,
        method = auc_method
      )
      if (auc_method == "NNE") {
        args$span <- nne_span
      }
      as.numeric(do.call(survivalROC::survivalROC, args)$AUC)
    }, error = function(e) {
      warning("AUC failed at t = ", tp, ": ", conditionMessage(e))
      NA_real_
    })
  }, numeric(1))
}

.coxpcapls_iAUC <- function(auc_curve) {
  mean(auc_curve, na.rm = TRUE)
}
