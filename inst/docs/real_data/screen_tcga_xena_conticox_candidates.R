project_dir <- "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls"
setwd(project_dir)

suppressPackageStartupMessages({
  library(matrixStats)
  library(survival)
})

method_files <- c(
  file.path(project_dir, "R", "pcapls_CR.R"),
  file.path(project_dir, "R", "DRPLS.R"),
  file.path(project_dir, "R", "DRPCAPLS.R"),
  file.path(project_dir, "R", "partial_cox.R"),
  file.path(project_dir, "R", "penalized_cox_baselines.R"),
  file.path(project_dir, "R", "train_test_validation.R")
)
for (method_file in method_files) {
  source(method_file)
}

xena_base_url <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download"

download_xena_file <- function(project, suffix, cache_dir) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(cache_dir, paste0(project, suffix))
  if (file.exists(out_file)) {
    return(out_file)
  }

  url <- paste0(xena_base_url, "/", project, suffix)
  message("Downloading ", url)
  utils::download.file(url, out_file, mode = "wb", quiet = FALSE)
  out_file
}

prepare_xena_tcga_project <- function(project,
                                      cache_dir = file.path(project_dir, "data_cache"),
                                      prepared_dir = file.path(project_dir, "data_cache")) {
  dir.create(prepared_dir, recursive = TRUE, showWarnings = FALSE)
  prepared_file <- file.path(prepared_dir, paste0(project, "_xena_prepared.RData"))

  if (file.exists(prepared_file)) {
    load(prepared_file)
    return(list(X = X, time = time, status = status, project = project))
  }

  expr_file <- download_xena_file(project, ".star_counts.tsv.gz", cache_dir)
  survival_file <- download_xena_file(project, ".survival.tsv.gz", cache_dir)

  expr <- utils::read.delim(gzfile(expr_file), check.names = FALSE)
  gene_id <- expr[[1]]
  expr_mat <- as.matrix(expr[-1])
  storage.mode(expr_mat) <- "numeric"
  rownames(expr_mat) <- gene_id
  colnames(expr_mat) <- names(expr)[-1]

  keep_col <- substr(colnames(expr_mat), 14, 15) == "01"
  expr_mat <- expr_mat[, keep_col, drop = FALSE]

  survival_df <- utils::read.delim(gzfile(survival_file), check.names = FALSE)
  survival_df$sample_short <- substr(survival_df$sample, 1, 15)
  expr_short <- substr(colnames(expr_mat), 1, 15)
  match_idx <- match(expr_short, survival_df$sample_short)
  keep_match <- !is.na(match_idx)
  expr_mat <- expr_mat[, keep_match, drop = FALSE]
  match_idx <- match_idx[keep_match]

  time <- as.numeric(survival_df$OS.time[match_idx])
  status <- as.integer(survival_df$OS[match_idx] != 0)

  keep_gene <- rowMeans(expr_mat > 1, na.rm = TRUE) >= 0.20
  expr_mat <- expr_mat[keep_gene, , drop = FALSE]
  madv <- matrixStats::rowMads(expr_mat, na.rm = TRUE)
  expr_mat <- expr_mat[is.finite(madv) & madv > 0, , drop = FALSE]

  X0 <- scale(t(expr_mat))
  ok <- is.finite(time) &
    is.finite(status) &
    time > 0 &
    rowSums(is.finite(X0)) == ncol(X0)

  X <- X0[ok, , drop = FALSE]
  time <- time[ok]
  status <- status[ok]

  save(X, time, status, project, file = prepared_file)
  list(X = X, time = time, status = status, project = project)
}

run_xena_tcga_project <- function(project,
                                  result_dir = file.path(project_dir, "candidate_results"),
                                  seed = 2026,
                                  test_prop = 0.3,
                                  k = 5,
                                  ncomp_candidates_DRPLS = 1:10,
                                  ncomp_candidates_pcox = 2:10,
                                  alpha_candidates = seq(0, 1, by = 0.1),
                                  lambda_candidates_ridgecox = 10 ^ seq(-4, 1, by = 1),
                                  alpha_candidates_enetcox = c(0.25, 0.5, 0.75),
                                  lambda_candidates_enetcox = 10 ^ seq(-5, 1, by = 1),
                                  p_thresh = 0.05,
                                  residual_type = "martingale",
                                  auc_method = "NNE") {
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  dat <- prepare_xena_tcga_project(project)
  X <- dat$X
  time <- dat$time
  status <- as.integer(dat$status != 0)

  cat("\n==", project, "Xena quick screen ==\n")
  cat("Patients:", nrow(X), "\n")
  cat("Genes before Cox preselection:", ncol(X), "\n")
  cat("Events:", sum(status == 1), "\n")
  cat("Censoring rate:", mean(status == 0), "\n")
  cat("Follow-up range:", paste(range(time, na.rm = TRUE), collapse = " - "), "\n")

  split <- .cvtest_make_split(status, test_prop = test_prop, seed = seed + 1L)
  train_idx <- split$train_idx
  test_idx <- split$test_idx

  cat("Cox preselection on training set...\n")
  feature_idx <- .cvtest_select_features(
    X = X[train_idx, , drop = FALSE],
    time = time[train_idx],
    status = status[train_idx],
    use_preselection = TRUE,
    p_thresh = p_thresh
  )
  cat("Genes after Cox preselection:", length(feature_idx), "\n")
  X_work <- X[, feature_idx, drop = FALSE]

  event_time <- time[status == 1]
  auc_time_grid <- as.numeric(stats::quantile(
    event_time,
    probs = seq(0.1, 0.8, length.out = 30),
    na.rm = TRUE
  ))
  auc_time_grid <- sort(unique(auc_time_grid))

  auc_curves_list_drpls <- val_drpls_cv_test(
    X_work, time, status,
    test_prop = test_prop,
    train_idx = train_idx,
    test_idx = test_idx,
    k = k,
    ncomp_candidates = ncomp_candidates_DRPLS,
    auc_time_grid = auc_time_grid,
    use_preselection = FALSE,
    seed = seed,
    auc_method = auc_method
  )
  cat("PLSDR test iAUC:", auc_curves_list_drpls$test_iAUC, "\n")

  auc_curves_list_drpcapls <- val_conticox_cv_test(
    X_work, time, status,
    test_prop = test_prop,
    train_idx = train_idx,
    test_idx = test_idx,
    k = k,
    ncomp_candidates = ncomp_candidates_DRPLS,
    alpha_candidates = alpha_candidates,
    auc_time_grid = auc_time_grid,
    use_preselection = FALSE,
    normalize_alpha_terms = TRUE,
    residual_type = residual_type,
    seed = seed,
    auc_method = auc_method
  )
  cat("ContiCox test iAUC:", auc_curves_list_drpcapls$test_iAUC, "\n")

  auc_curves_list_pcox <- val_partial_cox_cv_test(
    X_work, time, status,
    test_prop = test_prop,
    train_idx = train_idx,
    test_idx = test_idx,
    k = k,
    ncomp_candidates = ncomp_candidates_pcox,
    auc_time_grid = auc_time_grid,
    use_preselection = FALSE,
    seed = seed,
    auc_method = auc_method
  )
  cat("pCox test iAUC:", auc_curves_list_pcox$test_iAUC, "\n")

  auc_curves_list_ridgecox <- val_ridge_cox_cv_test(
    X_work, time, status,
    test_prop = test_prop,
    train_idx = train_idx,
    test_idx = test_idx,
    k = k,
    lambda_candidates = lambda_candidates_ridgecox,
    auc_time_grid = auc_time_grid,
    use_preselection = FALSE,
    seed = seed,
    auc_method = auc_method
  )
  cat("Ridge-Cox test iAUC:", auc_curves_list_ridgecox$test_iAUC, "\n")

  auc_curves_list_enetcox <- val_elastic_net_cox_cv_test(
    X_work, time, status,
    test_prop = test_prop,
    train_idx = train_idx,
    test_idx = test_idx,
    k = k,
    alpha_candidates = alpha_candidates_enetcox,
    lambda_candidates = lambda_candidates_enetcox,
    auc_time_grid = auc_time_grid,
    use_preselection = FALSE,
    seed = seed,
    auc_method = auc_method
  )
  cat("EN-Cox test iAUC:", auc_curves_list_enetcox$test_iAUC, "\n")

  summary <- data.frame(
    project = project,
    method = c("PLSDR", "ContiCox", "pCox", "Ridge-Cox", "EN-Cox"),
    test_iAUC = c(
      auc_curves_list_drpls$test_iAUC,
      auc_curves_list_drpcapls$test_iAUC,
      auc_curves_list_pcox$test_iAUC,
      auc_curves_list_ridgecox$test_iAUC,
      auc_curves_list_enetcox$test_iAUC
    ),
    stringsAsFactors = FALSE
  )
  summary <- summary[order(-summary$test_iAUC), ]
  print(summary, row.names = FALSE)

  tcga_auc_results <- list(
    dataset = project,
    source = "UCSC Xena GDC hub star_counts",
    auc_time_grid = auc_time_grid,
    train_idx = train_idx,
    test_idx = test_idx,
    test_prop = test_prop,
    feature_idx = feature_idx,
    PLSDR = auc_curves_list_drpls,
    ContiCox = auc_curves_list_drpcapls,
    pCox = auc_curves_list_pcox,
    RidgeCox = auc_curves_list_ridgecox,
    ElasticNetCox = auc_curves_list_enetcox
  )

  result_file <- file.path(result_dir, paste0(project, "_xena_validation_results.RData"))
  csv_file <- file.path(result_dir, paste0(project, "_xena_summary.csv"))
  save(
    auc_curves_list_drpls,
    auc_curves_list_drpcapls,
    auc_curves_list_pcox,
    auc_curves_list_ridgecox,
    auc_curves_list_enetcox,
    tcga_auc_results,
    summary,
    file = result_file
  )
  utils::write.csv(summary, csv_file, row.names = FALSE)
  invisible(summary)
}

args <- commandArgs(trailingOnly = TRUE)
projects <- if (length(args) == 0L) "TCGA-LGG" else args
for (project in projects) {
  run_xena_tcga_project(project)
}
