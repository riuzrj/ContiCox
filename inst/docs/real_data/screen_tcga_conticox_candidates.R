project_dir <- "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls"
setwd(project_dir)

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(edgeR)
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

pick <- function(df, candidates) {
  for (nm in candidates) {
    if (nm %in% names(df)) {
      return(df[[nm]])
    }
  }
  rep(NA, nrow(df))
}

prepare_tcga_project <- function(project, cache_dir = file.path(project_dir, "data_cache")) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_file <- file.path(cache_dir, paste0(project, "_prepared.RData"))

  if (file.exists(cache_file)) {
    load(cache_file)
    return(list(X = X, time = time, status = status, project = project))
  }

  query <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  GDCdownload(query, files.per.chunk = 40)
  se <- GDCprepare(query, summarizedExperiment = TRUE)

  counts <- SummarizedExperiment::assay(se)
  barcodes <- colnames(counts)
  keep_col <- substr(barcodes, 14, 15) == "01"
  counts <- counts[, keep_col, drop = FALSE]

  dge <- edgeR::DGEList(counts = counts)
  cpm0 <- edgeR::cpm(dge)
  keep_gene <- rowMeans(cpm0 > 1) >= 0.20
  dge <- dge[keep_gene, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  log_cpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  expr_mat <- t(log_cpm)

  meta <- as.data.frame(SummarizedExperiment::colData(se))
  meta <- meta[keep_col, , drop = FALSE]
  d2d <- suppressWarnings(as.numeric(pick(meta, c(
    "days_to_death", "days_to_death.x", "paper_days_to_death"
  ))))
  d2f <- suppressWarnings(as.numeric(pick(meta, c(
    "days_to_last_follow_up", "days_to_last_followup", "days_to_lastfollowup",
    "days_to_last_follow_up.x", "paper_days_to_last_followup"
  ))))
  vital_status <- tolower(as.character(pick(meta, c(
    "vital_status", "vital_status.x", "paper_vital_status"
  ))))

  time <- ifelse(!is.na(d2d), d2d, d2f)
  status <- ifelse(!is.na(d2d) & d2d > 0, 1,
                   ifelse(vital_status %in% c("dead", "deceased"), 1, 0))

  madv <- matrixStats::colMads(expr_mat, na.rm = TRUE)
  sel <- which(is.finite(madv) & madv > 0)
  X0 <- scale(expr_mat[, sel, drop = FALSE])
  ok <- is.finite(time) &
    is.finite(status) &
    time > 0 &
    rowSums(is.finite(X0)) == ncol(X0)

  X <- X0[ok, , drop = FALSE]
  time <- time[ok]
  status <- status[ok]
  save(X, time, status, project, file = cache_file)

  list(X = X, time = time, status = status, project = project)
}

run_tcga_project <- function(project,
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
                             use_preselection = TRUE,
                             residual_type = "martingale",
                             auc_method = "NNE") {
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  dat <- prepare_tcga_project(project)
  X <- dat$X
  time <- dat$time
  status <- as.integer(dat$status != 0)

  cat("\n==", project, "==\n")
  cat("Patients:", nrow(X), "\n")
  cat("Genes:", ncol(X), "\n")
  cat("Events:", sum(status == 1), "\n")
  cat("Censoring rate:", mean(status == 0), "\n")
  cat("Follow-up range:", paste(range(time, na.rm = TRUE), collapse = " - "), "\n")

  split <- .cvtest_make_split(status, test_prop = test_prop, seed = seed + 1L)
  train_idx <- split$train_idx
  test_idx <- split$test_idx

  event_time <- time[status == 1]
  auc_time_grid <- as.numeric(stats::quantile(
    event_time,
    probs = seq(0.1, 0.8, length.out = 30),
    na.rm = TRUE
  ))
  auc_time_grid <- sort(unique(auc_time_grid))

  auc_curves_list_drpls <- val_drpls_cv_test(
    X, time, status,
    test_prop = test_prop,
    train_idx = train_idx,
    test_idx = test_idx,
    k = k,
    ncomp_candidates = ncomp_candidates_DRPLS,
    auc_time_grid = auc_time_grid,
    use_preselection = use_preselection,
    seed = seed,
    auc_method = auc_method
  )
  cat("PLSDR test iAUC:", auc_curves_list_drpls$test_iAUC, "\n")

  auc_curves_list_drpcapls <- val_conticox_cv_test(
    X, time, status,
    test_prop = test_prop,
    train_idx = train_idx,
    test_idx = test_idx,
    k = k,
    ncomp_candidates = ncomp_candidates_DRPLS,
    alpha_candidates = alpha_candidates,
    auc_time_grid = auc_time_grid,
    use_preselection = use_preselection,
    normalize_alpha_terms = TRUE,
    residual_type = residual_type,
    seed = seed,
    auc_method = auc_method
  )
  cat("ContiCox test iAUC:", auc_curves_list_drpcapls$test_iAUC, "\n")

  auc_curves_list_pcox <- val_partial_cox_cv_test(
    X, time, status,
    test_prop = test_prop,
    train_idx = train_idx,
    test_idx = test_idx,
    k = k,
    ncomp_candidates = ncomp_candidates_pcox,
    auc_time_grid = auc_time_grid,
    use_preselection = use_preselection,
    seed = seed,
    auc_method = auc_method
  )
  cat("pCox test iAUC:", auc_curves_list_pcox$test_iAUC, "\n")

  auc_curves_list_ridgecox <- val_ridge_cox_cv_test(
    X, time, status,
    test_prop = test_prop,
    train_idx = train_idx,
    test_idx = test_idx,
    k = k,
    lambda_candidates = lambda_candidates_ridgecox,
    auc_time_grid = auc_time_grid,
    use_preselection = use_preselection,
    seed = seed,
    auc_method = auc_method
  )
  cat("Ridge-Cox test iAUC:", auc_curves_list_ridgecox$test_iAUC, "\n")

  auc_curves_list_enetcox <- val_elastic_net_cox_cv_test(
    X, time, status,
    test_prop = test_prop,
    train_idx = train_idx,
    test_idx = test_idx,
    k = k,
    alpha_candidates = alpha_candidates_enetcox,
    lambda_candidates = lambda_candidates_enetcox,
    auc_time_grid = auc_time_grid,
    use_preselection = use_preselection,
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
    auc_time_grid = auc_time_grid,
    train_idx = train_idx,
    test_idx = test_idx,
    test_prop = test_prop,
    PLSDR = auc_curves_list_drpls,
    ContiCox = auc_curves_list_drpcapls,
    pCox = auc_curves_list_pcox,
    RidgeCox = auc_curves_list_ridgecox,
    ElasticNetCox = auc_curves_list_enetcox
  )

  result_file <- file.path(result_dir, paste0(project, "_validation_results.RData"))
  csv_file <- file.path(result_dir, paste0(project, "_summary.csv"))
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
  run_tcga_project(project)
}
