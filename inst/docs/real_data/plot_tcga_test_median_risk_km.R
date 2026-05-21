# Plot Kaplan-Meier curves from held-out test risk scores.
#
# Usage:
#   1. Prepare X, time, and status in the current R session.
#   2. Source this file to refit cv-test models and draw the figure:
#        source("plot_tcga_test_median_risk_km.R")
#
# If cv-test result objects already exist, the script reuses them instead of
# refitting. Set FORCE_REFIT_TCGA_TEST_MEDIAN_RISK_KM=1 to force refitting.
#
# The script uses each method's test_marker and splits test patients by the
# median test-set risk score. Higher Cox linear predictor means higher risk.

project_dir <- "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls"
default_figure_dir <- file.path(project_dir, "figures")
default_output_base <- file.path(default_figure_dir, "tcga_test_median_risk_km")

.require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package `", pkg, "` is required.")
  }
}

.first_existing_object <- function(object_names, envir) {
  for (object_name in object_names) {
    if (exists(object_name, envir = envir, inherits = TRUE)) {
      return(get(object_name, envir = envir, inherits = TRUE))
    }
  }
  NULL
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

source_cvtest_method_files <- function(project_root = project_dir) {
  method_files <- c(
    file.path(project_root, "R", "auc_utils.R"),
    file.path(project_root, "R", "pcapls_CR.R"),
    file.path(project_root, "R", "DRPLS.R"),
    file.path(project_root, "R", "DRPCAPLS.R"),
    file.path(project_root, "R", "partial_cox.R"),
    file.path(project_root, "R", "penalized_cox_baselines.R"),
    file.path(project_root, "R", "train_test_validation.R")
  )
  missing_files <- method_files[!file.exists(method_files)]
  if (length(missing_files) > 0L) {
    stop("Cannot find method source files: ", paste(missing_files, collapse = ", "))
  }
  for (method_file in method_files) {
    source(method_file)
  }
  invisible(method_files)
}

prepare_survival_xy <- function(X, time, status) {
  X <- as.matrix(X)
  time <- as.numeric(time)
  status <- as.integer(status != 0L)

  if (nrow(X) != length(time) || length(time) != length(status)) {
    stop("nrow(X), length(time), and length(status) must be identical.")
  }

  ok <- is.finite(time) & time > 0 & is.finite(status) &
    rowSums(is.finite(X)) == ncol(X)

  if (sum(ok) < 10L) {
    stop("Too few valid samples after removing invalid survival or expression rows.")
  }
  if (any(!ok)) {
    message(sum(!ok), " samples removed before fitting because of invalid time/status/X.")
  }

  list(
    X = X[ok, , drop = FALSE],
    time = time[ok],
    status = status[ok],
    kept = which(ok)
  )
}

make_train_event_auc_grid <- function(time, status, train_idx,
                                      probs = seq(0.2, 0.8, length.out = 30)) {
  event_time <- time[train_idx][status[train_idx] == 1L]
  if (length(event_time) >= 3L) {
    grid <- stats::quantile(event_time, probs = probs, na.rm = TRUE, names = FALSE)
  } else {
    grid <- stats::quantile(time[train_idx], probs = probs, na.rm = TRUE, names = FALSE)
  }
  grid <- sort(unique(as.numeric(grid[is.finite(grid) & grid > 0])))
  if (length(grid) == 0L) {
    stop("Cannot construct auc_time_grid from training survival times.")
  }
  grid
}

fit_tcga_cvtest_results <- function(
    X,
    time,
    status,
    test_prop = 0.3,
    seed = 2026,
    k = 5,
    auc_method = c("IPCW", "NNE", "KM"),
    auc_time_grid = NULL,
    auc_time_probs = seq(0.2, 0.8, length.out = 30),
    ncomp_candidates_DRPLS = 1:10,
    ncomp_candidates_pcox = 2:10,
    alpha_candidates = seq(0, 1, by = 0.1),
    lambda_candidates_ridgecox = 10 ^ seq(-4, 1, by = 1),
    alpha_candidates_enetcox = c(0.25, 0.5, 0.75),
    lambda_candidates_enetcox = 10 ^ seq(-5, 1, by = 1),
    use_preselection = TRUE,
    p_thresh = 0.05,
    residual_type = "martingale",
    normalize_alpha_terms = TRUE,
    methods = c("PLSDR", "ContiCox", "pCox", "Ridge-Cox", "EN-Cox"),
    quiet = FALSE,
    source_methods = TRUE
) {
  auc_method <- match.arg(auc_method)
  methods <- unique(methods)

  if (isTRUE(source_methods)) {
    source_cvtest_method_files(project_root = project_dir)
  }

  dat <- prepare_survival_xy(X, time, status)
  X <- dat$X
  time <- dat$time
  status <- dat$status

  split <- .cvtest_make_split(status, test_prop = test_prop, seed = seed)
  train_idx <- split$train_idx
  test_idx <- split$test_idx

  if (is.null(auc_time_grid)) {
    auc_time_grid <- make_train_event_auc_grid(
      time = time,
      status = status,
      train_idx = train_idx,
      probs = auc_time_probs
    )
  } else {
    auc_time_grid <- sort(unique(as.numeric(auc_time_grid[is.finite(auc_time_grid) & auc_time_grid > 0])))
  }

  cat("Training samples:", length(train_idx), "\n")
  cat("Test samples:", length(test_idx), "\n")
  cat("Test events:", sum(status[test_idx] == 1L), "\n")
  cat("AUC method:", auc_method, "\n")

  results <- list()

  if ("PLSDR" %in% methods) {
    results[["PLSDR"]] <- .quiet_eval(val_drpls_cv_test(
      X, time, status,
      test_prop = test_prop,
      train_idx = train_idx,
      test_idx = test_idx,
      k = k,
      ncomp_candidates = ncomp_candidates_DRPLS,
      auc_time_grid = auc_time_grid,
      use_preselection = use_preselection,
      p_thresh = p_thresh,
      auc_method = auc_method,
      seed = seed
    ), quiet = quiet)
  }

  if ("ContiCox" %in% methods) {
    results[["ContiCox"]] <- .quiet_eval(val_conticox_cv_test(
      X, time, status,
      test_prop = test_prop,
      train_idx = train_idx,
      test_idx = test_idx,
      k = k,
      ncomp_candidates = ncomp_candidates_DRPLS,
      alpha_candidates = alpha_candidates,
      auc_time_grid = auc_time_grid,
      use_preselection = use_preselection,
      p_thresh = p_thresh,
      auc_method = auc_method,
      seed = seed,
      residual_type = residual_type,
      normalize_alpha_terms = normalize_alpha_terms
    ), quiet = quiet)
  }

  if ("pCox" %in% methods) {
    results[["pCox"]] <- .quiet_eval(val_partial_cox_cv_test(
      X, time, status,
      test_prop = test_prop,
      train_idx = train_idx,
      test_idx = test_idx,
      k = k,
      ncomp_candidates = ncomp_candidates_pcox,
      auc_time_grid = auc_time_grid,
      use_preselection = use_preselection,
      p_thresh = p_thresh,
      auc_method = auc_method,
      seed = seed
    ), quiet = quiet)
  }

  if ("Ridge-Cox" %in% methods) {
    results[["Ridge-Cox"]] <- .quiet_eval(val_ridge_cox_cv_test(
      X, time, status,
      test_prop = test_prop,
      train_idx = train_idx,
      test_idx = test_idx,
      k = k,
      lambda_candidates = lambda_candidates_ridgecox,
      auc_time_grid = auc_time_grid,
      use_preselection = use_preselection,
      p_thresh = p_thresh,
      auc_method = auc_method,
      seed = seed
    ), quiet = quiet)
  }

  if ("EN-Cox" %in% methods) {
    results[["EN-Cox"]] <- .quiet_eval(val_elastic_net_cox_cv_test(
      X, time, status,
      test_prop = test_prop,
      train_idx = train_idx,
      test_idx = test_idx,
      k = k,
      alpha_candidates = alpha_candidates_enetcox,
      lambda_candidates = lambda_candidates_enetcox,
      auc_time_grid = auc_time_grid,
      use_preselection = use_preselection,
      p_thresh = p_thresh,
      auc_method = auc_method,
      seed = seed
    ), quiet = quiet)
  }

  list(
    X = X,
    time = time,
    status = status,
    train_idx = train_idx,
    test_idx = test_idx,
    auc_time_grid = auc_time_grid,
    auc_method = auc_method,
    results = results,
    kept = dat$kept
  )
}

collect_tcga_cvtest_results <- function(envir = parent.frame()) {
  object_map <- c(
    PLSDR = "auc_curves_list_drpls",
    ContiCox = "auc_curves_list_drpcapls",
    pCox = "auc_curves_list_pcox",
    "Ridge-Cox" = "auc_curves_list_ridgecox",
    "EN-Cox" = "auc_curves_list_enetcox"
  )

  results <- list()
  for (method_name in names(object_map)) {
    object_name <- unname(object_map[method_name])
    if (exists(object_name, envir = envir, inherits = TRUE)) {
      results[[method_name]] <- get(object_name, envir = envir, inherits = TRUE)
    }
  }

  if (exists("tcga_auc_results", envir = envir, inherits = TRUE)) {
    tcga_auc_results <- get("tcga_auc_results", envir = envir, inherits = TRUE)
    for (method_name in names(object_map)) {
      if (method_name %in% names(tcga_auc_results)) {
        results[[method_name]] <- tcga_auc_results[[method_name]]
      }
    }
    if ("ElasticNetCox" %in% names(tcga_auc_results)) {
      results[["EN-Cox"]] <- tcga_auc_results[["ElasticNetCox"]]
    }
  }

  if (length(results) == 0L) {
    stop(
      "No cv-test result objects found. Run TCGA_HNSC.R first, or pass ",
      "`results = list(ContiCox = auc_curves_list_drpcapls, ...)`."
    )
  }

  results
}

make_test_risk_group_data <- function(result, time, status,
                                      time_scale = 365.25,
                                      high_risk_when = c("greater_equal", "greater")) {
  high_risk_when <- match.arg(high_risk_when)

  if (is.null(result$test_marker)) {
    stop("The result object does not contain `test_marker`; use a val_*_cv_test result.")
  }

  risk_score <- as.numeric(result$test_marker)
  test_idx <- result$test_idx

  if (!is.null(test_idx) && length(risk_score) == length(test_idx)) {
    time_test <- time[test_idx]
    status_test <- status[test_idx]
  } else if (length(risk_score) == length(time)) {
    time_test <- time
    status_test <- status
  } else {
    stop("Cannot align test_marker with time/status. Check result$test_idx.")
  }

  status_test <- as.integer(status_test != 0L)
  ok <- is.finite(time_test) & time_test > 0 &
    is.finite(status_test) & is.finite(risk_score)
  time_test <- as.numeric(time_test[ok])
  status_test <- status_test[ok]
  risk_score <- risk_score[ok]

  cutoff <- stats::median(risk_score, na.rm = TRUE)
  if (high_risk_when == "greater_equal") {
    group <- ifelse(risk_score >= cutoff, "high-risk patients", "low-risk patients")
  } else {
    group <- ifelse(risk_score > cutoff, "high-risk patients", "low-risk patients")
  }

  group <- factor(group, levels = c("low-risk patients", "high-risk patients"))
  if (length(unique(group[!is.na(group)])) < 2L) {
    stop("Median split produced only one group; check test risk scores.")
  }

  data.frame(
    time = time_test / time_scale,
    status = status_test,
    risk_score = risk_score,
    group = group,
    cutoff = cutoff
  )
}

logrank_p_value <- function(km_data) {
  .require_package("survival")
  fit_diff <- survival::survdiff(
    survival::Surv(time, status) ~ group,
    data = km_data
  )
  stats::pchisq(fit_diff$chisq, df = length(fit_diff$n) - 1L, lower.tail = FALSE)
}

format_logrank_p <- function(p_value, digits = 3) {
  if (!is.finite(p_value)) {
    return("p=NA")
  }
  if (p_value < 10 ^ (-digits)) {
    return(paste0("p<", formatC(10 ^ (-digits), format = "f", digits = digits)))
  }
  paste0("p=", signif(p_value, digits = digits))
}

plot_test_risk_km_panel <- function(result, time, status,
                                    method_label,
                                    panel_label = "(a)",
                                    time_scale = 365.25,
                                    xlim = NULL,
                                    ylim = c(0, 1),
                                    p_pos = NULL,
                                    high_risk_when = "greater_equal") {
  .require_package("survival")

  km_data <- make_test_risk_group_data(
    result = result,
    time = time,
    status = status,
    time_scale = time_scale,
    high_risk_when = high_risk_when
  )
  fit <- survival::survfit(survival::Surv(time, status) ~ group, data = km_data)
  p_value <- logrank_p_value(km_data)

  if (is.null(xlim)) {
    xlim <- c(0, max(km_data$time, na.rm = TRUE))
  }
  if (is.null(p_pos)) {
    p_pos <- c(xlim[1] + 0.43 * diff(xlim), ylim[1] + 0.68 * diff(ylim))
  }

  graphics::plot(
    fit,
    col = c("grey20", "grey20"),
    lty = c(1, 2),
    lwd = 1.15,
    mark.time = TRUE,
    xlim = xlim,
    ylim = ylim,
    xlab = "Survival in Years",
    ylab = "Survival Probability",
    conf.int = FALSE,
    main = "",
    bty = "o",
    las = 1
  )
  graphics::legend(
    "topright",
    legend = levels(km_data$group),
    col = c("grey20", "grey20"),
    lty = c(1, 2),
    lwd = 1.15,
    bty = "n",
    cex = 0.95
  )
  graphics::text(xlim[1] + 0.06 * diff(xlim), ylim[2] - 0.08 * diff(ylim),
                 labels = panel_label, font = 2, cex = 1.0)
  graphics::text(p_pos[1], p_pos[2], labels = format_logrank_p(p_value), cex = 0.95)
  graphics::mtext(method_label, side = 3, line = 0.2, cex = 0.85)

  invisible(list(data = km_data, fit = fit, p_value = p_value))
}

extract_test_auc_curve <- function(result) {
  if (!is.null(result$test_auc_curve)) {
    return(as.numeric(result$test_auc_curve))
  }
  if (!is.null(result$best_auc_curve)) {
    return(as.numeric(result$best_auc_curve))
  }
  stop("The result object does not contain test_auc_curve or best_auc_curve.")
}

plot_test_auc_panel <- function(results,
                                method_names = names(results),
                                panel_label = "(d)",
                                time_scale = 365.25,
                                xlim = NULL,
                                ylim = NULL,
                                ylab = "Area under the curve") {
  method_names <- method_names[method_names %in% names(results)]
  if (length(method_names) == 0L) {
    stop("No requested methods are available in `results`.")
  }

  auc_curves <- lapply(results[method_names], extract_test_auc_curve)
  time_grid <- results[[method_names[1L]]]$auc_time_grid
  if (is.null(time_grid)) {
    stop("Result object does not contain auc_time_grid.")
  }
  time_grid <- as.numeric(time_grid) / time_scale

  if (is.null(xlim)) {
    xlim <- range(time_grid, finite = TRUE)
  }
  if (is.null(ylim)) {
    auc_range <- range(unlist(auc_curves), finite = TRUE)
    ylim <- c(max(0.45, auc_range[1] - 0.03), min(1.00, auc_range[2] + 0.03))
  }

  line_types <- seq_len(length(method_names))
  line_cols <- grDevices::gray.colors(length(method_names), start = 0.15, end = 0.55)

  graphics::plot(
    time_grid,
    auc_curves[[1L]],
    type = "n",
    xlim = xlim,
    ylim = ylim,
    xlab = "time",
    ylab = ylab,
    main = "",
    bty = "o",
    las = 1
  )
  for (ii in seq_along(method_names)) {
    graphics::lines(
      time_grid,
      auc_curves[[ii]],
      col = line_cols[ii],
      lty = line_types[ii],
      lwd = 1.2
    )
  }
  graphics::legend(
    "topright",
    legend = method_names,
    col = line_cols,
    lty = line_types,
    lwd = 1.2,
    bty = "n",
    cex = 0.95
  )
  graphics::text(xlim[1] + 0.06 * diff(xlim), ylim[2] - 0.08 * diff(ylim),
                 labels = panel_label, font = 2, cex = 1.0)

  invisible(list(time = time_grid, auc_curves = auc_curves))
}

plot_tcga_test_median_risk_figure <- function(
    time = NULL,
    status = NULL,
    results = NULL,
    km_methods = c("PLSDR", "ContiCox", "pCox"),
    auc_methods = c("PLSDR", "ContiCox", "pCox"),
    time_scale = 365.25,
    output_base = default_output_base,
    save_pdf = TRUE,
    save_png = TRUE,
    png_res = 600,
    width = 7.2,
    height = 5.6,
    plot_to_current_device = interactive(),
    high_risk_when = "greater_equal"
) {
  caller <- parent.frame()
  if (is.null(time)) {
    time <- .first_existing_object(c("time", "HNSC_time", "GBM_time", "LUAD_time", "BRCA_time", "KIRC_time"), caller)
  }
  if (is.null(status)) {
    status <- .first_existing_object(c("status", "HNSC_status", "GBM_status", "LUAD_status", "BRCA_status", "KIRC_status"), caller)
  }
  if (is.null(time) || is.null(status)) {
    stop("Cannot find time/status. Pass them explicitly or run TCGA_HNSC.R first.")
  }
  if (is.null(results)) {
    results <- collect_tcga_cvtest_results(caller)
  }

  km_methods <- km_methods[km_methods %in% names(results)]
  auc_methods <- auc_methods[auc_methods %in% names(results)]
  if (length(km_methods) < 1L) {
    stop("No requested KM methods are available.")
  }
  if (length(auc_methods) < 1L) {
    auc_methods <- names(results)
  }

  plot_once <- function() {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mfrow = c(2, 2), mar = c(4.1, 4.3, 1.3, 1.0), oma = c(0, 0, 0, 0))

    panel_letters <- paste0("(", letters[seq_len(4L)], ")")
    km_outputs <- list()
    for (ii in seq_len(3L)) {
      if (ii <= length(km_methods)) {
        method_name <- km_methods[ii]
        km_outputs[[method_name]] <- plot_test_risk_km_panel(
          result = results[[method_name]],
          time = time,
          status = status,
          method_label = method_name,
          panel_label = panel_letters[ii],
          time_scale = time_scale,
          high_risk_when = high_risk_when
        )
      } else {
        graphics::plot.new()
      }
    }
    auc_output <- plot_test_auc_panel(
      results = results,
      method_names = auc_methods,
      panel_label = panel_letters[4L],
      time_scale = time_scale
    )
    invisible(list(km = km_outputs, auc = auc_output))
  }

  dir.create(dirname(output_base), recursive = TRUE, showWarnings = FALSE)

  plot_result <- NULL
  if (isTRUE(plot_to_current_device)) {
    plot_result <- plot_once()
  }
  if (isTRUE(save_pdf)) {
    grDevices::pdf(paste0(output_base, ".pdf"), width = width, height = height,
                   pointsize = 10, useDingbats = FALSE)
    plot_result <- plot_once()
    grDevices::dev.off()
  }
  if (isTRUE(save_png)) {
    grDevices::png(paste0(output_base, ".png"), width = width, height = height,
                   units = "in", res = png_res, pointsize = 10)
    plot_result <- plot_once()
    grDevices::dev.off()
  }
  if (is.null(plot_result)) {
    temp_pdf <- tempfile(fileext = ".pdf")
    grDevices::pdf(temp_pdf, width = width, height = height,
                   pointsize = 10, useDingbats = FALSE)
    plot_result <- plot_once()
    grDevices::dev.off()
    unlink(temp_pdf)
  }

  km_summary <- do.call(rbind, lapply(names(plot_result$km), function(method_name) {
    km_data <- plot_result$km[[method_name]]$data
    data.frame(
      Method = method_name,
      Test_N = nrow(km_data),
      Events = sum(km_data$status == 1L),
      Cutoff_median_test_risk = unique(km_data$cutoff)[1L],
      Low_risk_N = sum(km_data$group == "low-risk patients"),
      High_risk_N = sum(km_data$group == "high-risk patients"),
      Logrank_p = plot_result$km[[method_name]]$p_value,
      row.names = NULL
    )
  }))

  list(
    summary = km_summary,
    plot_result = plot_result,
    output_pdf = if (isTRUE(save_pdf)) paste0(output_base, ".pdf") else NA_character_,
    output_png = if (isTRUE(save_png)) paste0(output_base, ".png") else NA_character_
  )
}

run_tcga_test_median_risk_analysis <- function(
    X,
    time,
    status,
    output_base = default_output_base,
    save_fit_rds = TRUE,
    fit_rds = paste0(output_base, "_fit_results.rds"),
    km_methods = c("PLSDR", "ContiCox", "pCox"),
    auc_methods = c("PLSDR", "ContiCox", "pCox", "Ridge-Cox", "EN-Cox"),
    methods = unique(c(km_methods, auc_methods)),
    test_prop = 0.3,
    seed = 2026,
    k = 5,
    auc_method = c("IPCW", "NNE", "KM"),
    auc_time_grid = NULL,
    auc_time_probs = seq(0.2, 0.8, length.out = 30),
    ncomp_candidates_DRPLS = 1:10,
    ncomp_candidates_pcox = 2:10,
    alpha_candidates = seq(0, 1, by = 0.1),
    lambda_candidates_ridgecox = 10 ^ seq(-4, 1, by = 1),
    alpha_candidates_enetcox = c(0.25, 0.5, 0.75),
    lambda_candidates_enetcox = 10 ^ seq(-5, 1, by = 1),
    use_preselection = TRUE,
    p_thresh = 0.05,
    residual_type = "martingale",
    normalize_alpha_terms = TRUE,
    quiet = FALSE,
    save_pdf = TRUE,
    save_png = TRUE,
    png_res = 600,
    width = 7.2,
    height = 5.6,
    plot_to_current_device = interactive(),
    high_risk_when = "greater_equal"
) {
  auc_method <- match.arg(auc_method)

  fit <- fit_tcga_cvtest_results(
    X = X,
    time = time,
    status = status,
    test_prop = test_prop,
    seed = seed,
    k = k,
    auc_method = auc_method,
    auc_time_grid = auc_time_grid,
    auc_time_probs = auc_time_probs,
    ncomp_candidates_DRPLS = ncomp_candidates_DRPLS,
    ncomp_candidates_pcox = ncomp_candidates_pcox,
    alpha_candidates = alpha_candidates,
    lambda_candidates_ridgecox = lambda_candidates_ridgecox,
    alpha_candidates_enetcox = alpha_candidates_enetcox,
    lambda_candidates_enetcox = lambda_candidates_enetcox,
    use_preselection = use_preselection,
    p_thresh = p_thresh,
    residual_type = residual_type,
    normalize_alpha_terms = normalize_alpha_terms,
    methods = methods,
    quiet = quiet
  )

  plot <- plot_tcga_test_median_risk_figure(
    time = fit$time,
    status = fit$status,
    results = fit$results,
    km_methods = km_methods,
    auc_methods = auc_methods,
    output_base = output_base,
    save_pdf = save_pdf,
    save_png = save_png,
    png_res = png_res,
    width = width,
    height = height,
    plot_to_current_device = plot_to_current_device,
    high_risk_when = high_risk_when
  )

  out <- list(fit = fit, plot = plot)

  if (isTRUE(save_fit_rds)) {
    dir.create(dirname(fit_rds), recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, fit_rds)
    out$fit_rds <- fit_rds
  }

  out
}

if (identical(Sys.getenv("RUN_TCGA_TEST_MEDIAN_RISK_KM", unset = "1"), "1")) {
  has_xy <- exists("X", inherits = TRUE) &&
    exists("time", inherits = TRUE) &&
    exists("status", inherits = TRUE)
  has_existing_results <- exists("time", inherits = TRUE) && exists("status", inherits = TRUE) &&
    any(vapply(
      c("auc_curves_list_drpls", "auc_curves_list_drpcapls", "auc_curves_list_pcox",
        "auc_curves_list_ridgecox", "auc_curves_list_enetcox", "tcga_auc_results"),
      exists,
      logical(1),
      inherits = TRUE
    ))
  force_refit <- identical(Sys.getenv("FORCE_REFIT_TCGA_TEST_MEDIAN_RISK_KM", unset = "0"), "1")

  if (has_xy && (force_refit || !has_existing_results)) {
    tcga_test_median_risk_analysis <- run_tcga_test_median_risk_analysis(
      X = X,
      time = time,
      status = status
    )
    print(tcga_test_median_risk_analysis$plot$summary)
    message("Saved: ", tcga_test_median_risk_analysis$plot$output_pdf)
    message("Saved: ", tcga_test_median_risk_analysis$plot$output_png)
    message("Saved fit/results: ", tcga_test_median_risk_analysis$fit_rds)
  } else if (has_existing_results) {
    tcga_test_median_risk_plot <- plot_tcga_test_median_risk_figure()
    print(tcga_test_median_risk_plot$summary)
    message("Saved: ", tcga_test_median_risk_plot$output_pdf)
    message("Saved: ", tcga_test_median_risk_plot$output_png)
  } else {
    message(
      "Plot functions loaded. Create X, time, and status first, then call ",
      "run_tcga_test_median_risk_analysis(X, time, status)."
    )
  }
}

Sys.setenv(RUN_TCGA_TEST_MEDIAN_RISK_KM = "0")

X <- HNSC_X
time <- HNSC_time
status <- HNSC_status

out <- run_tcga_test_median_risk_analysis(
  X = X,
  time = time,
  status = status,
  
  methods = c("PLSDR", "ContiCox", "pCox", "Ridge-Cox", "EN-Cox"),
  km_methods = c("PLSDR", "ContiCox", "pCox"),
  auc_methods = c("PLSDR", "ContiCox", "pCox", "Ridge-Cox", "EN-Cox"),
  
  test_prop = 0.3,
  seed = 2026,
  k = 5,
  
  ncomp_candidates_DRPLS = 1:10,
  ncomp_candidates_pcox = 2:10,
  alpha_candidates = seq(0, 1, by = 0.1),
  lambda_candidates_ridgecox = 10 ^ seq(-4, 1, by = 1),
  alpha_candidates_enetcox = c(0.25, 0.5, 0.75),
  lambda_candidates_enetcox = 10 ^ seq(-5, 1, by = 1),
  
  use_preselection = TRUE,
  p_thresh = 0.05,
  
  residual_type = "martingale",
  normalize_alpha_terms = TRUE,
  
  auc_method = "NNE",
  
  output_base = "figures/TCGA_HNSC_test_median_risk_refit"
)

out$plot$summary