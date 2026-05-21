project_dir <- "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls"
setwd(project_dir)

script_lines <- readLines(file.path(project_dir, "simulation", "search_target_gap_settings.R"), warn = FALSE)
settings_start <- grep("^settings <- list", script_lines)[1L]
eval(parse(text = script_lines[seq_len(settings_start - 1L)]), envir = .GlobalEnv)

simulate_variance_signal <- function(
    n_samples = 80,
    n_total_genes = 2000,
    signal_genes = 80,
    signal_scale = 3.0,
    signal_rho = 0.9,
    decoy_blocks = 8,
    decoy_genes = 40,
    decoy_scale = 1.0,
    decoy_rho = 0.7,
    signal_strength = 0.8,
    base_rate = 0.1,
    censoring = 0.45
) {
  stopifnot(signal_genes + decoy_blocks * decoy_genes <= n_total_genes)

  X <- matrix(rnorm(n_samples * n_total_genes), nrow = n_samples)
  colnames(X) <- paste0("Gene_", seq_len(n_total_genes))

  z_signal <- rnorm(n_samples)
  current_col <- 1L
  for (j in seq_len(signal_genes)) {
    X[, current_col] <- signal_scale * (
      sqrt(signal_rho) * z_signal + sqrt(1 - signal_rho) * rnorm(n_samples)
    )
    current_col <- current_col + 1L
  }

  for (block in seq_len(decoy_blocks)) {
    z_decoy <- rnorm(n_samples)
    for (j in seq_len(decoy_genes)) {
      X[, current_col] <- decoy_scale * (
        sqrt(decoy_rho) * z_decoy + sqrt(1 - decoy_rho) * rnorm(n_samples)
      )
      current_col <- current_col + 1L
    }
  }

  X <- as.matrix(scale(X, center = TRUE, scale = FALSE))
  lp <- signal_strength * as.numeric(scale(z_signal))

  event_time <- -log(runif(n_samples)) / (base_rate * exp(lp))
  censor_rate <- estimate_exponential_censor_rate(event_time, censoring)
  censor_time <- rexp(n_samples, rate = censor_rate)

  list(
    X = X,
    surv = data.frame(
      Time = pmin(event_time, censor_time),
      Status = as.integer(event_time <= censor_time)
    ),
    latent = data.frame(signal = z_signal),
    linear_predictor = lp,
    params = list(
      scenario = "variance_signal_drpls_gap",
      n_samples = n_samples,
      n_total_genes = n_total_genes,
      signal_genes = signal_genes,
      signal_scale = signal_scale,
      signal_rho = signal_rho,
      decoy_blocks = decoy_blocks,
      decoy_genes = decoy_genes,
      decoy_scale = decoy_scale,
      decoy_rho = decoy_rho,
      signal_strength = signal_strength,
      censoring = censoring
    )
  )
}

simulate_variance_mixture <- function(
    n_samples = 80,
    n_total_genes = 2000,
    signal_genes = 60,
    supervised_genes = 10,
    signal_scale = 3.0,
    supervised_scale = 1.0,
    signal_rho = 0.9,
    supervised_rho = 0.35,
    signal_weight = 0.85,
    supervised_weight = 0.35,
    decoy_blocks = 10,
    decoy_genes = 40,
    decoy_scale = 1.0,
    decoy_rho = 0.8,
    signal_strength = 0.8,
    base_rate = 0.1,
    censoring = 0.45
) {
  stopifnot(signal_genes + supervised_genes + decoy_blocks * decoy_genes <= n_total_genes)

  X <- matrix(rnorm(n_samples * n_total_genes), nrow = n_samples)
  colnames(X) <- paste0("Gene_", seq_len(n_total_genes))

  z_signal <- rnorm(n_samples)
  z_supervised <- rnorm(n_samples)
  current_col <- 1L
  for (j in seq_len(signal_genes)) {
    X[, current_col] <- signal_scale * (
      sqrt(signal_rho) * z_signal + sqrt(1 - signal_rho) * rnorm(n_samples)
    )
    current_col <- current_col + 1L
  }
  for (j in seq_len(supervised_genes)) {
    X[, current_col] <- supervised_scale * (
      sqrt(supervised_rho) * z_supervised +
        sqrt(1 - supervised_rho) * rnorm(n_samples)
    )
    current_col <- current_col + 1L
  }
  for (block in seq_len(decoy_blocks)) {
    z_decoy <- rnorm(n_samples)
    for (j in seq_len(decoy_genes)) {
      X[, current_col] <- decoy_scale * (
        sqrt(decoy_rho) * z_decoy + sqrt(1 - decoy_rho) * rnorm(n_samples)
      )
      current_col <- current_col + 1L
    }
  }

  X <- as.matrix(scale(X, center = TRUE, scale = FALSE))
  raw_score <- signal_weight * as.numeric(scale(z_signal)) +
    supervised_weight * as.numeric(scale(z_supervised))
  lp <- signal_strength * as.numeric(scale(raw_score))

  event_time <- -log(runif(n_samples)) / (base_rate * exp(lp))
  censor_rate <- estimate_exponential_censor_rate(event_time, censoring)
  censor_time <- rexp(n_samples, rate = censor_rate)

  list(
    X = X,
    surv = data.frame(
      Time = pmin(event_time, censor_time),
      Status = as.integer(event_time <= censor_time)
    ),
    latent = data.frame(signal = z_signal, supervised = z_supervised),
    linear_predictor = lp,
    params = list(
      scenario = "variance_mixture_drpls_gap",
      n_samples = n_samples,
      n_total_genes = n_total_genes,
      signal_genes = signal_genes,
      supervised_genes = supervised_genes,
      signal_scale = signal_scale,
      supervised_scale = supervised_scale,
      signal_rho = signal_rho,
      supervised_rho = supervised_rho,
      signal_weight = signal_weight,
      supervised_weight = supervised_weight,
      decoy_blocks = decoy_blocks,
      decoy_genes = decoy_genes,
      decoy_scale = decoy_scale,
      decoy_rho = decoy_rho,
      signal_strength = signal_strength,
      censoring = censoring
    )
  )
}

settings <- list(
  make_setting(
    "var_signal_n60_p2000_scale3_g07_c50",
    "simulate_variance_signal",
    list(n_samples = 60, n_total_genes = 2000,
         signal_genes = 80, signal_scale = 3.0,
         signal_rho = 0.9, signal_strength = 0.7,
         censoring = 0.50),
    301L
  ),
  make_setting(
    "var_signal_n70_p2000_scale4_g07_c50",
    "simulate_variance_signal",
    list(n_samples = 70, n_total_genes = 2000,
         signal_genes = 80, signal_scale = 4.0,
         signal_rho = 0.92, signal_strength = 0.7,
         censoring = 0.50),
    302L
  ),
  make_setting(
    "var_signal_n80_p2000_scale3_g08_c45",
    "simulate_variance_signal",
    list(n_samples = 80, n_total_genes = 2000,
         signal_genes = 60, signal_scale = 3.0,
         signal_rho = 0.9, signal_strength = 0.8,
         censoring = 0.45),
    303L
  ),
  make_setting(
    "var_signal_n80_p3000_scale4_g08_c45",
    "simulate_variance_signal",
    list(n_samples = 80, n_total_genes = 3000,
         signal_genes = 80, signal_scale = 4.0,
         signal_rho = 0.92, signal_strength = 0.8,
         decoy_blocks = 12, decoy_genes = 50,
         censoring = 0.45),
    304L
  ),
  make_setting(
    "var_mix_n70_p2000_scale4_g08_c50",
    "simulate_variance_mixture",
    list(n_samples = 70, n_total_genes = 2000,
         signal_genes = 80, supervised_genes = 8,
         signal_scale = 4.0, supervised_scale = 1.0,
         signal_rho = 0.92, supervised_rho = 0.25,
         signal_weight = 0.85, supervised_weight = 0.35,
         signal_strength = 0.8, censoring = 0.50),
    305L
  ),
  make_setting(
    "var_mix_n80_p2500_scale5_g08_c45",
    "simulate_variance_mixture",
    list(n_samples = 80, n_total_genes = 2500,
         signal_genes = 80, supervised_genes = 10,
         signal_scale = 5.0, supervised_scale = 1.0,
         signal_rho = 0.94, supervised_rho = 0.30,
         signal_weight = 0.80, supervised_weight = 0.45,
         decoy_blocks = 12, decoy_genes = 45,
         signal_strength = 0.8, censoring = 0.45),
    306L
  )
)

summaries <- data.frame()
for (setting in settings) {
  cat("\n--- DRPLS-gap coarse screening:", setting$name, "---\n")
  flush.console()
  ans <- run_core_setting(setting, R = 5L, n_cores = 5L, seed = 15100)
  summaries <- rbind(summaries, ans)
  print(ans, digits = 4)
  flush.console()
}

cat("\n=== DRPLS-gap coarse screening summary ===\n")
summaries <- summaries[order(-summaries$diff_vs_DRPLS), ]
print(summaries, digits = 4)

best <- settings[[match(summaries$setting[1L], vapply(settings, `[[`, "", "name"))]]
cat("\n--- Confirming best DRPLS-gap setting with R = 10:", best$name, "---\n")
confirm <- run_core_setting(best, R = 10L, n_cores = 10L, seed = 15900)
print(confirm, digits = 4)

