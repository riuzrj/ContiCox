library(conticox)

set.seed(1)
n <- 60
p <- 12
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
lp <- 0.7 * X[, 1] - 0.5 * X[, 2]
event_time <- rexp(n, rate = 0.10 * exp(lp))
censor_time <- rexp(n, rate = 0.04)
time <- pmin(event_time, censor_time)
status <- as.integer(event_time <= censor_time)

grid <- as.numeric(stats::quantile(time, probs = c(0.25, 0.50, 0.75)))

cmp <- compare_survival_methods(
  X = X,
  time = time,
  status = status,
  methods = c("ContiCox", "DRPLS", "partialCox", "RidgeCox", "ElasticNetCox"),
  test_prop = 0.3,
  k = 2,
  ncomp_candidates = 1:2,
  conticox_alpha_candidates = c(0, 0.5),
  ridge_lambda_candidates = 10 ^ c(-3, -1),
  enet_alpha_candidates = c(0.25, 0.5),
  enet_lambda_candidates = 10 ^ c(-3, -1),
  auc_time_grid = grid,
  use_preselection = FALSE,
  auc_method = "NNE",
  seed = 1,
  stop_on_error = TRUE,
  verbose = FALSE
)

stopifnot(is.data.frame(cmp$summary))
stopifnot(nrow(cmp$summary) == 5L)
stopifnot(all(is.finite(cmp$summary$test_iAUC)))

auc <- auc_curve(time, status, marker = lp, auc_time_grid = grid, auc_method = "NNE")
stopifnot(length(auc) == length(grid))
stopifnot(is.finite(integrated_auc(auc)))
