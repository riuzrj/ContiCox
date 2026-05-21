library(conticox)

set.seed(123)
n <- 80
p <- 20
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
lp <- 0.8 * X[, 1] - 0.6 * X[, 2] + 0.3 * X[, 3]
event_time <- rexp(n, rate = 0.08 * exp(lp))
censor_time <- rexp(n, rate = 0.04)
time <- pmin(event_time, censor_time)
status <- as.integer(event_time <= censor_time)

res <- compare_survival_methods(
  X = X,
  time = time,
  status = status,
  methods = c("ContiCox", "DRPLS", "partialCox", "RidgeCox", "ElasticNetCox"),
  test_prop = 0.3,
  k = 2,
  ncomp_candidates = 1:2,
  conticox_alpha_candidates = c(0, 0.5, 1),
  ridge_lambda_candidates = 10 ^ c(-3, -1),
  enet_alpha_candidates = c(0.25, 0.5),
  enet_lambda_candidates = 10 ^ c(-3, -1),
  auc_method = "NNE",
  seed = 123
)

print(res$summary)
