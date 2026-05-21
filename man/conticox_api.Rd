\name{conticox}
\alias{conticox}
\alias{auc_curve}
\alias{compare_survival_methods}
\alias{conticox_cv}
\alias{conticox_cv_test}
\alias{conticox_fit}
\alias{enetcox_cv}
\alias{enetcox_cv_test}
\alias{integrated_auc}
\alias{partialcox_cv}
\alias{partialcox_cv_test}
\alias{partialcox_fit}
\alias{penalizedcox_cv}
\alias{penalizedcox_cv_test}
\alias{plsdr_cv}
\alias{plsdr_cv_test}
\alias{plsdr_fit}
\alias{ridgecox_cv}
\alias{ridgecox_cv_test}
\title{ContiCox survival prediction methods}
\description{
User-facing wrappers for the ContiCox continuum PCA-PLS Cox model and
competing survival prediction baselines used in the paper.
}
\usage{
conticox_fit(X, time, status, n_components = min(7, ncol(as.matrix(X))),
  alpha = 0.5, residual_type = "deviance",
  normalize_alpha_terms = TRUE, ...)

conticox_cv(..., auc_method = c("IPCW", "NNE", "KM"))

conticox_cv_test(..., auc_method = c("IPCW", "NNE", "KM"))

plsdr_fit(X, time, status, n_components = min(7, ncol(as.matrix(X))), ...)

plsdr_cv(..., auc_method = c("IPCW", "NNE", "KM"))

plsdr_cv_test(..., auc_method = c("IPCW", "NNE", "KM"))

partialcox_fit(X, time, status,
  n_components = min(7, ncol(as.matrix(X))), ...)

partialcox_cv(..., auc_method = c("IPCW", "NNE", "KM"))

partialcox_cv_test(..., auc_method = c("IPCW", "NNE", "KM"))

ridgecox_cv(...)

ridgecox_cv_test(...)

enetcox_cv(...)

enetcox_cv_test(...)

penalizedcox_cv(..., auc_method = c("IPCW", "NNE", "KM"))

penalizedcox_cv_test(..., auc_method = c("IPCW", "NNE", "KM"))

auc_curve(time, status, marker, auc_time_grid,
  auc_method = c("IPCW", "NNE", "KM"),
  ipcw_weighting = "marginal", nne_span = 0.30)

integrated_auc(auc_values)

compare_survival_methods(X, time, status,
  methods = c("ContiCox", "DRPLS", "partialCox", "RidgeCox", "ElasticNetCox"),
  test_prop = 0.3, train_idx = NULL, test_idx = NULL, k = 5,
  ncomp_candidates = 1:5,
  conticox_alpha_candidates = seq(0, 1, by = 0.1),
  ridge_lambda_candidates = 10^seq(-4, 1, by = 1),
  enet_alpha_candidates = c(0.25, 0.5, 0.75),
  enet_lambda_candidates = 10^seq(-4, 1, by = 1),
  auc_time_grid = NULL, use_preselection = FALSE, p_thresh = 0.05,
  max_features = NULL, auc_method = c("IPCW", "NNE", "KM"), seed = 123,
  residual_type = "deviance", normalize_alpha_terms = TRUE,
  standardize = TRUE, stop_on_error = FALSE, verbose = TRUE)
}
\arguments{
\item{X}{Numeric predictor matrix with rows as subjects and columns as predictors.}
\item{time}{Observed survival or follow-up time.}
\item{status}{Event indicator, coded as 1 for event and 0 for censoring.}
\item{n_components, ncomp_candidates}{Number of latent components or candidate component grid.}
\item{alpha, conticox_alpha_candidates}{Continuum parameter for the ContiCox method.}
\item{auc_method}{Time-dependent AUC estimator: \code{"IPCW"}, \code{"NNE"}, or \code{"KM"}.}
\item{marker}{Risk score or linear predictor used for AUC calculation.}
\item{auc_time_grid}{Evaluation time grid.}
\item{auc_values}{Vector of time-dependent AUC values.}
\item{methods}{Methods to compare.}
\item{test_prop, train_idx, test_idx}{Train/test split controls.}
\item{k}{Number of cross-validation folds.}
\item{p_thresh, use_preselection, max_features}{Optional feature screening controls.}
\item{ridge_lambda_candidates, enet_alpha_candidates, enet_lambda_candidates}{Penalized Cox tuning grids.}
\item{residual_type, normalize_alpha_terms, standardize}{Model fitting options.}
\item{seed}{Random seed for reproducible splits.}
\item{stop_on_error}{If \code{TRUE}, stop when any method fails.}
\item{verbose}{If \code{TRUE}, print method progress messages.}
\item{ipcw_weighting, nne_span}{Options passed to the AUC estimators.}
\item{...}{Additional arguments passed to the underlying method functions.}
}
\value{
Fit wrappers return fitted model objects. CV and CV/test wrappers return lists
with selected tuning parameters, AUC curves, and iAUC values.
\code{compare_survival_methods()} returns a list with \code{summary},
\code{results}, train/test indices, the AUC time grid, and the AUC method.
}
\examples{
set.seed(1)
n <- 40
p <- 8
X <- matrix(rnorm(n * p), n, p)
lp <- 0.5 * X[, 1] - 0.4 * X[, 2]
event_time <- rexp(n, rate = 0.10 * exp(lp))
censor_time <- rexp(n, rate = 0.04)
time <- pmin(event_time, censor_time)
status <- as.integer(event_time <= censor_time)

grid <- as.numeric(quantile(time, probs = c(0.3, 0.6)))
res <- compare_survival_methods(
  X, time, status,
  methods = c("ContiCox", "DRPLS"),
  k = 2,
  ncomp_candidates = 1:2,
  conticox_alpha_candidates = c(0, 0.5),
  auc_time_grid = grid,
  auc_method = "NNE",
  verbose = FALSE
)
res$summary
}
