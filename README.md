# conticox

`conticox` is a small R package that collects the survival prediction methods
used in the current paper code:

- `ContiCox`: continuum PCA-PLS Cox model
- `DRPLS`: PLS dimension reduction Cox baseline
- `partialCox`: partial Cox component baseline
- `RidgeCox`: ridge-penalized Cox baseline
- `ElasticNetCox`: elastic-net Cox baseline

The package keeps the original method implementations and adds a simpler API
for fitting, cross-validation, train/test validation, and time-dependent AUC
calculation.

## Install locally

From the project root:

```r
install.packages("conticoxpkg", repos = NULL, type = "source")
```

or from the command line:

```sh
R CMD INSTALL conticoxpkg
```

## Basic use

```r
library(conticox)

# X: n by p predictor matrix
# time: observed survival/follow-up time
# status: event indicator, 1 = event, 0 = censored

res <- compare_survival_methods(
  X = X,
  time = time,
  status = status,
  methods = c("ContiCox", "DRPLS", "partialCox", "RidgeCox", "ElasticNetCox"),
  test_prop = 0.3,
  k = 5,
  ncomp_candidates = 1:10,
  conticox_alpha_candidates = seq(0, 1, by = 0.1),
  auc_method = "NNE",
  seed = 123
)

res$summary
```

Use `auc_method = "IPCW"` when the paper analysis should use the IPCW
time-dependent AUC estimator. Use `auc_method = "NNE"` to reproduce the earlier
nearest-neighbor estimator based on `survivalROC`.

## Single-method wrappers

```r
fit <- conticox_fit(X, time, status, n_components = 3, alpha = 0.6)

cv_test <- conticox_cv_test(
  X = X,
  time = time,
  status = status,
  test_prop = 0.3,
  k = 5,
  ncomp_candidates = 1:10,
  alpha_candidates = seq(0, 1, by = 0.1),
  auc_method = "IPCW"
)

cv_test$test_iAUC
```

Available train/test wrappers:

- `conticox_cv_test()`
- `plsdr_cv_test()`
- `partialcox_cv_test()`
- `ridgecox_cv_test()`
- `enetcox_cv_test()`

## AUC utilities

```r
grid <- as.numeric(quantile(time, probs = seq(0.2, 0.8, length.out = 10)))
auc <- auc_curve(time, status, marker = risk_score, auc_time_grid = grid,
                 auc_method = "IPCW")
integrated_auc(auc)
```

## Paper scripts

Simulation and real-data scripts used for the paper are bundled as supporting
materials under:

```r
system.file("docs", package = "conticox")
```

They are stored in `inst/docs/simulation` and `inst/docs/real_data` in the
source package. These files are not loaded as package code.
