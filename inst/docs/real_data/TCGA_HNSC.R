
library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)
library(matrixStats)
library(dplyr)
library(survival)

## =========================
## 1) 下载 RNA-seq 计数 (STAR - Counts)
## =========================

options(timeout = 6000)                        # 拉长超时
dir.create("GDCdata", showWarnings = FALSE)


# HNSC
qry_expr <- GDCquery(
  project       = "TCGA-HNSC",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(qry_expr)
se <- GDCprepare(qry_expr, summarizedExperiment = TRUE)

counts <- SummarizedExperiment::assay(se)            # 行=基因(ENSEMBL)，列=样本条形码
barcodes <- colnames(counts)

## =========================
## 2) 筛选原发肿瘤样本 (01)
## =========================
sample_type <- substr(barcodes, 14, 15)              # "01"=Primary solid Tumor
keep_col <- sample_type == "01"
counts   <- counts[, keep_col, drop = FALSE]
barcodes <- barcodes[keep_col]

## =========================
## 3) 归一化：edgeR TMM → logCPM
## =========================
dge <- DGEList(counts = counts)

# 过滤低表达：至少在 20% 样本 CPM > 1（可按需调整）
cpm0 <- edgeR::cpm(dge)                               # 非 log，用于过滤判断
keep_gene <- rowMeans(cpm0 > 1) >= 0.20
dge <- dge[keep_gene, , keep.lib.sizes = FALSE]

# 重新计算 TMM 因子
dge <- calcNormFactors(dge, method = "TMM")

# 得到 log2-CPM（加 prior.count 避免 log(0)）
logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 1) # 基因×样本
expr_mat <- t(logCPM)                                  # 样本×基因 (n×p)

## =========================
## 4) 生存信息（临床）并与样本对齐
## =========================
## =========================
## 4) 生存信息：直接用 se 的 colData（与样本列一一对应）
## =========================
meta <- as.data.frame(SummarizedExperiment::colData(se))

# 只保留之前选过的原发肿瘤样本（keep_col）
meta <- meta[keep_col, , drop = FALSE]        # 与 counts[, keep_col] 同步

# 小工具：从多个候选字段里取第一个存在的列
pick <- function(df, candidates) {
  for (nm in candidates) if (nm %in% names(df)) return(df[[nm]])
  return(rep(NA, nrow(df)))
}

d2d <- suppressWarnings(as.numeric(pick(meta, c("days_to_death","days_to_death.x","paper_days_to_death"))))
d2f <- suppressWarnings(as.numeric(pick(meta, c("days_to_last_follow_up","days_to_last_followup",
                                                "days_to_last_follow_up.x","paper_days_to_last_followup"))))
vs  <- tolower(as.character(pick(meta, c("vital_status","vital_status.x","paper_vital_status"))))

HNSC_time   <- ifelse(!is.na(d2d), d2d, d2f)            # 生存时间（天）
HNSC_status <- ifelse(!is.na(d2d) & d2d > 0, 1,         # 事件：死亡=1；否则=0
                 ifelse(vs %in% c("dead","deceased"), 1, 0))


# 选高变基因并标准化，得到 X
madv <- matrixStats::colMads(expr_mat, na.rm = TRUE)
sel <- which(is.finite(madv) & madv > 0)
HNSC_X0 <- scale(expr_mat[, sel, drop = FALSE])  # 样本×基因

# 可选：去掉任何含 NA 的样本并同步 time/status
HNSC_ok <- is.finite(HNSC_time) &
  is.finite(HNSC_status) &
  HNSC_time > 0 &
  rowSums(is.finite(HNSC_X0)) == ncol(HNSC_X0)
HNSC_X <- HNSC_X0[HNSC_ok, , drop = FALSE]
HNSC_time <- HNSC_time[HNSC_ok]
HNSC_status <- HNSC_status[HNSC_ok]

# 快速自检（行数一致、顺序一致）
stopifnot(nrow(HNSC_X) == length(HNSC_time), length(HNSC_time) == length(HNSC_status))
nrow(HNSC_X)  # 样本数
ncol(HNSC_X)  # 特征数

cat("Number of patients:", nrow(HNSC_X), "\n")
cat("Number of events:", sum(HNSC_status == 1), "\n")
cat("Censoring rate:", mean(HNSC_status == 0), "\n")
cat("Final number of genes:", ncol(HNSC_X), "\n")
cat("Median follow-up:", median(HNSC_time, na.rm = TRUE), "\n")
cat("Follow-up range:", range(HNSC_time, na.rm = TRUE), "\n")
hist(HNSC_time[HNSC_status == 1],
     breaks = 30,
     main = "Distribution of Event Times",
     xlab = "Event time (days)",
     ylab = "Number of events",
     col = "gray80",
     border = "white")
sum(HNSC_status == 1 & HNSC_time >= 100 & HNSC_time <= 1000, na.rm = TRUE)



# GBM
qry_expr <- GDCquery(
  project       = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(qry_expr)
se <- GDCprepare(qry_expr, summarizedExperiment = TRUE)

counts <- SummarizedExperiment::assay(se)            # 行=基因(ENSEMBL)，列=样本条形码
barcodes <- colnames(counts)

sample_type <- substr(barcodes, 14, 15)
keep_col <- sample_type == "01"

counts <- counts[, keep_col, drop = FALSE]
barcodes <- barcodes[keep_col]


# 归一化
dge <- DGEList(counts = counts)

# 过滤低表达：至少在 20% 样本 CPM > 1（可按需调整）
cpm0 <- edgeR::cpm(dge)                               # 非 log，用于过滤判断
keep_gene <- rowMeans(cpm0 > 1) >= 0.20
dge <- dge[keep_gene, , keep.lib.sizes = FALSE]

# 重新计算 TMM 因子
dge <- calcNormFactors(dge, method = "TMM")

# 得到 log2-CPM（加 prior.count 避免 log(0)）
logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 1) # 基因×样本
expr_mat <- t(logCPM)     

# 生存信息
meta <- as.data.frame(SummarizedExperiment::colData(se))
meta <- meta[keep_col, , drop = FALSE]

pick <- function(df, candidates) {
  for (nm in candidates) if (nm %in% names(df)) return(df[[nm]])
  return(rep(NA, nrow(df)))
}

d2d <- as.numeric(pick(meta, c("days_to_death", "days_to_death.x", "paper_days_to_death")))
d2f <- as.numeric(pick(meta, c("days_to_last_follow_up", "days_to_last_followup",
                               "days_to_last_follow_up.x", "paper_days_to_last_followup")))
vs  <- tolower(as.character(pick(meta, c("vital_status","vital_status.x","paper_vital_status"))))

GBM_time <- ifelse(!is.na(d2d), d2d, d2f)
GBM_status <- ifelse(!is.na(d2d) & d2d > 0, 1,
                     ifelse(vs %in% c("dead","deceased"), 1, 0))

madv <- matrixStats::colMads(expr_mat, na.rm = TRUE)
sel <- which(is.finite(madv) & madv > 0)
GBM_X0 <- scale(expr_mat[, sel, drop = FALSE])

GBM_ok <- is.finite(GBM_time) &
  is.finite(GBM_status) &
  GBM_time > 0 &
  rowSums(is.finite(GBM_X0)) == ncol(GBM_X0)
GBM_X <- GBM_X0[GBM_ok, , drop = FALSE]
GBM_time <- GBM_time[GBM_ok]
GBM_status <- GBM_status[GBM_ok]

stopifnot(nrow(GBM_X) == length(GBM_time), length(GBM_time) == length(GBM_status))

cat("GBM number of patients:", nrow(GBM_X), "\n")
cat("GBM number of events:", sum(GBM_status == 1), "\n")
cat("GBM censoring rate:", mean(GBM_status == 0), "\n")
cat("GBM final number of genes:", ncol(GBM_X), "\n")
cat("GBM median follow-up:", median(GBM_time, na.rm = TRUE), "\n")
cat("GBM follow-up range:", range(GBM_time, na.rm = TRUE), "\n")


#LUAD
# ==== TCGA-LUAD ====
qry_expr <- GDCquery(
  project       = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(qry_expr, files.per.chunk = 40)
se <- GDCprepare(qry_expr, summarizedExperiment = TRUE)

counts   <- SummarizedExperiment::assay(se)
barcodes <- colnames(counts)

sample_type <- substr(barcodes, 14, 15)
keep_col    <- sample_type == "01"
counts      <- counts[, keep_col, drop = FALSE]
barcodes    <- barcodes[keep_col]

dge  <- DGEList(counts = counts)
cpm0 <- edgeR::cpm(dge)
keep_gene <- rowMeans(cpm0 > 1) >= 0.20
dge  <- dge[keep_gene, , keep.lib.sizes = FALSE]
dge  <- calcNormFactors(dge, method = "TMM")
logCPM   <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
expr_mat <- t(logCPM)

meta <- as.data.frame(SummarizedExperiment::colData(se))
meta <- meta[keep_col, , drop = FALSE]

d2d <- suppressWarnings(as.numeric(pick(meta, c("days_to_death","days_to_death.x","paper_days_to_death"))))
d2f <- suppressWarnings(as.numeric(pick(meta, c("days_to_last_follow_up","days_to_lastfollowup",
                                                "days_to_last_follow_up.x","paper_days_to_last_followup"))))
vs  <- tolower(as.character(pick(meta, c("vital_status","vital_status.x","paper_vital_status"))))

LUAD_time <- ifelse(!is.na(d2d), d2d, d2f)
LUAD_status <- ifelse(!is.na(d2d) & d2d > 0, 1,
                      ifelse(vs %in% c("dead","deceased"), 1, 0))

madv <- matrixStats::colMads(expr_mat, na.rm = TRUE)
sel <- which(is.finite(madv) & madv > 0)
LUAD_X0 <- scale(expr_mat[, sel, drop = FALSE])


LUAD_ok <- is.finite(LUAD_time) &
  is.finite(LUAD_status) &
  LUAD_time > 0 &
  rowSums(is.finite(LUAD_X0)) == ncol(LUAD_X0)
LUAD_X <- LUAD_X0[LUAD_ok, , drop = FALSE]
LUAD_time <- LUAD_time[LUAD_ok]
LUAD_status <- LUAD_status[LUAD_ok]

stopifnot(nrow(LUAD_X) == length(LUAD_time), length(LUAD_time) == length(LUAD_status))

cat("LUAD number of patients:", nrow(LUAD_X), "\n")
cat("LUAD number of events:", sum(LUAD_status == 1), "\n")
cat("LUAD censoring rate:", mean(LUAD_status == 0), "\n")
cat("LUAD final number of genes:", ncol(LUAD_X), "\n")
cat("LUAD median follow-up:", median(LUAD_time, na.rm = TRUE), "\n")
cat("LUAD follow-up range:", range(LUAD_time, na.rm = TRUE), "\n")
hist(LUAD_time[LUAD_status == 1],
     breaks = 30,
     main = "Distribution of Event Times",
     xlab = "Event time (days)",
     ylab = "Number of events",
     col = "gray80",
     border = "white")
sum(LUAD_status == 1 & LUAD_time >= 100 & LUAD_time <= 1000, na.rm = TRUE)


# BRCA
# ==== TCGA-BRCA ====
qry_expr <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(qry_expr)
se <- GDCprepare(qry_expr, summarizedExperiment = TRUE)

counts   <- SummarizedExperiment::assay(se)
barcodes <- colnames(counts)

sample_type <- substr(barcodes, 14, 15)
keep_col    <- sample_type == "01"
counts      <- counts[, keep_col, drop = FALSE]
barcodes    <- barcodes[keep_col]

dge  <- DGEList(counts = counts)
cpm0 <- edgeR::cpm(dge)
keep_gene <- rowMeans(cpm0 > 1) >= 0.20
dge  <- dge[keep_gene, , keep.lib.sizes = FALSE]
dge  <- calcNormFactors(dge, method = "TMM")
logCPM   <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
expr_mat <- t(logCPM)

meta <- as.data.frame(SummarizedExperiment::colData(se))
meta <- meta[keep_col, , drop = FALSE]

d2d <- suppressWarnings(as.numeric(pick(meta, c("days_to_death","days_to_death.x","paper_days_to_death"))))
d2f <- suppressWarnings(as.numeric(pick(meta, c("days_to_last_follow_up","days_to_lastfollowup",
                                                "days_to_last_follow_up.x","paper_days_to_last_followup"))))
vs  <- tolower(as.character(pick(meta, c("vital_status","vital_status.x","paper_vital_status"))))

BRCA_time <- ifelse(!is.na(d2d), d2d, d2f)
BRCA_status <- ifelse(!is.na(d2d) & d2d > 0, 1,
                      ifelse(vs %in% c("dead","deceased"), 1, 0))

madv <- matrixStats::colMads(expr_mat, na.rm = TRUE)
sel <- which(is.finite(madv) & madv > 0)
BRCA_X0 <- scale(expr_mat[, sel, drop = FALSE])

BRCA_ok <- is.finite(BRCA_time) &
  is.finite(BRCA_status) &
  BRCA_time > 0 &
  rowSums(is.finite(BRCA_X0)) == ncol(BRCA_X0)
BRCA_X <- BRCA_X0[BRCA_ok, , drop = FALSE]
BRCA_time <- BRCA_time[BRCA_ok]
BRCA_status <- BRCA_status[BRCA_ok]

stopifnot(nrow(BRCA_X) == length(BRCA_time), length(BRCA_time) == length(BRCA_status))

cat("BRCA number of patients:", nrow(BRCA_X), "\n")
cat("BRCA number of events:", sum(BRCA_status == 1), "\n")
cat("BRCA censoring rate:", mean(BRCA_status == 0), "\n")
cat("BRCA final number of genes:", ncol(BRCA_X), "\n")
cat("BRCA median follow-up:", median(BRCA_time, na.rm = TRUE), "\n")
cat("BRCA follow-up range:", range(BRCA_time, na.rm = TRUE), "\n")


# KIRC
# ==== TCGA-KIRC ====
qry_expr <- GDCquery(
  project       = "TCGA-KIRC",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(qry_expr, files.per.chunk = 40)  
se <- GDCprepare(qry_expr, summarizedExperiment = TRUE)

counts   <- SummarizedExperiment::assay(se)
barcodes <- colnames(counts)

sample_type <- substr(barcodes, 14, 15)
keep_col    <- sample_type == "01"
counts      <- counts[, keep_col, drop = FALSE]
barcodes    <- barcodes[keep_col]

dge  <- DGEList(counts = counts)
cpm0 <- edgeR::cpm(dge)
keep_gene <- rowMeans(cpm0 > 1) >= 0.20
dge  <- dge[keep_gene, , keep.lib.sizes = FALSE]
dge  <- calcNormFactors(dge, method = "TMM")
logCPM   <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
expr_mat <- t(logCPM)

meta <- as.data.frame(SummarizedExperiment::colData(se))
meta <- meta[keep_col, , drop = FALSE]

d2d <- suppressWarnings(as.numeric(pick(meta, c("days_to_death","days_to_death.x","paper_days_to_death"))))
d2f <- suppressWarnings(as.numeric(pick(meta, c("days_to_last_follow_up","days_to_lastfollowup",
                                                "days_to_last_follow_up.x","paper_days_to_last_followup"))))
vs  <- tolower(as.character(pick(meta, c("vital_status","vital_status.x","paper_vital_status"))))

KIRC_time <- ifelse(!is.na(d2d), d2d, d2f)
KIRC_status <- ifelse(!is.na(d2d) & d2d > 0, 1,
                      ifelse(vs %in% c("dead","deceased"), 1, 0))

madv <- matrixStats::colMads(expr_mat, na.rm = TRUE)
sel <- which(is.finite(madv) & madv > 0)
KIRC_X0 <- scale(expr_mat[, sel, drop = FALSE])

KIRC_ok <- is.finite(KIRC_time) &
  is.finite(KIRC_status) &
  KIRC_time > 0 &
  rowSums(is.finite(KIRC_X0)) == ncol(KIRC_X0)
KIRC_X <- KIRC_X0[KIRC_ok, , drop = FALSE]
KIRC_time <- KIRC_time[KIRC_ok]
KIRC_status <- KIRC_status[KIRC_ok]

stopifnot(nrow(KIRC_X) == length(KIRC_time), length(KIRC_time) == length(KIRC_status))

cat("KIRC number of patients:", nrow(KIRC_X), "\n")
cat("KIRC number of events:", sum(KIRC_status == 1), "\n")
cat("KIRC censoring rate:", mean(KIRC_status == 0), "\n")
cat("KIRC final number of genes:", ncol(KIRC_X), "\n")
cat("KIRC median follow-up:", median(KIRC_time, na.rm = TRUE), "\n")
cat("KIRC follow-up range:", range(KIRC_time, na.rm = TRUE), "\n")




project_dir <- "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls"
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

X <- HNSC_X
time <- HNSC_time
status <- HNSC_status

X <- GBM_X
time <- GBM_time
status <- GBM_status

X <- LUAD_X
time <- LUAD_time
status <- LUAD_status

X <- BRCA_X
time <- BRCA_time
status <- BRCA_status

X <- KIRC_X
time <- KIRC_time
status <- KIRC_status

set.seed(2026)
tcga_surv_ok <- is.finite(time) & is.finite(status) & time > 0 &
  rowSums(is.finite(as.matrix(X))) == ncol(X)
if (any(!tcga_surv_ok)) {
  message(sum(!tcga_surv_ok), " samples removed before CV because of invalid or non-positive survival time.")
  X <- X[tcga_surv_ok, , drop = FALSE]
  time <- time[tcga_surv_ok]
  status <- status[tcga_surv_ok]
}
n <- nrow(X)
test_prop_tcga <- 0.3
split_tcga <- .cvtest_make_split(status, test_prop = test_prop_tcga, seed = 2017)
train_idx_tcga <- split_tcga$train_idx
test_idx_tcga <- split_tcga$test_idx
cat("Training samples:", length(train_idx_tcga), "\n")
cat("Test samples:", length(test_idx_tcga), "\n")
cat("Test events:", sum(status[test_idx_tcga] == 1), "\n")


event_time <- time[train_idx_tcga][status[train_idx_tcga] == 1]

auc_time_grid_tcga <- as.numeric(quantile(
  event_time,
  probs = seq(0.2, 0.8, length.out = 30),
  na.rm = TRUE
))

auc_time_grid_tcga <- sort(unique(auc_time_grid_tcga))
auc_time_grid_tcga

ncomp_candidates_DRPLS <- 1:10
ncomp_candidates_pcox <- 2:10
alpha_candidates <- seq(0, 1, by = 0.1)
lambda_candidates_ridgecox <- 10 ^ seq(-4, 1, by = 1)
alpha_candidates_enetcox <- c(0.25, 0.5, 0.75)
lambda_candidates_enetcox <- 10 ^ seq(-5, 1, by = 1)
use_preselection <- TRUE

auc_curves_list_drpls <- val_drpls_cv_test(
  X, time, status,
  test_prop = test_prop_tcga,
  train_idx = train_idx_tcga,
  test_idx = test_idx_tcga,
  k = 5,
  ncomp_candidates = ncomp_candidates_DRPLS,
  auc_time_grid = auc_time_grid_tcga,
  use_preselection = use_preselection,
  seed = 1
)
print(auc_curves_list_drpls$best_ncomp)
print(auc_curves_list_drpls$best_cv_iAUC)
print(auc_curves_list_drpls$test_iAUC)

auc_curves_list_drpcapls <- val_conticox_cv_test(
  X, time, status,
  test_prop = test_prop_tcga,
  train_idx = train_idx_tcga,
  test_idx = test_idx_tcga,
  k = 5,
  ncomp_candidates = ncomp_candidates_DRPLS,
  alpha_candidates = alpha_candidates,
  auc_time_grid = auc_time_grid_tcga,
  use_preselection = use_preselection,
  normalize_alpha_terms = TRUE,
  residual_type = "martingale",
  seed = 1
)

print(auc_curves_list_drpcapls$best_ncomp)
print(auc_curves_list_drpcapls$best_alpha)
print(auc_curves_list_drpcapls$best_cv_iAUC)
print(auc_curves_list_drpcapls$test_iAUC)

auc_curves_list_pcox <- val_partial_cox_cv_test(
  X, time, status,
  test_prop = test_prop_tcga,
  train_idx = train_idx_tcga,
  test_idx = test_idx_tcga,
  k = 5,
  ncomp_candidates = ncomp_candidates_pcox,
  auc_time_grid = auc_time_grid_tcga,
  use_preselection = use_preselection,
  seed = 1
)
print(auc_curves_list_pcox$best_ncomp)
print(auc_curves_list_pcox$best_cv_iAUC)
print(auc_curves_list_pcox$test_iAUC)

auc_curves_list_ridgecox <- val_ridge_cox_cv_test(
  X, time, status,
  test_prop = test_prop_tcga,
  train_idx = train_idx_tcga,
  test_idx = test_idx_tcga,
  k = 5,
  lambda_candidates = lambda_candidates_ridgecox,
  auc_time_grid = auc_time_grid_tcga,
  use_preselection = use_preselection,
  seed = 2026
)
print(auc_curves_list_ridgecox$best_lambda)
print(auc_curves_list_ridgecox$best_cv_iAUC)
print(auc_curves_list_ridgecox$test_iAUC)

auc_curves_list_ridgecox$test_auc_curve
summary(auc_curves_list_ridgecox$test_marker)
sd(auc_curves_list_ridgecox$test_marker)
length(unique(round(auc_curves_list_ridgecox$test_marker, 8)))
auc_curves_list_ridgecox$best_lambda


auc_curves_list_enetcox <- val_elastic_net_cox_cv_test(
  X, time, status,
  test_prop = test_prop_tcga,
  train_idx = train_idx_tcga,
  test_idx = test_idx_tcga,
  k = 5,
  alpha_candidates = alpha_candidates_enetcox,
  lambda_candidates = lambda_candidates_enetcox,
  auc_time_grid = auc_time_grid_tcga,
  use_preselection = use_preselection,
  seed = 2026
)
print(auc_curves_list_enetcox$best_alpha)
print(auc_curves_list_enetcox$best_lambda)
print(auc_curves_list_enetcox$best_cv_iAUC)
print(auc_curves_list_enetcox$test_iAUC)

time_points <- auc_time_grid_tcga
auc_curve_results <- list(
  PLSDR = auc_curves_list_drpls,
  ContiCox = auc_curves_list_drpcapls,
  pCox = auc_curves_list_pcox,
  `Ridge-Cox` = auc_curves_list_ridgecox,
  `EN-Cox` = auc_curves_list_enetcox
)
auc_curve_values <- lapply(auc_curve_results, `[[`, "best_auc_curve")

plot_auc_curves <- function(
    main_title = "Predictive accuracy (AUC)",
    ylim_use = c(0.6, 0.9)
) {
  cols <- c(
    "PLSDR" = "#4C78A8",
    "ContiCox" = "#E45756",
    "pCox" = "#72B7B2",
    "Ridge-Cox" = "#B279A2",
    "EN-Cox" = "#79706E"
  )
  ltys <- c(
    "PLSDR" = 1,
    "ContiCox" = 1,
    "pCox" = 1,
    "Ridge-Cox" = 2,
    "EN-Cox" = 3
  )
  lwds <- c(
    "PLSDR" = 2.2,
    "ContiCox" = 2.8,
    "pCox" = 2.2,
    "Ridge-Cox" = 1.8,
    "EN-Cox" = 1.8
  )
  method_names <- names(auc_curve_values)
  
  plot(
    time_points,
    auc_curve_values[[1]],
    type = "n",
    ylim = ylim_use,
    xlab = "Time (days)",
    ylab = "Predictive accuracy (AUC)",
    main = main_title
  )
  abline(h = 0.5, lty = 3, col = "grey80")
  
  for (method in method_names) {
    lines(
      time_points,
      auc_curve_values[[method]],
      col = cols[method],
      lty = ltys[method],
      lwd = lwds[method]
    )
  }
  
  legend_obj <- legend(
    "topright",
    legend = method_names,
    col = cols[method_names],
    lty = ltys[method_names],
    lwd = lwds[method_names],
    bty = "n",
    bg = "white",
    cex = 0.95,
    seg.len = 2.6,
    x.intersp = 0.95,
    y.intersp = 1.05,
    inset = 0
  )
  usr <- par("usr")
  segments(legend_obj$rect$left, legend_obj$rect$top,
           legend_obj$rect$left, legend_obj$rect$top - legend_obj$rect$h,
           lwd = 1.15, xpd = FALSE)
  segments(legend_obj$rect$left, legend_obj$rect$top - legend_obj$rect$h,
           usr[2], legend_obj$rect$top - legend_obj$rect$h,
           lwd = 1.15, xpd = FALSE)
}

plot_auc_curves(main_title = "Predictive accuracy (AUC)")

tcga_dataset_name <- "TCGA_KIRC"
tcga_auc_results <- list(
  dataset = tcga_dataset_name,
  auc_time_grid = time_points,
  train_idx = train_idx_tcga,
  test_idx = test_idx_tcga,
  test_prop = test_prop_tcga,
  PLSDR = auc_curves_list_drpls,
  ContiCox = auc_curves_list_drpcapls,
  pCox = auc_curves_list_pcox,
  RidgeCox = auc_curves_list_ridgecox,
  ElasticNetCox = auc_curves_list_enetcox
)

save(
  auc_curves_list_drpls,
  auc_curves_list_drpcapls,
  auc_curves_list_pcox,
  auc_curves_list_ridgecox,
  auc_curves_list_enetcox,
  tcga_auc_results,
  file = file.path(project_dir, paste0(tcga_dataset_name, "_validation_results_BRCA.RData"))
)

load(file.path(project_dir, paste0("TCGA_GBM", "_validation_results_HNSC.RData")))

save_tcga_auc_curve <- function(
    file_base,
    width = 6.5,
    height = 4.2,
    res = 600,
    main_title = "Predictive accuracy (AUC)"
) {
  dir.create(dirname(file_base), recursive = TRUE, showWarnings = FALSE)
  
  pdf_file <- paste0(file_base, ".pdf")
  png_file <- paste0(file_base, ".png")
  
  grDevices::pdf(pdf_file, width = width, height = height,
                 pointsize = 10, useDingbats = FALSE)
  plot_auc_curves(main_title = main_title)
  grDevices::dev.off()
  
  grDevices::png(png_file, width = width, height = height,
                 units = "in", res = res, pointsize = 10)
  plot_auc_curves(main_title = main_title)
  grDevices::dev.off()
  
  invisible(list(pdf = pdf_file, png = png_file))
}

tcga_auc_curve_files <- save_tcga_auc_curve(
  file_base = file.path(project_dir, "results", paste0(tcga_dataset_name, "_auc_curve_0519"))
)
print(tcga_auc_curve_files)



# save data
combine_tcga_validation_results <- function(
    project_dir,
    pattern = "_validation_results\\.RData$",
    output_dir = file.path(project_dir, "results")
) {
  method_sources <- list(
    PLSDR = c("PLSDR", "DRPLS", "auc_curves_list_drpls"),
    ContiCox = c("ContiCox", "DRPCA_PLS", "auc_curves_list_drpcapls"),
    pCox = c("pCox", "partialCox", "auc_curves_list_pcox"),
    `Ridge-Cox` = c("RidgeCox", "auc_curves_list_ridgecox"),
    `EN-Cox` = c("ElasticNetCox", "auc_curves_list_enetcox"),
    PCGCox = c("PCGCox", "PCARegCox", "auc_curves_list_pca_pcox")
  )
  
  extract_iAUC <- function(x) {
    if (!is.null(x$best_iAUC)) return(as.numeric(x$best_iAUC))
    if (!is.null(x$metrics$iAUC)) return(as.numeric(x$metrics$iAUC))
    NA_real_
  }
  
  extract_method_summary <- function(dataset, method, x) {
    data.frame(
      dataset = dataset,
      method = method,
      best_iAUC = extract_iAUC(x),
      best_ncomp = if (!is.null(x$best_ncomp)) as.numeric(x$best_ncomp) else NA_real_,
      best_alpha = if (!is.null(x$best_alpha)) as.numeric(x$best_alpha) else NA_real_,
      best_lambda = if (!is.null(x$best_lambda)) as.numeric(x$best_lambda) else NA_real_,
      residual_type = if (!is.null(x$residual_type)) as.character(x$residual_type) else NA_character_,
      normalize_alpha_terms = if (!is.null(x$normalize_alpha_terms)) {
        as.character(x$normalize_alpha_terms)
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )
  }
  
  extract_auc_curve <- function(dataset, method, x, fallback_grid = NULL) {
    auc <- x$best_auc_curve
    if (is.null(auc)) return(NULL)
    
    time_grid <- x$auc_time_grid
    if (is.null(time_grid)) time_grid <- fallback_grid
    if (is.null(time_grid)) return(NULL)
    
    n_use <- min(length(time_grid), length(auc))
    data.frame(
      dataset = dataset,
      method = method,
      time = as.numeric(time_grid[seq_len(n_use)]),
      AUC = as.numeric(auc[seq_len(n_use)]),
      stringsAsFactors = FALSE
    )
  }
  
  validation_files <- list.files(project_dir, pattern = pattern, full.names = TRUE)
  if (length(validation_files) == 0L) {
    stop("No validation result files found in: ", project_dir)
  }
  
  summary_list <- list()
  curve_list <- list()
  
  for (file_i in validation_files) {
    env_i <- new.env(parent = emptyenv())
    load(file_i, envir = env_i)
    
    if (exists("tcga_auc_results", envir = env_i, inherits = FALSE)) {
      res_i <- get("tcga_auc_results", envir = env_i)
    } else {
      res_i <- list(
        dataset = sub(pattern, "", basename(file_i)),
        auc_time_grid = NULL
      )
    }
    
    dataset_i <- if (!is.null(res_i$dataset)) {
      as.character(res_i$dataset)
    } else {
      sub(pattern, "", basename(file_i))
    }
    
    for (method_i in names(method_sources)) {
      method_obj <- NULL
      for (source_name in method_sources[[method_i]]) {
        if (!is.null(res_i[[source_name]])) {
          method_obj <- res_i[[source_name]]
          break
        }
        if (exists(source_name, envir = env_i, inherits = FALSE)) {
          method_obj <- get(source_name, envir = env_i)
          break
        }
      }
      
      if (!is.null(method_obj)) {
        summary_list[[length(summary_list) + 1L]] <-
          extract_method_summary(dataset_i, method_i, method_obj)
        curve_i <- extract_auc_curve(
          dataset = dataset_i,
          method = method_i,
          x = method_obj,
          fallback_grid = res_i$auc_time_grid
        )
        if (!is.null(curve_i)) {
          curve_list[[length(curve_list) + 1L]] <- curve_i
        }
      }
    }
  }
  
  tcga_validation_summary <- do.call(rbind, summary_list)
  tcga_auc_curve_long <- if (length(curve_list) > 0L) {
    do.call(rbind, curve_list)
  } else {
    data.frame()
  }
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  summary_csv <- file.path(output_dir, "TCGA_validation_summary.csv")
  curve_csv <- file.path(output_dir, "TCGA_auc_curve_long.csv")
  rdata_file <- file.path(output_dir, "TCGA_validation_combined_results.RData")
  
  utils::write.csv(tcga_validation_summary, summary_csv, row.names = FALSE)
  utils::write.csv(tcga_auc_curve_long, curve_csv, row.names = FALSE)
  save(
    tcga_validation_summary,
    tcga_auc_curve_long,
    validation_files,
    file = rdata_file
  )
  
  list(
    summary = tcga_validation_summary,
    auc_curve_long = tcga_auc_curve_long,
    files = list(
      summary_csv = summary_csv,
      curve_csv = curve_csv,
      rdata = rdata_file
    )
  )
}

if (!exists("project_dir", inherits = FALSE)) {
  project_dir <- "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls"
}
tcga_combined_results <- combine_tcga_validation_results(project_dir)
print(tcga_combined_results$files)
print(tcga_combined_results$summary)
