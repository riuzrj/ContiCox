
# install.packages("/Users/ruijuanzhong/Downloads/IMvigor210CoreBiologies_1.0.1.tar.gz", repos = NULL,
#                  type  = "source")
library(IMvigor210CoreBiologies)
library(DESeq)
library(dplyr)
library(survival)
library(survminer)
library(pls) 

data(cds) 

# --- 1. 数据提取 ---
expr_data <- counts(cds)     # 基因表达矩阵 (行=基因, 列=样本)
pheno_data <- pData(cds)     # 临床表型数据

# 构建生存数据
surv_data <- pheno_data %>%
  transmute(
    sample_id = rownames(pheno_data),
    time      = as.numeric(os),      # 确保转为数值型
    status    = as.numeric(censOS)   # 确保转为数值型
  )

# 检查样本对齐 (这是一个好习惯)
if (!all(colnames(expr_data) == surv_data$sample_id)) {
  stop("Error: 表达矩阵列名与生存数据行名不匹配，请先对齐数据！")
}

# --- 2. 构建初始矩阵 (转置：行=样本，列=基因) ---
X_raw <- t(expr_data) 
time_raw <- surv_data$time
status_raw <- surv_data$status

# ⚠️ 严重错误修正提醒：
# 你原来的代码：X <- as.numeric(X) 
# 这会将整个矩阵拉直成一条极长的向量，导致丢失行/列结构！
# 正确做法是确保矩阵模式为 numeric，但保持矩阵结构：
if (!is.numeric(X_raw)) {
  mode(X_raw) <- "numeric"
}

# --- 3. 核心步骤：处理缺失值 (NA) ---
# 找出在 X, time, status 中完全没有缺失值的样本索引
# 大多数情况下 expr_data 是完整的，但 time/status 容易缺
keep_idx <- complete.cases(X_raw, time_raw, status_raw)

# 打印日志看看删除了多少数据
n_total <- length(time_raw)
n_keep  <- sum(keep_idx)
cat("原始样本量:", n_total, "\n")
cat("清洗后样本量:", n_keep, "\n")
cat("剔除含 NA 样本数:", n_total - n_keep, "\n")

# 应用筛选
X <- X_raw[keep_idx, , drop = FALSE]
time <- time_raw[keep_idx]
status <- status_raw[keep_idx]

# --- 4. 数据标准化 (可选，但推荐) ---
# 你的代码里提到了 X_final_noNA，这里我们用清洗后的 X 代替
# 对 X 进行中心化 (Centering)
X_center <- colMeans(X)
X <- sweep(X, 2, X_center, FUN = "-")


# --- 5. 生成 Folds ---
# 必须在剔除 NA 之后，重新计算 n，再生成 folds
n <- nrow(X) 

# 设置种子确保可复现
set.seed(123) 
common_folds <- sample(rep(1:5, length.out = n))


IMvigor_drpls <- val_DRPLS_cv(X, time, status, k = 5,
                           ncomp_candidates = 1:8, auc_time_grid = NULL,
                           use_preselection = FALSE, 
                           auc_method = "NNE")



IMvigor_drpcapls <- val_DRPCAPLS_cv(X, time, status, k = 5,
                                 ncomp_candidates = 1:10, auc_time_grid = NULL,
                                 alpha_candidates = seq(0, 1, by = 0.1),
                                 use_preselection = FALSE, p_thresh = 0.05,
                                 auc_method = "NNE")

IMvigor_pcox <- val_partial_cox_cv(X, time, status, k =5,
                                ncomp_candidates = 2:10, auc_time_grid = NULL,
                                use_preselection = FALSE, p_thresh = 0.05,
                                auc_method = "NNE")

IMvigor_pcapcox <- val_partial_cox_pca_cv(X, time, status, k = 5,
                                       ncomp_candidates = 2:5, auc_time_grid = NULL,
                                       alpha_candidates = seq(0, 1, by = 0.1),
                                       use_preselection = FALSE, p_thresh = 0.05,
                                       auc_method = "NNE")

IMvigor_wpca_pcox <- val_partial_cox_pca_weight_cv(X, time, status, k = 5,
                                                ncomp_candidates = 2:9, auc_time_grid = NULL,
                                                use_preselection = FALSE, p_thresh = 0.05,
                                                auc_method = "NNE")

# plot
par(mar = c(5, 4, 4, 6), xpd = TRUE)
time_points = IMvigor_drpls$auc_time_grid
plot(time_points, IMvigor_drpls$best_auc_curve, type = "l", col = "blue",
     xlab = "Time (days)", ylab = "AUC", ylim = c(0.3, 1), lwd = 2,
     main = "AUC curves for different methods")
lines(time_points, IMvigor_drpcapls$best_auc_curve, col = "red", lwd = 2)
lines(time_points, IMvigor_pcox$best_auc_curve, col = "green", lwd = 2)
lines(time_points, IMvigor_pcapcox$best_auc_curve, col = "purple", lwd = 2)
lines(time_points, IMvigor_wpca_pcox$best_auc_curve, col = "orange", lwd = 2)
legend("topright", inset = c(-0.4, 0), legend = c("DRPLS","DRPCAPLS",
                                                  "Partial Cox", 
                                                  "Partial PCA-Cox",
                                                  "Weighted PCA-Cox"),
       col = c("red", "blue", "green", "purple"), lwd = 1, cex = 0.5, bty = "n")



# ======= high-low risk ==========
# —— 1. 随机划分训练/测试集 ——  
set.seed(42)
n          <- nrow(X)
train_frac <- 0.7
train_idx  <- sample(seq_len(n), size = floor(train_frac * n))
test_idx   <- setdiff(seq_len(n), train_idx)

X_train    <- X[train_idx, , drop=FALSE]
time_train <- time[train_idx]
status_train <- status[train_idx]

X_test     <- X[test_idx, , drop=FALSE]
time_test  <- time[test_idx]
status_test  <- status[test_idx]


# ========== DRPLS =============
res <- cox_pls_dr(
  X            = X_train,
  time         = time_train,
  status       = status_train,
  residual_type= "deviance",
  scale_x      = TRUE,
  scale_time   = FALSE,
  n_components = 6,        
  return_all   = TRUE
)

# 提取训练时用到的缩放参数和投影矩阵 A
X_center        <- res$X_center   # 长度 = p
X_scale         <- res$X_scale    # 长度 = p
A               <- res$A          # p × ncomp
final_cox_model <- res$final_cox_model
pls_scores_train <- res$pls_scores  # 训练集上各成分得分

# —— 3. 在测试集上计算 PLS 得分 ——  
# 先中心化 & 缩放
X_test_scaled <- sweep(X_test, 2, X_center, "-")
X_test_scaled <- sweep(X_test_scaled,  2, X_scale, "/")
# 然后投影
pls_scores_test <- as.data.frame(as.matrix(X_test_scaled) %*% A)
colnames(pls_scores_test) <- colnames(pls_scores_train)

# —— 4. 计算线性预测值（风险得分） ——  
risk_train <- predict(final_cox_model,
                      newdata = pls_scores_train,
                      type    = "lp")
risk_test  <- predict(final_cox_model,
                      newdata = pls_scores_test,
                      type    = "lp")

# —— 5. 以训练集平均风险得分为阈值分组 ——  
cutoff      <- mean(risk_train)
group_test  <- ifelse(risk_test >= cutoff, "high-risk", "low-risk")

# —— 6. 构建测试集生存数据框 ——  
test_df <- data.frame(
  time   = time_test,
  status = status_test,
  group  = factor(group_test, levels = c("low-risk","high-risk"))
)

# —— 7. 拟合 KM 曲线 & 计算对数秩检验 ——  
fit_test <- survfit(Surv(time, status) ~ group, data = test_df)
logr     <- survdiff(Surv(time, status) ~ group, data = test_df)
p_val    <- 1 - pchisq(logr$chisq, df = length(logr$n) - 1)

# —— 8. 绘制生存曲线 ——  
ggsurvplot(
  fit_test,
  data        = test_df,
  pval        = sprintf("p = %.4f", p_val),
  legend.labs = c("low-risk patients","high-risk patients"),
  legend.title= NULL,
  palette     = c("black","darkgrey"),
  linetype    = c("solid","dashed"),
  censor.shape= 3,
  xlab        = "Survival in Years",
  ylab        = "Survival Probability",
  ggtheme     = theme_bw(base_size = 14)
)


# =========== DRPCAPLS =============
res <- cox_pca_pls_dr(
  X            = X_train,
  time         = time_train,
  status       = status_train,
  scale_time   = FALSE,   # 保留原始时间
  n_components = 1,       # 你想要的成分数
  alpha        = 0.1,     # 超参数 alpha
  return_all   = TRUE
)

pls_scores_train <- res$pls_scores
A                <- res$A
final_cox_model  <- res$final_cox_model

# —— 4. 测试集上计算 PLS 得分 ——  
pls_scores_test <- as.data.frame(as.matrix(X_test) %*% A)
colnames(pls_scores_test) <- colnames(pls_scores_train)

# —— 5. 计算风险得分 ——  
risk_train <- predict(final_cox_model,
                      newdata = pls_scores_train,
                      type    = "lp")
risk_test  <- predict(final_cox_model,
                      newdata = pls_scores_test,
                      type    = "lp")

# —— 6. 以训练集平均风险得分为阈值分组 ——  
cutoff     <- mean(risk_train)
group_test <- ifelse(risk_test >= cutoff, "high-risk", "low-risk")

# —— 7. 构建测试集生存数据框 ——  
test_df <- data.frame(
  time   = time_test,
  status = status_test,
  group  = factor(group_test, levels = c("low-risk","high-risk"))
)

# —— 8. 拟合 KM 曲线 & 对数秩检验 ——  
fit_km <- survfit(Surv(time, status) ~ group, data = test_df)
logr   <- survdiff(Surv(time, status) ~ group, data = test_df)
p_val  <- 1 - pchisq(logr$chisq, df = length(logr$n) - 1)

# —— 9. 绘制 Kaplan–Meier 曲线 ——  
ggsurvplot(
  fit_km,
  data        = test_df,
  pval        = sprintf("p = %.4f", p_val),
  legend.labs = c("low-risk patients","high-risk patients"),
  legend.title= NULL,
  palette     = c("black","darkgrey"),
  linetype    = c("solid","dashed"),
  censor.shape= 3,
  xlab        = "Survival in Years",
  ylab        = "Survival Probability",
  ggtheme     = theme_bw(base_size = 14)
)



# =========== Partial cox =============
res <- partial_cox(
  X            = X_train,
  time         = time_train,
  status       = status_train,
  n_components = 10      # 你选择的组件数
)

T_train         <- res$components        # 已是中心化后的分量得分
final_cox_model <- res$final_model

# —— 测试集直接投影（无需再中心化） ——  
T_test <- transform_partial_cox(
  X_new      = as.matrix(X_test),
  betas_list = res$betas_list,
  w_list     = res$w_list
)

# —— 计算风险得分 ——  
risk_train <- predict(final_cox_model,
                      newdata = as.data.frame(T_train),
                      type    = "lp")
risk_test  <- predict(final_cox_model,
                      newdata = as.data.frame(T_test),
                      type    = "lp")

# —— 用训练集平均 risk 做阈值分组 ——  
cutoff     <- mean(risk_train)
group_test <- ifelse(risk_test >= cutoff, "high-risk", "low-risk")

# —— 构造测试集生存表 ——  
test_df <- data.frame(
  time   = time_test,
  status = status_test,
  group  = factor(group_test, c("low-risk","high-risk"))
)

# —— 拟合 KM 曲线 & 计算 p 值 ——  
fit_km <- survfit(Surv(time, status) ~ group, data = test_df)
logr   <- survdiff(Surv(time, status) ~ group, data = test_df)
p_val  <- 1 - pchisq(logr$chisq, df = length(logr$n) - 1)

# —— 绘图 ——  
ggsurvplot(
  fit_km,
  data        = test_df,
  pval        = sprintf("p = %.4f", p_val),
  legend.labs = c("low-risk patients", "high-risk patients"),
  legend.title= NULL,
  palette     = c("black", "darkgrey"),
  linetype    = c("solid", "dashed"),
  censor.shape= 3,
  xlab        = "Survival Time",
  ylab        = "Survival Probability",
  ggtheme     = theme_bw(base_size = 14)
)


# ============ PCA Partial Cox =============
res <- partial_cox_pca(
  X            = X_train,
  time         = time_train,
  status       = status_train,
  n_components = 3,    # 选 5 个成分，按需修改
  alpha        = 0.3   # P‐PCA 中的 alpha 超参
)

# 提取训练时的 T_matrix 和最终 Cox 模型
T_train        <- res$components       # n_train × n_components
final_cox_mod  <- res$final_model

# —— 4. 在测试集上直接投影计算 T_test ——  
# （因为 X 已经在外面中心化过，无需再减均值）
T_test <- transform_partial_cox_pca(
  X_new        = X_test,
  betas_list   = res$betas_list,
  weights_list = res$weights_list
)

final_cox_mod <- if (is.null(res$final_model)) {
  message("Refitting final Cox model with coxph()")
  coxph(
    Surv(time_train, status_train) ~ .,
    data    = as.data.frame(T_train),
    control = coxph.control(iter.max = 100)
  )
} else {
  res$final_model
}


# —— 5. 计算线性预测值（risk score） ——  
risk_train <- predict(final_cox_mod,
                      newdata = as.data.frame(T_train),
                      type    = "lp")
risk_test  <- predict(final_cox_mod,
                      newdata = as.data.frame(T_test),
                      type    = "lp")

# —— 6. 以训练集平均 risk 为阈值对测试集分组 ——  
cutoff     <- mean(risk_train)
group_test <- ifelse(risk_test >= cutoff, "high-risk", "low-risk")
group_test <- factor(group_test, levels = c("low-risk","high-risk"))

# —— 7. 构建测试集生存数据框 ——  
test_df <- data.frame(
  time   = time_test,
  status = status_test,
  group  = group_test
)

# —— 8. 拟合 Kaplan–Meier 曲线 & 对数秩检验 ——  
fit_km <- survfit(Surv(time, status) ~ group, data = test_df)
logres <- survdiff(Surv(time, status) ~ group, data = test_df)
p_val  <- 1 - pchisq(logres$chisq, df = length(logres$n) - 1)

# —— 9. 绘制生存曲线 ——  
ggsurvplot(
  fit_km,
  data        = test_df,
  pval        = sprintf("p = %.4f", p_val),
  legend.labs = c("low-risk patients", "high-risk patients"),
  legend.title= NULL,
  palette     = c("black","darkgrey"),
  linetype    = c("solid","dashed"),
  censor.shape= 3,
  xlab        = "Survival Time",
  ylab        = "Survival Probability",
  ggtheme     = theme_bw(base_size = 14)
)



# ============ Pre PCA Partial Cox =============
use_pcr    <- TRUE 

if (use_pcr) {
  pca_obj   <- prcomp(X_train, center = TRUE, scale. = FALSE)
  scores_tr <- pca_obj$x
  # 在测试集上做同样的中心化 + 投影
  X_test_ctr <- sweep(X_test, 2, pca_obj$center, "-")
  scores_te  <- as.matrix(X_test_ctr) %*% pca_obj$rotation
} else {
  scores_tr <- as.matrix(X_train)
  scores_te <- as.matrix(X_test)
}

# —— 4. 在训练集上拟合 Partial Cox ——  
res <- partial_cox(
  X            = scores_tr,
  time         = time_train,
  status       = status_train,
  n_components = 6
)

# 提取训练时的组件矩阵和最终模型
T_train       <- res$components       # n_train × ncomp
final_cox_mod <- res$final_model      # coxph 对象

# —— 5. 在测试集上计算组件得分 ——  
T_test <- transform_partial_cox(
  X_new      = scores_te,
  betas_list = res$betas_list,
  w_list     = res$w_list
)

# —— 6. 计算线性预测值（risk score） ——  
risk_train <- predict(final_cox_mod,
                      newdata = as.data.frame(T_train),
                      type    = "lp")
risk_test  <- predict(final_cox_mod,
                      newdata = as.data.frame(T_test),
                      type    = "lp")
#risk_test

# —— 7. 以训练集平均风险为阈值对测试集分组 ——  
cutoff     <- mean(risk_train)
group_test <- factor(
  ifelse(risk_test >= cutoff, "high-risk", "low-risk"),
  levels = c("low-risk","high-risk")
)

# —— 8. 构建测试集生存数据框 ——  
# 构造原始的 test_df
test_df <- data.frame(
  time   = time_test,
  status = status_test,
  group  = group_test
)

# 过滤掉任何含 NA 的行
test_df_clean <- subset(test_df,
                        !is.na(time)   &
                          !is.na(status) &
                          !is.na(group))

# 检查下到底剩多少行
if (nrow(test_df_clean) == 0) {
  stop("After filtering, no non-missing observations remain in test_df.")
}

# 再来拟合生存曲线
fit_km <- survfit(Surv(time, status) ~ group, data = test_df_clean)

# 对数秩检验
logr  <- survdiff(Surv(time, status) ~ group, data = test_df_clean)
p_val <- 1 - pchisq(logr$chisq, df = length(logr$n) - 1)

# 绘图
ggsurvplot(
  fit_km,
  data        = test_df_clean,
  pval        = sprintf("p = %.4f", p_val),
  legend.labs = c("low-risk patients","high-risk patients"),
  legend.title= NULL,
  palette     = c("black","darkgrey"),
  linetype    = c("solid","dashed"),
  censor.shape= 3,
  xlab        = "Survival Time",
  ylab        = "Survival Probability",
  ggtheme     = theme_bw(base_size = 14)
)


# ============ PCAPLS Partial Cox =============
null_mod <- coxph(Surv(time_train, status_train) ~ 1)
Y_tr     <- residuals(null_mod, type = "deviance")

# 3.2 中心化 X_train
Xc_tr <- scale(X_train, center = TRUE, scale = FALSE)
cen   <- attr(Xc_tr, "scaled:center")

# 3.3 用 PCA-PLS 提取分量
alpha <- 0                  # 按需调整
ncomp <- 6                  # 你想要的组件数
pls_res <- acc_pca_pls_with_alpha(
  X     = Xc_tr,
  Y     = Y_tr,
  alpha = alpha,
  ncomp = ncomp
)

# pls_res <- pca_pls_with_alpha(
#   X     = Xc_tr,
#   Y     = Y_tr,
#   alpha = alpha,
#   ncomp = ncomp
# )

# 3.4 构造训练集的低维矩阵
X_tr_mod <- as.matrix(Xc_tr) %*% pls_res$Projection

# 3.5 拟合 Partial Cox
model <- partial_cox(
  X            = X_tr_mod,
  time         = time_train,
  status       = status_train,
  n_components = ncomp
)

# —— 4. 在测试集上投影 ——  
# 4.1 中心化测试集（跟训练集同样的中心值）
Xc_te    <- sweep(X_test, 2, cen, "-")
# 4.2 转成低维
X_te_mod <- as.matrix(Xc_te) %*% pls_res$Projection
# 4.3 得到测试集的组件得分
T_test   <- transform_partial_cox(
  X_new      = X_te_mod,
  betas_list = model$betas_list,
  w_list     = model$w_list
)

# —— 5. 计算风险得分（linear predictor） ——  
# 训练集得分（可选，用来确定 cutoff）
T_train    <- model$components
risk_train <- predict(model$final_model,
                      newdata = as.data.frame(T_train),
                      type    = "lp")
# 测试集得分
risk_test  <- predict(model$final_model,
                      newdata = as.data.frame(T_test),
                      type    = "lp")

# —— 6. 分组 ——  
cutoff     <- mean(risk_train, na.rm = TRUE)
group_test <- factor(
  ifelse(risk_test >= cutoff, "high-risk", "low-risk"),
  levels = c("low-risk","high-risk")
)

# —— 7. 构建测试集生存数据框并清洗 ——  
test_df <- data.frame(
  time   = time_test,
  status = status_test,
  group  = group_test
)
test_df <- subset(test_df,
                  !is.na(time)   &
                    !is.na(status) &
                    !is.na(group))

# —— 8. 拟合 KM 曲线 & 计算对数秩检验 ——  
fit_km <- survfit(Surv(time, status) ~ group, data = test_df)
logr   <- survdiff(Surv(time, status) ~ group, data = test_df)
p_val  <- 1 - pchisq(logr$chisq, df = length(logr$n) - 1)

# —— 9. 绘制 Kaplan–Meier 曲线 ——  
ggsurvplot(
  fit_km,
  data        = test_df,
  pval        = sprintf("p = %.4f", p_val),
  legend.labs = c("low-risk patients", "high-risk patients"),
  legend.title= NULL,
  palette     = c("black", "darkgrey"),
  linetype    = c("solid", "dashed"),
  censor.shape= 3,
  xlab        = "Survival Time",
  ylab        = "Survival Probability",
  ggtheme     = theme_bw(base_size = 14)
)







