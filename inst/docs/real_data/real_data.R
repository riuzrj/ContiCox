
library(survivalROC)
library(GEOquery)
# load("/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls/Xmicro.censure_compl_imp.RData")
# load("/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls/micro.censure.RData")

# 数据准备
gse <- getGEO("GSE60", GSEMatrix=TRUE)
exprSet <- exprs(gse[[3]])
pdata <- pData(gse[[3]])

# beer dataset
gse_beer <- getGEO("GSE68571", GSEMatrix=TRUE)

# get features
exprSet_beer <- exprs(gse_beer[[1]])

# get phenotype
pdata_beer <- pData(gse_beer[[1]])

X <- t(exprSet_beer)

# 1. 提取随访时间 (time)
time <- as.numeric(
  gsub("followup time \\(months\\): ", "", pdata_beer$characteristics_ch1.10)
)

# 2. 提取事件状态 (status)
status <- as.numeric(
  gsub("death \\(1=dead,0=alive\\): ", "", pdata_beer$characteristics_ch1.11)
)

# 3. 找到哪些行既有 time 又有 status (非NA)
idx_complete <- !is.na(time) & !is.na(status)

# 4. 删除含缺失值的行
# 如果你想保持 pdata_beer 与 time/status 同步，需要对所有相关对象进行筛选
pdata_beer_clean <- pdata_beer[idx_complete, ]
time_clean <- time[idx_complete]
status_clean <- status[idx_complete]
X_clean <- X[idx_complete, ]

X_final_noNA <- X_clean[, !apply(X_clean, 2, function(x) any(is.na(x)))]

# center X
X_center <- colMeans(X_final_noNA)
X <- sweep(X_final_noNA, 2, X_center)

# Metzeler
gse_Metzeler <- getGEO("GSE12417", GSEMatrix=TRUE)
exprSet_Rosenwald <- exprs(gse_Rosenwald[[1]])
pdata_Rosenwald <- pData(gse_Rosenwald[[1]])
X_Rosenwald <- t(exprSet_Rosenwald)

# ==== compare methods ====
# DRPLS
result_drpls <- cox_pls_dr(X = X, time = time_clean, status = status_clean, 
                           n_components = 5, return_all = TRUE)
summary(result_drpls$final_cox_model)

# DRPCAPLS
result_drpcapls <- cox_pca_pls_dr(X = X_final_noNA, time = time_clean, status = status_clean, 
                                  n_components = 9, return_all = TRUE, alpha = 0.9)
summary(result_drpcapls$final_cox_model)

result_drpcapls_select_alpha <- cox_pca_pls_dr_select_alpha(X = X_final_noNA, time = time_clean, status = status_clean, 
                                                            n_components = 7, return_all = TRUE)
result_drpcapls_select_alpha$results


# partial cox regression
result_pcox <- partial_cox(X = X_final_noNA, time = time_clean, status = status_clean, 
                           n_components = 9)
summary(result_pcox$final_model)

# partial cox regression with pca
result_pcox_pca <- partial_cox_pca(X = X_final_noNA, time = time_clean, status = status_clean, 
                                   n_components = 9, alpha = 0.9)
summary(result_pcox_pca$final_model)

result_pcox_pca_select_alpha <- partial_cox_pca_select_alpha(X = X_final_noNA, time = time_clean, status = status_clean, 
                                                             n_components = 9)
result_pcox_pca_select_alpha$results



# === predict c-indec as criteria ===
# train test split
set.seed(123)  # 为了结果可重复，建议设置随机种子
n <- nrow(X)          # 样本总数
train_size <- round(0.7 * n)     # 训练集数量

# 随机抽取训练集索引
train_index <- sample(seq_len(n), size = train_size)

# 训练集
pdata_train  <- pdata_beer_clean[train_index, ]
time_train   <- time_clean[train_index]
status_train <- status_clean[train_index]
X_train      <- X[train_index, ]

# 测试集
pdata_test  <- pdata_beer_clean[-train_index, ]
time_test   <- time_clean[-train_index]
status_test <- status_clean[-train_index]
X_test      <- X[-train_index, ]

dim(X_train)
dim(X_test)

# DRPLS
result_drpls <- cox_pls_dr(X = X_train, time = time_train, status = status_train, 
                           n_components = 5, return_all = TRUE)
T_test <- X_test %*% result_drpls$A
colnames(T_test) <- paste0("comp_", 1:ncol(T_test))
T_test <- as.data.frame(T_test)
marker_lp <- predict(result_drpls$final_cox_model, newdata = T_test, type = "lp")
marker_risk <- predict(result_drpls$final_cox_model, newdata = T_test, type = "risk")

timepoints <- unique(time_test)
timepoints <- timepoints[order(timepoints)]
AUCs <- numeric(length(timepoints))

for (i in seq_along(timepoints)) {
  timepoint <- timepoints[i]
  AUCs[i] <- survivalROC(Stime = time_test, status = status_test, marker = marker_lp, 
                                       method = "KM", predict.time = timepoint)$AUC
}

plot(
  timepoints, AUCs, 
  type = "l",        # 用线绘制
  lty = 2,           # 虚线
  lwd = 2,           # 线宽
  col = "blue",
  xlab = "Time (months)",
  ylab = "AUC",
  main = "Time-dependent AUC curve (Single Model)"
)
abline(h = 0.5, lty = 2, col = "gray")  # 可选：绘制 AUC=0.5 的参考线

# DRPCAPLS
result_drpcapls <- cox_pca_pls_dr_select_alpha(X = X_train, time = time_train, status = status_train, 
                                  n_components = 5, return_all = TRUE, alpha = 0.9)
T_test <- X_test %*% result_drpcapls$A
colnames(T_test) <- paste0("comp_", 1:ncol(T_test))
T_test <- as.data.frame(T_test)
marker_lp_drpcapls <- predict(result_drpcapls$best_model, newdata = T_test, type = "lp")
marker_risk_drpcapls <- predict(result_drpcapls$best_model, newdata = T_test, type = "risk")

result_drpcapls$results

timepoints <- unique(time_test)
timepoints <- timepoints[order(timepoints)]
AUCs_drpcapls <- numeric(length(timepoints))

for (i in seq_along(timepoints)) {
  timepoint <- timepoints[i]
  AUCs_drpcapls[i] <- survivalROC(Stime = time_test, status = status_test, marker = marker_lp_drpcapls, 
                                       method = "KM", predict.time = timepoint)$AUC
}

plot(
  timepoints, AUCs_drpcapls, 
  type = "l",        # 用线绘制
  lty = 2,           # 虚线
  lwd = 2,           # 线宽
  col = "blue",
  xlab = "Time (months)",
  ylab = "AUC",
  main = "Time-dependent AUC curve (Single Model)"
)
abline(h = 0.5, lty = 2, col = "gray")  # 可选：绘制 AUC=0.5 的参考线

# partial cox regression
result_pcox <- partial_cox(X = X_train, time = time_train, status = status_train, 
                           n_components = 5)
betas_list <- result_pcox$betas_list
w_list <- result_pcox$w_list

T_test <- transform_partial_cox(X_test, betas_list, w_list)
colnames(T_test) <- paste0("T", 1:ncol(T_test))
T_test <- as.data.frame(T_test)
marker_lp_pcox <- predict(result_pcox$final_model, newdata = T_test, type = "lp")
marker_risk_pcox <- predict(result_pcox$final_model, newdata = T_test, type = "risk")

timepoints <- unique(time_test)
timepoints <- timepoints[order(timepoints)]
AUCs_pcox <- numeric(length(timepoints))

for (i in seq_along(timepoints)) {
  timepoint <- timepoints[i]
  AUCs_pcox[i] <- survivalROC(Stime = time_test, status = status_test, marker = marker_lp_pcox, 
                                       method = "KM", predict.time = timepoint)$AUC
}

plot(
  timepoints, AUCs_pcox, 
  type = "l",        # 用线绘制
  lty = 2,           # 虚线
  lwd = 2,           # 线宽
  col = "blue",
  xlab = "Time (months)",
  ylab = "AUC",
  main = "Time-dependent AUC curve (Single Model)"
)
abline(h = 0.5, lty = 2, col = "gray") 

# partial cox pca
result_pcox_pca <- partial_cox_pca_select_alpha(X = X_train, time = time_train, status = status_train, 
                                   n_components = 5)
betas_list <- result_pcox_pca$betas_list
weights_list <- result_pcox_pca$weights_list

T_test <- transform_partial_cox_pca(X_test, betas_list, weights_list)
colnames(T_test) <- paste0("T", 1:ncol(T_test))
T_test <- as.data.frame(T_test)

marker_lp_pcox_pca <- predict(result_pcox_pca$best_model, newdata = T_test, type = "lp")
marker_risk_pcox_pca <- predict(result_pcox_pca$best_model, newdata = T_test, type = "risk")

timepoints <- unique(time_test)
timepoints <- timepoints[order(timepoints)]
AUCs_pcox_pca <- numeric(length(timepoints))

for (i in seq_along(timepoints)) {
  timepoint <- timepoints[i]
  AUCs_pcox_pca[i] <- survivalROC(Stime = time_test, status = status_test, marker = marker_lp_pcox_pca, 
                                       method = "KM", predict.time = timepoint)$AUC
}

plot(
  timepoints, AUCs_pcox_pca, 
  type = "l",        # 用线绘制
  lty = 2,           # 虚线
  lwd = 2,           # 线宽
  col = "blue",
  xlab = "Time (months)",
  ylab = "AUC",
  main = "Time-dependent AUC curve (Single Model)"
)
abline(h = 0.5, lty = 2, col = "gray")  # 可选：绘制 AUC=0.5 的参考线

# 将四个数据画在一张图上
op <- par(mar = c(5, 4, 4, 8) + 0.1)
plot(timepoints, AUCs, type = "l", lty = 1, lwd = 2, col = "blue",
     ylim = c(0,1),  # 确保纵轴从0到1
     xlab = "Time (months)", ylab = "AUC",
     main = "Time-dependent AUC for Four Models")
lines(timepoints, AUCs_drpcapls, type = "l", lty = 2, lwd = 2, col = "red")
lines(timepoints, AUCs_pcox, type = "l", lty = 3, lwd = 2, col = "green")
lines(timepoints, AUCs_pcox_pca, type = "l", lty = 4, lwd = 2, col = "purple")
abline(h = 0.5, lty = 2, col = "gray")  # 可选：绘制 AUC=0.5 的参考线
legend("topright",
       inset = c(-0.7, 0),  # 负值越大, 图例越往右移
       xpd = TRUE,           # 允许在图外绘制
       legend = c("Cox", "DR-PCAPLS", "Partial Cox", "Partial Cox PCA"),
       col = c("blue","red","green","purple"),
       lty = c(1,2,3,4),
       lwd = 2,
       cex = 0.9,            # 可调小图例文字
       seg.len = 1,          # 缩短图例线段
       bty = "n")            # 去掉边框



# === predict time-dependent AUC as criteria ===
# train test split
set.seed(123)  # 为了结果可重复，建议设置随机种子
n <- nrow(X)          # 样本总数
train_size <- round(0.7 * n)     # 训练集数量

# 随机抽取训练集索引
train_index <- sample(seq_len(n), size = train_size)

# 训练集
pdata_train  <- pdata_beer_clean[train_index, ]
time_train   <- time_clean[train_index]
status_train <- status_clean[train_index]
X_train      <- X[train_index, ]

# 测试集
pdata_test  <- pdata_beer_clean[-train_index, ]
time_test   <- time_clean[-train_index]
status_test <- status_clean[-train_index]
X_test      <- X[-train_index, ]

dim(X_train)
dim(X_test)

# T
timepoints <- unique(time_test)
timepoints <- timepoints[order(timepoints)]

# DRPLS
result_drpls <- val_DRPLS(X_train, time_train, status_train, 
                          X_test, time_test, status_test, ncomp_candidates = 1:5, 
                          auc_time_grid = timepoints)
best_AUCs_drpls <- result_drpls$best_auc_curve

# DRPCAPLS
result_drpcapls <- val_DRPCAPLS(X_train, time_train, status_train, 
                                X_test, time_test, status_test, ncomp_candidates = 1:10, 
                                alpha_candidates = seq(0, 1, 0.1),
                                auc_time_grid = timepoints)
best_AUCs_drpcapls <- result_drpcapls$best_auc_curve

# partial cox
result_partialcox <- val_partial_cox(X_train, time_train, status_train, 
                                     X_test, time_test, status_test, ncomp_candidates = 2:10, 
                                     auc_time_grid = timepoints)
best_AUCs_pcox <- result_partialcox$best_auc_curve

# partial cox pca
result_pcox_pca <- val_partial_cox_pca(X_train, time_train, status_train, 
                                       X_test, time_test, status_test, ncomp_candidates = 2:10, 
                                       alpha_candidates = seq(0, 1, 0.1),
                                       auc_time_grid = timepoints)
best_AUCs_pcox_pca <- result_pcox_pca$best_auc_curve

# plot
op <- par(mar = c(5, 4, 4, 8) + 0.1)
plot(timepoints, best_AUCs_drpls, type = "l", lty = 1, lwd = 2, col = "blue",
     ylim = c(0,1),  # 确保纵轴从0到1
     xlab = "Time", ylab = "AUC",
     main = "Time-dependent AUC for Four Models")
lines(timepoints, best_AUCs_drpcapls, type = "l", lty = 2, lwd = 2, col = "red")
lines(timepoints, best_AUCs_pcox, type = "l", lty = 3, lwd = 2, col = "green")
lines(timepoints, best_AUCs_pcox_pca, type = "l", lty = 4, lwd = 2, col = "purple")
abline(h = 0.5, lty = 2, col = "gray")  # 可选：绘制 AUC=0.5 的参考线
legend("topright",
       inset = c(-0.8, 0),  # 负值越大, 图例越往右移
       xpd = TRUE,           # 允许在图外绘制
       legend = c("DRCox", "DR-PCAPLS", "Partial Cox", "Partial Cox PCA"),
       col = c("blue","red","green","purple"),
       lty = c(1,2,3,4),
       lwd = 2,
       cex = 0.9,            # 可调小图例文字
       seg.len = 1,          # 缩短图例线段
       bty = "n")  





