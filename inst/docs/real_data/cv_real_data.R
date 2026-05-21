library(survivalROC)
library(GEOquery)
library(caret)

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

# breast cancer
gse_breast <- getGEO("GSE2034", GSEMatrix=FALSE)
# get features
exprSet_breast <- exprs(gse_breast[[1]])
# get phenotype
pdata_breast <- pData(gse_breast[[1]])


# ====== cross validation ======
set.seed(123)  

# 创建 5 折交叉验证的分组
folds <- createFolds(1:nrow(X), k = 5)
auc_curves_list_drpls <- list()    # 用于存储每一折的 AUC 曲线

# 检查数据
test_index <- folds[[1]]
train_index <- setdiff(1:nrow(X), test_index)
time_train  <- time_clean[train_index]
status_train <- status_clean[train_index]
X_train     <- X[train_index, ]
time_test  <- time_clean[test_index]
status_test <- status_clean[test_index]
X_test     <- X[test_index, ]
summary(time_test)
hist(time_test, breaks = 20, main = "Histogram of Time Test", xlab = "Time")

# DRPLS
for(i in 1:5) {
  # 定义测试集和训练集索引
  test_index <- folds[[i]]
  train_index <- setdiff(1:nrow(X), test_index)
  
  # 构建训练集
  # pdata_train <- pdata_beer_clean[train_index, ]
  time_train  <- time_clean[train_index]
  status_train <- status_clean[train_index]
  X_train     <- X[train_index, ]
  
  # 构建测试集
  # pdata_test <- pdata_beer_clean[test_index, ]
  time_test  <- time_clean[test_index]
  status_test <- status_clean[test_index]
  X_test     <- X[test_index, ]
  
  # 定义时间点：取测试集中所有唯一时间点，并排序
  timepoints <- sort(unique(time_test))
  
  # DRPLS 验证
  result_drpls <- val_DRPLS(X_train, time_train, status_train, 
                            X_test, time_test, status_test, 
                            ncomp_candidates = 1:5, 
                            auc_time_grid = timepoints)
  
  # if(any(is.na(result_drpls$best_auc_curve))) {
  #   warning(sprintf("第 %d 折 DRPLS 结果中含有 NA，跳过该折", i))
  #   next
  # }
  
  auc_curves_list_drpls[[length(auc_curves_list_drpls) + 1]] <- 
    result_drpls$best_auc_curve
}


auc_matrix_drpls <- do.call(rbind, auc_curves_list_drpls)
mean_auc_curve_drpls <- colMeans(auc_matrix_drpls, na.rm = TRUE)
print(mean_auc_curve_drpls)

# DRPPCAPLS
auc_curves_list_drpcapls <- list()

for(i in 1:5) {
  print(i)
  test_index <- folds[[i]]
  train_index <- setdiff(1:nrow(X), test_index)
  
  #pdata_train <- pdata_beer_clean[train_index, ]
  time_train  <- time_clean[train_index]
  status_train <- status_clean[train_index]
  X_train     <- X[train_index, ]
  
  #pdata_test <- pdata_beer_clean[test_index, ]
  time_test  <- time_clean[test_index]
  status_test <- status_clean[test_index]
  X_test     <- X[test_index, ]
  
  timepoints <- sort(unique(time_test))
  
  result_drpcapls <- val_DRPCAPLS(X_train, time_train, status_train, 
                                  X_test, time_test, status_test, 
                                  ncomp_candidates = 1:10, 
                                  alpha_candidates = seq(0, 1, 0.1),
                                  auc_time_grid = timepoints)
  
  # if(any(is.na(result_drpcapls$best_auc_curve))) {
  #   warning(sprintf("第 %d 折 DRPCAPLS 结果中含有 NA，跳过该折", i))
  #   next
  # }
  
  auc_curves_list_drpcapls[[length(auc_curves_list_drpcapls) + 1]] <- result_drpcapls$best_auc_curve
}

auc_matrix_drpcapls <- do.call(rbind, auc_curves_list_drpcapls)
mean_auc_curve_drpcapls <- colMeans(auc_matrix_drpcapls, na.rm = TRUE)
print(mean_auc_curve_drpcapls)

# partial cox 
auc_curves_list_pcox <- list()

for(i in 1:5) {
  test_index <- folds[[i]]
  train_index <- setdiff(1:nrow(X), test_index)
  
  #pdata_train <- pdata_beer_clean[train_index, ]
  time_train  <- time_clean[train_index]
  status_train <- status_clean[train_index]
  X_train     <- X[train_index, ]
  
  #pdata_test <- pdata_beer_clean[test_index, ]
  time_test  <- time_clean[test_index]
  status_test <- status_clean[test_index]
  X_test     <- X[test_index, ]
  
  timepoints <- sort(unique(time_test))
  
  result_pcox <- val_partial_cox(X_train, time_train, status_train, 
                                 X_test, time_test, status_test, 
                                 ncomp_candidates = 2:10, 
                                 auc_time_grid = timepoints)
  
  # if(any(is.na(result_pcox$best_auc_curve))) {
  #   warning(sprintf("第 %d 折 Partial Cox 结果中含有 NA，跳过该折", i))
  #   next
  # }
  
  auc_curves_list_pcox[[length(auc_curves_list_pcox) + 1]] <- result_pcox$best_auc_curve
}

auc_matrix_pcox <- do.call(rbind, auc_curves_list_pcox)
mean_auc_curve_pcox <- colMeans(auc_matrix_pcox, na.rm = TRUE)
print(mean_auc_curve_pcox)


# partial cox pca
auc_curves_list_pcox_pca <- list()

for(i in 1:5) {
  print(i)
  test_index <- folds[[i]]
  train_index <- setdiff(1:nrow(X), test_index)
  
  #pdata_train <- pdata_beer_clean[train_index, ]
  time_train  <- time_clean[train_index]
  status_train <- status_clean[train_index]
  X_train     <- X[train_index, ]
  
  #pdata_test <- pdata_beer_clean[test_index, ]
  time_test  <- time_clean[test_index]
  status_test <- status_clean[test_index]
  X_test     <- X[test_index, ]
  
  timepoints <- sort(unique(time_test))
  
  result_pcox_pca <- val_partial_cox_pca(X_train, time_train, status_train, 
                                         X_test, time_test, status_test, 
                                         ncomp_candidates = 2:10, 
                                         alpha_candidates = seq(0, 1, 0.1),
                                         auc_time_grid = timepoints)
  
  # if(any(is.na(result_pcox_pca$best_auc_curve))) {
  #   wasrning(sprintf("第 %d 折 Partial Cox PCA 结果中含有 NA，跳过该折", i))
  #   next
  # }
  
  auc_curves_list_pcox_pca[[length(auc_curves_list_pcox_pca) + 1]] <- result_pcox_pca$best_auc_curve
}

auc_matrix_pcox_pca <- do.call(rbind, auc_curves_list_pcox_pca)
mean_auc_curve_pcox_pca <- colMeans(auc_matrix_pcox_pca, na.rm = TRUE)
print(mean_auc_curve_pcox_pca)

# plot 
par(mar = c(5, 4, 4, 6), xpd = TRUE)
timepoints <- c(1:18)
plot(timepoints, mean_auc_curve_drpls, type = "l", col = "red", lwd = 2,
     xlab = "Time", ylab = "Mean AUC", ylim = c(0.4, 1), main = "AUC Curves")
lines(timepoints, mean_auc_curve_drpcapls, col = "blue", lwd = 2)
lines(timepoints, mean_auc_curve_pcox, col = "green", lwd = 2)
lines(timepoints, mean_auc_curve_pcox_pca, col = "purple", lwd = 2)
legend("topright",
       inset = c(-0.35, 0),  # 负的 inset 可以把图例移到绘图区域外
       legend = c("DRPLS", "DRPCAPLS", "Partial Cox", "Partial Cox PCA"),
       col = c("red", "blue", "green", "purple"),
       lwd = 2,
       cex = 0.6)  # 调小图例文字大小



