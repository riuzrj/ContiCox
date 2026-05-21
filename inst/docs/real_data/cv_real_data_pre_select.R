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

# DRPLS
beer_auc_curves_list_drpls <- val_DRPLS_cv(X, time_clean, status_clean, k = 5,
                                       ncomp_candidates = 1:5, auc_time_grid = NULL,
                                       use_preselection = FALSE, p_thresh = 0.05,
                                       auc_method = "NNE")
# DRPCAPLS
beer_auc_curves_list_drpcapls <- val_DRPCAPLS_cv(X, time_clean, status_clean, k = 5,
                                          ncomp_candidates = 1:10, 
                                          alpha_candidates = seq(0, 1, by = 0.1), 
                                          auc_time_grid = NULL,
                                          use_preselection = FALSE, p_thresh = 0.05,
                                          auc_method = "NNE")

# partial cox
beer_auc_curves_list_pcox <- val_partial_cox_cv(X, time_clean, status_clean, k = 5,
                                          ncomp_candidates = 2:10, 
                                          auc_time_grid = NULL,
                                          use_preselection = TRUE, p_thresh = 0.05,
                                          auc_method = "NNE")

# partial cox with pca
beer_auc_curves_list_pca_pcox <- val_partial_cox_pca_cv(X, time_clean, status_clean, k = 5,
                                          ncomp_candidates = 2:8, 
                                          auc_time_grid = NULL,
                                          alpha_candidates = seq(0, 1, by = 0.1),
                                          use_preselection = TRUE, p_thresh = 0.05,
                                          auc_method = "NNE")
save(beer_auc_curves_list_pca_pcox, file = "beer_auc_curves_list_pca_pcox.RData")
load(file = "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls/beer_auc_curves_list_pca_pcox.RData")

# PCA weight partial cox
beer_auc_curves_list_pca_weight_partial_cox <- val_partial_cox_pca_weight_cv(X, time_clean, status_clean, k = 5,
                                          ncomp_candidates = 2:10, 
                                          auc_time_grid = NULL,
                                          use_preselection = TRUE, p_thresh = 0.05,
                                          auc_method = "NNE")

# PC - PCR
beer_auc_curves_list_prepca_partial_cox <- val_prepca_partial_cox_cv(X, time_clean, status_clean, k = 5,
                                          ncomp_candidates = 2:9, 
                                          auc_time_grid = NULL,
                                          use_preselection = TRUE, p_thresh = 0.05,
                                          auc_method = "NNE", use_pcr = TRUE)

# PCAPLS - PCR
beer_auc_curves_list_pcapls_partial_cox <- val_pcapls_partial_cox_cv(X, time_clean, status_clean, k = 5,
                                          ncomp_candidates = 2:6, 
                                          auc_time_grid = NULL,
                                          use_preselection = TRUE, p_thresh = 0.05,
                                          auc_method = "NNE", alpha_candidates = seq(0, 1, by = 0.1))


# plot
par(mar = c(5, 4, 4, 6), xpd = TRUE)
time_points = beer_auc_curves_list_drpls$auc_time_grid
plot(time_points, beer_auc_curves_list_drpls$best_auc_curve, type = "l", col = "blue",
     xlab = "Time (months)", ylab = "AUC", ylim = c(0.6, 1), lwd = 2,
     main = "AUC curves for different methods")
lines(time_points, beer_auc_curves_list_drpcapls$best_auc_curve, col = "red", lwd = 2)
lines(time_points, beer_auc_curves_list_pcox$best_auc_curve, col = "green", lwd = 2)
lines(time_points, beer_auc_curves_list_pca_pcox$best_auc_curve, col = "purple", lwd = 2)
lines(time_points, beer_auc_curves_list_pca_weight_partial_cox$best_auc_curve, col = "black", lwd = 2)
lines(time_points, beer_auc_curves_list_prepca_partial_cox$best_auc_curve, col = "orange", lwd = 2)
lines(time_points, beer_auc_curves_list_pcapls_partial_cox$best_auc_curve, col = "gray", lwd = 2)
legend("topright", inset = c(-0.4, 0),legend = c("DRPLS", "DRPCAPLS", "Partial Cox", "PCA Partial Cox","wPCA Partial Cox",
                                                 "PrePCA-Partial_cox", "PCAPLS-Partial_cox"),
       col = c("red", "blue", "green", "purple", "black","orange", "gray"), lwd = 1, cex = 0.5, bty = "n")
# 只画前四个
legend("topright", inset = c(-0.4, 0),legend = c("DRPLS", "DRPCAPLS", "Partial Cox", "PCA Partial Cox"),
       col = c("red", "blue", "green", "purple"), lwd = 1, cex = 0.5, bty = "n")



# breast data
options(
  download.file.method.GEOquery = "curl",
  timeout = 600
)
library(GEOquery)
gse_breast <- getGEO("GSE2034", GSEMatrix=TRUE)
# get features
exprSet_breast <- exprs(gse_breast[[1]])
X <- t(exprSet_breast)
# get phenotype
pdata_breast <- pData(gse_breast[[1]])
# 1. 提取随访时间 (time)
clinical <- read.csv('/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls/breast.csv')
rownames(X) <- trimws(as.character(rownames(X)))
clinical$GEO.asscession.number <- trimws(as.character(clinical$GEO.asscession.number))
idx <- match(rownames(X), clinical$GEO.asscession.number)
clinical_order <- clinical[idx, ]
all(rownames(X) == clinical_order$GEO.asscession.number)
time_clean <- clinical_order$time.to.relapse.or.last.follow.up..months.
status_clean <- clinical_order$relapse..1.True.
sum(is.na(X))
X_center <- colMeans(X)
X <- sweep(X, 2, X_center)
X <- as.matrix(X)

# DRPLS
breast_auc_curves_list_drpls <- val_DRPLS_cv(X, time_clean, status_clean, k = 5,
                                       ncomp_candidates = 1:10, auc_time_grid = NULL,
                                       use_preselection = FALSE, p_thresh = 0.05,
                                       auc_method = "KM")

# DRPCAPLS
breast_auc_curves_list_drpcapls <- val_DRPCAPLS_cv(X, time_clean, status_clean, k = 5,
                                          ncomp_candidates = 1:10, 
                                          alpha_candidates = seq(0, 1, by = 0.1), 
                                          auc_time_grid = NULL,
                                          use_preselection = FALSE, p_thresh = 0.05,
                                          auc_method = "KM")

# partial cox
breast_auc_curves_list_pcox <- val_partial_cox_cv(X, time_clean, status_clean, k = 5,
                                          ncomp_candidates = 2:10, 
                                          auc_time_grid = NULL,
                                          use_preselection = FALSE, p_thresh = 0.05,
                                          auc_method = "KM")

# partial cox with pca
breast_auc_curves_list_pca_pcox <- val_partial_cox_pca_cv(X, time_clean, status_clean, k = 5,
                                          ncomp_candidates = 2:10, 
                                          auc_time_grid = NULL,
                                          use_preselection = TRUE, p_thresh = 0.05,
                                          auc_method = "NNE")
save(breast_auc_curves_list_pca_pcox, file = "breast_auc_curves_list_pca_pcox.RData")

# PCA weight partial cox
breast_auc_curves_list_pca_weight_partial_cox <- val_partial_cox_pca_weight_cv(X, time_clean, status_clean, k = 5,
                                          ncomp_candidates = 2:10, 
                                          auc_time_grid = NULL,
                                          use_preselection = TRUE, p_thresh = 0.05,
                                          auc_method = "NNE")

# PC - PCR
breast_auc_curves_list_prepca_partial_cox <- val_prepca_partial_cox_cv(X, time_clean, status_clean, k = 5,
                                          ncomp_candidates = 2:10, 
                                          auc_time_grid = NULL,
                                          use_preselection = TRUE, p_thresh = 0.05,
                                          auc_method = "NNE", use_pcr = TRUE)

# PCAPLS - PCR
breast_auc_curves_list_pcapls_partial_cox <- val_pcapls_partial_cox_cv(X, time_clean, status_clean, k = 5,
                                          ncomp_candidates = 2:10, 
                                          auc_time_grid = NULL,
                                          use_preselection = TRUE, p_thresh = 0.05,
                                          auc_method = "NNE", alpha_candidates = seq(0, 1, by = 0.1))

# plot
par(mar = c(5, 4, 4, 6), xpd = TRUE)
time_points = breast_auc_curves_list_drpls$auc_time_grid
plot(time_points, breast_auc_curves_list_drpls$best_auc_curve, type = "l", col = "blue",
     xlab = "Time (months)", ylab = "AUC", ylim = c(0.65, 1), lwd = 2,
     main = "AUC curves for different methods")
lines(time_points, breast_auc_curves_list_drpcapls$best_auc_curve, col = "red", lwd = 2)
lines(time_points, breast_auc_curves_list_pcox$best_auc_curve, col = "green", lwd = 2)
lines(time_points, breast_auc_curves_list_pca_pcox$best_auc_curve, col = "purple", lwd = 2)
lines(time_points, breast_auc_curves_list_pca_weight_partial_cox$best_auc_curve, col = "black", lwd = 2)
lines(time_points, breast_auc_curves_list_prepca_partial_cox$best_auc_curve, col = "orange", lwd = 2)
lines(time_points, breast_auc_curves_list_pcapls_partial_cox$best_auc_curve, col = "gray", lwd = 2)
legend("topright", inset = c(-0.4, 0),legend = c("DRPLS", "DRPCAPLS", "Partial Cox", "PCA Partial Cox", "wPCA Partial Cox",
                                                 "PrePCA-Partial_cox", "PCAPLS-Partial_cox"),
       col = c("red", "blue", "green", "purple", "black", "orange", "gray"), lwd = 1, cex = 0.5, bty = "n")



# AML
gse_AML <- getGEO("GSE12417", GSEMatrix=TRUE)

# get features
exprSet_AML <- exprs(gse_AML[[1]])
X <- t(exprSet_AML)

# get phenotype
pdata_AML <- pData(gse_AML[[1]])
# 1. 提取随访时间 (time)
time_AML <- as.numeric(sub(".*OS = (\\d+) days.*", "\\1", pdata_AML$characteristics_ch1))
# 2. 提取生存状态 (status)
status_AML <- as.numeric(sub(".*status \\(0=alive/1=dead\\): (\\d+).*", "\\1", pdata_AML$characteristics_ch1))

sum(is.na(X))
X_center <- colMeans(X)
X <- sweep(X, 2, X_center)
X <- as.matrix(X)

# DRPLS
AML_auc_curves_list_drpls <- val_DRPLS_cv(X, time_AML, status_AML, k = 5,
                                       ncomp_candidates = 1:6, auc_time_grid = NULL,
                                       use_preselection = TRUE, p_thresh = 0.05,
                                       auc_method = "NNE")
# DRPCAPLS
AML_auc_curves_list_drpcapls <- val_DRPCAPLS_cv(X, time_AML, status_AML, k = 5,
                                          ncomp_candidates = 1:10, 
                                          alpha_candidates = seq(0, 1, by = 0.1), 
                                          auc_time_grid = NULL,
                                          use_preselection = TRUE, p_thresh = 0.05,
                                          auc_method = "NNE")
# partial cox
AML_auc_curves_list_pcox <- val_partial_cox_cv(X, time_AML, status_AML, k = 5,
                                          ncomp_candidates = 2:8, 
                                          auc_time_grid = NULL,
                                          use_preselection = TRUE, p_thresh = 0.05,
                                          auc_method = "NNE")
# partial cox with pca
AML_auc_curves_list_pca_pcox <- val_partial_cox_pca_cv(X, time_AML, status_AML, k = 5,
                                          ncomp_candidates = 2:9, 
                                          auc_time_grid = NULL,
                                          use_preselection = TRUE, p_thresh = 0.05,
                                          auc_method = "NNE")
save(AML_auc_curves_list_pca_pcox, file = "AML_auc_curves_list_pca_pcox.RData")
load(file = "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls/AML_auc_curves_list_pca_pcox.RData")

# PCA weight partial cox
AML_auc_curves_list_pca_weight_partial_cox <- val_partial_cox_pca_weight_cv(X, time_AML, status_AML, k = 5,
                                          ncomp_candidates = 2:8, 
                                          auc_time_grid = NULL,
                                          use_preselection = TRUE, p_thresh = 0.05,
                                          auc_method = "NNE")
# PC - PCR
AML_auc_curves_list_prepca_partial_cox <- val_prepca_partial_cox_cv(X, time_AML, status_AML, k = 5,
                                          ncomp_candidates = 2:9, 
                                          auc_time_grid = NULL,
                                          use_preselection = TRUE, p_thresh = 0.05,
                                          auc_method = "NNE", use_pcr = TRUE)

# PCAPLS - PCR
AML_auc_curves_list_pcapls_partial_cox <- val_pcapls_partial_cox_cv(X, time_AML, status_AML, k = 5,
                                          ncomp_candidates = 2:6, 
                                          auc_time_grid = NULL,
                                          use_preselection = TRUE, p_thresh = 0.05,
                                          auc_method = "NNE", alpha_candidates = seq(0, 1, by = 0.1))

# plot
par(mar = c(5, 4, 4, 6), xpd = TRUE)
time_points = AML_auc_curves_list_drpls$auc_time_grid
plot(time_points, AML_auc_curves_list_drpls$best_auc_curve, type = "l", col = "blue",
     xlab = "Time (days)", ylab = "AUC", ylim = c(0.5, 1), lwd = 2,
     main = "AUC curves for different methods")
lines(time_points, AML_auc_curves_list_drpcapls$best_auc_curve, col = "red", lwd = 2)
lines(time_points, AML_auc_curves_list_pcox$best_auc_curve, col = "green", lwd = 2)
lines(time_points, AML_auc_curves_list_pca_pcox$best_auc_curve, col = "purple", lwd = 2)
lines(time_points, AML_auc_curves_list_pca_weight_partial_cox$best_auc_curve, col = "black", lwd = 2)
lines(time_points, AML_auc_curves_list_prepca_partial_cox$best_auc_curve, col = "orange", lwd = 2)
lines(time_points, AML_auc_curves_list_pcapls_partial_cox$best_auc_curve, col = "gray", lwd =1)
legend("topright", inset = c(-0.4, 0),legend = c("DRPLS", "DRPCAPLS", "Partial Cox", "PCA Partial Cox", "wPCA Partial Cox",
                                                 "PrePCA-Partial_cox", "PCAPLS-Partial_cox"),
       col = c("red", "blue", "green", "purple", "black","orange", "gray"), lwd = 1, cex = 0.5, bty = "n")




