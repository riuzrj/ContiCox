library(MASS) # 用于生成正态分布数据

# eigengene
simulate_one <- function(
    n_samples = 100,
    n_total_genes = 1000,
    n_modules = 4,
    n_genes_per_mod = 25,
    n_signal_modules = 2,
    r_min = 0.7,
    signal_strength = 1.0,
    signal_signs = NULL,
    signal_coef = c("module_average", "gene_unit"),
    standardize_lp = TRUE,
    base_rate = 0.1,
    censoring = 0.4
) {
  stopifnot(n_modules * n_genes_per_mod <= n_total_genes)
  stopifnot(n_signal_modules <= n_modules)
  stopifnot(r_min > 0, r_min < 1)
  stopifnot(censoring > 0, censoring < 1)
  signal_coef <- match.arg(signal_coef)
  
  if (is.null(signal_signs)) {
    signal_signs <- rep(1, n_signal_modules)
    if (n_signal_modules >= 2) {
      signal_signs[seq(2, n_signal_modules, by = 2)] <- -1
    }
  }
  if (length(signal_signs) != n_signal_modules) {
    stop("length(signal_signs) must equal n_signal_modules.")
  }
  
  n_noise_genes <- n_total_genes - (n_modules * n_genes_per_mod)
  datExpr <- matrix(NA, nrow = n_samples, ncol = n_total_genes)
  colnames(datExpr) <- paste0("Gene_", 1:n_total_genes)
  
  true_seeds <- matrix(NA, nrow = n_samples, ncol = n_modules)
  colnames(true_seeds) <- paste0("ME", seq_len(n_modules))
  current_col <- 1
  
  for (I in 1:n_modules) {
    seed_vector <- rnorm(n_samples)
    true_seeds[, I] <- seed_vector
    var_seed <- var(seed_vector)
    
    for (k in 1:n_genes_per_mod) {
      r_kI <- 1 - (k / n_genes_per_mod) * (1 - r_min)
      
      epsilon_k <- rnorm(n_samples)
      var_eps <- var(epsilon_k)
      
      a_k <- sqrt((var_seed / var_eps) * ((1 / (r_kI^2)) - 1))
      x_k <- seed_vector + a_k * epsilon_k
      
      datExpr[, current_col] <- x_k
      current_col <- current_col + 1
    }
  }
  
  if (n_noise_genes > 0L) {
    noise_data <- matrix(rnorm(n_samples * n_noise_genes),
                         nrow = n_samples, ncol = n_noise_genes)
    datExpr[, current_col:n_total_genes] <- noise_data
  }
  
  beta_true_base <- rep(0, n_total_genes)
  coef_scale <- if (identical(signal_coef, "gene_unit")) 1 else 1 / n_genes_per_mod
  for (mod in seq_len(n_signal_modules)) {
    idx_start <- (mod - 1L) * n_genes_per_mod + 1L
    idx_end <- mod * n_genes_per_mod
    beta_true_base[idx_start:idx_end] <- signal_signs[mod] * coef_scale
  }
  names(beta_true_base) <- colnames(datExpr)
  
  raw_score <- as.numeric(datExpr %*% beta_true_base)
  raw_sd <- stats::sd(raw_score)
  if (isTRUE(standardize_lp) && is.finite(raw_sd) && raw_sd > 0) {
    raw_score <- as.numeric(scale(raw_score))
    beta_true <- (signal_strength / raw_sd) * beta_true_base
    lp <- signal_strength * raw_score
  } else {
    beta_true <- signal_strength * beta_true_base
    lp <- signal_strength * raw_score
  }
  names(beta_true) <- colnames(datExpr)
  lambda <- base_rate * exp(lp)
  T_true <- rexp(n_samples, rate = lambda)
  
  censor_rate <- base_rate * censoring / (1 - censoring)
  C_time <- rexp(n_samples, rate = censor_rate)
  
  Time <- pmin(T_true, C_time)
  Status <- as.integer(T_true <= C_time)
  
  list(
    X = datExpr,
    seeds = true_seeds,
    beta_true_base = beta_true_base,
    beta_true = beta_true,
    linear_predictor = lp,
    surv = data.frame(Time = Time, Status = Status),
    params = list(
      scenario = "eigengene_module",
      n_samples = n_samples,
      n_total_genes = n_total_genes,
      n_modules = n_modules,
      n_genes_per_mod = n_genes_per_mod,
      n_signal_modules = n_signal_modules,
      r_min = r_min,
      signal_strength = signal_strength,
      signal_signs = signal_signs,
      signal_coef = signal_coef,
      standardize_lp = standardize_lp,
      base_rate = base_rate,
      censoring = censoring
    )
  )
}

estimate_exponential_censor_rate <- function(event_time, target_censoring) {
  stopifnot(target_censoring > 0, target_censoring < 1)
  event_time <- event_time[is.finite(event_time) & event_time > 0]
  if (length(event_time) < 1L) {
    stop("event_time must contain positive finite values.")
  }
  objective <- function(rate) {
    mean(1 - exp(-rate * event_time)) - target_censoring
  }
  upper <- 1 / stats::median(event_time)
  if (!is.finite(upper) || upper <= 0) {
    upper <- 1
  }
  while (objective(upper) < 0) {
    upper <- upper * 2
  }
  stats::uniroot(objective, lower = 0, upper = upper)$root
}

simulate_block_gaussian <- function(
    n_samples = 100,
    n_total_genes = 1000,
    n_blocks = 6,
    genes_per_block = 25,
    n_signal_blocks = 2,
    rho = 0.9,
    signal_strength = 1.0,
    signal_signs = c(1, -1),
    base_rate = 0.1,
    censoring = 0.4
) {
  stopifnot(n_blocks * genes_per_block <= n_total_genes)
  stopifnot(n_signal_blocks <= n_blocks)
  stopifnot(rho > 0, rho < 1)
  stopifnot(length(signal_signs) == n_signal_blocks)
  
  X <- matrix(rnorm(n_samples * n_total_genes), nrow = n_samples)
  colnames(X) <- paste0("Gene_", seq_len(n_total_genes))
  block_scores <- matrix(NA_real_, nrow = n_samples, ncol = n_blocks)
  
  current_col <- 1L
  for (block in seq_len(n_blocks)) {
    z <- rnorm(n_samples)
    block_scores[, block] <- z
    for (k in seq_len(genes_per_block)) {
      X[, current_col] <- sqrt(rho) * z + sqrt(1 - rho) * rnorm(n_samples)
      current_col <- current_col + 1L
    }
  }
  
  X <- as.matrix(scale(X, center = TRUE, scale = FALSE))
  raw_score <- rowSums(sweep(
    block_scores[, seq_len(n_signal_blocks), drop = FALSE],
    2,
    signal_signs,
    `*`
  ))
  raw_score <- as.numeric(scale(raw_score))
  lp <- signal_strength * raw_score
  
  event_time <- -log(runif(n_samples)) / (base_rate * exp(lp))
  censor_rate <- estimate_exponential_censor_rate(event_time, censoring)
  censor_time <- rexp(n_samples, rate = censor_rate)
  
  list(
    X = X,
    latent = block_scores,
    linear_predictor = lp,
    surv = data.frame(
      Time = pmin(event_time, censor_time),
      Status = as.integer(event_time <= censor_time)
    ),
    params = list(
      scenario = "block_gaussian_opposite_signal",
      n_samples = n_samples,
      n_total_genes = n_total_genes,
      n_blocks = n_blocks,
      genes_per_block = genes_per_block,
      n_signal_blocks = n_signal_blocks,
      rho = rho,
      signal_strength = signal_strength,
      signal_signs = signal_signs,
      base_rate = base_rate,
      censoring = censoring,
      censor_rate = censor_rate
    )
  )
}

simulate_variance_signal <- function(
    n_samples = 80,
    n_total_genes = 3000,
    signal_genes = 8,
    signal_scale = 8.0,
    signal_rho = 0.95,
    decoy_blocks = 12,
    decoy_genes = 50,
    decoy_scale = 1.0,
    decoy_rho = 0.7,
    signal_strength = 0.7,
    base_rate = 0.1,
    censoring = 0.45
) {
  stopifnot(signal_genes + decoy_blocks * decoy_genes <= n_total_genes)
  stopifnot(signal_rho > 0, signal_rho < 1)
  stopifnot(decoy_rho > 0, decoy_rho < 1)
  
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
    latent = data.frame(signal = z_signal),
    linear_predictor = lp,
    surv = data.frame(
      Time = pmin(event_time, censor_time),
      Status = as.integer(event_time <= censor_time)
    ),
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
      base_rate = base_rate,
      censoring = censoring,
      censor_rate = censor_rate
    )
  )
}

set.seed(123)

R <- 100
eigengene_module_grid <- expand.grid(
  n_total_genes = c(1000L, 2000L),
  n_modules = c(4L, 6L),
  r_min = c(0.5, 0.7, 0.9),
  signal_strength = c(0.5, 1.0, 1.5),
  KEEP.OUT.ATTRS = FALSE
)

sim_list_eigengene <- replicate(
  R,
  simulate_one(
    n_samples = 100,
    n_total_genes = 1000,
    n_modules = 6,
    n_genes_per_mod = 25,
    n_signal_modules = 2,
    r_min = 0.9,
    signal_strength = 1.0,
    signal_signs = c(1, -1),
    signal_coef = "module_average",
    standardize_lp = TRUE,
    censoring = 0.4
  ),
  simplify = FALSE
)
sim_list <- sim_list_eigengene


# cluster
simulate_cluster_complex <- function(
    n_samples      = 100,
    n_total_genes  = 1000,
    n_signal_block1 = 25,    # 信号块 1：与簇 + 生存相关
    n_signal_block2 = 25,    # 信号块 2：只与簇相关（生存不直接用）
    n_clusters     = 2,
    r_min          = 0.8,    # 块内相关性下界（越大共线性越强）
    delta1         = 1.0,    # 簇在信号块 1 上的均值差异尺度
    delta2         = 0.8,    # 簇在信号块 2 上的均值差异尺度
    beta_cluster   = NULL,   # 簇对 log-hazard 的系数向量，长度 = n_clusters
    beta_cont      = 0.8,    # 连续风险评分的系数
    base_rate      = 0.1,    # 基础事件率
    censor_rate    = 0.4     # 目标删失比例（近似）
) {
  stopifnot(n_signal_block1 + n_signal_block2 <= n_total_genes)
  
  # ---- 0. 病人簇标签：大致均匀分配 ----
  cluster <- factor(rep(seq_len(n_clusters), length.out = n_samples))
  cluster_id <- as.integer(cluster)
  
  # 若未指定 beta_cluster，自动生成一个从 -1 到 1 的梯度
  if (is.null(beta_cluster)) {
    beta_cluster <- seq(-1, 1, length.out = n_clusters)
  } else {
    if (length(beta_cluster) != n_clusters)
      stop("length(beta_cluster) 必须等于 n_clusters")
  }
  
  # ---- 1. 预分配表达矩阵 ----
  X <- matrix(NA, nrow = n_samples, ncol = n_total_genes)
  colnames(X) <- paste0("Gene_", seq_len(n_total_genes))
  
  # ---- 2. 为两个信号块生成簇特异的 latent seed ----
  # block 1：与簇 + 生存相关
  # block 2：与簇相关，但不直接进入生存模型
  
  # 每个簇的均值偏移
  mu1 <- seq(-delta1 / 2, delta1 / 2, length.out = n_clusters)
  mu2 <- seq(-delta2 / 2, delta2 / 2, length.out = n_clusters)
  
  # 病人级 latent seed
  z1 <- rnorm(n_samples)  # 基础
  z2 <- rnorm(n_samples)
  
  for (g in seq_len(n_clusters)) {
    idx_g <- which(cluster_id == g)
    z1[idx_g] <- z1[idx_g] + mu1[g]
    z2[idx_g] <- z2[idx_g] + mu2[g]
  }
  
  # ---- 3. 按 seed + 噪声 的方式生成块内强相关基因 ----
  make_block <- function(seed_vec, n_genes_in_block, start_col) {
    n <- length(seed_vec)
    var_seed <- var(seed_vec)
    cur_col <- start_col
    
    for (k in seq_len(n_genes_in_block)) {
      # 随 k 变化的相关系数，从 1 ~ r_min 之间
      r_k <- 1 - (k / n_genes_in_block) * (1 - r_min)
      
      eps_k <- rnorm(n)
      var_eps <- var(eps_k)
      
      a_k <- sqrt((var_seed / var_eps) * ((1 / (r_k^2)) - 1))
      x_k <- seed_vec + a_k * eps_k
      
      X[, cur_col] <<- x_k
      cur_col <- cur_col + 1
    }
    return(cur_col)
  }
  
  current_col <- 1
  # Block 1：生存相关块
  current_col <- make_block(z1, n_signal_block1, current_col)
  # Block 2：只与簇相关块
  current_col <- make_block(z2, n_signal_block2, current_col)
  
  # ---- 4. 噪声基因：完全无关，任意正态 ----
  if (current_col <= n_total_genes) {
    n_noise <- n_total_genes - current_col + 1
    X[, current_col:n_total_genes] <- matrix(
      rnorm(n_samples * n_noise),
      nrow = n_samples, ncol = n_noise
    )
  }
  
  # ---- 5. 构造生存时间 ----
  # 连续风险评分：用 block1 的 latent seed（或其标准化）
  z1_scaled <- as.numeric(scale(z1))
  
  lp <- beta_cluster[cluster_id] + beta_cont * z1_scaled
  lambda_event <- base_rate * exp(lp)
  T_true <- rexp(n_samples, rate = lambda_event)
  
  # 删失率近似控制：用指数删失率 λ_C ≈ base_rate * c/(1-c)
  lambda_c <- base_rate * censor_rate / (1 - censor_rate)
  C_time <- rexp(n_samples, rate = lambda_c)
  
  Time   <- pmin(T_true, C_time)
  Status <- as.integer(T_true <= C_time)
  
  list(
    X    = X,
    surv = data.frame(Time = Time, Status = Status),
    cluster = cluster,          # 病人簇
    latent = list(z1 = z1, z2 = z2),
    params = list(
      n_samples       = n_samples,
      n_total_genes   = n_total_genes,
      n_signal_block1 = n_signal_block1,
      n_signal_block2 = n_signal_block2,
      n_clusters      = n_clusters,
      r_min           = r_min,
      delta1          = delta1,
      delta2          = delta2,
      beta_cluster    = beta_cluster,
      beta_cont       = beta_cont,
      base_rate       = base_rate,
      censor_rate     = censor_rate
    )
  )
}

R <- 100
set.seed(2026)
sim_list <- replicate(
  R,
  simulate_cluster_complex(),
  simplify = FALSE
)




#factorial
simulate_factorial <- function(
    n_samples      = 100,
    n_total_genes  = 1000,
    n_groups       = 4,
    genes_per_group = 25,      # 每个 group 的“结构基因”数；4*25=100
    rho_within     = 0.7,      # 同一 group 内的相关系数
    n_signal_groups = 2,       # 前几个 group 与生存相关（这里是前 2 个 group）
    beta_signal    = 0.8,      # 信号基因对 log-hazard 的系数强度
    base_rate      = 0.1,      # 基础事件率（控制时间尺度）
    censor_rate    = 0.4       # 目标删失比例（近似）
) {
  ## 一些检查
  n_structured <- n_groups * genes_per_group
  if (n_structured > n_total_genes)
    stop("n_groups * genes_per_group 不能大于 n_total_genes")

  if (n_signal_groups > n_groups)
    stop("n_signal_groups 不能大于 n_groups")

  ## ---- 1. 构造相关矩阵 R（block-diagonal compound symmetry）----
  # 单个 group 的 compound symmetry 相关矩阵
  make_cs <- function(p, rho) {
    R <- matrix(rho, nrow = p, ncol = p)
    diag(R) <- 1
    R
  }

  R_blocks <- lapply(seq_len(n_groups), function(g) {
    make_cs(genes_per_group, rho_within)
  })
  # block-diagonal 拼接
  R_struct <- as.matrix(Matrix::bdiag(R_blocks))

  # 如有需要，可以在这里检查 R_struct 是否正定：
  # eigen(R_struct)$values

  ## ---- 2. PCA 分解得到 factor loading 矩阵 F ----
  # 这里用 eigen 分解：R = V Λ V^T，
  # 取 F = V Λ^{1/2}，则 F F^T = R
  eig <- eigen(R_struct, symmetric = TRUE)
  V   <- eig$vectors
  L   <- eig$values

  # 数值稳定处理：把极小的负特征值截成 0
  L[L < 0] <- 0
  F_mat <- V %*% diag(sqrt(L), nrow = length(L))

  ## ---- 3. 生成独立正态 X，再乘 F 得到有相关结构的 Z ----
  # X: k × n_samples 独立 N(0,1)
  X_latent <- matrix(rnorm(n_structured * n_samples),
                     nrow = n_structured, ncol = n_samples)

  # Z = F X（注意维度：F(k×k) * X(k×n) = Z(k×n)）
  Z_struct <- F_mat %*% X_latent        # k × n
  Z_struct <- t(Z_struct)               # 变成 n × k（样本在行）

  ## ---- 4. 拼接噪声基因 ----
  X <- matrix(NA, nrow = n_samples, ncol = n_total_genes)
  colnames(X) <- paste0("Gene_", seq_len(n_total_genes))

  # 前 n_structured 列是有结构的 4 groups
  X[, 1:n_structured] <- Z_struct

  # 剩余噪声基因（独立 N(0,1)）
  if (n_structured < n_total_genes) {
    n_noise <- n_total_genes - n_structured
    X[, (n_structured + 1):n_total_genes] <-
      matrix(rnorm(n_samples * n_noise),
             nrow = n_samples, ncol = n_noise)
  }

  ## ---- 5. 构造生存时间：只用前 n_signal_groups 个 group ----
  # 信号基因索引：前 n_signal_groups * genes_per_group 列
  n_signal_genes <- n_signal_groups * genes_per_group
  X_signal <- X[, 1:n_signal_genes, drop = FALSE]

  # 为了稳定，把信号基因先做标准化
  X_signal_sc <- scale(X_signal)

  # 简单起见，用这些信号基因的平均值作为风险评分
  risk_score <- rowMeans(X_signal_sc)

  # log-hazard：lp = β0 + β_signal * risk
  lp            <- beta_signal * risk_score
  lambda_event  <- base_rate * exp(lp)
  T_true        <- rexp(n_samples, rate = lambda_event)

  ## ---- 6. 指数删失，控制删失率约为 censor_rate ----
  # 在简单情形下，用 λ_C ≈ base_rate * c/(1-c) 控制边际删失
  lambda_c <- base_rate * censor_rate / (1 - censor_rate)
  C_time   <- rexp(n_samples, rate = lambda_c)

  Time   <- pmin(T_true, C_time)
  Status <- as.integer(T_true <= C_time)

  ## ---- 7. 返回 ----
  list(
    X    = X,
    surv = data.frame(Time = Time, Status = Status),
    params = list(
      n_samples       = n_samples,
      n_total_genes   = n_total_genes,
      n_groups        = n_groups,
      genes_per_group = genes_per_group,
      rho_within      = rho_within,
      n_signal_groups = n_signal_groups,
      beta_signal     = beta_signal,
      base_rate       = base_rate,
      censor_rate     = censor_rate
    )
  )
}

R <- 100
set.seed(2025)
sim_list_factorial <- replicate(
  R,
  simulate_factorial(),
  simplify = FALSE
)

set.seed(2027)
sim_list_block_gaussian <- replicate(
  R,
  simulate_block_gaussian(
    n_samples = 100,
    n_total_genes = 1000,
    n_blocks = 6,
    genes_per_block = 25,
    n_signal_blocks = 2,
    rho = 0.9,
    signal_strength = 1.0,
    signal_signs = c(1, -1),
    censoring = 0.4
  ),
  simplify = FALSE
)
sim_list <- sim_list_block_gaussian

R <- 100
set.seed(2028)
sim_list_variance_signal <- replicate(
  R,
  simulate_variance_signal(
    n_samples = 80,
    n_total_genes = 3000,
    signal_genes = 8,
    signal_scale = 8.0,
    signal_rho = 0.95,
    decoy_blocks = 12,
    decoy_genes = 50,
    signal_strength = 0.7,
    censoring = 0.45
  ),
  simplify = FALSE
)
sim_list_eigengene <- sim_list_variance_signal
sim_list <- sim_list_eigengene






# parallel
library(doParallel)
library(foreach)

method_files <- c(
  file.path("R", "pcapls_CR.R"),
  file.path("R", "DRPLS.R"),
  file.path("R", "DRPCAPLS.R"),
  file.path("R", "partial_cox.R"),
  file.path("R", "penalized_cox_baselines.R")
)
for (method_file in method_files) {
  if (file.exists(method_file)) {
    source(method_file)
  }
}

.quiet_eval <- function(expr, quiet = TRUE) {
  if (!isTRUE(quiet)) {
    return(eval.parent(substitute(expr)))
  }
  value <- NULL
  invisible(utils::capture.output(
    value <- suppressWarnings(suppressMessages(eval.parent(substitute(expr)))),
    type = "output"
  ))
  value
}

eval_iAUC_one_dataset <- function(
    X, time, status,
    k = 5,
    auc_time_grid = seq(1, 45, length.out = 30),
    ncomp_candidates_DRPLS = 1:10,
    ncomp_candidates_pcox  = 2:10,
    alpha_candidates       = seq(0, 1, by = 0.1),
    lambda_candidates_ridgecox = 10 ^ seq(-4, 1, by = 1),
    alpha_candidates_enetcox = c(0.25, 0.5, 0.75),
    lambda_candidates_enetcox = 10 ^ seq(-4, 1, by = 1),
    use_preselection       = TRUE,
    folds = NULL
) {
  stopifnot(nrow(X) == length(time), length(time) == length(status))
  n <- nrow(X)
  
  # 每个数据集使用同一份 folds 比较所有方法；顺序和并行入口会显式传入。
  if (is.null(folds)) {
    folds <- sample(rep(1:k, length.out = n))
  } else {
    if (length(folds) != n) {
      stop("length(folds) (", length(folds), ") != nrow(X) (", n, ")")
    }
    if (!all(folds %in% seq_len(k))) {
      stop("folds 中的取值应在 1:k 之间")
    }
  }
  
  get_iAUC <- function(res) {
    if (is.null(res) || is.null(res$best_iAUC)) return(NA_real_)
    val <- res$best_iAUC
    if (!is.finite(val)) return(NA_real_)
    val
  }
  get_fixed_alpha_iAUC <- function(res, alpha_value) {
    if (is.null(res) || is.null(res$results)) return(NA_real_)
    keep <- is.finite(res$results$iAUC) &
      abs(res$results$alpha - alpha_value) < sqrt(.Machine$double.eps)
    if (!any(keep)) return(NA_real_)
    max(res$results$iAUC[keep], na.rm = TRUE)
  }
  get_selected_alpha <- function(res) {
    if (is.null(res) || is.null(res$best_alpha)) return(NA_real_)
    val <- res$best_alpha
    if (!is.finite(val)) return(NA_real_)
    val
  }
  get_selected_ncomp <- function(res) {
    if (is.null(res) || is.null(res$best_ncomp)) return(NA_real_)
    val <- res$best_ncomp
    if (!is.finite(val)) return(NA_real_)
    val
  }
  
  res_drpls <- tryCatch(
    val_DRPLS_cv(
      X, time, status,
      k = k,
      ncomp_candidates = ncomp_candidates_DRPLS,
      auc_time_grid    = auc_time_grid,
      use_preselection = use_preselection,
      folds            = folds,
      seed             = NULL
    ),
    error = function(e) NULL
  )
  
  res_drpcapls <- tryCatch(
    val_DRPCAPLS_cv(
      X, time, status,
      k = k,
      ncomp_candidates = ncomp_candidates_DRPLS,
      alpha_candidates = sort(unique(c(0, 1, alpha_candidates))),
      auc_time_grid    = auc_time_grid,
      use_preselection = use_preselection,
      folds            = folds,
      seed             = NULL
    ),
    error = function(e) NULL
  )
  
  res_pcox <- tryCatch(
    val_partial_cox_cv(
      X, time, status,
      k = k,
      ncomp_candidates = ncomp_candidates_pcox,
      auc_time_grid    = auc_time_grid,
      use_preselection = use_preselection,
      folds            = folds,
      seed             = NULL
    ),
    error = function(e) NULL
  )

  res_ridgecox <- tryCatch(
    val_ridge_cox_cv(
      X, time, status,
      k = k,
      lambda_candidates = lambda_candidates_ridgecox,
      auc_time_grid    = auc_time_grid,
      use_preselection = use_preselection,
      folds            = folds,
      seed             = NULL
    ),
    error = function(e) {
      warning("RidgeCox failed: ", e$message)
      NULL
    }
  )

  res_enetcox <- tryCatch(
    val_elastic_net_cox_cv(
      X, time, status,
      k = k,
      alpha_candidates = alpha_candidates_enetcox,
      lambda_candidates = lambda_candidates_enetcox,
      auc_time_grid    = auc_time_grid,
      use_preselection = use_preselection,
      folds            = folds,
      seed             = NULL
    ),
    error = function(e) {
      warning("ElasticNetCox failed: ", e$message)
      NULL
    }
  )
  
  list(
    iAUC = c(
      DRPLS          = get_iAUC(res_drpls),
      CPPC_alpha0    = get_fixed_alpha_iAUC(res_drpcapls, 0),
      CPPC_alpha1    = get_fixed_alpha_iAUC(res_drpcapls, 1),
      DRPCA_PLS      = get_iAUC(res_drpcapls),
      partialCox     = get_iAUC(res_pcox),
      RidgeCox       = get_iAUC(res_ridgecox),
      ElasticNetCox  = get_iAUC(res_enetcox)
    ),
    selected_alpha = get_selected_alpha(res_drpcapls),
    selected_ncomp = get_selected_ncomp(res_drpcapls)
  )
}

mc_iAUC_seq <- function(
    sim_list,
    k = 5,
    auc_time_grid = seq(1, 45, length.out = 30),
    ncomp_candidates_DRPLS = 1:10,
    ncomp_candidates_pcox  = 2:10,
    alpha_candidates       = seq(0, 1, by = 0.1),
    lambda_candidates_ridgecox = 10 ^ seq(-4, 1, by = 1),
    alpha_candidates_enetcox = c(0.25, 0.5, 0.75),
    lambda_candidates_enetcox = 10 ^ seq(-4, 1, by = 1),
    use_preselection       = TRUE,
    seed = 123,
    quiet = TRUE
) {
  R <- length(sim_list)
  method_names <- c(
    "DRPLS", "CPPC_alpha0", "CPPC_alpha1", "DRPCA_PLS",
    "partialCox", "RidgeCox", "ElasticNetCox"
  )
  auc_mat <- matrix(NA_real_, nrow = R, ncol = length(method_names))
  selected_alpha <- rep(NA_real_, R)
  selected_ncomp <- rep(NA_real_, R)
  colnames(auc_mat) <- c(
    method_names
  )
  
  for (r in seq_len(R)) {
    dat <- sim_list[[r]]
    X <- dat$X
    time <- dat$surv$Time
    status <- dat$surv$Status
    if (!is.null(seed)) {
      set.seed(seed + r - 1L)
    }
    folds <- sample(rep(seq_len(k), length.out = nrow(X)))
    
    eval_res <- .quiet_eval(eval_iAUC_one_dataset(
      X, time, status,
      k = k,
      auc_time_grid = auc_time_grid,
      ncomp_candidates_DRPLS = ncomp_candidates_DRPLS,
      ncomp_candidates_pcox  = ncomp_candidates_pcox,
      alpha_candidates       = alpha_candidates,
      lambda_candidates_ridgecox = lambda_candidates_ridgecox,
      alpha_candidates_enetcox = alpha_candidates_enetcox,
      lambda_candidates_enetcox = lambda_candidates_enetcox,
      use_preselection       = use_preselection,
      folds                  = folds
    ), quiet = quiet)
    auc_mat[r, ] <- eval_res$iAUC[method_names]
    selected_alpha[r] <- eval_res$selected_alpha
    selected_ncomp[r] <- eval_res$selected_ncomp
  }
  
  mean_iAUC <- colMeans(auc_mat, na.rm = TRUE)
  sd_iAUC   <- apply(auc_mat, 2, sd, na.rm = TRUE)
  n_eff     <- colSums(is.finite(auc_mat))
  mcse_iAUC <- sd_iAUC / sqrt(n_eff)
  
  ci_mean <- cbind(
    CI_lower = mean_iAUC - 1.96 * mcse_iAUC,
    CI_upper = mean_iAUC + 1.96 * mcse_iAUC
  )
  
  list(draws = auc_mat, selected_alpha = selected_alpha,
       selected_ncomp = selected_ncomp, mean = mean_iAUC, sd = sd_iAUC,
       mcse = mcse_iAUC, ci_mean = ci_mean, n_eff = n_eff, R = R)
}

.format_mean_sd <- function(mean, sd, digits = 3) {
  ifelse(
    is.finite(mean) & is.finite(sd),
    paste0(formatC(mean, format = "f", digits = digits),
           " ± ",
           formatC(sd, format = "f", digits = digits)),
    NA_character_
  )
}

.format_p_value <- function(p, digits = 3) {
  ifelse(
    is.finite(p),
    ifelse(p < 10 ^ (-digits),
           paste0("<", formatC(10 ^ (-digits), format = "f", digits = digits)),
           formatC(p, format = "f", digits = digits)),
    NA_character_
  )
}

summarize_mc_iAUC <- function(mc_res,
                              proposed = "DRPCA_PLS",
                              baselines = NULL,
                              digits = 3) {
  draws <- as.data.frame(mc_res$draws)
  methods <- colnames(draws)
  if (is.null(baselines)) {
    baselines <- setdiff(methods, proposed)
  }
  
  mean_vals <- colMeans(mc_res$draws, na.rm = TRUE)
  sd_vals <- apply(mc_res$draws, 2, stats::sd, na.rm = TRUE)
  n_eff <- colSums(is.finite(mc_res$draws))
  mcse_vals <- sd_vals / sqrt(n_eff)
  ci_lower <- mean_vals - 1.96 * mcse_vals
  ci_upper <- mean_vals + 1.96 * mcse_vals
  
  summary_table <- data.frame(
    Method = methods,
    Mean = as.numeric(mean_vals[methods]),
    SD = as.numeric(sd_vals[methods]),
    Mean_SD = .format_mean_sd(mean_vals[methods], sd_vals[methods], digits),
    MCSE = as.numeric(mcse_vals[methods]),
    CI_lower = as.numeric(ci_lower[methods]),
    CI_upper = as.numeric(ci_upper[methods]),
    N_eff = as.integer(n_eff[methods]),
    row.names = NULL
  )
  
  paired_table <- do.call(rbind, lapply(baselines, function(base) {
    if (!(proposed %in% methods) || !(base %in% methods)) {
      return(NULL)
    }
    keep <- is.finite(draws[[proposed]]) & is.finite(draws[[base]])
    diff_vals <- draws[[proposed]][keep] - draws[[base]][keep]
    n_pair <- length(diff_vals)
    diff_mean <- if (n_pair > 0L) mean(diff_vals) else NA_real_
    diff_sd <- if (n_pair > 1L) stats::sd(diff_vals) else NA_real_
    diff_mcse <- if (n_pair > 1L) diff_sd / sqrt(n_pair) else NA_real_
    p_val <- if (n_pair > 1L && any(diff_vals != 0)) {
      tryCatch(
        stats::wilcox.test(draws[[proposed]][keep],
                           draws[[base]][keep],
                           paired = TRUE,
                           exact = FALSE)$p.value,
        error = function(e) NA_real_
      )
    } else {
      NA_real_
    }
    data.frame(
      Proposed = proposed,
      Baseline = base,
      N_pair = n_pair,
      Mean_diff = diff_mean,
      SD_diff = diff_sd,
      MCSE_diff = diff_mcse,
      CI_lower = diff_mean - 1.96 * diff_mcse,
      CI_upper = diff_mean + 1.96 * diff_mcse,
      Wilcoxon_p = p_val,
      Wilcoxon_p_fmt = .format_p_value(p_val, digits),
      row.names = NULL
    )
  }))
  
  if (is.null(paired_table)) {
    paired_table <- data.frame(
      Proposed = character(),
      Baseline = character(),
      N_pair = integer(),
      Mean_diff = numeric(),
      SD_diff = numeric(),
      MCSE_diff = numeric(),
      CI_lower = numeric(),
      CI_upper = numeric(),
      Wilcoxon_p = numeric(),
      Wilcoxon_p_fmt = character()
    )
  }
  
  list(
    summary = summary_table,
    paired_difference = paired_table
  )
}

plot_paired_difference <- function(mc_res,
                                   proposed = "DRPCA_PLS",
                                   baselines = c(
                                     "DRPLS", "partialCox",
                                     "RidgeCox", "ElasticNetCox"
                                   ),
                                   method_labels = c(
                                     DRPLS = "DRPLS",
                                     CPPC_alpha0 = "CPPC-Cox alpha=0",
                                     CPPC_alpha1 = "CPPC-Cox alpha=1",
                                     DRPCA_PLS = "ContiCox",
                                     partialCox = "partialCox",
                                     RidgeCox = "Ridge Cox",
                                     ElasticNetCox = "Elastic-net Cox"
                                   ),
                                   diff_limits = NULL,
                                   diff_breaks = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plot_paired_difference().")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package `tidyr` is required for plot_paired_difference().")
  }
  
  draws <- as.data.frame(mc_res$draws)
  if (!(proposed %in% colnames(draws))) {
    stop("proposed method not found in mc_res$draws: ", proposed)
  }
  baselines <- baselines[baselines %in% colnames(draws)]
  if (length(baselines) < 1L) {
    stop("No requested baselines were found in mc_res$draws.")
  }
  
  paired_table <- summarize_mc_iAUC(
    mc_res,
    proposed = proposed,
    baselines = baselines
  )$paired_difference
  paired_table$Baseline_label <- ifelse(
    paired_table$Baseline %in% names(method_labels),
    unname(method_labels[paired_table$Baseline]),
    paired_table$Baseline
  )
  baseline_order <- ifelse(
    baselines %in% names(method_labels),
    unname(method_labels[baselines]),
    baselines
  )
  paired_table$Baseline_label <- factor(
    paired_table$Baseline_label,
    levels = rev(baseline_order)
  )
  
  diff_wide <- as.data.frame(lapply(baselines, function(base) {
    draws[[proposed]] - draws[[base]]
  }))
  colnames(diff_wide) <- baselines
  diff_wide$rep <- seq_len(nrow(diff_wide))
  diff_long <- tidyr::pivot_longer(
    diff_wide,
    cols = -rep,
    names_to = "Baseline",
    values_to = "Difference"
  )
  diff_long <- diff_long[is.finite(diff_long$Difference), , drop = FALSE]
  diff_long$Baseline_label <- ifelse(
    diff_long$Baseline %in% names(method_labels),
    unname(method_labels[diff_long$Baseline]),
    diff_long$Baseline
  )
  diff_long$Baseline_label <- factor(
    diff_long$Baseline_label,
    levels = baseline_order
  )
  
  x_scale <- list()
  y_scale <- list()
  if (!is.null(diff_breaks)) {
    x_scale <- c(x_scale, list(ggplot2::scale_x_continuous(breaks = diff_breaks)))
    y_scale <- c(y_scale, list(ggplot2::scale_y_continuous(breaks = diff_breaks)))
  }
  if (!is.null(diff_limits)) {
    x_scale <- c(x_scale, list(ggplot2::coord_cartesian(xlim = diff_limits)))
    y_scale <- c(y_scale, list(ggplot2::coord_cartesian(ylim = diff_limits)))
  }
  
  forest <- ggplot2::ggplot(paired_table, ggplot2::aes(y = Baseline_label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_segment(
      ggplot2::aes(x = CI_lower, xend = CI_upper,
                   yend = Baseline_label),
      linewidth = 0.55
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = Mean_diff),
      shape = 21, size = 2.4, fill = "red", color = "black"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    x_scale +
    ggplot2::labs(
      x = "Paired iAUC difference: ContiCox - baseline",
      y = "Baseline"
    )
  
  violin <- ggplot2::ggplot(
    diff_long,
    ggplot2::aes(x = Baseline_label, y = Difference, fill = Baseline_label)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.35, color = NA, width = 0.55) +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA,
                          fill = "white", color = "black", linewidth = 0.35) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 21,
                          size = 2.2, fill = "red", color = "black") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 15, hjust = 1)) +
    y_scale +
    ggplot2::labs(
      x = "Baseline",
      y = "Paired iAUC difference: ContiCox - baseline"
    )
  
  list(
    summary = paired_table,
    differences = diff_long,
    forest = forest,
    violin = violin
  )
}

summarize_selected_alpha <- function(mc_res, digits = 3) {
  if (is.null(mc_res$selected_alpha)) {
    stop("mc_res does not contain selected_alpha.")
  }
  
  alpha <- mc_res$selected_alpha
  alpha <- alpha[is.finite(alpha)]
  if (length(alpha) < 1L) {
    empty_summary <- data.frame(
      N = 0L,
      Mean = NA_real_,
      SD = NA_real_,
      Mean_SD = NA_character_,
      Median = NA_real_,
      Q25 = NA_real_,
      Q75 = NA_real_,
      Prop_low = NA_real_,
      Prop_mid = NA_real_,
      Prop_high = NA_real_
    )
    return(list(summary = empty_summary, frequency = data.frame()))
  }
  
  summary_table <- data.frame(
    N = length(alpha),
    Mean = mean(alpha),
    SD = stats::sd(alpha),
    Mean_SD = .format_mean_sd(mean(alpha), stats::sd(alpha), digits),
    Median = stats::median(alpha),
    Q25 = as.numeric(stats::quantile(alpha, 0.25, names = FALSE)),
    Q75 = as.numeric(stats::quantile(alpha, 0.75, names = FALSE)),
    Prop_low = mean(alpha <= 0.3),
    Prop_mid = mean(alpha > 0.3 & alpha < 0.7),
    Prop_high = mean(alpha >= 0.7),
    row.names = NULL
  )
  
  alpha_levels <- sort(unique(c(seq(0, 1, by = 0.1), alpha)))
  frequency <- as.data.frame(table(
    factor(alpha, levels = alpha_levels),
    useNA = "no"
  ))
  colnames(frequency) <- c("selected_alpha", "Freq")
  frequency$selected_alpha <- as.numeric(as.character(frequency$selected_alpha))
  frequency$Prop <- frequency$Freq / sum(frequency$Freq)
  frequency <- frequency[frequency$Freq > 0, , drop = FALSE]
  rownames(frequency) <- NULL
  
  list(summary = summary_table, frequency = frequency)
}

scientific_iAUC_palette <- function(method_levels) {
  base_palette <- c(
    "DRPLS" = "#0072B2",
    "CPPC-Cox alpha=0" = "#56B4E9",
    "CPPC-Cox alpha=1" = "#E69F00",
    "CPPC-Cox CV" = "#009E73",
    "ContiCox" = "#009E73",
    "partialCox" = "#D55E00",
    "Ridge Cox" = "#CC79A7",
    "Elastic-net Cox" = "#999999"
  )
  missing_levels <- setdiff(method_levels, names(base_palette))
  if (length(missing_levels) > 0L) {
    base_palette <- c(
      base_palette,
      stats::setNames(
        grDevices::hcl.colors(length(missing_levels), palette = "Dark 3"),
        missing_levels
      )
    )
  }
  base_palette[method_levels]
}

plot_mc_iAUC <- function(mc_res,
                         method_labels = c(
                           DRPLS = "DRPLS",
                           CPPC_alpha0 = "CPPC-Cox alpha=0",
                           CPPC_alpha1 = "CPPC-Cox alpha=1",
                           DRPCA_PLS = "CPPC-Cox CV",
                           partialCox = "partialCox",
                           RidgeCox = "Ridge Cox",
                           ElasticNetCox = "Elastic-net Cox"
                         ),
                         method_order = names(method_labels),
                         y_limits = c(0.5, 1.0),
                         y_breaks = seq(0.5, 1.0, by = 0.1),
                         method_palette = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plot_mc_iAUC().")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package `tidyr` is required for plot_mc_iAUC().")
  }
  
  draws <- as.data.frame(mc_res$draws)
  draws$rep <- seq_len(nrow(draws))
  df_long <- tidyr::pivot_longer(
    draws,
    cols = -rep,
    names_to = "Method",
    values_to = "iAUC"
  )
  df_long <- df_long[is.finite(df_long$iAUC), , drop = FALSE]
  df_long$Method <- ifelse(
    df_long$Method %in% names(method_labels),
    unname(method_labels[df_long$Method]),
    df_long$Method
  )
  order_labels <- unname(method_labels[method_order[method_order %in% names(method_labels)]])
  df_long$Method <- factor(df_long$Method, levels = order_labels)
  y_scale <- list()
  if (!is.null(y_breaks)) {
    y_scale <- c(y_scale, list(ggplot2::scale_y_continuous(breaks = y_breaks)))
  }
  if (!is.null(y_limits)) {
    y_scale <- c(y_scale, list(ggplot2::coord_cartesian(ylim = y_limits)))
  }
  if (is.null(method_palette)) {
    method_palette <- scientific_iAUC_palette(levels(df_long$Method))
  } else {
    method_palette <- method_palette[levels(df_long$Method)]
  }
  fill_scale <- ggplot2::scale_fill_manual(values = method_palette, drop = FALSE)
  
  boxplot <- ggplot2::ggplot(df_long, ggplot2::aes(x = Method, y = iAUC, fill = Method)) +
    ggplot2::geom_boxplot(outlier.size = 0.8, width = 0.55) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 21,
                          size = 2.2, fill = "white", color = "black") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 15, hjust = 1)) +
    fill_scale +
    y_scale +
    ggplot2::labs(x = "Method", y = "iAUC")
  
  violin <- ggplot2::ggplot(df_long, ggplot2::aes(x = Method, y = iAUC, fill = Method)) +
    ggplot2::geom_violin(trim = FALSE, scale = "width", alpha = 0.35,
                         color = NA, width = 0.55) +
    ggplot2::geom_boxplot(width = 0.10, outlier.shape = NA,
                          fill = "white", color = "black", linewidth = 0.35) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 21,
                          size = 2.2, fill = "white", color = "black") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 15, hjust = 1)) +
    fill_scale +
    y_scale +
    ggplot2::labs(x = "Method", y = "iAUC")
  
  list(data = df_long, boxplot = boxplot, violin = violin)
}

plot_main_iAUC <- function(mc_res,
                           method_labels = c(
                             DRPLS = "PLSDR",
                             DRPCA_PLS = "ContiCox",
                             partialCox = "pCox",
                             RidgeCox = "Ridge-Cox",
                             ElasticNetCox = "EN-Cox"
                           ),
                           method_order = names(method_labels),
                           y_limits = c(0.5, 1.0),
                           y_breaks = seq(0.5, 1.0, by = 0.1),
                           method_palette = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plot_main_iAUC().")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package `tidyr` is required for plot_main_iAUC().")
  }
  
  draws <- as.data.frame(mc_res$draws)
  keep_methods <- method_order[method_order %in% colnames(draws)]
  draws <- draws[, keep_methods, drop = FALSE]
  draws$rep <- seq_len(nrow(draws))
  df_long <- tidyr::pivot_longer(
    draws,
    cols = -rep,
    names_to = "Method",
    values_to = "iAUC"
  )
  df_long <- df_long[is.finite(df_long$iAUC), , drop = FALSE]
  df_long$Method <- unname(method_labels[df_long$Method])
  order_labels <- unname(method_labels[keep_methods])
  df_long$Method <- factor(df_long$Method, levels = order_labels)
  y_scale <- list()
  if (!is.null(y_breaks)) {
    y_scale <- c(y_scale, list(ggplot2::scale_y_continuous(breaks = y_breaks)))
  }
  if (!is.null(y_limits)) {
    y_scale <- c(y_scale, list(ggplot2::coord_cartesian(ylim = y_limits)))
  }
  if (is.null(method_palette)) {
    method_palette <- c(
      "PLSDR" = "#4C78A8",
      "ContiCox" = "#E45756",
      "pCox" = "#72B7B2",
      "Ridge-Cox" = "#B279A2",
      "EN-Cox" = "#79706E"
    )[levels(df_long$Method)]
  } else {
    method_palette <- method_palette[levels(df_long$Method)]
  }
  fill_scale <- ggplot2::scale_fill_manual(values = method_palette, drop = FALSE)
  
  boxplot <- ggplot2::ggplot(df_long, ggplot2::aes(x = Method, y = iAUC, fill = Method)) +
    ggplot2::geom_boxplot(outlier.size = 0.8, width = 0.55) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 21,
                          size = 2.2, fill = "white", color = "black") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 15, hjust = 1)) +
    fill_scale +
    y_scale +
    ggplot2::labs(x = "Method", y = "iAUC")
  
  violin <- ggplot2::ggplot(df_long, ggplot2::aes(x = Method, y = iAUC, fill = Method)) +
    ggplot2::geom_violin(trim = FALSE, scale = "width", alpha = 0.35,
                         color = NA, width = 0.55) +
    ggplot2::geom_boxplot(width = 0.10, outlier.shape = NA,
                          fill = "white", color = "black", linewidth = 0.35) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 21,
                          size = 2.2, fill = "white", color = "black") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 15, hjust = 1)) +
    fill_scale +
    y_scale +
    ggplot2::labs(x = "Method", y = "iAUC")
  
  list(data = df_long, boxplot = boxplot, violin = violin)
}

plot_selected_alpha <- function(mc_res) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plot_selected_alpha().")
  }
  if (is.null(mc_res$selected_alpha)) {
    stop("mc_res does not contain selected_alpha.")
  }
  
  df_alpha <- data.frame(
    rep = seq_along(mc_res$selected_alpha),
    selected_alpha = mc_res$selected_alpha
  )
  df_alpha <- df_alpha[is.finite(df_alpha$selected_alpha), , drop = FALSE]
  
  ggplot2::ggplot(df_alpha, ggplot2::aes(x = selected_alpha)) +
    ggplot2::geom_bar(width = 0.08, fill = "#4C78A8", color = "white") +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, by = 0.1),
                                limits = c(-0.05, 1.05)) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(x = "Selected alpha", y = "Frequency")
}



mc_res_seq <- mc_iAUC_seq(sim_list[1])

mc_res_seq$draws

mc_res_seq <- mc_iAUC_seq(sim_list[1])
mc_res_seq$draws
mc_res_seq$mean
mc_res_seq_stats <- summarize_mc_iAUC(mc_res_seq)
mc_res_seq_stats$summary
mc_res_seq_stats$paired_difference
mc_res_seq_alpha_stats <- summarize_selected_alpha(mc_res_seq)
mc_res_seq_alpha_stats$summary
mc_res_seq_alpha_stats$frequency

mc_iAUC_par <- function(
    sim_list,
    k = 5,
    auc_time_grid = seq(1, 45, length.out = 30),
    ncomp_candidates_DRPLS = 1:10,
    ncomp_candidates_pcox  = 2:10,
    alpha_candidates       = seq(0, 1, by = 0.1),
    lambda_candidates_ridgecox = 10 ^ seq(-4, 1, by = 1),
    alpha_candidates_enetcox = c(0.25, 0.5, 0.75),
    lambda_candidates_enetcox = 10 ^ seq(-4, 1, by = 1),
    use_preselection       = TRUE,
    n_cores = max(1L, parallel::detectCores() - 2L),
    seed = 123,
    project_dir = NULL,
    quiet = TRUE
) {
  `%dopar%` <- foreach::`%dopar%`
  R <- length(sim_list)
  if (R < 1L) {
    stop("sim_list must contain at least one simulated dataset.")
  }
  
  n_cores <- max(1L, min(as.integer(n_cores), R))
  required_files <- c(
    file.path("R", "pcapls_CR.R"),
    file.path("R", "DRPLS.R"),
    file.path("R", "DRPCAPLS.R"),
    file.path("R", "partial_cox.R"),
    file.path("R", "penalized_cox_baselines.R")
  )
  
  if (is.null(project_dir)) {
    candidate_dirs <- unique(c(
      getwd(),
      "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls"
    ))
    has_all_files <- vapply(
      candidate_dirs,
      function(path) all(file.exists(file.path(path, required_files))),
      logical(1)
    )
    if (!any(has_all_files)) {
      stop(
        "Cannot locate method source files. Please call mc_iAUC_par(..., ",
        "project_dir = \"/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls\")."
      )
    }
    project_dir <- candidate_dirs[which(has_all_files)[1L]]
  }
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  missing_files <- required_files[!file.exists(file.path(project_dir, required_files))]
  if (length(missing_files) > 0L) {
    stop(
      "project_dir does not contain required files: ",
      paste(missing_files, collapse = ", ")
    )
  }
  
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)
  
  parallel::clusterCall(cl, function(project_dir, required_files, quiet) {
    load_worker <- function() {
      setwd(project_dir)
      missing_files <- required_files[!file.exists(required_files)]
      if (length(missing_files) > 0L) {
        stop("Worker cannot find files: ", paste(missing_files, collapse = ", "))
      }
      source(file.path("R", "pcapls_CR.R"))
      source(file.path("R", "DRPLS.R"))
      source(file.path("R", "DRPCAPLS.R"))
      source(file.path("R", "partial_cox.R"))
      source(file.path("R", "penalized_cox_baselines.R"))
      NULL
    }
    if (isTRUE(quiet)) {
      invisible(utils::capture.output(
        suppressWarnings(suppressMessages(load_worker())),
        type = "output"
      ))
    } else {
      load_worker()
    }
    NULL
  }, project_dir, required_files, quiet)
  
  auc_mat <- foreach::foreach(
    r = seq_len(R),
    .combine = rbind,
    .export = c("eval_iAUC_one_dataset", ".quiet_eval"),
    .packages = c(
      "survival", "pls", "survivalROC", "prodlim",
      "survcomp", "survAUC", "RSpectra", "glmnet"
    )
  ) %dopar% {
    if (!is.null(seed)) {
      set.seed(seed + r - 1L)
    }
    
    dat <- sim_list[[r]]
    X <- dat$X
    time <- dat$surv$Time
    status <- dat$surv$Status
    folds <- sample(rep(seq_len(k), length.out = nrow(X)))
    
    eval_res <- .quiet_eval(eval_iAUC_one_dataset(
      X, time, status,
      k = k,
      auc_time_grid = auc_time_grid,
      ncomp_candidates_DRPLS = ncomp_candidates_DRPLS,
      ncomp_candidates_pcox  = ncomp_candidates_pcox,
      alpha_candidates       = alpha_candidates,
      lambda_candidates_ridgecox = lambda_candidates_ridgecox,
      alpha_candidates_enetcox = alpha_candidates_enetcox,
      lambda_candidates_enetcox = lambda_candidates_enetcox,
      use_preselection       = use_preselection,
      folds                  = folds
    ), quiet = quiet)
    c(eval_res$iAUC,
      selected_alpha = eval_res$selected_alpha,
      selected_ncomp = eval_res$selected_ncomp)
  }
  
  method_names <- c(
    "DRPLS", "CPPC_alpha0", "CPPC_alpha1", "DRPCA_PLS",
    "partialCox", "RidgeCox", "ElasticNetCox"
  )
  worker_names <- c(method_names, "selected_alpha", "selected_ncomp")
  auc_mat <- if (is.null(dim(auc_mat))) {
    matrix(auc_mat, nrow = 1L, byrow = TRUE)
  } else {
    as.matrix(auc_mat)
  }
  if (is.null(colnames(auc_mat))) {
    colnames(auc_mat) <- worker_names
  }
  selected_alpha <- as.numeric(auc_mat[, "selected_alpha"])
  selected_ncomp <- as.numeric(auc_mat[, "selected_ncomp"])
  auc_mat <- auc_mat[, method_names, drop = FALSE]
  colnames(auc_mat) <- c(
    method_names
  )
  
  mean_iAUC <- colMeans(auc_mat, na.rm = TRUE)
  sd_iAUC   <- apply(auc_mat, 2, sd, na.rm = TRUE)
  n_eff     <- colSums(is.finite(auc_mat))
  mcse_iAUC <- sd_iAUC / sqrt(n_eff)
  
  ci_mean <- cbind(
    CI_lower = mean_iAUC - 1.96 * mcse_iAUC,
    CI_upper = mean_iAUC + 1.96 * mcse_iAUC
  )
  
  list(draws = auc_mat, selected_alpha = selected_alpha,
       selected_ncomp = selected_ncomp, mean = mean_iAUC, sd = sd_iAUC,
       mcse = mcse_iAUC, ci_mean = ci_mean, n_eff = n_eff, R = R)
}

mc_res_par <- mc_iAUC_par(
  sim_list[1:100],
  n_cores = 10,
  project_dir = "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls"
)

load("/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls/simulation/mc_res_par_0511_eigengene_opposite_modules.RData")

mc_res_par$draws
mc_res_par$mean
mc_res_par_stats <- summarize_mc_iAUC(mc_res_par)
mc_res_par_stats$summary
mc_res_par_stats$paired_difference
mc_res_par_plots <- plot_mc_iAUC(mc_res_par)
mc_res_par_plots$boxplot
mc_res_par_plots$violin
mc_res_main_plots <- plot_main_iAUC(
  mc_res_par,
  y_limits = c(0.5, 1.0),
  y_breaks = seq(0.5, 1.0, by = 0.1)
)
mc_res_main_plots$boxplot
mc_res_main_plots$violin
paired_difference_plots <- plot_paired_difference(mc_res_par)
paired_difference_plots$forest
paired_difference_plots$violin
selected_alpha_stats <- summarize_selected_alpha(mc_res_par)
selected_alpha_stats$summary
selected_alpha_stats$frequency
selected_alpha_plot <- plot_selected_alpha(mc_res_par)
selected_alpha_plot



save(mc_res_par, file = "/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls/simulation/mc_res_par_0511_eigengene_opposite_modules.RData")

load("/Users/ruijuanzhong/CR_PCA_PLS/cox_pcapls/mc_res_seq_highcorr.RData")


# plot
draws <- as.data.frame(mc_res_par$draws[, c(
  "DRPLS", "DRPCA_PLS", "partialCox", "RidgeCox", "ElasticNetCox"
), drop = FALSE])
library(tidyr)
library(dplyr)
library(ggplot2)

df_long <- draws %>%
  mutate(rep = row_number()) %>%
  pivot_longer(-rep, names_to = "Method", values_to = "iAUC")



ggplot(df_long, aes(Method, iAUC)) +
  geom_boxplot(outlier.size = 0.8) +
  theme_minimal() +
  labs(y = "iAUC")

sum_df <- data.frame(
  Method = names(mc_res_par$mean),
  Mean   = as.numeric(mc_res_par$mean),
  L      = mc_res_par$ci_mean[,1],
  U      = mc_res_par$ci_mean[,2]
)

ggplot(sum_df, aes(Method, Mean)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = L, ymax = U), width = 0.15) +
  theme_minimal() +
  labs(y = "Mean iAUC (95% MC CI)")

draws$diff_DRPCA_vs_pCox <- draws$DRPCA_PLS - draws$partialCox

ggplot(draws, aes(x = diff_DRPCA_vs_pCox)) +
  geom_histogram(bins = 20) +
  theme_minimal() +
  labs(x = "DRPCA_PLS - partialCox", y = "Count")




library(ggplot2)
library(dplyr)
library(tidyr)

draws <- as.data.frame(mc_res_par$draws)

df_long <- draws %>%
  mutate(rep = row_number()) %>%
  pivot_longer(-rep, names_to = "Method", values_to = "iAUC") %>%
  mutate(
    Method = recode(Method,
                    "DRPLS" = "DRPLS",
                    "DRPCA_PLS" = "ContiCox",
                    "partialCox" = "partialCox",
                    "RidgeCox" = "Ridge Cox",
                    "ElasticNetCox" = "Elastic-net Cox"
    )
  ) %>%
  filter(is.finite(iAUC))

mean_df <- df_long %>%
  group_by(Method) %>%
  summarise(mean_iAUC = mean(iAUC), .groups = "drop")

df_long$Method <- factor(df_long$Method,
                         levels = c("DRPLS", "ContiCox", "partialCox",
                                    "Ridge Cox", "Elastic-net Cox")
)
mean_df$Method <- factor(mean_df$Method, levels = levels(df_long$Method))

# 先算一个“自动留白”的 y 轴范围
library(ggplot2)

# 如果你还没固定顺序
df_long$Method <- factor(df_long$Method,
                         levels = c("DRPLS", "ContiCox", "partialCox",
                                    "Ridge Cox", "Elastic-net Cox")
)

p_target <- ggplot(df_long, aes(x = Method, y = iAUC, fill = Method)) +
  geom_violin(
    trim = FALSE,        # 关键 1：让顶部/底部有自然尾巴
    scale = "width",
    alpha = 0.35,
    color = NA,
    width = 0.42,
    adjust = 1.0
  ) +
  geom_boxplot(
    width = 0.08,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.35
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2.2,
    shape = 21,
    fill = "red",
    color = "black",
    stroke = 0.4
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +  # 关键 2：上端留白
  coord_cartesian(ylim = c(0.84, 0.96)) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 10, hjust = 0.5)
  ) +
  labs(x = "Method", y = "iAUC")

print(p_target)


