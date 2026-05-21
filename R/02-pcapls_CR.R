
# PCA-PLS with alpha

.pcapls_safe_positive <- function(x, eps = 1e-12) {
  x <- as.numeric(x)[1]
  if (!is.finite(x) || x <= eps) return(1)
  x
}

.pcapls_pca_norm <- function(X, eps = 1e-12) {
  val <- tryCatch(
    RSpectra::svds(X, k = 1, nu = 0, nv = 0)$d[1]^2,
    error = function(e) svd(X, nu = 0, nv = 0)$d[1]^2
  )
  .pcapls_safe_positive(val, eps = eps)
}

.pcapls_pls_norm <- function(X, Y, eps = 1e-12) {
  G <- crossprod(X, Y)
  val <- if (ncol(G) == 1L) {
    sum(G^2)
  } else {
    svd(G, nu = 0, nv = 0)$d[1]^2
  }
  .pcapls_safe_positive(val, eps = eps)
}

pca_pls_with_alpha <- function(X, Y, alpha, ncomp,
                               normalize_terms = FALSE,
                               norm_eps = 1e-12) {      
  n <- nrow(X)
  p <- ncol(X)
  if (is.vector(Y)) Y <- matrix(Y, ncol = 1)
  q <- ncol(Y)  # Y 
  
  # 
  T_scores <- matrix(0, n, ncomp)  #  X-scores
  P_loadings <- matrix(0, p, ncomp)  #  X-loadings
  R_loadings <- matrix(0, q, ncomp)  #  y-loadings ( q  ncomp )
  W_weights <- matrix(0, p, ncomp)  #  w  (p  ncomp )
  print("beginning PCA-PLS with alpha")
  
  for (i in 1:ncomp) {
    # Step 1:  H 
    print(paste("H matrix: ", i))
    if (isTRUE(normalize_terms)) {
      pca_scale <- .pcapls_pca_norm(X, eps = norm_eps)
      pls_scale <- .pcapls_pls_norm(X, Y, eps = norm_eps)
    } else {
      pca_scale <- 1
      pls_scale <- 1
    }
    H <- ((1 - alpha) / pca_scale) * crossprod(X) +
      (alpha / pls_scale) * crossprod(X, Y) %*% crossprod(Y, X)
    print("H matrix calculated.")
    #print(paste("H matrix: ", i))
    
    # Step 2:  w
    eigen_decomp <- eigs_sym(H, k = 1, which = "LM")
    print("Eigen decomposition done.")
    
    w <- eigen_decomp$vectors[, 1]  
    W_weights[, i] <- w  #  w 
    
    # Step 3:  X-scores: t = X %*% w
    t_score <- X %*% w
    
    # Step 4:  X-loadings: p = X' t / (t' t)
    p_loading <- crossprod(X, t_score) / as.numeric(crossprod(t_score))
    
    # Step 5:  y-loadings: r = Y' t / (t' t)
    r_loading <- crossprod(Y, t_score) / as.numeric(crossprod(t_score))  #  q x 1
    
    # Step 6: Deflate X  Y: X = X - t p', Y = Y - t r'
    X <- X - t_score %*% t(p_loading)
    Y <- Y - t_score %*% t(r_loading)  # Y  r_loading  q x 1
    
    # 
    T_scores[, i] <- t_score
    P_loadings[, i] <- p_loading
    R_loadings[, i] <- r_loading
  }
  
  # Step 7:  R
  #  P^T W
  PW <- crossprod(P_loadings, W_weights)
  
  #  (P^T W) 
  #PWinv <- solve(PW)
  PWinv <- backsolve(PW, diag(ncomp))
  
  #  R = W (P^T W)^{-1}
  R_proj <- W_weights %*% PWinv
  
  B <- R_proj %*% t(R_loadings)
  
  # 
  return(list(
    coefficients = B,
    loading.weights = W_weights,
    loadings = P_loadings,
    Yloadings = R_loadings,
    scores = T_scores,
    Projection = R_proj,
    normalize_terms = normalize_terms
  ))
}


# accelerate the computation by using matrix-free method
acc_pca_pls_with_alpha <- function(X, Y, alpha, ncomp,
                                   normalize_terms = FALSE,
                                   norm_eps = 1e-12) {      
  n <- nrow(X)
  p <- ncol(X)
  if (is.vector(Y)) Y <- matrix(Y, ncol = 1)
  q <- ncol(Y)  # Y 
  
  # 
  T_scores <- matrix(0, n, ncomp)  #  X-scores
  P_loadings <- matrix(0, p, ncomp)  #  X-loadings
  R_loadings <- matrix(0, q, ncomp)  #  y-loadings ( q  ncomp )
  W_weights <- matrix(0, p, ncomp)  #  w  (p  ncomp )
  pca_norms <- numeric(ncomp)
  pls_norms <- numeric(ncomp)
  
  for (i in 1:ncomp) {
    #print(paste("H matrix: ", i))
    if (isTRUE(normalize_terms)) {
      pca_norms[i] <- .pcapls_pca_norm(X, eps = norm_eps)
      pls_norms[i] <- .pcapls_pls_norm(X, Y, eps = norm_eps)
    } else {
      pca_norms[i] <- 1
      pls_norms[i] <- 1
    }
    extra_list <- list(
      X = X,
      Y = Y,
      alpha = alpha,
      pca_scale = pca_norms[i],
      pls_scale = pls_norms[i]
    )
    #print(str(extra_list))
    #  H*v H 
    Af <- function(v, args) {
      #  extra  X, Y, alpha
    
      X <- args$X
      Y <- args$Y
      alpha <- args$alpha
      pca_scale <- args$pca_scale
      pls_scale <- args$pls_scale

      # :  X %*% v
      tmp <- X %*% v

      # :  Y^T * (Xv) Y %*% (Y^T * (Xv))
      tmp2 <- ((1 - alpha) / pca_scale) * tmp +
        (alpha / pls_scale) * (Y %*% (crossprod(Y, tmp)))

      # :  X^T * tmp2 H*v
      return(crossprod(X, tmp2))
    }
    
    eigen_decomp <- tryCatch(
      eigs_sym(
        Af,
        k = 1,
        n = ncol(X),
        which = "LM",
        args = extra_list  # 
      ),
      error = function(e) {
        # RSpectra can fail on some nearly singular problems; fall back to
        # an explicit symmetric eigen decomposition for robustness.
        H <- ((1 - alpha) / pca_norms[i]) * crossprod(X) +
          (alpha / pls_norms[i]) * crossprod(X, Y) %*% crossprod(Y, X)
        H <- (H + t(H)) / 2
        eig <- eigen(H, symmetric = TRUE)
        list(
          values = eig$values[1],
          vectors = matrix(eig$vectors[, 1], ncol = 1)
        )
      }
    )


    # Step 2:  w
    # eigen_decomp <- eigs_sym(Af, k = 1, n = ncol(X), which = "LM",
    #                          extra = list(X = X, Y = Y, alpha = alpha))
  
    w <- eigen_decomp$vectors[, 1]  
    W_weights[, i] <- w  #  w 
    
    # Step 3:  X-scores: t = X %*% w
    t_score <- X %*% w
    
    # Step 4:  X-loadings: p = X' t / (t' t)
    p_loading <- crossprod(X, t_score) / as.numeric(crossprod(t_score))
    
    # Step 5:  y-loadings: r = Y' t / (t' t)
    r_loading <- crossprod(Y, t_score) / as.numeric(crossprod(t_score))  #  q x 1
    
    # Step 6: Deflate X  Y: X = X - t p', Y = Y - t r'
    X <- X - t_score %*% t(p_loading)
    Y <- Y - t_score %*% t(r_loading)  # Y  r_loading  q x 1
    
    # 
    T_scores[, i] <- t_score
    P_loadings[, i] <- p_loading
    R_loadings[, i] <- r_loading
  }
  
  # Step 7:  R
  #  P^T W
  PW <- crossprod(P_loadings, W_weights)
  
  #  (P^T W) 
  #PWinv <- solve(PW)
  PWinv <- backsolve(PW, diag(ncomp))
  
  #  R = W (P^T W)^{-1}
  R_proj <- W_weights %*% PWinv
  
  B <- R_proj %*% t(R_loadings)
  
  # 
  return(list(
    coefficients = B,
    loading.weights = W_weights,
    loadings = P_loadings,
    Yloadings = R_loadings,
    scores = T_scores,
    Projection = R_proj,
    normalize_terms = normalize_terms,
    pca_norms = pca_norms,
    pls_norms = pls_norms
  ))
}

acc_pca_pls_with_alpha_raw <- function(X, Y, alpha, ncomp) {
  acc_pca_pls_with_alpha(
    X = X,
    Y = Y,
    alpha = alpha,
    ncomp = ncomp,
    normalize_terms = FALSE
  )
}

acc_pca_pls_with_alpha_normalized <- function(X, Y, alpha, ncomp,
                                              norm_eps = 1e-12) {
  acc_pca_pls_with_alpha(
    X = X,
    Y = Y,
    alpha = alpha,
    ncomp = ncomp,
    normalize_terms = TRUE,
    norm_eps = norm_eps
  )
}
