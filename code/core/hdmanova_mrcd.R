hdmanova_mrcd <- function(X,
                          group,
                          eta = 0.75,
                          B = 199,
                          alpha = 0.05,
                          seed = NULL,
                          ridge_eps = 1e-8) {
  #------------------------------------------------------------
  # Weighted MRCD-based robust Wilks' Lambda test
  #
  # Implements the algorithm described in the manuscript:
  # 1) Group-wise MRCD estimates (t_k, Sigma_k)
  # 2) Pool centered data z_ki = x_ki - t_k
  # 3) Pooled MRCD estimates (delta_hat, C_hat)
  # 4) Preliminary centers mu_tilde_k = t_k + delta_hat
  # 5) Robust distances + normal-calibrated weights
  # 6) Weighted means, WR, BR, TR
  # 7) Lambda_R = |WR| / |TR|
  # 8) Permutation p-value with inclusive correction
  #
  # Input:
  #   X         : n x p numeric matrix
  #   group     : factor / vector of group labels
  #   eta       : MRCD trimming parameter
  #   B         : number of permutations
  #   alpha     : significance level used both in weighting and test decision
  #   seed      : optional seed for permutation reproducibility
  #   ridge_eps : tiny ridge added only if needed for numerical stability
  #
  # Output:
  #   list with statistic, p.value, WR, BR, weights, etc.
  #------------------------------------------------------------
  
  #-------------------------
  # checks
  #-------------------------

  
  if (is.data.frame(X)) X <- as.matrix(X)
  
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`X` must be a numeric matrix or coercible data.frame.")
  }
  if (anyNA(X)) {
    stop("`X` contains missing values.")
  }
  
  n <- nrow(X)
  p <- ncol(X)
  
  if (length(group) != n) {
    stop("Length of `group` must equal nrow(X).")
  }
  
  group <- as.factor(group)
  G <- nlevels(group)
  
  if (G < 2L) {
    stop("At least two groups are required.")
  }
  
  nk <- as.integer(table(group))
  if (any(nk < 2L)) {
    stop("Each group must contain at least 2 observations.")
  }
  
  if (!is.numeric(eta) || length(eta) != 1L || is.na(eta) || eta <= 0 || eta >= 1) {
    stop("`eta` must be a single number in (0, 1).")
  }
  
  if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 1) {
    stop("`B` must be a positive integer.")
  }
  B <- as.integer(B)
  
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).")
  }
  
  #-------------------------
  # package checks
  #-------------------------
  if (!requireNamespace("rrcov", quietly = TRUE)) {
    stop("Package 'rrcov' is required for CovMrcd(). Please install it first.")
  }
  
  if (!requireNamespace("rrcov", quietly = TRUE)) {
    stop("Package 'rrcov' is required. Run scripts/setup_truba_packages.R first.")
  }
  
  #-------------------------
  # internal helpers
  #-------------------------
  make_spd <- function(S, eps = ridge_eps) {
    S <- (S + t(S)) / 2
    ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    min_ev <- min(ev)
    
    if (!is.finite(min_ev)) {
      stop("Non-finite eigenvalue encountered while stabilizing a matrix.")
    }
    
    if (min_ev <= eps) {
      S <- S + diag(abs(min_ev) + eps, nrow(S))
    }
    (S + t(S)) / 2
  }
  
  logdet_spd <- function(S) {
    S <- make_spd(S)
    as.numeric(determinant(S, logarithm = TRUE)$modulus)
  }
  
  fit_mrcd <- function(Xmat, eta) {
    fit <- tryCatch(
      rrcov::CovMrcd(Xmat, alpha = eta),
      error = function(e) {
        stop("MRCD fit failed: ", conditionMessage(e))
      }
    )
    fit
  }
  
  extract_center <- function(fit) {
    if (exists("getCenter", where = asNamespace("rrcov"), inherits = FALSE)) {
      cen <- tryCatch(rrcov::getCenter(fit), error = function(e) NULL)
      if (!is.null(cen)) return(as.numeric(cen))
    }
    if (!is.null(fit@center)) return(as.numeric(fit@center))
    stop("Could not extract MRCD center.")
  }
  
  extract_cov <- function(fit) {
    if (exists("getCov", where = asNamespace("rrcov"), inherits = FALSE)) {
      cv <- tryCatch(rrcov::getCov(fit), error = function(e) NULL)
      if (!is.null(cv)) return(as.matrix(cv))
    }
    if (!is.null(fit@cov)) return(as.matrix(fit@cov))
    stop("Could not extract MRCD covariance.")
  }
  
  standardized_eigenvalues <- function(Z) {
    # Z is pooled centered data matrix
    # We standardize columns and use eigenvalues of the resulting correlation-type matrix.
    Zs <- scale(Z, center = TRUE, scale = FALSE)
    sds <- apply(Zs, 2L, stats::sd)
    
    # protect zero or nearly zero scales
    sds[!is.finite(sds) | sds < sqrt(.Machine$double.eps)] <- 1
    
    Zstd <- sweep(Zs, 2L, sds, "/")
    Rhat <- stats::cov(Zstd)
    Rhat <- make_spd(Rhat)
    lam <- eigen(Rhat, symmetric = TRUE, only.values = TRUE)$values
    lam[!is.finite(lam)] <- 0
    lam
  }
  
  compute_weights <- function(Xmat, mu_tilde, C_hat, lam, alpha) {
    C_hat <- make_spd(C_hat)
    C_inv <- solve(C_hat)
    
    centered <- sweep(Xmat, 2L, mu_tilde, "-")
    RD2 <- rowSums((centered %*% C_inv) * centered)
    
    zq <- qnorm(1 - alpha)
    
    tr2 <- sum(lam^2) - (p^2 / n)
    tr2 <- max(tr2, .Machine$double.eps)
    
    cpn <- 1 + sum(lam^2) / (p^(3/2))
    cpn <- max(cpn, 1)
    
    s_alpha <- 1 + (exp(-(zq^2) / 2) / ((1 - alpha) * sqrt(pi))) * sqrt(tr2 / p)
    
    Zscore <- (RD2 / s_alpha - p) / sqrt(2 * tr2 * cpn)
    
    w <- ifelse(Zscore <= zq, 1L, 0L)
    
    list(
      RD2 = as.numeric(RD2),
      Zscore = as.numeric(Zscore),
      weight = as.integer(w),
      tr2 = as.numeric(tr2),
      cpn = as.numeric(cpn),
      s_alpha = as.numeric(s_alpha)
    )
  }
  
  compute_lambda_once <- function(group_labels, return_details = TRUE) {
    group_labels <- as.factor(group_labels)
    split_idx <- split(seq_len(n), group_labels)
    nk_loc <- as.integer(lengths(split_idx))
    G_loc <- length(split_idx)
    
    #--------------------------------------------
    # Step 1: groupwise MRCD
    #--------------------------------------------
    t_list <- vector("list", G_loc)
    S_list <- vector("list", G_loc)
    
    for (k in seq_len(G_loc)) {
      Xk <- X[split_idx[[k]], , drop = FALSE]
      fit_k <- fit_mrcd(Xk, eta = eta)
      t_list[[k]] <- extract_center(fit_k)
      S_list[[k]] <- make_spd(extract_cov(fit_k))
    }
    
    #--------------------------------------------
    # Step 2: pooled centered data
    # z_ki = x_ki - t_k
    #--------------------------------------------
    Z <- X
    for (k in seq_len(G_loc)) {
      idx <- split_idx[[k]]
      Z[idx, ] <- sweep(X[idx, , drop = FALSE], 2L, t_list[[k]], "-")
    }
    
    #--------------------------------------------
    # Step 3: pooled MRCD on centered data
    #--------------------------------------------
    fit_pool <- fit_mrcd(Z, eta = eta)
    delta_hat <- extract_center(fit_pool)
    C_hat <- make_spd(extract_cov(fit_pool))
    
    #--------------------------------------------
    # Step 4: preliminary centers
    # mu_tilde_k = t_k + delta_hat
    #--------------------------------------------
    mu_tilde_list <- lapply(t_list, function(tk) tk + delta_hat)
    
    #--------------------------------------------
    # Step 5: compute robust distances and weights
    #--------------------------------------------
    lam <- standardized_eigenvalues(Z)
    
    weight_vec <- integer(n)
    RD2_vec <- numeric(n)
    Zscore_vec <- numeric(n)
    
    v_k <- integer(G_loc)
    mu_hat_list <- vector("list", G_loc)
    
    for (k in seq_len(G_loc)) {
      idx <- split_idx[[k]]
      Xk <- X[idx, , drop = FALSE]
      
      ww <- compute_weights(
        Xmat = Xk,
        mu_tilde = mu_tilde_list[[k]],
        C_hat = C_hat,
        lam = lam,
        alpha = alpha
      )
      
      weight_vec[idx] <- ww$weight
      RD2_vec[idx] <- ww$RD2
      Zscore_vec[idx] <- ww$Zscore
      
      v_k[k] <- sum(ww$weight)
      
      # fallback rule
      if (v_k[k] > 1L) {
        mu_hat_list[[k]] <- colSums(Xk * ww$weight) / v_k[k]
      } else {
        mu_hat_list[[k]] <- t_list[[k]]
        v_k[k] <- nk_loc[k]
      }
    }
    
    #--------------------------------------------
    # Step 6: overall weighted mean
    #--------------------------------------------
    total_v <- sum(v_k)
    mu_hat <- Reduce(`+`, Map(function(v, m) v * m, v_k, mu_hat_list)) / total_v
    
    #--------------------------------------------
    # Step 7: WR and BR
    # WR = sum (v_k - 1) Sigma_k
    # BR = sum v_k (mu_k - mu)(mu_k - mu)^T
    #--------------------------------------------
    WR <- matrix(0, nrow = p, ncol = p)
    BR <- matrix(0, nrow = p, ncol = p)
    
    for (k in seq_len(G_loc)) {
      WR <- WR + (v_k[k] - 1) * S_list[[k]]
      d <- mu_hat_list[[k]] - mu_hat
      BR <- BR + v_k[k] * tcrossprod(d)
    }
    
    WR <- make_spd(WR)
    TR <- make_spd(WR + BR)
    
    logdet_WR <- logdet_spd(WR)
    logdet_TR <- logdet_spd(TR)
    
    Lambda_R <- exp(logdet_WR - logdet_TR)
    
    if (!return_details) {
      return(Lambda_R)
    }
    
    list(
      statistic = as.numeric(Lambda_R),
      WR = WR,
      BR = BR,
      TR = TR,
      group_centers_mrcd = t_list,
      pooled_delta = delta_hat,
      pooled_scatter = C_hat,
      mu_tilde = mu_tilde_list,
      mu_hat = mu_hat_list,
      overall_mu_hat = mu_hat,
      weights = weight_vec,
      v_k = v_k,
      RD2 = RD2_vec,
      Zscore = Zscore_vec
    )
  }
  
  #-------------------------
  # observed statistic
  #-------------------------
  obs <- compute_lambda_once(group, return_details = TRUE)
  lambda_obs <- obs$statistic
  
  #-------------------------
  # permutation distribution
  # reject for small Lambda
  # inclusive p-value:
  # (1 + #{Lambda_b <= Lambda_obs}) / (B + 1)
  #-------------------------
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  lambda_perm <- numeric(B)
  
  for (b in seq_len(B)) {
    g_perm <- sample(group, size = n, replace = FALSE)
    lambda_perm[b] <- compute_lambda_once(g_perm, return_details = FALSE)
  }
  
  p_perm <- (1 + sum(lambda_perm <= lambda_obs)) / (B + 1)
  
  out <- list(
    method = "Weighted MRCD-based robust Wilks' Lambda test",
    statistic = as.numeric(lambda_obs),
    p.value = as.numeric(p_perm),
    perm_stats = lambda_perm,
    WR = obs$WR,
    BR = obs$BR,
    TR = obs$TR,
    weights = obs$weights,
    v_k = obs$v_k,
    RD2 = obs$RD2,
    Zscore = obs$Zscore,
    group_centers_mrcd = obs$group_centers_mrcd,
    pooled_delta = obs$pooled_delta,
    pooled_scatter = obs$pooled_scatter,
    mu_tilde = obs$mu_tilde,
    mu_hat = obs$mu_hat,
    overall_mu_hat = obs$overall_mu_hat,
    n = n,
    p = p,
    G = G,
    nk = nk,
    eta = eta,
    B = B,
    alpha = alpha,
    call = match.call()
  )
  
  class(out) <- "hdmanova_mrcd"
  out
}