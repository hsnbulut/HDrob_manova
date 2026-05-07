hdmanova_cheng <- function(X,
                           group,
                           B = 199,
                           alpha = 0.05,
                           standardized = TRUE) {
  #------------------------------------------------------------
  # Cheng et al. (2024)-style geometric median-based
  # high-dimensional MANOVA test with bootstrap calibration.
  #
  # This implementation is adapted from the user-provided
  # HD_MANOVA_GM_functions.R code, but only retains the
  # geometric median test and wraps it into a stable interface.
  #
  # Input:
  #   X            : n x p numeric matrix
  #   group        : factor / vector of group labels
  #   B            : number of bootstrap iterations
  #   alpha        : significance level
  #   standardized : if TRUE, use standardized GM statistic
  #
  # Output:
  #   list(method, statistic, p.value, bootstrap_stats, raw, call)
  #------------------------------------------------------------
  
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }
  
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
  K <- nlevels(group)
  
  if (K < 2L) {
    stop("At least two groups are required.")
  }
  
  if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 1) {
    stop("`B` must be a positive integer.")
  }
  B <- as.integer(B)
  
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).")
  }
  
  if (!is.logical(standardized) || length(standardized) != 1L || is.na(standardized)) {
    stop("`standardized` must be TRUE or FALSE.")
  }
  
  #------------------------------------------------------------
  # package checks
  #------------------------------------------------------------
  if (!requireNamespace("Gmedian", quietly = TRUE)) {
    stop("Package 'Gmedian' is required. Run scripts/setup_truba_packages.R first.")
  }
  if (!requireNamespace("extraDistr", quietly = TRUE)) {
    stop("Package 'extraDistr' is required. Run scripts/setup_truba_packages.R first.")
  }
  
  #------------------------------------------------------------
  # helpers
  #------------------------------------------------------------
  split_dat <- split.data.frame(as.data.frame(X), group)
  dat <- lapply(split_dat, as.matrix)
  N <- vapply(dat, nrow, integer(1))
  p_dim <- ncol(X)
  
  if (any(N < 2L)) {
    stop("Each group must contain at least 2 observations.")
  }
  
  geom_median <- function(M) {
    as.numeric(Gmedian::Weiszfeld(M)$median)
  }
  
  safe_diag_sqrt <- function(M) {
    d <- diag(M)
    d[!is.finite(d) | d <= sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
    sqrt(d)
  }
  
  compute_group_medians <- function(dat_list) {
    lapply(dat_list, geom_median)
  }
  
  compute_group_median_cov <- function(Xk, med_k) {
    # Direct adaptation of the user-provided estimator
    Xk_central <- Xk - rep(1, nrow(Xk)) %o% med_k
    norm_Xk <- sqrt(diag(Xk_central %*% t(Xk_central)))
    
    # protect against division by zero
    norm_Xk[!is.finite(norm_Xk) | norm_Xk < sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
    
    c_Xk <- mean(1 / norm_Xk)
    B_Xk <- (t(Xk_central) %*% diag(1 / norm_Xk)) %*%
      t(t(Xk_central) %*% diag(1 / norm_Xk)) / nrow(Xk)
    
    B_Xk / (c_Xk^2)
  }
  
  compute_unstd_stat <- function(median_list, N) {
    K <- length(median_list)
    dif_matrix <- matrix(0, K, K)
    
    for (k1 in 1:(K - 1)) {
      for (k2 in (k1 + 1):K) {
        tmp <- abs(median_list[[k1]] - median_list[[k2]])
        dif_matrix[k1, k2] <- max(sqrt(N[k1] * N[k2] / (N[k1] + N[k2])) * tmp)
      }
    }
    max(dif_matrix)
  }
  
  compute_std_stat <- function(median_list, median_sd_est, N) {
    K <- length(median_list)
    dif_matrix_std <- matrix(0, K, K)
    
    for (k1 in 1:(K - 1)) {
      for (k2 in (k1 + 1):K) {
        Bxy <- (N[k2] * median_sd_est[[k1]] + N[k1] * median_sd_est[[k2]]) / (N[k1] + N[k2])
        bxy <- safe_diag_sqrt(Bxy)
        tmp <- abs(median_list[[k1]] - median_list[[k2]]) / bxy
        dif_matrix_std[k1, k2] <- max(sqrt(N[k1] * N[k2] / (N[k1] + N[k2])) * tmp)
      }
    }
    max(dif_matrix_std)
  }
  
  #------------------------------------------------------------
  # observed statistic
  #------------------------------------------------------------
  median_list <- compute_group_medians(dat)
  
  median_sd_est <- vector("list", K)
  for (k in seq_len(K)) {
    median_sd_est[[k]] <- compute_group_median_cov(dat[[k]], median_list[[k]])
  }
  
  test_stat_median <- compute_unstd_stat(median_list, N)
  test_stat_median_std <- compute_std_stat(median_list, median_sd_est, N)
  
  #------------------------------------------------------------
  # bootstrap
  #------------------------------------------------------------
  res_median_boot <- numeric(B)
  res_median_std_boot <- numeric(B)
  
  for (b in seq_len(B)) {
    dat_boot <- vector("list", K)
    median_boot_list <- vector("list", K)
    
    for (k in seq_len(K)) {
      Z_Xk <- extraDistr::rsign(N[k]) %o% rep(1, p_dim)
      dat_boot[[k]] <- (dat[[k]] - rep(1, N[k]) %o% median_list[[k]]) * Z_Xk
      median_boot_list[[k]] <- geom_median(dat_boot[[k]])
    }
    
    res_median_boot[b] <- compute_unstd_stat(median_boot_list, N)
    
    median_sd_est_boot <- vector("list", K)
    for (k in seq_len(K)) {
      median_sd_est_boot[[k]] <- compute_group_median_cov(dat_boot[[k]], median_boot_list[[k]])
    }
    
    res_median_std_boot[b] <- compute_std_stat(median_boot_list, median_sd_est_boot, N)
  }
  
  statistic <- if (standardized) test_stat_median_std else test_stat_median
  boot_stats <- if (standardized) res_median_std_boot else res_median_boot
  
  # right-tailed inclusive bootstrap p-value
  p_value <- (1 + sum(boot_stats >= statistic)) / (B + 1)
  
  raw <- list(
    test_stat_median = test_stat_median,
    test_stat_median_std = test_stat_median_std,
    res_median_boot = res_median_boot,
    res_median_std_boot = res_median_std_boot,
    median_list = median_list,
    median_sd_est = median_sd_est
  )
  
  out <- list(
    method = if (standardized) {
      "Cheng et al. (2024) geometric median-based HD-MANOVA test (standardized)"
    } else {
      "Cheng et al. (2024) geometric median-based HD-MANOVA test (unstandardized)"
    },
    statistic = as.numeric(statistic),
    p.value = as.numeric(p_value),
    bootstrap_stats = boot_stats,
    raw = raw,
    standardized = standardized,
    B = B,
    alpha = alpha,
    n = n,
    p = p,
    K = K,
    N = N,
    call = match.call()
  )
  
  class(out) <- "hdmanova_cheng"
  out
}