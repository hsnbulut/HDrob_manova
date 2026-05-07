hdmanova_schott <- function(X, group) {
  #------------------------------------------------------------
  # Schott (2007) high-dimensional one-way MANOVA test
  #
  # Input:
  #   X     : numeric matrix/data.frame, n x p
  #   group : group labels of length n
  #
  # Output:
  #   list with statistic, z.value, p.value, B, W, etc.
  #------------------------------------------------------------
  
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }
  
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`X` must be a numeric matrix or coercible data.frame.")
  }
  
  if (anyNA(X)) {
    stop("`X` contains missing values. Please remove or impute them first.")
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
  
  df_between <- G - 1L
  df_within  <- n - G
  
  if (df_within <= 1L) {
    stop("Need n - G > 1 for Schott's variance estimator.")
  }
  
  #------------------------------------------------------------
  # Group means and overall mean
  #------------------------------------------------------------
  overall_mean <- colMeans(X)
  
  group_means <- t(vapply(
    split.data.frame(as.data.frame(X), group),
    FUN = function(dfk) colMeans(as.matrix(dfk)),
    FUN.VALUE = numeric(p)
  ))
  
  #------------------------------------------------------------
  # Between-group scatter matrix B
  # B = sum_k n_k (xbar_k - xbar)(xbar_k - xbar)^T
  #------------------------------------------------------------
  B <- matrix(0, nrow = p, ncol = p)
  for (k in seq_len(G)) {
    d <- group_means[k, ] - overall_mean
    B <- B + nk[k] * tcrossprod(d)
  }
  
  #------------------------------------------------------------
  # Within-group scatter matrix W
  # W = sum_k sum_i (x_ki - xbar_k)(x_ki - xbar_k)^T
  #------------------------------------------------------------
  W <- matrix(0, nrow = p, ncol = p)
  split_X <- split.data.frame(as.data.frame(X), group)
  
  for (k in seq_len(G)) {
    Xk <- as.matrix(split_X[[k]])
    centered_k <- sweep(Xk, 2L, group_means[k, ], FUN = "-")
    W <- W + crossprod(centered_k)
  }
  
  trB <- sum(diag(B))
  trW <- sum(diag(W))
  trW2 <- sum(W * W)  # trace(W %*% W)
  
  #------------------------------------------------------------
  # Schott test statistic
  #------------------------------------------------------------
  T_S <- (1 / sqrt(n - 1)) * ((trB / df_between) - (trW / df_within))
  
  a_hat_num <- trW2 - (trW^2 / df_within)
  a_hat_den <- (df_within + 2) * (df_within - 1)
  a_hat <- a_hat_num / a_hat_den
  
  if (!is.finite(a_hat) || a_hat <= 0) {
    stop("Estimated variance component a_hat is non-positive or not finite.")
  }
  
  varhat_T_S <- 2 * a_hat / (df_between * df_within)
  
  if (!is.finite(varhat_T_S) || varhat_T_S <= 0) {
    stop("Estimated variance of T_S is non-positive or not finite.")
  }
  
  Z_S <- T_S / sqrt(varhat_T_S)
  
  # Right-tailed test
  p_value <- 1 - pnorm(Z_S)
  
  out <- list(
    method = "Schott (2007) high-dimensional one-way MANOVA test",
    statistic = unname(T_S),
    z.value = unname(Z_S),
    p.value = unname(p_value),
    a_hat = unname(a_hat),
    varhat = unname(varhat_T_S),
    B = B,
    W = W,
    trB = unname(trB),
    trW = unname(trW),
    n = n,
    p = p,
    G = G,
    nk = nk,
    df_between = df_between,
    df_within = df_within,
    call = match.call()
  )
  
  class(out) <- "hdmanova_schott"
  out
}