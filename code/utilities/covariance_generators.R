make_sigma <- function(p, structure = c("identity", "ar1"), rho = 0.5) {
  structure <- match.arg(structure)
  
  if (!is.numeric(p) || length(p) != 1L || is.na(p) || p < 1) {
    stop("`p` must be a positive integer.")
  }
  p <- as.integer(p)
  
  if (structure == "identity") {
    return(diag(p))
  }
  
  if (!is.numeric(rho) || length(rho) != 1L || is.na(rho) || abs(rho) >= 1) {
    stop("For AR(1), `rho` must be a single numeric value with |rho| < 1.")
  }
  
  idx <- seq_len(p)
  Sigma <- rho ^ abs(outer(idx, idx, "-"))
  Sigma
}