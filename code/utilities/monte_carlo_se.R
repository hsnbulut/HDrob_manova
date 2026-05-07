mcse_binomial <- function(rate, R) {
  if (!is.numeric(rate) || length(rate) != 1L || is.na(rate) || rate < 0 || rate > 1) {
    stop("`rate` must be a single number in [0, 1].")
  }
  if (!is.numeric(R) || length(R) != 1L || is.na(R) || R <= 0) {
    stop("`R` must be a positive integer.")
  }
  
  R <- as.integer(R)
  sqrt(rate * (1 - rate) / R)
}