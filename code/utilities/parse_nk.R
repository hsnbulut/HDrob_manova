parse_nk <- function(x) {
  if (length(x) != 1L || !is.character(x)) {
    stop("`x` must be a single character string such as '20,20,20'.")
  }
  
  parts <- strsplit(x, ",", fixed = TRUE)[[1L]]
  parts <- trimws(parts)
  
  if (length(parts) == 0L) {
    stop("No group sizes found in `x`.")
  }
  
  nk <- suppressWarnings(as.integer(parts))
  
  if (anyNA(nk)) {
    stop("All group sizes in `x` must be valid integers.")
  }
  
  if (any(nk <= 0L)) {
    stop("All group sizes must be positive.")
  }
  
  nk
}