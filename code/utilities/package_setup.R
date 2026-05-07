ensure_package <- function(pkg,
                           lib = Sys.getenv("R_LIBS_USER"),
                           repos = "https://cloud.r-project.org",
                           quietly = TRUE) {
  if (!is.character(pkg) || length(pkg) != 1L || !nzchar(pkg)) {
    stop("`pkg` must be a single non-empty character string.")
  }
  
  if (!nzchar(lib)) {
    lib <- path.expand("~/R_libs")
  }
  
  if (!dir.exists(lib)) {
    dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  }
  
  currently_available <- requireNamespace(pkg, quietly = quietly)
  
  if (!currently_available) {
    message(sprintf(
      "Package '%s' is not available. Attempting installation into: %s",
      pkg, lib
    ))
    
    tryCatch(
      utils::install.packages(pkg, lib = lib, repos = repos, quiet = quietly),
      error = function(e) {
        stop(sprintf(
          "Automatic installation of package '%s' failed.\nOriginal error: %s",
          pkg, conditionMessage(e)
        ))
      }
    )
  }
  available_after_install <- requireNamespace(pkg, quietly = quietly, lib.loc = lib)

  if (!available_after_install) {
    stop(sprintf(
      "Package '%s' is still unavailable after installation attempt.",
      pkg
    ))
  }
  
  invisible(TRUE)
}

ensure_required_packages <- function() {
  required <- c(
    "rrcov",
    "MASS",
    "data.table",
    "optparse",
    "here",
    "Gmedian",
    "extraDistr"
  )
  
  for (pkg in required) {
    ensure_package(pkg)
  }
  
  invisible(TRUE)
}

check_required_packages <- function() {
  required <- c(
    "rrcov",
    "MASS",
    "data.table",
    "optparse",
    "here",
    "Gmedian",
    "extraDistr"
  )
  
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  
  if (length(missing) > 0) {
    stop(
      paste(
        "Missing required packages:",
        paste(missing, collapse = ", "),
        "\nRun scripts/setup_truba_packages.R on login node first, or install them manually."
      )
    )
  }
  
  invisible(TRUE)
}