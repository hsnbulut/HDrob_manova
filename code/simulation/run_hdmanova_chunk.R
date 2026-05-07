suppressPackageStartupMessages({
  library(data.table)
  library(MASS)
  library(optparse)
  library(here)
  library(parallel)
})

source(here::here("code", "utilities", "package_setup.R"))
check_required_packages()

source(here::here("code", "utilities", "parse_nk.R"))
source(here::here("code", "utilities", "covariance_generators.R"))
source(here::here("code", "utilities", "monte_carlo_se.R"))

source(here::here("code", "core", "hdmanova_mrcd.R"))
source(here::here("code", "competitors", "hdmanova_schott.R"))
source(here::here("code", "competitors", "hdmanova_cheng.R"))

option_list <- list(
  make_option("--grid", type = "character", dest = "grid"),
  make_option("--study", type = "character", dest = "study"),
  make_option("--task-id", type = "integer", dest = "task_id"),
  make_option("--outdir", type = "character", dest = "outdir"),
  make_option("--R", type = "integer", default = 1000L, dest = "R"),
  make_option("--R-block", type = "integer", default = 100L, dest = "R_block"),
  make_option("--B", type = "integer", default = 199L, dest = "B"),
  make_option("--alpha", type = "double", default = 0.05, dest = "alpha"),
  make_option("--seed", type = "integer", default = 12345L, dest = "seed"),
  make_option("--eta", type = "double", default = 0.75, dest = "eta"),
  make_option("--delta-rob", type = "double", default = 4.0, dest = "delta_rob"),
  make_option("--rho", type = "double", default = 0.5, dest = "rho"),
  make_option("--mc-cores", type = "integer", default = 1L, dest = "mc_cores")
)

opt <- parse_args(OptionParser(option_list = option_list))

cat("Parsed options:\n")
print(opt)

if (is.null(opt$grid) || is.null(opt$study) || is.null(opt$outdir) || is.null(opt$task_id)) {
  stop("Arguments --grid, --study, --task-id and --outdir are required.")
}

if (opt$R %% opt$R_block != 0L) {
  stop("`R` must be divisible by `R_block`.")
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

#------------------------------------------------------------
# Data generators
#------------------------------------------------------------
rmvn <- function(n, mu, Sigma) {
  MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
}

generate_type1_data <- function(nk, p, Sigma) {
  G <- length(nk)
  X_list <- vector("list", G)
  g <- integer(sum(nk))
  pos <- 1L
  
  for (k in seq_len(G)) {
    Xk <- rmvn(nk[k], rep(0, p), Sigma)
    X_list[[k]] <- Xk
    g[pos:(pos + nk[k] - 1L)] <- k
    pos <- pos + nk[k]
  }
  
  list(X = do.call(rbind, X_list), group = factor(g))
}

generate_power_data <- function(nk, p, Sigma, d) {
  G <- length(nk)
  X_list <- vector("list", G)
  g <- integer(sum(nk))
  pos <- 1L
  
  mu_list <- vector("list", G)
  mu_list[[1L]] <- rep(0, p)
  
  if (G == 2L) {
    mu_list[[2L]] <- rep(d, p)
  } else if (G == 3L) {
    mu_list[[2L]] <- rep(d, p)
    mu_list[[3L]] <- rep(2 * d, p)
  } else {
    stop("Only G = 2 or G = 3 is supported.")
  }
  
  for (k in seq_len(G)) {
    Xk <- rmvn(nk[k], mu_list[[k]], Sigma)
    X_list[[k]] <- Xk
    g[pos:(pos + nk[k] - 1L)] <- k
    pos <- pos + nk[k]
  }
  
  list(X = do.call(rbind, X_list), group = factor(g))
}

generate_robustness_data <- function(nk, p, Sigma, eps, delta) {
  G <- length(nk)
  q <- floor(sqrt(p))
  mu_pos <- c(rep(delta, q), rep(0, p - q))
  mu_neg <- c(rep(-delta, q), rep(0, p - q))
  mu_zero <- rep(0, p)
  Sigma_out <- 0.25 * Sigma
  
  X_list <- vector("list", G)
  g <- integer(sum(nk))
  pos <- 1L
  
  for (k in seq_len(G)) {
    n_k <- nk[k]
    n_out <- floor(eps * n_k)
    n_clean <- n_k - n_out
    
    X_clean <- rmvn(n_clean, mu_zero, Sigma)
    
    if (G == 2L) {
      mu_out <- if (k == 1L) mu_pos else mu_neg
    } else if (G == 3L) {
      mu_out <- if (k == 1L) mu_pos else if (k == 2L) mu_neg else mu_zero
    } else {
      stop("Only G = 2 or G = 3 is supported.")
    }
    
    if (n_out > 0L) {
      X_out <- rmvn(n_out, mu_out, Sigma_out)
      Xk <- rbind(X_clean, X_out)
    } else {
      Xk <- X_clean
    }
    
    if (nrow(Xk) > 1L) {
      Xk <- Xk[sample.int(nrow(Xk)), , drop = FALSE]
    }
    
    X_list[[k]] <- Xk
    g[pos:(pos + n_k - 1L)] <- k
    pos <- pos + n_k
  }
  
  list(X = do.call(rbind, X_list), group = factor(g))
}

#------------------------------------------------------------
# Safe extraction
#------------------------------------------------------------
safe_pvalue <- function(obj) {
  if (is.list(obj) && !is.null(obj$p.value)) {
    return(as.numeric(obj$p.value))
  }
  NA_real_
}

#------------------------------------------------------------
# One replication
#------------------------------------------------------------
run_one_rep <- function(row, rep_seed, opt) {
  set.seed(rep_seed)
  
  nk <- parse_nk(row$nk)
  p <- as.integer(row$p)
  Sigma <- make_sigma(p = p, structure = row$cov_type, rho = opt$rho)

  # For robustness scenarios, use the contamination-adaptive MRCD trimming rule
  # discussed in the manuscript: eta = 1 - 1.5 * eps.
  # This preserves the existing fixed eta behavior for type-I and power studies.
  
  # eta_current <- opt$eta
  # if (identical(as.character(row$study), "robustness")) {
  #   eta_current <- 1 - 1.5 * as.numeric(row$eps)
  #   eta_current <- max(0.55, min(0.95, eta_current))
  # }
  
  # For robustness scenarios, use the contamination-adaptive MRCD trimming rule
  # discussed in the manuscript: eta = 1 - 2 * eps.
  # This preserves the existing fixed eta behavior for type-I and power studies.
  
  eta_current <- opt$eta
  if (identical(as.character(row$study), "robustness")) {
    eta_current <- 1 - 2.0 * as.numeric(row$eps)
    eta_current <- max(0.55, min(0.95, eta_current))
  }
  
  dat <- switch(
    row$study,
    type1 = generate_type1_data(nk, p, Sigma),
    power = generate_power_data(nk, p, Sigma, d = row$signal_d),
    robustness = generate_robustness_data(nk, p, Sigma, eps = row$eps, delta = opt$delta_rob),
    stop("Unknown study type.")
  )
  
  res_prop <- tryCatch(
    hdmanova_mrcd(dat$X, dat$group, eta = eta_current, B = opt$B, alpha = opt$alpha),
    error = function(e) list(p.value = NA_real_, error = conditionMessage(e))
  )
  
  res_schott <- tryCatch(
    hdmanova_schott(dat$X, dat$group),
    error = function(e) list(p.value = NA_real_, error = conditionMessage(e))
  )
  
  res_cheng <- tryCatch(
    hdmanova_cheng(dat$X, dat$group, B = opt$B, alpha = opt$alpha, standardized = TRUE),
    error = function(e) list(p.value = NA_real_, error = conditionMessage(e))
  )
  
  data.table(
    rep_id = NA_integer_,
    p_proposed = safe_pvalue(res_prop),
    p_schott = safe_pvalue(res_schott),
    p_cheng = safe_pvalue(res_cheng)
  )
}

#------------------------------------------------------------
# Load grid and select exactly one scenario
#------------------------------------------------------------
grid <- fread(opt$grid)
if (!"study" %in% names(grid)) {
  grid[, study := opt$study]
}

if (opt$task_id < 1L || opt$task_id > nrow(grid)) {
  stop(sprintf("task-id %s is outside the grid range 1..%s", opt$task_id, nrow(grid)))
}

row <- grid[opt$task_id]

cat(sprintf("Running task-id %d for scenario %s\n", opt$task_id, row$scenario_id))

#------------------------------------------------------------
# Blocked simulation
#------------------------------------------------------------
n_blocks <- opt$R / opt$R_block

scenario_stub <- row$scenario_id

if (row$study == "power") {
  d_tag <- sprintf("d%02d", as.integer(round(100 * row$signal_d)))
  scenario_stub <- paste0(scenario_stub, "_", d_tag)
}

if (row$study == "robustness") {
  eps_tag <- sprintf("eps%02d", as.integer(round(100 * row$eps)))
  scenario_stub <- paste0(scenario_stub, "_", eps_tag)
}

for (block_id in seq_len(n_blocks)) {
  start_rep <- (block_id - 1L) * opt$R_block + 1L
  end_rep <- block_id * opt$R_block
  
  block_file <- file.path(
    opt$outdir,
    sprintf("%s_block%02d.rds", scenario_stub, block_id)
  )
  
  if (file.exists(block_file)) {
    cat(sprintf("Skipping existing block file: %s\n", block_file))
    next
  }
  
  block_seeds <- opt$seed + seq.int(from = start_rep, to = end_rep) + opt$task_id * 1000000L
  
  if (opt$mc_cores > 1L) {
    block_res <- mclapply(
      seq_along(block_seeds),
      function(i) {
        ans <- run_one_rep(row, block_seeds[i], opt)
        ans$rep_id <- start_rep + i - 1L
        ans
      },
      mc.cores = opt$mc_cores
    )
  } else {
    block_res <- lapply(
      seq_along(block_seeds),
      function(i) {
        ans <- run_one_rep(row, block_seeds[i], opt)
        ans$rep_id <- start_rep + i - 1L
        ans
      }
    )
  }
  
  block_dt <- rbindlist(block_res, fill = TRUE)
  
  save_obj <- list(
    meta = as.list(row),
    study = row$study,
    scenario_id = row$scenario_id,
    block_id = block_id,
    start_rep = start_rep,
    end_rep = end_rep,
    R_block = opt$R_block,
    B = opt$B,
    alpha = opt$alpha,
    # eta = if (row$study == "robustness") max(0.55, min(0.95, 1 - 1.5 * as.numeric(row$eps))) else opt$eta,
    # eta_rule = if (row$study == "robustness") "adaptive: eta = 1 - 1.5 * eps" else "fixed",
    eta = if (row$study == "robustness") max(0.55, min(0.95, 1 - 2.0 * as.numeric(row$eps))) else opt$eta,
    eta_rule = if (row$study == "robustness") "adaptive: eta = 1 - 2.0 * eps" else "fixed",
    delta_rob = opt$delta_rob,
    contamination_note = if (row$study == "robustness") "H0 true for clean data; structured outliers injected into selected groups" else NA_character_,
    results = block_dt
  )
  
  saveRDS(save_obj, block_file)
  cat(sprintf("Saved block %02d: %s\n", block_id, block_file))
}