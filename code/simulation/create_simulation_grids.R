suppressPackageStartupMessages({
  library(data.table)
})

source("code/utilities/parse_nk.R")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

make_base_scenarios <- function() {
  out <- list()
  s_id <- 1L
  
  add_scenario <- function(G, nk, p_vals, design_label, out, s_id) {
    cov_types <- c("identity", "ar1")
    
    for (p in p_vals) {
      for (cov_type in cov_types) {
        out[[length(out) + 1L]] <- data.table(
          scenario_id = sprintf("S%02d", s_id),
          G = G,
          nk = paste(nk, collapse = ","),
          p = as.integer(p),
          cov_type = cov_type,
          design = design_label
        )
        s_id <- s_id + 1L
      }
    }
    
    list(out = out, s_id = s_id)
  }
  
  # Balanced
  tmp <- add_scenario(2L, c(20L, 20L), c(30L, 50L), "balanced", out, s_id)
  out <- tmp$out; s_id <- tmp$s_id
  
  tmp <- add_scenario(2L, c(30L, 30L), c(50L, 100L), "balanced", out, s_id)
  out <- tmp$out; s_id <- tmp$s_id
  
  tmp <- add_scenario(2L, c(50L, 50L), c(50L, 100L), "balanced", out, s_id)
  out <- tmp$out; s_id <- tmp$s_id
  
  tmp <- add_scenario(3L, c(20L, 20L, 20L), c(30L, 50L), "balanced", out, s_id)
  out <- tmp$out; s_id <- tmp$s_id
  
  tmp <- add_scenario(3L, c(30L, 30L, 30L), c(50L, 100L), "balanced", out, s_id)
  out <- tmp$out; s_id <- tmp$s_id
  
  tmp <- add_scenario(3L, c(50L, 50L, 50L), c(50L, 100L), "balanced", out, s_id)
  out <- tmp$out; s_id <- tmp$s_id
  
  # Unbalanced
  tmp <- add_scenario(2L, c(30L, 20L), c(30L, 50L), "unbalanced", out, s_id)
  out <- tmp$out; s_id <- tmp$s_id
  
  tmp <- add_scenario(2L, c(50L, 30L), c(50L, 100L), "unbalanced", out, s_id)
  out <- tmp$out; s_id <- tmp$s_id
  
  tmp <- add_scenario(3L, c(30L, 20L, 20L), c(30L, 50L), "unbalanced", out, s_id)
  out <- tmp$out; s_id <- tmp$s_id
  
  tmp <- add_scenario(3L, c(50L, 30L, 20L), c(50L, 100L), "unbalanced", out, s_id)
  out <- tmp$out; s_id <- tmp$s_id
  
  rbindlist(out)
}

make_type1_grid <- function(base_grid) {
  grid <- copy(base_grid)
  grid[, study := "type1"]
  grid[, signal_d := NA_real_]
  grid[, eps := NA_real_]
  grid[]
}

make_power_grid <- function(base_grid, d_vals = c(0.2, 0.4, 0.6)) {
  rbindlist(lapply(d_vals, function(d) {
    g <- copy(base_grid)
    g[, study := "power"]
    g[, signal_d := d]
    g[, eps := NA_real_]
    g
  }))
}

make_robustness_grid <- function(base_grid, eps_vals = c(0.10, 0.20)) {
  rbindlist(lapply(eps_vals, function(eps) {
    g <- copy(base_grid)
    g[, study := "robustness"]
    g[, signal_d := NA_real_]
    g[, eps := eps]
    g
  }))
}

base_grid <- make_base_scenarios()

type1_grid <- make_type1_grid(base_grid)
power_grid <- make_power_grid(base_grid, d_vals = c(0.2, 0.4, 0.6))
robustness_grid <- make_robustness_grid(base_grid, eps_vals = c(0.10, 0.20))

fwrite(type1_grid, "results/tables/grid_type1.csv")
fwrite(power_grid, "results/tables/grid_power.csv")
fwrite(robustness_grid, "results/tables/grid_robustness.csv")

cat("Type-I grid rows:", nrow(type1_grid), "\n")
cat("Power grid rows:", nrow(power_grid), "\n")
cat("Robustness grid rows:", nrow(robustness_grid), "\n")