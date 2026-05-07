## =========================================================
## TYPE-I RESULTS: combine + tables + p-value plots
## For HDrob_manova project
## =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

## ---------------------------------------------------------
## User settings
## ---------------------------------------------------------
alpha_nominal <- 0.05
R_total <- 1000
blocks_per_scenario <- 10

raw_dir   <- file.path("results", "raw_chunks", "type1")
table_dir <- file.path("results", "tables")
fig_dir   <- file.path("results", "figures", "type1")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,   recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------
## Helper: Monte Carlo standard error
## ---------------------------------------------------------
mc_se <- function(phat, R) {
  sqrt(phat * (1 - phat) / R)
}

## ---------------------------------------------------------
## Helper: read one block safely
## ---------------------------------------------------------
read_one_block <- function(f) {
  obj <- readRDS(f)
  
  required_top <- c(
    "meta", "study", "scenario_id", "block_id",
    "start_rep", "end_rep", "results"
  )
  miss_top <- setdiff(required_top, names(obj))
  if (length(miss_top) > 0) {
    stop(sprintf(
      "File %s is missing top-level fields: %s",
      basename(f), paste(miss_top, collapse = ", ")
    ))
  }
  
  dt <- as.data.table(obj$results)
  
  ## expected columns in results
  ## from your structure: 100 obs x 4 vars
  ## likely rep_id, p_proposed, p_schott, p_cheng
  expected_cols <- c("rep_id", "p_proposed", "p_schott", "p_cheng")
  miss_cols <- setdiff(expected_cols, names(dt))
  if (length(miss_cols) > 0) {
    stop(sprintf(
      "File %s is missing result columns: %s\nAvailable columns: %s",
      basename(f),
      paste(miss_cols, collapse = ", "),
      paste(names(dt), collapse = ", ")
    ))
  }
  
  ## attach metadata
  dt[, scenario_id := obj$scenario_id]
  dt[, block_id    := obj$block_id]
  dt[, start_rep   := obj$start_rep]
  dt[, end_rep     := obj$end_rep]
  dt[, source_file := basename(f)]
  
  ## flatten meta
  meta <- obj$meta
  if (is.list(meta) && length(meta) > 0) {
    for (nm in names(meta)) {
      val <- meta[[nm]]
      if (length(val) == 1L) {
        dt[, (nm) := val]
      }
    }
  }
  
  dt[]
}

## ---------------------------------------------------------
## 1) Read all block files
## ---------------------------------------------------------
files <- list.files(raw_dir, pattern = "\\.rds$", full.names = TRUE)

if (length(files) == 0L) {
  stop("No .rds files found in: ", raw_dir)
}

cat("Number of block files found:", length(files), "\n")

all_dt <- rbindlist(lapply(files, read_one_block), fill = TRUE)

cat("Total rows after combine:", nrow(all_dt), "\n")
cat("Columns in combined data:\n")
print(names(all_dt))

## ---------------------------------------------------------
## 2) Basic integrity checks
## ---------------------------------------------------------
## Expected:
## 40 scenarios x 10 blocks x 100 reps = 40000 rows
## but we do not hard-stop if different; we report.

## block count per scenario
block_check <- all_dt[
  , .(n_blocks = uniqueN(block_id)),
  by = scenario_id
][order(scenario_id)]

## rep count per scenario
rep_check <- all_dt[
  , .(
    n_reps = .N,
    min_rep = min(rep_id, na.rm = TRUE),
    max_rep = max(rep_id, na.rm = TRUE)
  ),
  by = scenario_id
][order(scenario_id)]

## merged completeness check
scenario_check <- merge(block_check, rep_check, by = "scenario_id", all = TRUE)
scenario_check[, expected_blocks := blocks_per_scenario]
scenario_check[, expected_reps   := R_total]
scenario_check[, block_complete  := (n_blocks == expected_blocks)]
scenario_check[, rep_complete    := (n_reps   == expected_reps)]

fwrite(
  scenario_check,
  file.path(table_dir, "type1_scenario_completeness_check.csv")
)

cat("\nScenario completeness summary:\n")
print(scenario_check)

## ---------------------------------------------------------
## 3) Rejection indicators
## ---------------------------------------------------------
all_dt[, reject_proposed := as.integer(p_proposed < alpha_nominal)]
all_dt[, reject_schott   := as.integer(p_schott   < alpha_nominal)]
all_dt[, reject_cheng    := as.integer(p_cheng    < alpha_nominal)]

## ---------------------------------------------------------
## 4) Scenario-level empirical size table
## ---------------------------------------------------------
## Keep key scenario descriptors if present
candidate_meta_cols <- c("scenario_id", "G", "nk", "p", "cov_type", "study")
group_cols <- candidate_meta_cols[candidate_meta_cols %in% names(all_dt)]

type1_table <- all_dt[
  ,
  .(
    R_effective_proposed = sum(!is.na(p_proposed)),
    R_effective_schott   = sum(!is.na(p_schott)),
    R_effective_cheng    = sum(!is.na(p_cheng)),
    
    proposed_size = mean(reject_proposed, na.rm = TRUE),
    schott_size   = mean(reject_schott,   na.rm = TRUE),
    cheng_size    = mean(reject_cheng,    na.rm = TRUE)
  ),
  by = group_cols
][order(scenario_id)]

type1_table[, proposed_mcse := mc_se(proposed_size, R_effective_proposed)]
type1_table[, schott_mcse   := mc_se(schott_size,   R_effective_schott)]
type1_table[, cheng_mcse    := mc_se(cheng_size,    R_effective_cheng)]

## absolute deviation from nominal alpha
type1_table[, proposed_absdev := abs(proposed_size - alpha_nominal)]
type1_table[, schott_absdev   := abs(schott_size   - alpha_nominal)]
type1_table[, cheng_absdev    := abs(cheng_size    - alpha_nominal)]

fwrite(
  type1_table,
  file.path(table_dir, "type1_empirical_size_table.csv")
)

cat("\nScenario-level type1 table:\n")
print(type1_table)

## ---------------------------------------------------------
## 5) A LaTeX-ready compact table (rounded)
## ---------------------------------------------------------
latex_tab <- copy(type1_table)

num_cols <- c(
  "proposed_size", "proposed_mcse",
  "schott_size",   "schott_mcse",
  "cheng_size",    "cheng_mcse",
  "proposed_absdev", "schott_absdev", "cheng_absdev"
)
for (cc in intersect(num_cols, names(latex_tab))) {
  latex_tab[, (cc) := round(get(cc), 4)]
}

fwrite(
  latex_tab,
  file.path(table_dir, "type1_empirical_size_table_rounded.csv")
)

## ---------------------------------------------------------
## 6) Overall summary across all scenarios
## ---------------------------------------------------------
overall_summary <- data.table(
  method = c("Proposed", "Schott", "Cheng"),
  mean_size = c(
    mean(type1_table$proposed_size, na.rm = TRUE),
    mean(type1_table$schott_size,   na.rm = TRUE),
    mean(type1_table$cheng_size,    na.rm = TRUE)
  ),
  mean_abs_dev = c(
    mean(type1_table$proposed_absdev, na.rm = TRUE),
    mean(type1_table$schott_absdev,   na.rm = TRUE),
    mean(type1_table$cheng_absdev,    na.rm = TRUE)
  ),
  max_abs_dev = c(
    max(type1_table$proposed_absdev, na.rm = TRUE),
    max(type1_table$schott_absdev,   na.rm = TRUE),
    max(type1_table$cheng_absdev,    na.rm = TRUE)
  )
)

fwrite(
  overall_summary,
  file.path(table_dir, "type1_overall_summary.csv")
)

cat("\nOverall summary:\n")
print(overall_summary)

## ---------------------------------------------------------
## 7) Prepare data for p-value plots
## ---------------------------------------------------------
## Davidson-MacKinnon style:
## Under correct calibration, ordered p-values should track Uniform(0,1).

make_pplot_dt <- function(pvals, scenario_id, method_name) {
  pvals <- pvals[is.finite(pvals)]
  n <- length(pvals)
  if (n == 0L) return(NULL)
  
  p_sorted <- sort(pvals)
  ## expected order stats approximation
  u_exp <- ((1:n) - 0.5) / n
  
  data.table(
    scenario_id = scenario_id,
    method = method_name,
    order_id = 1:n,
    expected = u_exp,
    observed = p_sorted,
    discrepancy = p_sorted - u_exp
  )
}

pplot_dt <- rbindlist(
  lapply(split(all_dt, by = "scenario_id"), function(dd) {
    sid <- dd$scenario_id[1]
    rbindlist(list(
      make_pplot_dt(dd$p_proposed, sid, "Proposed"),
      make_pplot_dt(dd$p_schott,   sid, "Schott"),
      make_pplot_dt(dd$p_cheng,    sid, "Cheng")
    ), fill = TRUE)
  }),
  fill = TRUE
)

## Save combined p-value data if needed
fwrite(
  pplot_dt,
  file.path(table_dir, "type1_pvalue_plot_data.csv")
)

## ---------------------------------------------------------
## 8) Plot function: P-value plot
## ---------------------------------------------------------
plot_pvalue <- function(dt, sid, outdir) {
  dsub <- dt[scenario_id == sid]
  
  g <- ggplot(dsub, aes(x = expected, y = observed, color = method)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, linewidth = 0.6) +
    geom_line(linewidth = 0.7) +
    coord_cartesian(xlim = c(0, 0.20), ylim = c(0, 0.20)) +
    labs(
      title = paste("P-value plot -", sid),
      x = "Expected uniform order statistics",
      y = "Observed ordered p-values"
    ) +
    theme_bw(base_size = 12)
  
  ggsave(
    filename = file.path(outdir, paste0("pvalue_plot_", sid, ".png")),
    plot = g,
    width = 7, height = 5, dpi = 300
  )
}
## ---------------------------------------------------------
## 9) Plot function: P-value discrepancy plot
## ---------------------------------------------------------
plot_discrepancy <- function(dt, sid, outdir) {
  dsub <- dt[scenario_id == sid]
  
  g <- ggplot(dsub, aes(x = expected, y = discrepancy, color = method)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.6) +
    geom_line(linewidth = 0.7) +
    coord_cartesian(xlim = c(0, 0.20)) +
    labs(
      title = paste("P-value discrepancy plot -", sid),
      x = "Expected uniform order statistics",
      y = "Observed - expected"
    ) +
    theme_bw(base_size = 12)
  
  ggsave(
    filename = file.path(outdir, paste0("pvalue_discrepancy_", sid, ".png")),
    plot = g,
    width = 7, height = 5, dpi = 300
  )
}
## ---------------------------------------------------------
## 10) Create plots for all scenarios
## ---------------------------------------------------------
scenario_ids <- sort(unique(all_dt$scenario_id))

for (sid in scenario_ids) {
  plot_pvalue(pplot_dt, sid, fig_dir)
  plot_discrepancy(pplot_dt, sid, fig_dir)
}

cat("\nAll scenario p-value plots saved to:\n", fig_dir, "\n")

## ---------------------------------------------------------
## 11) Optional: combined faceted plots for quick inspection
## ---------------------------------------------------------
g_all_p <- ggplot(pplot_dt, aes(x = expected, y = observed, color = method)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, linewidth = 0.4) +
  geom_line(linewidth = 0.5) +
  coord_cartesian(xlim = c(0, 0.20), ylim = c(0, 0.20)) +
  facet_wrap(~ scenario_id, ncol = 4) +
  labs(
    title = "P-value plots for all Type-I scenarios",
    x = "Expected uniform order statistics",
    y = "Observed ordered p-values"
  ) +
  theme_bw(base_size = 10)

ggsave(
  filename = file.path(fig_dir, "pvalue_plots_all_scenarios.png"),
  plot = g_all_p,
  width = 16, height = 20, dpi = 300
)

g_all_d <- ggplot(pplot_dt, aes(x = expected, y = discrepancy, color = method)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ scenario_id, ncol = 4) +
  labs(
    title = "P-value discrepancy plots for all Type-I scenarios",
    x = "Expected uniform order statistics",
    y = "Observed - expected"
  ) +
  theme_bw(base_size = 10)

ggsave(
  filename = file.path(fig_dir, "pvalue_discrepancy_all_scenarios.png"),
  plot = g_all_d,
  width = 16, height = 20, dpi = 300
)

## ---------------------------------------------------------
## 12) Optional: method comparison plot of empirical sizes
## ---------------------------------------------------------
size_long <- rbindlist(list(
  type1_table[, .(scenario_id, method = "Proposed", size = proposed_size)],
  type1_table[, .(scenario_id, method = "Schott",   size = schott_size)],
  type1_table[, .(scenario_id, method = "Cheng",    size = cheng_size)]
))

g_size <- ggplot(size_long, aes(x = scenario_id, y = size, color = method, group = method)) +
  geom_hline(yintercept = alpha_nominal, linetype = 2, linewidth = 0.6) +
  geom_point(size = 1.8) +
  geom_line(linewidth = 0.7) +
  labs(
    title = "Empirical Type-I error rates across scenarios",
    x = "Scenario",
    y = "Empirical size"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave(
  filename = file.path(fig_dir, "type1_empirical_sizes.png"),
  plot = g_size,
  width = 12, height = 6, dpi = 300
)

## ---------------------------------------------------------
## 13) Console summary
## ---------------------------------------------------------
cat("\n====================================================\n")
cat("TYPE-I combine completed.\n")
cat("Raw blocks read: ", length(files), "\n", sep = "")
cat("Scenarios found: ", length(unique(all_dt$scenario_id)), "\n", sep = "")
cat("Figures saved in: ", fig_dir, "\n", sep = "")
cat("Tables saved in:  ", table_dir, "\n", sep = "")
cat("====================================================\n")