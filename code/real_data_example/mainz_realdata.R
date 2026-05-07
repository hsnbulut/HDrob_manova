###############################################################################
# File: code/real_data_example/mainz_realdata.R
# Purpose:
#   Final real-data application for the proposed MRCD-based robust HD-MANOVA test.
#
# Dataset:
#   breastCancerMAINZ gene expression data.
#
# Design:
#   1) Load the breastCancerMAINZ dataset.
#   2) Select the largest tumor-grade group.
#   3) Create three pseudo-groups from this single real group:
#        G1 = 40, G2 = 40, G3 = 40.
#   4) Select p = 300 high-MAD genes within the selected group.
#   5) Standardize variables using the selected group only.
#   6) Run the three methods on clean pseudo-groups:
#        - Schott
#        - Cheng
#        - proposed MRCD
#   7) Apply structured block contamination:
#        - last 12 observations of G1: +5 shift on first 120 selected genes
#        - last 12 observations of G2: -5 shift on first 120 selected genes
#        - G3 unchanged
#   8) Re-run the three methods.
#
# Outputs:
#   results/tables/mainz_realdata_final_results.csv
#   results/tables/mainz_realdata_group_shift_summary.csv
#   results/tables/mainz_realdata_contamination_design.csv
#   results/raw_chunks/realdata/mainz_realdata_final_objects.rds
#   results/figures/mainz_realdata_final_pvalues.png
###############################################################################

rm(list = ls())

# -------------------------------------------------------------------------
# 0) User settings
# -------------------------------------------------------------------------
set.seed(2026)

n_per_group <- 40
p_use <- 300

eta_mrcd <- 0.75
alpha_level <- 0.05

B_boot <- 499
B_perm <- 499

contam_m <- 12
contam_block_size <- 120
contam_shift <- 5

split_seed <- 2026

# -------------------------------------------------------------------------
# 1) Packages
# -------------------------------------------------------------------------
cran_pkgs <- c(
  "here", "ggplot2", "dplyr", "rrcov", "Gmedian", "extraDistr"
)

bioc_pkgs <- c(
  "Biobase", "breastCancerMAINZ"
)

to_install_cran <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install_cran) > 0) {
  install.packages(to_install_cran)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

to_install_bioc <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install_bioc) > 0) {
  BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)
}

invisible(lapply(c(cran_pkgs, bioc_pkgs), library, character.only = TRUE))

# -------------------------------------------------------------------------
# 2) Output folders
# -------------------------------------------------------------------------
dir.create(here::here("results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("results", "raw_chunks", "realdata"), recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# 3) Source methods
# -------------------------------------------------------------------------
source(here::here("code", "competitors", "hdmanova_schott.R"))
source(here::here("code", "competitors", "hdmanova_cheng.R"))
source(here::here("code", "core", "hdmanova_mrcd.R"))

# -------------------------------------------------------------------------
# 4) Load breastCancerMAINZ ExpressionSet
# -------------------------------------------------------------------------
pkg_data <- data(package = "breastCancerMAINZ")$results[, "Item"]

if ("mainz" %in% pkg_data) {
  data_obj_name <- "mainz"
} else {
  mainz_like <- grep("mainz", pkg_data, ignore.case = TRUE, value = TRUE)
  if (length(mainz_like) == 0) {
    stop("Could not find a MAINZ data object in breastCancerMAINZ package.")
  }
  data_obj_name <- mainz_like[1]
}

env_data <- new.env()
data(list = data_obj_name, package = "breastCancerMAINZ", envir = env_data)
eset <- env_data[[data_obj_name]]

if (!inherits(eset, "ExpressionSet")) {
  stop("Loaded object is not an ExpressionSet. Please inspect breastCancerMAINZ data objects.")
}

expr_mat <- Biobase::exprs(eset)  # genes x samples
pd <- Biobase::pData(eset)

X_full <- t(expr_mat)             # samples x genes

keep_cols <- apply(X_full, 2, function(z) all(is.finite(z)))
X_full <- X_full[, keep_cols, drop = FALSE]

cat("\nLoaded breastCancerMAINZ object:", data_obj_name, "\n")
cat("Full data: n =", nrow(X_full), ", p =", ncol(X_full), "\n")

# -------------------------------------------------------------------------
# 5) Detect tumor-grade column
# -------------------------------------------------------------------------
detect_grade_column <- function(pd) {
  nms <- names(pd)
  
  preferred <- grep("grade|grading|hist|tumou?r", nms, ignore.case = TRUE, value = TRUE)
  
  score_col <- function(v) {
    vv <- as.character(v)
    vv <- vv[!is.na(vv) & vv != ""]
    if (length(vv) == 0) return(-Inf)
    
    tb <- table(vv)
    u <- length(tb)
    maxn <- max(tb)
    
    if (u >= 3 && u <= 5) return(maxn)
    return(-Inf)
  }
  
  candidates <- unique(c(preferred, nms))
  scores <- sapply(candidates, function(nm) score_col(pd[[nm]]))
  
  candidates <- candidates[is.finite(scores)]
  scores <- scores[is.finite(scores)]
  
  if (length(candidates) == 0) {
    cat("\nAvailable phenotype columns:\n")
    print(names(pd))
    stop("Could not detect a plausible tumor-grade column.")
  }
  
  candidates[which.max(scores)]
}

grade_col <- detect_grade_column(pd)
grade_raw <- as.character(pd[[grade_col]])

valid_grade <- !is.na(grade_raw) & grade_raw != ""
grade <- factor(grade_raw[valid_grade])
X_full <- X_full[valid_grade, , drop = FALSE]

cat("\nDetected tumor-grade column:", grade_col, "\n")
cat("Grade table:\n")
print(table(grade))

# -------------------------------------------------------------------------
# 6) Select largest tumor-grade group
# -------------------------------------------------------------------------
grade_tab <- table(grade)
base_grade <- names(which.max(grade_tab))
base_idx <- which(grade == base_grade)

if (length(base_idx) < 3 * n_per_group) {
  stop(
    "Largest grade group has only ", length(base_idx),
    " observations. Need at least ", 3 * n_per_group, "."
  )
}

X_base_all <- X_full[base_idx, , drop = FALSE]

cat("\nSelected base tumor-grade group:", base_grade, "\n")
cat("Available observations in base group:", nrow(X_base_all), "\n")

# -------------------------------------------------------------------------
# 7) Select high-MAD genes within the selected grade group
# -------------------------------------------------------------------------
gene_mad <- apply(X_base_all, 2, mad)
gene_sd  <- apply(X_base_all, 2, sd)

keep_gene <- is.finite(gene_mad) & is.finite(gene_sd) & (gene_sd > 1e-8)

X_base_all <- X_base_all[, keep_gene, drop = FALSE]
gene_mad <- gene_mad[keep_gene]

p_use <- min(p_use, ncol(X_base_all))
ord <- order(gene_mad, decreasing = TRUE)
sel_cols <- ord[1:p_use]

X_base_all <- X_base_all[, sel_cols, drop = FALSE]

# Standardize using the selected base group only
mu_ref <- colMeans(X_base_all)
sd_ref <- apply(X_base_all, 2, sd)
sd_ref[!is.finite(sd_ref) | sd_ref < 1e-8] <- 1

X_base_all <- sweep(X_base_all, 2, mu_ref, "-")
X_base_all <- sweep(X_base_all, 2, sd_ref, "/")

cat("\nSelected analysis matrix from one real tumor-grade group:\n")
cat("n =", nrow(X_base_all), ", p =", ncol(X_base_all), "\n")

# -------------------------------------------------------------------------
# 8) Helper functions
# -------------------------------------------------------------------------
extract_p <- function(obj) {
  if (is.null(obj) || !is.list(obj)) return(NA_real_)
  
  for (nm in c("p.value", "p_value", "pval", "p")) {
    if (nm %in% names(obj)) return(as.numeric(obj[[nm]]))
  }
  
  NA_real_
}

extract_stat <- function(obj) {
  if (is.null(obj) || !is.list(obj)) return(NA_real_)
  
  for (nm in c("statistic", "test_stat", "lambda", "T", "LR")) {
    if (nm %in% names(obj)) return(as.numeric(obj[[nm]]))
  }
  
  NA_real_
}

run_safe <- function(expr) {
  out <- try(expr, silent = TRUE)
  
  if (inherits(out, "try-error")) {
    list(ok = FALSE, result = NULL, error = as.character(out))
  } else {
    list(ok = TRUE, result = out, error = NA_character_)
  }
}

run_all_methods <- function(X, group, B_boot = 1999, B_perm = 499, alpha = 0.05) {
  schott <- run_safe(
    hdmanova_schott(X = X, group = group)
  )
  
  cheng <- run_safe(
    hdmanova_cheng(
      X = X,
      group = group,
      B = B_boot,
      alpha = alpha,
      standardized = TRUE
    )
  )
  
  mrcd <- run_safe(
    hdmanova_mrcd(
      X = X,
      group = group,
      eta = eta_mrcd,
      B = B_perm,
      alpha = alpha,
      seed = 20260427
    )
  )
  
  data.frame(
    method = c("Schott", "Cheng", "Proposed MRCD"),
    ok = c(schott$ok, cheng$ok, mrcd$ok),
    statistic = c(
      extract_stat(schott$result),
      extract_stat(cheng$result),
      extract_stat(mrcd$result)
    ),
    p_value = c(
      extract_p(schott$result),
      extract_p(cheng$result),
      extract_p(mrcd$result)
    ),
    error = c(
      schott$error,
      cheng$error,
      mrcd$error
    ),
    stringsAsFactors = FALSE
  )
}

make_pseudogroups <- function(X_base_all, seed, n_per_group = 40) {
  set.seed(seed)
  
  idx <- sample(seq_len(nrow(X_base_all)), size = 3 * n_per_group, replace = FALSE)
  
  g1_idx <- idx[1:n_per_group]
  g2_idx <- idx[(n_per_group + 1):(2 * n_per_group)]
  g3_idx <- idx[(2 * n_per_group + 1):(3 * n_per_group)]
  
  X <- rbind(
    X_base_all[g1_idx, , drop = FALSE],
    X_base_all[g2_idx, , drop = FALSE],
    X_base_all[g3_idx, , drop = FALSE]
  )
  
  group <- factor(c(
    rep("G1", n_per_group),
    rep("G2", n_per_group),
    rep("G3", n_per_group)
  ), levels = c("G1", "G2", "G3"))
  
  list(
    X = X,
    group = group,
    g1_idx = g1_idx,
    g2_idx = g2_idx,
    g3_idx = g3_idx,
    seed = seed
  )
}

apply_structured_contamination <- function(X, group,
                                           m = 12,
                                           block_size = 120,
                                           shift = 5) {
  Xc <- X
  
  idx_g1 <- which(group == "G1")
  idx_g2 <- which(group == "G2")
  
  rows_g1 <- tail(idx_g1, m)
  rows_g2 <- tail(idx_g2, m)
  
  cols_blk <- seq_len(min(block_size, ncol(Xc)))
  
  Xc[rows_g1, cols_blk] <- Xc[rows_g1, cols_blk] + shift
  Xc[rows_g2, cols_blk] <- Xc[rows_g2, cols_blk] - shift
  
  list(
    X = Xc,
    rows_g1 = rows_g1,
    rows_g2 = rows_g2,
    cols_blk = cols_blk,
    m = m,
    block_size = length(cols_blk),
    shift = shift
  )
}

pairwise_group_shift_summary <- function(X, group) {
  levs <- levels(group)
  
  gm <- rowsum(X, group) / as.vector(table(group))
  gm <- as.matrix(gm)
  
  out <- list()
  idx <- 1L
  
  for (i in 1:(length(levs) - 1)) {
    for (j in (i + 1):length(levs)) {
      g1 <- levs[i]
      g2 <- levs[j]
      
      m1 <- gm[g1, ]
      m2 <- gm[g2, ]
      
      pooled_sd <- apply(
        rbind(X[group == g1, , drop = FALSE], X[group == g2, , drop = FALSE]),
        2,
        sd
      )
      
      pooled_sd[!is.finite(pooled_sd) | pooled_sd < 1e-8] <- 1e-8
      
      d_raw <- m1 - m2
      d_std <- d_raw / pooled_sd
      
      out[[idx]] <- data.frame(
        group1 = g1,
        group2 = g2,
        corr_group_means = suppressWarnings(cor(m1, m2)),
        l2_distance = sqrt(sum(d_raw^2)),
        median_abs_mean_diff = median(abs(d_raw)),
        median_abs_std_mean_diff = median(abs(d_std)),
        max_abs_std_mean_diff = max(abs(d_std)),
        stringsAsFactors = FALSE
      )
      
      idx <- idx + 1L
    }
  }
  
  dplyr::bind_rows(out)
}

# -------------------------------------------------------------------------
# 9) Create clean pseudo-groups
# -------------------------------------------------------------------------
split_obj <- make_pseudogroups(
  X_base_all = X_base_all,
  seed = split_seed,
  n_per_group = n_per_group
)

X_before <- split_obj$X
group_before <- split_obj$group

cat("\nPseudo-groups before contamination:\n")
print(table(group_before))

# -------------------------------------------------------------------------
# 10) Run methods before contamination
# -------------------------------------------------------------------------
cat("\nRunning methods before contamination...\n")

before_results <- run_all_methods(
  X = X_before,
  group = group_before,
  B_boot = B_boot,
  B_perm = B_perm,
  alpha = alpha_level
) %>%
  dplyr::mutate(scenario = "Before contamination")

print(before_results)

# -------------------------------------------------------------------------
# 11) Apply structured block contamination
# -------------------------------------------------------------------------
contam_obj <- apply_structured_contamination(
  X = X_before,
  group = group_before,
  m = contam_m,
  block_size = contam_block_size,
  shift = contam_shift
)

X_after <- contam_obj$X
group_after <- group_before

cat("\nStructured contamination applied:\n")
cat("G1 shifted rows:", paste(contam_obj$rows_g1, collapse = ", "), "\n")
cat("G2 shifted rows:", paste(contam_obj$rows_g2, collapse = ", "), "\n")
cat("Number of shifted variables:", contam_obj$block_size, "\n")
cat("Shift magnitude:", contam_obj$shift, "\n")

# -------------------------------------------------------------------------
# 12) Run methods after contamination
# -------------------------------------------------------------------------
cat("\nRunning methods after contamination...\n")

after_results <- run_all_methods(
  X = X_after,
  group = group_after,
  B_boot = B_boot,
  B_perm = B_perm,
  alpha = alpha_level
) %>%
  dplyr::mutate(scenario = "After contamination")

print(after_results)

# -------------------------------------------------------------------------
# 13) Summaries
# -------------------------------------------------------------------------
results_final <- dplyr::bind_rows(before_results, after_results) %>%
  dplyr::select(scenario, method, ok, statistic, p_value, error)

shift_summary <- dplyr::bind_rows(
  pairwise_group_shift_summary(X_before, group_before) %>%
    dplyr::mutate(scenario = "Before contamination"),
  pairwise_group_shift_summary(X_after, group_after) %>%
    dplyr::mutate(scenario = "After contamination")
) %>%
  dplyr::select(
    scenario, group1, group2,
    corr_group_means,
    l2_distance,
    median_abs_mean_diff,
    median_abs_std_mean_diff,
    max_abs_std_mean_diff
  )

contamination_design <- data.frame(
  dataset = "breastCancerMAINZ",
  data_object = data_obj_name,
  base_grade_group = base_grade,
  n_per_group = n_per_group,
  p_selected = p_use,
  pseudo_group_seed = split_seed,
  contaminated_groups = "G1 and G2",
  unchanged_group = "G3",
  contaminated_observations_per_group = contam_m,
  contaminated_variables = contam_block_size,
  shift_G1 = contam_shift,
  shift_G2 = -contam_shift,
  eta_mrcd = eta_mrcd,
  B_boot_cheng = B_boot,
  B_perm_mrcd = B_perm,
  stringsAsFactors = FALSE
)

# -------------------------------------------------------------------------
# 14) Save outputs
# -------------------------------------------------------------------------
write.csv(
  results_final,
  here::here("results", "tables", "mainz_realdata_final_results.csv"),
  row.names = FALSE
)

write.csv(
  shift_summary,
  here::here("results", "tables", "mainz_realdata_group_shift_summary.csv"),
  row.names = FALSE
)

write.csv(
  contamination_design,
  here::here("results", "tables", "mainz_realdata_contamination_design.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    dataset = "breastCancerMAINZ",
    data_object = data_obj_name,
    phenotype_data = pd,
    grade_column = grade_col,
    base_grade = base_grade,
    selected_gene_indices = sel_cols,
    split_seed = split_seed,
    n_per_group = n_per_group,
    p_use = p_use,
    X_before = X_before,
    X_after = X_after,
    group_before = group_before,
    group_after = group_after,
    contamination = contam_obj,
    contamination_design = contamination_design,
    results = results_final,
    shift_summary = shift_summary
  ),
  here::here("results", "raw_chunks", "realdata", "mainz_realdata_final_objects.rds")
)

# -------------------------------------------------------------------------
# 15) Plot p-values
# -------------------------------------------------------------------------
plot_df <- results_final %>%
  dplyr::mutate(
    scenario = factor(
      scenario,
      levels = c("Before contamination", "After contamination")
    ),
    method = factor(
      method,
      levels = c("Proposed MRCD", "Schott", "Cheng")
    )
  )

g <- ggplot(plot_df, aes(x = method, y = p_value, fill = scenario)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_hline(yintercept = 0.05, linetype = 2, linewidth = 0.7) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Real-data-based high-dimensional MANOVA example",
    subtitle = "Clean pseudo-groups vs structured block-contaminated pseudo-groups",
    x = NULL,
    y = "p-value"
  )

ggsave(
  filename = here::here("results", "figures", "mainz_realdata_final_pvalues.png"),
  plot = g,
  width = 8,
  height = 5,
  dpi = 300
)

# -------------------------------------------------------------------------
# 16) Console output
# -------------------------------------------------------------------------
cat("\nSaved files:\n")
cat(" - results/tables/mainz_realdata_final_results.csv\n")
cat(" - results/tables/mainz_realdata_group_shift_summary.csv\n")
cat(" - results/tables/mainz_realdata_contamination_design.csv\n")
cat(" - results/raw_chunks/realdata/mainz_realdata_final_objects.rds\n")
cat(" - results/figures/mainz_realdata_final_pvalues.png\n")

cat("\nFinal p-values:\n")
print(results_final)

cat("\nContamination design:\n")
print(contamination_design)

cat("\nDone.\n")