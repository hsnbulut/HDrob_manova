suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--indir", type = "character", dest = "indir"),
  make_option("--outfile", type = "character", dest = "outfile"),
  make_option("--alpha", type = "double", default = 0.05, dest = "alpha")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$indir) || is.null(opt$outfile)) {
  stop("Arguments --indir and --outfile are required.")
}

files <- list.files(opt$indir, pattern = "\\.rds$", full.names = TRUE)

if (length(files) == 0L) {
  stop("No .rds block files found in ", opt$indir)
}

read_one <- function(f) {
  obj <- readRDS(f)
  dt <- as.data.table(obj$results)
  
  meta <- as.data.table(obj$meta)
  for (nm in names(meta)) {
    dt[[nm]] <- meta[[nm]][1]
  }
  
  dt[, block_id := obj$block_id]
  dt[, start_rep := obj$start_rep]
  dt[, end_rep := obj$end_rep]
  dt
}

all_res <- rbindlist(lapply(files, read_one), fill = TRUE)

summarise_one <- function(pvals, alpha = 0.05) {
  ok <- is.finite(pvals)
  if (!any(ok)) {
    return(list(rate = NA_real_, mcse = NA_real_, n_eff = 0L))
  }
  rate <- mean(pvals[ok] < alpha)
  mcse <- sqrt(rate * (1 - rate) / sum(ok))
  list(rate = rate, mcse = mcse, n_eff = sum(ok))
}

group_cols <- c("study", "scenario_id", "G", "nk", "p", "cov_type", "design", "signal_d", "eps")
group_cols <- group_cols[group_cols %in% names(all_res)]

summary_dt <- all_res[, {
  s1 <- summarise_one(p_proposed, alpha = opt$alpha)
  s2 <- summarise_one(p_schott, alpha = opt$alpha)
  s3 <- summarise_one(p_cheng, alpha = opt$alpha)
  
  .(
    proposed = s1$rate,
    proposed_mcse = s1$mcse,
    schott = s2$rate,
    schott_mcse = s2$mcse,
    cheng = s3$rate,
    cheng_mcse = s3$mcse,
    n_eff_proposed = s1$n_eff,
    n_eff_schott = s2$n_eff,
    n_eff_cheng = s3$n_eff,
    n_blocks = uniqueN(block_id),
    n_reps = .N
  )
}, by = group_cols]

dir.create(dirname(opt$outfile), recursive = TRUE, showWarnings = FALSE)
fwrite(summary_dt, opt$outfile)

cat("Combined scenarios:", nrow(summary_dt), "\n")
cat("Saved to:", opt$outfile, "\n")