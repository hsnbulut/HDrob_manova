# ============================================================
# create_size_power_plots.R
# Create size-power curves for all scenarios and d levels
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

raw_dir <- "results/raw_chunks/power"
out_dir <- "results/figures/power/size_power"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

alpha_grid <- seq(0, 1, length.out = 201)

d_info <- data.table(
  dtag = c("d20", "d40", "d60"),
  dval = c(0.2, 0.4, 0.6)
)

make_curve <- function(scenario_id, dtag, dval) {
  
  files <- list.files(
    raw_dir,
    pattern = paste0("^", scenario_id, "_", dtag, "_block[0-9]+\\.rds$"),
    full.names = TRUE
  )
  
  if (length(files) != 10) {
    warning(paste0(scenario_id, "_", dtag, ": expected 10 files, found ", length(files)))
  }
  
  A <- rbindlist(lapply(files, function(f) readRDS(f)$results), fill = TRUE)
  
  curve_dt <- rbindlist(list(
    data.table(
      alpha = alpha_grid,
      rejection_rate = sapply(alpha_grid, function(a) mean(A$p_proposed <= a, na.rm = TRUE)),
      method = "Proposed"
    ),
    data.table(
      alpha = alpha_grid,
      rejection_rate = sapply(alpha_grid, function(a) mean(A$p_schott <= a, na.rm = TRUE)),
      method = "Schott (2007)"
    ),
    data.table(
      alpha = alpha_grid,
      rejection_rate = sapply(alpha_grid, function(a) mean(A$p_cheng <= a, na.rm = TRUE)),
      method = "Cheng-GM"
    )
  ))
  
  curve_dt[, scenario_id := scenario_id]
  curve_dt[, dtag := dtag]
  curve_dt[, signal_d := dval]
  
  curve_dt
}

all_curves <- list()

for (ii in seq_len(nrow(d_info))) {
  
  dtag <- d_info$dtag[ii]
  dval <- d_info$dval[ii]
  
  for (s in sprintf("S%02d", 1:40)) {
    
    cat("Creating:", s, dtag, "\n")
    
    curve_dt <- make_curve(s, dtag, dval)
    all_curves[[length(all_curves) + 1]] <- curve_dt
    
    p <- ggplot(curve_dt, aes(x = alpha, y = rejection_rate, color = method)) +
      geom_line(linewidth = 1) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      labs(
        title = paste0("Size-power curve: ", s, ", d = ", dval),
        x = expression(alpha),
        y = "Empirical rejection rate",
        color = "Method"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "bottom",
        panel.grid.minor = element_blank()
      )
    
    ggsave(
      filename = file.path(out_dir, paste0("size_power_", s, "_", dtag, ".png")),
      plot = p,
      width = 6,
      height = 4.5,
      dpi = 300
    )
    
    ggsave(
      filename = file.path(out_dir, paste0("size_power_", s, "_", dtag, ".pdf")),
      plot = p,
      width = 6,
      height = 4.5
    )
  }
}

all_curves_dt <- rbindlist(all_curves, fill = TRUE)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

fwrite(
  all_curves_dt,
  "results/tables/size_power_curve_data_all.csv"
)

cat("\nDONE.\n")
cat("Size-power plots saved to:", out_dir, "\n")
cat("Expected PNG files: 120\n")
cat("Expected PDF files: 120\n")
cat("Curve data saved to: results/tables/size_power_curve_data_all.csv\n")