# ============================================================
# summarize_power_outputs.R
# Creates CSV summaries, LaTeX tables, and size-power plots
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

raw_dir <- "results/raw_chunks/power"
table_dir <- "results/tables"
fig_dir <- "results/figures/power"
sp_dir <- file.path(fig_dir, "size_power")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sp_dir, recursive = TRUE, showWarnings = FALSE)

alpha_level <- 0.05
d_tags <- c("d20", "d40", "d60")
d_vals <- c(d20 = 0.2, d40 = 0.4, d60 = 0.6)

mcse <- function(phat, R) {
  sqrt(phat * (1 - phat) / R)
}

# ------------------------------------------------------------
# 1. Combine raw chunks to scenario-level summaries
# ------------------------------------------------------------
combine_one_d <- function(dtag) {
  
  files <- list.files(
    raw_dir,
    pattern = paste0("_", dtag, "_block[0-9]+\\.rds$"),
    full.names = TRUE
  )
  
  if (length(files) == 0) {
    warning("No files found for ", dtag)
    return(NULL)
  }
  
  all_dt <- rbindlist(lapply(files, function(f) {
    obj <- readRDS(f)
    dt <- as.data.table(obj$results)
    meta <- as.list(obj$meta)
    
    dt[, scenario_id := meta$scenario_id]
    dt[, G := meta$G]
    dt[, nk := meta$nk]
    dt[, p := meta$p]
    dt[, cov_type := meta$cov_type]
    dt[, design := meta$design]
    dt[, signal_d := meta$signal_d]
    dt[, block_id := obj$block_id]
    dt
  }), fill = TRUE)
  
  summary_dt <- all_dt[
    ,
    {
      pp <- mean(p_proposed < alpha_level, na.rm = TRUE)
      ps <- mean(p_schott   < alpha_level, na.rm = TRUE)
      pc <- mean(p_cheng    < alpha_level, na.rm = TRUE)
      
      .(
        R = .N,
        proposed_power = pp,
        proposed_mcse  = mcse(pp, sum(!is.na(p_proposed))),
        schott_power   = ps,
        schott_mcse    = mcse(ps, sum(!is.na(p_schott))),
        cheng_power    = pc,
        cheng_mcse     = mcse(pc, sum(!is.na(p_cheng))),
        na_proposed    = sum(is.na(p_proposed)),
        na_schott      = sum(is.na(p_schott)),
        na_cheng       = sum(is.na(p_cheng))
      )
    },
    by = .(scenario_id, G, nk, p, cov_type, design, signal_d)
  ]
  
  summary_dt[, scenario_num := as.integer(gsub("S", "", scenario_id))]
  setorder(summary_dt, scenario_num)
  summary_dt[, scenario_num := NULL]
  
  out_csv <- file.path(table_dir, paste0("power_summary_", dtag, ".csv"))
  fwrite(summary_dt, out_csv)
  
  cat("Written:", out_csv, "\n")
  cat("Scenarios:", nrow(summary_dt), " | Total R:", sum(summary_dt$R), "\n")
  
  incomplete <- summary_dt[R != 1000]
  if (nrow(incomplete) > 0) {
    cat("WARNING incomplete scenarios for", dtag, "\n")
    print(incomplete[, .(scenario_id, R)])
  }
  
  summary_dt
}

tabs <- rbindlist(lapply(d_tags, combine_one_d), fill = TRUE)
fwrite(tabs, file.path(table_dir, "power_summary_all.csv"))
cat("Written:", file.path(table_dir, "power_summary_all.csv"), "\n")

# ------------------------------------------------------------
# 2. Paper-ready CSV tables
# ------------------------------------------------------------
make_paper_csv <- function(dtag) {
  dval <- d_vals[[dtag]]
  dt <- tabs[abs(signal_d - dval) < 1e-12]
  
  paper_dt <- dt[
    ,
    .(
      Scenario = scenario_id,
      G = G,
      `Group Sizes` = paste0("(", nk, ")"),
      p = p,
      Structure = ifelse(cov_type == "identity", "I", "AR(1)"),
      `Proposed Power` = round(proposed_power, 3),
      `Proposed MCSE`  = round(proposed_mcse, 4),
      `Schott Power`   = round(schott_power, 3),
      `Schott MCSE`    = round(schott_mcse, 4),
      `Cheng-GM Power` = round(cheng_power, 3),
      `Cheng-GM MCSE`  = round(cheng_mcse, 4)
    )
  ]
  
  out_csv <- file.path(table_dir, paste0("table_power_", dtag, "_paper.csv"))
  fwrite(paper_dt, out_csv)
  cat("Written:", out_csv, "\n")
  
  paper_dt
}

paper_tabs <- lapply(d_tags, make_paper_csv)

# ------------------------------------------------------------
# 3. LaTeX tables matching manuscript format
# ------------------------------------------------------------
write_latex_table <- function(paper_dt, dtag) {
  
  out_tex <- file.path(table_dir, paste0("table_power_", dtag, ".tex"))
  
  lines <- c(
    "\\begin{table}[ht]",
    "\\centering",
    paste0("\\caption{Empirical power and Monte Carlo standard errors under dense alternatives ($d = ",
           d_vals[[dtag]], "$)}"),
    paste0("\\label{tab:power_", dtag, "}"),
    "\\scriptsize",
    "\\begin{tabular}{ccccccccccc}",
    "\\hline",
    "Scenario & G & $(n_1,\\ldots,n_G)$ & $p$ & Structure & \\multicolumn{2}{c}{Proposed} & \\multicolumn{2}{c}{Schott (2007)} & \\multicolumn{2}{c}{Cheng-GM} \\\\",
    " & & & & & Power & MCSE & Power & MCSE & Power & MCSE \\\\",
    "\\hline"
  )
  
  for (i in seq_len(nrow(paper_dt))) {
    r <- paper_dt[i]
    lines <- c(lines, paste0(
      r$Scenario, " & ",
      r$G, " & ",
      r$`Group Sizes`, " & ",
      r$p, " & ",
      r$Structure, " & ",
      sprintf("%.3f", r$`Proposed Power`), " & ",
      sprintf("%.4f", r$`Proposed MCSE`), " & ",
      sprintf("%.3f", r$`Schott Power`), " & ",
      sprintf("%.4f", r$`Schott MCSE`), " & ",
      sprintf("%.3f", r$`Cheng-GM Power`), " & ",
      sprintf("%.4f", r$`Cheng-GM MCSE`), " \\\\"
    ))
  }
  
  lines <- c(lines, "\\hline", "\\end{tabular}", "\\end{table}")
  
  writeLines(lines, out_tex)
  cat("Written:", out_tex, "\n")
}

for (i in seq_along(d_tags)) {
  write_latex_table(paper_tabs[[i]], d_tags[i])
}

# ------------------------------------------------------------
# 4. Average power CSV and figure
# ------------------------------------------------------------
avg_power <- tabs[
  ,
  .(
    proposed_power = mean(proposed_power),
    schott_power   = mean(schott_power),
    cheng_power    = mean(cheng_power)
  ),
  by = signal_d
]

fwrite(avg_power, file.path(table_dir, "power_average_by_d.csv"))
cat("Written:", file.path(table_dir, "power_average_by_d.csv"), "\n")

avg_long <- melt(
  avg_power,
  id.vars = "signal_d",
  variable.name = "method",
  value.name = "power"
)

avg_long[, method := factor(
  method,
  levels = c("proposed_power", "schott_power", "cheng_power"),
  labels = c("Proposed", "Schott (2007)", "Cheng-GM")
)]

p_avg <- ggplot(avg_long, aes(x = signal_d, y = power, color = method, group = method)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = c(0.2, 0.4, 0.6)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Signal strength (d)", y = "Average empirical power", color = "Method") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "power_curve_average.png"), p_avg, width = 7, height = 4.5, dpi = 300)
ggsave(file.path(fig_dir, "power_curve_average.pdf"), p_avg, width = 7, height = 4.5)

# ------------------------------------------------------------
# 5. Size-power curve data and plots for all scenarios
# ------------------------------------------------------------
alpha_grid <- seq(0, 1, length.out = 201)

make_size_power_one <- function(scenario, dtag) {
  
  files <- list.files(
    raw_dir,
    pattern = paste0("^", scenario, "_", dtag, "_block[0-9]+\\.rds$"),
    full.names = TRUE
  )
  
  if (length(files) != 10) {
    warning("Expected 10 files but found ", length(files), " for ", scenario, "_", dtag)
  }
  
  A <- rbindlist(lapply(files, function(f) readRDS(f)$results), fill = TRUE)
  
  out <- rbindlist(list(
    data.table(alpha = alpha_grid, rejection_rate = sapply(alpha_grid, function(a) mean(A$p_proposed <= a, na.rm = TRUE)), method = "Proposed"),
    data.table(alpha = alpha_grid, rejection_rate = sapply(alpha_grid, function(a) mean(A$p_schott   <= a, na.rm = TRUE)), method = "Schott (2007)"),
    data.table(alpha = alpha_grid, rejection_rate = sapply(alpha_grid, function(a) mean(A$p_cheng    <= a, na.rm = TRUE)), method = "Cheng-GM")
  ))
  
  out[, scenario_id := scenario]
  out[, signal_d := d_vals[[dtag]]]
  out[, dtag := dtag]
  
  out
}

all_curve_dt <- rbindlist(
  lapply(d_tags, function(dtag) {
    rbindlist(lapply(sprintf("S%02d", 1:40), function(s) {
      make_size_power_one(s, dtag)
    }))
  }),
  fill = TRUE
)

fwrite(all_curve_dt, file.path(table_dir, "size_power_curve_data_all.csv"))
cat("Written:", file.path(table_dir, "size_power_curve_data_all.csv"), "\n")

for (dtag in d_tags) {
  for (s in sprintf("S%02d", 1:40)) {
    
    curve_dt <- all_curve_dt[scenario_id == s & dtag == !!dtag]
    
    p <- ggplot(curve_dt, aes(x = alpha, y = rejection_rate, color = method)) +
      geom_line(linewidth = 1) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      labs(
        title = paste0("Size-power curve: ", s, ", d = ", d_vals[[dtag]]),
        x = expression(alpha),
        y = "Empirical rejection rate",
        color = "Method"
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom", panel.grid.minor = element_blank())
    
    ggsave(file.path(sp_dir, paste0("size_power_", s, "_", dtag, ".png")), p, width = 6, height = 4.5, dpi = 300)
    ggsave(file.path(sp_dir, paste0("size_power_", s, "_", dtag, ".pdf")), p, width = 6, height = 4.5)
  }
}

cat("\nDONE.\n")
cat("CSV outputs are in:", table_dir, "\n")
cat("Figures are in:", fig_dir, "\n")