library(data.table)

tab <- fread("results/tables/power/power_summary_all.csv")

format_table <- function(dval) {
  dt <- tab[signal_d == dval]
  
  dt[, scenario := scenario_id]
  
  dt[, proposed := sprintf("%.3f", proposed_power)]
  dt[, schott   := sprintf("%.3f", schott_power)]
  dt[, cheng    := sprintf("%.3f", cheng_power)]
  
  dt[, .(scenario, G, nk, p, cov_type, proposed, schott, cheng)]
}

write_latex <- function(dt, file) {
  
  cat("\\begin{tabular}{cccccccc}\n", file=file)
  cat("\\hline\n", file=file, append=TRUE)
  cat("Scenario & G & nk & p & Cov & Proposed & Schott & Cheng \\\\\n", file=file, append=TRUE)
  cat("\\hline\n", file=file, append=TRUE)
  
  for(i in 1:nrow(dt)) {
    row <- dt[i]
    line <- paste(
      row$scenario, "&",
      row$G, "&",
      row$nk, "&",
      row$p, "&",
      row$cov_type, "&",
      row$proposed, "&",
      row$schott, "&",
      row$cheng, "\\\\"
    )
    cat(line, "\n", file=file, append=TRUE)
  }
  
  cat("\\hline\n\\end{tabular}\n", file=file, append=TRUE)
}

write_latex(format_table(0.2), "results/tables/table_power_d20.tex")
write_latex(format_table(0.4), "results/tables/table_power_d40.tex")
write_latex(format_table(0.6), "results/tables/table_power_d60.tex")

cat("Latex tables created.\n")