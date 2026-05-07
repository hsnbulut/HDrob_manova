# ============================================================
# Robustness p-value plots - corrected version
# ============================================================

library(ggplot2)
library(dplyr)

dir.create("results/figures/robustness/pvalue_plots/all",
           recursive = TRUE, showWarnings = FALSE)

dat <- read.csv("results/tables/robustness_main_raw_long.csv")

# Scenario labels
dat$scenario_id <- as.character(dat$scenario_id)
dat$scenario_id <- sub("^S0+", "S", dat$scenario_id)

# Method order
dat$method <- factor(
  dat$method,
  levels = c("Proposed", "Schott", "Cheng-GM")
)

make_pvalue_data <- function(data_sub) {
  data_sub %>%
    group_by(method) %>%
    arrange(p_value, .by_group = TRUE) %>%
    mutate(
      i = row_number(),
      n = n(),
      expected = i / (n + 1),
      observed = p_value
    ) %>%
    ungroup()
}

make_pvalue_plot <- function(scen, eps_value) {
  
  tmp <- dat %>%
    filter(scenario_id == scen, abs(eps - eps_value) < 1e-8)
  
  if (nrow(tmp) == 0) {
    stop("No data found for scenario = ", scen, " and eps = ", eps_value)
  }
  
  pp <- make_pvalue_data(tmp)
  
  ggplot(pp, aes(x = expected, y = observed, color = method)) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_step(linewidth = 0.7) +
    coord_cartesian(xlim = c(0, 0.20), ylim = c(0, 0.20)) +
    labs(
      title = paste0("P-value plot - ", scen, ", eps = ", eps_value),
      x = "Expected uniform order statistics",
      y = "Observed ordered p-values",
      color = "Method"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    )
}


scenarios <- sort(
  unique(dat$scenario_id),
  index.return = FALSE
)

scenarios <- scenarios[order(as.numeric(sub("S", "", scenarios)))]

eps_values <- sort(unique(dat$eps))

for (eps_value in eps_values) {
  for (scen in scenarios) {
    
    p <- make_pvalue_plot(scen, eps_value)
    
    eps_tag <- ifelse(abs(eps_value - 0.10) < 1e-8, "eps10", "eps20")
    
    ggsave(
      filename = paste0(
        "results/figures/robustness/pvalue_plots/all/",
        "robustness_pvalue_", scen, "_", eps_tag, ".png"
      ),
      plot = p,
      width = 6.2,
      height = 4.8,
      dpi = 300
    )
  }
}

pdf(
  "results/figures/robustness/pvalue_plots/robustness_pvalue_plots_all.pdf",
  width = 6.2,
  height = 4.8
)

for (eps_value in eps_values) {
  for (scen in scenarios) {
    print(make_pvalue_plot(scen, eps_value))
  }
}

dev.off()

selected <- data.frame(
  scenario_id = c("S9", "S16", "S9", "S16"),
  eps = c(0.10, 0.10, 0.20, 0.20)
)

rep_dat <- dat %>%
  inner_join(selected, by = c("scenario_id", "eps")) %>%
  mutate(
    panel_label = paste0(scenario_id, ", \u03b5 = ", eps)
  )

rep_pp <- rep_dat %>%
  group_by(panel_label, method) %>%
  arrange(p_value, .by_group = TRUE) %>%
  mutate(
    i = row_number(),
    n = n(),
    expected = i / (n + 1),
    observed = p_value
  ) %>%
  ungroup()

p_rep <- ggplot(rep_pp, aes(x = expected, y = observed, color = method)) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_step(linewidth = 0.7) +
  coord_cartesian(xlim = c(0, 0.20), ylim = c(0, 0.20)) +
  facet_wrap(~ panel_label, ncol = 2) +
  labs(
    x = "Expected uniform order statistics",
    y = "Observed ordered p-values",
    color = "Method"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

ggsave(
  "results/figures/robustness/robustness_pvalue_representative.pdf",
  p_rep,
  width = 7.2,
  height = 5.8
)

ggsave(
  "results/figures/robustness/robustness_pvalue_representative.png",
  p_rep,
  width = 7.2,
  height = 5.8,
  dpi = 300
)

p_rep