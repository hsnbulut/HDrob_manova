library(data.table)
library(ggplot2)

tab <- fread("results/tables/power_summary_all.csv")

avg <- tab[, .(
  proposed = mean(proposed_power),
  schott   = mean(schott_power),
  cheng    = mean(cheng_power)
), by=signal_d]

avg_long <- melt(avg, id="signal_d")

p <- ggplot(avg_long,
            aes(signal_d,value,color=variable)) +
  geom_line(size=1.2) +
  geom_point(size=3) +
  labs(x="d",y="Power") +
  theme_minimal()

ggsave("results/figures/power_curve.png",p,6,4)