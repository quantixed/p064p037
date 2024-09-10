library(ggplot2)
library(cowplot)
library(patchwork)

## Script ----

prots <- data.frame(
  n = c(601, 3622),
  label = c("INV", "Total")
)

p1 <- ggplot(prots) +
  geom_col(aes(x = label, y = n, fill = label), position = "dodge") +
  scale_fill_manual(values = c("#009988", "#bbbbbb")) +
  labs(x = "", y = "Proteins") +
  lims(y = c(0,4000)) +
  theme_cowplot(9) +
  theme(legend.position = "none")

secreted <- data.frame(
  n = c(63 / 601, 259 / 3622),
  label = c("INV", "Total")
)

p2 <- ggplot(secreted) +
  geom_col(aes(x = label, y = n, fill = label), position = "dodge") +
  scale_fill_manual(values = c("#009988", "#bbbbbb")) +
  labs(x = "Secreted", y = "Fraction") +
  scale_y_continuous(breaks = seq(0,0.12,0.02)) +
  theme_cowplot(9) +
  theme(legend.position = "none")

tm <- data.frame(
  n = c(168 / 601, 433 / 3622),
  label = c("INV", "Total")
)

p3 <- ggplot(tm) +
  geom_col(aes(x = label, y = n, fill = label), position = "dodge") +
  scale_fill_manual(values = c("#009988", "#bbbbbb")) +
  labs(x = "TMDs", y = "Fraction") +
  lims(y = c(0,0.3)) +
  theme_cowplot(9) +
  theme(legend.position = "none")

q1 <- p1 | p2 | p3
ggsave("Output/Plots/bars.pdf", q1, width = 86, height = 43, units = "mm")
 
