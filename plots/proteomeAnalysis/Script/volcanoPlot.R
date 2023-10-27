# the goal of this script is to visualise Peyton's dataset
library(ggplot2)
library(ggrepel)

# functions ----

make_shortnames_unique <- function(x) {
  y <- x
  
  for(i in 1:length(y)) {
    z <- unlist(strsplit(y[i],";"))
    if(length(z) == 1) {
      next
    }
    for(j in 1:length(z)) {
      if(length(unique(y)) < length(unique(append(y,z[j])))) {
        y[i] <- z[j]
        break
      }
    }
  }
  
  return(y)
}

# script ----
pe <- read.delim("Data/rankTable_INVvsControlPE.txt")

# some names to be disambiguated, e.g. RAB8A;RAB8B from RAB8A
pe$uSHORTNAME <- make_shortnames_unique(pe$so_SHORTNAME)
# label the proteins of interest
pe$interest <- ifelse(pe$so_colorWave == 3, pe$uSHORTNAME, "")

# zoomed out view (no protein labels)
p1 <- ggplot(pe, aes(x = log2(so_ratioWave), y = -log10(so_allTWave))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey") +
  geom_point(aes(colour = as.factor(so_colorWave)), size = 1, shape = 16, alpha = 0.5) +
  scale_colour_manual(values = c("0" = "#808080", "1" = "#ff8080", "2" = "#8080ff", "3" = "#ff80ff")) +
  scale_x_continuous(limits = c(-11,11), breaks = seq(-10,10,2)) +
  scale_y_continuous(limits = c(0,20), breaks = seq(0,20,2), labels = function(i) 10^-i) +
  xlab(bquote(INV~-~Control~(Log[2]))) +
  ylab("P-value") +
  theme_classic(8) +
  theme(legend.position = "none")

ggsave("Output/Plots/vp_pe.pdf", p1, width = 80, height = 100, units = "mm", family = "Helvetica")

# zoom in (with protein labels)
p2 <- ggplot(pe, aes(x = log2(so_ratioWave), y = -log10(so_allTWave), label = interest)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey") +
  geom_point(aes(colour = as.factor(so_colorWave)), size = 1, shape = 16, alpha = 0.5) +
  scale_colour_manual(values = c("0" = "#808080", "1" = "#ff8080", "2" = "#8080ff", "3" = "#ff80ff")) +
  geom_text_repel(size = 1.5, max.overlaps = 15) +
  scale_x_continuous(limits = c(1,8), breaks = seq(1,8,1)) +
  scale_y_continuous(limits = c(0,12), breaks = seq(0,12,2), labels = function(i) 10^-i) +
  xlab(bquote(INV~-~Control~(Log[2]))) +
  ylab("P-value") +
  theme_classic(8) +
  theme(legend.position = "none")

ggsave("Output/Plots/vp_pe2.pdf", p2, width = 80, height = 100, units = "mm", family = "Helvetica")
