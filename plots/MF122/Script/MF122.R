library(dplyr)
library(tidyr)
library(ggplot2)
library(ggforce)
library(cowplot)

## Functions ----

print_stats <- function(df, comb1, comb2) {
  query <- paste0(comb1, "-", comb2)
  pval <- df[query, "p adj"]
  if (is.na(pval)) {
    query <- paste0(comb2, "-", comb1)
    pval <- df[query, "p adj"]
  }
  return(pval)
}

## Script ----


# list all files in the Data directory
# each group of experiment files must be in a subdirectory
filelist <- list.files("Data", recursive = TRUE, full.names = TRUE)
# load all files and rbind into big dataframe
# add columns to each file to identify the file and expt
df <- do.call(rbind, lapply(filelist, function(x) {
  temp <- read.csv(x)
  temp$file <- basename(x)
  temp$expt <- dirname(x)
  temp}))
# file column has name of the file, name is of the form foo_bar_1.csv,
# extract foo and bar into two columns
df$condA <- sapply(strsplit(df$file, "_"), "[", 1)
df$condB <- sapply(strsplit(df$file, "_"), "[", 2)
df$cell <- sapply(strsplit(df$file, "_"), "[", 3)
# remove fed+Baf as we will not plot this condition
df <- df[df$condB != "fed+Baf", ]
# remove puncta with less than 200 pixels
df <- df[df$PIXEL_COUNT > 200, ]
# scale pixels to Volume..micron.3.
df$Volume..micron.3. <- df$PIXEL_COUNT * 0.0392857^3

# combine condA and condB in column called condAB
df$condAB <- paste(df$condA, df$condB, sep = "_")
# make a column for labels
proper <- c("siGL2_fed+DMSO" = "Fed - BafA1\nsiCtrl",
            "siGL2_fed+Baf" = "Fed + BafA1\nsiCtrl",
            "siGL2_starved+DMSO" = "Starved - BafA1\nsiCtrl",
            "siGL2_starved+Baf" = "Starved + BafA1\nsiCtrl",
            "siTPD54_fed+DMSO" = "Fed - BafA1\nsiTPD54",
            "siTPD54_fed+Baf" = "Fed + BafA1\nsiTPD54",
            "siTPD54_starved+DMSO" = "Starved - BafA1\nsiTPD54",
            "siTPD54_starved+Baf" = "Starved + BafA1\nsiTPD54")
df$condAB <- factor(df$condAB, levels = names(proper), labels = proper)

# per cell summary
summary_df <- df %>%
  group_by(expt, condA, condB, cell) %>%
  summarise(medianvol = median(Volume..micron.3.),
            meanvol = mean(Volume..micron.3.),
            n = n(),
            totalvol = sum(Volume..micron.3.))
summary_df$condAB <- paste(summary_df$condA, summary_df$condB, sep = "_")
summary_df$condAB <- factor(summary_df$condAB,
                            levels = names(proper), labels = proper)

# per expt summary
expt_summary_df <- summary_df %>%
  group_by(expt, condA, condB) %>%
  summarise(overall_mean = mean(meanvol),
            overall_median = mean(medianvol),
            overall_n = mean(n),
            overall_total = mean(totalvol))
expt_summary_df$condAB <- paste(expt_summary_df$condA,
                                expt_summary_df$condB, sep = "_")
expt_summary_df$condAB <- factor(expt_summary_df$condAB,
                                 levels = names(proper), labels = proper)

ggplot() +
  geom_sina(data = summary_df, aes(x = condAB, y = meanvol,
                                   colour = expt, shape = condA),
            alpha = 0.5, position = "auto", size = 0.8, maxwidth = 0.3) +
  geom_point(data = expt_summary_df, aes(x = condAB, y = overall_mean,
                                         fill = expt),
             shape = 22, size = 1.5, stroke = 0.5, alpha = 0.7) +
  scale_color_manual(values = c("#4477aa", "#ccbb44", "#ee6677", "#000000")) +
  scale_fill_manual(values = c("#4477aa", "#ccbb44", "#ee6677", "#000000")) +
  scale_shape_manual(values = c(1, 16)) +
  labs(x = "", y = "Mean volume (µm³)") +
  lims(y = c(0, NA)) +
  theme_cowplot(9) +
  theme(legend.position = "none")
ggsave("Output/Plots/meanvol.pdf", width = 88, height = 50, units = "mm")

ggplot() +
  geom_sina(data = summary_df, aes(x = condAB, y = n,
                                   colour = expt, shape = condA),
            alpha = 0.5, position = "auto", size = 0.8, maxwidth = 0.3) +
  geom_point(data = expt_summary_df, aes(x = condAB, y = overall_n,
                                         fill = expt),
             shape = 22, size = 1.5, stroke = 0.5, alpha = 0.7) +
  scale_color_manual(values = c("#4477aa", "#ccbb44", "#ee6677", "#000000")) +
  scale_fill_manual(values = c("#4477aa", "#ccbb44", "#ee6677", "#000000")) +
  scale_shape_manual(values = c(1, 16)) +
  labs(x = "", y = "LC3 puncta per cell") +
  lims(y = c(0, NA)) +
  theme_cowplot(9) +
  theme(legend.position = "none")
ggsave("Output/Plots/n.pdf", width = 88, height = 50, units = "mm")


## Stats ----

lm_meanvol <- aov(overall_mean ~ condA * condB, data = expt_summary_df)
summary(lm_meanvol)
stats_vol <- TukeyHSD(lm_meanvol, conf.level = 0.95)

lm_count <- aov(overall_n ~ condA * condB, data = expt_summary_df)
summary(lm_count)
stats_n <- TukeyHSD(lm_count, conf.level = 0.95)

stats_vol <- as.data.frame(stats_vol$`condA:condB`)
stats_n <- as.data.frame(stats_n$`condA:condB`)

# print the p adj values of interest
print_stats(stats_vol, "siGL2:starved+DMSO", "siGL2:fed+DMSO")
print_stats(stats_vol, "siTPD54:starved+DMSO", "siGL2:fed+DMSO")
print_stats(stats_n, "siGL2:starved+DMSO", "siGL2:starved+Baf")
print_stats(stats_n, "siTPD54:starved+DMSO", "siTPD54:starved+Baf")

# difference in overall_mean for each expt compare "Starved - BafA1\nsiCtrl" and  "Starved - BafA1\nsiTPD54"
expt_summary_df %>%
  ungroup() %>%
  select(condAB, expt, overall_mean) %>%
  pivot_wider(names_from = condAB, values_from = overall_mean) %>%
  mutate(diff = (`Starved - BafA1\nsiCtrl` - `Starved - BafA1\nsiTPD54`) / `Starved - BafA1\nsiCtrl` * 100)
