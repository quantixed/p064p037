## Script to plot ATG9 on Golgi data
## Author: Mary Fesenko

library(ggplot2)
library(ggforce)
library(cowplot)

# list all files in the Data directory
# each group of experiment files must be in a subdirectory
filelist <- list.files("Data", recursive = TRUE, full.names = TRUE)
# load all files and rbind into big dataframe
# add columns to each file to identify the file and expt
df <- do.call(rbind, lapply(filelist, function(x) {
  temp <- read.csv(x);
  temp$file <- basename(x);
  temp$expt <- dirname(x);
  temp}))
# file column has name of the file, name is of the form foo_bar_1.csv, extract foo and bar into two columns
df$condA <- sapply(strsplit(df$file, "_"), "[", 1)
df$condB <- sapply(strsplit(df$file, "_"), "[", 2)
df$cell <- sapply(strsplit(df$file, "_"), "[", 3)

# When you use "[" as an argument, it refers to the extraction function [, 
# if its first argument is set to X, the X-ths element of each list would be extracted. 
# so, "[", 1 would extract the first element in a list

df[df == "starv"] <- "starved"
df[df == "siGL2"] <- "siCtrl"

## Graphing ----

# calculate experiment means
summary_df <- aggregate(Mean ~ condA + condB + expt, data = df, FUN = mean)

# combine condA and condB into one column
df$condAB <- paste(df$condA, df$condB, sep = "_")
summary_df$condAB <- paste(summary_df$condA, summary_df$condB, sep = "_")

# Plotting
ggplot() +
  geom_sina(data = df, aes(x = condAB, y = Mean, colour = expt, shape = condA), alpha = 0.5, position = "auto", size = 0.8, maxwidth = 0.3) +
  geom_point(data = summary_df, aes(x = condAB, y = Mean, fill = expt), shape = 22, size = 1.5, stroke = 0.5, alpha = 0.7) +
  scale_color_manual(values = c("#4477aa", "#ccbb44", "#ee6677", "#000000")) +
  scale_fill_manual(values = c("#4477aa", "#ccbb44", "#ee6677", "#000000")) +
  scale_shape_manual(values = c(1, 16)) +
  labs(x = "", y = "ATG9A immunofluorescence (A.U.)") +
  lims(y = c(0, NA)) +
  theme_cowplot(9) +
  theme(legend.position = "none")
ggsave("Output/Plots/ATG9onGolgi_plot.pdf", width = 88, height = 50, units = "mm")

lm <- aov(Mean ~ condA * condB, data = summary_df)
summary(lm)
TukeyHSD(lm, conf.level = 0.95)
