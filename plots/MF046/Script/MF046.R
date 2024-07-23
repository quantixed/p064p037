library(ggplot2)
library(ggforce)
library(dplyr)
library(tidyr)

# function to load in the data into a big data frame
# relies on files being called `allResults_MF046.csv` etc
# can be of form `foo_expid.csv` the code picks up `expid` from this filename
# uses column names as defined by Mary in her analysis script.
# this is the Mito & Cyto version of the ATG9 and TPD54 analysis.


## Functions ----
process_data <- function() {
  all_files <- list.files("Data/mitocyto", pattern = "*.csv", full.names = TRUE)
  df_all <- read.csv(all_files[1], header = TRUE)
  df_all <- remove_imagej_rows_if_there(df_all)
  fname <- tools::file_path_sans_ext(basename(all_files[1]))
  fname <- gsub("^.*?_", "", fname)
  df_all$expt <- fname
  for (filename in all_files[-1]) {
    df_temp <- read.csv(filename)
    fname <- tools::file_path_sans_ext(basename(filename))
    fname <- gsub("^.*?_", "", fname)
    df_temp$expt <- fname
    df_temp <- remove_imagej_rows_if_there(df_temp)
    df_all <- rbind(df_all, df_temp)
  }
  return(df_all)
}

remove_imagej_rows_if_there <- function(df) {
  testIf <- names(df)
  if (testIf[1] == "X") {
    df$X <- NULL
  }
  
  return(df)
}

getRatio <- function(df, s, label) {
  subdf <- df[df$id == s & df$Channel == label, ]
  ratio <- subdf$Mean.BG[
    grep("^post", subdf$Condition)[1]] /
    subdf$Mean.BG[grep("^pre", subdf$Condition)[1]]
  
  return(ratio)
}

calculate_ratios <- function(df, uid) {
  newdf <- data.frame(
    uid = uid,
    rTPD54 = rep(0, length(uid)),
    rATG9 = rep(0, length(uid)),
    condA = rep("", length(uid)),
    condB = rep("", length(uid)),
    expt = rep("", length(uid))
  )
  for (i in seq_along(uid)) {
    newdf$rTPD54[i] <- getRatio(df, uid[i], "TPD54")
    newdf$rATG9[i] <- getRatio(df, uid[i], "ATG9")
    
    newdf$condA[i] <- unlist(strsplit(uid[i], "_"))[2]
    newdf$condB[i] <- unlist(strsplit(uid[i], "_"))[3]
    newdf$expt[i] <- unlist(strsplit(uid[i], "_"))[1]
  }
  
  newdf$condA <- as.factor(newdf$condA)
  newdf$condA <- factor(newdf$condA, levels = c("WT", "Mut", "Cl35"))
  newdf$condB <- as.factor(newdf$condB)
  newdf$expt <- as.factor(newdf$expt)
  
  return(newdf)
}

calculate_loss_ratios <- function(df) {
  newdf <- data.frame(
    rTPD54 = rep(0, length(uid)),
    rATG9 = rep(0, length(uid)),
    condA = rep("", length(uid)),
    condB = rep("", length(uid)),
    expt = rep("", length(uid))
  )
  for (i in seq_along(uid)) {
    newdf$rTPD54[i] <- getRatio(df, uid[i], "TPD54")
    newdf$rATG9[i] <- getRatio(df, uid[i], "ATG9")
    
    newdf$condA[i] <- unlist(strsplit(uid[i], "_"))[2]
    newdf$condB[i] <- unlist(strsplit(uid[i], "_"))[3]
    newdf$expt[i] <- unlist(strsplit(uid[i], "_"))[1]
  }
  
  newdf$condA <- as.factor(newdf$condA)
  newdf$condA <- factor(newdf$condA, levels = c("WT", "Mut", "Cl35"))
  newdf$condB <- as.factor(newdf$condB)
  newdf$expt <- as.factor(newdf$expt)
  
  return(newdf)
}

# function for labelling x-axis in a more meaningful way
meaningfully_x_label <- function(p) {
  # p is ggplot
  p <- p + scale_x_discrete(labels = c(
    "WT:cont" = "TPD54\nFKBP-mCherry",
    "WT:test" = "TPD54\nATG9A-FKBP-mCherry",
    "Mut:cont" = "R159E\nFKBP-mCherry",
    "Mut:test" = "R159E\nATG9A-FKBP-mCherry",
    "Cl35:cont" = "EndoTPD54\nFKBP-mCherry",
    "Cl35:test" = "EndoTPD54\nATG9A-FKBP-mCherry"
  ))
  return(p)
}

## Script ----

data <- process_data()
# Condition can have trailing /
data$Condition <- gsub("/$", "", data$Condition)
# PrePost extracted
data$PrePost <- sapply(strsplit(data$Condition, "_"), "[[", 1)
# extract conditions
data$CondA <- sapply(strsplit(data$Condition, "_"), "[[", 2)
data$CondB <- sapply(strsplit(data$Condition, "_"), "[[", 3)
# generate id with Condition(pre/post) and Channel left as variables
# data$id <- paste(data$expt, data$Condition, data$File, data$cell, sep = "_")
data$id <- paste(data$expt, data$CondA, data$CondB, data$File, sep = "_")
# take unique list of each cell that was image
uid <- unique(data$id)
# now we need to separate mito and cyto measurements
# they are designated 0 (mito) or 1 (cyto) in the Slice column
mito_data <- data[data$Slice == 0, ]
cyto_data <- data[data$Slice == 1, ]
# do the postOverPre calculations
mito_summary_df <- calculate_ratios(mito_data, uid)
cyto_summary_df <- calculate_ratios(cyto_data, uid)

# filter for ATG9 actually rerouting - relies on the rows being the same!
cyto_summary_df <- cyto_summary_df[mito_summary_df$rATG9 > 1, ]
mito_summary_df <- mito_summary_df[mito_summary_df$rATG9 > 1, ]

condition_names <- c(
  "cont" = "FKBP-mCherry",
  "test" = "ATG9A-FKBP-mCherry",
  "ATG9" = "mCherry",
  "TPD54" = "TPD54"
)

# cyto data can be used to compare post as a group with pre as a group per expt
# we take the average of each channel per group per expt and then find the
# ratio of minus/plus treatment
cytoloss <- cyto_data %>%
  group_by(Channel, expt, PrePost, CondA, CondB) %>%
  summarise(mean = mean(Mean.BG)) %>%
  pivot_wider(names_from = PrePost, values_from = mean)
cytoloss$ratio <- cytoloss$postRapalog / cytoloss$preRapalog
cytoloss$CondA <- factor(cytoloss$CondA, levels = c("WT", "Mut", "Cl35"))

mitocytocomp <- data %>%
  pivot_wider(values_from = "Mean.BG", id_cols = c("id", "Condition", "PrePost", "CondA", "CondB", "expt"), names_from = c("Channel", "Slice"))
mitocytocomp$TPD54 <- mitocytocomp$TPD54_0 / (mitocytocomp$TPD54_0 + mitocytocomp$TPD54_1)
mitocytocomp$ATG9 <- mitocytocomp$ATG9_0 / (mitocytocomp$ATG9_0 + mitocytocomp$ATG9_1)
mitocytocomp <- mitocytocomp %>%
  filter(ATG9 > 0.5) %>%
  select(-c(TPD54_0, TPD54_1, ATG9_0, ATG9_1)) %>%
  pivot_longer(cols = c("TPD54", "ATG9"), names_to = "Channel", values_to = "mmc")

sp_mmc <- mitocytocomp %>%
  group_by(CondA, PrePost, Channel, expt) %>%
  summarise(
    n = n(),
    mmc = mean(mmc)
  )

###########################################
## generate plots
## colour/shape scheme is
## shape: open control, closed treatment
## colours: WT, Mut, Cl35 = "#117733", "#999933", "#44aa99"
## colours: Cont, WT, Mut = "#88ccee", "#117733", "#999933"
## colours: superplot:" #4477aa", "#ccbb44", "#ee6677"
## This next part will now generate the exact same plots as using mitoAnalysis.R
### and the data from the original ImageJ script

p1 <- ggplot(mito_summary_df, aes(x = condA:condB, y = rTPD54)) +
  geom_hline(yintercept = 1, linetype = "dashed", col = "grey") +
  geom_sina(aes(colour = condA, shape = condB)) +
  scale_colour_manual(values = c("#117733", "#999933", "#44aa99")) +
  scale_shape_manual(values = c(1, 16)) +
  stat_summary(geom = "point", fun = "mean", col = "black", size = 2) +
  ylim(c(0, NA)) +
  labs(x = "", y = "TPD54 (Post/Pre)") +
  theme_bw(10) +
  theme(legend.position = "none")

p2 <- ggplot(mito_summary_df, aes(x = condA:condB, y = rATG9)) +
  geom_hline(yintercept = 1, linetype = "dashed", col = "grey") +
  geom_sina(aes(colour = condA, shape = condB)) +
  scale_colour_manual(values = c("#117733", "#999933", "#44aa99")) +
  scale_shape_manual(values = c(1, 16)) +
  stat_summary(geom = "point", fun = "mean", col = "black", size = 2) +
  ylim(c(0, NA)) +
  labs(x = "", y = "mCherry (Post/Pre)") +
  theme_bw(10) +
  theme(legend.position = "none")

p3 <- ggplot(mito_summary_df, aes(x = rATG9, y = rTPD54)) +
  geom_hline(yintercept = 1, linetype = "dashed", col = "grey") +
  geom_vline(xintercept = 1, linetype = "dashed", col = "grey") +
  geom_point(aes(colour = condA, shape = condB), alpha = 0.5) +
  geom_smooth(method = lm,
              se = FALSE, linetype = "dashed",
              alpha = 0.75, aes(colour = condA)) +
  facet_wrap(. ~ condB, labeller = as_labeller(condition_names)) +
  scale_colour_manual(values = c("#117733", "#999933", "#44aa99")) +
  scale_shape_manual(values = c(1, 16)) +
  labs(x = "mCherry (Post/Pre)", y = "TPD54 (Post/Pre)") +
  lims(x = c(0.5, 3), y = c(0.5, 2.5)) +
  theme_bw(10) +
  coord_fixed() +
  theme(legend.position = "none")

p1 <- meaningfully_x_label(p1)
p2 <- meaningfully_x_label(p2)

ggsave("Output/Plots/tpd54_ratio.png", p1,
       width = 7, height = 4, units = "in", dpi = 300)
ggsave("Output/Plots/atg9_ratio.png", p2,
       width = 7, height = 4, units = "in", dpi = 300)
ggsave("Output/Plots/tpd54vatg9_ratio.png", p3,
       width = 8, height = 4, units = "in", dpi = 300)

# Now make super plot versions
# create the mean of ratios per experiment by grouping variables
sp_df <- group_by(mito_summary_df, condA, condB, expt)
sp_summary <- summarise(sp_df,
                        n = n(),
                        rTPD54 = mean(rTPD54),
                        rATG9 = mean(rATG9)
)

# superplot
p4 <- ggplot() +
  geom_hline(yintercept = 1, linetype = "dashed", col = "grey") +
  geom_sina(
    data = mito_summary_df,
    aes(x = condA:condB, y = rTPD54, colour = expt),
    alpha = 0.5, position = "auto", size = 0.8, maxwidth = 0.4,
    shape = 16
  ) +
  geom_point(
    data = sp_summary,
    aes(x = condA:condB, y = rTPD54, fill = expt),
    shape = 22, size = 1.5, stroke = 0.5, alpha = 0.7
  ) +
  scale_color_manual(values = c("#4477aa", "#ccbb44", "#ee6677", "#000000")) +
  scale_fill_manual(values = c("#4477aa", "#ccbb44", "#ee6677", "#000000")) +
  scale_y_continuous(
    breaks = c(0, 1, 2),
    minor_breaks = c(0.5, 1.5), limits = c(0, NA)
  ) +
  labs(x = "", y = expression(Mito ~ TPD54 ~ (F[post] / F[pre]))) +
  theme_bw(9) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

p5 <- ggplot() +
  geom_hline(yintercept = 1, linetype = "dashed", col = "grey") +
  geom_sina(
    data = mito_summary_df,
    aes(x = condA:condB, y = rATG9, colour = expt),
    alpha = 0.5, position = "auto", size = 0.8, maxwidth = 0.4,
    shape = 16
  ) +
  geom_point(
    data = sp_summary, aes(x = condA:condB, y = rATG9, fill = expt),
    shape = 22, size = 1.5, stroke = 0.5, alpha = 0.7
  ) +
  scale_color_manual(values = c("#4477aa", "#ccbb44", "#ee6677", "#000000")) +
  scale_fill_manual(values = c("#4477aa", "#ccbb44", "#ee6677", "#000000")) +
  ylim(c(0, NA)) +
  labs(x = "", y = expression(Mito ~ mCherry ~ (F[post] / F[pre]))) +
  theme_bw(9) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

p4 <- meaningfully_x_label(p4)
p5 <- meaningfully_x_label(p5)

ggsave("Output/Plots/tpd54_ratio_sp.pdf", p4,
       width = 86, height = 90, units = "mm"
)
ggsave("Output/Plots/atg9_ratio_sp.pdf", p5,
       width = 86, height = 90, units = "mm"
)


#
# Stats

TPD54_model <- aov(rTPD54 ~ condA * condB, data = sp_summary)
summary(TPD54_model)
TukeyHSD(TPD54_model, conf.level = 0.95)
