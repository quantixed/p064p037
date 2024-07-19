library(ggplot2)
library(ggforce)
library(dplyr)
library(tidyr)
library(reshape2)

# function to load in the data into a big data frame
# relies on files being called `allResults_MF046.csv` etc
# can be of form `foo_expid.csv` the code picks up `expid` from this filename
# uses column names as defined by Mary in her analysis script.
# this is the Mito & Cyto version of the ATG9 and TPD54 analysis.

# this script is for MF053 and subsequent repeats of fixed cell expts
# prepost is actually control and rapalog (because it's fixed cell)
# GFP-FKBP, GFP-FKBP-TPD54 and R159E are rerouted, ATG9A is stained
# so we have mito/cyto for each channel per cell,
# but for post/pre i.e. rapa/ctrl we can only do this for the population

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

# this function is probably unnecessary since we could probably do it with pivot
# here we take the mito / cyto ratio for each cell (each channel)
# option 0 is mito / cyto, option 1 is mito / (mito + cyto)
calculate_mitocyto_ratios <- function(df, opt) {
  # get the combinations that we want ratios for *includes expt*
  combinations <- group_by(df, Channel, Filename,
                           Condition, PrePost, CondA, expt) %>%
    summarise()
  # make column for ratios
  combinations$Ratio <- rep(0, nrow(combinations))

  for (i in seq_len(nrow(combinations))) {
    mitoSignal <- df$Mean.BG[df$Slice == 0 &
                               df$Channel == combinations$Channel[i] &
                               df$Filename == combinations$Filename[i] &
                               df$Condition == combinations$Condition[i] &
                               df$expt == combinations$expt[i]]
    cytoSignal <- df$Mean.BG[df$Slice == 1 &
                               df$Channel == combinations$Channel[i] &
                               df$Filename == combinations$Filename[i] &
                               df$Condition == combinations$Condition[i] &
                               df$expt == combinations$expt[i]]
    if (opt == 0) {
      combinations$Ratio[i] <- mitoSignal[1] / cytoSignal[1]
    } else {
      combinations$Ratio[i] <- mitoSignal[1] / (mitoSignal[1] + cytoSignal[1])
    }
  }
  combinations$CondA <- as.factor(combinations$CondA)
  combinations$CondA <- factor(combinations$CondA,
                               levels = c("Cont", "WT", "Mut"))
  combinations$Channel <- as.factor(combinations$Channel)

  return(combinations)
}

# function for labelling x-axis in a more meaningful way
meaningfully_x_label <- function(p) {
  # p is ggplot
  p <- p + scale_x_discrete(labels = c(
    "ATG9:Cont" = "ATG9A\nGFP-FKBP",
    "ATG9:WT" = "ATG9A\nGFP-FKBP-TPD54",
    "ATG9:Mut" = "ATG9A\nGFP-FKBP-TPD54(R159E)",
    "TPD54:Cont" = "GFP-FKBP",
    "TPD54:WT" = "GFP-FKBP-TPD54",
    "TPD54:Mut" = "GFP-FKBP-TPD54(R159E)"
  ))
  return(p)
}

## Variables ----

slice_names <- c(
  "0" = "Mitochondria",
  "1" = "Cytoplasm",
  "ATG9" = "ATG9A",
  "TPD54" = "GFP"
)

## Script ----

data <- process_data()
# Condition can have trailing /
data$Condition <- gsub("/$", "", data$Condition)
# PrePost extracted
data$PrePost <- sapply(strsplit(data$Condition, "_"), "[[", 1)
# extract conditions
data$CondA <- sapply(strsplit(data$Condition, "_"), "[[", 2)

# summary of data

# we need a wide version for the scatter plot
# we use opt 1 to get values for ratio on [0,1] scale
data_wide_0 <- calculate_mitocyto_ratios(data, opt = 0)
data_wide_1 <- calculate_mitocyto_ratios(data, opt = 1)

# we need get rid of cells where GFP mito / cyto was 0.5 or less
# equivalent is 1 or less of raw mito / cyto
# this is done from both groups (control and rapa)
# ATG9 value for that cell is also removed
data_wide_0 <- data_wide_0 %>%
  pivot_wider(names_from = Channel, values_from = Ratio) %>%
  filter(TPD54 > 1) %>%
  pivot_longer(cols = c("ATG9", "TPD54"),
               names_to = "Channel", values_to = "Ratio")
data_wide_1 <- data_wide_1 %>%
  pivot_wider(names_from = Channel, values_from = Ratio) %>%
  filter(TPD54 > 0.5) %>%
  pivot_longer(cols = c("ATG9", "TPD54"),
               names_to = "Channel", values_to = "Ratio")

# cyto data can be used to compare post as a group with pre as a group per expt
# cytoloss <- data %>%
#   filter(Slice == 1) %>%
#   group_by(Channel,expt,PrePost,CondA) %>%
#   summarise(mean = mean(Mean.BG)) %>%
#   pivot_wider(names_from = PrePost, values_from = mean)
# cytoloss$ratio <- cytoloss$postRapalog / cytoloss$preRapalog
# cytoloss$CondA <- factor(cytoloss$CondA, levels = c("Cont", "WT", "Mut"))

###########################################
## generate plots
## colour/shape scheme is
## shape: open control, closed treatment
## colours: WT, Mut, Cl35 = #117733", "#999933", "#44aa99"
## colours: Cont, WT, Mut = #88ccee", "#11773", "#999933"
## colours: superplot: #4477aa", "#ccbb44", "#ee6677"

p1 <- ggplot(data, aes(
  x = factor(CondA, levels = c("Cont", "WT", "Mut")),
  y = Mean.BG,
  shape = factor(PrePost, levels = c("preRapalog", "postRapalog")),
  colour = factor(CondA, levels = c("Cont", "WT", "Mut"))
)) +
  geom_sina(alpha = 0.75) +
  facet_wrap(Slice ~ Channel,
             labeller = as_labeller(slice_names), scales = "free") +
  scale_shape_manual(values = c(1, 16)) +
  scale_color_manual(values = c("#88ccee", "#117733", "#999933")) +
  labs(x = "", y = "Signal") +
  theme_bw(10) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c(
    "Cont" = "GFP-FKBP",
    "WT" = "GFP-FKBP-TPD54",
    "Mut" = "GFP-FKBP-TPD54(R159E)"
  ))

ggsave("Output/Plots/meanbg_raw.png", p1,
       width = 8, height = 6, units = "in", dpi = 300)

p2 <- ggplot(
  data_wide_0,
  aes(
    x = factor(CondA, levels = c("Cont", "WT", "Mut")),
    y = Ratio,
    shape = factor(PrePost, levels = c("preRapalog", "postRapalog")),
    colour = factor(CondA, levels = c("Cont", "WT", "Mut"))
  )
) +
  geom_hline(yintercept = 1, linetype = "dashed", col = "grey") +
  geom_sina(alpha = 0.75) +
  facet_wrap(. ~ Channel, labeller = as_labeller(slice_names)) +
  scale_shape_manual(values = c(1, 16)) +
  scale_color_manual(values = c("#88ccee", "#117733", "#999933")) +
  scale_y_continuous(limits = c(0.5, 64), trans = "log2") +
  labs(x = "", y = expression(Ratio ~ (F[mito] / F[cyto]))) +
  theme_bw(10) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c(
    "Cont" = "GFP-FKBP",
    "WT" = "GFP-FKBP-TPD54",
    "Mut" = "GFP-FKBP-TPD54(R159E)"
  ))

ggsave("Output/Plots/mitocytoratio_0.png", p2,
       width = 9, height = 5, units = "in", dpi = 300)

p3 <- ggplot(
  data_wide_1,
  aes(
    x = factor(CondA, levels = c("Cont", "WT", "Mut")),
    y = Ratio,
    shape = factor(PrePost, levels = c("preRapalog", "postRapalog")),
    colour = factor(CondA, levels = c("Cont", "WT", "Mut"))
  )
) +
  geom_hline(yintercept = 0.5, linetype = "dashed", col = "grey") +
  geom_sina(alpha = 0.75) +
  facet_wrap(. ~ Channel, labeller = as_labeller(slice_names)) +
  scale_shape_manual(values = c(1, 16)) +
  scale_color_manual(values = c("#88ccee", "#117733", "#999933")) +
  lims(y = c(0, 1)) +
  labs(x = "", y = expression(Mito ~ Fraction ~ (F[mito] / F[total]))) +
  theme_bw(10) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c(
    "Cont" = "GFP-FKBP",
    "WT" = "GFP-FKBP-TPD54",
    "Mut" = "GFP-FKBP-TPD54(R159E)"
  ))

ggsave("Output/Plots/mitocytoratio_1.png", p3,
       width = 9, height = 5, units = "in", dpi = 300)

p4 <- data_wide_0 %>%
  pivot_wider(
    names_from = "Channel",
    values_from = "Ratio"
  ) %>%
  ggplot(aes(
    x = TPD54,
    y = ATG9,
    shape = factor(PrePost, levels = c("preRapalog", "postRapalog")),
    colour = factor(CondA, levels = c("Cont", "WT", "Mut"))
  )) +
  geom_hline(yintercept = 1, linetype = "dashed", col = "grey") +
  geom_vline(xintercept = 1, linetype = "dashed", col = "grey") +
  geom_point(alpha = 0.5) +
  geom_smooth(method = lm, se = FALSE, linetype = "dashed", alpha = 0.75) +
  scale_shape_manual(values = c(1, 16)) +
  scale_color_manual(values = c("#88ccee", "#117733", "#999933")) +
  facet_wrap(. ~ factor(PrePost, levels = c("preRapalog", "postRapalog"),
                        labels = c("Control", "Rapalog"))) +
  scale_y_continuous(limits = c(0.125, 128), trans = "log2") +
  scale_x_continuous(limits = c(0.125, 128), trans = "log2") +
  labs(
    x = expression(GFP ~ Ratio ~ (F[mito] / F[cyto])),
    y = expression(ATG9A ~ Ratio ~ (F[mito] / F[cyto]))
  ) +
  theme_bw(10) +
  coord_fixed() +
  theme(legend.position = "none")

ggsave("Output/Plots/mitocytoratio_scatter_0.png", p4,
       width = 10, height = 8, units = "in", dpi = 300)

p5 <- data_wide_1 %>%
  pivot_wider(
    names_from = "Channel",
    values_from = "Ratio"
  ) %>%
  ggplot(aes(
    x = TPD54,
    y = ATG9,
    shape = factor(PrePost, levels = c("preRapalog", "postRapalog")),
    colour = factor(CondA, levels = c("Cont", "WT", "Mut"))
  )) +
  geom_hline(yintercept = 0.5, linetype = "dashed", col = "grey") +
  geom_vline(xintercept = 0.5, linetype = "dashed", col = "grey") +
  geom_point(alpha = 0.5) +
  geom_smooth(method = lm, se = FALSE, linetype = "dashed", alpha = 0.75) +
  scale_shape_manual(values = c(1, 16)) +
  scale_color_manual(values = c("#88ccee", "#117733", "#999933")) +
  facet_wrap(. ~ factor(PrePost, levels = c("preRapalog", "postRapalog"),
                        labels = c("Control", "Rapalog"))) +
  lims(x = c(0.4, 1), y = c(0.4, 1)) +
  labs(
    x = expression(GFP ~ Ratio ~ (F[mito] / F[total])),
    y = expression(ATG9A ~ Ratio ~ (F[mito] / F[total]))
  ) +
  theme_bw(10) +
  coord_fixed() +
  theme(legend.position = "none")

ggsave("Output/Plots/mitocytoratio_scatter_1.png", p5,
       width = 10, height = 8, units = "in", dpi = 300)

# Now make super plot versions

# create the mean of ratios per experiment by grouping variables
sp_summary_0 <- data_wide_0 %>%
  group_by(CondA, PrePost, Channel, expt) %>%
  summarise(
    n = n(),
    Ratio = mean(Ratio)
  )
sp_summary_1 <- data_wide_1 %>%
  group_by(CondA, PrePost, Channel, expt) %>%
  summarise(
    n = n(),
    Ratio = mean(Ratio)
  )

# superplot
p9 <- ggplot() +
  geom_hline(yintercept = 1, linetype = "dashed", col = "grey") +
  geom_sina(data = data_wide_0, aes(
    x = factor(CondA,
               levels = c("Cont", "WT", "Mut")):
      factor(PrePost, levels = c("preRapalog", "postRapalog")),
    y = Ratio,
    shape = factor(PrePost, levels = c("preRapalog", "postRapalog")),
    colour = expt
  ), alpha = 0.5, position = "auto") +
  geom_point(data = sp_summary_0, aes(
    x = factor(CondA,
               levels = c("Cont", "WT", "Mut")):
      factor(PrePost, levels = c("preRapalog", "postRapalog")),
    y = Ratio,
    fill = expt
  ), shape = 22, size = 2, stroke = 0.5, alpha = 0.7) +
  scale_shape_manual(values = c(1, 16)) +
  scale_color_manual(values = c("#4477aa", "#ccbb44", "#ee6677", "#000000")) +
  scale_fill_manual(values = c("#4477aa", "#ccbb44", "#ee6677", "#000000")) +
  facet_wrap(. ~ Channel, labeller = as_labeller(slice_names)) +
  scale_y_continuous(limits = c(0.5, 64), trans = "log2") +
  labs(x = "", y = expression(Ratio ~ (F[mito] / F[cyto]))) +
  theme_bw(9) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c(
    "Cont:preRapalog" = "-\nGFP-FKBP",
    "Cont:postRapalog" = "+\nGFP-FKBP",
    "WT:preRapalog" = "-\nGFP-FKBP-TPD54",
    "WT:postRapalog" = "+\nGFP-FKBP-TPD54",
    "Mut:preRapalog" = "-\nGFP-FKBP-TPD54(R159E)",
    "Mut:postRapalog" = "+\nGFP-FKBP-TPD54(R159E)"
  ))

ggsave("Output/Plots/mitocytoratio_sp_0.png", p9,
       width = 10, height = 6, units = "in", dpi = 300)

p10 <- ggplot() +
  geom_hline(yintercept = 0.5, linetype = "dashed", col = "grey") +
  geom_sina(data = data_wide_1, aes(
    x = factor(CondA, levels = c("Cont", "WT", "Mut")):
      factor(PrePost, levels = c("preRapalog", "postRapalog")),
    y = Ratio,
    shape = factor(PrePost, levels = c("preRapalog", "postRapalog")),
    colour = expt
  ), alpha = 0.5, position = "auto") +
  geom_point(data = sp_summary_1, aes(
    x = factor(CondA, levels = c("Cont", "WT", "Mut")):
      factor(PrePost, levels = c("preRapalog", "postRapalog")),
    y = Ratio,
    fill = expt
  ), shape = 22, size = 2, stroke = 0.5, alpha = 0.7) +
  scale_shape_manual(values = c(1, 16)) +
  scale_color_manual(values = c("#4477aa", "#ccbb44", "#ee6677", "#000000")) +
  scale_fill_manual(values = c("#4477aa", "#ccbb44", "#ee6677", "#000000")) +
  facet_wrap(. ~ Channel, labeller = as_labeller(slice_names)) +
  labs(x = "", y = expression(Ratio ~ (F[mito] / F[total]))) +
  theme_bw(9) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c(
    "Cont:preRapalog" = "-\nGFP-FKBP",
    "Cont:postRapalog" = "+\nGFP-FKBP",
    "WT:preRapalog" = "-\nGFP-FKBP-TPD54",
    "WT:postRapalog" = "+\nGFP-FKBP-TPD54",
    "Mut:preRapalog" = "-\nGFP-FKBP-TPD54(R159E)",
    "Mut:postRapalog" = "+\nGFP-FKBP-TPD54(R159E)"
  ))

ggsave("Output/Plots/mitocytoratio_sp_1.png", p10,
       width = 10, height = 6, units = "in", dpi = 300)


# this is to look at "raw" data - not saved
p11 <- ggplot(data, aes(
  x = factor(CondA, levels = c("Cont", "WT", "Mut")),
  y = Mean.BG,
  shape = factor(PrePost, levels = c("preRapalog", "postRapalog")),
  colour = factor(expt)
)) +
  geom_sina(alpha = 0.75) +
  facet_wrap(Slice ~ Channel, labeller = as_labeller(slice_names),
             scales = "free") +
  scale_shape_manual(values = c(1, 16)) +
  labs(x = "", y = "Signal") +
  theme_bw(10) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c(
    "Cont" = "GFP-FKBP",
    "WT" = "GFP-FKBP-TPD54",
    "Mut" = "GFP-FKBP-TPD54(R159E)"
  ))

cytoloss <- sp_summary_1 %>%
  select(-n) %>%
  pivot_wider(
    names_from = c("Channel", "PrePost"),
    values_from = "Ratio"
  )
cytoloss$ATG9_cytoloss <- ((1 - cytoloss$ATG9_postRapalog) -
                             (1 - cytoloss$ATG9_preRapalog)) /
  (1 - cytoloss$ATG9_preRapalog) * 100
cytoloss$TPD54_cytoloss <- ((1 - cytoloss$TPD54_postRapalog) -
                              (1 - cytoloss$TPD54_preRapalog)) /
  (1 - cytoloss$TPD54_preRapalog) * 100

p11 <- cytoloss %>%
  pivot_longer(cols = c("ATG9_cytoloss", "TPD54_cytoloss"),
               names_to = "Channel", values_to = "Loss") %>%
  ggplot() +
  geom_sina(aes(x = CondA, y = Loss, colour = expt), size = 3.5, alpha = 0.9) +
  scale_color_manual(values = c("#4477aa", "#ccbb44", "#ee6677", "#000000")) +
  facet_wrap(. ~ Channel) +
  labs(x = "", y = "Cytoplasmic Loss (%)") +
  theme_bw(10) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c(
    "Cont" = "GFP-FKBP",
    "WT" = "GFP-FKBP-TPD54",
    "Mut" = "GFP-FKBP-TPD54(R159E)"
  ))

ggsave("Output/Plots/cytoloss.png", p11,
       width = 8, height = 2.2, units = "in", dpi = 300)

ggplot() +
  geom_hline(yintercept = 0.5, linetype = "dashed", col = "grey") +
  geom_sina(data = data_wide_1 %>% filter(Channel == "ATG9"), aes(
    x = factor(CondA, levels = c("Cont", "WT", "Mut")):
      factor(PrePost, levels = c("preRapalog", "postRapalog")),
    y = Ratio,
    shape = factor(PrePost, levels = c("preRapalog", "postRapalog")),
    colour = expt
  ), alpha = 0.5, position = "auto") +
  geom_point(data = sp_summary_1 %>% filter(Channel == "ATG9"), aes(
    x = factor(CondA, levels = c("Cont", "WT", "Mut")):
      factor(PrePost, levels = c("preRapalog", "postRapalog")),
    y = Ratio,
    fill = expt
  ), shape = 22, size = 2, stroke = 0.5, alpha = 0.7) +
  scale_shape_manual(values = c(1, 16)) +
  scale_color_manual(values = c("#4477aa", "#ccbb44", "#ee6677", "#000000")) +
  scale_fill_manual(values = c("#4477aa", "#ccbb44", "#ee6677", "#000000")) +
  labs(x = "", y = expression(Ratio ~ (F[mito] / F[total]))) +
  theme_bw(11) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c(
    "Cont:preRapalog" = "-\nGFP-FKBP",
    "Cont:postRapalog" = "+\nGFP-FKBP",
    "WT:preRapalog" = "-\nGFP-FKBP-TPD54",
    "WT:postRapalog" = "+\nGFP-FKBP-TPD54",
    "Mut:preRapalog" = "-\nGFP-FKBP-TPD54(R159E)",
    "Mut:postRapalog" = "+\nGFP-FKBP-TPD54(R159E)"
  ))
ggsave("Output/Plots/mitocytoratio_sp_2.png", width = 12,
       height = 6, units = "in", dpi = 300)

#
# Stats
stat_summary <- sp_summary_0 %>%
  pivot_wider(names_from = Channel, values_from = Ratio)

TPD54_model <- aov(TPD54 ~ CondA * PrePost, data = stat_summary)
summary(TPD54_model)
TukeyHSD(TPD54_model, conf.level = 0.95)

ATG9_model <- aov(ATG9 ~ CondA * PrePost, data = stat_summary)
summary(ATG9_model)
TukeyHSD(ATG9_model, conf.level = 0.95)

## Figure ----

p12 <- ggplot() +
  geom_hline(yintercept = 0.5, linetype = "dashed", col = "grey") +
  geom_sina(
    data = data_wide_1,
    aes(
      x = factor(CondA,
                 levels = c("Cont", "WT", "Mut")
      ):
        factor(PrePost,
               levels = c("preRapalog", "postRapalog")
        ),
      y = Ratio,
      shape = factor(PrePost, levels = c("preRapalog", "postRapalog")),
      colour = expt
    ),
    alpha = 0.5, position = "auto", size = 0.8,
    maxwidth = 0.4
  ) +
  geom_point(
    data = sp_summary_1,
    aes(
      x = factor(CondA,
                 levels = c("Cont", "WT", "Mut")
      ):
        factor(PrePost, levels = c("preRapalog", "postRapalog")),
      y = Ratio, fill = expt
    ),
    shape = 22, size = 1.5, stroke = 0.5, alpha = 0.7
  ) +
  scale_shape_manual(values = c(1, 16)) +
  scale_color_manual(values = c("#4477aa", "#ccbb44", "#ee6677", "#000000")) +
  scale_fill_manual(values = c("#4477aa", "#ccbb44", "#ee6677", "#000000")) +
  facet_wrap(. ~ factor(Channel, levels = c("TPD54", "ATG9")),
             labeller = as_labeller(slice_names), nrow = 2) +
  labs(x = "", y = expression(Ratio ~ (F[mito] / F[total]))) +
  lims(y = c(0.5, 1)) +
  theme_bw(9) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c(
    "Cont:preRapalog" = "-\nGFP-FKBP",
    "Cont:postRapalog" = "+\nGFP-FKBP",
    "WT:preRapalog" = "-\nGFP-FKBP-TPD54",
    "WT:postRapalog" = "+\nGFP-FKBP-TPD54",
    "Mut:preRapalog" = "-\nGFP-FKBP-TPD54(R159E)",
    "Mut:postRapalog" = "+\nGFP-FKBP-TPD54(R159E)"
  ))

ggsave("Output/Plots/mitototalratio_sp.pdf", p12,
       width = 6.5, height = 10.5, units = "cm")
