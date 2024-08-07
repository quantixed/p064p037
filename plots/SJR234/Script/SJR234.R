library(ggplot2)
library(ggforce)
library(dplyr)
library(tidyr)
library(cowplot)

# function to load in the data into a big data frame
# uses column names as defined by Mary in her analysis script.
# condition A is the cargo to be tested
# condition B is the FKBP tagged protein that is rerouted

## Static variables ----
# these are the conditions that we want to compare
condition_b_names <- c(
  "emptyv" = "mCherry",
  "atg9" = "ATG9A",
  "tpdwt" = "WT",
  "tpdmut" = "R159E"
)

condition_a_names <- c(
  "bif1" = "SH3GLB1",
  "arfip2" = "ARFIP2",
  "daglb" = "DAGLB",
  "pi4k2a" = "PI4K2A",
  "pi4k3b" = "PI4KB"
)


## Functions ----
# function to load in the data into a big data frame
process_data <- function() {
  all_files <- list.files("Data", pattern = "*.csv", full.names = TRUE)
  df_all <- data.frame()
  for (filename in all_files) {
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

factorise_condition <- function(v, condition_names) {
  # v is a vector of conditions
  # condition_names is a named vector of condition names
  # returns a factorised vector
  v <- factor(v, levels = names(condition_names))
  levels(v) <- condition_names
  return(v)
}

# function for labelling x-axis in a more meaningful way
meaningfully_x_label <- function(p) {
  # p is ggplot
  # make a named vector using condition_a_names and condition_b_names
  # repeat condition_a_names twice

  left <- paste0(rep(condition_a_names, times = 1,
                     each = length(condition_b_names)),
                 ":", rep(condition_b_names, length(condition_a_names)))
  right <- paste0(rep(condition_a_names, times = 1,
                      each = length(condition_b_names)),
                  "\n", rep(condition_b_names, length(condition_a_names)))
  named_vector <- setNames(right, left)

  p <- p + scale_x_discrete(labels = named_vector)
  return(p)
}

print_stats <- function(df, comb1, comb2) {
  query <- paste0(comb1, "-", comb2)
  pval <- stats[query, "p adj"]
  if (is.na(pval)) {
    query <- paste0(comb2, "-", comb1)
    pval <- stats[query, "p adj"]
  }
  return(pval)
}

## Script ----
# load in the data
data <- process_data()
# Condition can have trailing /
data$Condition <- gsub("/$", "", data$Condition)
# PrePost extracted
data$PrePost <- sapply(strsplit(data$Condition, "_"), "[[", 1)
# extract conditions
data$CondA <- tolower(sapply(strsplit(data$Condition, "_"), "[[", 3))
data$CondB <- tolower(sapply(strsplit(data$Condition, "_"), "[[", 2))

# generate id with Condition(pre/post) and Channel left as variables
# data$id <- paste(data$expt, data$Condition, data$File, data$cell, sep = "_")
data$id <- paste(data$expt, data$CondA, data$CondB, data$File, sep = "_x_")

# now we need to pivot the data and take the ratios
wide <- data %>%
  select(!c(Mean, BG, Area)) %>%
  pivot_wider(
    id_cols = c(id, Slice),
    names_from = c(Channel, PrePost),
    values_from = Mean.BG
  )
# now we need to calculate the ratios
wide$rFKBPTagPOI <- wide$FKBPTagPOI_postRapalog / wide$FKBPTagPOI_preRapalog
wide$rATG9Cargo <- wide$ATG9Cargo_postRapalog / wide$ATG9Cargo_preRapalog

# here we will remove the non-responders, i.e. the ones where rFKBPTagPOI
# is less than 1 but do it per cell.
# So if one cell is a non-responder then we will remove it from the analysis
wide <- wide %>%
  select(-c(ATG9Cargo_postRapalog,
            ATG9Cargo_preRapalog,
            FKBPTagPOI_postRapalog,
            FKBPTagPOI_preRapalog)) %>%
  pivot_wider(
    id_cols = id,
    names_from = Slice,
    values_from = c(rFKBPTagPOI, rATG9Cargo)
  ) %>%
  filter(rFKBPTagPOI_0 > 1) %>%
  # pivot_longer so that we have two columns of data
  # rFKBPTagPOI and rATG9Cargo with a column called Slice
  pivot_longer(
    cols = c(rFKBPTagPOI_0, rFKBPTagPOI_1, rATG9Cargo_0, rATG9Cargo_1),
    names_to = "Slice",
    values_to = "ratio"
  ) %>%
  # pivot_wider so that we have two columns of data
  # rFKBPTagPOI and rATG9Cargo with a column called Slice
  separate(col = Slice, c("temp2", "Slice")) %>%
  pivot_wider(
    names_from = temp2,
    values_from = ratio
  )


# pull out the condition names
wide$expt <- tolower(sapply(strsplit(wide$id, "_x_"), "[[", 1))
wide$condA <- sapply(strsplit(wide$id, "_x_"), "[[", 2)
wide$condB <- sapply(strsplit(wide$id, "_x_"), "[[", 3)
wide$condA <- factorise_condition(wide$condA, condition_a_names)
wide$condB <- factorise_condition(wide$condB, condition_b_names)
# we need to figure out which replicate is which because we may have three
# experimental replicates for each condition but they may be spread over many
# different experiment identifiers so that we can't just use the expt column
# if we assume that condA, exp combos have all the conditions in them then we
# can use that to figure out which replicate is which
reps <- wide %>%
  group_by(condA, expt) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(condA) %>%
  mutate(rep = as.character(seq_len(n()))) %>%
  ungroup() %>%
  select(condA, rep, expt)
wide <- left_join(wide, reps, by = c("condA", "expt"))

# get rid of any rows with NA in them
wide <- wide[complete.cases(wide), ]

# plots ----

# Now make super plot versions
# create the mean of ratios per experiment by grouping variables
sp_mito_summary <- wide %>%
  filter(Slice == 0) %>%
  group_by(condA, condB, rep) %>%
  summarise(
    n = n(),
    rFKBPTagPOI = mean(rFKBPTagPOI),
    rATG9Cargo = mean(rATG9Cargo)
  )

p11 <- ggplot() +
  geom_hline(yintercept = 1, linetype = "dashed", col = "grey") +
  geom_sina(
    data = wide %>% filter(Slice == 0),
    aes(x = condB, y = rATG9Cargo, colour = rep),
    alpha = 0.5, position = "auto", size = 0.8, maxwidth = 0.4,
    shape = 16
  ) +
  geom_point(data = sp_mito_summary,
             aes(x = condB, y = rATG9Cargo, fill = rep),
             shape = 22, size = 1.5, stroke = 0.5, alpha = 0.7) +
  scale_color_manual(values = c("#4477aa", "#ccbb44", "#ee6677", "#000000")) +
  scale_fill_manual(values = c("#4477aa", "#ccbb44", "#ee6677", "#000000")) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "", y = expression(Mitochondrial ~ GFP ~ (F[post] / F[pre]))) +
  facet_wrap(. ~ condA, ncol = 1) +
  theme_cowplot(9) +
  theme(legend.position = "none")

ggsave("Output/Plots/cargo_mitoratio_sp.pdf", p11,
       width = 51, height = 151, units = "mm")


# Stats ----

cat("Summary statistics\nCondA is cargo\nCondB is tag\n\n")

fkbp_model <- aov(rFKBPTagPOI ~ condA * condB, data = sp_mito_summary)
summary(fkbp_model)
TukeyHSD(fkbp_model, conf.level = 0.95)

cargo_model <- aov(rATG9Cargo ~ condA * condB, data = sp_mito_summary)
summary(cargo_model)
stats <- TukeyHSD(cargo_model, conf.level = 0.95)
stats <- as.data.frame(stats$`condA:condB`)
# print the p adj values of interest
print_stats(stats, "SH3GLB1:mCherry", "SH3GLB1:ATG9A")
print_stats(stats, "ARFIP2:mCherry", "ARFIP2:ATG9A")
print_stats(stats, "DAGLB:mCherry", "DAGLB:ATG9A")
print_stats(stats, "PI4K2A:mCherry", "PI4K2A:ATG9A")
print_stats(stats, "PI4KB:mCherry", "PI4KB:ATG9A")
print_stats(stats, "SH3GLB1:R159E", "SH3GLB1:WT")
print_stats(stats, "ARFIP2:R159E", "ARFIP2:WT")
print_stats(stats, "DAGLB:R159E", "DAGLB:WT")
print_stats(stats, "PI4K2A:R159E", "PI4K2A:WT")
print_stats(stats, "PI4KB:R159E", "PI4KB:WT")

