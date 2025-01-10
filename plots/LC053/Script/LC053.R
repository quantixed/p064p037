library(dplyr)
library(tidyr)
library(ggplot2)
library(ggforce)
library(cowplot)

condition_a_names <- c(
  "Fed_DMSO" = "Fed + DMSO",
  "Fed_Baf" = "Fed + Baf",
  "Starved_DMSO" = "Starved + DMSO",
  "Starved_Baf" = "Starved + Baf"
)

condition_b_names <- c(
  "Control" = "siCtrl",
  "TPD54" = "siTPD54"
)

factorise_condition <- function(v, condition_names) {
  # v is a vector of conditions
  # condition_names is a named vector of condition names
  # returns a factorised vector
  v <- factor(v, levels = names(condition_names))
  levels(v) <- condition_names
  return(v)
}


# load in the data & fix labels
data <- read.csv("Data/Quantified_LC3.csv")
data$Starve <- sapply(strsplit(data$Condition, "_"), "[", 1)
data$Baf <- sapply(strsplit(data$Condition, "_"), "[", 2)
data$CondA <- paste(data$Starve, data$Baf, sep = "_")

# get ratios of LC3:Vinculin
ratio <- data %>% 
  group_by(Date, CondA, siRNA) %>%
  mutate(Ratio = Intensity[Protein=="LC3-II"]/Intensity[Protein=="Vinculin"])
ratio <- ratio[-which(ratio$Protein == "Vinculin"), ]  

# convert Date to character
ratio$Date <- as.character(ratio$Date) 

# factorise conditions
ratio$CondA <- factorise_condition(ratio$CondA, condition_a_names)
ratio$siRNA <- factorise_condition(ratio$siRNA, condition_b_names)


#Find mean & standard error of each conditions
mean_data <- ratio %>%
  group_by(Condition) %>%
  summarise(mean_Ratio = mean(Ratio),
            SE = sd(Ratio) / sqrt(length(Ratio)))

mean_data$Condition <- as.character(mean_data$Condition)
mean_data$Starve <- sapply(strsplit(mean_data$Condition, "_"), "[", 1)
mean_data$Baf <- sapply(strsplit(mean_data$Condition, "_"), "[", 2)
mean_data$siRNA <- sapply(strsplit(mean_data$Condition,"_"), "[",3)
mean_data$CondA <- paste(mean_data$Starve, mean_data$Baf, sep = "_")
# factorise conditions
mean_data$CondA <- factorise_condition(mean_data$CondA, condition_a_names)
mean_data$siRNA <- factorise_condition(mean_data$siRNA, condition_b_names)

#Re-order conditions
mean_data$CondA <- factor(mean_data$CondA, levels = c("Fed + DMSO", "Fed + Baf", 
                                                "Starved + DMSO", "Starved + Baf"))

ggplot(data = mean_data, aes(x = siRNA:CondA, y = mean_Ratio)) +
  geom_point(data = ratio,
             aes(x= siRNA:CondA, y=Ratio, colour = Date, shape = siRNA)) +
  scale_shape_manual(values = c("siCtrl" = 1, "siTPD54" = 16)) + # Custom Shapes
  scale_colour_manual(values = c("#4477aa", "#ccbb44", "#ee6677")) +
  geom_errorbar(aes(ymin = mean_Ratio - SE, ymax = mean_Ratio + SE), width = 0, 
                linewidth = 0.2) +
  geom_crossbar(aes(ymin = mean_Ratio, ymax = mean_Ratio), width = 0.8,
                linewidth = 0.4) +
  ylim(c(0, NA)) +  # Adjust the y-axis limits if needed
  labs(x = "", y = "LC3-II:Vinculin") +  # Adjust axis labels
  theme_cowplot(9) +
  theme(legend.position = "none",
        axis.text.x=element_blank())
ggsave("Output/Plots/LC3.pdf", width = 58, height = 38,
       units = "mm", bg = "white")

# Normalise to control
# to do this we will normalise all the siCtrl ratios to the mean_Ratio for Fed_DMSO_Control
# and normalise all the siTPD54 ratios to the mean_Ratio for Fed_DMSO_TPD54

# get the mean_Ratio for Fed_DMSO_Control
control_ratio <- mean_data$mean_Ratio[mean_data$CondA == "Fed + DMSO" & mean_data$siRNA == "siCtrl"]
# get the mean_Ratio for Fed_DMSO_TPD54
tpd54_ratio <- mean_data$mean_Ratio[mean_data$CondA == "Fed + DMSO" & mean_data$siRNA == "siTPD54"]
# normalise the ratios
ratio$Normalised_Ratio <- ratio$Ratio
ratio$Normalised_Ratio[ratio$siRNA == "siCtrl"] <- ratio$Ratio[ratio$siRNA == "siCtrl"] / control_ratio
ratio$Normalised_Ratio[ratio$siRNA == "siTPD54"] <- ratio$Ratio[ratio$siRNA == "siTPD54"] / tpd54_ratio

mean_norm_data <- ratio %>%
  group_by(Condition) %>%
  summarise(mean_norm_ratio = mean(Normalised_Ratio),
            SE = sd(Normalised_Ratio) / sqrt(length(Normalised_Ratio)))
mean_data$mean_norm_ratio <- mean_norm_data$mean_norm_ratio
mean_data$norm_se <- mean_norm_data$SE

# remove rows corresponding to Baf in Baf column
mean_data_short <- mean_data[!grepl("Baf", mean_data$Baf),]
ratio_short <- ratio[!grepl("Baf", ratio$Baf),]

# plot the normalised data
ggplot(data = mean_data_short, aes(x = siRNA:CondA, y = mean_norm_ratio)) +
  geom_point(data = ratio_short, aes(y = Normalised_Ratio, colour = Date, shape = siRNA)) +
  scale_shape_manual(values = c("siCtrl" = 1, "siTPD54" = 16)) + # Custom Shapes
  scale_colour_manual(values = c("#4477aa", "#ccbb44", "#ee6677")) +
  geom_errorbar(aes(ymin = mean_norm_ratio - norm_se, ymax = mean_norm_ratio + norm_se), width = 0, 
                linewidth = 0.2) +
  geom_crossbar(data  = mean_data_short, aes(ymin = mean_norm_ratio, ymax = mean_norm_ratio), width = 0.8,
                linewidth = 0.4) +
  ylim(c(0, NA)) +  # Adjust the y-axis limits if needed
  labs(x = "", y = "LC3-II (Normalized)") +  # Adjust axis labels
  theme_cowplot(9) +
  theme(legend.position = "none",
        axis.text.x=element_blank())
ggsave("Output/Plots/LC3_norm.pdf", width = 32, height = 38,
       units = "mm", bg = "white")
