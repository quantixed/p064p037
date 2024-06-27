library(ggplot2)
library(dplyr)
library(ggrepel)

# read in the cleaned data
ATG <- read.csv("Output/Data/ATG_cleaned.txt")

# load in the INV data
INV <- read.delim("Output/Data/combined_inv_proteome_with_intensities.txt")
INV <- INV[, c("so_SHORTNAME", "so_ratioWave")]
names(INV) <- c("Gene.names", "INV")

# match if a gene is in the INV dataset
ATG$cat <- ifelse(ATG$Gene.names %in% INV$Gene.names, "INV", "ATG")
# get list of uids that are in the INV dataset
# get a list of INV matches - this is needed because after collapse of ATG, we might not have the INV names anymore
INV_matches <- data.frame(labelnames = ATG$Gene.names[ATG$cat == "INV"],
                          uid = ATG$uid[ATG$cat == "INV"])
# collapse back to a single uid, discard duplicates
ATG <- ATG %>% distinct(uid, .keep_all = TRUE)
# relabel. If uid is in the INV_uids list, set cat to INV
ATG$cat <- ifelse(ATG$uid %in% INV_matches$uid, "INV", "ATG")
# merge to add column which has teh labename of any uid that is in the INV dataset
ATG <- merge(ATG, INV_matches, by = "uid", all.x = TRUE)

# crap <- read.csv("Data/1719326562_gp.txt", sep = "\t")
# # filter for NA in Ave.SC column
# crap <- crap[!is.na(crap$Ave.SC), ]
# # select the Mapped.Gene.Symbol and Ave.SC columns
# crap <- crap[, c("Mapped.Gene.Symbol", "Ave.SC")]
# # merge the crap data with the ATG data
# ATG <- merge(ATG, crap, by.x = "Gene.names", by.y = "Mapped.Gene.Symbol", all.x = TRUE)

ATG %>% 
  # filter(Ave.SC < 5) %>% 
  ggplot(aes(x = logFC_ES_vs_FM_ATG9A.Immunoisolation_intensity_filtered, y = LFQ, shape = cat, colour = cat)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  geom_point(size = 1) +
  scale_color_manual(values = c("ATG" = "#7f7f7f7f", "INV" = "#ab66f0")) +
  scale_shape_manual(values = c("ATG" = 1, "INV" = 16)) +
  geom_text_repel(data = subset(ATG, cat == "INV"), aes(label = labelnames), size = 1.5, max.overlaps = 10, box.padding = 0.2) +
  scale_y_log10(limits = c(1e6,1e10)) +
  lims(x = c(-2.5, 2.5)) +
  labs(x = "ES - FM (Log2)", y = "Intensity") +
  theme_classic(9) +
  theme(legend.position = "none")

ggsave("Output/Plots/ATG_INV_label.pdf", width = 8, height = 8, bg = "white", units = "cm")
