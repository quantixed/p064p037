library(plyr)
library(ggplot2)
library(dplyr)

## Functions ----
# this function will separate values with semi colons and give them new rows, then sum duplicate entries
df_to_single_protein <- function(df, Genelist) {
  df$uid <- as.character(1:nrow(df))
  build <- data.frame()
  for(i in 1:length(df[,Genelist])) {
    brick <- data.frame()
    gene <- unlist(strsplit(df[i,Genelist], ";"))
    for(j in 1:length(gene)){
      brick <- rbind(brick, df[i, ])
      brick[j, Genelist] <- gene[j]
    } 
    build <- rbind(build,brick)
  }
  
  separated_with_uid <- build[, c(Genelist, "uid")]
  # remove rows with duplicate Genelist values (keep first)
  separated_with_uid <- separated_with_uid[!duplicated(separated_with_uid[,Genelist]), ]
  # sum duplicate entries in build
  build <- ddply(build, Genelist, numcolwise(sum))
  build <- merge(build, separated_with_uid, by = Genelist, all.x = TRUE)
  
  return(build)
}

## Script ----

# read in the original data
INV <- read.delim("Output/Data/combined_inv_proteome_with_intensities.txt")
SLMV <- read.delim("Data/SLMV_MassSpec.txt")
ATG <- read.delim("Data/ATG9a proteomics.txt")

# INV dataset cleaning.
INV <- INV[, c("so_SHORTNAME", "intensity")]

# SLMV dataset cleaning.
SLMV <- SLMV[, c("Gene.Names", "Peptides")]
names(SLMV)[names(SLMV) == "Gene.Names"] <- "Gene.names"
# Gene.names are separated by a space, separate them by a semicolon instead
SLMV$Gene.names <- gsub(" ", ";", SLMV$Gene.names)

# ATG9 dataset cleaning. Gets raw LFQ values from log10
ATG$LFQ <- sapply(ATG$Intensity_ES_vs_FM_ATG9A.Immunoisolation_intensity_filtered, function(df){ 10^df})
names(ATG)[names(ATG) == "Gene_names"] <- "Gene.names"
ATG$FC10 <- ATG$Intensity_ES_vs_FM_ATG9A.Immunoisolation_intensity_filtered - ATG$Intensity_Mix_CTRL.IgM.Immunoisolation
# if Intensity_Mix_CTRL.IgM.Immunoisolation is 0 then FC10 is 1
ATG$FC10[ATG$Intensity_Mix_CTRL.IgM.Immunoisolation == 0] <- 1
# drop rows with NA
ATG <- ATG[!is.na(ATG$FC10), ]
# drop rows with blank Gene.names
ATG <- ATG[ATG$Gene.names != "", ]
# drop rows with logFC_ES_vs_FM_ATG9A.Immunoisolation_intensity_filtered < -3
ATG <- ATG[ATG$logFC_ES_vs_FM_ATG9A.Immunoisolation_intensity_filtered > -3, ]
# drop rows with FC10 < 0.6
ATG <- ATG[ATG$FC10 > 0.6, ]

# before simplifying ATG, we will save the cleaned data
ATG_cleaned <- df_to_single_protein(ATG,"Gene.names")
write.csv(ATG_cleaned, "Output/Data/ATG_cleaned.txt", row.names = FALSE)

# simplify the datasets
names(INV) <- c("Gene.names", "INV")
names(SLMV) <- c("Gene.names", "SLMV")
ATG <- ATG %>% 
  select(Gene.names, LFQ)
names(ATG) <- c("Gene.names", "ATG")

# Separates values with semi colons and gives them new rows, then sums duplicate entries
# INV <- df_to_single_protein(INV,"Gene.names") # not needed for INV
SLMV <- df_to_single_protein(SLMV,"Gene.names")
ATG <- df_to_single_protein(ATG,"Gene.names")

write.csv(INV, "Output/Data/INV.txt", row.names = FALSE)
write.csv(SLMV, "Output/Data/SLMV.txt", row.names = FALSE)
write.csv(ATG, "Output/Data/ATG.txt", row.names = FALSE)

