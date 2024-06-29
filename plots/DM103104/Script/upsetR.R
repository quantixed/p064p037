library(dplyr)
library(UpSetR)
library(eulerr)

## Functions ----

make_shortnames_unique <- function(x) {
  y <- x

  for (i in seq_along(y)) {
    z <- unlist(strsplit(y[i], ";"))
    if (length(z) == 1) {
      next
    }
    for (j in seq_along(z)) {
      if (length(unique(y)) < length(unique(append(y, z[j])))) {
        y[i] <- z[j]
        break
      }
    }
  }

  return(y)
}

classify_upset_groups <- function(listIn) {
  all2 <- data.frame(protein = unique(unlist(listIn)))
  all1 <- lapply(listIn, function(x) {
    data.frame(protein = x)
  }) %>%
    bind_rows(.id = "path")

  all_int <- lapply(all2$protein, function(x) {
    # pull the name of the intersections
    intersection <- all1 %>%
      dplyr::filter(protein == x) %>%
      arrange(path) %>%
      pull("path") %>%
      paste0(collapse = "|")

    # build the dataframe
    data.frame(protein = x, int = intersection)
  }) %>%
    bind_rows()

  # all_int %>%
  #   group_by(int) %>%
  #   summarise(n = n()) %>%
  #   arrange(desc(n))

  return(all_int)
}

make_upset_table <- function(x) {
  y <- x %>%
    group_by(int) %>%
    summarise(n = n()) %>%
    arrange(desc(n))
  return(y)
}

# script ----
dm_df <- read.delim("Data/rankTable_WTvsControl.txt")
pe_df <- read.delim("Data/rankTable_INVvsControl.txt")
re_df <- read.delim("Data/rankTable_R159EvsControl.txt")

dm_df$so_SHORTNAME <- make_shortnames_unique(dm_df$so_SHORTNAME)
pe_df$so_SHORTNAME <- make_shortnames_unique(pe_df$so_SHORTNAME)
re_df$so_SHORTNAME <- make_shortnames_unique(re_df$so_SHORTNAME)

dm_df <- dm_df[dm_df$so_colorWave == 3, ]
pe_df <- pe_df[pe_df$so_colorWave == 3, ]
re_df <- re_df[re_df$so_colorWave == 3, ]

list_all <- list(
  DM = dm_df$so_SHORTNAME,
  PE = pe_df$so_SHORTNAME,
  RE = re_df$so_SHORTNAME
)

upset(fromList(list_all), order.by = "freq")

all_int_all <- classify_upset_groups(list_all)

# look at any proteins that were enriched in R159E
all_int_all[grepl("RE", all_int_all$int), ]

# we need to construct a table which shows each protein in this dataset
# name, shortname, uniprot ID, fold change, p-value
# precedence should be DM, then PE, then RE for the data
# the order should be by fold change,
# RE proteins at the end, with TPD52L2 exempted

# remove RE only
all_int_all <- all_int_all[all_int_all$int != "RE", ]
# first combine with dm_df
combine_df <- merge(all_int_all, dm_df,
                    by.x = "protein", by.y = "so_SHORTNAME",
                    all.x = TRUE, sort = FALSE)
# we now have complete rows and a bunch of NA rows that weren't matched,
# take complete half
df1 <- combine_df[!is.na(combine_df$so_ratioWave), ]
# missing values are:
missing_df <- all_int_all[is.na(combine_df$so_ratioWave), ]
df2 <- merge(missing_df, pe_df,
             by.x = "protein", by.y = "so_SHORTNAME",
             all.x = TRUE, sort = FALSE)
combine_df <- rbind(df1, df2)
# clean up other columns
combine_df$so_NAME <- make_shortnames_unique(combine_df$so_NAME)
combine_df$so_PID <- make_shortnames_unique(combine_df$so_PID)
# order by fold change
combine_df <- combine_df[order(combine_df$so_ratioWave, decreasing = TRUE), ]
# we currently have these columns
# keep <- c("protein","int","so_NAME","so_PID","so_productWave",
#           "so_colorWave","so_allTWave","so_ratioWave","so_keyW")
# place in desired order
keep <- c("protein", "so_NAME", "so_PID", "so_ratioWave", "so_allTWave", "int")
combine_df <- combine_df[, keep]
rownames(combine_df) <- 1:nrow(combine_df)
names(combine_df) <- c("so_SHORTNAME", "so_NAME", "so_PID",
                       "so_ratioWave", "so_allTWave", "int")
# save this first
write.table(combine_df, "Output/Data/combined_proteome.txt", sep = "\t", row.names = FALSE)

# now format for nice output
combine_df <- cbind(1:nrow(combine_df), combine_df)
names(combine_df) <- c("Rank", "Gene name", "Protein name", "Protein ID",
                       "Fold change", "P value", "Dataset")
# Intersection is things like DM|PE we want WT + INV
combine_df$Dataset <- gsub("|", " + ", combine_df$Dataset, fixed = TRUE)
combine_df$Dataset <- gsub("DM", "WT", combine_df$Dataset)
combine_df$Dataset <- gsub("RE", "R159E", combine_df$Dataset)
combine_df$Dataset <- gsub("PE", "INV", combine_df$Dataset)
# we now have combine_df as our final output for presentation as a Table
write.csv(format(combine_df, digits = 2),
          "Output/Data/combined_proteome.csv", row.names = FALSE)

# we require the background for comparisons - we do not take RE dataset
# because it is not background for complete_df
dm_df <- read.delim("Data/rankTable_WTvsControl.txt")
pe_df <- read.delim("Data/rankTable_INVvsControl.txt")
dm_df$so_SHORTNAME <- make_shortnames_unique(dm_df$so_SHORTNAME)
pe_df$so_SHORTNAME <- make_shortnames_unique(pe_df$so_SHORTNAME)
# clean up other columns
dm_df$so_NAME <- make_shortnames_unique(dm_df$so_NAME)
dm_df$so_PID <- make_shortnames_unique(dm_df$so_PID)
pe_df$so_NAME <- make_shortnames_unique(pe_df$so_NAME)
pe_df$so_PID <- make_shortnames_unique(pe_df$so_PID)
# we'll combine these
universe_df <- rbind(dm_df, pe_df)
# place in desired order
keep <- c("so_SHORTNAME", "so_NAME", "so_PID", "so_ratioWave", "so_allTWave")
universe_df <- universe_df[, keep]
universe_df <- universe_df[!duplicated(universe_df$so_SHORTNAME), ]
rownames(universe_df) <- 1:nrow(universe_df)
# save this first
write.table(universe_df, "Output/Data/universe_proteome.txt",
            sep = "\t", row.names = FALSE)
