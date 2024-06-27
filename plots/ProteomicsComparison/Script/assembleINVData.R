# script to combine both datasets into one data frame

# functions ----
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

# script ----

# load the list of enriched proteins to allow us to select
# the proteins of interest
# this is the Igor-style output from DM103104
enriched_proteins <- read.delim("Data/combined_proteome.txt", header = TRUE)
# now load the data (not ranked table) output from Igor
# for both DM and PE experiments
dm_df <- read.delim("Data/withMeans_WTvsControl.txt", header = TRUE)
pe_df <- read.delim("Data/withMeans_INVvsControl.txt", header = TRUE)
# maybe we don't need the combined data frame after all
# order by descending productWave and
# then make sure shortnames and PIDs are unique
dm_df <- dm_df[order(dm_df$productWave, decreasing = TRUE), ]
pe_df <- pe_df[order(pe_df$productWave, decreasing = TRUE), ]
dm_df$SHORTNAME <- make_shortnames_unique(dm_df$SHORTNAME)
pe_df$HORTNAME <- make_shortnames_unique(pe_df$SHORTNAME)
dm_df$PID <- make_shortnames_unique(dm_df$PID)
pe_df$PID <- make_shortnames_unique(pe_df$PID)
# we will also filter here
dm_df <- dm_df[dm_df$colorWave == 3, ]
pe_df <- pe_df[pe_df$colorWave == 3, ]

# now we need to put the LFQ intensity data back in place with the original data
temp <-  merge(enriched_proteins, dm_df,
               by.x = "so_PID", by.y = "PID", all.x = TRUE, sort = FALSE)
# add the meanCond1 column to original enriched_proteins data frame
enriched_proteins$dm_int <- 2^temp$meanCond1
temp <-  merge(enriched_proteins, pe_df,
               by.x = "so_PID", by.y = "PID", all.x = TRUE, sort = FALSE)
enriched_proteins$pe_int <- 2^temp$meanCond1
# this is not perfect due to intensity differences between the two experiments
enriched_proteins$intensity <- ifelse(is.na(enriched_proteins$dm_int),
                                      enriched_proteins$pe_int,
                                      ifelse(enriched_proteins$dm_int >
                                               enriched_proteins$pe_int,
                                             enriched_proteins$dm_int,
                                             enriched_proteins$pe_int))
# now save the data frame to Output/Data/
write.table(enriched_proteins,
            file = "Output/Data/combined_inv_proteome_with_intensities.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
