# this script works well and does what we want
# uses v18, and we are missing protein classifications

# if (!requireNamespace("BiocManager")) install.packages("BiocManager")
# BiocManager::install("PANTHER.db")
# if (!requireNamespace("AnnotationHub")) BiocManager::install("AnnotationHub")
# library(AnnotationHub)
# ah <- AnnotationHub()
# query(ah, "PANTHER.db")[[1]]

library(PANTHER.db)
library(treemap)

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

build_classifications <- function(df) {
  uids <- unique(df$UNIPROT)
  parent_df <- data.frame()
  for (uid in uids) {
    # extract rows that match uid
    subdf <- subset(df, UNIPROT == uid)
    if (nrow(subdf) == 1) {
      result <- pull_solo_ancestors(uid = uid,
                                    id = subdf$CLASS_ID[1],
                                    term = subdf$CLASS_TERM[1])
    } else {
      rightrow <- determine_correct_row(subdf)
      result <- pull_solo_ancestors(uid = uid,
                                    id = subdf$CLASS_ID[rightrow],
                                    term = subdf$CLASS_TERM[rightrow])
    }
    parent_df <- rbind(parent_df, result)
  }

  return(parent_df)
}

determine_correct_row <- function(df) {
  best <- 0L
  # we have 2 or more rows
  for (i in seq_len(nrow(df))) {
    # we want the row with child i.e. two ancestors (top and parent)
    ancestors <- traverseClassTree(PANTHER.db, df$CLASS_ID[i],
                                   scope = "ANCESTOR")
    if (length(ancestors) == 2) {
      return(i)
    }
    if (length(ancestors) == 1) {
      best <- i
    }
    # we do not deal with the possibility that we only find
    # three or more ancestors for all rows
    # we just return first row in that case
  }
  best <- max(1L, best)
  return(best)
}

pull_solo_ancestors <- function(uid, id, term) {
  if (is.na(id) || nchar(id) == 0) {
    # deal with a blank protein class id
    df <- data.frame(UNIPROT = uid,
                     class_id = id,
                     CLASS_ID = id,
                     CLASS_TERM = term,
                     UCLASS_TERM = NA,
                     UCLASS_ID = NA)
    return(df)
  }
  ancestors <- traverseClassTree(PANTHER.db, id, scope = "ANCESTOR")

  if (length(ancestors) > 2) {
    classterm <- select(PANTHER.db, ancestors[3], "CLASS_TERM", "CLASS_ID")
    uclassterm <- select(PANTHER.db, ancestors[2], "CLASS_TERM", "CLASS_ID")
    df <- data.frame(UNIPROT = uid,
                     class_id = id,
                     CLASS_ID = ancestors[2],
                     CLASS_TERM = classterm$CLASS_TERM[1],
                     UCLASS_ID = ancestors[3],
                     UCLASS_TERM = uclassterm$CLASS_TERM[1])
  }
  # next we deal with the following cases, but term may be blank
  if (is.na(term) || nchar(term) == 0) {
    orig <- select(PANTHER.db, id, "CLASS_TERM", "CLASS_ID")
    term <- orig$CLASS_TERM[1]
  }
  if (length(ancestors) == 1) {
    df <- data.frame(UNIPROT = uid,
                     class_id = id,
                     CLASS_ID = id,
                     CLASS_TERM = term,
                     UCLASS_ID = id,
                     UCLASS_TERM = term)
  }
  if (length(ancestors) == 2) {
    details <- select(PANTHER.db, ancestors, "CLASS_TERM", "CLASS_ID")
    df <- data.frame(UNIPROT = uid,
                     class_id = id,
                     CLASS_ID = id,
                     CLASS_TERM = term,
                     UCLASS_ID = details$CLASS_ID[1],
                     UCLASS_TERM = details$CLASS_TERM[1])
  }

  return(df)
}

process_and_make_treemap <- function(path, label = "", en = NULL) {
  if (!is.null(path)) {
    mydata <- read.delim(path)
    en <- mydata[mydata$so_colorWave == 3, ]
    # we need a unique list of protein IDs
    en$so_PID <- make_shortnames_unique(en$so_PID)
    en_up <- unique(en$so_PID)
  } else {
    if (is.null(en)) {
      cat("Need unique list of protein IDs!\n")
      return(NULL)
    }
    # we need a unique list of protein IDs
    en$so_PID <- make_shortnames_unique(en$so_PID)
    en_up <- unique(en$so_PID)
  }

  # retrieve all class id/terms, and the blanks
  mylist <- select(PANTHER.db, en_up, columns = c("CLASS_ID", "CLASS_TERM"),
                   keytype = "UNIPROT", jointype = "left")
  # merge with shortname so that we can patch obviously missing classifications
  left <- subset(en, select = c(so_PID, so_SHORTNAME))
  all_df <- merge(mylist, left, by.x = "UNIPROT", by.y = "so_PID",
                  all = TRUE, sort = FALSE)
  # patch some missing classes
  patch <- data.frame(so_SHORTNAME = c("TLN1", "VCL", "^FLOT",
                                       "^TPD5", "VAT1", "^SAR1",
                                       "^RAB", "^SSR", "^RTN",
                                       "L1CAM",  "NRP1", "COG5",
                                       "FYCO1", "ERLIN1", "SLC30A1",
                                       "ARHGEF2", "GNAS", "SQSTM1",
                                       "LAMTOR1", "SFT2D1", "ARL8A",
                                       "CKAP4", "WDFY1", "ATG9A",
                                       "SLC3A2", "LMBRD2", "AGTRAP"),
                      CLASS_ID = c("PC00041", "PC00041", "PC00151",
                                   "PC00151", "PC00258", "PC00151",
                                   "PC00208", "PC00151", "PC00151",
                                   "PC00125", "PC00233", "PC00151",
                                   "PC00151", "PC00151", "PC00227",
                                   "PC00095", "PC00095",
                                   "PC00151", "PC00151", "PC00151",
                                   "PC00151", "PC00151",
                                   "PC00151", "PC00151", "PC00227",
                                   "PC00197", "PC00197"))

  for (i in seq_len(nrow(all_df))) {
    if (!is.na(all_df$CLASS_ID[i])) {
      next
    }
    for (j in seq_len(nrow(patch))) {
      if (grepl(patch$so_SHORTNAME[j], all_df$so_SHORTNAME[i])) {
        all_df$CLASS_ID[i] <- patch$CLASS_ID[j]
        break
      }
    }
  }

  # we require the child and parent pair (we do not want grandchildren)
  parent_df <- build_classifications(all_df)
  # now merge back
  final_df <- merge(en, parent_df, by.x = "so_PID", by.y = "UNIPROT",
                    all = TRUE, sort = FALSE)

  # list out the unclassified proteins
  # unclassified <- subset(final_df, is.na(UCLASS_TERM))

  # any proteins with UCLASS_ID of NA are unclassified
  final_df$UCLASS_TERM <- ifelse(is.na(final_df$UCLASS_ID),
                                 "Unclassified", final_df$UCLASS_TERM)
  final_df$CLASS_TERM <- ifelse(is.na(final_df$UCLASS_ID),
                                "Unclassified", final_df$CLASS_TERM)

  # save this output
  write.csv(final_df, paste0("Output/Data/", label, "panther.csv"))

  sum_df <- final_df %>%
    dplyr::group_by(UCLASS_TERM, CLASS_TERM) %>%
    dplyr::summarise(n = dplyr::n())

  png(paste0("Output/Plots/", label, "Treemap.png"),
      width = 1200, height = 800)
  treemap(sum_df,
          index = c("UCLASS_TERM", "CLASS_TERM"),
          vSize = "n",
          type = "index",
          overlap.labels = 0.5,
          title = "")
  dev.off()

  return(sum_df)
}

pthOrganisms(PANTHER.db) <- "HUMAN"

# harmonise the treemaps

pe <- process_and_make_treemap(path = "Data/rankTable_INVvsControlPE.txt",
                               label = "PE")
dm <- process_and_make_treemap(path = "Data/rankTable_WTvsControlDM.txt",
                               label = "DM")
both <- process_and_make_treemap(path = "Data/rankTable_INVvsControlBoth.txt",
                                 label = "both")

# process combined data
mydata <- read.delim("Data/combined_proteome.txt")
combined <- process_and_make_treemap(path = NULL,
                                     label = "combined", en = mydata)

# merge. First merge give n.x and n.y
all <- merge(pe, dm, by = c("UCLASS_TERM", "CLASS_TERM"), all = TRUE)
names(all) <- c("UCLASS_TERM", "CLASS_TERM", "n.pe", "n.dm")
# next merge adds n
all <- merge(all, both, by = c("UCLASS_TERM", "CLASS_TERM"), all = TRUE)
names(all) <- c("UCLASS_TERM", "CLASS_TERM", "n.pe ", "n.dm", "n.both")
all <- merge(all, combined, by = c("UCLASS_TERM", "CLASS_TERM"), all = TRUE)
# replace NA with 0
all[is.na(all)] <- 0
names(all) <- c("UCLASS_TERM", "CLASS_TERM",
                "n.pe", "n.dm", "n.both", "n.combined")

make_a_treemap <- function(df, label = "", N = NULL) {
  df$m <- rowSums(df[, 3:5])
  pdf(paste0("Output/Plots/", label, "Treemap.pdf"), width = 13.5, height = 7)
  treemap(df,
          index = c("UCLASS_TERM", "CLASS_TERM"),
          vSize = N,
          type = "index",
          fontcolor.labels = c("#000000", "#222222"),
          fontsize.labels = c(16, 12),
          fontface.labels = c(2, 1),
          bg.labels = c("transparent"),
          align.labels = list(
            c("left", "top"),
            c("right", "bottom")
          ),
          overlap.labels = 0.5,
          title = "",
          palette = "Set3",
          drop.unused.levels = FALSE)
  dev.off()
}

make_a_treemap(all, "pe2", N = "n.pe")
make_a_treemap(all, "dm2", N = "n.dm")
make_a_treemap(all, "both2", N = "n.both")
make_a_treemap(all, "combined2", N = "n.combined")

# load the combined panther data frame back in
combined_df <- read.csv("Output/Data/combinedpanther.csv")
combined_df <- combined_df[order(
  combined_df$UCLASS_TERM, combined_df$CLASS_TERM), ]
# keep so_SHORTNAME, so_NAME, so_PID, UCLASS_TERM and CLASS_TERM
combined_df <- combined_df[, c("so_SHORTNAME",
                               "so_NAME", "so_PID",
                               "UCLASS_TERM", "CLASS_TERM")]
# change names to "Gene name", "Protein Name", "UniProt ID", "Class", "Subclass"
names(combined_df) <- c("Protein",
                        "Protein name", "Protein ID", "Class", "Subclass")
# save this output
write.csv(combined_df, "Output/Data/combinedpanther_sort.csv")
nrow(combined_df[combined_df$UCLASS_TERM == "Unclassified", ])
# # output PDF after downsizing by 50%:
# # Uncategorized category is 48.519 x 58.907 mm
# # so unit are per protein is in mm2
# (48.519 * 58.907) / 92
# # draw a box to represent 1 protein should be one side of
# sqrt((48.519 * 58.907) / 92 * 1)
