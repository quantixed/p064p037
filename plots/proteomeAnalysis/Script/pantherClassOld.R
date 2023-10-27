# this script works well and does what we want
# uses v18, and we are missing protein classifications

# BiocManager::install("PANTHER.db")
library(PANTHER.db)
library(treemap)

# functions ----

pull_ancestors <- function(term) {
  ancestors <- traverseClassTree(PANTHER.db, term, scope="ANCESTOR")
  if(length(ancestors) > 2) {
    df <- data.frame(CLASS_ID = term,
                     UCLASS_TERM = NA,
                     UCLASS_ID = NA)
  }
  if(length(ancestors) == 1) {
    df <- data.frame(CLASS_ID = term,
                     UCLASS_TERM = "Protein Class",
                     UCLASS_ID = "P00000")
  }
  if (length(ancestors) == 2) {
    parents <- traverseClassTree(PANTHER.db, term, scope="PARENT")
    parents <- select(PANTHER.db, ancestors, "CLASS_TERM", "CLASS_ID")
    df <- data.frame(CLASS_ID = term,
                     UCLASS_TERM = parents$CLASS_TERM[1],
                     UCLASS_ID = parents$CLASS_ID[1])
  }
  
  return(df)
}

make_shortnames_unique <- function(x) {
  y <- x
  
  for(i in 1:length(y)) {
    z <- unlist(strsplit(y[i],";"))
    if(length(z) == 1) {
      next
    }
    for(j in 1:length(z)) {
      if(length(unique(y)) < length(unique(append(y,z[j])))) {
        y[i] <- z[j]
        break
      }
    }
  }
  
  return(y)
}

process_and_make_treemap <- function(path,label = "") {
  mydata <- read.delim(path)
  en <- mydata[mydata$so_colorWave == 3, ]
  # en_up <- sapply(strsplit(en$so_PID,";"), `[`, 1)
  # we need a unique list of protein IDs
  en_up <- make_shortnames_unique(en$so_PID)
  
  # Families is too low level so we will omit
  mylist <- select(PANTHER.db, en_up, columns = c("CLASS_ID","CLASS_TERM"), keytype = "UNIPROT")
  # we can add jointype = "left" to the previous line to get the unmatched ones too
  
  ## Protein class *the top-most level* is the grandparent.
  # we require the parent and child pair (we do not want grandchildren) for the treemap
  
  parent_df <- data.frame()
  for(i in 1:nrow(mylist)) {
    new_df <- pull_ancestors(mylist$CLASS_ID[i])
    parent_df <- rbind(parent_df, new_df)
  }
  # there are duplicated rows so let's eliminate those
  parent_df <- parent_df[!duplicated(parent_df),]
  # now merge so we get matches per original class ID
  all_df <- merge(mylist, parent_df, by = "CLASS_ID", all.x = TRUE, sort = FALSE)
  # NAs should be a grandchild
  all_df <- all_df[!is.na(all_df$UCLASS_ID),]
  # parent-child combos should be everything where UCLASS_ID does not equal P00000
  upper_df <- all_df[all_df$UCLASS_ID == "P00000",]
  lower_df <- all_df[all_df$UCLASS_ID != "P00000",]
  final_df <- merge(upper_df, lower_df, by = "UNIPROT", all = TRUE, sort = FALSE)
  # these are all the mapped terms
  final_df <- 
    # save this output
    write.csv(final_df, paste0("Output/Data/",label,"panther.csv"))
  
  sum_df <- final_df %>%
    dplyr::group_by(CLASS_TERM.x,CLASS_TERM.y) %>% 
    dplyr::summarise(n = dplyr::n())
  
  png(paste0("Output/Plots/",label,"Treemap.png"), width = 1200, height = 800)
  treemap(sum_df,
          index = c("CLASS_TERM.x","CLASS_TERM.y"),
          vSize = "n",
          type = "index",
          overlap.labels = 0.5,
          title="")
  dev.off()
  
  return(sum_df)
}

# pthOrganisms(PANTHER.db) <- "HUMAN"

# harmonise the treemaps

pe <- process_and_make_treemap("Data/rankTable_INVvsControlPE.txt","PE")
dm <- process_and_make_treemap("Data/rankTable_WTvsControlDM.txt","DM")
both <- process_and_make_treemap("Data/rankTable_INVvsControlBoth.txt","both")

# make childless terms have the child name of x
pe$CLASS_TERM.y <- ifelse(is.na(pe$CLASS_TERM.y), pe$CLASS_TERM.x, pe$CLASS_TERM.y)
dm$CLASS_TERM.y <- ifelse(is.na(dm$CLASS_TERM.y), dm$CLASS_TERM.x, dm$CLASS_TERM.y)
both$CLASS_TERM.y <- ifelse(is.na(both$CLASS_TERM.y), both$CLASS_TERM.x, both$CLASS_TERM.y)

# merge. First merge give n.x and n.y
all <- merge(pe, dm, by = c("CLASS_TERM.x","CLASS_TERM.y"), all = TRUE)
# next merge adds n
all <- merge(all, both, by = c("CLASS_TERM.x","CLASS_TERM.y"), all = TRUE)
# replace NA with 0
all[is.na(all)] <- 0
names(all) <- c("CLASS_TERM.x","CLASS_TERM.y","n.pe","n.dm","n.both")

make_a_treemap <- function(df,label = "", N = NULL) {
  df$m <- rowSums(df[,3:5])
  png(paste0("Output/Plots/",label,"Treemap.png"), width = 1200, height = 800)
  treemap(df,
          index = c("CLASS_TERM.x","CLASS_TERM.y"),
          vSize = N,
          type = "index",
          align.labels=list(
            c("left", "top"), 
            c("center", "center")
          ),
          overlap.labels = 0.5,
          title="",
          palette = "Set3",
          drop.unused.levels = FALSE)
  dev.off()
}

make_a_treemap(all,"pe2", N = "n.pe")
make_a_treemap(all, "dm2", N = "n.dm")
make_a_treemap(all, "both2", N = "n.both")
