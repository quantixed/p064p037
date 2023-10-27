# Import of hits and disambiguating is fine.
# I can either map to v16 with PANTHER.db which leaves 150 proteins unclassified in a 513 dataset
# or I can import from the web interface using v17. This matches all but 100 from the same dataset.
# or I can do the same thing with rbioapi which also queries v17 and gives the same result
# probelm is that the result is difficult to parse from rbioapi and in both cases requires mapping child/parent terms
# this can PROBABLY be done with PANTHER.db but if there any changes between v16 and v17 then I will make a mistake
# best thing to do would be to use an updated form of PANTHER.db. To do that I will need to wait for Jurgen to send me the updated package
# I could check if rbioapi is easy to use for mapping


remotes::install_github("moosa-r/rbioapi",force = TRUE)
BiocManager::install("PANTHER.db")
library(rbioapi)
library(PANTHER.db)
library(treemap)

# functions ----

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


pull_ancestors <- function(id, term) {
  if(nchar(id) == 0) {
    # deal with a blank protein class id
    df <- data.frame(class_id = id,
                     CLASS_ID = id,
                     CLASS_TERM = term,
                     UCLASS_TERM = NA,
                     UCLASS_ID = NA)
    return(df)
  }
  ancestors <- traverseClassTree(PANTHER.db, id, scope="ANCESTOR")
  
  if(length(ancestors) > 2) {
    classterm <- select(PANTHER.db, ancestors[3], "CLASS_TERM", "CLASS_ID")
    uclassterm <- select(PANTHER.db, ancestors[2], "CLASS_TERM", "CLASS_ID")
    df <- data.frame(class_id = id,
                     CLASS_ID = ancestors[2],
                     CLASS_TERM = classterm$CLASS_TERM[1],
                     UCLASS_ID = ancestors[3],
                     UCLASS_TERM = uclassterm$CLASS_TERM[1])
    
  }
  # next we deal with the following cases, but term may be blank
  if(nchar(term) == 0) {
    orig <- select(PANTHER.db, id, "CLASS_TERM", "CLASS_ID")
    term <- orig$CLASS_TERM[1]
  }
  if(length(ancestors) == 1) {
    df <- data.frame(class_id = id,
                     CLASS_ID = id,
                     CLASS_TERM = term,
                     UCLASS_ID = id,
                     UCLASS_TERM = term)
  }
  if (length(ancestors) == 2) {
    details <- select(PANTHER.db, ancestors, "CLASS_TERM", "CLASS_ID")
    df <- data.frame(class_id = id,
                     CLASS_ID = id,
                     CLASS_TERM = term,
                     UCLASS_ID = details$CLASS_ID[1],
                     UCLASS_TERM = details$CLASS_TERM[1])
  }
  
  return(df)
}

## function for the workflow

make_panther_treemap <- function(path, label) {
  mydata <- read.delim(path)
  en <- mydata[mydata$so_colorWave == 3, ]
  # we need a unique list of protein IDs
  en_up <- en$up <- make_shortnames_unique(en$so_PID)
  
  result <- rba_panther_mapping(genes = en_up, organism = 9606)
  
  ## this loop will extract the PANTHER classes (v17)
  pc_df <- data.frame()
  for(i in 1:length(result[[2]]$gene)) {
    temp <- result[[2]]$gene[[i]]
    gene <- temp$accession
    pc <- pcid <- ""
    for(j in 1:length(temp$annotation_type_list$annotation_data_type)) {
      if(is.atomic(temp$annotation_type_list$annotation_data_type[[j]])) {
        break
      }
      if(!is.null(temp$annotation_type_list$annotation_data_type[[j]]) &&
         !is.null(temp$annotation_type_list$annotation_data_type[[j]]$content) &&
         temp$annotation_type_list$annotation_data_type[[j]]$content == "ANNOT_TYPE_ID_PANTHER_PC") {
        pc <- temp$annotation_type_list$annotation_data_type[[j]]$annotation_list$annotation$name
        pcid <- temp$annotation_type_list$annotation_data_type[[j]]$annotation_list$annotation$id
        break
      }
    }
    tempdf <- data.frame(accession = gene,
                         class_id = pcid,
                         class_name = pc)
    pc_df <- rbind(pc_df,tempdf)
  }
  
  pc_df$UNIPROT <- sapply(strsplit(pc_df$accession, "\\|"), function(x) x[3], simplify = FALSE)
  pc_df$UNIPROT <- sub("UniProtKB=","", pc_df$UNIPROT)
  
  all_df <- merge(en, pc_df, by.x = "up", by.y = "UNIPROT", all.x = TRUE)
  all_df[is.na(all_df)] <- ""
  # patch some missing classes
  patch <- data.frame(so_SHORTNAME = c("TLN1","VCL","^FLOT","^TPD5","VAT1","^SAR1","^RAB","^SSR","^RTN"),
                      class_id = c("PC00041","PC00041","PC00151","PC00151","PC00258","PC00151","PC00208","PC00151","PC00151"))
  
  for(i in 1:nrow(all_df)) {
    if(all_df$class_id[i] != "") {
      next
    }
    for(j in 1:nrow(patch)) {
      if(grepl(patch$so_SHORTNAME[j],all_df$so_SHORTNAME[i])) {
        all_df$class_id[i] <- patch$class_id[j]
        break
      }
    }
  }
  
  # we require the child and parent pair
  pc_ids <- unique(all_df[c("class_id", "class_name")])
  
  parent_df <- data.frame()
  for(i in 1:nrow(pc_ids)) {
    new_df <- pull_ancestors(pc_ids[i,1],pc_ids[i,2])
    parent_df <- rbind(parent_df, new_df)
  }
  
  # there are duplicated rows so let's eliminate those
  parent_df <- parent_df[!duplicated(parent_df),]
  # now merge so we get matches per original class ID
  final_df <- merge(all_df, parent_df, by.x = "class_id", by.y = "class_id", all.x = TRUE, sort = FALSE)
  # any proteins with UCLASS_ID of NA are unclassified
  final_df$UCLASS_TERM <- ifelse(is.na(final_df$UCLASS_ID),"Unclassified",final_df$UCLASS_TERM)
  final_df$CLASS_TERM <- ifelse(is.na(final_df$UCLASS_ID),"Unclassified",final_df$CLASS_TERM)
  
  # save this output
  write.csv(final_df, paste0("Output/Data/",label,"panther.csv"))

  sum_df <- final_df %>%
    dplyr::group_by(UCLASS_TERM,CLASS_TERM) %>% 
    dplyr::summarise(n = dplyr::n())
    
  png(paste0("Output/Plots/", label,"Treemap.png"), width = 1200, height = 800)
    treemap(sum_df,
            index = c("UCLASS_TERM","CLASS_TERM"),
            vSize = "n",
            type = "index",
            overlap.labels = 0.5,
            title="",
            align.labels=list(
              c("left", "top"), 
              c("center", "center")
            ),
            palette = "Set3",
            drop.unused.levels = FALSE)
  dev.off()
  
  return(sum_df)
}

make_harmony_treemap <- function(df,label = "", N = NULL) {
  df$m <- rowSums(df[,3:5])
  png(paste0("Output/Plots/",label,"Treemap.png"), width = 1200, height = 800)
  treemap(df,
          index = c("UCLASS_TERM","CLASS_TERM"),
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


pthOrganisms(PANTHER.db) <- "HUMAN"
pe <- make_panther_treemap("Data/rankTable_INVvsControlPE.txt","pe")
dm <- make_panther_treemap("Data/rankTable_WTvsControlDM.txt", "dm")
both <- make_panther_treemap("Data/rankTable_INVvsControlBoth.txt", "both")

# merge. First merge give n.x and n.y
all <- merge(pe, dm, by = "CLASS_TERM", all = TRUE)
# next merge adds n
all <- merge(all, both, by = "CLASS_TERM", all = TRUE)
# replace NA with 0
all[is.na(all)] <- 0
all <- all[, c("CLASS_TERM","UCLASS_TERM","n.x","n.y","n")]
all <- all[, names(all) %in% c("CLASS_TERM","UCLASS_TERM","n.x","n.y","n")]
names(all) <- c("CLASS_TERM","UCLASS_TERM","n.pe","n.dm","n.both")

make_harmony_treemap(all,"pe2", N = "n.pe")
make_harmony_treemap(all, "dm2", N = "n.dm")
make_harmony_treemap(all, "both2", N = "n.both")
