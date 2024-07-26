library(ggplot2)
library(ggforce)
library(dplyr)
library(tidyr)
library(reshape2)

## Functions ----

process_data <- function() {
  all_files <- list.files("Data", pattern = "*.csv", full.names = TRUE)
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

process_time <- function() {
  all_dirs <- list.dirs("Data/timestamps")
  switch <- 1
  for (dirname in all_dirs) {
    all_files <- list.files(dirname, pattern = "*.txt", full.names = TRUE)
    # text files are a single column of data with no header
    for (filename in all_files) {
      times <- scan(file = filename, quiet = TRUE)
      # we want to keep every 4th row and put into data frame
      df_temp <- data.frame(tt = times[seq(1, length(times), 4)])
      # we need to create a column that shows the Frame numbers
      df_temp$Frame <- seq_len(nrow(df_temp))
      fname <- tools::file_path_sans_ext(basename(filename))
      # edit filename to match the image file that was analysed
      fname <- gsub("_deltaT", "-Registered", fname)
      dname <- basename(dirname(filename))
      df_temp$id <- paste(dname, fname, sep = "_")
      if (switch == 1) {
        df_all <- df_temp
        switch <- 0
      } else {
        df_all <- rbind(df_all, df_temp)
      }
    }
  }

  return(df_all)
}

## Script ----

data <- process_data()
# condition extracted
data$Condition <- sapply(strsplit(data$Filename, "_"), "[[", 2)
data$Condition <- tolower(data$Condition)
# generate unique id for each cell
data$id <- paste(data$expt, data$Filename, sep = "_")
# to do background subtraction
# we need to match up frames for ROI 1,2,3 with ROI 0
data_wide <- data %>%
  select(!Area) %>%
  pivot_wider(names_from = ROI, values_from = Mean)
data_wide$cell <- data_wide$"1" - data_wide$"0"
data_wide$mito <- data_wide$"2" - data_wide$"0"
data_wide$cyto <- data_wide$"3" - data_wide$"0"
data_wide$mito_cell <- data_wide$mito / data_wide$cell
data_wide$cyto_cell <- data_wide$cyto / data_wide$cell
data_wide$mito_cyto <- data_wide$mito / data_wide$cyto

# timestamp data
tstamps <- process_time()
# add time to data. To do this we need an id/frame key
data_wide$key <- paste(data_wide$id, data_wide$Frame, sep = "_")
tstamps$key <- paste(tstamps$id, tstamps$Frame, sep = "_")
data_wide <- plyr::join(data_wide, tstamps)

# we need Area for ROI = 2 or 3
temp_df <- data %>%
  filter(ROI > 1) %>%
  select(!Mean) %>%
  pivot_wider(names_from = ROI, values_from = Area)
temp_df$key <- paste(temp_df$id, temp_df$Frame, temp_df$Channel, sep = "_")
temp_df <- data.frame(
  key = temp_df$key,
  area_mito = temp_df$"2",
  area_cyto = temp_df$"3"
)
data_wide$key <- paste(data_wide$id,
  data_wide$Frame, data_wide$Channel,
  sep = "_"
)
data_wide <- plyr::join(data_wide, temp_df)
# clean up a bit
data_wide[, c("key", "0", "1", "2", "3")] <- list(NULL)

# calculate intden
data_wide$intden_mito <- data_wide$mito * data_wide$area_mito
data_wide$intden_cyto <- data_wide$cyto * data_wide$area_cyto
data_wide$frac_mito <- data_wide$intden_mito /
  (data_wide$intden_mito + data_wide$intden_cyto)
data_wide$frac_mean_mito <- data_wide$mito /
  (data_wide$mito + data_wide$cyto)

write.csv(data_wide, "Output/Data/data_wide.csv", row.names = FALSE)
