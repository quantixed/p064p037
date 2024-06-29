if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager") 
}
# BiocManager::install("biomaRt")
library(httr)
library(stringr)
library(ggplot2)
library(biomaRt)
library(dplyr)
library(tidyr)
library(cowplot)

## FUNCTIONS ----

isJobReady <- function(jobId) {
  pollingInterval <- 5
  nTries <- 20
  for (i in 1:nTries) {
    url <- paste("https://rest.uniprot.org/idmapping/status/", jobId, sep = "")
    r <- GET(url = url, accept_json())
    status <- content(r, as = "parsed")
    if (!is.null(status[["results"]]) || !is.null(status[["failedIds"]])) {
      return(TRUE)
    }
    if (!is.null(status[["messages"]])) {
      print(status[["messages"]])
      return(FALSE)
    }
    Sys.sleep(pollingInterval)
  }
  return(FALSE)
}

retrieveUniprotInfo <- function(x) {
  files <- list(
    ids = paste0(x, collapse = ", "),
    from = "Gene_Name",
    to = "UniProtKB-Swiss-Prot",
    taxId = "9606"
  )
  r <- POST(url = "https://rest.uniprot.org/idmapping/run", body = files,
            encode = "multipart", accept_json())
  submission <- content(r, as = "parsed", encoding = "UTF-8")

  if (isJobReady(submission[["jobId"]])) {
    url <- paste("https://rest.uniprot.org/idmapping/details/",
                 submission[["jobId"]], sep = "")
    r <- GET(url = url, accept_json())
    details <- content(r, as = "parsed", encoding = "UTF-8")
    url <- details[["redirectURL"]]
    url <- paste(url, "?fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Cft_transmem%2Clength%2Ccc_function%2Ccc_subcellular_location%2Cgo_p%2Cgo_c&format=tsv&size=500", sep = "")
    r <- GET(url = url, accept_json())
    resultsTable <- read.table(text = content(r),
                               sep = "\t", header = TRUE, fill = TRUE,
                               quote = "")

    while (exists("link", where = r$headers)) {
      url <- str_extract(r$headers$link[1], "(?<=<).*(?=>)")
      r <- GET(url = url, accept_json())
      resultsTable <- rbind(resultsTable,
                            read.table(text = content(r),
                                       sep = "\t", header = TRUE, fill = TRUE,
                                       quote = ""))
    }
  }
  return(resultsTable)
}

incidence_of_tm <- function(df1, df2,
                            df1label = "df1", df2label = "df2",
                            tmopt = FALSE, lab = "") {
  # incidence of transmembrane domain-containing proteins
  cat("\nIncidence of TM domain containing proteins\n\n")
  temp1 <- df1 %>%
    summarise(count = n(), count1 = sum(tms > 0))
  temp2 <- df2 %>%
    summarise(count = n(), count1 = sum(tms > 0))
  tmdata <- rbind(
    c(temp1[1, 2], temp1[1, 1] - temp1[1, 2]),
    c(temp2[1, 2], temp2[1, 1] - temp2[1, 2])
  )
  rownames(tmdata) <- c(df1label, df2label)
  colnames(tmdata) <- c("TM", "noTM")
  print(tmdata)
  print(chisq.test(tmdata), similate.p.value = TRUE)

  # incidence of number of TMDs in TMD-containing population
  cat("Frequency of TM domain incidences\n\n")
  temp3 <- df1 %>%
    group_by(tms) %>%
    summarise(count = n())
  temp4 <- df2 %>%
    group_by(tms) %>%
    summarise(count = n())
  temp5 <- merge(temp3, temp4, by = "tms", all.x = TRUE)
  temp5$count.y <- ifelse(is.na(temp5$count.y) == TRUE, 0, temp5$count.y)
  # get the "names"
  n <- temp5$tms
  # transpose the data columns (excluding "name" column)
  temp6 <- as.data.frame(t(temp5[, -1]))
  colnames(temp6) <- n
  # remove the 0 tm column?
  if (tmopt == TRUE) {
    tmincidence <- temp6[, -1]
  } else {
    tmincidence <- temp6
  }
  rownames(tmincidence) <- c(df1label, df2label)
  # chi squared test
  print(tmincidence)
  print(chisq.test(tmincidence, simulate.p.value = TRUE))

  # now do posthoc test and return p-value results
  pval <- posthoc(tmincidence)
  print(pval)

  # plotting
  propsdf <- temp5
  names(propsdf) <- c("tms", df1label, df2label)
  propsdf[[df1label]] <- propsdf[[df1label]] / sum(propsdf[[df1label]])
  propsdf[[df2label]] <- propsdf[[df2label]] / sum(propsdf[[df2label]])

  dflabels <- c(df1label, df2label)
  plotdf <- propsdf %>% pivot_longer(
    cols = all_of(dflabels),
    names_to = "group",
    values_to = "freq"
  )
  ggplot(plotdf, aes(x = tms, y = freq, fill = group)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("#009988", "#bbbbbb")) +
    labs(x = "Transmembrane domains", y = "Fraction") +
    scale_x_continuous(limits = c(0.5, max(plotdf$tms) + 0.5),
                       breaks = seq(from = 1, to = max(plotdf$tms), by = 2)) +
    theme_cowplot(9) +
    theme(legend.position = "none")
  plotpath <- paste0("Output/Plots/", lab, "tm_all", df1label, df2label, ".pdf")
  ggsave(plotpath, width = 9, height = 4, units = "cm", dpi = 300)

  plotdf %>%
    filter(tms > 0) %>%
    ggplot(aes(x = tms, y = freq, fill = group)) +
    scale_fill_manual(values = c("#009988", "#bbbbbb")) +
    geom_col(position = "dodge") +
    labs(x = "Transmembrane domains", y = "Fraction") +
    scale_x_continuous(limits = c(0.5, max(plotdf$tms) + 0.5),
                       breaks = seq(from = 1, to = max(plotdf$tms), by = 2)) +
    theme_cowplot(9) +
    theme(legend.position = "none")
  plotpath <- paste0("Output/Plots/", lab, "tm_gt1", df1label, df2label, ".pdf")
  ggsave(plotpath, width = 9, height = 4, units = "cm", dpi = 300)
}

posthoc <- function(df) {
  mnames <- colnames(df)
  nc <- ncol(df)
  m <- as.matrix(df)
  m_sum <- rowSums(m)
  v <- NULL
  for (i in 1:nc) {
    n <- rbind(
      c(m[1, i], m_sum[1] - m[1, i]),
      c(m[2, i], m_sum[2] - m[2, i])
    )
    v <- rbind(v, chisq.test(n, simulate.p.value = TRUE)$p.value)
  }
  result <- data.frame(
    tms = mnames,
    p.value = v
  )
  alpha <- 0.05 / (nc - 1)
  result$bonferroni <- ifelse(result$p.value < alpha, "reject", "accept")
  # Bonferroni-Holm - order by p-value
  result <- result[order(result$p.value, decreasing = FALSE), ]
  result$bh <- ifelse(result$p.value < 1 &
                        result$p.value < (0.05 / (nc - (row(result) - 1))),
                      "reject", "accept")[, 1]
  result <- result[order(as.numeric(result$tms), decreasing = FALSE), ]
  return(result)
}

incidence_of_secreted <- function(df1, df2,
                                  df1label = "df1", df2label = "df2") {
  # incidence of transmembrane domain-containing proteins
  cat("\nIncidence of secreted proteins\n\n")

  v1 <- str_count(df1$Subcellular.location..CC., "Secrete")
  v2 <- str_count(df2$Subcellular.location..CC., "Secrete")

  sdata <- rbind(
    c(sum(v1), nrow(df1) - sum(v1)),
    c(sum(v2), nrow(df2) - sum(v2))
  )
  rownames(sdata) <- c(df1label, df2label)
  colnames(sdata) <- c("Secreted", "NotSecreted")
  print(sdata)
  print(chisq.test(sdata), simulate.p.value = TRUE)
}

process_output_dataset <- function(path1, path2 = NULL, label) {
  if (is.null(path2)) {
    # means that path1 is universe with hits marked with colorWave = 3
    # import our data
    universe <- read.delim(path1)
  } else {
    # import our data
    hits <- read.delim(path1)
    universe <- read.delim(path2)
    # check that all unique hits are in universe
    uhits <- unique(hits$so_SHORTNAME)
    uuniverse <- unique(universe$so_SHORTNAME)
    if (length(uhits) != length(uhits %in% uuniverse)) {
      cat("Some hits are not found in the full dataset (universe)\n")
    }
  }
  # we need to use the first instance if there are semicolon separated values
  universe$so_SHORTNAME <- sapply(universe$so_SHORTNAME,
                                  function(x) strsplit(x, ";")[[1]][1])
  universe$so_PID <- sapply(universe$so_PID,
                            function(x) strsplit(x, ";")[[1]][1])
  hits$so_SHORTNAME <- sapply(hits$so_SHORTNAME,
                              function(x) strsplit(x, ";")[[1]][1])
  hits$so_PID <- sapply(hits$so_PID,
                        function(x) strsplit(x, ";")[[1]][1])
  # retrieve the uniprot data
  uni_results <- retrieveUniprotInfo(universe$so_SHORTNAME)
  # merge our data with the results from uniprot
  df <- merge(universe, uni_results, by.x = "so_PID",
              by.y = "Entry", all.x = TRUE, all.y = FALSE)
  # df <- merge(universe, uni_results, by.x = c("so_SHORTNAME", "so_PID"),
  #              by.y = c("From", "Entry"), all.x = TRUE, all.y = FALSE)
  # get rid of ones where we found nothing
  df <- df[!is.na(df$Entry), ]
  # count the number of transmembrane domains in each entry
  df$tms <- str_count(df$Transmembrane, "TRANSMEM")
  # at this point we can have duplicated entries where one has
  # 0 TM domains and the other has more we will reject the duplicate(s)
  # with the lowest value, because duplicated takes the first value
  # we can sort first before pruning duplicates
  df <- df[order(df$so_SHORTNAME, -df$tms), ]
  df <- df[!duplicated(df$so_SHORTNAME), ]
  # isolate our "hits"
  if (is.null(path2)) {
    hits <- df[df$so_colorWave == 3, ]
  } else {
    hits <- df[df$so_SHORTNAME %in% hits$so_SHORTNAME, ]
  }

  # generate results ~ INV vs whole IP
  cat(paste0("*** ", label, " ***\n"))
  incidence_of_tm(df1 = df, df2 = hits,
                  df1label = "IP", df2label = "INV", lab = label)
  incidence_of_secreted(df1 = df, df2 = hits,
                        df1label = "IP", df2label = "INV")

  # save data note that the fold change and p-value relates to universe
  # so only use output for Uniprot data
  write.csv(df, "Output/Data/universe_proteome_tm.csv", row.names = FALSE)
  write.csv(hits, "Output/Data/combined_proteome_tm.csv", row.names = FALSE)
}

## SCRIPT ----

con <- file("Output/Data/tmCalc.txt")
sink(con, append = TRUE)
process_output_dataset(path1 = "Data/combined_proteome.txt",
                       path2 = "Data/universe_proteome.txt",
                       label = "combined")
sink()
