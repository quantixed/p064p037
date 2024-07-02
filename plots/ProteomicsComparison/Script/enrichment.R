# enrichment analysis

## Functions ----

enrichment_test <- function(en1, total1,
                            en2, total2,
                            label1, label2, title = "Comparison") {
  df <- data.frame(
    "Positive" = c(en1, en2),
    "Negative" = c(total1 - en1, total2 - en2)
  )

  rownames(df) <- c(label1, label2)
  cat(title, "\n")
  print(df)
  cat(paste(label2, "positive (%):", en2 / total2 * 100, "\n",
            label1, "positive (%):", en1 / total1 * 100, "\n"))
  print(chisq.test(df), similate.p.value = TRUE)
}

## Script ----

hits <- read.delim("Data/combined_proteome.txt")
universe <- read.delim("Data/universe_proteome.txt")
# there are some entries in here that contain a semicolon i.e. more than one
# gene name, but this doesn't affect the enrichment test
SLMV <- read.csv("Output/Data/SLMV.txt")
ATG <- read.csv("Output/Data/ATG.txt")
CCV <- read.delim("Data/CCV.txt")

enrichment_test(
  en1 = sum(universe$so_SHORTNAME %in% ATG$Gene.names, na.rm = TRUE),
  total1 = nrow(universe),
  en2 = sum(hits$so_SHORTNAME %in% ATG$Gene.names, na.rm = TRUE),
  total2 = nrow(hits),
  label1 = "IP",
  label2 = "INV",
  title = "INVs and ATG dataset"
)

enrichment_test(
  en1 = sum(universe$so_SHORTNAME %in% SLMV$Gene.names, na.rm = TRUE),
  total1 = nrow(universe),
  en2 = sum(hits$so_SHORTNAME %in% SLMV$Gene.names, na.rm = TRUE),
  total2 = nrow(hits),
  label1 = "IP",
  label2 = "INV",
  title = "INVs and SLMV dataset"
)

enrichment_test(
  en1 = sum(universe$so_SHORTNAME %in% CCV$Genes, na.rm = TRUE),
  total1 = nrow(universe),
  en2 = sum(hits$so_SHORTNAME %in% CCV$Genes, na.rm = TRUE),
  total2 = nrow(hits),
  label1 = "IP",
  label2 = "INV",
  title = "INVs and CCV dataset"
)
