# enrichment analysis
library(ggplot2)
library(cowplot)

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

pie_maker <- function(x1, cols, x2 = NULL, label = "pie") {
  df <- data.frame(
    pvalues = x1,
    pcolours = as.factor(seq_along(cols))
  )
  cs <- cumsum(x1)
  df$ymin <- c(0, cs[seq_along(cs) - 1])
  df$ymax <- cs
  if (is.null(x2)) {
    df2 <- df
  } else {
    cs <- cumsum(x2)
    df2 <- data.frame(
      pvalues = x2,
      ymin = c(0, cs[seq_along(cs) - 1]),
      ymax = cs
    )
  }
  
  ggplot() +
    geom_rect(data = df, aes(fill = pcolours,
                             ymin = ymin, ymax = ymax, xmin = 0, xmax = 3)) +
    scale_fill_manual(values = cols) +
    geom_rect(data = df2, aes(ymin = ymin, ymax = ymax,
                              xmin = 0, xmax = 3),
              alpha = 0, colour = "black") +
    xlim(c(0, 3)) +
    theme_minimal() +
    theme(aspect.ratio = 1) +
    coord_polar(theta = "y")
  ggsave(paste0("Output/Plots/", label, ".pdf"),
         width = 5, height = 3, units = "cm", dpi = 300)
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


slmv <- data.frame(
  n = c(33 / 602, 67 / 3656),
  label = c("INV", "Total")
)

atg9 <- data.frame(
  n = c(143 / 602, 392 / 3656),
  label = c("INV", "Total")
)

ccv <- data.frame(
  n = c(9 / 602, 61 / 3656),
  label = c("INV", "Total")
)

# rbind these three into one data frame with a column to indicate the origin
prots <- rbind(
  cbind(slmv, df = "SLMV"),
  cbind(atg9, df = "ATG9"),
  cbind(ccv, df = "CCV")
)
prots$df <- factor(prots$df, levels = c("CCV", "ATG9", "SLMV"))
prots$label <- factor(prots$label, levels = c("Total", "INV"))

p1 <- ggplot(prots) +
  geom_col(aes(y = df, x = n, fill = label), position = "dodge") +
  scale_fill_manual(values = c("#bbbbbb","#009988")) +
  labs(x = "Fraction", y = "") +
  lims(x = c(0, 0.3)) +
  theme_cowplot(9) +
  theme(legend.position = "top")
p1
ggsave("Output/Plots/bars.pdf", p1, width = 79, height = 80, units = "mm")
