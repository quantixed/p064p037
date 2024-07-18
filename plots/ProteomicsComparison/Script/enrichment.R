# enrichment analysis
library(ggplot2)

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

# ATG9 incidence in INVs
pie_maker(x1 = c(143, 459),
          cols = c("#44aa99", "#117733"), label = "pie_atg9_incidence_inv")
# ATG9 incidence in all and in INVs
pie_maker(x1 = c(392 - 143, 143, 459, 3264 - 459),
          cols = c("#bbbbbb", "#44aa99", "#117733", "#dddddd"),
          x2 = c(392, 3264), label = "pie_atg9_incidence_ip")
# SLMV incidence in INVs
pie_maker(x1 = c(33, 569),
          cols = c("#44aa99", "#117733"), label = "pie_SLMV_incidence_inv")
# SLMV incidence in all and in INVs
pie_maker(x1 = c(67 - 33, 33, 569, 3589 - 569),
          cols = c("#bbbbbb", "#44aa99", "#117733", "#dddddd"),
          x2 = c(67, 3589), label = "pie_SLMV_incidence_ip")
# CCV incidence in INVs
pie_maker(x1 = c(9, 593),
          cols = c("#44aa99", "#117733"), label = "pie_CCV_incidence_inv")
# CCV incidence in all and in INVs
pie_maker(x1 = c(61 - 9, 9, 593, 3595 - 593),
          cols = c("#bbbbbb", "#44aa99", "#117733", "#dddddd"),
          x2 = c(61, 3595), label = "pie_CCV_incidence_ip")