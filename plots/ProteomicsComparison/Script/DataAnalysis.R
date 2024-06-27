library(eulerr)
library(ggplot2)
library(ggrepel)
library(UpSetR)

# need to read in the fold-change data for INVs
INV <- read.delim("Output/Data/combined_inv_proteome_with_intensities.txt")
INV <- INV[, c("so_SHORTNAME", "so_ratioWave")]
names(INV) <- c("Gene.names", "INV")

# Read in cleaned data files.
# Should have a col of gene names and LFQ values (or other quantitative measure)
SLMV <- read.csv("Output/Data/SLMV.txt")
ATG <- read.csv("Output/Data/ATG.txt")

names(SLMV)[names(SLMV) == "uid"] <- "SLMV_uid"
names(ATG)[names(ATG) == "uid"] <- "ATG_uid"

# Merge on gene names
IS <- merge(INV, SLMV)
IA <- merge(INV, ATG)
AS <- merge(ATG, SLMV)
ISA <- merge(IS, ATG)
# if NA in any, we cannot compare
IS <- IS[complete.cases(IS), ]
IA <- IA[complete.cases(IA), ]
AS <- AS[complete.cases(AS), ]
ISA <- ISA[complete.cases(ISA), ]

# change 0 to random value with mean of 100 sd of
ISA[ISA == 0] <- rnorm(n = length(ISA[ISA == 0]), mean = 100, sd = 10)

# correlation calculation
cor.test(x = log10(IA$INV), y = IA$ATG, method = "pearson")
# Pearson's product-moment correlation
#
# data:  log10(IA$INV) and IA$ATG
# t = 2.5867, df = 141, p-value = 0.0107
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.05046165 0.36426871
# sample estimates:
#       cor
# 0.2128472

INVvsSLMV_Plot <- ggplot(IS, aes(x = INV, y = SLMV, label = Gene.names)) +
  geom_point(size = 0.5) +
  geom_text_repel(size = 1.5, max.overlaps = 10, segment.color = "#7f7f7f7f") +
  scale_x_continuous(name = "INV", trans = "log2") +
  scale_y_continuous(name = "SLMV", trans = "log10") +
  theme_classic(9)

ggsave("Output/Plots/INVvsSLMV_Plot.pdf", INVvsSLMV_Plot,
  width = 8, height = 8,
  bg = "white", units = "cm"
)

INVvsATG_Plot <- ggplot(IA, aes(x = INV, y = ATG, label = Gene.names)) +
  geom_point(size = 0.5) +
  geom_text_repel(size = 1.5, max.overlaps = 10, segment.color = "#7f7f7f7f") +
  scale_x_continuous(name = "INV", trans = "log2") +
  scale_y_continuous(name = "ATG9A", trans = "log10", limits = c(1e+6, 1e+10)) +
  theme_classic(9)

ggsave("Output/Plots/INVvsATG_Plot.pdf", INVvsATG_Plot,
  width = 8, height = 8,
  bg = "white", units = "cm"
)

ATGvsSLMV_Plot <- ggplot(AS, aes(x = ATG, y = SLMV, label = Gene.names)) +
  geom_point(size = 0.5) +
  geom_text_repel(size = 1.5, max.overlaps = 10, segment.color = "#7f7f7f7f") +
  scale_x_continuous(name = "ATG9A", trans = "log10", limits = c(1e+6, 1e+10)) +
  scale_y_continuous(name = "SLMV") +
  theme_classic(9)

ggsave("Output/Plots/ATGvsSLMV_Plot.pdf", ATGvsSLMV_Plot,
  width = 8, height = 8,
  bg = "white", units = "cm"
)

## Euler plot

make_euler_2 <- function(A, B, AB, names = c("A", "B")) {
  set1 <- c(
    A = A - AB,
    B = B - AB,
    "A&B" = AB
  )
  # change names
  name <- paste(names, collapse = "&")
  names(set1) <- c(names, name)
  fit1 <- euler(set1, shape = "ellipse")
  p <- plot(fit1, quantities = TRUE)
  name <- paste(names, collapse = "_vs_")
  ggsave(paste0("Output/Plots/", name, ".pdf"), p,
    width = 3.5, height = 3.5, units = "cm"
  )
}

make_euler_2(
  A = nrow(INV),
  B = length(unique(SLMV$SLMV_uid)),
  AB = length(unique(IS$SLMV_uid)),
  names = c("INV", "SLMV")
)
make_euler_2(
  A = nrow(INV),
  B = length(unique(ATG$ATG_uid)),
  AB = length(unique(IA$ATG_uid)),
  names = c("INV", "ATG9")
)

make_euler_3 <- function(A, B, C, AB, BC, AC, ABC, names = c("A", "B", "C")) {
  set1 <- c(
    A = A - AB - AC + ABC,
    B = B - AB - BC + ABC,
    C = C - BC - AC + ABC,
    "A&B" = AB - ABC,
    "A&C" = AC - ABC,
    "B&C" = BC - ABC,
    "A&B&C" = ABC
  )
  # change names
  name <- paste(names, collapse = "&")
  names(set1[1:3]) <- names
  fit1 <- euler(set1, shape = "ellipse")
  p <- plot(fit1, quantities = TRUE)
  name <- paste(names, collapse = "_vs_")
  ggsave(paste0("Output/Plots/", name, ".pdf"), p,
    width = 4, height = 4, units = "cm"
  )
}

make_euler_3(
  A = nrow(INV),
  B = length(unique(SLMV$SLMV_uid)),
  C = length(unique(ATG$ATG_uid)),
  AB = length(unique(IS$SLMV_uid)),
  BC = length(unique(AS$SLMV_uid)),
  AC = length(unique(IA$ATG_uid)),
  ABC = min(
    length(unique(ISA$ATG_uid)),
    length(unique(ISA$SLMV_uid))
  ),
  names = c("INV", "SLMV", "ATG")
)
