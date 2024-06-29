# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

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

# use the faux-Igor output from R - combined data
data_all <- read.delim("Data/universe_proteome.txt", header = TRUE)
hits <- read.delim("Data/combined_proteome.txt", header = TRUE)

# # we have entries like "ANXA2;ANXA2P2" so...
data_all$UNIPROT <- make_shortnames_unique(data_all$so_PID)
hits$UNIPROT <- make_shortnames_unique(hits$so_PID)
entrez_all <- bitr(data_all$UNIPROT,
                   fromType = "UNIPROT", toType = "ENTREZID",
                   OrgDb = "org.Hs.eg.db")
df <- merge(data_all, entrez_all, by = "UNIPROT", all = TRUE)
# remove na
df <- df[!is.na(df$ENTREZID), ]
# separate out the INV proteins
df_en <- df[df$UNIPROT %in% hits$UNIPROT, ]

make_ego_plots <- function(obj, s) {
  l <- obj@result[["Description"]][1:50]
  # exclude chromatin, chromosome terms
  l <- l[!grepl("chrom", l)]
  # exclude duplicate MF term
  l <- l[!grepl("containing anhydrides", l)]
  p1 <- goplot(obj)
  ggsave(filename = paste0("Output/Plots/", "cp_go_", s, ".pdf"), plot = p1,
         width = 20, height = 20, units = "cm", dpi = 300, bg = "white")
  p2 <- dotplot(obj, showCategory = l[1:15], font.size = 12)
  p2 <- p2 + scale_size_area(max_size = 6)
  ggsave(filename = paste0("Output/Plots/", "cp_dot_", s, ".pdf"), plot = p2,
         width = 16, height = 12, units = "cm")
  p3 <- cnetplot(obj, showCategory = l[1:10],
                 cex.params = list(category_label = 0.6, gene_label = 0.4),
                 max.overlaps = 12)
  ggsave(filename = paste0("Output/Plots/", "cp_cnet_", s, ".pdf"), plot = p3,
         width = 25, height = 19, units = "cm", dpi = 300, bg = "white")
}

ego <- enrichGO(gene = unique(df_en$ENTREZID),
                universe = unique(df$ENTREZID),
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego2 <- enrichGO(gene = unique(df_en$ENTREZID),
                 universe = unique(df$ENTREZID),
                 OrgDb = org.Hs.eg.db,
                 keyType = "ENTREZID",
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05,
                 readable = TRUE)

ego3 <- enrichGO(gene = unique(df_en$ENTREZID),
                 universe = unique(df$ENTREZID),
                 OrgDb = org.Hs.eg.db,
                 keyType = "ENTREZID",
                 ont = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05,
                 readable = TRUE)

make_ego_plots(ego, "CC_enrichedOverAll")
make_ego_plots(ego2, "BP_enrichedOverAll")
make_ego_plots(ego3, "MF_enrichedOverAll")

# write.csv(df_en$ENTREZID, file = "Output/Data/df_en.csv", row.names = FALSE)

# KEGG

kk <- enrichKEGG(gene = unique(df_en$UNIPROT),
                 universe = unique(df$UNIPROT),
                 organism  = "hsa",
                 keyType = "uniprot",
                 pvalueCutoff = 0.05)
dotplot(kk, showCategory = 20, font.size = 8)
cnetplot(kk, showCategory = 10,
         cex.params = list(category_label = 0.6,
                           gene_label = 0.4),
         max.overlaps = 12)
