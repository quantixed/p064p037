BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# use violin plot output from Igor
data_all <- read.delim("Data/rankTable_INVvsControl.txt", header = TRUE)

# # we have entries like "ANXA2;ANXA2P2" so...
# data_all$SYMBOL <- gsub("\"","",data_all$V1)
# data_all$SYMBOL <- sapply(strsplit(data_all$SYMBOL,";"), `[`, 1)
# entrez_all <- bitr(data_all$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# df <- merge(data_all,entrez_all, by = "SYMBOL", all = TRUE)
# Mapping SMBOL gave a higher failure rate than using Uniprot IDs
data_all$UNIPROT <- sapply(strsplit(data_all$so_PID,";"), `[`, 1)
entrez_all <- bitr(data_all$UNIPROT, fromType="UNIPROT", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df <- merge(data_all, entrez_all, by = "UNIPROT", all = TRUE)
# remove na
df <- df[!is.na(df$ENTREZID),]
df_en <- df[df$so_colorWave == 3, ]

make_ego_plots <- function(obj,s) {
  p1 <- goplot(obj)
  ggsave(filename = paste0("Output/Plots/","go_",s,".pdf"), plot = p1, width = 20, height = 20, units = "cm", dpi = 300, bg = "white")
  p2 <- dotplot(obj, showCategory = 20, font.size = 8)
  ggsave(filename = paste0("Output/Plots/","dot_",s,".pdf"), plot = p2, width = 17, height = 17, units = "cm", dpi = 300, bg = "white")
  p3 <- cnetplot(obj, showCategory = 10, cex.params = list(category_label = 0.6, gene_label = 0.4), max.overlaps = 12)
  ggsave(filename = paste0("Output/Plots/","cnet_",s,".pdf"), plot = p3, width = 25, height = 19, units = "cm", dpi = 300, bg = "white")
}

# ggo <- groupGO(gene = unique(df_en$ENTREZID),
#                OrgDb = org.Hs.eg.db,
#                ont = "BP",
#                level = 3,
#                readable = TRUE)
# ggo@result  <- ggo@result[order(ggo@result$Count, decreasing = TRUE),]
# barplot(ggo, drop = TRUE, showCategory = 20)

ego <- enrichGO(gene = unique(df_en$ENTREZID),
                universe = unique(df$ENTREZID),
                OrgDb = org.Hs.eg.db,
                keyType = 'ENTREZID',
                ont = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego2 <- enrichGO(gene = unique(df_en$ENTREZID),
                 universe = unique(df$ENTREZID),
                 OrgDb = org.Hs.eg.db,
                 keyType = 'ENTREZID',
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05,
                 readable = TRUE)

# ego3 <- enrichGO(gene = unique(df_en$ENTREZID),
#                  OrgDb = org.Hs.eg.db,
#                  keyType = 'ENTREZID',
#                  ont = "CC",
#                  pAdjustMethod = "BH",
#                  pvalueCutoff = 0.01,
#                  qvalueCutoff = 0.05,
#                  readable = TRUE)
# 
# ego4 <- enrichGO(gene = unique(df_en$ENTREZID),
#                  OrgDb = org.Hs.eg.db,
#                  keyType = 'ENTREZID',
#                  ont = "BP",
#                  pAdjustMethod = "BH",
#                  pvalueCutoff = 0.01,
#                  qvalueCutoff = 0.05,
#                  readable = TRUE)

make_ego_plots(ego,"CC_enrichedOverAll")
make_ego_plots(ego2,"BP_enrichedOverAll")
# make_ego_plots(ego3,"CC_enrichedOverUni")
# make_ego_plots(ego4,"BP_enrichedOverUni")

# in https://pubmed.ncbi.nlm.nih.gov/37248224/ they use 1.2 fold enrichment and p < 0.05
df_enlo <- df[(df$so_ratioWave > 1.2) & (df$so_allTWave <= 0.05), ]

write.csv(df_enlo$ENTREZID, file = "Output/Data/df_enlo.csv", row.names = FALSE)
write.csv(df_en$ENTREZID, file = "Output/Data/df_en.csv", row.names = FALSE)

# KEGG

kk <- enrichKEGG(gene = unique(df_en$ENTREZID),
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
dotplot(kk, showCategory = 20, font.size=8)
cnetplot(kk, showCategory = 5, cex.params = list(category_label = 0.6, gene_label = 0.4), max.overlaps = 12)

kk2 <- enrichKEGG(gene = unique(df_en$ENTREZID),
                  universe = unique(df$ENTREZID),
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)
dotplot(kk2, showCategory = 20, font.size=8)
cnetplot(kk2, showCategory = 5, cex.params = list(category_label = 0.6, gene_label = 0.4), max.overlaps = 12)

kk3 <- enrichKEGG(gene = unique(df_en$UNIPROT),
                  universe = unique(df$UNIPROT),
                  organism     = 'hsa',
                  keyType = "uniprot",
                  pvalueCutoff = 0.05)
dotplot(kk3, showCategory = 20, font.size=8)
cnetplot(kk3, showCategory = 10, cex.params = list(category_label = 0.6, gene_label = 0.4), max.overlaps = 12)
