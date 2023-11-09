# the goal of this script is to visualise Daniel's dataset
library(ggplot2)
library(ggrepel)

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


# script ----
wt_ct <- read.delim("Data/rankTable_WTvsControl.txt")
wt_re <- read.delim("Data/rankTable_WTvsR159E.txt")
re_ct <- read.delim("Data/rankTable_R159EvsControl.txt")

all_df <- merge(wt_ct,wt_re, by = "so_SHORTNAME", sort = FALSE)
all_df <- merge(all_df,re_ct, by = "so_SHORTNAME", sort = FALSE)
# after this names are messed up
keep <- c("so_SHORTNAME","so_NAME.x","so_PID.x","so_productWave.x","so_colorWave.x","so_allTWave.x","so_ratioWave.x",
          "so_productWave.y","so_colorWave.y","so_allTWave.y","so_ratioWave.y",
          "so_productWave","so_colorWave","so_allTWave","so_ratioWave")
all_df <- all_df[,names(all_df) %in% keep]
names(all_df) <- c("SHORTNAME","NAME","PID","productWave_wt_ct","colorWave_wt_ct","allTWave_wt_ct","ratioWave_wt_ct",
                   "productWave_wt_re","colorWave_wt_re","allTWave_wt_re","ratioWave_wt_re",
                   "productWave_re_ct","colorWave_re_ct","allTWave_re_ct","ratioWave_re_ct")
# some names to be disambiguated, e.g. RAB8A;RAB8B from RAB8A
all_df$uSHORTNAME <- make_shortnames_unique(all_df$SHORTNAME)

# coolors <- c("#403075",
#              "#822b66",
#              "#aa3939",
#              "#2c4770",
#              "#808080",
#              "#aa1839",
#              "#2d882d",
#              "#83a136",
#              "#aa9739")

coolors <- rev(c("#4122a6","#712771","#f11919",
             "#2c5270","#808080","#ae7d3c",
             "#14c114","#95a83a","#f1cc19"))

# bit0 is significant in 1st set, bit1 is significant in 2nd
all_df$pcat1 <- ifelse(all_df$allTWave_wt_ct < 0.05,1,0) + ifelse(all_df$allTWave_wt_re < 0.05,2,0)
# now we have 3 levels (up, mid, down) for two coordinates, 9 levels, 0-8.
all_df$encat1 <- ifelse(all_df$ratioWave_wt_ct > 2,2,0) +
  ifelse((all_df$ratioWave_wt_ct < 2) & (all_df$ratioWave_wt_ct > 0.5),1,0) +
  ifelse(all_df$ratioWave_wt_re > 2,6,0) +
  ifelse((all_df$ratioWave_wt_re < 2) & (all_df$ratioWave_wt_re > 0.5),3,0)
# labels for only the proteins of interest
all_df$label1 <- ifelse((all_df$pcat1 == 3) & (all_df$encat1 == 8), all_df$uSHORTNAME, "")


p1 <- ggplot(all_df, aes(x = log2(ratioWave_wt_ct), y = log2(ratioWave_wt_re), label = label1)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") +
  geom_point(aes(colour = as.factor(encat1), size = ceiling(pcat1 / 4) / 10, shape = as.factor(pcat1)), alpha = 0.75) +
  geom_text_repel(size = 1.8, max.overlaps = 25) +
  scale_shape_manual(values = c(1,0,5,16)) +
  scale_color_manual(values = coolors) +
  scale_x_continuous(limits = c(-2.5,10), breaks=seq(-2,10,2)) +
  scale_y_continuous(limits = c(-2.5,8), breaks=seq(-2,8,2)) +
  coord_fixed() +
  labs(x = bquote(WT~-~Control~(Log[2])), y = bquote(WT~-~R159E~(Log[2]))) +
  theme_classic(8) +
  scale_size(range = c(0.5,1.5)) +
  theme(legend.position = "none")

ggsave("Output/Plots/WTvREvsWTvCT.pdf",p1, width = 80, height = 80, units = "mm")
ggsave("Output/Plots/WTvREvsWTvCTKey.pdf",p1+theme(legend.position = "right"), width = 120, height = 120, units = "mm")


# bit0 is significant in 1st set, bit1 is significant in 2nd
all_df$pcat2 <- ifelse(all_df$allTWave_wt_ct < 0.05,1,0) + ifelse(all_df$allTWave_re_ct < 0.05,2,0)
# now we have 3 levels (up, mid, down) for two coordinates, 9 levels, 0-8.
all_df$encat2 <- ifelse(all_df$ratioWave_wt_ct > 2,2,0) +
  ifelse((all_df$ratioWave_wt_ct < 2) & (all_df$ratioWave_wt_ct > 0.5),1,0) +
  ifelse(all_df$ratioWave_re_ct > 2,6,0) +
  ifelse((all_df$ratioWave_re_ct < 2) & (all_df$ratioWave_re_ct > 0.5),3,0)
# labels for only the proteins of interest
all_df$label2 <- ifelse((all_df$pcat2 == 3) & (all_df$encat2 == 8), all_df$uSHORTNAME, "")


p2 <- ggplot(all_df, aes(x = log2(ratioWave_wt_ct), y = log2(ratioWave_re_ct), label = label2)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") +
  geom_point(aes(colour = as.factor(encat2), size = ceiling(pcat2 / 4) / 10, shape = as.factor(pcat2)), alpha = 0.75) +
  geom_text_repel(size = 1.8, max.overlaps = 25) +
  scale_shape_manual(values = c(1,0,5,16)) +
  scale_color_manual(values = coolors) +
  scale_x_continuous(limits = c(-2.5,10), breaks=seq(-2,10,2)) +
  scale_y_continuous(limits = c(-2.5,8), breaks=seq(-2,8,2)) +
  coord_fixed() +
  labs(x = bquote(WT~-~Control~(Log[2])), y = bquote(R159E~-~Control~(Log[2]))) +
  theme_classic(8) +
  scale_size(range = c(0.5,1.5)) +
  theme(legend.position = "none")

ggsave("Output/Plots/WTvREvsREvCT.pdf",p2, width = 80, height = 80, units = "mm")

p3 <- ggplot(all_df, aes(x = log2(ratioWave_wt_ct), y = log2(ratioWave_re_ct), label = label1)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") +
  geom_point(aes(colour = as.factor(encat2), size = ceiling(pcat2 / 4) / 10, shape = as.factor(pcat2)), alpha = 0.75) +
  geom_text_repel(size = 1.8, max.overlaps = 25) +
  scale_shape_manual(values = c(1,0,5,16)) +
  scale_color_manual(values = coolors) +
  scale_x_continuous(limits = c(-2.5,10), breaks=seq(-2,10,2)) +
  scale_y_continuous(limits = c(-2.5,8), breaks=seq(-2,8,2)) +
  coord_fixed() +
  labs(x = bquote(WT~-~Control~(Log[2])), y = bquote(R159E~-~Control~(Log[2]))) +
  theme_classic(8) +
  scale_size(range = c(0.5,1.5)) +
  theme(legend.position = "none")

ggsave("Output/Plots/WTvREvsREvCTalt.pdf",p3, width = 80, height = 80, units = "mm")

# volcano plots ----

vp_label <- function(x) {
  ifelse(
    x == 0, "1",
    parse(text = gsub("1e", "10^", scales::label_scientific()(10^-x)))
  )
}


make_volcano_plot <- function(df, vs, aa, bb, poi = NULL) {
  xw <- sym(paste0("ratioWave_",vs))
  yw <- sym(paste0("allTWave_",vs))
  zw <- sym(paste0("colorWave_",vs))
  # try to label all proteins of interest
  if(is.null(poi)) {
    df$interest <- ifelse(df[[zw]] == 3, df$uSHORTNAME, "")
  } else {
    df$interest <- ifelse(df[[zw]] == 3 & df$uSHORTNAME %in% poi, df$uSHORTNAME, "")
  }
  
  p <- ggplot(df, aes(x = log2(!!xw), y = -log10(!!yw), label = interest)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey") +
    geom_point(aes(colour = as.factor(!!zw)), size = 1, shape = 16, alpha = 0.5) +
    scale_colour_manual(values = c("0" = "#808080", "1" = "#ff8080", "2" = "#8080ff", "3" = "#ff80ff")) +
    geom_text_repel(size = 1.5, max.overlaps = 15, segment.color = "#7f7f7f7f") +
    scale_x_continuous(limits = c(-3,10), breaks = seq(-2,10,2)) +
    scale_y_continuous(limits = c(0,12), breaks = seq(0,12,2), labels = vp_label) +
    xlab(bquote(.(aa)~-~.(bb)~(Log[2]))) +
    ylab("P-value") +
    theme_classic(8) +
    theme(legend.position = "none")
  ggsave(paste0("Output/Plots/vp_",vs,".pdf"), p, width = 48, height = 80, units = "mm", family = "Helvetica")
  return(p)
}

toMatch <- c("^TPD5","^RAB","^VAMP","ATG9","SCAMP","DAGLB","IGF2R","PI4K2A","SERINC","SYNGR", "SYPL")
pois <- unique(grep(paste(toMatch,collapse="|"), all_df$uSHORTNAME, value = TRUE))

p4 <- make_volcano_plot(all_df,"wt_ct","WT","Control", poi = pois)
p5 <- make_volcano_plot(all_df,"wt_re","WT","R159E", poi = pois)
p6 <- make_volcano_plot(all_df,"re_ct","R159E","Control", poi = pois)
  
