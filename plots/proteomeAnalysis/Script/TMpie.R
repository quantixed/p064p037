library(ggplot2)
pie_maker <-  function(x1, cols, x2 = NULL, label = "pie") {
  df <- data.frame(pvalues = x1,
                   pcolours = as.factor(1:length(cols)))
  cs <- cumsum(x1)
  df$ymin <- c(0,cs[1:length(cs)-1])
  df$ymax <- cs
  if(is.null(x2)) {
    df2 <- df
  } else {
    cs <- cumsum(x2)
    df2 <- data.frame(pvalues = x2,
                      ymin = c(0,cs[1:length(cs)-1]),
                      ymax = cs)
  }
  
  ggplot() +
    geom_rect(data = df, aes(fill = pcolours, ymin = ymin, ymax = ymax, xmin = 0, xmax = 3)) +
    scale_fill_manual(values = cols) +
    geom_rect(data = df2, aes(ymin = ymin, ymax = ymax, xmin = 0, xmax = 3), alpha = 0, colour = "black") +
    xlim(c(0,3)) +
    theme_minimal() +
    theme(aspect.ratio = 1) +
    coord_polar(theta = "y")
  ggsave(paste0("Output/Plots/", label,".pdf"), width = 5, height = 3, units = "cm", dpi = 300)
}


# pie_maker(x1 = c(634,3768-634), cols = c("#00a652","#ed8200"), label = "pie_coverage")
# pie_maker(x1 = c(176,458), cols = c("#008441","#22b56a"), label = "pie_tm_incidence_inv")
# # pie_maker(x1 = c(643,3125), cols = c("#a5003a","#e02a6b"), label = "pie_tm_incidence_ip")
# pie_maker(x1 = c(467,176,458,2667), cols = c("#bd6800","#008441","#22b56a","#ffa230"), x2 = c(643,3125), label = "pie_tm_incidence_ip")
# pie_maker(x1 = c(63,571), cols = c("#008441","#22b56a"), label = "pie_secreted_inv")
# # pie_maker(x1 = c(265,3503), cols = c("#bd6800","#ffa230"), label = "pie_secreted_ip")
# pie_maker(x1 = c(202,63,571,2932), cols = c("#bd6800","#008441","#22b56a","#ffa230"), x2 = c(265,3503), label = "pie_secreted_ip")

pie_maker(x1 = c(601,3631-601), cols = c("#009988","#bbbbbb"), label = "pie_coverage")
pie_maker(x1 = c(173,428), cols = c("#44aa99","#117733"), label = "pie_tm_incidence_inv")
pie_maker(x1 = c(636-173,173,428,2995-428), cols = c("#bbbbbb","#44aa99","#117733","#dddddd"), x2 = c(636,2995), label = "pie_tm_incidence_ip")
pie_maker(x1 = c(61,540), cols = c("#44aa99","#117733"), label = "pie_secreted_inv")
pie_maker(x1 = c(256-61,61,540,3375-540), cols = c("#bbbbbb","#44aa99","#117733","#dddddd"), x2 = c(256,3375), label = "pie_secreted_ip")
