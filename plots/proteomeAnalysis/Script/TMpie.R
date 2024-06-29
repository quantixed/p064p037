library(ggplot2)

## Functions ----

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

# Values here are hardcoded and taken from tmCalc.txt

# pie to show INV as a fraction of total
pie_maker(x1 = c(601, 3622 - 601),
          cols = c("#009988", "#bbbbbb"), label = "pie_coverage")
# TM incidence in INVs
pie_maker(x1 = c(168, 433),
          cols = c("#44aa99", "#117733"), label = "pie_tm_incidence_inv")
# TM incidence in all and in INVs
pie_maker(x1 = c(617 - 168, 168, 433, 3005 - 433),
          cols = c("#bbbbbb", "#44aa99", "#117733", "#dddddd"),
          x2 = c(636, 2995), label = "pie_tm_incidence_ip")
# secreted incidence in INVs
pie_maker(x1 = c(63, 538),
          cols = c("#44aa99", "#117733"), label = "pie_secreted_inv")
# secreted incidence in all and in INVs
pie_maker(x1 = c(259 - 63, 63, 538, 3363 - 538),
          cols = c("#bbbbbb", "#44aa99", "#117733", "#dddddd"),
          x2 = c(256, 3375), label = "pie_secreted_ip")
