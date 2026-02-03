fetch_geo_supp <- function(gse, destdir = "data") {
  dir.create(destdir, showWarnings = FALSE, recursive = TRUE)
  GEOquery::getGEOSuppFiles(GEO = gse, baseDir = destdir, makeDirectory = TRUE)
  invisible(file.path(destdir, gse))
}


kable_head <- function(x, n = 5, caption = NULL) {
  knitr::kable(utils::head(x, n), caption = caption)
}

plot_box <- function(mat, main = "", ylab = "log2(counts+1)") {
  df <- as.data.frame(mat)
  df_long <- df |>
    mutate(gene = rownames(df)) |>
    pivot_longer(-gene, names_to = "sample", values_to = "value") |>
    mutate(value = log2(value + 1))
  
  ggplot(df_long, aes(x = sample, y = value)) +
    geom_boxplot(outlier.size = 0.2) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = main, x = "Sample", y = ylab)
}

plot_density <- function(mat, main = "") {
  df <- as.data.frame(mat)
  df_long <- df |>
    mutate(gene = rownames(df)) |>
    pivot_longer(-gene, names_to = "sample", values_to = "value") |>
    mutate(value = log2(value + 1))
  
  ggplot(df_long, aes(x = value, colour = sample)) +
    geom_density() +
    theme_bw() +
    labs(title = main, x = "log2(counts+1)", y = "Density") +
    guides(colour = "none")
}
