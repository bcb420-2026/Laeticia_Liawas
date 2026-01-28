
fetch_geo_supp <- function(gse, destdir = "data") {
  dir.create(destdir, showWarnings = FALSE, recursive = TRUE)
  GEOquery::getGEOSuppFiles(GEO = gse, baseDir = destdir, makeDirectory = TRUE)
  invisible(file.path(destdir, gse))
}
