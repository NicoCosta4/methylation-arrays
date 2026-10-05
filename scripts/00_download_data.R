# Download the raw IDAT files of GSE245203 from GEO into data/idats/.

suppressPackageStartupMessages(library(GEOquery))

gse <- "GSE245203"
out_dir <- "data/idats"
expected_idats <- 16  # 8 samples x 2 channels (Grn + Red)

tar_path <- file.path("data", gse, paste0(gse, "_RAW.tar"))
if (!file.exists(tar_path)) {
  downloaded <- getGEOSuppFiles(gse, baseDir = "data", filter_regex = "_RAW\\.tar$")
  # GEOquery reports "No supplemental files found" both when the series has none
  # and when NCBI refuses the connection (e.g. HTTP 403 for a blocked IP).
  if (is.null(downloaded) || !file.exists(tar_path)) {
    stop("Could not download ", basename(tar_path), " from GEO. If the series does ",
         "have a RAW.tar, NCBI may be refusing connections from your network; test with:\n",
         "  curl -I https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/", gse, "/suppl/")
  }
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
untar(tar_path, exdir = out_dir)

n_idats <- length(list.files(out_dir, pattern = "\\.idat\\.gz$"))
if (n_idats != expected_idats) {
  stop("Expected ", expected_idats, " IDAT files, found ", n_idats)
}
message("Done: ", n_idats, " IDAT files in ", out_dir)
