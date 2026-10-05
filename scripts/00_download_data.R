# Download the raw IDAT files of GSE245203 from GEO into data/idats/.

suppressPackageStartupMessages(library(GEOquery))

gse <- "GSE245203"
out_dir <- "data/idats"
expected_idats <- 16  # 8 samples x 2 channels (Grn + Red)

tar_path <- file.path("data", gse, paste0(gse, "_RAW.tar"))
if (!file.exists(tar_path)) {
  getGEOSuppFiles(gse, baseDir = "data", filter_regex = "_RAW\\.tar$")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
untar(tar_path, exdir = out_dir)

n_idats <- length(list.files(out_dir, pattern = "\\.idat\\.gz$"))
if (n_idats != expected_idats) {
  stop("Expected ", expected_idats, " IDAT files, found ", n_idats)
}
message("Done: ", n_idats, " IDAT files in ", out_dir)
