# Download the raw IDAT files of GSE245203 from GEO into data/idats/.

suppressPackageStartupMessages(library(GEOquery))

gse <- "GSE245203"
out_dir <- "data/idats"
expected_idats <- 16  # 8 samples x 2 channels (Grn + Red)
max_tries <- 3

tar_path <- file.path("data", gse, paste0(gse, "_RAW.tar"))

# A download interrupted by a dropped connection or Ctrl+C leaves a truncated
# tar behind; listing its contents fails, so we treat it as missing.
tar_is_complete <- function(path) {
  file.exists(path) &&
    !inherits(tryCatch(untar(path, list = TRUE), error = identity, warning = identity), "condition")
}

# getGEOSuppFiles() creates data/<GSE>/ but not its parent: create it on a fresh clone
dir.create(dirname(tar_path), recursive = TRUE, showWarnings = FALSE)

if (file.exists(tar_path) && !tar_is_complete(tar_path)) {
  message("Removing incomplete download: ", tar_path)
  unlink(tar_path)
}

attempt <- 0
while (!tar_is_complete(tar_path) && attempt < max_tries) {
  attempt <- attempt + 1
  message("Downloading ", basename(tar_path), " (attempt ", attempt, "/", max_tries, ")")
  # GEOquery reports "No supplemental files found" both when the series has none
  # and when NCBI refuses the connection (e.g. HTTP 403 from one of its servers).
  tryCatch(
    getGEOSuppFiles(gse, baseDir = "data", filter_regex = "_RAW\\.tar$"),
    error = function(e) message("  failed: ", conditionMessage(e))
  )
  if (!tar_is_complete(tar_path)) unlink(tar_path)
}

if (!tar_is_complete(tar_path)) {
  stop("Could not download ", basename(tar_path), " from GEO after ", max_tries, " attempts. ",
       "If the series does have a RAW.tar, the NCBI server refused or dropped the connection; ",
       "retry later or test with:\n",
       "  curl -I https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/", gse, "/suppl/")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
untar(tar_path, exdir = out_dir)

n_idats <- length(list.files(out_dir, pattern = "\\.idat\\.gz$"))
if (n_idats != expected_idats) {
  stop("Expected ", expected_idats, " IDAT files, found ", n_idats)
}
message("Done: ", n_idats, " IDAT files in ", out_dir)
