# Step 2 — Normalization (stratified quantile normalization, minfi::preprocessQuantile)
# and probe filtering. Every filter is logged so the README can report how many
# probes each one removed.

source("scripts/_common.R")
suppressPackageStartupMessages({
  library(minfi)
  library(limma)
  library(ggplot2)
})

qc     <- readRDS(results_path("rds", "01_qc.rds"))
rgSet  <- qc$rgSet
detP   <- qc$detP
groups <- rgSet$Sample_Group

mSetRaw <- preprocessRaw(rgSet)
# preprocessQuantile estimates each sample's sex to normalize chrX/chrY separately.
# When all samples share one sex (here, a single cell line) it warns about an
# "inconsistency"; harmless, as sex-chromosome probes are removed below.
mSet    <- preprocessQuantile(rgSet)

# ---- Beta-value densities before and after normalization ---------------------
# A random subset of probes keeps the plot fast; the shape is identical.
set.seed(1)
probes <- sample(intersect(featureNames(mSet), featureNames(mSetRaw)), 50000)

density_data <- function(betas, label) {
  data.frame(
    beta   = as.vector(betas[probes, ]),
    Sample = rep(colnames(betas), each = length(probes)),
    Group  = rep(groups, each = length(probes)),
    Data   = label
  )
}
# offset = 100 (Illumina's default) avoids 0/0 betas for probes with no signal
dens <- rbind(density_data(getBeta(mSetRaw, offset = 100), "Raw"),
              density_data(getBeta(mSet), "Normalized"))
dens$Data <- factor(dens$Data, levels = c("Raw", "Normalized"))

p_dens <- ggplot(dens, aes(beta, group = Sample, colour = Group)) +
  geom_density(linewidth = 0.4) +
  facet_wrap(~Data) +
  labs(x = "Beta value", y = "Density", colour = NULL) +
  theme_minimal()
ggsave(results_path("figures", "02_beta_densities.png"), p_dens, width = 9, height = 4, dpi = 150)

# ---- Probe filtering ----------------------------------------------------------
filter_log <- data.frame(filter = "Probes after normalization", removed = NA_integer_,
                         remaining = nrow(mSet))

apply_filter <- function(mSet, keep, name) {
  filter_log <<- rbind(filter_log, data.frame(filter = name, removed = sum(!keep),
                                              remaining = sum(keep)))
  mSet[keep, ]
}

detP <- detP[match(featureNames(mSet), rownames(detP)), ]
mSet <- apply_filter(mSet, rowSums(detP < cfg$qc$probe_detp_max) == ncol(mSet),
                     sprintf("Detection p >= %s in any sample", cfg$qc$probe_detp_max))

mSet <- apply_filter(mSet, startsWith(featureNames(mSet), "cg"), "Non-CpG probes (ch.*)")

if (cfg$filtering$remove_sex_chromosomes) {
  mSet <- apply_filter(mSet, !as.character(seqnames(mSet)) %in% c("chrX", "chrY"),
                       "Sex chromosomes (chrX, chrY)")
}

if (cfg$filtering$remove_snp_probes) {
  no_snps <- featureNames(dropLociWithSnps(mSet))
  mSet <- apply_filter(mSet, featureNames(mSet) %in% no_snps,
                       "SNP at CpG or single-base extension site")
}

if (cfg$filtering$remove_cross_reactive) {
  if (cfg$array != "EPIC") {
    stop("reference/crosshyb_probes.txt covers 450K/EPIC v1 only; ",
         "set filtering$remove_cross_reactive to false for ", cfg$array)
  }
  crosshyb <- readLines("reference/crosshyb_probes.txt")
  mSet <- apply_filter(mSet, !featureNames(mSet) %in% crosshyb, "Cross-reactive probes")
}

write.csv(filter_log, results_path("tables", "02_probe_filtering.csv"), row.names = FALSE)
print(filter_log, row.names = FALSE)

# ---- Outputs for differential methylation --------------------------------------
beta <- getBeta(mSet)
M    <- getM(mSet)

# EPICv2 contains replicate probes for some CpGs (e.g. cg00000029_TC21);
# average them so every CpG appears once.
ann <- as.data.frame(getAnnotation(mSet))
ann <- ann[, c("chr", "pos", "Relation_to_Island", "UCSC_RefGene_Name", "UCSC_RefGene_Group")]

if (cfg$array == "EPICv2") {
  cpg_id <- sub("_[A-Z0-9]+$", "", rownames(M))
  M      <- avereps(M, ID = cpg_id)
  beta   <- avereps(beta, ID = cpg_id)
  ann    <- ann[!duplicated(cpg_id), ]
  rownames(ann) <- cpg_id[!duplicated(cpg_id)]
}

saveRDS(list(beta = beta, M = M, pheno = as.data.frame(colData(mSet)), annotation = ann),
        results_path("rds", "02_normalized.rds"))
message("Normalization done: ", nrow(M), " CpGs x ", ncol(M), " samples")
