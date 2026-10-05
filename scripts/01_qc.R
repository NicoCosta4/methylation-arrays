# Step 1 — Quality control: read IDATs, compute detection p-values and drop
# samples whose mean detection p-value is above the threshold in config.yml.

source("scripts/_common.R")
suppressPackageStartupMessages({
  library(minfi)
  library(ggplot2)
})

ss <- read.csv(cfg$paths$samplesheet, stringsAsFactors = FALSE)
ss$Basename <- file.path(cfg$paths$idat_dir, ss$Basename)

rgSet <- read.metharray.exp(targets = ss, verbose = TRUE)
rgSet@annotation <- array_annotation(cfg$array)
sampleNames(rgSet) <- ss$Sample_Name

detP <- detectionP(rgSet)

qc <- data.frame(
  Sample_Name       = colnames(detP),
  Sample_Group      = ss$Sample_Group,
  mean_detP         = colMeans(detP),
  failed_probes_pct = 100 * colMeans(detP > cfg$qc$probe_detp_max)
)
qc$pass <- qc$mean_detP < cfg$qc$sample_detp_max
write.csv(qc, results_path("tables", "01_qc_samples.csv"), row.names = FALSE)

p_detp <- ggplot(qc, aes(Sample_Name, mean_detP, fill = Sample_Group)) +
  geom_col() +
  geom_hline(yintercept = cfg$qc$sample_detp_max, colour = "red", linetype = "dashed") +
  labs(x = NULL, y = "Mean detection p-value", fill = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(results_path("figures", "01_detection_pvalues.png"), p_detp, width = 7, height = 4, dpi = 150)

# minfi intensity QC: samples below the dashed line ((mMed + uMed) / 2 < 10.5)
# have low overall signal and are usually failed hybridizations.
intensity <- as.data.frame(getQC(preprocessRaw(rgSet)))
intensity$Sample_Name  <- rownames(intensity)
intensity$Sample_Group <- ss$Sample_Group

p_int <- ggplot(intensity, aes(mMed, uMed, colour = Sample_Group, label = Sample_Name)) +
  geom_abline(intercept = 21, slope = -1, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3) +
  geom_text(vjust = -1, size = 3, show.legend = FALSE) +
  labs(x = "Median methylated intensity (log2)", y = "Median unmethylated intensity (log2)", colour = NULL) +
  theme_minimal()
ggsave(results_path("figures", "01_median_intensities.png"), p_int, width = 6, height = 5, dpi = 150)

if (any(!qc$pass)) {
  message("Removing samples that fail QC: ", paste(qc$Sample_Name[!qc$pass], collapse = ", "))
}
rgSet <- rgSet[, qc$pass]
detP  <- detP[, qc$pass]

saveRDS(list(rgSet = rgSet, detP = detP), results_path("rds", "01_qc.rds"))
message("QC done: ", sum(qc$pass), "/", nrow(qc), " samples kept")
