# Step 4 — Summary figures: PCA, volcano plot and heatmap of significant DMPs.

source("scripts/_common.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexHeatmap)
})

d     <- readRDS(results_path("rds", "02_normalized.rds"))
dmp   <- readRDS(results_path("rds", "03_dmp_all.rds"))
pheno <- d$pheno
group_col <- cfg$dmp$group_column
block_col <- cfg$dmp$block_column

group_colours <- setNames(c("#E4572E", "#2E86AB"), sort(unique(pheno[[group_col]])))
# Okabe-Ito palette: distinguishable for colour-blind readers
block_levels  <- sort(unique(pheno[[block_col]]))
block_colours <- setNames(c("#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#F0E442",
                            "#0072B2", "#D55E00", "#999999")[seq_along(block_levels)],
                          block_levels)

# ---- PCA on the most variable CpGs --------------------------------------------
top_var <- order(apply(d$M, 1, var), decreasing = TRUE)[1:min(10000, nrow(d$M))]
pca     <- prcomp(t(d$M[top_var, ]), scale. = FALSE)
var_pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

pca_df <- data.frame(pca$x[, 1:2], Sample = rownames(pca$x),
                     Group = pheno[[group_col]], Mouse = pheno[[block_col]])

p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = Group, shape = Mouse, label = Sample)) +
  geom_point(size = 3.5) +
  geom_text(vjust = -1.2, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = group_colours) +
  labs(x = paste0("PC1 (", var_pct[1], "%)"), y = paste0("PC2 (", var_pct[2], "%)"),
       colour = NULL, title = "PCA — 10,000 most variable CpGs (M-values)") +
  theme_minimal()
ggsave(results_path("figures", "04_pca.png"), p_pca, width = 7, height = 5, dpi = 150)

# ---- Volcano plot --------------------------------------------------------------
dmp$status <- ifelse(dmp$significant, dmp$direction, "not significant")

p_volcano <- ggplot(dmp, aes(delta_beta, -log10(p_value), colour = status)) +
  geom_point(size = 0.4, alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1) * cfg$dmp$min_abs_delta_beta, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c(hyper = "#E4572E", hypo = "#2E86AB", "not significant" = "grey75")) +
  labs(x = paste0("Delta beta (", cfg$dmp$contrast, ")"), y = expression(-log[10](p)), colour = NULL,
       title = paste0(sum(dmp$significant), " DMPs (FDR < ", cfg$dmp$fdr_max,
                      ", |delta beta| >= ", cfg$dmp$min_abs_delta_beta, ")")) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_minimal()
ggsave(results_path("figures", "04_volcano.png"), p_volcano, width = 7, height = 5, dpi = 150)

# ---- Heatmap of significant DMPs ----------------------------------------------
sig_cpgs <- dmp$CpG[dmp$significant]
if (length(sig_cpgs) >= 2) {
  column_annotation <- HeatmapAnnotation(
    Group = pheno[[group_col]],
    Mouse = pheno[[block_col]],
    col   = setNames(list(group_colours, block_colours), c("Group", "Mouse"))
  )
  ht <- Heatmap(
    d$beta[sig_cpgs, ],
    name = "Beta",
    col = circlize::colorRamp2(c(0, 0.5, 1), c("#2E86AB", "white", "#E4572E")),
    top_annotation = column_annotation,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    clustering_distance_columns = "euclidean",
    column_title = paste(length(sig_cpgs), "significant DMPs")
  )
  png(results_path("figures", "04_heatmap.png"), width = 7, height = 7, units = "in", res = 150)
  draw(ht)
  invisible(dev.off())
}

message("Figures written to ", file.path(cfg$paths$results_dir, "figures"))
