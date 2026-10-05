# Step 3 — Differentially methylated positions (DMPs) with limma.
# Statistics are computed on M-values (better variance properties); effect sizes
# are reported as delta-beta, which is easier to interpret. Samples from the same
# mouse are correlated, so the mouse enters the model as a blocking factor.

source("scripts/_common.R")
suppressPackageStartupMessages(library(limma))

d     <- readRDS(results_path("rds", "02_normalized.rds"))
group <- factor(d$pheno[[cfg$dmp$group_column]])
block <- d$pheno[[cfg$dmp$block_column]]

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
contrast <- makeContrasts(contrasts = cfg$dmp$contrast, levels = design)

# Only a single consensus correlation is needed, so estimate it on a random subset
# of CpGs: same value in a fraction of the time of using all ~800k probes.
set.seed(1)
cor_rows <- sample(nrow(d$M), min(50000, nrow(d$M)))
corfit   <- duplicateCorrelation(d$M[cor_rows, ], design, block = block)
message("Within-mouse correlation (consensus): ", round(corfit$consensus.correlation, 3))

fit <- lmFit(d$M, design, block = block, correlation = corfit$consensus.correlation)
fit <- eBayes(contrasts.fit(fit, contrast))
res <- topTable(fit, number = Inf, sort.by = "none")

# Delta-beta: the same contrast applied to the per-group mean beta values
group_means <- sapply(levels(group), function(g) rowMeans(d$beta[, group == g, drop = FALSE]))
res$delta_beta <- as.vector(group_means %*% contrast[levels(group), 1])

# Promoter / gene body / intergenic, from the UCSC RefGene annotation
gene_region <- function(refgene_group) {
  parts <- strsplit(ifelse(is.na(refgene_group), "", refgene_group), ";")
  promoter_terms <- c("TSS1500", "TSS200", "5'UTR", "5UTR", "1stExon", "exon_1")
  vapply(parts, function(p) {
    if (length(p) == 0) return("intergenic")
    if (any(p %in% promoter_terms)) "promoter" else "body"
  }, character(1))
}

# RefGene lists one gene name per transcript ("RBP1;RBP1;RBP1"); keep unique names
unique_genes <- function(x) {
  vapply(strsplit(ifelse(is.na(x), "", x), ";"),
         function(g) paste(unique(g), collapse = ";"), character(1))
}

ann <- d$annotation[rownames(res), ]
dmp <- data.frame(
  CpG             = rownames(res),
  chr             = ann$chr,
  pos             = ann$pos,
  gene            = unique_genes(ann$UCSC_RefGene_Name),
  gene_region     = gene_region(ann$UCSC_RefGene_Group),
  island_relation = ann$Relation_to_Island,
  round(group_means, 3),
  delta_beta      = round(res$delta_beta, 3),
  logFC_M         = round(res$logFC, 3),
  t               = round(res$t, 3),
  p_value         = signif(res$P.Value, 3),
  fdr             = signif(res$adj.P.Val, 3)
)
dmp <- dmp[order(dmp$p_value), ]

dmp$significant <- dmp$fdr < cfg$dmp$fdr_max & abs(dmp$delta_beta) >= cfg$dmp$min_abs_delta_beta
dmp$direction   <- ifelse(dmp$delta_beta > 0, "hyper", "hypo")
sig <- dmp[dmp$significant, ]

summary_table <- as.data.frame.matrix(table(sig$gene_region, sig$direction))
write.csv(summary_table, results_path("tables", "03_dmp_summary.csv"))
write.csv(sig, results_path("tables", "03_dmp_significant.csv"), row.names = FALSE)
saveRDS(dmp, results_path("rds", "03_dmp_all.rds"))

message(nrow(sig), " DMPs (FDR < ", cfg$dmp$fdr_max, ", |delta-beta| >= ",
        cfg$dmp$min_abs_delta_beta, "): ", sum(sig$direction == "hyper"), " hyper, ",
        sum(sig$direction == "hypo"), " hypo in ", cfg$dmp$contrast)
