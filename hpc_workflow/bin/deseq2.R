#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DESeq2))

parse_args <- function(args) {
  result <- list()
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (i == length(args)) stop(paste("Missing value for", args[[i]]))
    result[[key]] <- args[[i + 1]]
    i <- i + 2
  }
  result
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("counts", "samplesheet", "outdir", "design", "contrast", "min-count", "min-samples", "alpha")
missing <- required[!required %in% names(args)]
if (length(missing) > 0) stop(paste("Missing arguments:", paste(missing, collapse = ", ")))

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
counts_df <- read.delim(args$counts, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(counts_df) < 3) stop("Count matrix must contain a feature column and at least two samples.")
if (anyDuplicated(counts_df[[1]]) > 0) stop("Count matrix has duplicate feature IDs.")
rownames(counts_df) <- counts_df[[1]]
counts_df <- counts_df[, -1, drop = FALSE]
counts <- as.matrix(counts_df)
storage.mode(counts) <- "numeric"
if (any(!is.finite(counts)) || any(counts < 0)) stop("Counts must be finite and non-negative.")
counts <- round(counts)
storage.mode(counts) <- "integer"

metadata <- read.csv(args$samplesheet, check.names = FALSE, stringsAsFactors = FALSE)
if (!"sample" %in% colnames(metadata)) stop("Sample sheet requires a sample column.")
if (anyDuplicated(metadata$sample) > 0) stop("Sample sheet has duplicate sample IDs.")
rownames(metadata) <- metadata$sample
missing_samples <- setdiff(colnames(counts), rownames(metadata))
if (length(missing_samples) > 0) stop(paste("Samples absent from metadata:", paste(missing_samples, collapse = ", ")))
metadata <- metadata[colnames(counts), , drop = FALSE]

design_formula <- as.formula(args$design)
design_variables <- all.vars(design_formula)
missing_variables <- setdiff(design_variables, colnames(metadata))
if (length(missing_variables) > 0) stop(paste("Design variables absent from sample sheet:", paste(missing_variables, collapse = ", ")))
for (variable in design_variables) {
  if (is.character(metadata[[variable]])) metadata[[variable]] <- factor(metadata[[variable]])
  if (anyNA(metadata[[variable]])) stop(paste("Missing metadata in design variable", variable))
}

min_count <- as.integer(args[["min-count"]])
min_samples <- as.integer(args[["min-samples"]])
keep <- rowSums(counts >= min_count) >= min_samples
if (!any(keep)) stop("No features passed the count filter. Reduce --min_count or --min_samples.")
filtered_counts <- counts[keep, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(countData = filtered_counts, colData = metadata, design = design_formula)
dds <- DESeq(dds, quiet = TRUE)

contrast <- strsplit(args$contrast, ",", fixed = TRUE)[[1]]
if (length(contrast) != 3) stop("--contrast must be variable,numerator,denominator")
contrast <- trimws(contrast)
if (!contrast[[1]] %in% colnames(metadata)) stop("Contrast variable is absent from the sample sheet.")
levels_available <- levels(factor(metadata[[contrast[[1]]]]))
if (!all(contrast[2:3] %in% levels_available)) {
  stop(paste("Contrast levels must be in:", paste(levels_available, collapse = ", ")))
}

alpha <- as.numeric(args$alpha)
res <- results(dds, contrast = contrast, alpha = alpha)
result_df <- as.data.frame(res)
result_df$feature_id <- rownames(result_df)
result_df$significant <- !is.na(result_df$padj) & result_df$padj <= alpha
result_df$regulation <- ifelse(is.na(result_df$log2FoldChange), "NA", ifelse(result_df$log2FoldChange >= 0, "up", "down"))
result_df <- result_df[order(result_df$padj, result_df$pvalue, na.last = TRUE), ]
result_df <- result_df[, c("feature_id", setdiff(colnames(result_df), "feature_id"))]
write.table(result_df, file.path(args$outdir, "differential_expression.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

normalized <- counts(dds, normalized = TRUE)
write.table(cbind(feature_id = rownames(normalized), normalized), file.path(args$outdir, "normalized_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample = colnames(dds), size_factor = sizeFactors(dds)), file.path(args$outdir, "size_factors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
vst_counts <- assay(vsd)
write.table(cbind(feature_id = rownames(vst_counts), vst_counts), file.path(args$outdir, "vst_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(args$outdir, "MA_plot.pdf"), width = 7, height = 6)
plotMA(res, alpha = alpha, ylim = c(-5, 5))
dev.off()

if (ncol(vst_counts) >= 3) {
  pca <- prcomp(t(vst_counts), scale. = FALSE)
  percent <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
  pdf(file.path(args$outdir, "PCA_plot.pdf"), width = 7, height = 6)
  plot(pca$x[, 1], pca$x[, 2], xlab = sprintf("PC1 (%.1f%%)", percent[1]), ylab = sprintf("PC2 (%.1f%%)", percent[2]), pch = 19)
  text(pca$x[, 1], pca$x[, 2], labels = rownames(pca$x), pos = 3, cex = 0.7)
  dev.off()
}

sample_distance <- as.matrix(dist(t(vst_counts)))
write.table(cbind(sample = rownames(sample_distance), sample_distance), file.path(args$outdir, "sample_distances.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

parameters <- data.frame(
  parameter = c("design", "contrast", "min_count", "min_samples", "alpha", "features_input", "features_tested"),
  value = c(args$design, args$contrast, min_count, min_samples, alpha, nrow(counts), nrow(filtered_counts))
)
write.table(parameters, file.path(args$outdir, "analysis_parameters.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

sink(file.path(args$outdir, "session_info.txt"))
sessionInfo()
sink()
