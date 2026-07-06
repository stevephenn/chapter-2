# =============================================================================
# 01_MG_MT_organism_DESeq2_clean.R
# Reproducible organism-level DESeq2 analysis for MG and MT.
# Gene filter: >=5 counts in >=2 samples.
# Organism filter: >=10 total counts across 16 Low/High samples.
# DESeq2 contrast: High vs Low; positive log2FC = enriched in High/HME.
# =============================================================================

# ---- Packages ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2", ask = FALSE, update = FALSE)
for (pkg in c("data.table", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(ggplot2)
})

# ---- Configuration ----
project_dir <- "/Users/stephencourtney/Library/Mobile Documents/com~apple~CloudDocs/analysis_folder/"
setwd(project_dir)

file_path_MG <- "counts-HL_combined2_featureCounts.txt"
file_path_MT <- "hungate_alignment_primarybam_all.Fc_copy_2.txt"  # amend only if renamed

out_dir <- file.path(project_dir, "results_clean")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_gene_count <- 5
min_gene_samples <- 2
min_org_total <- 10
alpha <- 0.05

# MG sample IDs
low_MG <- c("SRR873595","SRR873596","SRR873597","SRR873598",
            "SRR873599","SRR873600","SRR873601","SRR873602")
high_MG <- c("SRR1206671","SRR873604","SRR873605","SRR873606",
             "SRR873607","SRR873608","SRR873609","SRR873610")

# MT sample IDs
low_MT <- c("SRR873450","SRR873451","SRR873452","SRR873453",
            "SRR873454","SRR873455","SRR873456","SRR873457")
high_MT <- c("SRR1206249","SRR873459","SRR873460","SRR873461",
             "SRR873462","SRR873463","SRR873464","SRR873465")

# ---- Helper functions ----
read_featurecounts <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  lines <- readLines(path, warn = FALSE)
  header_i <- grep("^Geneid\\t", lines)[1]
  if (is.na(header_i)) stop("Could not find featureCounts Geneid header in: ", path)

  fc <- data.table::fread(path, skip = header_i - 1, header = TRUE,
                          sep = "\\t", data.table = FALSE, check.names = FALSE)
  req <- c("Geneid","Chr","Start","End","Strand","Length")
  if (!all(req %in% names(fc))) stop("Missing required featureCounts columns.")
  fc
}

extract_srr <- function(x) {
  m <- regexpr("SRR[0-9]+", x)
  ifelse(m > 0, regmatches(x, m), NA_character_)
}

collapse_duplicate_sample_columns <- function(mat) {
  base <- sub("\\.[0-9]+$", "", colnames(mat))
  out <- t(rowsum(t(mat), group = base, reorder = FALSE))
  storage.mode(out) <- "integer"
  out
}

calc_org_tpm <- function(org_counts, org_length_bp) {
  org_length_kb <- org_length_bp / 1000
  if (any(!is.finite(org_length_kb)) || any(org_length_kb <= 0)) {
    stop("Invalid effective organism length during TPM calculation.")
  }
  rpk <- sweep(org_counts, 1, org_length_kb, "/")
  denom <- colSums(rpk)
  if (any(denom <= 0)) stop("At least one sample has zero RPK total.")
  sweep(rpk, 2, denom, "/") * 1e6
}

run_org_deseq2 <- function(file_path, low_ids, high_ids, label,
                            min_gene_count = 5, min_gene_samples = 2,
                            min_org_total = 10, alpha = 0.05) {

  cat("\\n--- Running ", label, " ---\\n", sep = "")
  fc <- read_featurecounts(file_path)
  annotation_cols <- c("Geneid","Chr","Start","End","Strand","Length")

  # Rename count columns to SRR identifiers and build a gene x sample matrix.
  sample_cols <- setdiff(names(fc), annotation_cols)
  srr <- extract_srr(sample_cols)
  if (anyNA(srr)) stop(label, ": unable to extract SRR ID from one or more sample headers.")
  names(fc) <- c(annotation_cols, make.unique(srr))

  gene_counts <- as.matrix(fc[, -(1:6), drop = FALSE])
  storage.mode(gene_counts) <- "integer"
  rownames(gene_counts) <- fc$Geneid
  gene_counts <- collapse_duplicate_sample_columns(gene_counts)

  # Retain exactly the 8 Low and 8 High samples in declared order.
  keep_samples <- c(low_ids, high_ids)
  missing <- setdiff(keep_samples, colnames(gene_counts))
  if (length(missing) > 0) stop(label, ": missing sample IDs: ", paste(missing, collapse = ", "))
  gene_counts_LH <- gene_counts[, keep_samples, drop = FALSE]

  condition <- ifelse(colnames(gene_counts_LH) %in% low_ids, "Low", "High")
  coldata <- data.frame(condition = factor(condition, levels = c("Low", "High")),
                        row.names = colnames(gene_counts_LH))
  stopifnot(identical(rownames(coldata), colnames(gene_counts_LH)))

  # Organism is encoded in Chr before '#'.
  organism <- trimws(sub("#.*", "", fc$Chr))
  gene_length_bp <- as.numeric(fc$Length)
  stopifnot(length(organism) == nrow(gene_counts_LH),
            length(gene_length_bp) == nrow(gene_counts_LH))
  if (anyNA(organism) || any(organism == "")) stop(label, ": missing organism assignment.")
  if (anyNA(gene_length_bp) || any(gene_length_bp <= 0)) stop(label, ": invalid gene length.")

  # Gene filter: >=5 counts in >=2 samples.
  keep_gene <- rowSums(gene_counts_LH >= min_gene_count) >= min_gene_samples
  if (!any(keep_gene)) stop(label, ": no genes passed the 5-in-2 gene filter.")

  gene_counts_filt <- gene_counts_LH[keep_gene, , drop = FALSE]
  organism_filt <- organism[keep_gene]
  lengths_filt <- gene_length_bp[keep_gene]

  # Organism aggregation and >=10 total-count filter.
  org_counts_pre10 <- rowsum(gene_counts_filt, group = organism_filt, reorder = TRUE)
  storage.mode(org_counts_pre10) <- "integer"
  org_lengths_pre10 <- tapply(lengths_filt, organism_filt, sum)
  org_lengths_pre10 <- org_lengths_pre10[rownames(org_counts_pre10)]

  keep_org <- rowSums(org_counts_pre10) >= min_org_total
  org_counts <- org_counts_pre10[keep_org, , drop = FALSE]
  org_lengths <- org_lengths_pre10[rownames(org_counts)]
  storage.mode(org_counts) <- "integer"
  stopifnot(identical(rownames(coldata), colnames(org_counts)))

  # DESeq2. DESeq2 handles BH-adjustment internally in results().
  dds <- DESeqDataSetFromMatrix(countData = org_counts, colData = coldata,
                                design = ~ condition)
  dds <- DESeq(dds, sfType = "poscounts", quiet = TRUE)
  res <- results(dds, contrast = c("condition", "High", "Low"), alpha = alpha)

  res_tbl <- as.data.frame(res)
  res_tbl$organism <- rownames(res_tbl)
  res_tbl <- res_tbl[, c("organism", setdiff(names(res_tbl), "organism"))]
  res_tbl <- res_tbl[order(is.na(res_tbl$padj), res_tbl$padj, res_tbl$pvalue), ]

  tpm <- calc_org_tpm(org_counts, org_lengths)

  qc <- data.frame(
    sample = colnames(org_counts),
    condition = coldata$condition,
    organism_library_size = colSums(org_counts),
    deseq2_size_factor = sizeFactors(dds),
    row.names = NULL
  )

  cat("Genes before filter: ", nrow(gene_counts_LH), "\\n", sep = "")
  cat("Genes after 5-in-2 filter: ", nrow(gene_counts_filt), "\\n", sep = "")
  cat("Organisms after >=10 total-count filter: ", nrow(org_counts), "\\n", sep = "")
  cat("Significant organisms (padj < 0.05): ",
      sum(!is.na(res_tbl$padj) & res_tbl$padj < alpha), "\\n", sep = "")

  list(label = label, coldata = coldata, org_counts = org_counts,
       org_lengths = org_lengths, tpm = tpm, dds = dds,
       results = res, results_table = res_tbl, qc = qc)
}

save_outputs <- function(x, out_dir) {
  prefix <- x$label

  write.table(x$results_table,
              file.path(out_dir, paste0(prefix, "_organism_DESeq2_High_vs_Low.tsv")),
              sep = "\\t", quote = FALSE, row.names = FALSE)

  write.table(x$tpm,
              file.path(out_dir, paste0(prefix, "_organism_TPM.tsv")),
              sep = "\\t", quote = FALSE, col.names = NA)

  write.table(x$qc,
              file.path(out_dir, paste0(prefix, "_QC.tsv")),
              sep = "\\t", quote = FALSE, row.names = FALSE)

  # PCA and dispersion QC figure.
  pdf(file.path(out_dir, paste0(prefix, "_DESeq2_QC.pdf")), width = 8, height = 6)
  vsd <- varianceStabilizingTransformation(x$dds, blind = FALSE)
  print(plotPCA(vsd, intgroup = "condition") + ggtitle(paste0(prefix, " VST PCA")))
  plotDispEsts(x$dds, main = paste0(prefix, " dispersion estimates"))
  dev.off()

  saveRDS(list(coldata = x$coldata, org_counts = x$org_counts,
               org_lengths = x$org_lengths, tpm = x$tpm, dds = x$dds,
               results = x$results, results_table = x$results_table, qc = x$qc),
          file.path(out_dir, paste0(prefix, "_final_objects.rds")))
}

# ---- Run MG and MT ----------------------------------------------------------
MG_final <- run_org_deseq2(file_path_MG, low_MG, high_MG, "MG",
                            min_gene_count, min_gene_samples, min_org_total, alpha)
MT_final <- run_org_deseq2(file_path_MT, low_MT, high_MT, "MT",
                            min_gene_count, min_gene_samples, min_org_total, alpha)

save_outputs(MG_final, out_dir)
save_outputs(MT_final, out_dir)

# ---- Cross-omic TPM correlations --------------------------------------------
# MG and MT samples are unpaired. Therefore correlate across shared organisms.
shared_org <- intersect(rownames(MG_final$tpm), rownames(MT_final$tpm))
MG_tpm <- MG_final$tpm[shared_org, , drop = FALSE]
MT_tpm <- MT_final$tpm[shared_org, , drop = FALSE]

mg_low_cols <- rownames(MG_final$coldata)[MG_final$coldata$condition == "Low"]
mg_high_cols <- rownames(MG_final$coldata)[MG_final$coldata$condition == "High"]
mt_low_cols <- rownames(MT_final$coldata)[MT_final$coldata$condition == "Low"]
mt_high_cols <- rownames(MT_final$coldata)[MT_final$coldata$condition == "High"]

mg_mean <- rowMeans(log2(MG_tpm + 1))
mt_mean <- rowMeans(log2(MT_tpm + 1))
mg_effect <- rowMeans(log2(MG_tpm[, mg_high_cols, drop = FALSE] + 1)) -
             rowMeans(log2(MG_tpm[, mg_low_cols, drop = FALSE] + 1))
mt_effect <- rowMeans(log2(MT_tpm[, mt_high_cols, drop = FALSE] + 1)) -
             rowMeans(log2(MT_tpm[, mt_low_cols, drop = FALSE] + 1))

ct_mean <- cor.test(mg_mean, mt_mean, method = "spearman", exact = FALSE)
ct_effect <- cor.test(mg_effect, mt_effect, method = "spearman", exact = FALSE)

cor_summary <- data.frame(
  comparison = c("Overall mean log2(TPM + 1)", "High minus Low log2(TPM + 1) effect"),
  n_shared_organisms = length(shared_org),
  spearman_rho = c(unname(ct_mean$estimate), unname(ct_effect$estimate)),
  p_value = c(ct_mean$p.value, ct_effect$p.value)
)

write.table(cor_summary, file.path(out_dir, "MG_MT_TPM_Spearman_correlations.tsv"),
            sep = "\\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

cat("\\nComplete. Outputs written to: ", out_dir, "\\n", sep = "")
