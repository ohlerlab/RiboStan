#!/usr/bin/env Rscript
# Compare translation-supported Ribostan uORFs vs RiboTaper periodic ORFs.
args <- commandArgs(trailingOnly = TRUE)
workdir <- if (length(args) >= 1) args[[1]] else "."
setwd(workdir)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
})

rs <- readRDS("ribostan_uorfs.rds")
if (nrow(rs) == 0) {
  warning("No periodic Ribostan uORFs in ribostan_uorfs.rds")
}

rt_files <- list.files(workdir, pattern = "sORFs_sign|ORFs_sign", recursive = TRUE, full.names = TRUE)
rt_files <- rt_files[!grepl("sign_filtered_multi$|sign_filtered_cds_multi$|/ORFs_all$", rt_files)]
rt_files <- rt_files[file.info(rt_files)$size > 200]
message("RiboTaper candidate files: ", length(rt_files))

read_rt_table <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(NULL)
  df <- tryCatch(
    read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = ""),
    error = function(e) NULL
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)
  if (all(is.na(df[1, ]))) return(NULL)
  df
}

prio <- rt_files[grepl("best_periodicity/ORFs_sign_notfiltered_multi$", rt_files)]
if (length(prio) == 0) {
  prio <- rt_files[grepl("ORFs_sign_notfiltered_multi$", rt_files)]
}
if (length(prio) == 0) {
  prio <- rt_files[grepl("ORFs_sign_nocds_nofilter$", rt_files)]
}
if (length(prio) == 0) prio <- rt_files

rt_list <- lapply(prio, read_rt_table)
rt_list <- rt_list[!vapply(rt_list, is.null, logical(1))]

if (length(rt_list) == 0) {
  cat("ERROR: no usable RiboTaper ORF tables under ", workdir, "\n", sep = "")
  quit(status = 1)
}

rt <- do.call(rbind, rt_list)
for (f in unique(prio)) message("using: ", f)

# RiboTaper ORF_length is track length in nucleotides minus 1, not codons.
rt_df <- data.frame(
  transcript_id = rt$transcript_id,
  orf_id_tr = rt$ORF_id_tr,
  strand = rt$strand,
  width_nt = as.integer(rt$ORF_length) + 1L,
  width_codons = as.integer(rt$ORF_length) / 3L,
  p_value = if ("ORF_pval_multi_ribo" %in% names(rt)) rt$ORF_pval_multi_ribo else NA_real_,
  periodic = if ("ORF_pval_multi_ribo" %in% names(rt)) rt$ORF_pval_multi_ribo < 0.05 else FALSE,
  category = if ("type" %in% names(rt)) rt$type else "unknown",
  stringsAsFactors = FALSE
)
rt_df <- rt_df[!is.na(rt_df$transcript_id) & rt_df$transcript_id != "NA", , drop = FALSE]

# Translation-supported only: periodic signal and RiboTaper minimum track length (>5 nt).
rt_u <- rt_df[rt_df$periodic & rt_df$width_nt > 5L, , drop = FALSE]
if (any(grepl("uORF", rt_df$category, ignore.case = TRUE), na.rm = TRUE)) {
  rt_u <- rt_u[grepl("uORF", rt_u$category, ignore.case = TRUE), , drop = FALSE]
}

saveRDS(rt_u, "ribotaper_uorfs.rds")
write.table(rt_u, "ribotaper_uorfs.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

n_candidates <- NA_integer_
cand_path <- "ribostan_uorfs_candidates.tsv"
if (file.exists(cand_path)) {
  n_candidates <- nrow(read.delim(cand_path, stringsAsFactors = FALSE))
}

rs$key <- paste(rs$transcript_id, rs$start, rs$end, sep = "|")
rt$key <- paste(rt_u$transcript_id, rt_u$orf_id_tr, sep = "|")

both <- intersect(rs$key, rt$key)
rs_only <- setdiff(rs$key, rt$key)
rt_only <- setdiff(rt$key, rs$key)

summary_lines <- c(
  "=== Translation-supported uORF comparison (chr22 vignette data) ===",
  "",
  "Both sides require triplet-periodicity evidence (Ribostan q<0.05;",
  "RiboTaper multitaper p<0.05). Sequence-only candidates are excluded.",
  "",
  sprintf("Ribostan sequence candidates (diagnostic only): %s",
          if (is.na(n_candidates)) "n/a" else as.character(n_candidates)),
  sprintf("Ribostan periodic uORFs (q < 0.05): %d", nrow(rs)),
  sprintf("RiboTaper periodic ORFs (p < 0.05, width > 5 nt): %d", nrow(rt_u)),
  sprintf("Exact match (tx + coords): %d", length(both)),
  sprintf("Ribostan periodic only: %d", length(rs_only)),
  sprintf("RiboTaper periodic only: %d", length(rt_only)),
  "",
  "--- Ribostan periodic uORF lengths (after stop-codon strip) ---",
  sprintf("  width 3bp (1 codon): %d", sum(rs$width_nt == 3)),
  sprintf("  width 6bp (2 codons): %d", sum(rs$width_nt == 6)),
  sprintf("  width <= 9bp: %d", sum(rs$width_nt <= 9)),
  "",
  "--- RiboTaper periodic ORF length distribution (nt) ---",
  paste(capture.output(summary(rt_u$width_nt)), collapse = "\n"),
  ""
)

overlap_pairs <- function(a, b) {
  if (nrow(a) == 0 || nrow(b) == 0) return(0L)
  if (!all(c("seqnames", "start", "end") %in% names(a))) return(0L)
  gr_a <- GRanges(a$seqnames, IRanges(a$start, a$end), strand = a$strand)
  gr_b <- GRanges(b$seqnames, IRanges(b$start, b$end))
  m <- findOverlaps(gr_a, gr_b, minoverlap = 3L)
  length(unique(queryHits(m)))
}

summary_lines <- c(summary_lines,
  sprintf("Ribostan periodic uORFs with >=3bp genomic overlap to RiboTaper: %d / %d",
          overlap_pairs(rs, rt_u), nrow(rs)),
  "",
  "=== Interpretation notes ===",
  "1. Ribostan defaults: ATG-only scan, minimumLength=3 (>=6 bp after stop strip).",
  "2. Comparison uses periodicity-filtered Ribostan uORFs, not raw sequence candidates.",
  "3. RiboTaper side uses periodic ORFs with track length > 5 nt.",
  "4. Sparse chr22 RPFs and shared ribo/RNA BAM may limit sensitivity on both sides.",
  "5. Coordinate systems differ; treat exact-match counts as strict and overlap as directional."
)

writeLines(summary_lines, "comparison_summary.txt")
cat(paste(summary_lines, collapse = "\n"), "\n")
