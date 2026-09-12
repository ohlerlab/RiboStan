#!/usr/bin/env Rscript
# Detailed comparison of translation-supported Ribostan vs RiboTaper ORFs.
args <- commandArgs(trailingOnly = TRUE)
workdir <- if (length(args) >= 1) args[[1]] else "."
setwd(workdir)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
})

rs <- read.delim("ribostan_uorfs.tsv", stringsAsFactors = FALSE)
if (!"set" %in% names(rs)) rs$set <- "periodic"

rt_path <- "ORFs_CCDS/best_periodicity/ORFs_sign_notfiltered_multi"
if (!file.exists(rt_path)) {
  stop("Missing ", rt_path, " — run RiboTaper ORF finder first.")
}
rt <- read.delim(rt_path, sep = "\t", stringsAsFactors = FALSE, quote = "")

parse_orf_gen <- function(x) {
  parts <- strsplit(x, "_", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) < 3) return(NA_character_)
    p[[1]]
  }, character(1))
}

parse_gen_start <- function(x) {
  parts <- strsplit(x, "_", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) < 3) return(NA_integer_)
    as.integer(p[[2]])
  }, integer(1))
}

parse_gen_end <- function(x) {
  parts <- strsplit(x, "_", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) < 3) return(NA_integer_)
    as.integer(p[[3]])
  }, integer(1))
}

rt$seqnames <- parse_orf_gen(rt$ORF_id_gen)
rt$start <- parse_gen_start(rt$ORF_id_gen)
rt$end <- parse_gen_end(rt$ORF_id_gen)
rt <- rt[!is.na(rt$start) & !is.na(rt$end), , drop = FALSE]

# ORF_length is track length minus 1 (nucleotides), not codons.
rt$width_nt <- as.integer(rt$ORF_length) + 1L
rt$width_codons <- (rt$width_nt - 1L) / 3L
rt$periodic <- rt$ORF_pval_multi_ribo < 0.05

# Translation-supported only on both sides.
rt_periodic <- rt[rt$periodic & rt$width_nt > 5L, , drop = FALSE]

n_candidates <- NA_integer_
cand_path <- "ribostan_uorfs_candidates.tsv"
if (file.exists(cand_path)) {
  n_candidates <- nrow(read.delim(cand_path, stringsAsFactors = FALSE))
}

rs_gr <- GRanges(
  rs$seqnames, IRanges(rs$start, rs$end), strand = rs$strand,
  orf_id = rs$orf_id
)
rt_gr <- GRanges(
  rt_periodic$seqnames, IRanges(rt_periodic$start, rt_periodic$end),
  strand = rt_periodic$strand, orf_id = rt_periodic$ORF_id_tr
)

m <- findOverlaps(rs_gr, rt_gr, minoverlap = 9L)
rs_hit <- unique(queryHits(m))
rt_hit <- unique(subjectHits(m))

reciprocal <- function(min_frac = 0.5) {
  hits <- findOverlaps(rs_gr, rt_gr, minoverlap = 3L)
  if (length(hits) == 0) return(integer(0))
  q <- queryHits(hits)
  s <- subjectHits(hits)
  ok <- mapply(function(qi, si) {
    ov <- pmin(end(rs_gr[qi]), end(rt_gr[si])) - pmax(start(rs_gr[qi]), start(rt_gr[si])) + 1L
    wq <- width(rs_gr[qi])
    ws <- width(rt_gr[si])
    (ov / wq >= min_frac) | (ov / ws >= min_frac)
  }, q, s)
  unique(q[ok])
}

rs_recip <- reciprocal(0.5)

lines <- c(
  "=== Ribostan vs RiboTaper (translation-supported, genomic overlap, chr22) ===",
  "",
  sprintf("Ribostan sequence candidates (not compared): %s",
          if (is.na(n_candidates)) "n/a" else as.character(n_candidates)),
  sprintf("Ribostan periodic uORFs (q < 0.05): %d", nrow(rs)),
  sprintf("  3bp (1 codon after stop strip): %d", sum(rs$width_nt == 3)),
  sprintf("  6bp (2 codons): %d", sum(rs$width_nt == 6)),
  sprintf("  <= 9bp total: %d", sum(rs$width_nt <= 9)),
  sprintf("  overlaps annotated CDS: %d", sum(rs$overlaps_cds)),
  "",
  sprintf("RiboTaper periodic ORFs (p<0.05, width>5 nt): %d", nrow(rt_periodic)),
  sprintf("  <= 9 nt: %d", sum(rt_periodic$width_nt <= 9)),
  sprintf("  median length (nt): %.0f", median(rt_periodic$width_nt)),
  "",
  "--- Overlap (genomic, same strand) ---",
  sprintf("Ribostan periodic uORFs with >=9bp overlap: %d / %d (%.1f%%)",
          length(rs_hit), nrow(rs), if (nrow(rs) > 0) 100 * length(rs_hit) / nrow(rs) else 0),
  sprintf("Ribostan periodic uORFs with >=50%% reciprocal overlap: %d / %d (%.1f%%)",
          length(rs_recip), nrow(rs), if (nrow(rs) > 0) 100 * length(rs_recip) / nrow(rs) else 0),
  sprintf("RiboTaper periodic ORFs hit: %d / %d",
          length(rt_hit), nrow(rt_periodic)),
  "",
  "--- RiboTaper filtered uORFs (sORFs_sign_filtered_cds) ---",
  sprintf("  count: %d",
          max(0L, nrow(read.delim("ORFs_CCDS/best_periodicity/sORFs_sign_filtered_cds",
                                  sep = "\t", stringsAsFactors = FALSE)) - 1L)),
  "",
  "=== Interpretation ===",
  "• Both tools are compared on translation-supported ORFs only.",
  "• Ribostan sequence candidates (minimumLength=3, >=6 bp post-strip) are exported",
  "  separately and are not used in overlap statistics here.",
  "• RiboTaper requires triplet periodicity in RPF data (multitaper p<0.05) and",
  "  minimum track length > 5 nt.",
  "• Sparse chr22 RPFs limit sensitivity; low overlap may reflect data, not method bugs."
)

if (length(rs_recip) > 0) {
  lines <- c(lines, "", "--- Example reciprocal overlaps (first 5) ---")
  for (i in head(rs_recip, 5)) {
    h <- findOverlaps(rs_gr[i], rt_gr, minoverlap = 3L)
    if (length(h) == 0) next
    j <- subjectHits(h)[1]
    lines <- c(lines, sprintf(
      "  RS %s (%s:%d-%d %dbp) <-> RT %s (%s:%d-%d %d nt)",
      rs$orf_id[i], rs$seqnames[i], rs$start[i], rs$end[i], rs$width_nt[i],
      rt_periodic$ORF_id_tr[j], rt_periodic$seqnames[j],
      rt_periodic$start[j], rt_periodic$end[j], rt_periodic$width_nt[j]
    ))
  }
}

writeLines(lines, "comparison_detailed.txt")
cat(paste(lines, collapse = "\n"), "\n")

if (length(m) > 0) {
  ov_df <- data.frame(
    ribostan_orf = rs$orf_id[queryHits(m)],
    ribostan_start = rs$start[queryHits(m)],
    ribostan_end = rs$end[queryHits(m)],
    ribostan_width = rs$width_nt[queryHits(m)],
    ribotaper_orf = rt_periodic$ORF_id_tr[subjectHits(m)],
    ribotaper_start = rt_periodic$start[subjectHits(m)],
    ribotaper_end = rt_periodic$end[subjectHits(m)],
    ribotaper_width_nt = rt_periodic$width_nt[subjectHits(m)],
    stringsAsFactors = FALSE
  )
  write.table(ov_df, "comparison_overlaps.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
  message("Wrote comparison_overlaps.tsv (", nrow(ov_df), " pairs)")
}
