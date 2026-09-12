#!/usr/bin/env Rscript
# Per-uORF diagnostic: for each RiboTaper HEK293 uORF, why does/doesn't
# Ribostan recover it; for each Ribostan periodic non-CDS-overlap call, what
# is it if it isn't in RiboTaper's filtered uORF set. Uses only the already
# written outputs of a single compare_hek293_prep.R / compare_hek293_report.R
# run (no rerun of the expensive periodicity test).
args <- commandArgs(trailingOnly = TRUE)
workdir <- if (length(args) >= 1) args[[1]] else "."
suffix  <- if (length(args) >= 2) args[[2]] else Sys.getenv("HEK293_OUT_SUFFIX", "ribotaper_matched")
setwd(workdir)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
})

# Must match compare_hek293_prep.R / compare_hek293_report.R out_path().
p <- function(base) sub("\\.", paste0(suffix, "."), paste0(base, ".tsv"))

cand     <- read.delim(p("ribostan_hek293_candidates_full"), stringsAsFactors = FALSE)
at       <- read.delim(p("ribostan_hek293_all_tested"), stringsAsFactors = FALSE)
periodic <- read.delim(p("ribostan_hek293_periodic"), stringsAsFactors = FALSE)
rt       <- read.delim(p("ribotaper_hek293_uorfs"), stringsAsFactors = FALSE)
matched  <- read.delim(sub("\\.", paste0(suffix, "."),
                          "hek293_comparison_overlaps_nonCDSoverlap.tsv"),
                        stringsAsFactors = FALSE)

stopifnot(all(cand$orf_id == at$orf_id))
full <- cbind(
  cand[, c("orf_id", "transcript_id", "seqnames", "start", "end", "strand",
           "width_nt", "overlaps_cds")],
  at[, c("p_value", "q_value", "n_psites", "pct_inframe")]
)
# survivor of the nested-uORF selection <=> got a p_value written back
full$is_survivor <- !is.na(full$p_value)
full$stop_key <- paste0(full$transcript_id, ":", ifelse(full$strand == "-", full$start, full$end))

gr_from_df <- function(df, id_col = "orf_id", extra = NULL) {
  gr <- GRanges(seqnames = df$seqnames, ranges = IRanges(df$start, df$end),
                strand = df$strand)
  mcols(gr)$orf_id <- df[[id_col]]
  for (e in extra) mcols(gr)[[e]] <- df[[e]]
  gr
}

cand_gr <- gr_from_df(full, extra = c("p_value", "n_psites", "pct_inframe",
                                       "is_survivor", "stop_key", "transcript_id"))
rt_gr   <- gr_from_df(rt, "orf_id_tr", extra = "transcript_id")

## ---------------------------------------------------------------
## Analysis 1: why does Ribostan miss each RiboTaper filtered uORF?
## RiboTaper's BED splits multi-exon uORFs into one row per exon segment
## (908 rows / 763 unique orf_id_tr / 701 transcripts), so classification is
## done per unique orf_id_tr, pooling overlap evidence across all of its rows.
## ---------------------------------------------------------------
hits <- findOverlaps(rt_gr, cand_gr, minoverlap = 9L)
recovered_rt <- unique(matched$ribotaper_orf)

rt_ids <- unique(rt$orf_id_tr)
row_idx_by_id <- split(seq_len(nrow(rt)), rt$orf_id_tr)

classify_miss <- function(id) {
  row_i <- row_idx_by_id[[id]]
  ov <- unique(subjectHits(hits)[queryHits(hits) %in% row_i])
  tx <- rt$transcript_id[row_i[1]]
  if (length(ov) == 0L) {
    if (!(tx %in% full$transcript_id)) {
      return("A1_transcript_absent")
    }
    return("A2_transcript_present_no_candidate_at_locus")
  }
  rows <- full[ov, , drop = FALSE]
  if (!any(rows$is_survivor)) {
    if (all(is.na(rows$n_psites) | rows$n_psites < 2)) {
      return("B1_candidate_exists_too_few_reads_to_test")
    }
    return("B2_candidate_exists_lost_nested_selection_to_sibling")
  }
  best <- rows[rows$is_survivor, , drop = FALSE]
  best <- best[order(best$p_value), , drop = FALSE][1, ]
  if (best$p_value >= 0.05) return("C1_tested_not_significant")
  if (is.na(best$n_psites) || best$n_psites < 6) return("C2_significant_but_lt6_psites")
  if (is.na(best$pct_inframe) || best$pct_inframe <= 0.5) return("C3_significant_but_failed_inframe_gate")
  "D_passed_all_filters"
}

rt_id_class <- vapply(rt_ids, function(id) {
  if (id %in% recovered_rt) "RECOVERED" else classify_miss(id)
}, character(1))
rt_id_df <- data.frame(orf_id_tr = rt_ids, class = rt_id_class, stringsAsFactors = FALSE)
rt_id_df$transcript_id <- rt$transcript_id[match(rt_id_df$orf_id_tr, rt$orf_id_tr)]

cat("=== Analysis 1: fate of each unique RiboTaper filtered uORF (n=", length(rt_ids),
    " unique uORFs; ", nrow(rt), " BED rows after multi-exon splitting) ===\n", sep = "")
print(table(rt_id_df$class))
cat("\n(D_passed_all_filters should be ~0; any >0 rows are boundary/dedup edge cases)\n\n")

for (cl in setdiff(unique(rt_id_df$class), "RECOVERED")) {
  idx <- which(rt_id_df$class == cl)
  ex <- head(rt_id_df$orf_id_tr[idx], 6)
  cat(sprintf("-- %s (n=%d) examples: %s\n", cl, length(idx), paste(ex, collapse = ", ")))
}

write.table(rt_id_df, sub("\\.", paste0(suffix, "."), "diag_ribotaper_fate.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

## ---------------------------------------------------------------
## Analysis 2: what are the Ribostan periodic (non-CDS-overlap) calls
## that don't match any RiboTaper filtered uORF?
## ---------------------------------------------------------------
periodic_noncds <- periodic[!is.na(periodic$overlaps_cds) & !periodic$overlaps_cds, , drop = FALSE]
matched_rs_ids <- unique(matched$ribostan_orf)
excess <- periodic_noncds[!(periodic_noncds$orf_id %in% matched_rs_ids), , drop = FALSE]

all_bed <- file.path(workdir, "GSM1887643_HEK293_all_ORFs.bed.gz")
con <- gzfile(all_bed, open = "rt")
raw <- read.delim(con, header = FALSE, stringsAsFactors = FALSE, comment.char = "", quote = "")
close(con)
names(raw) <- c("seqnames", "start0", "end0", "name", "p_sites", "strand")[seq_len(ncol(raw))]
meta <- strsplit(raw$name, ";", fixed = TRUE)
raw$orf_id_tr <- vapply(meta, `[`, character(1), 1)
raw$category  <- vapply(meta, function(x) if (length(x) >= 2) x[[2]] else NA_character_, character(1))
raw$start <- raw$start0 + 1L
raw$end   <- raw$end0
allrt_gr <- GRanges(seqnames = raw$seqnames, ranges = IRanges(raw$start, raw$end),
                    strand = raw$strand)
mcols(allrt_gr)$category <- raw$category
mcols(allrt_gr)$orf_id_tr <- raw$orf_id_tr

excess_gr <- gr_from_df(excess)
hits2 <- findOverlaps(excess_gr, allrt_gr, minoverlap = 9L)
excess$class <- "E4_no_ribotaper_overlap_at_all"
if (length(hits2) > 0) {
  cat_per_query <- split(mcols(allrt_gr)$category[subjectHits(hits2)], queryHits(hits2))
  best_cat <- vapply(cat_per_query, function(cs) {
    # all_ORFs.bed has 922 "uORF" rows vs 908 in the final filtered BED -
    # these 14 extra rows are RiboTaper's own raw uORF candidates that its
    # ORF-selection/QC step rejected. A hit here means Ribostan flagged a
    # locus RiboTaper *also* considered as a uORF candidate but excluded
    # from its own final list, not a locus RiboTaper never looked at.
    if ("uORF" %in% cs) return("E0_overlaps_ribotaper_rejected_uorf_candidate")
    if ("ORFs_ccds" %in% cs) return("E1_overlaps_ORFs_ccds")
    if (any(cs %in% c("dORF", "Overl_dORF"))) return("E2_overlaps_dORF")
    if (any(cs %in% c("ncORFS", "nonccds_coding_ORFs"))) return("E3_overlaps_other_noncoding_ORF")
    "E5_overlaps_uncategorized"
  }, character(1))
  qidx <- as.integer(names(cat_per_query))
  excess$class[qidx] <- best_cat
}

cat("\n=== Analysis 2: fate of Ribostan periodic non-CDS-overlap calls NOT matching a RiboTaper uORF (n=",
    nrow(excess), " of ", nrow(periodic_noncds), " total) ===\n", sep = "")
print(table(excess$class))
cat("\n")

# For the fully novel bucket, check call strength (are these weak/borderline
# or robust calls) and whether they sit on a transcript RiboTaper called
# translated at all (any category) vs a transcript RiboTaper never touched.
novel <- excess[excess$class == "E4_no_ribotaper_overlap_at_all", , drop = FALSE]
novel$rt_touches_transcript <- novel$transcript_id %in% raw$transcript_id
cat(sprintf(
  "Of the %d fully novel calls: %d (%.0f%%) are on a transcript RiboTaper found\n",
  nrow(novel), sum(novel$rt_touches_transcript),
  100 * mean(novel$rt_touches_transcript)))
cat(sprintf("  ANY signal on (different position/category); %d (%.0f%%) are on transcripts RiboTaper never called at all.\n",
            sum(!novel$rt_touches_transcript), 100 * mean(!novel$rt_touches_transcript)))

at_novel <- at[match(novel$orf_id, at$orf_id), ]
cat(sprintf(
  "Novel-call strength: median n_psites = %s, median p_value = %.4g, median pct_inframe = %.2f\n",
  median(at_novel$n_psites, na.rm = TRUE),
  median(at_novel$p_value, na.rm = TRUE),
  median(at_novel$pct_inframe, na.rm = TRUE)
))
cat(sprintf(
  "  (for comparison, matched/recovered calls: median n_psites = %s, median p_value = %.4g)\n",
  median(at$n_psites[match(matched_rs_ids, at$orf_id)], na.rm = TRUE),
  median(at$p_value[match(matched_rs_ids, at$orf_id)], na.rm = TRUE)
))

write.table(excess, sub("\\.", paste0(suffix, "."), "diag_ribostan_excess.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nWrote diag_ribotaper_fate_", suffix, ".tsv and diag_ribostan_excess_", suffix, ".tsv\n", sep = "")
