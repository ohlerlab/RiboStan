#!/usr/bin/env Rscript
# Fill Ribostan p-values into an existing compare_ribotaper_psite_vectors TSV.
args <- commandArgs(trailingOnly = TRUE)
workdir <- if (length(args) >= 1) args[[1]] else "."
suffix <- if (length(args) >= 2) args[[2]] else "ribotaper_matched"
cmp_tag <- if (length(args) >= 3) args[[3]] else "matched_genome"
setwd(workdir)

suppressPackageStartupMessages(library(GenomicRanges))

p <- function(base) sub("\\.", paste0(suffix, "."), paste0(base, ".tsv"))
cmp_path <- sub("\\.", paste0(suffix, "."),
                paste0("compare_ribotaper_psite_vectors_", cmp_tag, ".tsv"))
cmp <- read.delim(cmp_path, stringsAsFactors = FALSE)
rt <- read.delim(p("ribotaper_hek293_uorfs"), stringsAsFactors = FALSE)
if (!"start" %in% names(rt)) {
  rt$start <- as.integer(rt$start0) + 1L
  rt$end <- as.integer(rt$end0)
}
cand <- read.delim(p("ribostan_hek293_candidates_full"), stringsAsFactors = FALSE)
at <- read.delim(p("ribostan_hek293_all_tested"), stringsAsFactors = FALSE)
cand <- cand[, setdiff(names(cand), c("p_value", "n_psites", "pct_inframe")), drop = FALSE]
full <- merge(cand, at[, c("orf_id", "p_value", "n_psites", "pct_inframe")],
              by = "orf_id", all.x = TRUE)
full$is_survivor <- !is.na(full$p_value)

rt_gr <- GRanges(rt$seqnames, IRanges(rt$start, rt$end), strand = rt$strand)
mcols(rt_gr)$orf_id_tr <- rt$orf_id_tr
cand_gr <- GRanges(full$seqnames, IRanges(full$start, full$end), strand = full$strand)
mcols(cand_gr) <- full
hits <- findOverlaps(rt_gr, cand_gr, minoverlap = 9L)
row_idx <- split(seq_len(nrow(rt)), rt$orf_id_tr)

best_rs <- function(orf_id_tr) {
  qi <- row_idx[[orf_id_tr]]
  empty <- list(p = NA_real_, n = NA_real_, fi = NA_real_, orf_id = NA_character_)
  if (is.null(qi)) return(empty)
  ov <- unique(subjectHits(hits)[queryHits(hits) %in% qi])
  if (!length(ov)) return(empty)
  rows <- full[ov, , drop = FALSE]
  surv <- rows[rows$is_survivor, , drop = FALSE]
  if (!nrow(surv)) return(empty)
  best <- surv[which.min(surv$p_value), ]
  list(p = as.numeric(best$p_value), n = as.numeric(best$n_psites),
       fi = as.numeric(best$pct_inframe), orf_id = as.character(best$orf_id))
}

for (i in seq_len(nrow(cmp))) {
  rs <- best_rs(cmp$orf_id_tr[i])
  cmp$p_ribostan[i] <- rs$p
  cmp$n_ribostan[i] <- rs$n
  cmp$fi_ribostan[i] <- rs$fi
  cmp$rs_orf_id[i] <- rs$orf_id
}
write.table(cmp, cmp_path, sep = "\t", row.names = FALSE, quote = FALSE)

summarize <- function(label) {
  x <- cmp[cmp$fate == label, , drop = FALSE]
  cat(sprintf("\n=== %s (n=%d) ===\n", label, nrow(x)))
  with_vec <- x[!is.na(x$p_ribotaper_psites), , drop = FALSE]
  cat(sprintf("  RT vector: median p=%.4g  pass(p<0.05)=%d/%d  median n=%s\n",
              median(with_vec$p_ribotaper_psites),
              sum(with_vec$p_ribotaper_psites < 0.05), nrow(with_vec),
              median(with_vec$rt_n_psites)))
  both <- with_vec[!is.na(with_vec$p_ribostan), , drop = FALSE]
  cat(sprintf("  Ribostan: median p=%.4g  pass(p<0.05)=%d/%d\n",
              median(both$p_ribostan), sum(both$p_ribostan < 0.05), nrow(both)))
  cat(sprintf("  RT pass & RS fail: %d; RT fail & RS pass: %d; both pass: %d; both fail: %d\n",
              sum(both$p_ribotaper_psites < 0.05 & both$p_ribostan >= 0.05),
              sum(both$p_ribotaper_psites >= 0.05 & both$p_ribostan < 0.05),
              sum(both$p_ribotaper_psites < 0.05 & both$p_ribostan < 0.05),
              sum(both$p_ribotaper_psites >= 0.05 & both$p_ribostan >= 0.05)))
  cat(sprintf("  median published p_sites=%s; median rebuilt n=%s\n",
              median(x$rt_published_p_sites), median(x$rt_n_psites)))
}
summarize("C1_tested_not_significant")
summarize("RECOVERED")
message("Updated ", cmp_path)
