#!/usr/bin/env Rscript
# Compare multitaper p-values on RiboTaper-exact P-site vectors vs Ribostan.
# Default: transcript-space P-sites (how RiboTaper tests uORFs). Set
# HEK293_RT_PSITE_MODE=genome for genomic bedtools coverage (legacy).
args <- commandArgs(trailingOnly = TRUE)
workdir <- if (length(args) >= 1) args[[1]] else "."
suffix <- if (length(args) >= 2) args[[2]] else Sys.getenv("HEK293_OUT_SUFFIX", "ribotaper_matched")
setwd(workdir)

suppressPackageStartupMessages({
  library(Ribostan)
  library(GenomicRanges)
  library(IRanges)
})

p <- function(base) sub("\\.", paste0(suffix, "."), paste0(base, ".tsv"))

ensure_genomic_coords <- function(df) {
  if (!"start" %in% names(df) && "start0" %in% names(df)) {
    df$start <- as.integer(df$start0) + 1L
    df$end <- as.integer(df$end0)
  }
  df
}

fate <- read.delim(sub("\\.", paste0(suffix, "."), "diag_ribotaper_fate.tsv"),
                   stringsAsFactors = FALSE)
rt <- ensure_genomic_coords(read.delim(p("ribotaper_hek293_uorfs"), stringsAsFactors = FALSE))
cand <- read.delim(p("ribostan_hek293_candidates_full"), stringsAsFactors = FALSE)
at <- read.delim(p("ribostan_hek293_all_tested"), stringsAsFactors = FALSE)
# candidates_full may carry placeholder NA stats; prefer all_tested values.
stat_cols <- c("p_value", "n_psites", "pct_inframe")
cand <- cand[, setdiff(names(cand), stat_cols), drop = FALSE]
full <- merge(cand, at[, c("orf_id", stat_cols)], by = "orf_id", all.x = TRUE)
stopifnot("p_value" %in% names(full), mean(is.na(full$p_value)) < 1)
full$is_survivor <- !is.na(full$p_value)

c1_ids <- fate$orf_id_tr[fate$class == "C1_tested_not_significant"]
rec_ids <- fate$orf_id_tr[fate$class == "RECOVERED"]
eval_ids <- unique(c(fate$orf_id_tr[fate$class == "C1_tested_not_significant"],
                     head(rec_ids, 100L)))

mode <- tolower(Sys.getenv("HEK293_RT_PSITE_MODE", "transcript"))
psites_bed <- Sys.getenv(
  "HEK293_RIBOTAPER_PSITES",
  if (mode == "transcript") {
    file.path("ribotaper_psites_transcript_published_bam", "P_sites_all.bed6")
  } else {
    file.path("ribotaper_psites_published_bam", "P_sites_all.clean.bed6")
  }
)
stopifnot(file.exists(psites_bed))
message("Mode=", mode, "; evaluating ", length(c1_ids), " C1 + ",
        min(100L, length(rec_ids)), " RECOVERED using ", psites_bed)

run_ftest <- function(vec) {
  if (is.null(vec) || length(vec) < 2L || sum(vec) <= 1L) {
    return(c(Fstat = NA, p.value = NA, spec_coef = NA))
  }
  Ribostan:::ftestvect(as.numeric(vec))
}

# --- Ribostan overlap lookup (once) -------------------------------------------
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
  if (length(ov) == 0L) return(empty)
  rows <- full[ov, , drop = FALSE]
  surv <- rows[rows$is_survivor, , drop = FALSE]
  if (nrow(surv) == 0L) return(empty)
  best <- surv[which.min(surv$p_value), ]
  list(
    p = as.numeric(best$p_value),
    n = as.numeric(best$n_psites),
    fi = as.numeric(best$pct_inframe),
    orf_id = as.character(best$orf_id)
  )
}

fate_of <- function(id) {
  if (id %in% c1_ids) "C1_tested_not_significant" else "RECOVERED"
}

# --- Transcript-space vectors (RiboTaper uORF testing) ------------------------
load_tx_psite_index <- function(bed_path) {
  message("Loading transcript P-sites from ", bed_path, " ...")
  ps <- read.delim(bed_path, header = FALSE, stringsAsFactors = FALSE,
                   col.names = c("transcript_id", "start0", "end0", "name", "score", "strand"))
  ps <- ps[nzchar(ps$transcript_id), , drop = FALSE]
  ps$pos <- as.integer(ps$start0) + 1L
  # Aggregate duplicate P-site calls at the same transcript position.
  agg <- aggregate(rep(1L, nrow(ps)), by = list(tx = ps$transcript_id, pos = ps$pos), FUN = sum)
  names(agg)[3L] <- "depth"
  agg
}

build_tx_vec <- function(orf_id_tr, rt_df, ps_agg) {
  rows <- rt_df[rt_df$orf_id_tr == orf_id_tr, , drop = FALSE]
  if (nrow(rows) == 0L || any(is.na(rows$tx_start)) || any(is.na(rows$tx_end))) {
    return(NULL)
  }
  tx <- rows$transcript_id[1L]
  t0 <- min(rows$tx_start)
  t1 <- max(rows$tx_end)
  width <- t1 - t0 + 1L
  if (width < 2L) return(NULL)
  vec <- integer(width)
  sub <- ps_agg[ps_agg$tx == tx & ps_agg$pos >= t0 & ps_agg$pos <= t1, , drop = FALSE]
  if (nrow(sub) == 0L) return(vec)
  idx <- sub$pos - t0 + 1L
  vec[idx] <- vec[idx] + as.integer(sub$depth)
  vec
}

# --- Genomic vectors via GenomicRanges (no bedtools) --------------------------
build_genome_vec_fn <- function(eval_ids, rt_df, psites_bed_path) {
  message("Loading genomic P-sites from ", psites_bed_path, " ...")
  # Stream only chr* rows with 6 fields to avoid malformed bedtools leftovers.
  tmp <- tempfile(fileext = ".bed")
  status <- system(sprintf(
    "awk 'BEGIN{FS=OFS=\"\\t\"} NF==6 && $1 ~ /^chr/ {print $1,$2,$3,\".\",0,$6}' %s > %s",
    shQuote(psites_bed_path), shQuote(tmp)
  ))
  if (status != 0L) stop("awk filter of P-sites failed")
  ps <- read.delim(tmp, header = FALSE, stringsAsFactors = FALSE,
                   col.names = c("seqnames", "start0", "end0", "name", "score", "strand"))
  unlink(tmp)
  ps_gr <- GRanges(
    ps$seqnames,
    IRanges(as.integer(ps$start0) + 1L, as.integer(ps$end0)),
    strand = ps$strand
  )
  message("Loaded ", length(ps_gr), " genomic P-sites")

  rt_eval <- rt_df[rt_df$orf_id_tr %in% eval_ids, , drop = FALSE]
  seg_i <- ave(seq_len(nrow(rt_eval)), rt_eval$orf_id_tr, FUN = seq_along)
  rt_eval$seg <- seg_i
  seg_gr <- GRanges(
    rt_eval$seqnames,
    IRanges(rt_eval$start, rt_eval$end),
    strand = rt_eval$strand
  )
  mcols(seg_gr)$orf_id_tr <- rt_eval$orf_id_tr
  mcols(seg_gr)$seg <- rt_eval$seg

  ov <- findOverlaps(ps_gr, seg_gr, ignore.strand = FALSE)
  message("P-site overlaps with eval segments: ", length(ov))

  qh <- queryHits(ov)
  sh <- subjectHits(ov)
  hit_orf <- mcols(seg_gr)$orf_id_tr[sh]
  hit_seg <- mcols(seg_gr)$seg[sh]
  hit_strand <- as.character(strand(seg_gr)[sh])
  seg_start <- start(seg_gr)[sh]
  seg_end <- end(seg_gr)[sh]
  ps_pos <- start(ps_gr)[qh]
  pos_5to3 <- ifelse(hit_strand == "+", ps_pos - seg_start + 1L, seg_end - ps_pos + 1L)

  hits_df <- data.frame(
    orf_id_tr = hit_orf,
    seg = hit_seg,
    pos = pos_5to3,
    stringsAsFactors = FALSE
  )

  seg_meta <- unique(rt_eval[, c("orf_id_tr", "seg", "start", "end", "strand")])
  seg_meta$width <- seg_meta$end - seg_meta$start + 1L

  function(orf_id) {
    meta <- seg_meta[seg_meta$orf_id_tr == orf_id, , drop = FALSE]
    if (nrow(meta) == 0L) return(NULL)
    st <- meta$strand[1]
    segs <- sort(meta$seg)
    if (st == "-") segs <- rev(segs)
    chunks <- lapply(segs, function(s) {
      w <- meta$width[meta$seg == s][1]
      vec <- integer(w)
      hh <- hits_df[hits_df$orf_id_tr == orf_id & hits_df$seg == s, , drop = FALSE]
      if (nrow(hh)) {
        ok <- hh$pos >= 1L & hh$pos <= w
        if (any(ok)) {
          vec <- as.integer(tabulate(hh$pos[ok], nbins = w))
        }
      }
      vec
    })
    unlist(chunks, use.names = FALSE)
  }
}

if (mode == "transcript") {
  ps_agg <- load_tx_psite_index(psites_bed)
  build_vec <- function(orf_id) build_tx_vec(orf_id, rt, ps_agg)
} else {
  build_vec <- build_genome_vec_fn(eval_ids, rt, psites_bed)
}

message("Computing multitaper per ORF ...")
rows <- lapply(eval_ids, function(id) {
  v <- build_vec(id)
  t <- run_ftest(v)
  rs <- best_rs(id)
  data.frame(
    orf_id_tr = id,
    fate = fate_of(id),
    rt_n_psites = if (is.null(v)) 0L else sum(v),
    rt_published_p_sites = max(rt$p_sites[rt$orf_id_tr == id], na.rm = TRUE),
    p_ribotaper_psites = unname(t["p.value"]),
    p_ribostan = rs$p,
    rs_orf_id = rs$orf_id,
    n_ribostan = rs$n,
    fi_ribostan = rs$fi,
    stringsAsFactors = FALSE
  )
})
res <- do.call(rbind, rows)

out_tag <- Sys.getenv("HEK293_COMPARE_OUT_TAG", "")
out_base <- if (nzchar(out_tag)) {
  paste0("compare_ribotaper_psite_vectors_", out_tag, ".tsv")
} else {
  "compare_ribotaper_psite_vectors.tsv"
}
out <- sub("\\.", paste0(suffix, "."), out_base)
write.table(res, out, sep = "\t", row.names = FALSE, quote = FALSE)

summarize <- function(label) {
  x <- res[res$fate == label, , drop = FALSE]
  cat(sprintf("\n=== %s (n=%d) ===\n", label, nrow(x)))
  if (nrow(x) == 0L) return(invisible(NULL))
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
}
summarize("C1_tested_not_significant")
summarize("RECOVERED")
message("Wrote ", out)
