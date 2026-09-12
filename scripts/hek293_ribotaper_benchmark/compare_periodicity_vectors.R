#!/usr/bin/env Rscript
# Rebuild per-uORF P-site coverage vectors on matched inputs under several
# P-site constructions, run the same multitaper test (ftestvect) on each,
# and compare p-values to the stored pipeline result + RiboTaper fate classes.
args <- commandArgs(trailingOnly = TRUE)
workdir <- if (length(args) >= 1) args[[1]] else "."
suffix <- if (length(args) >= 2) args[[2]] else Sys.getenv("HEK293_OUT_SUFFIX", "ribotaper_matched")
setwd(workdir)

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args_all[grep("--file=", args_all)])
bench_dir <- if (length(script_path)) dirname(script_path) else "."
source(file.path(bench_dir, "hek293_gtf_config.R"), local = TRUE)

suppressPackageStartupMessages({
  library(Ribostan)
  library(GenomicRanges)
  library(IRanges)
  library(GenomicFeatures)
  library(dplyr)
})

out_path <- function(name) sub("\\.", paste0(suffix, "."), name)
pbase <- function(base) sub("\\.", paste0(suffix, "."), paste0(base, ".tsv"))

message("Loading matched inputs ...")
gtf_cfg <- resolve_hek293_gtf(workdir)
anno <- load_annotation(gtf_cfg$path, "GRCh37.p13.genome.fa", add_uorfs = TRUE)
rpfs <- readRDS(Sys.getenv("HEK293_PROJECTED_RDS", "rpfs_genome_projected_rrnafilt_matched.rds"))
offsets_df <- read.table(out_path("ribostan_offsets.tsv"), header = TRUE, stringsAsFactors = FALSE)
at <- read.delim(pbase("ribostan_hek293_all_tested"), stringsAsFactors = FALSE)
cand <- read.delim(pbase("ribostan_hek293_candidates_full"), stringsAsFactors = FALSE)
rt <- read.delim(pbase("ribotaper_hek293_uorfs"), stringsAsFactors = FALSE)
fate <- read.delim(sub("\\.", paste0(suffix, "."), "diag_ribotaper_fate.tsv"),
                   stringsAsFactors = FALSE)

# --- P-site builders ---------------------------------------------------------
message("Building P-site sets (matched prep + with phase-shiftdf) ...")
apply_phase_shiftdf <- as.logical(Sys.getenv("HEK293_APPLY_PHASE_SHIFT", "FALSE"))
ps_pipeline <- get_psite_gr(rpfs, offsets_df, anno, apply_phase_shiftdf = apply_phase_shiftdf)
ps_with_shiftdf <- get_psite_gr(rpfs, offsets_df, anno, apply_phase_shiftdf = TRUE)

orf_coverage <- function(psites, orf_id, anno) {
  orf <- anno$trspacecds[orf_id]
  if (length(orf) == 0L) return(NULL)
  sub <- psites[as.character(psites$orf) == orf_id]
  if (length(sub) == 0L) return(NULL)
  mapped <- mapToTranscripts(sub, orf, ignore.strand = TRUE)
  mapped <- mapped[sub$orf[mapped$xHits] == orf_id]
  if (length(mapped) == 0L) return(NULL)
  cov <- coverage(mapped)
  vec <- as.numeric(unlist(cov, use.names = FALSE))
  if (length(vec) == 0L || sum(vec > 0) <= 1L) return(NULL)
  vec
}

run_ftest <- function(vec) {
  if (is.null(vec)) return(c(Fstat = NA, p.value = NA, spec_coef = NA))
  Ribostan:::ftestvect(vec)
}

pct_inframe <- function(vec) {
  if (is.null(vec)) return(NA_real_)
  total <- sum(vec)
  if (total == 0) return(NA_real_)
  sum(vec[seq(1, length(vec), by = 3)]) / total
}

# --- Map RiboTaper uORFs -> best overlapping Ribostan candidate ---------------
full <- cbind(
  cand[, c("orf_id", "transcript_id", "seqnames", "start", "end", "strand")],
  at[, c("p_value", "n_psites", "pct_inframe")]
)
full$is_survivor <- !is.na(full$p_value)
cand_gr <- GRanges(full$seqnames, IRanges(full$start, full$end), strand = full$strand)
mcols(cand_gr) <- full
rt_gr <- GRanges(rt$seqnames, IRanges(rt$start, rt$end), strand = rt$strand)
mcols(rt_gr)$orf_id_tr <- rt$orf_id_tr
hits <- findOverlaps(rt_gr, cand_gr, minoverlap = 9L)

rt_ids <- unique(rt$orf_id_tr)
row_idx <- split(seq_len(nrow(rt)), rt$orf_id_tr)

best_candidate <- function(rt_id) {
  qi <- row_idx[[rt_id]]
  ov <- unique(subjectHits(hits)[queryHits(hits) %in% qi])
  if (length(ov) == 0L) return(NA_character_)
  rows <- full[ov, , drop = FALSE]
  surv <- rows[rows$is_survivor, , drop = FALSE]
  if (nrow(surv) == 0L) return(NA_character_)
  surv$orf_id[which.min(surv$p_value)]
}

message("Evaluating ", length(rt_ids), " unique RiboTaper uORFs ...")
rows <- lapply(rt_ids, function(rt_id) {
  cand_id <- best_candidate(rt_id)
  cls <- fate$class[match(rt_id, fate$orf_id_tr)]
  if (is.na(cand_id)) {
    return(data.frame(
      orf_id_tr = rt_id, rs_orf_id = NA_character_, fate = cls,
      rt_p_sites = max(rt$p_sites[rt$orf_id_tr == rt_id], na.rm = TRUE),
      stored_p = NA, stored_n = NA, stored_fi = NA,
      p_pipeline = NA, p_with_shiftdf = NA,
      n_pipeline = NA, fi_pipeline = NA,
      stringsAsFactors = FALSE
    ))
  }
  stored <- full[full$orf_id == cand_id, ]
  v_pipe <- orf_coverage(ps_pipeline, cand_id, anno)
  v_shift <- orf_coverage(ps_with_shiftdf, cand_id, anno)
  t_pipe <- run_ftest(v_pipe)
  t_shift <- run_ftest(v_shift)
  data.frame(
    orf_id_tr = rt_id,
    rs_orf_id = cand_id,
    fate = cls,
    rt_p_sites = max(rt$p_sites[rt$orf_id_tr == rt_id], na.rm = TRUE),
    stored_p = stored$p_value[1],
    stored_n = stored$n_psites[1],
    stored_fi = stored$pct_inframe[1],
    p_pipeline = t_pipe["p.value"],
    p_with_shiftdf = t_shift["p.value"],
    n_pipeline = sum(v_pipe, na.rm = TRUE),
    fi_pipeline = pct_inframe(v_pipe),
    stringsAsFactors = FALSE
  )
})
res <- bind_rows(rows)

# Sanity: recomputed pipeline p should match stored p for survivors
both <- res[!is.na(res$stored_p) & !is.na(res$p_pipeline), ]
if (nrow(both) > 0) {
  both$log10_p_diff <- abs(log10(both$p_pipeline) - log10(both$stored_p))
  message(sprintf(
    "Sanity stored vs recomputed pipeline p: median |log10 diff| = %.3g (max %.3g, n=%d)",
    median(both$log10_p_diff, na.rm = TRUE),
    max(both$log10_p_diff, na.rm = TRUE),
    nrow(both)
  ))
}

summarize_class <- function(df, label, pcol) {
  x <- df[!is.na(df[[pcol]]), ]
  cat("\n=== ", label, " (pcol=", pcol, ", n=", nrow(x), ") ===\n", sep = "")
  for (cl in c("RECOVERED", "C1_tested_not_significant", "A1_transcript_absent")) {
    m <- x$fate == cl
    if (!any(m)) next
    cat(sprintf("  %s: median p=%.3g  pass(p<0.05)=%d/%d  median fi=%.2f\n",
                cl, median(x[[pcol]][m], na.rm = TRUE),
                sum(x[[pcol]][m] < 0.05, na.rm = TRUE), sum(m),
                median(x$fi_pipeline[m], na.rm = TRUE)))
  }
}

cat("\n=== Per-vector multitaper on matched genome inputs ===\n")
cat("Modes: pipeline = prep get_psite_gr (apply_phase_shiftdf=",
    apply_phase_shiftdf, "); with_shiftdf = Ribostan default\n", sep = "")
summarize_class(res, "Matched prep psites", "p_pipeline")
summarize_class(res, "With CDS phase-shiftdf", "p_with_shiftdf")

c1 <- res[res$fate == "C1_tested_not_significant" & !is.na(res$p_pipeline), ]
if (nrow(c1) > 0) {
  c1$flipped_by_noshift <- c1$p_with_shiftdf >= 0.05 & c1$p_pipeline < 0.05
  cat(sprintf(
    "\nC1 loci where dropping phase-shiftdf flips p<0.05: %d / %d (%.1f%%)\n",
    sum(c1$flipped_by_noshift, na.rm = TRUE), nrow(c1),
    100 * mean(c1$flipped_by_noshift, na.rm = TRUE)
  ))
  c1$flipped_worse <- c1$p_pipeline < 0.05 & c1$p_with_shiftdf >= 0.05
  cat(sprintf(
    "C1 loci where phase-shiftdf HELPS (with_shiftdf pass, prep fail): %d\n",
    sum(c1$flipped_worse, na.rm = TRUE)
  ))
}

out_tsv <- sub("\\.", paste0(suffix, "."), "compare_periodicity_vectors.tsv")
write.table(res, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
message("Wrote ", out_tsv)
