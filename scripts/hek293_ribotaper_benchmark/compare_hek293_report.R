#!/usr/bin/env Rscript
# Compare Ribostan periodic uORFs vs RiboTaper HEK293 published uORFs (GSE73136).
args <- commandArgs(trailingOnly = TRUE)
workdir <- if (length(args) >= 1) args[[1]] else "."
setwd(workdir)

out_suffix <- Sys.getenv("HEK293_OUT_SUFFIX", "")
out_path <- function(name) sub("\\.", paste0(out_suffix, "."), name)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
})

rs_all <- read.delim(out_path("ribostan_hek293_periodic.tsv"), stringsAsFactors = FALSE)
rt <- read.delim(out_path("ribotaper_hek293_uorfs.tsv"), stringsAsFactors = FALSE)

settings_path <- out_path("ribostan_hek293_settings.txt")
run_settings <- if (file.exists(settings_path)) {
  s <- readLines(settings_path)
  setNames(sub("^[^=]+=", "", s), sub("=.*$", "", s))
} else c(fdr_correct = "TRUE", nested_rule = "density", multimap = "all",
         min_frame_frac = "0", offset_method = "cds_incl", multimap_qc = "FALSE",
         gtf_name = Sys.getenv("HEK293_GTF_NAME", "gencode.v19.ccds_appris.gtf"))
if (is.na(run_settings["min_frame_frac"])) run_settings["min_frame_frac"] <- "0"
if (is.na(run_settings["offset_method"])) run_settings["offset_method"] <- "cds_incl"
if (is.na(run_settings["multimap_qc"])) run_settings["multimap_qc"] <- "FALSE"
if (is.na(run_settings["read_source"])) run_settings["read_source"] <- "transcriptome_bam"
sig_label <- if (isTRUE(as.logical(run_settings[["fdr_correct"]]))) {
  "q < 0.05 (BH-FDR corrected)"
} else "p < 0.05 (raw, uncorrected)"

gr_from_df <- function(df, id_col = "orf_id") {
  GRanges(
    seqnames = df$seqnames,
    ranges = IRanges(start = df$start, end = df$end),
    strand = df$strand,
    orf_id = df[[id_col]]
  )
}

rt_gr <- gr_from_df(rt, "orf_id_tr")

# RiboTaper's own category rule excludes CDS-overlap from its "uORF" label:
# "uORFs and dORFs as ORFs in CCDS genes non overlapping with any CDS exon"
# (Calviello et al. 2016 Methods). Ribostan's uORF finder does not apply
# this exclusion, so a like-for-like comparison needs the non-overlapping
# subset (overlaps_cds == FALSE) as well as the full periodic set.
has_overlap_col <- "overlaps_cds" %in% names(rs_all)
rs_nonoverlap <- if (has_overlap_col) {
  rs_all[!is.na(rs_all$overlaps_cds) & !rs_all$overlaps_cds, , drop = FALSE]
} else rs_all[0, , drop = FALSE]

compare_set <- function(rs, rt, rt_gr, label) {
  if (nrow(rs) == 0L) {
    return(list(
      label = label, n = 0L, n_hit = 0L, n_recip = 0L, n_rt_hit = 0L,
      pct_hit = 0, pct_recip = 0, pct_rt_hit = 0, median_width = NA_real_,
      overlaps = data.frame()
    ))
  }
  rs_gr <- gr_from_df(rs, "orf_id")
  hits <- findOverlaps(rs_gr, rt_gr, minoverlap = 9L)
  rs_hit <- unique(queryHits(hits))
  rt_hit <- unique(subjectHits(hits))

  reciprocal <- function(min_frac = 0.5) {
    h <- findOverlaps(rs_gr, rt_gr, minoverlap = 3L)
    if (length(h) == 0) return(integer(0))
    q <- queryHits(h)
    s <- subjectHits(h)
    ok <- mapply(function(qi, si) {
      ov <- pmin(end(rs_gr[qi]), end(rt_gr[si])) -
        pmax(start(rs_gr[qi]), start(rt_gr[si])) + 1L
      wq <- width(rs_gr[qi])
      ws <- width(rt_gr[si])
      (ov / wq >= min_frac) | (ov / ws >= min_frac)
    }, q, s)
    unique(q[ok])
  }
  rs_recip <- reciprocal(0.5)

  ov_df <- if (length(hits) > 0) {
    data.frame(
      ribostan_orf = rs$orf_id[queryHits(hits)],
      ribostan_width = rs$width_nt[queryHits(hits)],
      ribotaper_orf = rt$orf_id_tr[subjectHits(hits)],
      ribotaper_width = rt$width_nt[subjectHits(hits)],
      ribotaper_p_sites = rt$p_sites[subjectHits(hits)],
      stringsAsFactors = FALSE
    )
  } else data.frame()

  list(
    label = label, n = nrow(rs), n_hit = length(rs_hit), n_recip = length(rs_recip),
    n_rt_hit = length(rt_hit),
    pct_hit = 100 * length(rs_hit) / nrow(rs),
    pct_recip = 100 * length(rs_recip) / nrow(rs),
    pct_rt_hit = 100 * length(rt_hit) / nrow(rt),
    median_width = median(rs$width_nt),
    overlaps = ov_df
  )
}

res_all <- compare_set(rs_all, rt, rt_gr, "all periodic uORFs")
res_nonoverlap <- compare_set(rs_nonoverlap, rt, rt_gr, "non-CDS-overlapping periodic uORFs only")

fmt_result <- function(res) {
  c(
    sprintf("[%s] n = %d, median width (nt) = %s", res$label, res$n,
            if (is.na(res$median_width)) "n/a" else sprintf("%.0f", res$median_width)),
    sprintf("  >=9bp overlap to a RiboTaper uORF: %d / %d (%.1f%%)",
            res$n_hit, res$n, res$pct_hit),
    sprintf("  >=50%% reciprocal overlap: %d / %d (%.1f%%)",
            res$n_recip, res$n, res$pct_recip),
    sprintf("  RiboTaper uORFs recovered: %d / %d (%.1f%%)",
            res$n_rt_hit, nrow(rt), res$pct_rt_hit)
  )
}

lines <- c(
  "=== Ribostan vs RiboTaper HEK293 uORF comparison ===",
  "Reference: Calviello et al. 2016 (GSE73136 / GSM1887643)",
  "Ribo-seq: SRR2433794 (~32M reads, 27-30 nt footprints)",
  "Annotation: Gencode v19 (hg19/GRCh37), matching original RiboTaper analysis",
  "",
  sprintf("RiboTaper published uORFs (filtered BED): %d", nrow(rt)),
  sprintf("  unique transcripts: %d", length(unique(rt$transcript_id))),
  sprintf("  median width (nt): %.0f", median(rt$width_nt)),
  "",
  sprintf("Ribostan sequence candidates (pre-filter): %s",
          if (file.exists(out_path("ribostan_hek293_candidates.tsv"))) {
            as.character(nrow(read.delim(out_path("ribostan_hek293_candidates.tsv"),
                                         stringsAsFactors = FALSE)))
          } else if (file.exists(out_path("ribostan_hek293_candidate_count.txt"))) {
            readLines(out_path("ribostan_hek293_candidate_count.txt"), n = 1L)
          } else "n/a"),
  sprintf("Ribostan translation-supported uORFs (%s, n_psites >= 6): %d",
          sig_label, nrow(rs_all)),
  if (has_overlap_col) {
    sprintf("  of which overlap the annotated CDS (overlaps_cds = TRUE): %d",
            sum(!is.na(rs_all$overlaps_cds) & rs_all$overlaps_cds))
  } else "  (overlaps_cds not available in this run's output)",
  "",
  "--- Agreement: all periodic uORFs (Ribostan's definition) ---",
  fmt_result(res_all),
  "",
  if (has_overlap_col) c(
    "--- Agreement: RiboTaper-matched definition (excludes uORFs overlapping annotated CDS) ---",
    fmt_result(res_nonoverlap)
  ) else "(overlaps_cds column absent - RiboTaper-matched breakdown skipped)",
  "",
  "Notes:",
  "- RiboTaper uORFs are translation-supported (periodic RPF + ORF reconstruction).",
  "- RiboTaper's own category rule excludes uORFs that overlap any annotated CDS",
  "  exon (Calviello et al. 2016 Methods); those get bucketed as ORFs_ccds instead.",
  "  Ribostan's uORF finder does not currently apply this exclusion, so the",
  "  'RiboTaper-matched' block above is the fairer like-for-like comparison.",
  "- RiboTaper ORF_id_tr uses transcript coordinates; genomic BED used for overlap.",
  sprintf("- Annotation: %s%s",
          run_settings[["gtf_name"]],
          if (!is.na(run_settings["gtf_source"]) && nzchar(run_settings["gtf_source"])) {
            sprintf(" (RiboTaper filter on %s: use_ccds=%s use_appris=%s)",
                    run_settings["gtf_source"],
                    run_settings["gtf_use_ccds"],
                    run_settings["gtf_use_appris"])
          } else if (!is.na(run_settings["gtf_use_ccds"])) {
            ""
          } else " (explicit HEK293_GTF_NAME override)"),
  sprintf("- Significance test: %s", sig_label),
  sprintf("- Nested-uORF selection rule: %s",
          if (identical(run_settings[["nested_rule"]], "ribotaper")) {
            "RiboTaper's upstream-ATG-with-signal rule"
          } else "density-based preselection"),
  sprintf("- Multi-mapping read handling: %s", run_settings[["multimap"]]),
  sprintf("- In-frame-majority requirement (min_frame_frac): %s",
          if (identical(run_settings[["min_frame_frac"]], "0")) {
            "disabled"
          } else paste0("> ", run_settings[["min_frame_frac"]],
                        " (RiboTaper's Phase_P_sites_frame/P_sites_sum gate)")),
  sprintf("- P-site offset calibration: %s",
          if (identical(run_settings[["offset_method"]], "metagene_peak")) {
            "RiboTaper's own metagene-peak method (get_metagene_peak_offsets)"
          } else "Ribostan's default CDS-inclusion-maximising method (cds_incl)"),
  sprintf("- Multi-mapping/coverage post-hoc QC (flag_multimapping_uorfs): %s",
          if (isTRUE(as.logical(run_settings[["multimap_qc"]]))) {
            "applied (RiboTaper's pct_covered_onlymulti/pct_covered_all < 0.3 AND pct_covered_all > 0.3)"
          } else "not applied"),
  sprintf("- Read alignment source: %s",
          if (identical(run_settings[["read_source"]], "genome_primary")) {
            "RiboTaper-style: genome-aligned STAR + primary-alignment-only (samtools -F 0x100), projected to transcript space"
          } else "Ribostan's default: direct transcriptome-aligned BAM")
)

writeLines(lines, out_path("hek293_comparison_summary.txt"))
cat(paste(lines, collapse = "\n"), "\n")

if (nrow(res_all$overlaps) > 0) {
  write.table(res_all$overlaps, out_path("hek293_comparison_overlaps.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}
if (nrow(res_nonoverlap$overlaps) > 0) {
  write.table(res_nonoverlap$overlaps, out_path("hek293_comparison_overlaps_nonCDSoverlap.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}
