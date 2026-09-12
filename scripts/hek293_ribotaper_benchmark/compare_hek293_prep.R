#!/usr/bin/env Rscript
# Run Ribostan on HEK293 Ribo-seq (SRR2433794) for comparison to RiboTaper uORFs.
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
bench_dir <- if (length(script_path)) dirname(script_path) else "."
workdir <- {
  tra <- commandArgs(trailingOnly = TRUE)
  if (length(tra) >= 1) tra[[1]] else "."
}
setwd(workdir)

source(file.path(bench_dir, "load_ribotaper_uorfs.R"), local = TRUE)
source(file.path(bench_dir, "hek293_gtf_config.R"), local = TRUE)

suppressPackageStartupMessages({
  library(Ribostan)
  library(GenomicRanges)
  library(IRanges)
  library(magrittr)
  library(dplyr)
})

out_suffix <- Sys.getenv("HEK293_OUT_SUFFIX", "")
out_path <- function(name) file.path(workdir, sub("\\.", paste0(out_suffix, "."), name))

gtf_cfg <- resolve_hek293_gtf(workdir)
gtf <- gtf_cfg$path
gtf_name <- gtf_cfg$name
fafile <- file.path(workdir, "GRCh37.p13.genome.fa")
bam <- file.path(workdir, "HEK293_ribo_transcriptome.bam")
bed_gz <- file.path(workdir, "GSM1887643_HEK293_filtered_ORFs.bed.gz")

# Methodology knobs for aligning with RiboTaper's own published pipeline
# (Calviello et al. 2016): raw-p vs BH-FDR significance, RiboTaper's
# upstream-ATG-with-signal nested-uORF rule, multi-mapping read handling,
# RiboTaper's own in-frame-majority requirement for a "periodic" call, and
# RiboTaper's own metagene-peak P-site offset calibration procedure. See
# R/periodicity_tests.R, R/get_ritpms.R and R/get_offsets.R for the
# underlying options.
fdr_correct <- as.logical(Sys.getenv("HEK293_FDR_CORRECT", "FALSE"))
nested_rule <- Sys.getenv("HEK293_NESTED_RULE", "ribotaper")
multimap <- Sys.getenv("HEK293_MULTIMAP", "all")
min_frame_frac <- as.numeric(Sys.getenv("HEK293_MIN_FRAME_FRAC", "0.5"))
offset_method <- Sys.getenv("HEK293_OFFSET_METHOD", "metagene_peak")
n_cores <- as.integer(Sys.getenv("HEK293_N_CORES", "1"))
multimap_qc <- as.logical(Sys.getenv("HEK293_MULTIMAP_QC", "FALSE"))
apply_phase_shiftdf <- as.logical(Sys.getenv("HEK293_APPLY_PHASE_SHIFT", "FALSE"))
# "transcriptome_bam" (default): Ribostan's normal path, get_readgr() on a
# transcriptome-aligned BAM. "genome_primary": RiboTaper-style alignment -
# reads pre-projected (convert_genome_to_txspace.R) from a genome-aligned,
# primary-alignment-only BAM (STAR + samtools -F 0x100), matching how
# GSM1887643 itself was processed, then broadcast to every transcript whose
# exons they land in (equivalent to multimap="all" in transcript space).
read_source <- Sys.getenv("HEK293_READ_SOURCE", "genome_primary")
projected_rds_name <- Sys.getenv("HEK293_PROJECTED_RDS",
                                 "rpfs_genome_projected_rrnafilt_matched.rds")
stopifnot(!is.na(fdr_correct), nested_rule %in% c("density", "ribotaper"),
          multimap %in% c("all", "primary", "unique"),
          !is.na(min_frame_frac), min_frame_frac >= 0, min_frame_frac < 1,
          offset_method %in% c("cds_incl", "metagene_peak"),
          !is.na(n_cores), n_cores >= 1, !is.na(multimap_qc),
          !is.na(apply_phase_shiftdf),
          read_source %in% c("transcriptome_bam", "genome_primary"))
message("Settings: fdr_correct=", fdr_correct, " nested_rule=", nested_rule,
        " multimap=", multimap, " min_frame_frac=", min_frame_frac,
        " offset_method=", offset_method, " n_cores=", n_cores,
        " multimap_qc=", multimap_qc, " apply_phase_shiftdf=", apply_phase_shiftdf,
        " read_source=", read_source,
        " gtf=", gtf_name,
        if (!is.na(gtf_cfg$use_ccds)) {
          paste0(" (source=", gtf_cfg$source,
                 " use_ccds=", gtf_cfg$use_ccds,
                 " use_appris=", gtf_cfg$use_appris, ")")
        } else {
          " (explicit HEK293_GTF_NAME)"
        })
writeLines(
  c(paste0("fdr_correct=", fdr_correct), paste0("nested_rule=", nested_rule),
    paste0("multimap=", multimap), paste0("min_frame_frac=", min_frame_frac),
    paste0("offset_method=", offset_method), paste0("multimap_qc=", multimap_qc),
    paste0("apply_phase_shiftdf=", apply_phase_shiftdf),
    paste0("read_source=", read_source),
    paste0("projected_rds_name=", projected_rds_name),
    paste0("gtf_name=", gtf_name),
    paste0("gtf_source=", ifelse(is.na(gtf_cfg$source), "", gtf_cfg$source)),
    paste0("gtf_use_ccds=", ifelse(is.na(gtf_cfg$use_ccds), "", gtf_cfg$use_ccds)),
    paste0("gtf_use_appris=", ifelse(is.na(gtf_cfg$use_appris), "", gtf_cfg$use_appris))),
  out_path("ribostan_hek293_settings.txt")
)

stopifnot(file.exists(gtf), file.exists(fafile), file.exists(bed_gz))
if (read_source == "transcriptome_bam") {
  stopifnot(file.exists(bam))
}

message("Loading Gencode v19 annotation + uORF candidates ...")
anno <- load_annotation(gtf, fafile, add_uorfs = TRUE)

uorf_to_df <- function(anno_obj, ids, label) {
  ids <- intersect(ids, names(anno_obj$cdsgrl))
  if (length(ids) == 0L) {
    return(data.frame(
      set = character(), orf_id = character(), transcript_id = character(),
      seqnames = character(), start = integer(), end = integer(),
      strand = character(), tx_start = integer(), tx_end = integer(),
      width_nt = integer(), width_codons = integer(),
      p_value = numeric(), q_value = numeric(), n_psites = integer(),
      pct_inframe = numeric(),
      overlaps_cds = logical(), out_of_phase = logical(),
      stringsAsFactors = FALSE
    ))
  }
  # cdsgrl holds real genomic coordinates (chromosome, possibly split across
  # exons for spliced ORFs); trspacecds holds transcript-relative coordinates
  # (its "seqnames" is the transcript ID, not a chromosome). Both are needed:
  # genomic for comparing against RiboTaper's BED intervals, transcript-space
  # for comparing against RiboTaper's ORF_id_tr-encoded tx_start/tx_end.
  cds_grl <- anno_obj$cdsgrl[ids]
  tr_orfs <- anno_obj$trspacecds[ids]
  mcol_num <- function(gr, col) {
    vals <- mcols(gr)[[col]]
    if (is.null(vals)) rep(NA_real_, length(gr)) else as.numeric(vals)
  }
  mcol_int <- function(gr, col) {
    vals <- mcols(gr)[[col]]
    if (is.null(vals)) rep(NA_integer_, length(gr)) else as.integer(vals)
  }
  mcol_lgl <- function(gr, col) {
    vals <- mcols(gr)[[col]]
    if (is.null(vals)) rep(NA, length(gr)) else as.logical(vals)
  }
  data.frame(
    set = label,
    orf_id = ids,
    transcript_id = sub("_[0-9]+$", "", ids),
    seqnames = vapply(cds_grl, function(g) as.character(GenomicRanges::seqnames(g)[1]), character(1)),
    start = vapply(cds_grl, function(g) min(GenomicRanges::start(g)), integer(1)),
    end = vapply(cds_grl, function(g) max(GenomicRanges::end(g)), integer(1)),
    strand = vapply(cds_grl, function(g) as.character(GenomicRanges::strand(g)[1]), character(1)),
    tx_start = as.integer(GenomicRanges::start(tr_orfs)),
    tx_end = as.integer(GenomicRanges::end(tr_orfs)),
    width_nt = as.integer(GenomicRanges::width(tr_orfs)),
    width_codons = (as.integer(GenomicRanges::width(tr_orfs)) - 1L) / 3L,
    p_value = mcol_num(tr_orfs, "p.value"),
    q_value = mcol_num(tr_orfs, "q.value"),
    n_psites = mcol_int(tr_orfs, "n_psites"),
    # Fraction of P-sites in the ORF's own reading frame (RiboTaper's
    # Phase_P_sites_frame / P_sites_sum gate) - see min_frame_frac below.
    pct_inframe = mcol_num(tr_orfs, "pct_inframe"),
    # Matches RiboTaper's own uORF/ORFs_ccds split: "uORFs ... non
    # overlapping with any CDS exon". overlaps_cds is computed by
    # load_annotation()'s .annotate_uorf_overlaps() and is FALSE/NA for
    # non-uORF ORFs, so it is only meaningful for uORF rows.
    overlaps_cds = mcol_lgl(tr_orfs, "overlaps_cds"),
    out_of_phase = mcol_lgl(tr_orfs, "out_of_phase"),
    stringsAsFactors = FALSE
  )
}

candidate_ids <- names(anno$uORF)[anno$uORF]
writeLines(as.character(length(candidate_ids)), out_path("ribostan_hek293_candidate_count.txt"))
message("Ribostan sequence candidates: ", length(candidate_ids))
candidates_df <- uorf_to_df(anno, candidate_ids, "sequence_candidate")
write.table(candidates_df, out_path("ribostan_hek293_candidates_full.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

if (read_source == "genome_primary") {
  rpfs_rds <- file.path(workdir, Sys.getenv("HEK293_PROJECTED_RDS", "rpfs_genome_projected.rds"))
  stopifnot(file.exists(rpfs_rds))
  message("Loading RiboTaper-style genome-primary reads pre-projected to ",
          "transcript space (", rpfs_rds, ") ...")
  rpfs <- readRDS(rpfs_rds)
} else {
  message("Reading Ribo-seq BAM (multimap = ", multimap, ") ...")
  rpfs <- get_readgr(bam, anno, multimap = multimap)
}
message("RPFs: ", length(rpfs), " ; read lengths: ",
        paste(sort(unique(rpfs$readlen)), collapse = ", "))

message("Estimating P-site offsets (method = ", offset_method, ") ...")
offsets_df <- get_offsets(rpfs, anno, method = offset_method)
write.table(offsets_df, out_path("ribostan_offsets.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

psites <- get_psite_gr(rpfs, offsets_df, anno, apply_phase_shiftdf = apply_phase_shiftdf)
message("P-sites: ", length(psites))

message("Periodicity testing all uORF candidates (remove = FALSE, keeps full stats) ...")
alpha <- 0.05
min_psites <- 6L
tested_anno <- periodicity_filter_uORFs(psites, anno, remove = FALSE,
                                         alpha = alpha, min_psites = min_psites,
                                         n_cores = n_cores, nested_rule = nested_rule,
                                         fdr_correct = fdr_correct,
                                         min_frame_frac = min_frame_frac)

# Full per-candidate stats table: every sequence-candidate uORF, whether
# tested-and-passed, tested-and-failed, pruned as a redundant nested
# candidate (p/q NA), or never reached the multitaper test at all. This is
# what lets us diagnose *why* a specific RiboTaper-published uORF wasn't
# recovered by Ribostan, rather than just observing that it wasn't.
all_tested_df <- uorf_to_df(tested_anno, candidate_ids, "sequence_candidate_tested")
write.table(all_tested_df, out_path("ribostan_hek293_all_tested.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Same pass/fail logic as periodicity_filter_uORFs(remove = TRUE), applied
# here to the already-tested stats so the (slow) multitaper test does not
# need to be run twice. sig_col follows fdr_correct: RiboTaper's own methods
# paper thresholds on the raw per-ORF p-value (p < 0.05), not a genome-wide
# BH-corrected q-value, so fdr_correct=FALSE reproduces that directly.
# min_frame_frac mirrors RiboTaper's own in-frame-majority gate (disabled
# by default; set HEK293_MIN_FRAME_FRAC=0.5 to match RiboTaper exactly).
sig_col <- if (fdr_correct) "q_value" else "p_value"
periodic_ids <- all_tested_df$orf_id[
  !is.na(all_tested_df[[sig_col]]) & all_tested_df[[sig_col]] < alpha &
    !is.na(all_tested_df$n_psites) & all_tested_df$n_psites >= min_psites &
    (is.na(all_tested_df$pct_inframe) | all_tested_df$pct_inframe > min_frame_frac)
]

# RiboTaper's own post-hoc multi-mapping/coverage QC (Calviello et al. 2016):
# drops ORFs whose full-read-footprint coverage is dominated by
# multi-mapping reads and/or too sparse overall. Uses `rpfs` (full aligned
# read spans), not `psites` (single-nt P-sites) - RiboTaper's own
# pct_region_covered_ribo comes from bedtools coverage on whole reads, so
# using P-sites here would understate coverage and over-filter. Only
# meaningful with multimap != "unique" upstream (otherwise no
# multi-mapping reads survive to evaluate).
if (multimap_qc) {
  message("Applying RiboTaper's multi-mapping/coverage QC (threshold=0.3, min_coverage=0.3) ...")
  mm_qc <- flag_multimapping_uorfs(rpfs, tested_anno, periodic_ids,
                                    threshold = 0.3, min_coverage = 0.3)
  write.table(mm_qc, out_path("ribostan_hek293_multimap_qc.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  message("Multi-mapping/coverage QC: ", sum(mm_qc$flagged), " / ",
          nrow(mm_qc), " periodic uORFs flagged and dropped")
  periodic_ids <- mm_qc$orf_id[!mm_qc$flagged]
}
periodic_df <- uorf_to_df(tested_anno, periodic_ids, "periodic")

write.table(periodic_df, out_path("ribostan_hek293_periodic.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

rt_uorfs <- load_ribotaper_hek293_uorfs(bed_gz)
write.table(rt_uorfs, out_path("ribotaper_hek293_uorfs.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(rt_uorfs, out_path("ribotaper_hek293_uorfs.rds"))

message("Ribostan candidates: ", length(candidate_ids))
message("Ribostan translation-supported (", sig_col, "<", alpha,
        ", n_psites>=", min_psites, "): ", nrow(periodic_df))
message("RiboTaper published uORFs: ", nrow(rt_uorfs))
