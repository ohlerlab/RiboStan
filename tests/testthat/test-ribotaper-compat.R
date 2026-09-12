# Tests for features added to align Ribostan's periodicity filtering more
# closely with RiboTaper's own published methodology (Calviello et al.
# 2016): raw (uncorrected) p-value filtering, RiboTaper's own
# upstream-ATG-with-signal rule for nested uORFs, and its multi-mapping
# read QC exclusion. See R/periodicity_tests.R and R/get_ritpms.R.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})


# ---------------------------------------------------------------------------
# select_ribotaper_nested_uorfs: RiboTaper's "most upstream ATG with >5
# P-sites (>50% in-frame) between it and its closest downstream neighbour"
# rule, applied to already-tested candidates sharing a stop codon.
# ---------------------------------------------------------------------------

#' Build a minimal trspacecds-shaped GRanges for synthetic nested-uORF tests
#' @keywords internal
.mk_trspacecds <- function(tx, starts, ends, ids) {
  gr <- GenomicRanges::GRanges(
    seqnames = tx,
    ranges   = IRanges::IRanges(start = starts, end = ends),
    strand   = "+"
  )
  names(gr) <- ids
  gr
}

test_that("select_ribotaper_nested_uorfs picks the upstream ATG when its gap has signal", {
  # Three candidates on tx1 sharing stop=30, spaced 9nt apart (in-frame):
  # c1 (start=1, not significant), c2 (start=10, significant),
  # c3 (start=19, significant). c2 should win because the gap [10,18] to
  # its neighbour c3 has >5 psites, >50% in-frame (positions 10,13,16).
  orfs <- .mk_trspacecds("tx1", c(1L, 10L, 19L), c(30L, 30L, 30L),
                        c("c1", "c2", "c3"))
  anno <- list(trspacecds = orfs)

  tested_df <- data.frame(
    orf_id = c("c1", "c2", "c3"),
    p.value = c(0.20, 0.01, 0.02),
    stringsAsFactors = FALSE
  )

  # In-frame positions relative to gap_start=10 are 10, 13, 16 (step 3).
  # Give them 3 psites each (9 in-frame) plus 2 out-of-frame psites so the
  # gap has 11 psites, 9/11 > 50% in-frame, clearing the ">5 psites" bar.
  gap_pos <- c(rep(c(10L, 13L, 16L), 3), 11L, 12L)
  psites <- GenomicRanges::GRanges(
    seqnames = "tx1",
    ranges   = IRanges::IRanges(start = gap_pos, width = 1L),
    strand   = "+"
  )
  names(psites) <- paste0("read", seq_along(psites))

  winners <- select_ribotaper_nested_uorfs(psites, anno,
                                           c("c1", "c2", "c3"), tested_df)
  expect_equal(winners, "c2")
})


test_that("select_ribotaper_nested_uorfs falls back to lowest p-value when no gap qualifies", {
  # Same layout, but no psites at all in the c2->c3 gap: RiboTaper falls
  # back to the single most-significant candidate in the whole group.
  orfs <- .mk_trspacecds("tx1", c(1L, 10L, 19L), c(30L, 30L, 30L),
                        c("c1", "c2", "c3"))
  anno <- list(trspacecds = orfs)

  tested_df <- data.frame(
    orf_id = c("c1", "c2", "c3"),
    p.value = c(0.20, 0.01, 0.001),
    stringsAsFactors = FALSE
  )
  # No psites anywhere - the gap check for c2 always fails (n_gap == 0).
  psites <- GenomicRanges::GRanges(
    seqnames = character(0),
    ranges   = IRanges::IRanges(start = integer(0), width = integer(0))
  )
  names(psites) <- character(0)

  winners <- select_ribotaper_nested_uorfs(psites, anno,
                                           c("c1", "c2", "c3"), tested_df)
  expect_equal(winners, "c3")
})


test_that("select_ribotaper_nested_uorfs falls back when fewer than 2 candidates are significant", {
  orfs <- .mk_trspacecds("tx1", c(1L, 10L), c(30L, 30L), c("c1", "c2"))
  anno <- list(trspacecds = orfs)

  # Only c2 clears p < 0.05, so there's nothing to compare a gap against;
  # RiboTaper keeps the single significant (== best) candidate.
  tested_df <- data.frame(
    orf_id = c("c1", "c2"),
    p.value = c(0.20, 0.03),
    stringsAsFactors = FALSE
  )
  psites <- GenomicRanges::GRanges(
    seqnames = character(0),
    ranges   = IRanges::IRanges(start = integer(0), width = integer(0))
  )
  names(psites) <- character(0)

  winners <- select_ribotaper_nested_uorfs(psites, anno,
                                           c("c1", "c2"), tested_df)
  expect_equal(winners, "c2")
})


test_that("select_ribotaper_nested_uorfs(min_frame_frac) nulls out candidates that fail the frame gate", {
  # Same 3-candidate layout as the "gap has signal" test above, but c2 -
  # the candidate that would otherwise win via the gap-signal rule - is
  # tagged as majority out-of-frame (pct_inframe = 0.3). RiboTaper would
  # never have assigned it a usable p-value in the first place, so with
  # min_frame_frac = 0.5 it must lose its spot to c3, the next-best
  # candidate by raw p-value among the remaining (in-frame) group members.
  orfs <- .mk_trspacecds("tx1", c(1L, 10L, 19L), c(30L, 30L, 30L),
                        c("c1", "c2", "c3"))
  anno <- list(trspacecds = orfs)

  tested_df <- data.frame(
    orf_id = c("c1", "c2", "c3"),
    p.value = c(0.20, 0.01, 0.02),
    pct_inframe = c(0.9, 0.3, 0.9),
    stringsAsFactors = FALSE
  )
  gap_pos <- c(rep(c(10L, 13L, 16L), 3), 11L, 12L)
  psites <- GenomicRanges::GRanges(
    seqnames = "tx1",
    ranges   = IRanges::IRanges(start = gap_pos, width = 1L),
    strand   = "+"
  )
  names(psites) <- paste0("read", seq_along(psites))

  # Without the gate, c2 wins exactly as in the ungated test above.
  winners_default <- select_ribotaper_nested_uorfs(psites, anno,
                                                    c("c1", "c2", "c3"),
                                                    tested_df)
  expect_equal(winners_default, "c2")

  # With the gate, c2's p-value is nulled out, leaving only c3 significant
  # (< 2 significant candidates), so the group falls back to the lowest
  # raw p-value *among candidates that still have one* - c3.
  winners_gated <- select_ribotaper_nested_uorfs(psites, anno,
                                                 c("c1", "c2", "c3"),
                                                 tested_df,
                                                 min_frame_frac = 0.5)
  expect_equal(winners_gated, "c3")
})


test_that("select_ribotaper_nested_uorfs excludes candidates below the P-site gate", {
  orfs <- .mk_trspacecds("tx1", c(1L, 10L), c(30L, 30L), c("c1", "c2"))
  anno <- list(trspacecds = orfs)
  tested_df <- data.frame(
    orf_id = c("c1", "c2"),
    p.value = c(0.01, 0.001),
    n_psites = c(8L, 5L),
    stringsAsFactors = FALSE
  )
  psites <- GenomicRanges::GRanges(
    seqnames = character(0),
    ranges = IRanges::IRanges(start = integer(0), width = integer(0))
  )
  names(psites) <- character(0)

  winners <- select_ribotaper_nested_uorfs(
    psites,
    anno,
    c("c1", "c2"),
    tested_df,
    min_psites = 6L
  )
  expect_equal(winners, "c1")
})


test_that("select_ribotaper_nested_uorfs passes through singleton groups untouched", {
  orfs <- .mk_trspacecds(c("tx1", "tx2"), c(1L, 5L), c(30L, 40L),
                        c("solo1", "solo2"))
  anno <- list(trspacecds = orfs)
  tested_df <- data.frame(
    orf_id = c("solo1", "solo2"),
    p.value = c(0.9, NA_real_),
    stringsAsFactors = FALSE
  )
  psites <- GenomicRanges::GRanges(
    seqnames = character(0),
    ranges   = IRanges::IRanges(start = integer(0), width = integer(0))
  )
  names(psites) <- character(0)

  winners <- select_ribotaper_nested_uorfs(psites, anno,
                                           c("solo1", "solo2"), tested_df)
  expect_setequal(winners, c("solo1", "solo2"))
})


# ---------------------------------------------------------------------------
# periodicity_filter_uORFs: fdr_correct toggles q.value vs raw p.value as
# the significance column used for the remove = TRUE decision.
# ---------------------------------------------------------------------------
test_that("periodicity_filter_uORFs(fdr_correct = FALSE) filters on raw p.value", {
  data(chr22_anno)
  data(rpfs)
  data(offsets_df)
  psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)

  filt_raw <- suppressMessages(
    periodicity_filter_uORFs(psites, chr22_anno, remove = FALSE,
                             alpha = 0.05, n_cores = 1, fdr_correct = FALSE)
  )
  m <- S4Vectors::mcols(filt_raw$trspacecds)
  uorf_ids <- names(filt_raw$uORF)[filt_raw$uORF]
  tested <- uorf_ids[!is.na(m[uorf_ids, "p.value"])]
  skip_if(length(tested) < 1L, "No tested uORFs in fixture data")

  # Raw p-value is never stricter than the BH-corrected q-value, so the
  # raw-p-selected passing set must be a superset of the q-value-selected
  # passing set for the same alpha.
  passes_p <- tested[m[tested, "p.value"] < 0.05]
  passes_q <- tested[!is.na(m[tested, "q.value"]) & m[tested, "q.value"] < 0.05]
  expect_true(all(passes_q %in% passes_p))

  filt_hi <- suppressMessages(
    periodicity_filter_uORFs(psites, chr22_anno, remove = TRUE,
                             alpha = 0.05, n_cores = 1, fdr_correct = TRUE)
  )
  filt_lo <- suppressMessages(
    periodicity_filter_uORFs(psites, chr22_anno, remove = TRUE,
                             alpha = 0.05, n_cores = 1, fdr_correct = FALSE)
  )
  # Raw-p filtering is always at least as permissive as BH-corrected
  # filtering at the same alpha.
  expect_lte(sum(filt_hi$uORF), sum(filt_lo$uORF))
})


# ---------------------------------------------------------------------------
# ftest_orfs: pct_inframe, matching RiboTaper's Phase_P_sites_frame /
# P_sites_sum gate (CCDS_orf_finder.R / NONCCDS_orf_finder.R).
# ---------------------------------------------------------------------------
test_that("ftest_orfs computes pct_inframe as the fraction of P-sites in the ORF's own frame", {
  orf <- GenomicRanges::GRanges(seqnames = "tx1",
                                ranges = IRanges::IRanges(start = 1, end = 30),
                                strand = "+")
  names(orf) <- "orfA"
  anno <- list(trspacecds = orf, cdsgrl = orf)

  # In-frame (relative positions 1, 4, 7, 10 - frame 0, since the ORF's own
  # start maps to relative position 1): 4 psites. Out-of-frame (relative
  # positions 2, 3): 2 psites. pct_inframe = 4 / 6.
  positions <- c(1L, 4L, 7L, 10L, 2L, 3L)
  psites <- GenomicRanges::GRanges(
    seqnames = "tx1",
    ranges = IRanges::IRanges(start = positions, width = 1L),
    strand = "+",
    orf = "orfA"
  )
  names(psites) <- paste0("read", seq_along(psites))

  res <- suppressMessages(ftest_orfs(psites, anno, n_cores = 1))
  expect_equal(res$pct_inframe[res$orf_id == "orfA"], 4 / 6, tolerance = 1e-8)
})


# ---------------------------------------------------------------------------
# periodicity_filter_uORFs: min_frame_frac requires a majority of an ORF's
# P-sites to sit in its own reading frame, matching RiboTaper's gate.
# ---------------------------------------------------------------------------
test_that("periodicity_filter_uORFs attaches pct_inframe as an mcol", {
  data(chr22_anno)
  data(rpfs)
  data(offsets_df)
  psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)

  filt <- suppressMessages(
    periodicity_filter_uORFs(psites, chr22_anno, remove = FALSE, n_cores = 1)
  )
  m <- S4Vectors::mcols(filt$trspacecds)
  expect_true("pct_inframe" %in% colnames(m))
  uorf_ids <- names(filt$uORF)[filt$uORF]
  tested <- uorf_ids[!is.na(m[uorf_ids, "p.value"])]
  skip_if(length(tested) < 1L, "No tested uORFs in fixture data")
  expect_true(all(m[tested, "pct_inframe"] >= 0 & m[tested, "pct_inframe"] <= 1))
})


test_that("periodicity_filter_uORFs(min_frame_frac = 0.5) is at least as strict as the default", {
  data(chr22_anno)
  data(rpfs)
  data(offsets_df)
  psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)

  filt_default <- suppressMessages(
    periodicity_filter_uORFs(psites, chr22_anno, remove = TRUE,
                             alpha = 0.05, n_cores = 1, fdr_correct = FALSE)
  )
  filt_gated <- suppressMessages(
    periodicity_filter_uORFs(psites, chr22_anno, remove = TRUE,
                             alpha = 0.05, n_cores = 1, fdr_correct = FALSE,
                             min_frame_frac = 0.5)
  )
  expect_true(all(names(filt_gated$uORF)[filt_gated$uORF] %in%
                    names(filt_default$uORF)[filt_default$uORF]))
})


test_that("periodicity_filter_uORFs defaults to the RiboTaper nested rule", {
  expect_identical(
    eval(formals(periodicity_filter_uORFs)$nested_rule),
    c("ribotaper", "density")
  )

  data(chr22_anno)
  data(rpfs)
  data(offsets_df)
  psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)

  filt <- suppressMessages(
    periodicity_filter_uORFs(psites, chr22_anno, remove = TRUE,
                             alpha = 0.05, n_cores = 1,
                             nested_rule = "ribotaper", fdr_correct = FALSE)
  )
  expect_true(isTRUE(filt$uORF_periodicity_filtered))
  expect_true(sum(filt$uORF) <= sum(chr22_anno$uORF))
})


# ---------------------------------------------------------------------------
# flag_multimapping_uorfs: RiboTaper's multi-mapping-only coverage QC rule
# (pct_covered_onlymulti / pct_covered_all >= threshold -> flagged), computed
# from full read-footprint coverage (bedtools-coverage-style), not P-sites.
# ---------------------------------------------------------------------------
test_that("flag_multimapping_uorfs flags ORFs covered mostly by multi-mapping reads", {
  orfs <- .mk_trspacecds(c("tx1", "tx2"), c(1L, 1L), c(20L, 20L),
                        c("orfA", "orfB"))
  anno <- list(trspacecds = orfs)

  # orfA: two adjacent 5nt reads cover positions 1-10, all NH > 1
  # (multi-mapping only) -> flagged. orfB: same footprint on positions
  # 1-10, all NH == 1 (unique) -> not flagged.
  reads <- GenomicRanges::GRanges(
    seqnames = c("tx1", "tx1", "tx2", "tx2"),
    ranges = IRanges::IRanges(start = c(1L, 6L, 1L, 6L), width = 5L),
    strand = "+",
    NH = c(2L, 2L, 1L, 1L)
  )
  names(reads) <- paste0("read", seq_along(reads))

  res <- flag_multimapping_uorfs(reads, anno, c("orfA", "orfB"),
                                 threshold = 0.3)
  rownames(res) <- res$orf_id

  expect_equal(res["orfA", "pct_covered_all"], 0.5, tolerance = 1e-8)
  expect_equal(res["orfA", "pct_covered_onlymulti"], 0.5, tolerance = 1e-8)
  expect_true(res["orfA", "flagged"])
  expect_equal(res["orfB", "pct_covered_onlymulti"], 0, tolerance = 1e-8)
  expect_false(res["orfB", "flagged"])
})


test_that("flag_multimapping_uorfs also flags ORFs with insufficient coverage breadth", {
  orfs <- .mk_trspacecds(c("tx1", "tx2"), c(1L, 1L), c(20L, 20L),
                        c("orfA", "orfB"))
  anno <- list(trspacecds = orfs)

  # orfA: a single 2nt read covers only 2/20 positions (10% breadth), NH==1
  # -> multimap ratio is 0 (would pass threshold alone) but breadth is too
  # low. orfB: a 10nt read covers 10/20 positions (50% breadth), NH==1 ->
  # passes both.
  reads <- GenomicRanges::GRanges(
    seqnames = c("tx1", "tx2"),
    ranges = IRanges::IRanges(start = c(1L, 1L), width = c(2L, 10L)),
    strand = "+",
    NH = 1L
  )
  names(reads) <- paste0("read", seq_along(reads))

  res <- flag_multimapping_uorfs(reads, anno, c("orfA", "orfB"),
                                 threshold = 0.3, min_coverage = 0.3)
  rownames(res) <- res$orf_id

  expect_true(res["orfA", "low_coverage"])
  expect_true(res["orfA", "flagged"])
  expect_false(res["orfB", "low_coverage"])
  expect_false(res["orfB", "flagged"])
})


test_that("flag_multimapping_uorfs handles ORFs with zero coverage without erroring", {
  orfs <- .mk_trspacecds("tx1", 1L, 20L, "orfA")
  anno <- list(trspacecds = orfs)
  reads <- GenomicRanges::GRanges(
    seqnames = character(0),
    ranges   = IRanges::IRanges(start = integer(0), width = integer(0)),
    NH  = integer(0)
  )
  names(reads) <- character(0)

  res <- flag_multimapping_uorfs(reads, anno, "orfA")
  expect_equal(res$pct_covered_all, 0)
  expect_false(res$flagged)
})


# ---------------------------------------------------------------------------
# get_psite_gr: apply_phase_shiftdf=FALSE skips Ribostan's CDS phase reranking
# (RiboTaper uses metagene offset only).
# ---------------------------------------------------------------------------
test_that("get_psite_gr(apply_phase_shiftdf = FALSE) skips CDS phase reranking", {
  orfs <- GenomicRanges::GRanges(
    seqnames = "tx1",
    ranges = IRanges::IRanges(1L, 60L),
    strand = "+"
  )
  names(orfs) <- "cds1"
  GenomeInfoDb::seqlengths(orfs) <- 60L
  anno <- list(trspacecds = orfs, cdsgrl = orfs)

  starts <- c(2L, 5L, 8L, 11L, 1L, 3L, 6L)
  rpfs <- GenomicRanges::GRanges(
    seqnames = rep("tx1", length(starts)),
    ranges = IRanges::IRanges(starts, width = 1L),
    strand = "+",
    readlen = 1L
  )
  names(rpfs) <- paste0("read", seq_along(rpfs))
  GenomeInfoDb::seqlengths(rpfs) <- 60L
  offsets <- data.frame(
    readlen = rep(1L, 3L),
    phase = 0:2,
    p_offset = 0L
  )

  ps_with <- get_psite_gr(rpfs, offsets, anno, apply_phase_shiftdf = TRUE)
  ps_without <- get_psite_gr(rpfs, offsets, anno, apply_phase_shiftdf = FALSE)
  expect_equal(length(ps_with), length(ps_without))
  expect_false(all(start(ps_with) == start(ps_without)))
})

test_that("get_psite_gr reranks secondary phases when periodicity improves", {
  orfs <- GenomicRanges::GRanges(
    seqnames = "tx1",
    ranges = IRanges::IRanges(1L, 60L),
    strand = "+"
  )
  names(orfs) <- "cds1"
  GenomeInfoDb::seqlengths(orfs) <- 60L
  anno <- list(trspacecds = orfs, cdsgrl = orfs)

  starts_first <- c(1L, 4L, 7L, 10L, 2L, 5L, 3L)
  starts_second <- c(1L, 4L, 7L, 10L, 2L, 3L, 6L)
  starts <- c(starts_first, starts_second)
  rpfs <- GenomicRanges::GRanges(
    seqnames = rep("tx1", length(starts)),
    ranges = IRanges::IRanges(starts, width = 1L),
    strand = "+",
    readlen = rep(1:2, each = 7L)
  )
  names(rpfs) <- paste0("read", seq_along(rpfs))
  GenomeInfoDb::seqlengths(rpfs) <- 60L
  offsets <- data.frame(
    readlen = rep(1:2, each = 3L),
    phase = rep(0:2, 2L),
    p_offset = 0L
  )

  gated <- get_psite_gr(rpfs, offsets, anno, apply_phase_shiftdf = TRUE)
  unshifted <- get_psite_gr(rpfs, offsets, anno, apply_phase_shiftdf = FALSE)

  unshifted_frame <- (start(unshifted) - 1L) %% 3L
  gated_frame <- (start(gated) - 1L) %% 3L
  expect_equal(as.integer(table(factor(unshifted_frame, levels = 0:2))),
               c(8L, 3L, 3L))
  expect_equal(as.integer(table(factor(gated_frame, levels = 0:2))),
               c(8L, 4L, 2L))
})


# ---------------------------------------------------------------------------
# read_ribobam / get_readgr: multimap = "all" | "unique" | "primary"
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# get_metagene_peak_offsets: RiboTaper's own "read the peak off the metagene
# plot" offset-selection procedure.
# ---------------------------------------------------------------------------
test_that("get_metagene_peak_offsets finds the modal 5'-end-to-start-codon distance per read length", {
  anno <- list(cdsstarts = c(tx1 = 101L, tx2 = 201L))

  # tx1, readlen 28: 5' ends cluster at distance -12 (start = 101 - 12 = 89),
  # with weaker decoy peaks at other distances.
  starts_28 <- c(rep(89L, 8), rep(90L, 2), rep(80L, 1))
  # tx2, readlen 30: cluster at distance -15 (start = 201 - 15 = 186).
  starts_30 <- c(rep(186L, 6), rep(187L, 1))

  rpfs <- GenomicRanges::GRanges(
    seqnames = c(rep("tx1", length(starts_28)), rep("tx2", length(starts_30))),
    ranges = IRanges::IRanges(start = c(starts_28, starts_30), width = 1L),
    strand = "+",
    readlen = c(rep(28, length(starts_28)), rep(30, length(starts_30)))
  )

  offsets_df <- get_metagene_peak_offsets(rpfs, anno)
  expect_setequal(unique(offsets_df$phase), 0:2)

  off28 <- unique(offsets_df$p_offset[offsets_df$readlen == 28])
  off30 <- unique(offsets_df$p_offset[offsets_df$readlen == 30])
  expect_equal(off28, 12)
  expect_equal(off30, 15)
  # Same offset reported for every phase - matches get_psite_gr()'s
  # (readlen, phase) lookup schema and RiboTaper's own phase-independent
  # (single-per-read-length) offset choice.
  expect_equal(length(off28), 1L)
  expect_equal(length(off30), 1L)
})


test_that("get_readgr attaches an NH mcol and 'unique' mode keeps only NH == 1 reads", {
  data(chr22_anno)
  testbam <- system.file("extdata", "nchr22.bam", package = "Ribostan",
                        mustWork = TRUE)

  rpfs_all <- get_readgr(testbam, chr22_anno, multimap = "all")
  expect_true("NH" %in% colnames(GenomicRanges::mcols(rpfs_all)))
  expect_true(all(rpfs_all$NH >= 1L))

  rpfs_unique <- get_readgr(testbam, chr22_anno, multimap = "unique")
  expect_true(all(rpfs_unique$NH == 1L))
  # "unique" can only ever drop reads relative to "all", never add any.
  expect_lte(length(rpfs_unique), length(rpfs_all))

  rpfs_primary <- get_readgr(testbam, chr22_anno, multimap = "primary")
  expect_lte(length(rpfs_primary), length(rpfs_all))
})
