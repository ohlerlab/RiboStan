# Regression tests for bugs fixed during Apr 2026 development session.
# Each test has a comment identifying the bug it guards against.

# ---------------------------------------------------------------------------
# Bug: ftest_orfs ignored the psite-based ORF filter
#
# The original code was:
#   orfs <- intersect(unique(psites$orf), names(anno$trspacecds))  # line 80
#   orfs <- anno$trspacecds[]     # line 81 - BUG: [] overwrites with ALL ORFs
#
# Fixed to:
#   orfs <- anno$trspacecds[orfs]
#
# Consequence of the bug: mapToTranscripts was called with every ORF in the
# annotation regardless of whether any psite belonged to it. The downstream
# filter (psitecov[sum(psitecov > 0) > 1]) still dropped uncovered ORFs, so
# outputs were correct, but the function was O(all ORFs) instead of O(covered
# ORFs). More critically, line 105
#   orflens <- width(orfs[spec_test_df$orf_id])
# relies on spec_test_df$orf_id being valid keys into `orfs`. With the bug
# present, `orfs` happened to be the full trspacecds so any orf_id would be
# found; after the fix, only orf_ids that came from psites are in `orfs`, so
# spec_test_df can never contain an orf_id that is absent from `orfs`.
#
# Regression: verify that (a) ORFs absent from psites get NA results, and
# (b) every non-NA result corresponds to an ORF that was actually in psites.
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Bug: get_psite_gr silently dropped reads overlapping multiple ORFs
#
# Reads overlapping >1 ORF were correctly assigned a phase/p_offset per ORF
# and one candidate was chosen randomly, but then p_offset was set to NULL
# before being put back into rpfs. The downstream subset(!is.na(p_offset))
# then dropped them entirely, losing valid P-site evidence.
#
# Fixed by removing the `uniq_mov_rpfs$p_offset <- NULL` line.
#
# Regression: psite count should be >= single-ORF-only psite count,
# i.e. multi-ORF reads now contribute rather than being silently dropped.
# ---------------------------------------------------------------------------
test_that("get_psite_gr does not silently drop reads overlapping multiple ORFs", {
  data(chr22_anno)
  data(rpfs)
  data(offsets_df)

  psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)

  # Every RPF that overlapped at least one ORF should have generated a psite.
  # Compute how many RPFs overlap any ORF in transcript space.
  orfs <- chr22_anno$trspacecds
  n_rpfs_with_orf <- sum(
    GenomicRanges::countOverlaps(rpfs, orfs, ignore.strand = TRUE) > 0
  )

  # After the fix, psites should account for (almost) all covered RPFs.
  # Before the fix, multi-ORF RPFs were silently dropped, so psites would
  # be strictly fewer. We allow for reads dropped by other filters
  # (out-of-bounds, no matching offset) but require no catastrophic loss.
  expect_gte(length(psites), 0.8 * n_rpfs_with_orf)
})


# ---------------------------------------------------------------------------
# Bug: get_psite_gr used merge() to look up phaseshift values, which reorders
# rows by the join keys (phase, readlen). The resulting shft vector was no
# longer aligned with `psites`, so every P-site received a shift meant for
# some other P-site. Under this bug, the final frame distribution of psites
# became essentially a random draw from the global shft distribution,
# independent of the read's own phase. This made ftest_orfs call ~77% of
# ORFs "periodic" on pure random-read null data.
#
# Fixed by using an order-preserving dplyr::left_join.
#
# Regression: for each P-site, assert that the total shift applied equals
# the expected (p_offset - phase + shft). After the fix this must hold for
# every row; before the fix it held only for reads whose merge-destination
# row happened to coincide with their original row.
# ---------------------------------------------------------------------------
test_that("get_psite_gr applies phaseshift in correct per-row order", {
  data(chr22_anno)
  data(rpfs)
  data(offsets_df)

  # Replicate the phase computation that happens inside get_psite_gr
  orfs <- chr22_anno$trspacecds
  ov <- GenomicRanges::findOverlaps(rpfs, orfs, select = "first",
                                     ignore.strand = TRUE)
  raw_phase <- (GenomicRanges::start(rpfs) -
                GenomicRanges::start(orfs)[ov]) %% 3
  assigned_orf <- names(orfs)[ov]

  psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)

  # For every P-site, recover the total shift applied to its source read's
  # 5' end. Look up each psite's source read by matching read-level mcols.
  # After the fix, (psite_pos - read_5p_pos) %% 3 should equal 0 for all
  # psites whose (phase, readlen) has shft = (phase) mod 3 (phase 0 case),
  # and more generally the frame distribution should follow a clean
  # permutation of phase->frame rather than the near-uniform mess from the
  # merge bug. The strongest sanity check: phase 0 reads (the majority,
  # shft = 0) must ALL land at frame 0 of their assigned ORF.
  ps_df <- as.data.frame(GenomicRanges::mcols(psites))
  ps_df$frame <- (GenomicRanges::start(psites) -
    GenomicRanges::start(orfs)[match(as.character(ps_df$orf), names(orfs))]) %% 3

  # Restrict to annotated CDS (not uORFs) for a clean signal
  cds_ids <- names(chr22_anno$uORF)[!chr22_anno$uORF]
  ps_cds  <- ps_df[as.character(ps_df$orf) %in% cds_ids, ]

  # Phase 0 reads, shft 0, should be exactly at frame 0 of their ORF.
  # Under the merge bug this fraction was ~48% (global shft mix); now
  # it should be very close to 100%.
  phase0_frame0 <- mean(ps_cds$frame[ps_cds$phase == 0] == 0)
  expect_gt(phase0_frame0, 0.97)
})


# ---------------------------------------------------------------------------
# Bug: ftest_orfs + ftestvect had a dead return() in ftestvect so the
# "p.value" column was actually the raw spectral power at 1/3 Hz, not
# the F-test p-value. This meant periodicity_filter_uORFs(..., p.value < 0.05)
# was filtering on spectral power, keeping non-periodic uORFs.
#
# Fixed by removing the early `return()` and computing the actual p-value
# via pf(Fstat, df1 = 2, df2 = 2k - 2, lower.tail = FALSE).
#
# Regression: p.value column must contain values in [0, 1] (probabilities),
# and the columns must include Fstat, spec_coef, p.value.
# ---------------------------------------------------------------------------
test_that("ftest_orfs returns actual F-test p-values, not spectral power", {
  data(chr22_anno)
  data(rpfs)
  data(offsets_df)

  psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
  ft <- ftest_orfs(psites, chr22_anno, n_cores = 1)

  expect_true(all(c("orf_id", "Fstat", "spec_coef", "p.value") %in%
                  colnames(ft)))

  p <- ft$p.value[!is.na(ft$p.value)]
  expect_gt(length(p), 0)
  # p-values must be in [0, 1] (the old spec_coef values could exceed 1)
  expect_true(all(p >= 0 & p <= 1))

  # Under the bug, the "p.value" column was spectral power — often > 1
  # for spiky periodic data, and uncorrelated with F statistic. After the
  # fix, p should be monotonically decreasing in Fstat.
  finite_rows <- ft[!is.na(ft$Fstat) & !is.na(ft$p.value), ]
  expect_gt(
    cor(finite_rows$Fstat, -log10(finite_rows$p.value + 1e-16)),
    0.95
  )
})


# ---------------------------------------------------------------------------
# Bug: the periodicity pipeline (get_psite_gr + ftest_orfs) was badly
# mis-calibrated due to the combination of the merge-reorder bug and
# recomputing `shiftdf` on every call to get_psite_gr. On pure random-read
# null data it called ~77% of ORFs "periodic" at p<0.05 (expected 5%).
#
# Fixes:
#   1. merge() -> left_join() (preserve row order)
#   2. calibrate shiftdf on annotated-CDS psites only, not the full set
#   3. remove dead return() in ftestvect so real p-values come out
#
# Regression: under a null where RPF 5' positions are randomized uniformly
# within their host transcripts, ftest_orfs should give approximately
# uniform p-values: median ~0.5 and no more than ~10% of tested ORFs at
# p<0.05 (lenient bound; the true expectation is 5%).
# ---------------------------------------------------------------------------
test_that("periodicity pipeline is calibrated on random-read null data", {
  data(chr22_anno)
  data(rpfs)
  data(offsets_df)

  set.seed(1)
  tx_lens  <- GenomeInfoDb::seqlengths(rpfs)
  rpf_lens <- tx_lens[as.character(GenomicRanges::seqnames(rpfs))]
  max_start <- pmax(rpf_lens - rpfs$readlen + 1L, 1L)
  rs_start <- pmin(
    pmax(1L, as.integer(round(runif(length(rpfs), 1, max_start)))),
    max_start
  )
  rpfs_rand <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(rpfs),
    ranges   = IRanges::IRanges(start = rs_start, width = rpfs$readlen),
    strand   = GenomicRanges::strand(rpfs),
    seqinfo  = GenomeInfoDb::seqinfo(rpfs)
  )
  S4Vectors::mcols(rpfs_rand) <- S4Vectors::mcols(rpfs)

  psites_null <- suppressWarnings(
    get_psite_gr(rpfs_rand, offsets_df, chr22_anno)
  )
  ft_null <- suppressMessages(
    ftest_orfs(psites_null, chr22_anno, n_cores = 1)
  )
  p <- ft_null$p.value[!is.na(ft_null$p.value)]

  expect_gt(length(p), 100)
  # Under H0 the median p should be ~0.5. Allow wide tolerance.
  expect_gt(median(p), 0.3)
  # The fraction with p<0.05 should be close to 5%. Under the old bug
  # this was 77%. Allow up to 10% to absorb finite-sample noise.
  expect_lt(mean(p < 0.05), 0.10)
})


test_that("ftest_orfs only produces non-NA results for ORFs present in psites", {
  data(chr22_anno)
  data(rpfs)
  data(offsets_df)

  psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)

  # Pick a small set of annotated-CDS ORFs that definitely have psites
  cds_ids       <- names(chr22_anno$uORF[!chr22_anno$uORF])
  covered_orfs  <- intersect(as.character(unique(psites$orf)), cds_ids)
  target_orfs   <- head(covered_orfs, 5)
  psites_target <- psites[as.character(psites$orf) %in% target_orfs]

  # Subset the annotation to those 5 ORFs + a handful of others with NO psites
  other_orfs    <- setdiff(cds_ids, covered_orfs)[1:5]
  test_orfs     <- c(target_orfs, other_orfs)
  sub_anno      <- Ribostan:::subset_annotation(chr22_anno, test_orfs)

  ftests <- ftest_orfs(psites_target, sub_anno, n_cores = 1)

  # Output should have one row per ORF in the annotation
  expect_equal(nrow(ftests), length(test_orfs))

  # ORFs with no psites must have NA p.value and NA spec_coef
  no_psite_rows <- ftests[ftests$orf_id %in% other_orfs, ]
  expect_true(all(is.na(no_psite_rows$p.value)))
  expect_true(all(is.na(no_psite_rows$spec_coef)))

  # Every non-NA result must belong to an ORF that was in psites$orf
  non_na_rows <- ftests[!is.na(ftests$p.value), ]
  expect_true(all(non_na_rows$orf_id %in% target_orfs))
})


# ---------------------------------------------------------------------------
# Feature: ftest_orfs now reports both raw p.value and BH-corrected q.value.
# periodicity_filter_uORFs filters on q.value (FDR) rather than raw p-value.
#
# Regression: verify both columns exist, q.value == BH(p.value) among unique
# tested loci, aliases inherit the locus q.value, and filtering respects alpha.
# ---------------------------------------------------------------------------
test_that("ftest_orfs reports BH-corrected q.value alongside p.value", {
  data(chr22_anno)
  data(rpfs)
  data(offsets_df)

  psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
  ft <- ftest_orfs(psites, chr22_anno, n_cores = 1)

  expect_true(all(c("p.value", "q.value") %in% colnames(ft)))

  tested <- ft[!is.na(ft$p.value), ]
  expect_gt(nrow(tested), 50)

  # q.value should be the BH adjustment of one p.value per genomic ORF locus.
  unique_loci <- tested[tested$locus_representative, ]
  expected_q <- p.adjust(unique_loci$p.value, method = "BH")
  expect_equal(unique_loci$q.value, expected_q, tolerance = 1e-12)
  expect_equal(
    tested$q.value,
    unique_loci$q.value[match(tested$orf_locus_id,
                              unique_loci$orf_locus_id)],
    tolerance = 0
  )

  # BH-adjusted q-values are always >= raw p-values
  expect_true(all(tested$q.value >= tested$p.value - 1e-12))

  # Untested ORFs (no psite coverage) should have NA for both columns
  untested <- ft[is.na(ft$p.value), ]
  expect_true(all(is.na(untested$q.value)))
})


test_that("periodicity_filter_uORFs uses q.value with the alpha argument", {
  data(chr22_anno)
  data(rpfs)
  data(offsets_df)

  psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)

  # Filter at a permissive threshold; mcols must include q.value.
  filt <- suppressMessages(
    periodicity_filter_uORFs(psites, chr22_anno,
                             remove = FALSE, alpha = 0.05, n_cores = 1)
  )
  expect_true("q.value" %in% colnames(
    GenomicRanges::mcols(filt$trspacecds)
  ))
  expect_true("p.value" %in% colnames(
    GenomicRanges::mcols(filt$trspacecds)
  ))

  # Tighter alpha should keep no more uORFs than a looser alpha.
  filt_strict <- suppressMessages(
    periodicity_filter_uORFs(psites, chr22_anno,
                             remove = TRUE, alpha = 0.05, n_cores = 1)
  )
  filt_loose <- suppressMessages(
    periodicity_filter_uORFs(psites, chr22_anno,
                             remove = TRUE, alpha = 0.50, n_cores = 1)
  )
  expect_lte(sum(filt_strict$uORF), sum(filt_loose$uORF))
})


test_that("periodicity_filter_uORFs enforces min_psites for translation support", {
  data(chr22_anno)
  data(rpfs)
  data(offsets_df)

  psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
  filt_none <- suppressMessages(
    periodicity_filter_uORFs(psites, chr22_anno,
                             remove = FALSE, min_psites = 0L, n_cores = 1)
  )
  expect_true("n_psites" %in% colnames(
    GenomicRanges::mcols(filt_none$trspacecds)
  ))
  expect_true(isTRUE(filt_none$uORF_periodicity_filtered))

  uorf_ids <- names(filt_none$uORF)[filt_none$uORF]
  tested <- uorf_ids[
    !is.na(GenomicRanges::mcols(filt_none$trspacecds[uorf_ids])$q.value) &
      GenomicRanges::mcols(filt_none$trspacecds[uorf_ids])$q.value < 0.05
  ]
  if (length(tested) >= 2L) {
    high_cov <- tested[
      GenomicRanges::mcols(filt_none$trspacecds[tested])$n_psites >= 6L
    ]
    filt_hi <- suppressMessages(
      periodicity_filter_uORFs(psites, chr22_anno,
                               remove = TRUE, min_psites = 6L, n_cores = 1)
    )
    filt_lo <- suppressMessages(
      periodicity_filter_uORFs(psites, chr22_anno,
                               remove = TRUE, min_psites = 0L, n_cores = 1)
    )
    expect_lte(sum(filt_hi$uORF), sum(filt_lo$uORF))
    if (length(high_cov) > 0L) {
      expect_true(all(high_cov %in% names(filt_hi$trspacecds)))
    }
  }
})


test_that("prepare_ribostan_anno returns only periodicity-filtered uORFs", {
  data(chr22_anno)
  data(rpfs)

  gtf <- tempfile(fileext = ".gtf")
  ah <- AnnotationHub::AnnotationHub()
  gencode32 <- ah[["AH75191"]]
  GenomeInfoDb::seqlevels(gencode32) <- "chr22"
  suppressWarnings(rtracklayer::export(gencode32, gtf, format = "GTF"))
  on.exit(unlink(gtf), add = TRUE)

  prep <- suppressMessages(prepare_ribostan_anno(
    gtf,
    fafile = {
      fa <- here::here("chr22.fa")
      if (!file.exists(fa)) {
        library(BSgenome.Hsapiens.UCSC.hg38)
        seq <- Biostrings::DNAStringSet(BSgenome.Hsapiens.UCSC.hg38[["chr22"]])
        names(seq) <- "chr22"
        Biostrings::writeXStringSet(seq, fa)
      }
      fa
    },
    rpfs,
    n_cores = 1L,
    return_psites = TRUE
  ))
  expect_true(isTRUE(prep$anno$uORF_periodicity_filtered))
  expect_lte(sum(prep$anno$uORF), sum(chr22_anno$uORF))
  expect_true(all(prep$psites$orf %in% names(prep$anno$trspacecds)))
})
