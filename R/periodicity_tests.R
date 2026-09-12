#' @importFrom multitaper dpss spec.mtm dropFreqs
NULL

#' Test a numeric vector for periodicity using the multitaper package
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param psit A list of SimpleRleLists - RPF coverage split by readlength
#' @param k - The number of slepians to use for the test
#' @param bw - The bandwidth to use for the test
#'
#' @details This function pads short vectors up to 50, and performs a multitaper test for 3bp
#' peridicity on it's input, returning a spectral coefficient and p-value
#' @return a numeric vector with the spectral coefficient at 0.333... and the pvalue for the test

ftestvect <- function(psit, k = 24, bw = 12) {
  psit <- as.numeric(psit)
  slepians_values <- dpss(
    n = length(psit) %>% ifelse(. < 25, 50, .),
    k = k, nw = bw
  )

  if (length(psit) < 25) {
    remain <- 50 - length(psit)
    halfrmn <- as.integer(remain / 2)
    psit <- c(rep(0, halfrmn), psit, rep(0, remain %% 2 + halfrmn))
  }
  padding <- if (length(psit) < 1024 / 2) 1024 else "default"

  resSpec1 <- spec.mtm(as.ts(psit),
    k = k, nw = bw, nFFT = padding,
    centreWithSlepians = TRUE, Ftest = TRUE,
    maxAdaptiveIterations = 100, returnZeroFreq = FALSE,
    plot = FALSE, dpssIN = slepians_values
  )

  closestfreqind <- which.min(abs(resSpec1$freq - (1 / 3)))
  Fstat_3nt <- resSpec1$mtm$Ftest[closestfreqind]
  spect_3nt <- resSpec1$spec[closestfreqind]
  # spec.mtm's Ftest has F(2, 2k - 2) distribution under H0 of no
  # harmonic component.
  pval <- pf(q = Fstat_3nt, df1 = 2, df2 = 2 * k - 2, lower.tail = FALSE)
  c(Fstat = Fstat_3nt, spec_coef = spect_3nt, p.value = pval)
}

.orf_locus_key <- function(cdsgrl) {
  if (is(cdsgrl, "GRanges")) {
    gr_names <- names(cdsgrl)
    cdsgrl <- GenomicRanges::GRangesList(
      lapply(seq_along(cdsgrl), function(i) cdsgrl[i])
    )
    names(cdsgrl) <- gr_names
  }
  keys <- vapply(seq_along(cdsgrl), function(i) {
    gr <- cdsgrl[[i]]
    strand_value <- as.character(GenomicRanges::strand(gr))[1L]
    order_index <- if (strand_value == "-") {
      order(GenomicRanges::start(gr), decreasing = TRUE)
    } else {
      order(GenomicRanges::start(gr))
    }
    gr <- gr[order_index]
    paste(
      paste0(
        as.character(GenomicRanges::seqnames(gr)), ":",
        GenomicRanges::start(gr), "-", GenomicRanges::end(gr), ":",
        as.character(GenomicRanges::strand(gr))
      ),
      collapse = ";"
    )
  }, character(1))
  names(keys) <- names(cdsgrl)
  keys
}


#' Run Multitaper Tests on All ORFs
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param psites A GRanges object containing psites
#' @param anno annotation object
#' @param n_cores number of cores to use
#'
#' @details Candidates with exactly the same ordered genomic exon blocks and
#' strand are one ORF locus. P-site records fanned out across aliases are
#' deduplicated by read name and locus-relative position, the locus is tested
#' once, BH correction is applied across unique tested loci, and the resulting
#' statistics are propagated to every transcript alias.
#' @return A data frame with one row per ORF annotation. In addition to the
#' periodicity statistics, it includes \code{orf_locus_id},
#' \code{locus_representative}, \code{n_locus_aliases}, and
#' \code{compatible_orf_ids} so transcript membership is preserved while locus
#' calls can be counted once.
#' @export
#' @examples
#' data(chr22_anno)
#' data(rpfs)
#' data(offsets_df)
#' psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
#' ftests <- ftest_orfs(psites %>% head(10000), chr22_anno, 
#'   n_cores=1)

ftest_orfs <- function(psites, anno, n_cores=1) {
  all_orf_ids <- names(anno$cdsgrl)
  locus_by_orf <- .orf_locus_key(anno$cdsgrl)
  aliases_by_locus <- split(all_orf_ids, locus_by_orf[all_orf_ids])
  representative_by_locus <- vapply(
    aliases_by_locus, `[[`, character(1), 1L
  )

  covered_orf_ids <- intersect(
    as.character(unique(psites$orf)),
    names(anno$trspacecds)
  )
  orfs <- anno$trspacecds[covered_orf_ids]
  mapped <- GenomicFeatures::mapToTranscripts(
    psites, orfs, ignore.strand = TRUE
  )
  mapped <- mapped[
    psites$orf[mapped$xHits] == names(orfs)[mapped$transcriptsHits]
  ]

  mapped_orf <- names(orfs)[mapped$transcriptsHits]
  mapped_locus <- unname(locus_by_orf[mapped_orf])
  read_id <- names(psites)[mapped$xHits]
  if (is.null(read_id)) {
    read_id <- paste0("psite_", mapped$xHits)
  } else {
    missing_id <- is.na(read_id) | !nzchar(read_id)
    read_id[missing_id] <- paste0("psite_", mapped$xHits[missing_id])
  }

  # The same physical P-site record fanned out to several transcript aliases
  # contributes once to their shared genomic ORF locus. Keep distinct
  # positions for genuinely multimapping reads.
  mapped_position <- GenomicRanges::start(mapped)
  read_locus_key <- paste(
    mapped_locus, read_id, mapped_position, sep = "\r"
  )
  keep <- !duplicated(read_locus_key)
  mapped_locus <- mapped_locus[keep]
  mapped_position <- mapped_position[keep]

  position_by_locus <- split(mapped_position, mapped_locus)
  tested_loci <- intersect(names(position_by_locus),
                            names(representative_by_locus))
  psitecov <- lapply(tested_loci, function(locus_id) {
    representative <- representative_by_locus[[locus_id]]
    as.numeric(tabulate(
      position_by_locus[[locus_id]],
      nbins = GenomicRanges::width(anno$trspacecds[representative])
    ))
  })
  names(psitecov) <- tested_loci
  psitecov <- psitecov[
    vapply(psitecov, function(v) sum(v > 0) > 1L, logical(1))
  ]

  # Fraction of each candidate's P-sites that fall in its own reading frame
  # (relative positions 1, 4, 7, ... of the window, since mapToTranscripts
  # puts each ORF's own start at position 1). Mirrors RiboTaper's
  # Phase_P_sites_frame / P_sites_sum gate (CCDS_orf_finder.R /
  # NONCCDS_orf_finder.R): a candidate is only "periodic" there if this
  # exceeds 0.5, since the multitaper F-test alone detects 3-nt
  # periodicity at any phase, not specifically the candidate's own frame.
  pct_inframe <- vapply(psitecov, function(v) {
    v <- as.numeric(v)
    total <- sum(v)
    if (total == 0) return(NA_real_)
    sum(v[seq(1, length(v), by = 3)]) / total
  }, numeric(1))
  # now run multitaper tests on our data using multiple
  # cores if available
  message("running multitaper tests, this will be slow for a full dataset...")
  spec_tests <- psitecov %>% 
    parallel::mclapply(F = ftestvect, mc.cores = n_cores)
  if (length(spec_tests) == 0L) {
    locus_test_df <- data.frame(
      orf_locus_id = character(),
      Fstat = numeric(),
      spec_coef = numeric(),
      p.value = numeric(),
      q.value = numeric(),
      n_psites = integer(),
      pct_inframe = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    locus_test_df <- as.data.frame(t(simplify2array(spec_tests)))
    locus_test_df$orf_locus_id <- rownames(locus_test_df)
    rownames(locus_test_df) <- NULL
    locus_test_df$pct_inframe <-
      unname(pct_inframe[locus_test_df$orf_locus_id])
    locus_test_df$n_psites <- as.integer(vapply(
      psitecov[locus_test_df$orf_locus_id], sum, numeric(1)
    ))
    representative_ids <- unname(
      representative_by_locus[locus_test_df$orf_locus_id]
    )
    locus_widths <- GenomicRanges::width(
      anno$trspacecds[representative_ids]
    )
    locus_test_df$spec_coef <-
      sqrt(locus_test_df$spec_coef) / locus_widths
    # Correct across unique genomic ORF loci, not transcript aliases.
    locus_test_df$q.value <- p.adjust(
      locus_test_df$p.value, method = "BH"
    )
    locus_test_df <- locus_test_df[
      , c("orf_locus_id", "Fstat", "spec_coef", "p.value", "q.value",
          "n_psites", "pct_inframe")
    ]
  }

  compatible_orf_ids <- vapply(
    aliases_by_locus, paste, character(1), collapse = ","
  )
  alias_df <- data.frame(
    orf_id = all_orf_ids,
    orf_locus_id = unname(locus_by_orf[all_orf_ids]),
    stringsAsFactors = FALSE
  )
  alias_df$locus_representative <- alias_df$orf_id ==
    unname(representative_by_locus[alias_df$orf_locus_id])
  alias_df$n_locus_aliases <- as.integer(
    lengths(aliases_by_locus)[alias_df$orf_locus_id]
  )
  alias_df$compatible_orf_ids <- unname(
    compatible_orf_ids[alias_df$orf_locus_id]
  )
  left_join(alias_df, locus_test_df, by = "orf_locus_id")
}


#' Select one uORF per nested (shared-stop) group by P-site read density
#'
#' Alternative start codons upstream of a shared stop (e.g. two in-frame
#' ATGs reading into the same stop) generate uORF candidates that are
#' nested/redundant rather than independent: their P-site coverage is the
#' same signal seen through different windows. Testing all of them for
#' periodicity would inflate the multiple-testing burden with correlated
#' hits, so exactly one candidate per \code{(transcript, stop)} group is
#' kept - the one with the highest P-site read density (P-sites per
#' nucleotide of ORF length).
#'
#' This uses raw P-site counts rather than the multitaper spectral
#' coefficient, so the redundant candidates can be dropped \emph{before}
#' the (slow) periodicity test is run on them at all, rather than after.
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param psites GRanges of psites (output of \code{get_psite_gr}), with an
#'   \code{orf} mcol.
#' @param anno Annotation object with \code{trspacecds}.
#' @param orf_ids Character vector of candidate uORF ids to prune among.
#'   ORFs that do not share a \code{(transcript, stop)} group with any
#'   other id in \code{orf_ids} are always kept.
#' @return \code{orf_ids}, restricted to the highest-read-density member of
#'   each nested group. Order is not preserved.
#' @export
#' @examples
#' data(chr22_anno)
#' data(rpfs)
#' data(offsets_df)
#' psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
#' uorf_ids <- names(chr22_anno$uORF)[chr22_anno$uORF]
#' survivors <- prune_nested_uorfs(psites, chr22_anno, uorf_ids)
prune_nested_uorfs <- function(psites, anno, orf_ids) {
  orf_ids <- unique(orf_ids)
  if (length(orf_ids) == 0L) return(orf_ids)

  orfs <- anno$trspacecds[orf_ids]
  tx  <- sub("_[0-9]+$", "", orf_ids)
  stp <- GenomicRanges::end(orfs)
  grp <- paste0(tx, "@", stp)

  read_counts <- tabulate(
    match(as.character(psites$orf), orf_ids),
    nbins = length(orf_ids)
  )
  density <- read_counts / width(orfs)

  # Within each (tx, stop) group, keep exactly one winner: the ORF with
  # the highest read density. Ties (e.g. no P-site coverage at all) are
  # broken by the first occurrence so pruning is always decisive.
  idx_in_grp <- seq_along(orf_ids)
  winner_idx <- as.integer(stats::ave(
    idx_in_grp, grp,
    FUN = function(i) {
      v <- density[i]
      if (all(v == 0)) i[1] else i[which.max(v)]
    }
  ))
  orf_ids[idx_in_grp == winner_idx]
}


#' Select one uORF per nested (shared-stop) group by RiboTaper's own rule
#'
#' RiboTaper (Calviello et al. 2016) resolves multiple in-frame ATGs sharing
#' a stop codon differently to \code{\link{prune_nested_uorfs}}: every
#' candidate is periodicity-tested first (raw p-value, no FDR correction),
#' then, among the candidates significant at \code{alpha}, it walks from the
#' most upstream ATG and picks the first one for which the gap to its
#' closest downstream neighbour ATG has more than 5 P-sites, over 50% of
#' them in-frame - i.e. there is direct evidence that translation actually
#' initiates that far upstream, rather than leaky scanning into the next
#' ATG. If no candidate's gap clears that bar (or fewer than two candidates
#' are significant), the single candidate with the lowest raw p-value in
#' the group is kept instead. This matches the rule and its exact
#' thresholds (\code{P_sites_next_start > 5}, \code{pct_inframe > 0.5},
#' \code{p < 0.05} for group membership) in RiboTaper's
#' \code{NONCCDS_orf_finder.R} / \code{CCDS_orf_finder.R} (\code{max_sit}
#' selection, "max_P_sites" method).
#'
#' Unlike \code{\link{prune_nested_uorfs}}, this requires every candidate in
#' \code{orf_ids} to already have a periodicity p-value (i.e. it must be
#' called \emph{after} \code{\link{ftest_orfs}} has tested all of them,
#' not just a pre-selected representative) - RiboTaper's rule is inherently
#' a post-hoc choice among already-tested candidates.
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param psites GRanges of psites (output of \code{get_psite_gr}). Used
#'   here only for raw P-site positions (deduplicated so a single physical
#'   P-site attributed to several overlapping ORFs by
#'   \code{\link{get_psite_gr}} is not double-counted); the \code{orf} mcol
#'   is not used, since the "gap" window is defined purely by transcript
#'   position.
#' @param anno Annotation object with \code{trspacecds}.
#' @param orf_ids Character vector of candidate uORF ids to select among.
#' @param tested_df Data frame with columns \code{orf_id} and \code{p.value}
#'   (raw, uncorrected), covering every id in \code{orf_ids} - typically
#'   \code{ftest_orfs(psites, subset_annotation(anno, orf_ids))}. If it also
#'   has a \code{pct_inframe} column (as \code{\link{ftest_orfs}} now
#'   returns), that is folded into group-membership eligibility too - see
#'   \code{min_frame_frac}.
#' @param alpha Raw p-value threshold used both for group-membership
#'   eligibility and as RiboTaper's own hardcoded significance cutoff
#'   (default \code{0.05}, matching RiboTaper).
#' @param min_frame_frac Minimum fraction of a candidate's P-sites that must
#'   fall in its own reading frame for its p-value to count at all (default
#'   \code{0}, i.e. no extra requirement). Candidates at or below this
#'   fraction are treated as untested (p-value nulled out) for both the
#'   significance walk and the lowest-p-value fallback - matching
#'   RiboTaper, which never assigns \code{ORF_pval_multi_ribo} to such a
#'   candidate in the first place (see \code{\link{ftest_orfs}}'s
#'   \code{pct_inframe}). Set \code{min_frame_frac = 0.5} to match
#'   RiboTaper exactly. Ignored if \code{tested_df} has no
#'   \code{pct_inframe} column.
#' @param min_psites Minimum total P-sites required for a candidate to
#'   participate in nested selection. Candidates below this threshold are
#'   excluded before choosing a shared-stop winner. Default \code{0L}.
#' @return \code{orf_ids}, restricted to one representative per
#'   \code{(transcript, stop)} group. Order is not preserved.
#' @export
#' @examples
#' data(chr22_anno)
#' data(rpfs)
#' data(offsets_df)
#' psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
#' uorf_ids <- names(chr22_anno$uORF)[chr22_anno$uORF]
#' tested <- suppressMessages(ftest_orfs(psites,
#'   subset_annotation(chr22_anno, uorf_ids), n_cores = 1))
#' survivors <- select_ribotaper_nested_uorfs(psites, chr22_anno, uorf_ids,
#'   tested)
select_ribotaper_nested_uorfs <- function(psites, anno, orf_ids, tested_df,
                                          alpha = 0.05, min_frame_frac = 0,
                                          min_psites = 0L) {
  orf_ids <- unique(orf_ids)
  if (length(orf_ids) == 0L) return(orf_ids)

  orfs <- anno$trspacecds[orf_ids]
  tx  <- as.character(GenomicRanges::seqnames(orfs))
  stp <- GenomicRanges::end(orfs)
  st  <- GenomicRanges::start(orfs)
  grp <- paste0(tx, "@", stp)
  pval <- tested_df$p.value[match(orf_ids, tested_df$orf_id)]
  frame_frac <- if ("pct_inframe" %in% names(tested_df)) {
    tested_df$pct_inframe[match(orf_ids, tested_df$orf_id)]
  } else {
    rep(NA_real_, length(orf_ids))
  }

  # Raw, de-duplicated psites: one row per physical (read, position) pair
  # regardless of how many ORFs get_psite_gr() attributed it to, so gap
  # counts reflect actual P-site density rather than ORF-overlap fan-out.
  psite_key <- paste(as.character(GenomicRanges::seqnames(psites)),
                     GenomicRanges::start(psites), names(psites))
  raw_psites <- psites[!duplicated(psite_key)]
  pos_by_tx <- split(GenomicRanges::start(raw_psites),
                     as.character(GenomicRanges::seqnames(raw_psites)))

  # Candidates whose own P-sites are majority out-of-frame never get a
  # usable p-value in RiboTaper (its ORF_pval_multi_ribo stays NA and the
  # row is dropped before stop-groups are even formed), so treat them the
  # same way here: null the p-value out entirely, rather than only
  # excluding them from the "significant" subset. That way they are
  # skipped by both the significance walk *and* the best_overall fallback
  # below (`which.min` ignores NAs), matching RiboTaper's actual ordering
  # of operations.
  pval[!is.na(frame_frac) & frame_frac <= min_frame_frac] <- NA_real_
  if ("n_psites" %in% names(tested_df)) {
    n_psites <- tested_df$n_psites[match(orf_ids, tested_df$orf_id)]
    pval[!is.na(n_psites) & n_psites < min_psites] <- NA_real_
  }

  df <- data.frame(orf_id = orf_ids, tx = tx, start = st, grp = grp,
                   p.value = pval, stringsAsFactors = FALSE)

  grp_sizes <- table(df$grp)
  singleton_grps <- names(grp_sizes)[grp_sizes == 1L]
  winners <- df$orf_id[df$grp %in% singleton_grps]

  multi_df <- df[!df$grp %in% singleton_grps, , drop = FALSE]
  if (nrow(multi_df) == 0L) return(winners)

  for (g in split(multi_df, multi_df$grp)) {
    sig <- g[!is.na(g$p.value) & g$p.value < alpha, , drop = FALSE]
    best_overall <- g$orf_id[which.min(g$p.value)]
    if (length(best_overall) == 0L) best_overall <- g$orf_id[1L]
    if (nrow(sig) < 2L) {
      winners <- c(winners, best_overall)
      next
    }
    sig <- sig[order(sig$start), ]
    winner <- NA_character_
    tx_pos <- pos_by_tx[[sig$tx[1L]]]
    for (i in seq_len(nrow(sig) - 1L)) {
      gap_start <- sig$start[i]
      gap_end   <- sig$start[i + 1L] - 1L
      if (gap_end < gap_start || is.null(tx_pos)) next
      in_gap <- tx_pos[tx_pos >= gap_start & tx_pos <= gap_end]
      n_gap <- length(in_gap)
      if (n_gap > 5L) {
        pct_inframe <- sum((in_gap - gap_start) %% 3L == 0L) / n_gap
        if (pct_inframe > 0.5) {
          winner <- sig$orf_id[i]
          break
        }
      }
    }
    if (is.na(winner)) winner <- best_overall
    winners <- c(winners, winner)
  }
  winners
}


#' Run Multitaper tests on a set of ORFs
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param psites GRanges object with psite information
#' @param anno  annotation object
#' @param remove whether to remove uORFs that fail the translation-support
#'   criteria (default \code{TRUE}). When \code{FALSE}, all sequence candidates
#'   are kept and periodicity statistics are attached as mcols only.
#' @param alpha significance threshold; applied to \code{q.value} when
#'   \code{fdr_correct = TRUE} (default), or to the raw \code{p.value} when
#'   \code{fdr_correct = FALSE}. Default 0.05.
#' @param min_psites minimum number of P-sites overlapping a uORF for it to
#'   count as translation-supported (default \code{6}, matching RiboTaper's
#'   \code{P_sites_sum > 5} filter). Set to \code{0L} to disable.
#' @param n_cores number of cores to use
#' @param nested_rule how to pick one representative per nested (shared-stop)
#'   group of uORF candidates:
#'   \describe{
#'     \item{\code{"density"}}{\code{\link{prune_nested_uorfs}}
#'       picks the highest-P-site-density candidate \emph{before} the
#'       (slow) F-test runs, so only one candidate per group is ever
#'       tested. Cheaper, but not what RiboTaper does.}
#'     \item{\code{"ribotaper"} (default)}{every distinct candidate locus in every group
#'       is F-tested first, then \code{\link{select_ribotaper_nested_uorfs}} applies
#'       RiboTaper's own upstream-ATG-with-signal rule post-hoc. Matches
#'       RiboTaper exactly but is more expensive (tests every distinct ATG
#'       locus, not one representative per group) and, because
#'       \code{\link{ftest_orfs}}'s BH correction runs across that larger tested
#'       set, is normally
#'       paired with \code{fdr_correct = FALSE} (RiboTaper itself does not
#'       FDR-correct).}
#'   }
#' @param fdr_correct Whether the `remove = TRUE` / significance decision
#'   uses the Benjamini-Hochberg \code{q.value} (default \code{TRUE}) or the
#'   raw, uncorrected \code{p.value} (\code{FALSE}). RiboTaper's own
#'   published threshold is a raw \code{p < 0.05} per candidate with no
#'   multiple-testing correction; Ribostan defaults to BH-corrected q-values
#'   because it tests many more candidates genome-wide. Both p.value and
#'   q.value are always attached as mcols regardless of this setting.
#' @param min_frame_frac Minimum fraction of a uORF's P-sites that must fall
#'   in its own reading frame (default \code{0}, i.e. disabled/backward
#'   compatible). RiboTaper requires this fraction to exceed \code{0.5}
#'   before a candidate's periodicity p-value counts at all (see
#'   \code{\link{ftest_orfs}}'s \code{pct_inframe}) - the multitaper F-test
#'   alone detects 3-nt periodicity regardless of which of the 3 possible
#'   phases carries it, so without this gate a candidate can be "periodic"
#'   while most of its signal is actually in a different frame than its own
#'   start codon implies. Set to \code{0.5} to match RiboTaper exactly; also
#'   passed to \code{\link{select_ribotaper_nested_uorfs}} when
#'   \code{nested_rule = "ribotaper"} so group-membership eligibility uses
#'   the same gate.
#'
#' @details
#' Ribostan uORF calls require periodic RPF signal, analogous to RiboTaper:
#' sequence candidates from \code{\link{load_annotation}} are not treated as
#' bona fide uORFs until this function (or \code{\link{prepare_ribostan_anno}})
#' is run on P-sites derived from ribo-seq data.
#'
#' Nested/redundant uORF candidates that share a stop position with another
#' candidate on the same transcript are collapsed to one representative per
#' group, by \code{nested_rule}. The F-test is then applied to every
#' surviving uORF ("density") or to every raw candidate ("ribotaper"). Both
#' the raw p-value and the Benjamini-Hochberg q-value (FDR-corrected across
#' the unique tested genomic ORF loci) are propagated to every compatible
#' transcript alias and attached to `anno$trspacecds` as mcols, along with
#' \code{n_psites}, \code{pct_inframe}, \code{orf_locus_id},
#' \code{locus_representative}, \code{n_locus_aliases}, and
#' \code{compatible_orf_ids}. When `remove = TRUE`, uORFs that
#' fail the significance test (per \code{fdr_correct}), fall below
#' \code{min_psites} or \code{min_frame_frac}, or were pruned/untested are
#' dropped from the annotation.
#' @return an annotation object with
#'   locus-level Fstat/spec_coef/p.value/q.value/n_psites/pct_inframe and locus
#'   identity metadata attached per transcript alias,
#'   \code{uORF_periodicity_filtered = TRUE}, and if `remove = TRUE` only
#'   translation-supported uORF aliases retained.
#' @export
#' @examples
#' data(chr22_anno)
#' data(rpfs)
#' data(offsets_df)
#' psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
#' filteredanno <- periodicity_filter_uORFs(psites, chr22_anno)

periodicity_filter_uORFs <- function(psites, anno, remove = TRUE,
                                     alpha = 0.05, min_psites = 6L,
                                     n_cores = 1,
                                     nested_rule = c("ribotaper", "density"),
                                     fdr_correct = TRUE,
                                     min_frame_frac = 0) {
  stopifnot(!is.null(anno$uORF))
  stopifnot(is.numeric(min_psites), length(min_psites) == 1L, min_psites >= 0)
  nested_rule <- match.arg(nested_rule)
  uORFs <- anno$uORF
  uORFs <- unique(names(uORFs[uORFs]))
  locus_by_uorf <- .orf_locus_key(anno$cdsgrl[uORFs])
  aliases_by_locus <- split(uORFs, locus_by_uorf[uORFs])
  representative_by_locus <- vapply(
    aliases_by_locus, `[[`, character(1), 1L
  )
  compatible_orf_ids <- vapply(
    aliases_by_locus, paste, character(1), collapse = ","
  )
  uorf_psite_counts <- stats::setNames(
    tabulate(match(as.character(psites$orf), uORFs), nbins = length(uORFs)),
    uORFs
  )
  if (nested_rule == "density") {
    survivors <- prune_nested_uorfs(psites, anno, uORFs)
    ftestdf <- ftest_orfs(psites, subset_annotation(anno, survivors),
                          n_cores = n_cores)
  } else {
    # "ribotaper": every raw candidate must be tested before a
    # representative can be chosen, since the selection rule itself
    # depends on each candidate's own p-value.
    ftestdf_all <- ftest_orfs(psites, subset_annotation(anno, uORFs),
                              n_cores = n_cores)
    survivors <- select_ribotaper_nested_uorfs(psites, anno, uORFs,
                                               ftestdf_all, alpha = alpha,
                                               min_frame_frac = min_frame_frac,
                                               min_psites = min_psites)
    ftestdf <- ftestdf_all[ftestdf_all$orf_id %in% survivors, ]
  }
  # Initialise per-ORF metadata columns if missing. ftest_orfs now returns
  # Fstat, spec_coef, p.value, q.value, pct_inframe; historically it
  # returned only spec_coef and p.value.
  for (col in c("Fstat", "spec_coef", "p.value", "q.value", "pct_inframe")) {
    if (is.null(mcols(anno$trspacecds)[[col]]))
      mcols(anno$trspacecds)[[col]] <- NA_real_
  }
  if (is.null(mcols(anno$trspacecds)$n_psites))
    mcols(anno$trspacecds)$n_psites <- NA_integer_
  if (is.null(mcols(anno$trspacecds)$orf_locus_id))
    mcols(anno$trspacecds)$orf_locus_id <- NA_character_
  if (is.null(mcols(anno$trspacecds)$locus_representative))
    mcols(anno$trspacecds)$locus_representative <- FALSE
  if (is.null(mcols(anno$trspacecds)$n_locus_aliases))
    mcols(anno$trspacecds)$n_locus_aliases <- NA_integer_
  if (is.null(mcols(anno$trspacecds)$compatible_orf_ids))
    mcols(anno$trspacecds)$compatible_orf_ids <- NA_character_
  uorf_idx <- match(uORFs, names(anno$trspacecds))
  uorf_loci <- unname(locus_by_uorf[uORFs])
  mcols(anno$trspacecds)$orf_locus_id[uorf_idx] <- uorf_loci
  mcols(anno$trspacecds)$locus_representative[uorf_idx] <-
    uORFs == unname(representative_by_locus[uorf_loci])
  mcols(anno$trspacecds)$n_locus_aliases[uorf_idx] <- as.integer(
    lengths(aliases_by_locus)[uorf_loci]
  )
  mcols(anno$trspacecds)$compatible_orf_ids[uorf_idx] <- unname(
    compatible_orf_ids[uorf_loci]
  )
  n_psite_idx <- match(names(uorf_psite_counts), names(anno$trspacecds))
  mcols(anno$trspacecds)$n_psites[n_psite_idx] <- as.integer(uorf_psite_counts)
  ftest_idx <- match(ftestdf$orf_id, names(anno$trspacecds))
  for (col in c("Fstat", "spec_coef", "p.value", "q.value", "pct_inframe")) {
    mcols(anno$trspacecds)[[col]][ftest_idx] <- ftestdf[[col]]
  }
  mcols(anno$trspacecds)$n_psites[ftest_idx] <- ftestdf$n_psites
  if (remove) {
    sig_col <- if (fdr_correct) "q.value" else "p.value"
    periodic_uORFs <- ftestdf %>%
      filter(!is.na(.data[[sig_col]]) & .data[[sig_col]] < alpha &
               (is.na(.data$pct_inframe) | .data$pct_inframe > min_frame_frac)) %>%
      .$orf_id
    if (min_psites > 0L) {
      tested_psite_counts <- stats::setNames(
        ftestdf$n_psites, ftestdf$orf_id
      )
      periodic_uORFs <- periodic_uORFs[
        tested_psite_counts[periodic_uORFs] >= min_psites
      ]
    }
    non_periodic_uORFs <- setdiff(uORFs, periodic_uORFs)
    orfs_to_keep <- setdiff(names(anno$trspacecds), non_periodic_uORFs)
    anno <- subset_annotation(anno, orfs_to_keep)
  }
  anno$uORF_periodicity_filtered <- TRUE
  anno
}


#' Flag uORFs whose read coverage is dominated by multi-mapping reads
#'
#' Reproduces RiboTaper's own post-hoc quality-control step (Calviello
#' et al. 2016): after periodicity filtering, an ORF is excluded if (a) too
#' much of its coverage comes only from multi-mapping reads, rather than
#' from uniquely-mapping ones, or (b) its overall coverage breadth is too
#' sparse to trust. RiboTaper applies both conditions to its final
#' periodic-and-non-CDS-overlapping candidate set (see
#' \code{CCDS_orf_finder.R} / \code{NONCCDS_orf_finder.R}):
#' \code{all_sORFs_CCDS_periodic_nocds_filtered_multi <- }
#' \code{  all_sORFs_CCDS_periodic_nocds[(pct_covered_onlymulti_ribo /}
#' \code{  pct_region_covered_ribo) < 0.3, ]}, then further restricted to
#' \code{pct_region_covered_ribo > 0.3}.
#'
#' Critically, RiboTaper's \code{pct_region_covered_ribo} /
#' \code{pct_covered_onlymulti_ribo} come from \code{bedtools coverage -s
#' -split} run on the \emph{full aligned read footprint} in
#' \code{RIBO_best.bam} (all alignments) versus \code{RIBO_unique.bam}
#' (\code{MAPQ >= 50} only) - see \code{analyze_multi_clust.bash} /
#' \code{include_multi_nomerge.R}, whose output columns are literally named
#' \code{pct_region_covered_ribo} and \code{pct_covered_onlymulti_ribo}.
#' This is full-read-footprint coverage breadth, \emph{not} P-site coverage:
#' a single RPF spans ~28-30nt, so it marks many nucleotides covered at
#' once, unlike a P-site call which marks exactly one. Passing P-site-only
#' GRanges here would systematically understate coverage and make
#' \code{min_coverage} far stricter than RiboTaper's own threshold (an
#' earlier version of this function did exactly that).
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param reads GRanges of aligned read footprints spanning their full
#'   width (e.g. \code{\link{get_readgr}} called \emph{without}
#'   \code{offsets_df}, so reads are not collapsed to single-nucleotide
#'   P-sites), with an \code{NH} mcol. Must be in the same coordinate space
#'   as \code{anno$trspacecds} (transcript space). Call with reads read
#'   using \code{multimap != "unique"} (otherwise no multi-mapping reads
#'   survive to evaluate).
#' @param anno Annotation object with \code{trspacecds}.
#' @param orf_ids Character vector of ORF ids to evaluate.
#' @param threshold Numeric in [0, 1]; ORFs at or above this
#'   multi-mapping-only coverage fraction are flagged (default \code{0.3},
#'   matching RiboTaper).
#' @param min_coverage Numeric in [0, 1]; ORFs at or below this overall
#'   coverage breadth (\code{pct_covered_all}) are flagged (default
#'   \code{0}, i.e. disabled - only the zero-coverage edge case is caught).
#'   Set \code{min_coverage = 0.3} to match RiboTaper's own
#'   \code{pct_region_covered_ribo > 0.3} requirement exactly.
#' @return A data frame with one row per member of \code{orf_ids} and
#'   columns \code{orf_id}, \code{pct_covered_all}, \code{pct_covered_unique},
#'   \code{pct_covered_onlymulti}, \code{low_coverage} (logical;
#'   \code{pct_covered_all <= min_coverage}), \code{flagged} (logical;
#'   \code{TRUE} means RiboTaper would have excluded this ORF for
#'   multi-mapping dependence and/or insufficient coverage breadth).
#' @export
#' @examples
#' data(chr22_anno)
#' data(rpfs)
#' uorf_ids <- names(chr22_anno$uORF)[chr22_anno$uORF]
#' mm <- flag_multimapping_uorfs(rpfs, chr22_anno, uorf_ids)
flag_multimapping_uorfs <- function(reads, anno, orf_ids, threshold = 0.3,
                                     min_coverage = 0) {
  stopifnot("NH" %in% colnames(S4Vectors::mcols(reads)))
  orf_ids <- unique(orf_ids)
  orfs <- anno$trspacecds[orf_ids]
  widths <- GenomicRanges::width(orfs)

  ov <- GenomicRanges::findOverlaps(reads, orfs, ignore.strand = TRUE)
  if (length(ov) == 0L) {
    return(data.frame(
      orf_id = orf_ids, pct_covered_all = 0, pct_covered_unique = 0,
      pct_covered_onlymulti = 0, low_coverage = TRUE,
      flagged = FALSE, stringsAsFactors = FALSE
    ))
  }
  q <- S4Vectors::queryHits(ov)
  orf_idx <- S4Vectors::subjectHits(ov)
  nh <- S4Vectors::mcols(reads)$NH[q]
  # Each read's full aligned span, clipped to the ORF it overlaps -
  # matching bedtools coverage's per-base accounting rather than marking
  # only a single (e.g. start) position per read.
  rel_start <- pmax(GenomicRanges::start(reads)[q] -
                      GenomicRanges::start(orfs)[orf_idx] + 1L, 1L)
  rel_end   <- pmin(GenomicRanges::end(reads)[q] -
                      GenomicRanges::start(orfs)[orf_idx] + 1L, widths[orf_idx])

  # Fraction of each ORF's nucleotides covered by >= 1 read footprint,
  # restricted to rows where `mask` is TRUE (coverage *breadth*, not read
  # count) - computed via interval union (IRanges::reduce) per ORF.
  cover_frac <- function(mask) {
    keep <- mask & (rel_end >= rel_start)
    if (!any(keep)) return(rep(0, length(orf_ids)))
    ir <- IRanges::IRanges(start = rel_start[keep], end = rel_end[keep])
    # factor(..., levels = seq_along(orf_ids)) guarantees split() emits one
    # (possibly empty) group per ORF, in orf_ids order, so the vapply result
    # already lines up positionally - no separate re-indexing needed.
    grp <- factor(orf_idx[keep], levels = seq_along(orf_ids))
    red <- IRanges::reduce(split(ir, grp))
    covered <- vapply(red, function(x) sum(IRanges::width(x)), numeric(1))
    covered / widths
  }
  pct_all    <- cover_frac(rep(TRUE, length(q)))
  pct_unique <- cover_frac(nh == 1L)
  pct_onlymulti <- pct_all - pct_unique
  low_coverage <- pct_all <= min_coverage

  data.frame(
    orf_id = orf_ids,
    pct_covered_all = pct_all,
    pct_covered_unique = pct_unique,
    pct_covered_onlymulti = pct_onlymulti,
    low_coverage = low_coverage,
    flagged = (pct_all > 0) &
      (((pct_onlymulti / pct_all) >= threshold) | low_coverage),
    stringsAsFactors = FALSE
  )
}
