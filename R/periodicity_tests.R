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


#' Run Multitaper Tests on All ORFs
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param psites A GRanges object containing psites
#' @param anno annotation object
#' @param n_cores number of cores to use
#'
#' @details This function applies a multitaper test to
#' @return a numeric vector with the spectral coefficient at 0.333...
#' and the pvalue for the test
#' @export
#' @examples
#' data(chr22_anno)
#' data(rpfs)
#' data(offsets_df)
#' psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
#' ftests <- ftest_orfs(psites %>% head(10000), chr22_anno, 
#'   n_cores=1)

ftest_orfs <- function(psites, anno, n_cores=1) {
  #
  orfs <- intersect(unique(psites$orf), names(anno$trspacecds))
  orfs <- anno$trspacecds[orfs]
  psitecov <- psites %>%
    {
      x <- .
      out <- GenomicFeatures::mapToTranscripts(., orfs, ignore.strand = TRUE)
      out <- out[x$orf[out$xHits] == names(orfs)[out$transcriptsHits]]
      coverage(out)
    }
  #
  psitecov <- psitecov[sum(psitecov > 0) > 1]
  # now run multitaper tests on our data using multiple
  # cores if available
  message("running multitaper tests, this will be slow for a full dataset...")
  spec_tests <- psitecov %>% 
    parallel::mclapply(F = ftestvect, mc.cores = n_cores)
  # now format the output
  # ftestvect returns a named vector (Fstat, spec_coef, p.value).
  spec_test_df <- spec_tests %>%
    simplify2array() %>%
    t() %>%
    as.data.frame()
  spec_test_df$orf_id <- rownames(spec_test_df)
  rownames(spec_test_df) <- NULL
  # Length-normalised spectral coefficient, preserved from legacy output.
  orflens <- width(orfs[spec_test_df$orf_id])
  spec_test_df$spec_coef <- sqrt(spec_test_df$spec_coef) / orflens
  # Benjamini-Hochberg FDR control across the tested ORFs. Applied here
  # (before NA-padding to all anno$cdsgrl rows) so only ORFs actually
  # evaluated contribute to the correction.
  spec_test_df$q.value <- p.adjust(spec_test_df$p.value, method = "BH")
  spec_test_df <- spec_test_df[, c("orf_id", "Fstat", "spec_coef",
                                   "p.value", "q.value")]
  # put in NA values for things we couldn't test
  testdf <- tibble(orf_id = names(anno$cdsgrl)) %>%
    left_join(spec_test_df, by = "orf_id")
  testdf
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



#' Run Multitaper tests on a set of ORFs
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param psites GRanges object with psite information
#' @param anno  annotation object
#' @param remove  whether to remove non-periodic uORFs
#' @param alpha significance threshold applied to `q.value` (BH-corrected);
#'   default 0.05
#' @param n_cores  number of cores to use
#'
#' @details
#' Before the (slow) multitaper F-test is run, nested/redundant uORF
#' candidates that share a stop position with another candidate on the
#' same transcript are pruned down to one representative per group, by
#' P-site read density (see \code{\link{prune_nested_uorfs}}). The
#' F-test is then applied to every surviving uORF. Both the raw p-value
#' and the Benjamini-Hochberg q-value (FDR-corrected across the tested
#' uORFs) are attached to `anno$trspacecds` as mcols. When
#' `remove = TRUE`, uORFs with `q.value >= alpha` (including pruned or
#' otherwise untested candidates) are dropped from the annotation.
#' @return an annotation object with Fstat/spec_coef/p.value/q.value
#'   attached per-ORF, and if `remove = TRUE` the non-periodic uORFs
#'   filtered out.
#' @export
#' @examples
#' data(chr22_anno)
#' data(rpfs)
#' data(offsets_df)
#' psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
#' filteredanno <- periodicity_filter_uORFs(psites, chr22_anno)

periodicity_filter_uORFs <- function(psites, anno, remove = TRUE,
                                     alpha = 0.05, n_cores = 1) {
  stopifnot(!is.null(anno$uORF))
  uORFs <- anno$uORF
  uORFs <- unique(names(uORFs[uORFs]))
  survivors <- prune_nested_uorfs(psites, anno, uORFs)
  ftestdf <- ftest_orfs(psites, subset_annotation(anno, survivors),
                        n_cores = n_cores)
  # Initialise per-ORF metadata columns if missing. ftest_orfs now returns
  # Fstat, spec_coef, p.value, q.value; historically it returned only
  # spec_coef and p.value.
  for (col in c("Fstat", "spec_coef", "p.value", "q.value")) {
    if (is.null(mcols(anno$trspacecds)[[col]]))
      mcols(anno$trspacecds)[[col]] <- NA_real_
  }
  mcols(anno$trspacecds[ftestdf$orf_id]) <- ftestdf %>% select(-"orf_id")
  if (remove) {
    periodic_uORFs <- ftestdf %>%
      filter(!is.na(.data$q.value) & .data$q.value < alpha) %>%
      .$orf_id
    non_periodic_uORFs <- setdiff(uORFs, periodic_uORFs)
    orfs_to_keep <- setdiff(names(anno$trspacecds), non_periodic_uORFs)
    anno <- subset_annotation(anno, orfs_to_keep)
  }
  anno
}
