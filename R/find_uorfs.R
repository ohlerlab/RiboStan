#' Find upstream ORFs (uORFs) from 5' UTR annotation
#'
#' Self-contained replacement for \code{ORFik::findUORFs}. Finds ORFs within
#' 5' UTR regions, optionally extending the search space into the downstream
#' CDS and filtering out ORFs that coincide with the annotated CDS.
#'
#' Procedure:
#' 1. Build search space: 5' UTRs, extended into the CDS if \code{cds} is
#'    supplied (so that ORFs bridging the 5' UTR / CDS boundary are found).
#' 2. Extract transcript sequences from the FASTA file.
#' 3. Scan each sequence for in-frame start → stop codon pairs.
#' 4. Map ORF positions from transcript coordinates back to genomic coordinates.
#' 5. Filter out ORFs that are fully within, share a start with, share a stop
#'    with, or have a start inside any annotated CDS element.
#'
#' @param fiveUTRs \code{GRangesList} of 5' UTR exons, one element per
#'   transcript.
#' @param fa Character path to an indexed genome FASTA file, or a
#'   \code{FaFile} object.
#' @param startCodon Character vector; start codons to search for. Defaults
#'   to \code{"ATG"} only, matching the annotated-CDS convention. Callers
#'   may opt in to non-canonical starts (e.g. \code{c("ATG","CTG","GTG")})
#'   but should note this substantially increases the number of candidate
#'   uORFs found.
#' @param stopCodons Character vector; stop codons to recognise
#'   (default \code{c("TAA", "TAG", "TGA")}).
#' @param longestORF Logical; if \code{TRUE}, keep only the longest ORF sharing
#'   each stop codon position.
#' @param minimumLength Integer; minimum body length in nucleotides between the
#'   end of the start codon and the beginning of the stop codon.  The default
#'   (\code{3}) requires at least one in-frame codon between start and stop,
#'   matching RiboTaper's minimum track length (\code{length > 5} nt).
#' @param cds \code{GRangesList} of CDS exons.  When supplied, the CDS is
#'   appended to the 5' UTR search space and detected ORFs that overlap the
#'   annotated CDS boundary (fully within it, sharing its start, sharing
#'   its stop, or starting inside it) are removed. See
#'   \code{\link{.filter_uorfs}}.
#' @return A \code{GRangesList} of uORFs in genomic coordinates, named
#'   \code{"<transcript_id>_<N>"}. Each element carries a metadata column
#'   \code{start_codon} giving the start codon that was matched.
#' @keywords internal
find_uorfs <- function(fiveUTRs, fa,
                       startCodon    = "ATG",
                       stopCodons    = c("TAA", "TAG", "TGA"),
                       longestORF    = FALSE,
                       minimumLength = 3L,
                       cds           = NULL) {

  ## 1. Build search space ------------------------------------------------
  search_space <- fiveUTRs
  if (!is.null(cds) && length(cds) > 0L) {
    search_space <- .append_cds_to_utrs(fiveUTRs, cds)
  }

  ## 2. Extract transcript sequences ---------------------------------------
  fa_file <- if (is.character(fa)) Rsamtools::FaFile(fa) else fa
  seqs <- toupper(as.character(
    GenomicFeatures::extractTranscriptSeqs(fa_file, search_space)
  ))

  ## 3. Find ORFs in transcript coordinates (via C++) ----------------------
  start_codons <- toupper(as.character(startCodon))
  raw <- find_orfs_cpp(
    seqs         = seqs,
    start_codons = start_codons,
    stop_codons  = toupper(stopCodons),
    min_body     = as.integer(minimumLength)
  )

  if (length(raw$starts) == 0L) return(GenomicRanges::GRangesList())

  # Optionally retain only the longest ORF per (sequence, stop-site) pair
  if (isTRUE(longestORF)) {
    widths  <- raw$ends - raw$starts
    key     <- paste0(raw$indices, "_", raw$ends)
    keep    <- unlist(
      tapply(seq_along(key), key, function(i) i[which.max(widths[i])]),
      use.names = FALSE
    )
    raw <- lapply(raw, `[`, keep)
  }

  start_codon_hit <- start_codons[raw$start_codon_idx]

  ## 4. Map to genomic coordinates -----------------------------------------
  # tx_idx: which element of search_space each ORF belongs to (1-based)
  tx_idx     <- raw$indices
  all_orf_ir <- IRanges::IRanges(start = raw$starts, end = raw$ends)

  has_orfs   <- logical(length(search_space))
  has_orfs[unique(tx_idx)] <- TRUE
  tx_lengths <- tabulate(tx_idx, nbins = length(search_space))[unique(tx_idx)]

  # pmapFromTranscripts maps orf_gr[i] through search_space[tx_idx[i]]
  orf_gr      <- GenomicRanges::GRanges(
    seqnames = names(search_space)[tx_idx],
    ranges   = all_orf_ir
  )
  genomic_grl <- GenomicFeatures::pmapFromTranscripts(orf_gr,
                                                       search_space[tx_idx])
  # pmapFromTranscripts() returns one row per exon of the *parent
  # transcript* (search_space[tx_idx[i]]), not just the exon(s) the ORF
  # actually overlaps: non-overlapping exons come back as zero-width
  # filler ranges flagged hit=FALSE. Left unfiltered, every uORF's
  # GRangesList element balloons to the size of its transcript's full
  # exon count (a major source of memory/CPU blowup for highly-exonic
  # genes) and the spurious filler ranges break the 5'->3' contiguity
  # that load_annotation()'s downstream transcript-space mapping
  # (get_trspace_cds()) requires.
  hit <- genomic_grl@unlistData$hit
  genomic_grl <- genomic_grl[IRanges::relist(hit, genomic_grl)]
  # Every surviving row is a hit by construction now, so the column is
  # redundant going forward; drop it rather than leak it into cdsgrl.
  genomic_grl@unlistData$hit <- NULL

  ## 5. Name ORFs: <txname>_<rank> -----------------------------------------
  # Sequential rank within each transcript
  tx_names_rep <- names(search_space)[tx_idx]
  rank         <- unlist(
    tapply(seq_along(tx_idx), tx_idx, seq_along),
    use.names = FALSE
  )
  names(genomic_grl)  <- paste0(tx_names_rep, "_", rank)
  seqinfo(genomic_grl) <- seqinfo(fiveUTRs)
  # Record which start codon produced each uORF. Attach at both the GRL
  # level (one value per uORF) and the inner GRanges level (one value
  # per exon, same value for all exons of a uORF) so the column survives
  # the unlist() / split() cycle in load_annotation.
  S4Vectors::mcols(genomic_grl)$start_codon <- start_codon_hit
  S4Vectors::mcols(genomic_grl@unlistData)$start_codon <-
    rep(start_codon_hit, S4Vectors::elementNROWS(genomic_grl))

  ## 6. Filter ORFs that overlap the annotated CDS boundary ----------------
  if (!is.null(cds) && length(genomic_grl) > 0L) {
    genomic_grl <- .filter_uorfs(genomic_grl, cds)
  }

  genomic_grl
}


# ── Internal helpers ──────────────────────────────────────────────────────────

#' Scan a single nucleotide string for ORFs
#'
#' @param seq Character; uppercase DNA sequence.
#' @param start_codon Character; start codon (default \code{"ATG"}).
#' @param stop_codons Character vector; stop codons.
#' @param min_length Integer; minimum body length in nucleotides.
#' @param longest_orf Logical; return only longest ORF per stop site.
#' @return \code{IRanges} of ORF coordinates (1-based; both start and stop
#'   codons included).
#' @keywords internal
.find_orfs_in_seq <- function(seq,
                               start_codon = "ATG",
                               stop_codons = c("TAA", "TAG", "TGA"),
                               min_length  = 0L,
                               longest_orf = FALSE) {
  n       <- nchar(seq)
  starts  <- gregexpr(start_codon, seq, fixed = TRUE)[[1L]]
  if (starts[[1L]] == -1L) return(IRanges::IRanges())

  orf_starts <- integer(0)
  orf_ends   <- integer(0)

  for (s in as.integer(starts)) {
    pos <- s + 3L
    while (pos + 2L <= n) {
      codon <- substr(seq, pos, pos + 2L)
      if (codon %in% stop_codons) {
        # body = nt between end of start and start of stop = pos - s - 3
        if ((pos - s - 3L) >= min_length) {
          orf_starts <- c(orf_starts, s)
          orf_ends   <- c(orf_ends, pos + 2L)
        }
        break
      }
      pos <- pos + 3L
    }
  }

  if (length(orf_starts) == 0L) return(IRanges::IRanges())

  if (longest_orf) {
    widths <- orf_ends - orf_starts
    keep   <- unlist(
      tapply(seq_along(orf_ends), orf_ends,
             function(i) i[which.max(widths[i])]),
      use.names = FALSE
    )
    orf_starts <- orf_starts[keep]
    orf_ends   <- orf_ends[keep]
  }

  IRanges::IRanges(start = orf_starts, end = orf_ends)
}


#' Append CDS exons to each 5' UTR to extend the ORF search space
#'
#' Allows detection of ORFs that start in the 5' UTR and extend into the CDS.
#' Only transcripts present in both \code{fiveUTRs} and \code{cds} are
#' extended; others are returned unchanged.
#'
#' @param fiveUTRs GRangesList of 5' UTR exons.
#' @param cds GRangesList of CDS exons.
#' @return GRangesList with CDS appended (and reduced) for matching transcripts.
#' @keywords internal
.append_cds_to_utrs <- function(fiveUTRs, cds) {
  shared <- intersect(names(fiveUTRs), names(cds))
  if (length(shared) == 0L) return(fiveUTRs)

  # Combine element-wise (equivalent to ORFik::pc)
  combined <- GenomicRanges::GRangesList(
    mapply(c, fiveUTRs[shared], cds[shared], SIMPLIFY = FALSE)
  )
  names(combined) <- shared

  reduced <- GenomicRanges::reduce(combined, ignore.strand = FALSE)
  names(reduced)   <- shared
  seqinfo(reduced) <- seqinfo(fiveUTRs)
  # reduce() always returns each element's ranges sorted by ascending
  # genomic start, which is the wrong order for minus-strand transcripts
  # (5'->3' order there means descending genomic start). Downstream
  # consumers (pmapFromTranscripts inside find_uorfs, and resize_grl in
  # load_annotation) require exons in proper 5'->3' order.
  reduced <- sort_grl_st(reduced)

  only_utrs <- fiveUTRs[setdiff(names(fiveUTRs), shared)]
  # reduce() drops all inner (exon-level) mcols from 'reduced', while
  # 'only_utrs' still carries the exon_id/exon_name/exon_rank columns
  # inherited from the TxDb. Combining the two coerces exon_rank into a
  # logical NA vector (rather than leaving it absent/NULL), which
  # extractTranscriptSeqs()'s internal exon_rank check then rejects as
  # "not numeric". None of that per-exon annotation is used downstream
  # (only ranges/strand/seqnames matter for sequence extraction), so
  # drop it from both pieces before combining.
  S4Vectors::mcols(only_utrs, level = "within") <- NULL
  c(reduced, only_utrs)
}


#' Remove uORFs that coincide with annotated CDS boundaries
#'
#' Four checks are applied in sequence:
#' \enumerate{
#'   \item ORF is fully within a CDS element.
#'   \item ORF stop site equals a CDS stop site.
#'   \item ORF start site equals a CDS start site.
#'   \item ORF start site falls inside a CDS element.
#' }
#'
#' @param uorfs GRangesList of candidate uORFs.
#' @param cds GRangesList of annotated CDS.
#' @return Filtered GRangesList.
#' @keywords internal
.filter_uorfs <- function(uorfs, cds) {
  if (length(uorfs) == 0L) return(uorfs)

  drop <- function(grl, hits) {
    idx <- unique(S4Vectors::queryHits(hits))
    if (length(idx) > 0L) grl[-idx] else grl
  }

  uorfs <- drop(uorfs,
    GenomicRanges::findOverlaps(uorfs, cds, type = "within"))
  if (length(uorfs) == 0L) return(uorfs)

  uorfs <- drop(uorfs,
    GenomicRanges::findOverlaps(
      .stop_sites_gr(uorfs), .stop_sites_gr(cds), type = "within"))
  if (length(uorfs) == 0L) return(uorfs)

  uorfs <- drop(uorfs,
    GenomicRanges::findOverlaps(
      .start_sites_gr(uorfs), .start_sites_gr(cds), type = "within"))
  if (length(uorfs) == 0L) return(uorfs)

  drop(uorfs,
    GenomicRanges::findOverlaps(.start_sites_gr(uorfs), cds, type = "within"))
}


#' Width-1 GRanges at the genomic start (5'-most base) of each GRL element
#' @keywords internal
.start_sites_gr <- function(grl) {
  ul     <- unlist(grl, use.names = FALSE)
  is_neg <- as.character(GenomicRanges::strand(ul)) == "-"
  # + strand: start of the leftmost exon; - strand: end of the rightmost exon
  coords <- ifelse(is_neg,
                   GenomicRanges::end(ul),
                   GenomicRanges::start(ul))
  grp    <- factor(rep(seq_along(grl), S4Vectors::elementNROWS(grl)),
                   levels = seq_along(grl))
  sel    <- unlist(
    tapply(seq_along(ul), grp, function(g) {
      if (is_neg[g[[1L]]]) g[which.max(coords[g])]
      else                  g[which.min(coords[g])]
    }),
    use.names = FALSE
  )
  GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(ul)[sel],
    ranges   = IRanges::IRanges(start = coords[sel], width = 1L),
    strand   = GenomicRanges::strand(ul)[sel]
  )
}


#' Width-1 GRanges at the genomic stop (3'-most base) of each GRL element
#' @keywords internal
.stop_sites_gr <- function(grl) {
  ul     <- unlist(grl, use.names = FALSE)
  is_neg <- as.character(GenomicRanges::strand(ul)) == "-"
  # + strand: end of rightmost exon; - strand: start of leftmost exon
  coords <- ifelse(is_neg,
                   GenomicRanges::start(ul),
                   GenomicRanges::end(ul))
  grp    <- factor(rep(seq_along(grl), S4Vectors::elementNROWS(grl)),
                   levels = seq_along(grl))
  sel    <- unlist(
    tapply(seq_along(ul), grp, function(g) {
      if (is_neg[g[[1L]]]) g[which.min(coords[g])]
      else                  g[which.max(coords[g])]
    }),
    use.names = FALSE
  )
  GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(ul)[sel],
    ranges   = IRanges::IRanges(start = coords[sel], width = 1L),
    strand   = GenomicRanges::strand(ul)[sel]
  )
}
