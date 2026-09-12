#' Add 'correct offset' column to a set of reads
#'
#' Given rpfs, and named vectors of starts and stops, this adds the 'correct'
#' offset to the reads to use with determining p-site offsets
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param rpfs GRanges object with RPFs
#' @param anno annotation object
#'
#' @return the rpf granges but wth a 'cor_offset' column


add_cor_offset <- function(rpfs, anno) {
  cdsstarts <- anno$cdsstarts %>% unlist()
  cds_prestop_st <- anno$cds_prestop_st %>% unlist()

  rpfs <- rpfs %>%
    methods::as("GRanges") %>%
    addphase(cdsstarts) %>%
    {
      .$startoffset <- cdsstarts[as.vector(seqnames(.))] + .$phase - GenomicRanges::start(.)
      .$endoffset <- cds_prestop_st[as.vector(seqnames(.))] + .$phase - GenomicRanges::start(.)
      .$cor_offset <- ifelse(
        (.$startoffset < (.$readlen - 5)) & (.$startoffset > 5),
        .$startoffset,
        ifelse(
          (.$endoffset < (.$readlen - 5)) & (.$endoffset > 5),
          .$endoffset,
          NA
        )
      )
      .
    }
  rpfs
}



#' Get Offsets by maximizing p-sites within the CDS
#'
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param rpfs GRanges object with RPFs
#' @param anno annotation object
#'
#' @return a data frame with the best offset per read/phase

# TODO document the columns output by this function
# TODO maybe also shortedn this by removing arguments?
# TODO remove the hardcoded numbers

get_incl_max_offsets <- function(rpfs, anno) {
  stopifnot("readlen" %in% colnames(mcols(rpfs)))
  #
  rpfs <- rpfs %>% add_cor_offset(anno)
  #
  rpfs$readlen <- as.numeric(rpfs$readlen)
  if (!"score" %in% colnames(mcols(rpfs))) rpfs$score <- 1
  # data frame with numbers of reads gained at a specific offset
  gaindf <- rpfs %>%
    subset((startoffset < readlen - 5) & (startoffset > 5)) %>%
    mcols(.) %>%
    .[, c("startoffset", "score", "readlen", "phase")] %>%
    as.data.frame() %>%
    group_by(.data$phase, .data$readlen, .data$startoffset) %>%
    tally(wt = .data$score) %>%
    arrange(.data$startoffset) %>%
    group_by(.data$phase, .data$readlen) %>%
    mutate(gain = cumsum(.data$n)) %>%
    select(p_offset = 'startoffset', 'phase', 'gain')
  # data frame with numbers of reads lost at a specific offset
  lossdf <- rpfs %>%
    subset((.$endoffset < (readlen - 3)) & (.$endoffset > 3)) %>%
    mcols(.) %>%
    .[, c("endoffset", "score", "readlen", "phase")] %>%
    as.data.frame() %>%
    group_by(.data$phase, .data$readlen, .data$endoffset) %>%
    tally(wt = .data$score) %>%
    arrange(.data$endoffset) %>%
    group_by(.data$phase, .data$readlen) %>%
    mutate(loss = cumsum(.data$n)) %>%
    select('p_offset' = 'endoffset', 'phase', 'loss') %>%
    mutate('p_offset' = .data$p_offset + 3)
  # for a given offset, our gain is the start overlappers whose
  # offset is that or less, minus the stop guys whose offset is
  # less than that
  #
  netdf <- gaindf %>%
    full_join(lossdf, by = c("readlen", "p_offset", "phase")) %>%
    mutate(net = 
      replace_na(.data$gain, 0) - replace_na(.data$loss, 0)) %>%
    group_by(.data$readlen) %>%
    arrange(.data$readlen, .data$p_offset) %>%
    filter(!is.na(.data$net))
  #
  readlens <- rpfs$readlen %>% unique()
  #
  allpos <-
    lapply(readlens, function(rl) {
      expand.grid(
        readlen = rl,
        phase = 0:2,
        p_offset = 2:(rl - 3)
      )
    }) %>%
    bind_rows() %>%
    filter((.data$p_offset %% 3) == 0)
  netdf <- allpos %>% left_join(netdf,
    by = c("readlen", "p_offset", "phase")
  )
  netdf <- netdf %>%
    mutate(net = replace_na(.data$net, 0)) %>%
    mutate(twind = .data$net + lag(.data$net) + lead(.data$net))
  netdf <- netdf %>%
    group_by(.data$readlen) %>%
    filter(max(.data$net) > 32)
  best_offsets <- netdf %>%
    filter(!is.na(.data$twind)) %>%
    group_by(.data$readlen) %>%
    dplyr::slice((which.max(.data$twind) - 1):(which.max(.data$twind) + 1))
  #
  uniqueoffsetsfound <- best_offsets %>%
    group_by(.data$readlen, .data$phase) %>%
    filter(n() > 1) %>%
    nrow() %>%
    `==`(0)
  stopifnot(uniqueoffsetsfound)
  #
  best_offsets <- best_offsets%>%
    select(-'gain',-'loss',-'net',-'twind')
  best_offsets
}



#' Estimate P-site offsets via RiboTaper's metagene-peak method
#'
#' Reproduces the offset-selection procedure RiboTaper's own pipeline uses
#' (\code{create_metaplots.bash} + \code{metag.R}): build a metagene
#' histogram of the distance from each read's 5' end to the nearest
#' annotated start codon, separately per read length, and take the offset at
#' the histogram's peak (mode) - i.e. the single "obvious spike" a human
#' reads off the plotted PDF, rather than an inclusion-maximising search.
#' Unlike \code{\link{get_incl_max_offsets}} (Ribostan's default), this
#' resolves any integer nt offset and is not quantized to multiples of 3.
#'
#' RiboTaper itself picks one offset per read length only, with no
#' dependence on an individual read's ORF-relative frame phase (its
#' \code{P_sites_RNA_sites_calc.bash} takes a single offset per read
#' length). The same offset is reported here for phase 0/1/2 purely to
#' match \code{\link{get_psite_gr}}'s \code{(readlen, phase)} lookup schema;
#' this mirrors \code{\link{get_incl_max_offsets}}'s own output, which is
#' also phase-invariant in practice.
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param rpfs GRanges object with RPFs (as from \code{\link{get_readgr}}),
#'   with a \code{readlen} mcol.
#' @param anno annotation object (used for \code{anno$cdsstarts})
#' @param search_range Integer vector of length 2: the window of raw
#'   5'-end-to-start-codon distances (nt) to search for the peak in,
#'   matching RiboTaper's plausible P-site range (default \code{c(-35,-5)},
#'   i.e. the read's 5' end 5-35 nt upstream of the start codon).
#'
#' @return a data frame with columns \code{readlen}, \code{phase},
#'   \code{p_offset} - directly usable as \code{offsets_df} in
#'   \code{\link{get_psite_gr}}.
#' @export
#' @examples
#' data(chr22_anno)
#' data(rpfs)
#' offsets_df <- get_metagene_peak_offsets(rpfs, chr22_anno)
get_metagene_peak_offsets <- function(rpfs, anno, search_range = c(-35, -5)) {
  stopifnot("readlen" %in% colnames(mcols(rpfs)))
  cdsstarts <- anno$cdsstarts %>% unlist()
  starts <- cdsstarts[as.character(GenomicRanges::seqnames(rpfs))]
  distance <- GenomicRanges::start(rpfs) - starts
  keep <- !is.na(distance) & distance >= search_range[1] & distance <= search_range[2]
  df <- data.frame(
    readlen = as.numeric(rpfs$readlen[keep]),
    distance = distance[keep]
  )
  peaks <- df %>%
    dplyr::count(.data$readlen, .data$distance, name = "n") %>%
    dplyr::group_by(.data$readlen) %>%
    dplyr::slice_max(order_by = .data$n, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(readlen = .data$readlen, p_offset = -.data$distance)
  expand.grid(readlen = peaks$readlen, phase = 0:2) %>%
    dplyr::inner_join(peaks, by = "readlen") %>%
    dplyr::arrange(.data$readlen, .data$phase)
}


#' Get Offsets by maximizing p-sites within the CDS
#'
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param input GRanges object with RPFs or a bam file with RPFs
#' @param anno annotation object
#' @param method Method for p-site determination. \code{"cds_incl"}
#'   (default) maximises P-sites falling within annotated CDS regions
#'   (Ribostan's own algorithm). \code{"metagene_peak"} reproduces
#'   RiboTaper's own metagene-peak procedure instead - see
#'   \code{\link{get_metagene_peak_offsets}}.
#'
#' @return a dataframe with p-site offsets per readlength
#' @examples
#' cov <- get_readgr(ribobam, anno)
#' offsets_df <- get_offsets(cov, anno)
#' @export
#' @examples
#' data(chr22_anno)
#' data(rpfs)
#' data(offsets_df)
#' testoffsets_df <- get_offsets(rpfs, chr22_anno)
#'
# TODO Add reference to a-p-site paper
# input = ribobam
get_offsets <- function(input, anno, method = c("cds_incl", "metagene_peak")) {
  method <- match.arg(method)
  if (is.character(input) & (length(input) == 1)) {
    message("loading footprint data")
    rpfs <- get_readgr(input, anno)
  } else {
    stopifnot(is(input, "GRanges"))
    rpfs <- input
  }
  #
  if (method == "cds_incl") {
    best_offsets <- get_incl_max_offsets(rpfs, anno)
  } else if (method == "metagene_peak") {
    best_offsets <- get_metagene_peak_offsets(rpfs, anno)
  } else {
    stop("not yet implemented")
  }
  best_offsets
}
