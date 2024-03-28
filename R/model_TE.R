#' Generate a data frame of orf_id,position,codon
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param orfnames - the ORFs to get sequence data for
#' @param anno - an annotation object
#' @return returns a data frame with columns orf_id, pos, codon

get_codposdf <- function(orfnames, anno) {
  #
  cdsgrl <- anno$cdsgrl
  #
  bestcdsseq <- cdsgrl[orfnames] %>%
    sort_grl_st() %>%
    GenomicFeatures::extractTranscriptSeqs(x = anno$fafileob, .)
  #
  codposdf <- lapply(bestcdsseq, function(cdsseq) {
    codonmat <- codons(cdsseq) %>%
      {
        data.frame(pos = .@ranges@start, codon = as.data.frame(.)$seq)
      }
  })
  codposdf[[1]]%>%head
  codposdf <- bind_rows(codposdf, .id = "orf_id")
  colnames(codposdf) <- c("orf_id", "pos", "codon")
  codposdf
}


#' Get a coverage RleList of psite coverage over a set of ORFs
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param rpfs - granges object with positions of footprints, in transcript
#' space
#' @param offsets_df - a data frame with numeric columns readlen,phase,offset
#' @param anno - an annotation object with at least one ORF per transcript
#' @return a numeric SimpleRleList for each ORF with p-site density


# TODO make sure this can handle uORFs, dORFs.
get_psitecov <- function(rpfs, offsets_df, anno) {
  offsetcols <- c("readlen", "phase", "p_offset")
  p_offsets <- rpfs %>%
    # subset(readlen %in% offsets_df$readlen)%>%
    addphase(anno$cdsstarts) %>%
    {
      as.data.frame(mcols(.)[, c("readlen", "phase")])
    } %>%
    left_join(
      offsets_df %>% select(one_of(offsetcols)),
      by = c("readlen", "phase")
    ) %>%
    .$p_offset
  #
  p_off_na <- is.na(p_offsets)
  psitecov <- rpfs[!p_off_na] %>%
    resize(1, "start") %>%
    GenomicRanges::shift(p_offsets[!p_off_na]) %>%
    coverage()
  psitecov[anno$trspacecds[names(psitecov)]]
}
#


#' create a data frame with a and p site codons, and rpf counts per
#' transcript/position
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param psite_cov - a data frame of RPF footprints
#' @param anno - an annotation object
#' @return a data frame with columns orf_id, p_codon, a_codon, count

get_sitedf <- function(psite_cov, anno) {
  #
  stopifnot(length(psite_cov)>0)
  stopifnot(length(anno$cdsgrl)>0)
  sitedf <- psite_cov %>%
    lapply(. %>% matrix(nrow = 3) %>% colSums()) %>%
    stack()
  colnames(sitedf) <- c("count", "orf_id")
  #
  codposdf <- get_codposdf(names(psite_cov), anno)
  #
  stopifnot(
    all(unique(codposdf$orf_id) == unique(sitedf$orf_id))
  )
  #
  sitedf$p_codon <- codposdf$codon
  sitedf$a_codon <- codposdf %>%
    split(., .$orf_id) %>%
    lapply(function(x) lead(x$codon)) %>%
    unlist()
  sitedf %>% select('orf_id', 'p_codon', 'a_codon', 'count')
}

#' Estimate per codon occupancies
#'
#' Uses sampled RPF coverage gr object, offsets data frame, and an annotation
#' object
#' to get estimates of per-codon occupancy (which are proportional to
#' dwell times)
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param psites - a data frame of inferred RPF p-sites
#' @param offsets_df - a data frame with numeric columns readlen,phase,offset
#' @param anno - an annotation object
#' @param n_genes - how many genes to use for the model - most highly expressed
#'  used first
#' @param method - which method to use to estimate the codon occupancies,
#' @param covthresh - how many reads to require in an ORF for it to be used
#' defaults to RUST
#' @return a data frame with columns codon,position,estimate,upper,lower,
#' p.value
#' @examples
#' data(chr22_anno)
#' data(rpfs)
#' data(offsets_df)
#' psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
#' rust_codon_occ_df <- get_codon_occs(psites, offsets_df, chr22_anno,
#'   n_genes = 200, method = "RUST"
#' )
#' @export

get_codon_occs <- function(psites, offsets_df,
                           anno, n_genes = 1000, covthresh=10, method = "linear") {
  highcov <- psites$orf %>%
    table() %>%
    .[. > covthresh]
  stopifnot(length(highcov)>n_genes)
  toporfs <- highcov %>%
    sort() %>%
    names %>%
    intersect(names(anno$trspacecds))%>%
    tail(n_genes)
  stopifnot(length(toporfs)==n_genes)
  unq_orfs <- anno$trspacecds[toporfs]
  psitecov <- psites %>%
    subset(orf %in% toporfs) %>%
    {
      x <- .
      out <- GenomicFeatures::mapToTranscripts(., unq_orfs, ignore.strand = TRUE)
      out <- out[x$orf[out$xHits] == names(unq_orfs)[out$transcriptsHits]]
      coverage(out)
    }
  psitecov <- psitecov[toporfs]
  #
  sitedf <- get_sitedf(psitecov, anno)
  #
  sitedf <- sitedf %>%
    group_by(.data$orf_id) %>%
    mutate(trmean = mean(.data$count))
  sitedf <- sitedf %>% filter(.data$trmean != 0)
  #
  if (method == "linear") {
    nbfit <- glm.nb(
      data = sitedf,
      formula = count ~ 0 + offset(log(trmean)) + p_codon + a_codon
    )
    #
    codon_occ_df <- broom::tidy(nbfit) %>%
      as.data.frame() %>%
      mutate(codon = .data$term %>% str_replace("[pa]_codon", "")) %>%
      mutate(position = .data$term %>% str_extract("^[pa]_codon")) %>%
      mutate(lower = .data$estimate - 1.96 * .data$std.error) %>%
      mutate(upper = .data$estimate + 1.96 * .data$std.error)
    codon_occ_df
  } else if (method == "RUST_glm") {
    sitedf <- sitedf %>% mutate(rust = count > trmean, rtrmean = mean(rust))
    rustfit <- glm(
      data = sitedf,
      formula = rust ~ 0 + offset(log(rtrmean)) + p_codon + a_codon,
      family = "binomial"
    )
    codon_occ_df <- broom::tidy(rustfit) %>%
      as.data.frame() %>%
      mutate(codon = .data$term %>% str_replace("[pa]_codon", "")) %>%
      mutate(position = .data$term %>% str_extract("^[pa]_codon")) %>%
      mutate(lower = .data$estimate - 1.96 * .data$std.error) %>%
      mutate(upper = .data$estimate + 1.96 * .data$std.error)
  } else if (method == "RUST") {
    sitedf <- sitedf %>% mutate(rust = .data$count > .data$trmean,
        rtrmean = mean(rust))
    p_ests <- sitedf %>%
      group_by(.data$p_codon) %>%
      summarise(estimate = sum(.data$rust) / sum(.data$trmean)) %>%
      setNames(c("codon", "estimate")) %>%
      mutate(position = "p_codon")
    a_ests <- sitedf %>%
      group_by(.data$p_codon) %>%
      summarise(estimate = sum(.data$rust) / sum(.data$trmean)) %>%
      setNames(c("codon", "estimate")) %>%
      mutate(position = "a_codon")
    codon_occ_df <- bind_rows(p_ests, a_ests)
    codon_occ_df <- codon_occ_df %>% 
      mutate(upper = NA, lower = NA, p.value = NA)
    codon_occ_df
  } else {
    p_ests <- sitedf %>%
      group_by(.data$p_codon) %>%
      summarise(estimate = sum(.data$count) / sum(.data$trmean)) %>%
      setNames(c("codon", "estimate")) %>%
      mutate(position = "p_codon")
    a_ests <- sitedf %>%
      group_by(.data$p_codon) %>%
      summarise(estimate = sum(.data$count) / sum(.data$trmean)) %>%
      setNames(c("codon", "estimate")) %>%
      mutate(position = "a_codon")
    codon_occ_df <- bind_rows(p_ests, a_ests)
    codon_occ_df <- mutate(codon_occ_df, upper = NA, lower = NA, p.value = NA)
    codon_occ_df
  }
  codon_occ_df %>%
    ungroup() %>%
    select('codon', 'position', 'estimate', 'upper', 'lower', 'p.value')
}



#' Get the mean dwell time over all translated codons for a set of transcripts
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param anno - an annotation object
#' @param codon_model - a dwell-time data frame with columns position, codon,
#' estimate
#' @return a data frame with columns orf_id, mean_occ
#'
#' @details Use estimates of per codon dwell time to get the mean dwell time
#' @examples
#'
#'
#' data(chr22_anno)
#' data(rpfs)
#' data(offsets_df)
#' psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
#' rust_codon_occ_df <- get_codon_occs(psites, offsets_df, chr22_anno,
#'   n_genes = 1000, method = "RUST"
#' )
#' orf_elong <- get_orf_elong(chr22_anno, rust_codon_occ_df)
#' gn_elong <- gene_level_elong(orf_elong, ritpms, anno)
#' @export
get_orf_elong <- function(anno, codon_model) {
  #
  cdsgrl <- anno$cdsgrl
  #
  fafileob <- anno$fafileob
  Rsamtools::indexFa(fafileob$path)
  #
  elong_orfs <- names(cdsgrl)
  trsplit <- elong_orfs %>% split(floor(seq_along(.) / 1000))
  tr_elong <- lapply(trsplit, function(elong_orfs) {
    #
    bestcdsseq <- cdsgrl[elong_orfs] %>%
      sort_grl_st() %>%
      GenomicFeatures::extractTranscriptSeqs(x = fafileob, .)

    codonfreqdf <- bestcdsseq %>%
      oligonucleotideFrequency(3, step = 3) %>%
      {
        rownames(.) <- (names(bestcdsseq))
        .
      } %>%
      {
        . <- . / rowSums(.)
        .
      } %>%
      as.data.frame()
    codonfreqdf$orf_id <- rownames(codonfreqdf)
    rownames(codonfreqdf) <- NULL
    codonfreqdf <- codonfreqdf %>%
      tidyr::pivot_longer(-orf_id, names_to = "codon", values_to = "freq")

    ipos <- "a_codon"
    positions <- unique(codon_model$position)
    poselongs <- lapply(positions, function(ipos) {
      poselong <- codonfreqdf %>%
        left_join(
          codon_model %>% filter(position == ipos) %>%
            select('codon', 'estimate'),
          by = "codon"
        ) %>%
        group_by(orf_id) %>%
        summarise(estimate = 
          weighted.mean(.data$estimate, .data$freq, na.rm = TRUE))
      poselong
    })
    tr_elong <- poselongs %>%
      Reduce(f = function(x, y) left_join(x, y, by = "orf_id")) %>%
      {
        data.frame(orf_id = .[[1]], mean_occ = rowMeans(.[, -1]))
      }
    tr_elong
  })
  tr_elong <- bind_rows(tr_elong)
  tr_elong
}

#' Aggregate transcript elongation scores to gene level
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param tr_elong elongation rates per transcript, data frame with columns
#' orf_id, mean_occ
#' @param ritpms a named vector of ritpms per transcript
#' @param anno an annotation object
#' @param exclude_uORFs exclude uORFs from the gene level total;Defaults to TRUE
#' @return a data frame with columns mean_occ and gene_id
#'
#' @details This function aggregates transript level ritpms to gene level,with
#' the mean elongation rate over all the
#' genes transcripts weighted by transcript abundance
#' @seealso \code{\link{get_cds_reads}}, \code{\link{get_readlens}}
#' @export


gene_level_elong <- function(tr_elong, ritpms, anno, exclude_uORFs = TRUE) {
  trgiddf <- anno$trgiddf
  if (exclude_uORFs) {
    non_uORFs <- trgiddf %>%
      subset(!uORF) %>%
      .$orf_id
    tr_elong <- tr_elong %>% subset(orf_id %in% non_uORFs)
  }
  #
  gn_elong <- tr_elong %>%
    left_join(trgiddf, by = "orf_id") %>%
    group_by(.data$gene_id) %>%
    left_join(tibble::enframe(ritpms, "orf_id", "ritpm"),
      by = "orf_id"
    ) %>%
    ungroup() %>%
    mutate(ritpm = 
      replace_na(.data$ritpm, min(.data$ritpm, na.rm = TRUE) * 0.01)) %>%
    group_by(.data$gene_id) %>%
    summarise(mean_occ = weighted.mean(mean_occ, ritpm, na.rm = TRUE))
  gn_elong %>%
    ungroup() %>%
    select('gene_id', 'mean_occ')
}
