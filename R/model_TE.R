#' Estimate per codon occupancies
#'
#' Uses sampled RPF coverage gr object, offsets data frame, and an annotation
#' object
#' to get estimates of per-codon occupancy (which should be proportional to
#' dwell times)
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param anno - an annotation object
#' @param offsets_df - a data frame with numeric columns readlen,phase,offset
#' @param anno - an annotation object
#' @return a data frame with columns transcript_id, sum_occ
#'
#' @details Use estimates of per codon dwell time to get the mean dwell time
#' over transcrips

get_codposdf <- function(ribocovtrs, anno, fafile) {
  #
  fafileob <- Rsamtools::FaFile(fafile)
  #
  cdsgrl <- anno %>%
    subset(type == "CDS") %>%
    GenomicRanges::split(., .$transcript_id)
  #
  bestcdsseq <- cdsgrl[ribocovtrs] %>%
    sort_grl_st() %>%
    GenomicFeatures::extractTranscriptSeqs(x = fafileob, .)
  #
  codposdf <- lapply(bestcdsseq, function(cdsseq) {
    codonmat <- codons(cdsseq) %>%
      {
        cbind(pos = .@ranges@start, as.data.frame(.))
      } %>%
      identity()
  })
  codposdf %<>% bind_rows(.id = "transcript_id")
  codposdf
}

#' get a coverage RleList of psite coverage over a set of ORFs
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
  offsetcols <- c('readlen', 'phase', 'offset')
  offsets <- rpfs %>%
    # subset(readlen %in% offsets_df$readlen)%>%
    addphase(anno$cdsstarts) %>%
    {
      as.data.frame(mcols(.)[, c("readlen", "phase")])
    } %>%
    left_join(
      offsets_df %>% select(one_of(offsetcols)),
      by = c("readlen", "phase")
    ) %>%
    .$offset
  #
  psitecov <- rpfs[!is.na(offsets)] %>%
    resize(1, "start") %>%
    GenomicRanges::shift(offsets %>% keep(Negate(is.na))) %>%
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
#' @param rpf_cov - a data frame of RPF footprints
#' @param anno - an annotation object
#' @param fafile - genomic fasta file
#' @return 	a data frame with columns transcript_id, p_codon, a_codon, count

get_sitedf <- function(rpf_cov, anno, fafile) {
  #
  sitedf <- psitecov %>%
    lapply(. %>% matrix(nrow = 3) %>% colSums()) %>%
    stack() %>%
    set_colnames(c("count", "transcript_id"))
  #
  codposdf <- get_codposdf(names(psitecov), anno$anno, fafile)
  #
  stopifnot(
    all(unique(codposdf$transcript_id) == unique(sitedf$transcript_id))
  )
  #
  sitedf$p_codon <- codposdf$x
  sitedf$a_codon <- codposdf %>%
    GenomicRanges::split(., .$transcript_id) %>%
    map(~ lead(.$x)) %>%
    unlist()
  sitedf %>% select(transcript_id, p_codon, a_codon, count)
}

#' Estimate per codon occupancies
#'
#' Uses sampled RPF coverage gr object, offsets data frame, and an annotation
#' object
#' to get estimates of per-codon occupancy (which should be proportional to
#' dwell times)
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param rpfs - a data frame of RPF footprints
#' @param offsets_df - a data frame with numeric columns readlen,phase,offset
#' @param anno - an annotation object
#' @param n_genes - how many genes to use for the model - most highly expressed
#'  used first
#' @param method - which method to use to estimate the codon occupancies,
#' defaults to RUST
#' @return 	a data frame with columns codon,position,estimate,upper,lower,
#' p.value
#'
#' @export


get_codon_occs <- function(rpfs, offsets_df,
                           anno, n_genes = 1000, method = "linear") {
  #
  allpsitecov <- get_psitecov(rpfs, offsets_df, anno)
  #
  #
  psitecov <- allpsitecov[allpsitecov %>%
    sum() %>%
    `<`(10)]
  toptrs <- psitecov %>%
    mean() %>%
    sort() %>%
    tail(n_genes) %>%
    names()
  #
  psitecov <- psitecov[toptrs]
  #

  sitedf <- get_sitedf(psitecov, anno, fafile)
  library(MASS)
  #
  # message('faking data')
  # sitedf$count = 100
  # sitedf$count = sitedf$count * codonstrengths[sitedf$p_codon]
  # sitedf$count%<>%rpois(length(.),.)
  sitedf %<>% group_by(transcript_id) %>% mutate(trmean = mean(count))
  sitedf <- sitedf %<>% filter(trmean != 0)

  sitedf

  #
  if (method == "linear") {
    nbfit <- glm.nb(
      data = sitedf,
      formula = count ~ 0 + offset(log(trmean)) + p_codon + a_codon
    )
    #
    codon_occ_df <- broom::tidy(nbfit) %>%
      as.data.frame() %>%
      mutate(codon = term %>% str_replace("[pa]_codon", "")) %>%
      mutate(position = term %>% str_extract("^[pa]_codon")) %>%
      mutate(lower = estimate - 1.96 * std.error) %>%
      mutate(upper = estimate + 1.96 * std.error)
    codon_occ_df
    # codon_occ_df%>%filter(position=='p_codon')%>%
    # 	mutate(comp=codonstrengths[codon])%>%
    # 	{txtplot(.$comp,exp(.$estimate))}
  } else if (method == "RUST_glm") {
    sitedf %<>% mutate(rust = count > trmean, rtrmean = mean(rust))
    rustfit <- glm(
      data = sitedf,
      formula = rust ~ 0 + offset(log(rtrmean)) + p_codon + a_codon,
      family = "binomial"
    )
    rcodon_occ_df <- broom::tidy(rustfit) %>%
      as.data.frame() %>%
      mutate(codon = term %>% str_replace("[pa]_codon", "")) %>%
      mutate(position = term %>% str_extract("^[pa]_codon")) %>%
      mutate(lower = estimate - 1.96 * std.error) %>%
      mutate(upper = estimate + 1.96 * std.error)
    codon_occ_df <- rcodon_occ_df
    # codon_occ_df%>%filter(position=='p_codon')%>%
    # 	mutate(comp=codonstrengths[codon])%>%
    # 	{txtplot::txtplot(.$comp,.$estimate)}
  } else if (method == "RUST") {
    sitedf %<>% mutate(rust = count > trmean, rtrmean = mean(rust))
    p_ests <- sitedf %>%
      group_by(p_codon) %>%
      summarise(estimate = sum(rust) / sum(trmean)) %>%
      setNames(c("codon", "estimate")) %>%
      mutate(position = "p_codon")
    a_ests <- sitedf %>%
      group_by(p_codon) %>%
      summarise(estimate = sum(rust) / sum(trmean)) %>%
      setNames(c("codon", "estimate")) %>%
      mutate(position = "a_codon")
    codon_occ_df <- bind_rows(p_ests, a_ests)
    codon_occ_df %<>% mutate(upper = NA, lower = NA, p.value = NA)
    codon_occ_df
  } else {
    p_ests <- sitedf %>%
      group_by(p_codon) %>%
      summarise(estimate = sum(count) / sum(trmean)) %>%
      setNames(c("codon", "estimate")) %>%
      mutate(position = "p_codon")
    a_ests <- sitedf %>%
      group_by(p_codon) %>%
      summarise(estimate = sum(count) / sum(trmean)) %>%
      setNames(c("codon", "estimate")) %>%
      mutate(position = "a_codon")
    codon_occ_df <- bind_rows(p_ests, a_ests)
    codon_occ_df %<>% mutate(upper = NA, lower = NA, p.value = NA)
    codon_occ_df
  }
  codon_occ_df %>%
    ungroup() %>%
    select(codon, position, estimate, upper, lower, p.value)
}



#' Get the mean dwell time over all translated codons for a set of transcripts
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param anno - an annotation object
#' @param fafile - a fasta file with sequence for that anno object
#' @param codon_model - a dwell-time data frame with columns position, codon,
#' estimate
#' @return a data frame with columns transcript_id, sum_occ
#'
#' @details Use estimates of per codon dwell time to get the mean dwell time
#' over transcrips

# TODO - update this to be more fleixible in what a codon model is
get_predicted_codon_occs <- function(anno, fafile, codon_model) {
  #
  cdsgrl <- anno$anno %>%
    subset(type == "CDS") %>%
    GenomicRanges::split(., .$transcript_id)
  #
  fafileob <- Rsamtools::FaFile(fafile)
  #
  elongtrs <- unique(anno$anno$transcript_id)
  trsplit <- elongtrs %>% split(floor(seq_along(.) / 1000))
  tr_elong <- map_df(trsplit, function(elongtrs) {
    #
    bestcdsseq <- cdsgrl[elongtrs] %>%
      sort_grl_st() %>%
      GenomicFeatures::extractTranscriptSeqs(x = fafileob, .)

    codonfreqdf <- bestcdsseq %>%
      oligonucleotideFrequency(3, step = 3) %>%
      set_rownames(names(bestcdsseq)) %>%
      {
        . <- . / rowSums(.)
        .
      } %>%
      as.data.frame() %>%
      rownames_to_column("transcript_id") %>%
      pivot_longer(-transcript_id, names_to = "codon", values_to = "freq")

    ipos <- "a_codon"
    positions <- unique(codon_model$position)
    poselongs <- lapply(positions, function(ipos) {
      poselong <- codonfreqdf %>%
        left_join(
          codon_model %>% filter(position == ipos) %>%
            select(codon, estimate),
          by = "codon"
        ) %>%
        group_by(transcript_id) %>%
        summarise(estimate = weighted.mean(estimate, freq, na.rm = T))
      poselong
    })
    tr_elong <- poselongs %>%
      purrr::reduce(left_join, by = "transcript_id") %>%
      {
        tibble(transcript_id = .[[1]], sum_occ = rowSums(.[, -1]))
      }
    tr_elong
  })
  tr_elong
}

#' Aggregate transcript elongation scores to gene level
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param tr_elong elongation rates per transcript, data frame with columns
#' transcript_id, sum_occ
#' @param tpms a named vector of tpms per transcript
#' @param anno an annotation object
#' @return a data frame with columns sum_occ and gene_id
#'
#' @details This function aggregates transript level tpms to gene level,with
#' the mean elongation rate over all the
#' genes transcripts weighted by transcript abundance
#' @seealso \code{\link{get_cds_reads}}, \code{\link{get_readlens}}


gene_level_elong <- function(tr_elong, tpms, anno) {
  trgiddf <- anno$anno %>%
    mcols() %>%
    .[, c("transcript_id", "gene_id")] %>%
    as.data.frame() %>%
    distinct()
  #
  gn_elong <- tr_elong %>%
    left_join(trgiddf, by = "transcript_id") %>%
    group_by(gene_id) %>%
    left_join(tibble::enframe(tpms, "transcript_id", "tpm"), 
              by = "transcript_id") %>%
    ungroup() %>%
    mutate(tpm = replace_na(tpm, min(tpm, na.rm = T) * 0.01)) %>%
    group_by(gene_id) %>%
    summarise(sum_occ = weighted.mean(sum_occ, tpm, na.rm = T))
  gn_elong %>%
    ungroup() %>%
    select(gene_id, sum_occ)
}
