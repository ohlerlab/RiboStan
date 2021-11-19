#' Get A Granges object with locations of each codon in the given cds
#'
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param anno An annotation object
#' @param offsets_df - a data frame with numeric columns readlen,phase,offset
#' @param anno - an annotation object
#' @param n_wind_l_ext The 5' extent of the window around the codon to return
#' Defaults to: 45
#' @param n_wind_r_ext The 3' extent of the window around the codon to return
#' Defaults to: 9
#' @param n_start_buff The number of bp at the start of each cds to exclude
#' Defaults to 60
#' @param n_end_buff The number of bp at the end of each cds to exclude
#' Defaults to 60
#' @return a GRanges object with the codon identity in the names, and the
#'
#' @details Use estimates of per codon dwell time to get the mean dwell time
#' over transcrips

# TODO stop missusing the name slot, maybe move some of this out so we don't
# have so many arguments

get_cds_codons <- function(anno, fafileob,
                           n_wind_l_ext = 45,
                           n_wind_r_ext = 9,
                           n_start_buff = 60,
                           n_end_buff = 60) {
  allcodons <- getGeneticCode()
  allcodons <- allcodons %>% setNames(names(allcodons))
  cds_codons <- lapply(allcodons, {
    exonsgrl <- anno$exonsgrl
    trspacecds <- anno$trspacecds
    # define inner cds - cds but with start and end trimmed off.
    innercds <- trspacecds %>%
      subset(width > (3 + n_start_buff + n_end_buff)) %>%
      resize(width(.) - n_start_buff, "end") %>%
      resize(width(.) - n_end_buff, "start")
    exonseq <- exonsgrl %>% extractTranscriptSeqs(x = fafileob)
    # get matches to that codon (any frame)
    codmatches <- vmatchPattern(pattern = codon, exonseq) # exclude the start ccodon
    matchgr <- codmatches %>%
      unlist() %>%
      GRanges(names(.), .)
    # select only the in frame codons
    matchcdsstarts <- start(trspacecds[as.vector(seqnames(matchgr))])
    matchgr$cdspos <- start(matchgr) - matchcdsstarts
    matchgr %<>% subset(cdspos %% 3 == 0)
    # put seqlengths on this object
    seqlengths(matchgr) <- exonsgrl %>%
      width() %>%
      .[seqlevels(matchgr)] %>%
      sum()
    # use only matches in the inner cds
    matchgr <- matchgr %>% subsetByOverlaps(innercds)
    # expand the windows around these codons
    codmatchwindows <- matchgr %>%
      resize(n_wind_l_ext, "end") %>%
      resize(n_wind_r_ext, "start")
    # require these to be inside our cds
    codmatchwindows <- codmatchwindows[!is_out_of_bounds(codmatchwindows)]
    codmatchwindows %<>% subsetByOverlaps(innercds)
    innercdspos <- innercds %>%
      {
        strand(.) <- "+"
        .
      }
    # lift them to inner cds coordinates
    icds_codons <- allcodlist %>% mapToTranscripts(innercdspos)
    mcols(icds_codons) <- NULL
    icds_codons
  })
  cds_codons <- cds_codons %>%
    GRangesList() %>%
    unlist()
  cds_codons$codon <- names(cds_codons) %>%
    str_split("\\.") %>%
    map_chr(1)
  names(cds_codons) <- NULL
  cds_codons
}

#' Get A Granges object with locations of each codon in the given cds
#'
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param rpf_cov A list of SimpleRleLists - RPF coverage split by readlength
#' @param cds_codons - GRanges object containing windows in which to get rust scores
#'
#' @details Use estimates of per codon dwell time to get the mean dwell time over transcrips
#' @return A data frame with expected and actual rust score per position in window, codon, readlength

get_cov_rust_scores <- function(rpf_cov, cds_codons) {
  cdssampfpcov %>% map_df(.id = "readlen", function(cdsrlfpcov) {
    rustcdsrlfpcov <- cdsrlfpcov
    rustcdsrlfpcov <- rustcdsrlfpcov > mean(rustcdsrlfpcov)
    nz_trs <- any(rustcdsrlfpcov) %>% names(.)[.]
    #
    cds_codons_nz <- cds_codons %>% subset(seqnames %in% nz_trs)
    codtrs <- seqnames(cds_codons_nz) %>% as.character()
    cat(".")
    # calculate ro vals
    ro_cl <- rustcdsrlfpcov[cds_codons_nz] %>%
      split(cds_codons_nz$codon) %>%
      lapply(as.matrix) %>%
      map(colMeans)
    # also get evals
    tr_rust_evals <- rustcdsrlfpcov %>% mean()
    re_c <- tr_rust_evals[codtrs] %>%
      split(cds_codons_nz$codon) %>%
      map_dbl(mean)
    ro_cl <- ro_cl %>% map_df(.id = "codon", tibble::enframe, 
                              "position", "ro_cl")
    re_c <- tibble::enframe(re_c, "codon", "re_c")
    ro_cl %>% left_join(re_c, by = "codon")
  })
}

#' Get A Granges object with locations of each codon in the given cds
#'
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param rpf_cov A list of SimpleRleLists - RPF coverage split by readlength
#' @param cds_codons - GRanges object containing windows in which to get rust
#' scores
#'
#' @details Use estimates of per codon dwell time to get the mean dwell time
#' over transcrips
#' @return A data frame with expected and actual rust score per position in
#' window, codon, readlength

get_sample_profs <- function(covgrs, cds_codons, n_wind_l_ext = 45) {
  stopifnot(!is.null(names(covgrs)))
  cdsfpcovlist <- lapply(covgrs, function(covgr) {
    coverage(split(resize(covgr, 1, "start"), covgr$readlen))
  })
  rust_roel <- mclapply(
    mc.cores = 4, SIMPLIFY = F, cdsfpcovlist, cds_codons = cds_codons,
    get_cov_rust_scores
  )
  # https://www.nature.com/articles/ncomms12915#Sec10 see equation 3
  # of RUST paper
  frustprofilelist <- rust_roel %>%
    # frustprofilelist <- cn_norm%>%
    map_df(.id = "sample", . %>% bind_rows(.id = "readlen"))
  #
  frustprofilelist %<>% mutate(position = position - 1 - (n_wind_l_ext))
  frustprofilelist %<>% filter(!codon %in% c("TAG", "TAA", "TGA"))
  frustprofilelist %<>% mutate(count = ro_cl / re_c)
  frustprofilelist$nreadlen <- frustprofilelist$readlen %>% as.numeric()
  # frustprofilelist$readlen%<>%str_replace('^(\\d)','rl\\1')
  frustprofilelist
}


#' Calculate KL-divergence from a data frame of rust expected and actual scores
#'
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param frustprofilelist
#'
#' @return A data frame with expected and actual rust score per position in
#' window, codon, readlength

get_kl_df <- function(frustprofilelist) {
  frustprofilelist %>%
    group_by(sample, nreadlen, position) %>%
    mutate(ro_cl = ro_cl / sum(ro_cl), re_c = re_c / sum(re_c)) %>%
    summarise(KL = sum(ro_cl * log2(ro_cl / re_c)))
}

#' Get the most frequent element from a numeric vector
#'
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param x a vector
#'
#' @return the single most common value in the vector

most_freq <- function(x) {
  x %>%
    table() %>%
    sort() %>%
    names() %>%
    as.numeric() %>%
    tail(1)
}


#' Get a dataframe with p-site offsets from a datafram eof KL divergence
#'
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param frustprofilelist
#'
#' @return a dataframe with p-site offsets per readlength

select_offsets <- function(kl_df) {
  if (method == "a_max") {
    kl_offsets <- kl_df %>%
      group_by(nreadlen, position) %>%
      summarise(sumKL = sum(KL)) %>%
      filter(position < -6) %>%
      filter(position > -(nreadlen - 6)) %>%
      slice(which.max(sumKL)) %>%
      mutate(p_offset = position + 3)
  }
  kl_offsets %>% select(nreadlen, position, p_offset)
}



#' Plot KL divergence in selected windows
#'
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param kl_df  dataframe with KL divergence per readlength, position in window
#'  around codon
#' @param kl_offsets  dataframe with information on p-sites
#' @param selreadlens  integer vector of readlengths to plot
#'
#' @return a dataframe with p-site offsets per readlength


plot_kl_dv <- function(kl_df, kl_offsets, selreadlens = NULL) {
  stopifnot(c("position", "nreadlen", "KL", "sample") %in% colnames(kl_df))
  stopifnot(c("nreadlen", "sample", "p_offset") %in% colnames(kl_offsets))
  #
  kl_offsets2plot <- kl_offsets %>% filter(nreadlen %in% selreadlens)
  #
  kl_df %>%
    filter(position < -3) %>%
    filter(position > -(nreadlen - 6)) %>%
    # separate(sample,c('fraction','genotype','rep'),remove=F)%>%
    filter(nreadlen %in% selreadlens) %>%
    {
      qplot(data = ., x = position, y = KL) +
        theme_bw() +
        facet_grid(nreadlen ~ sample) +
        scale_y_continuous("RUST KL divergence") +
        scale_x_continuous("5 read position relative to codon ") +
        geom_vline(
          data = kl_offsets2plot, aes(xintercept = p_offset - 3),
          color = I("green"), linetype = 2
        ) +
        geom_vline(
          data = kl_offsets2plot, aes(xintercept = p_offset),
          color = I("blue"), linetype = 2
        ) +
        ggtitle("RUST KL divergence vs position")
    }
}
