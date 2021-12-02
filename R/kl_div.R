#' Get A Granges object with locations of each codon in the given cds
#'
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param anno An annotation object
#' @param offsets_df - a data frame with numeric columns readlen,phase,offset
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

# TODO this function is dumb...
# TODO stop missusing the name slot, maybe move some of this out so we don't
# have so many arguments

get_cds_codons <- function(anno,
                           n_wind_l_ext = 45,
                           n_wind_r_ext = 9,
                           n_start_buff = 60,
                           n_end_buff = 60) {
	stopifnot((n_start_buff%%3)==0)
	stopifnot((n_end_buff%%3)==0)
	fafileob <- anno$fafileob
	longtrs <- anno$longtrs
	longorfs <- anno$trgiddf$orf_id[match(longtrs,anno$trgiddf$transcript_id)]
	message('getting codon positions...')
	exonsgrl <- anno$exonsgrl[longtrs]
	trspacecds <- anno$trspacecds[longorfs]
	names(trspacecds)<-longtrs
	#sequence of relevant exons
	exonseq <- exonsgrl %>% 
		sort_grl_st%>%
		GenomicFeatures::extractTranscriptSeqs(x = fafileob)
	#data frame with codon and position info
	codposdf <- get_codposdf(longorfs, anno)
	#as a gr
	codposdf <- GRanges(codposdf$orf_id,IRanges(codposdf$pos),codon=codposdf$codon)
	codongr <- codposdf %>% GenomicFeatures::mapFromTranscripts(trspacecds)
	codongr$codon <- codposdf$codon[codongr$xHits]
	codongr<-resize(codongr, 3, 'start')
	#add in seqinfo
	seqlevels(codongr)<-names(exonseq)
	seqinfo(codongr)<-Seqinfo(names(exonseq),nchar(exonseq))
	strand(codongr)<-'+'
	# expand the windows around these codons
  codongrbak <- codongr
	codongr <- codongr[(start(codongr)-1)>=n_wind_l_ext]
  enddists <- (seqlengths(codongr)[as.character(seqnames(codongr))]-
    end(codongr))
  codongr <- codongr[enddists>=n_wind_r_ext]
  codmatchwindows <- codongr %>%
  	resize(n_wind_l_ext, "end") %>%
  	resize(width(.)+n_wind_r_ext, "start")
	# require these to be inside our cds
  # use only matches in the inner cds
  innercds <- trspacecds %>%
    subset(width > (3 + n_start_buff + n_end_buff)) %>%
    resize(width(.) - n_start_buff, "end") %>%
    resize(width(.) - n_end_buff, "start")
  codongr <- codongr %>% subsetByOverlaps(innercds)
	codmatchwindows <- codmatchwindows[!is_out_of_bounds(codmatchwindows)]
	# define inner cds - cds but with start and end trimmed off.
	codmatchwindows <- codmatchwindows%>%subsetByOverlaps(innercds)
	# lift them to inner cds coordinates
	cds_codons <- codmatchwindows %>% 
		GenomicFeatures::mapToTranscripts(innercds)
	cds_codons$codon <- codmatchwindows$codon[cds_codons$xHits]
	cds_codons$xHits <- NULL
	cds_codons$transcriptsHits <- NULL
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
  fp_profiles <- rpf_cov %>% lapply(function(cdsrlfpcov) {
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
      lapply(colMeans)
    # also get evals
    tr_rust_evals <- rustcdsrlfpcov %>% mean()
    re_c <- tr_rust_evals[codtrs] %>%
      split(cds_codons_nz$codon) %>%
      vapply(mean,1.0)
    ro_cl <- ro_cl %>% lapply(tibble::enframe, 
                              "position", "ro_cl")
    ro_cl <- bind_rows(ro_cl, .id = "codon")
    re_c <- tibble::enframe(re_c, "codon", "re_c")
    ro_cl %>% left_join(re_c, by = "codon")
  })
  fp_profiles <- bind_rows(.id = "readlen", fp_profiles)
  fp_profiles
}

#' Get A Granges object with locations of each codon in the given cds
#'
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param rpf_covs A list of SimpleRleLists - RPF coverage 
#' @param anno - an annotation object
#' scores
#'
#' @details gets profiles of read 5' end occurence around codons
#' @return A data frame with expected and actual rust score per position in
#' window, codon, readlength

get_metacodon_profs <- function(covgrs, anno, n_wind_l_ext = 45) {
  are_psites <- covgrs%>%vapply(function(x) 'orf'%in%colnames(mcols(x)),TRUE)
  stopifnot(!any(are_psites)) 
  cds_codons <- get_cds_codons(anno)
  
  stopifnot(!is.null(names(covgrs)))
  cdsfpcovlist <- lapply(covgrs, function(covgr) {
  	rlsplitcov <- split(resize(covgr, 1, "start"), covgr$readlen)
    lapply(rlsplitcov, coverage)
  })
  rust_roel <- mclapply(
    mc.cores = detectCores(), cdsfpcovlist, cds_codons = cds_codons,
    get_cov_rust_scores
  )
  # https://www.nature.com/articles/ncomms12915#Sec10 see equation 3
  # of RUST paper
  metacodondf <- bind_rows(rust_roel, .id = "sample")
  #
  metacodondf <- mutate(metacodondf, position = position - 1 - (n_wind_l_ext))
  metacodondf <- filter(metacodondf, !codon %in% c("TAG", "TAA", "TGA"))
  metacodondf <- mutate(metacodondf, count = ro_cl / re_c)
  metacodondf$nreadlen <- metacodondf$readlen %>% as.numeric()
  metacodondf
}

#' Calculate KL-divergence from a data frame of rust expected and actual scores
#'
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param metacodondf a data frame with metacodon info
#' @examples
#'   data(chr22_anno)
#'  data(rpfs)
#'  data(offsets_df)
#'  covgrs = list(sample1=rpfs)
#'  #note this doesn't work that well on a small subset
#'  kl_df<-get_kl_df(metacodondf, chr22_anno)
#' 
#' @return A data frame with expected and actual rust score per position in
#' window, codon, readlength

get_kl_df <- function(metacodondf, anno) {
  metacodondf %>%
    group_by(sample, nreadlen, position) %>%
    mutate(ro_cl = ro_cl / sum(ro_cl), re_c = re_c / sum(re_c)) %>%
    summarise(KL = sum(ro_cl * log2(ro_cl / re_c)),.groups='keep')
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
#' @param metacodondf A data frame with expected and actual rust score per position in
#' window, codon, readlength
#' @param method means by which to select best offfset currently only possible value 
#' is a_max - assumes the maximum KL divergence occurs under the A site
#'
#' @return a dataframe with p-site offsets per readlength, sample
#' @examples
#'   data(chr22_anno)
#'  data(rpfs)
#'  data(offsets_df)
#'  covgrs = list(sample1=rpfs)
#'  #note this doesn't work that well on a small subset
#'  kl_df<-get_kl_df(covgrs, chr22_anno)#' 
#'  kl_offsets <- select_offsets(kl_df)
#' @export


select_offsets <- function(kl_df, method='a_max') {
  if (method == "a_max") {
    kl_offsets <- kl_df %>%
      group_by(sample, nreadlen, position) %>%
      summarise(sumKL = sum(KL),.groups='keep') %>%
      filter(position < -6) %>%
      filter(position > -(nreadlen - 6)) %>%
      group_by(sample, nreadlen) %>%
      dplyr::slice(which.max(sumKL)) %>%
      mutate(p_offset = -(position + 3))
  }
  kl_offsets %>% ungroup%>%select(sample, nreadlen, p_offset)
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
#' @export

plot_kl_dv <- function(kl_df, kl_offsets, selreadlens = NULL) {
  stopifnot(c("position", "nreadlen", "KL", "sample") %in% colnames(kl_df))
  stopifnot(c("nreadlen", "sample", "p_offset") %in% colnames(kl_offsets))
  #
  if(!is.null(selreadlens)) {
    kl_offsets <- kl_offsets %>% filter(nreadlen %in% selreadlens)
  }
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
          data = kl_offsets, aes(xintercept = p_offset - 3),
          color = I("green"), linetype = 2
        ) +
        geom_vline(
          data = kl_offsets, aes(xintercept = p_offset),
          color = I("blue"), linetype = 2
        ) +
        ggtitle("RUST KL divergence vs position")
    }
}



#' Create a dataframe with a_site, p_site and e site occupancies
#'
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param metacodondf A data frame with expected and actual rust score per position in
#' window, codon, readlength
#' @param kl_offsets A dataframe with p-site offsets per readlength, sample
#'
#' @return a dataframe occupancies for e,p and a sites for each sample
#' @examples
#' data(chr22_anno)
#' data(rpfs)
#' data(offsets_df)
#' data(metacodondf)
#' #note this doesn't work that well on a small subset
#' kl_df<-get_kl_df(metacodondf, chr22_anno)
#' kl_offsets <- select_offsets(kl_df)
#' allcodondt <- export_codon_dts(metacodondf, kl_offsets)
#' @export

export_codon_dts <- function(metacodondf, kl_offsets){
  posseldf = bind_rows(
    kl_offsets%>%mutate(position=-p_offset+3,site='e_site')%>%select(nreadlen,position,site),
    kl_offsets%>%mutate(position=-p_offset,site='p_site')%>%select(nreadlen,position,site),
    kl_offsets%>%mutate(position=-p_offset-3,site='a_site')%>%select(nreadlen,position,site),
    kl_offsets%>%mutate(position=-p_offset-6,site='a_p3_site')%>%select(nreadlen,position,site)
  )
  #
  allcodondt = metacodondf%>%
    inner_join(posseldf, by=c('position','nreadlen'))%>%
    group_by(sample,site,codon)%>%
    summarise(rust = sum(ro_cl)/sum(re_c))%>%
    dplyr::rename('RUST_score'='rust')%>%
    tidyr::pivot_wider(names_from='site',values_from='RUST_score')
   #
  allcodondt
}