id <- function(cov) BiocGenerics::match(cov, unique(cov))


#' Add 'phase' column to a set of reads, given a vector of cds starts
#'
#' Given rpfs, and named vectors of starts and stops, this adds the 'phase'
#' to the reads
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param rpfs GRanges object with RPFs
#' @param cdsstarts named vector of cds starts
#'
#' @return the rpf granges but wth a 'phase' column

# TODO update this so it works with uORFs

addphase <- function(gr, cdsstarts) {
  cdsstarts <- cdsstarts[as.vector(seqnames(gr))]
  gr$phase <- unlist((GenomicRanges::start(gr) - cdsstarts) %% 3)
  gr
}

#' Select those read lengths which make up 95% of the reads
#'
#' This function takes in a vector of the lengths of each read, and determines
#' which read lengths should be included in order to filter out the 2.5%
#' shortest/longest reads.
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param readlens Numeric; numeric vector
#' @return a numeric integer vector, such that selecting reads with these
#' lengths will preserve at least 95%
#' of reads

get_readlens <- function(readlens) {
  readlens %>%
    as.numeric() %>%
    table() %>%
    cumsum() %>%
    {
      . / max(.)
    } %>%
    Filter(f = function(x) x > 0.025) %>%
    {
      Filter(f = function(x) 1 - x > 0.025, .)
    } %>%
    names() %>%
    as.numeric()
}


#' Read a bam file containing Ribosomal Footprints
#'
#' This function
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param ribobam String; full path to html report file.
#' @param strip_seqnames  whether the function should remove all text after the
#' first '|', useful if aligning to gencode fastas#
#' Defaults to \code{TRUE}
#' @return This function returns a granges object with the read names in the
#' names slot and a metadata column
#' denoting readlength.
#'
#' @details This function reads in a bam file containing ribosomal footprints.
#' It uses only reads without splice sites and reads which align to the positive
#'  strand, as it's designed to work on transcriptomic alignments.
#'
#' @seealso \code{\link{get_cds_reads}}, \code{\link{get_readlens}}

read_ribobam <- function(ribobam, which, strip_seqnames = TRUE) {
  #
  flags <- Rsamtools::scanBamFlag(
    isUnmappedQuery = FALSE,
    isMinusStrand = FALSE
  )
  bparam <- Rsamtools::ScanBamParam(simpleCigar = TRUE, flags)
  ribogr <- GenomicAlignments::readGAlignments(
    ribobam,
    use.names = TRUE, param = bparam
  )
  #
  readlens <- get_readlens(GenomicAlignments::qwidth(ribogr))
  wfilt <- GenomicAlignments::qwidth(ribogr) <= max(readlens) &
    (min(readlens) <= GenomicAlignments::qwidth(ribogr))
  ribogr <- ribogr[wfilt]
  # this is the number get_bamdf gets in py 4628644
  mcols(ribogr)$readlen <- GenomicAlignments::qwidth(ribogr)
  ribogr <- as(ribogr, "GenomicRanges")
  # strip seqnames for when bam files contains a large fasta header
  if (strip_seqnames) {
    seqlevels(ribogr) <- str_replace(seqlevels(ribogr), "\\|.*", "")
  }
  # name the reads with integers as they appear in the sorted object
  ribogr <- sort(ribogr)
  names(ribogr) <- id(names(ribogr))
  ribogr
}

#' Filter an RPF GR for overlap with coding sequences
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param cov GRanges; A granges object with
#' @param strip_seqnames  whether the function should remove all text after the
#' first '|', useful if aligning to gencode fastas Defaults to \code{TRUE}
#' @return This function returns a granges object with the read names in the
#' names slot and a metadata column
#' denoting readlength.
#'
#' @details This function reads in a bam file containing ribosomal footprints.
#' It uses only reads without splice sites
#' and reads which align to the positive strand, as it's designed to work on
#' transcriptomic alignments.
#'
#' @seealso \code{\link{get_cds_reads}}, \code{\link{get_readlens}}

get_cds_reads <- function(cov, anno) {
  trspacecds <- anno$trspacecds

  # load genomic or transcriptomic bam
  genomicbam <- seqnames(cov) %>%
    head(1) %>%
    as.vector() %>%
    str_detect("chr")
  if (genomicbam) {
    cov <- cov %>%
      resize(1) %>%
      GenomicFeatures::mapToTranscripts()
    cov$readlen <- mcols(cov)$readlen[cov$xHits]
    cov$name <- names(cov)[cov$xHits]
    cov <- subset(cov, between(readlen, min(readlens), max(readlens)))
    cov <- resize(cov, 1, "start")
    cov <- sort(cov)
  } else {
    cov <- cov
    cov <- sort(cov)
  }

  sharedseqnames <- unique(seqnames(trspacecds)) %>%
    unlist() %>%
    intersect(unique(seqnames(cov)))
  cov <- cov %>% keepSeqlevels(sharedseqnames, pruning.mode = "coarse")
  seqlevels(cov) <- seqlevels(trspacecds)
  seqinfo(cov) <- seqinfo(trspacecds)
  cov <- cov %>% IRanges::subsetByOverlaps(trspacecds)
}



#' Merge the seqlevels of two Granges objects
#' This function checks two GRanges objects have compatible ranges, and then
#' adds in to gr1 what's missing from gr2.
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param gr1 - granges object
#' @param gr2 - granges object
#' @return gr1 but with expanded seqinfo

mergeseqlevels <- function(gr1, gr2) {
  unseqs <- union(seqlevels(gr1), seqlevels(gr2))
  shared <- intersect(seqlevels(gr1), seqlevels(gr2))
  s1unique <- setdiff(seqlevels(gr1), shared)
  stopifnot(
    all(seqinfo(gr1)[shared]@seqlengths == seqinfo(gr2)[shared]@seqlengths)
  )
  mergeseqinfo <- bind_rows(
    as.data.frame(seqinfo(gr1)[s1unique]),
    as.data.frame(seqinfo(gr2))
  )
  mergeseqinfo <- as(mergeseqinfo, "Seqinfo")
  seqlevels(gr1) <- seqlevels(mergeseqinfo)
  seqinfo(gr1) <- mergeseqinfo
  gr1
}


#' convert a granges object of footprints to a psites object
#'
#' This function applies offsets to RPF data, attributing to each RPF
#' an offset specific a readlength and phase. Where phase is ambigous,
#' because more than one ORF overlaps the RPF, the RPF has both offsets
#' applied, and in rare cases where more than one possible psite location
#' applies, one is randomly selected.
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param rpfs - granges object with positions of footprints, in transcript
#' space
#' @param offsets_df - a data frame with numeric columns readlen,phase,offset
#' @param anno - an annotation object with at least one ORF per transcript,
#' the function uses the 'uORF' slot if it's available
#' @return a GRanges object of width-1 psites with columns readlen,
#' phase and offset
#' @examples
#' data(chr22_anno)
#' data(rpfs)
#' data(offsets_df)
#' psites <- get_psite_gr(rpfs, offsets_df, anno)
#' @export
get_psite_gr <- function(rpfs, offsets_df, anno) {
  orfs <- c(anno$trspacecds)
  #
  rpfs <- mergeseqlevels(rpfs, orfs)
  orfs <- mergeseqlevels(orfs, rpfs)
  #
  orfov_ind <- findOverlaps(rpfs, orfs, select = "first", ignore.strand = TRUE)
  rpfs$orf <- S4Vectors::Rle(names(orfs)[orfov_ind])
  rpfs$phase <- GenomicRanges::start(rpfs) - GenomicRanges::start(orfs)[orfov_ind]
  rpfs$phase <- rpfs$phase %% 3
  #
  # now get offsets for all
  offsetcols <- c("readlen", "phase", "p_offset")
  rpfs$p_offset <- rpfs %>%
    {
      as.data.frame(mcols(.)[, c("readlen", "phase")])
    } %>%
    left_join(
      offsets_df %>% select(one_of(offsetcols)),
      by = c("readlen", "phase")
    ) %>%
    .$p_offset
  rpfs <- rpfs %>% subset(!is.na(p_offset))
  orfov <- GenomicRanges::countOverlaps(rpfs, orfs, ignore.strand = TRUE)
  is_mov_rpf <- orfov > 1
  # get phase for multi orf alignments
  mov_rpfs <- rpfs[is_mov_rpf]
  mov_rpfs_ov <- orfov[is_mov_rpf]
  mov_rpfs_ind <- findOverlaps(mov_rpfs, orfs, ignore.strand = TRUE)
  mov_rpfs <- rep(mov_rpfs, mov_rpfs_ov)
  mov_rpfs$orf <- names(orfs)[subjectHits(mov_rpfs_ind)]
  mov_rpfs$phase <- GenomicRanges::start(mov_rpfs) - GenomicRanges::start(orfs)[subjectHits(mov_rpfs_ind)]
  mov_rpfs$phase <- mov_rpfs$phase %% 3
  # get offsets for our multi orf alignments
  mov_rpfs$p_offset <- mov_rpfs %>%
    {
      as.data.frame(mcols(.)[, c("readlen", "phase")])
    } %>%
    left_join(
      offsets_df %>% select(one_of(offsetcols)),
      by = c("readlen", "phase")
    ) %>%
    .$p_offset
  #
  # for the remaining ambigous ones choose randomly
  uniq_mov_rpfs <- mov_rpfs %>%
    resize(1, "start") %>%
    shift(.$p_offset) %>%
    IRanges::subsetByOverlaps(orfs) %>%
    sample() %>%
    {
      .[match(unique(names(.)), names(.))]
    } %>%
    .[order(as.numeric(names(.)))] %>%
    shift(-.$p_offset) %>%
    sort()
  # add in phase to the ambiguous ones
  uniq_mov_rpfs$p_offset <- NULL
  rpfs[names(uniq_mov_rpfs)] <- uniq_mov_rpfs

  psites <- rpfs %>%
    subset(!is.na(p_offset)) %>%
    resize(1, "start") %>%
    shift(., .$p_offset)
  psites <- psites[!is_out_of_bounds(psites)]
  # phaseshift the psites

  rl_phs_ <- mcols(psites)[, c("phase", "readlen")] %>%
    as.data.frame() %>%
    group_by(.data$phase, .data$readlen) %>%
    tally()
  shiftdf <- rl_phs_ %>%
    group_by(.data$readlen) %>%
    mutate(shft = rank(-.data$n) - 1) %>%
    as.data.frame()

  phaseshifts <- merge(
    mcols(psites)[, c("phase", "readlen")],
    shiftdf,
    all.x = TRUE
  )
  #
  psites <- psites %>%
    GenomicRanges::shift(-.$phase) %>%
    GenomicRanges::shift(phaseshifts$shft)
  psites
}






#' Read and filter a bam file of RPFs
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param ribobam GRanges; A bam file with RPFs
#' @param anno An annotation object
#' @param startstop option to select only reads that overlapt he start/stop
#' for offset determination,Defaults to FALSE
#' @param strip_seqnames  whether the function should remove all text after the
#' first '|', useful if aligning to gencode fastas Defaults to \code{TRUE}
#' @return This function returns a granges object with the read names in the
#' names slot and a metadata column
#' denoting readlength.
#'
#' @details This function reads in a bam file containing ribosomal footprints.
#' It uses only reads without splice sites and reads which align to the
#' positive strand, as it's designed to work on transcriptomic alignments.
#'
#' @seealso \code{\link{get_cds_reads}}, \code{\link{get_readlens}}
#' @examples
#' data(chr22_anno)
#' testbam <- system.file("extdata", "chr22.bam", package = "Ribostan", mustWork = TRUE)
#' rpfs <- get_readgr(testbam, chr22_anno)
#' @export

get_readgr <- function(ribobam, anno, offsets_df = NULL, startstop = FALSE, strip_seqnames = TRUE) {
  bamseqnames <- seqinfo(Rsamtools::BamFile(ribobam))@seqnames
  if (strip_seqnames) bamseqnames <- str_replace(bamseqnames, "\\|.*", "")
  orf_trs <- unlist(unique(seqnames(anno$trspacecds)))
  stopifnot(
    mean(orf_trs %in% bamseqnames) > .5
  )
  if (startstop) {
    seltrs <- intersect(names(anno$trspacecds), bamseqnames)
    which <- c(
      anno$trspacecds[seltrs] %>% resize(1, "start"),
      anno$trspacecds[seltrs] %>% resize(1, "end")
    )
    strand(which) <- "+"
  } else {
    which <- NULL
  }
  cov <- read_ribobam(ribobam, which)
  if (!is.null(offsets_df)) {
    cov <- get_psite_gr(cov, offsets_df, anno)
  } else {
    cov <- get_cds_reads(cov, anno)
  }
  cov
}


#' get_read_spmat
#'
#' get optimal ritpms using stan
#'
#' @param psites GRanges object of
#' @return A matrix of the infile

get_read_spmat <- function(psites, anno) {
  stopifnot(length(psites) > 0)
  orfs <- c(anno$trspacecds)
  #
  psites <- mergeseqlevels(psites, orfs)
  orfs <- mergeseqlevels(orfs, psites)
  orfs <- IRanges::subsetByOverlaps(orfs, psites)
  #
  spmat <- Matrix::sparseMatrix(
    i = names(psites) %>% id(),
    j = psites$orf %>% id(),
    x = 1
  )
  colnames(spmat) <- psites$orf %>% unique()
  rownames(spmat) <- names(psites) %>% unique()
  spmat <- spmat / Matrix::rowSums(spmat)
  spmat
}



#' Use a sparse mapping matrix to optimize ritpms
#'
#' get optimal ritpms using stan
#'
#' @param spmat a sparse numeric matrix
#' @param anno an annotation object
#' @param iternum how many iterations to optimize TPMs for;
#' Defaults to 500
#' @param verbose whether to show optimization messages from stan;
#' Defaults to FALSE
#' @return A matrix of the infile

optimize_ritpms <- function(spmat, anno, iternum = 500, verbose = FALSE) {
  trlens <- anno$trspacecds %>%
  width() %>%
  unlist() %>%
  setNames(names(anno$trspacecds))
  # now let's try the whole shebang in rstan
  setdiff(colnames(spmat),names(trlens))
  stopifnot(colnames(spmat)%in%names(trlens))
  sptrlens <- trlens[colnames(spmat)]
  fdata <- list(nonorm_trlen = sptrlens)
  fdata <- c(fdata, spmat %>% rstan::extract_sparse_parts(.))
  fdata$trlen <- fdata$nonorm_trlen %>%
  {
    . / sum(.)
  }
  fdata$TR <- spmat %>% ncol()
  fdata$R <- spmat %>% nrow()
  fdata$V <- fdata$w %>% length()
  fdata$Ulen <- fdata$u %>% length()
  fdata$classweights <- rep(1, fdata$R)
  init <- list(ritpm = spmat %>%
  {
    Matrix::colSums(.)
    } %>% `/`(fdata$trlen) %>%
    {
      . / sum(.)
      })
  #
  modelcode <- "
  data {
    int TR;// number of TRs
    int R;// number of reads
    int V;
    int Ulen;
    vector [V] w ;
    int  v [V];
    int   u [Ulen];
    vector [R] classweights;
  }
  parameters {
    simplex [TR] n;
  }
  model {
    target += log(
    csr_matrix_times_vector(R, TR, w, v, u, n)
    ).* classweights;
  }
  "
  eqritpm_mod <- rstan::stan_model(model_code = modelcode)
  message('optimizing...')
  opt <- rstan::optimizing(
    eqritpm_mod,
    data = fdata,
    init = init,
    verbose = verbose,
    iter = iternum
    )
  opt$seqnames <- colnames(spmat)
  opt$trlen <- fdata$nonorm_trlen[opt$seqnames]
  opt
}



#' sample_cols_spmat
#'
#' This function takes in a sparse matrix, and then samples a single value
#' from each column, with sampling weights within each column equal to
#' the columns values
#'
#' @param spmat a sparse numeric matrix describing RPF's multimapping
#' @param return_mat whether to a matrix for selecting the elements, 
#'                  as opposed to the summary data frame;defaults to TRUE
#' @return A matrix of the


sample_cols_spmat <- function(spmat, return_mat = TRUE) {
  #
  # so I think p gives us element which is the final one for each
  matsum <- Matrix::summary(spmat)
  cs <- cumsum(spmat@x)
  colsums <- Matrix::colSums(spmat)
  pts <- spmat@p
  prevcs <- c(0, cs[pts %>% tail(-1)])[matsum$j]
  cs <- cs - prevcs
  xold <- spmat@x
  spmat@x <- cs
  spmatnm <- Matrix::t(Matrix::t(spmat) / colsums)
  nvals <- spmat@x %>% length()
  passrand <- runif(nvals) < spmatnm@x
  spmat@x <- xold
  spmat@x <- spmat@x * passrand
  pass_summ <- Matrix::summary(spmat) %>% as.data.frame()
  pass_summ <- pass_summ[pass_summ$x != 0, ]
  pass_summ <- pass_summ[diff(c(0, pass_summ$j)) > 0, ]
  # returnmatrix or the summary
  if (return_mat) {
    outmat <- Matrix::sparseMatrix(
      i = pass_summ$i, j = pass_summ$j, x = pass_summ$x,
      dims = dim(spmat)
    )
    stopifnot(dim(outmat) == dim(spmat))
    outmat
  } else {
    pass_summ
  }
}

################################################################################

#' get_ritpms
#'
#' This function takes in a sparse matrix, and then samples a single value
#' from each column, with sampling weights within each column equal to
#' the columns values
#' @param psites A GRanges object containing RPFs
#' @param anno An annotation object with a gene-transcript table
#' @return a vector of normalized footprint densities
#' @examples
#' data(chr22_anno)
#' data(rpfs)
#' data(offsets_df)
#' data(ms_df)
#'
#' psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
#' # now, use Stan to estimate normalized p-site densities for our data
#' ritpms <- get_ritpms(psites, chr22_anno)
#' @export

get_ritpms <- function(psites, anno) {
  psites <- psites%>%subset(orf %in% names(anno$trspacecds))
  # quantify_orfs
  spmat <- get_read_spmat(psites, anno)
  #
  ritpm_opt <- optimize_ritpms(spmat, anno, iternum = 100)
  #
  if (ritpm_opt$par %>% names() %>% str_detect("^ritpm\\[\\d+\\]$") %>% any()) {
    ritpmpars <- ritpm_opt$par %>%
      names() %>%
      str_subset("^ritpm\\[\\d+\\]$")
    ritpms <- ritpm_opt$par[ritpmpars]
    names(ritpms) <- ritpm_opt$seqnames
    ritpms <- ritpms * 1e6
    ritpms
  } else {
    ritpmpars <- ritpm_opt$par %>%
      names() %>%
      str_subset("^n\\[\\d+\\]$")
    ritpms <- ritpm_opt$par[ritpmpars]
    names(ritpms) <- ritpm_opt$seqnames
    ritpms <- ritpms / ritpm_opt$trlen
    ritpms <- ritpms / sum(ritpms)
    ritpms <- ritpms * 1e6
    ritpms
  }
}

#' sample_cov_gr
#'
#' This function takes in a coverage GR and takes one out of each multimap
#' weighting according to the ritpms
#'
#' @param psites A GRanges object containing psites
#' @param anno An annotation object with a gene-transcript table
#' @param ritpms a vector of ribosoome densities
#' @return a vector of normalized footprint densities
#' @export

sample_cov_gr <- function(psites, anno, ritpms) {
  spmat <- get_read_spmat(psites, anno)
  #

  spmat <- Matrix::t(spmat) * ritpms[colnames(spmat)]
  matsample <- sample_cols_spmat(spmat,return_mat = FALSE)

  #
  iddf <- tibble(
    rind = seq_along(psites),
    j = names(psites) %>% id(),
    i = as.numeric(id(psites$orf))
  )
  #
  rinds <- iddf %>%
    inner_join(matsample, by = c("j", "i")) %>%
    .$rind
  #
  sampcov <- psites[rinds]
  #
  sampcov
}


#' Aggregate transcript level RiboPMs to the gene-level
#' Ignores uORF expression
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param ripms - named vector of ribosome densities
#' @param anno An annotation object with a gene-transcript table
#' @return a data frame with columns gene_id, expr
#' @details
#'
#' @examples
#' data(chr22_anno)
#' data(rpfs)
#' data(offsets_df)
#' data(ms_df)
#'
#' psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
#' # now, use Stan to estimate normalized p-site densities for our data
#' ritpms <- get_ritpms(psites, chr22_anno)
#' # and get these at the gene level (ignoring uORFs)
#' gritpms <- gene_level_expr(ritpms, chr22_anno)
#' @export

gene_level_expr <- function(ripms, anno) {
  trgiddf <- anno$trgiddf
  if('uORF'%in%colnames(trgiddf)){
    trgiddf <- trgiddf %>% subset(!uORF)
  }
  #
  gn_expr <- left_join(
    trgiddf,
    tibble::enframe(ripms, "orf_id", "ritpm"),
    by = "orf_id"
  )
  #
  gn_expr <- gn_expr %>%
    group_by(.data$gene_id) %>%
    summarise(expr = sum(replace_na(.data$ritpm, 0)))
  gn_expr %>% select('gene_id', 'expr')
}


#' Given a bam file, and transcript fasta, output a file of ribosome densities
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param ribobam A bam file with RPFs
#' @param ribofasta A gencode style fasta to which RPFs were aligned
#' @param outfile the file to which the RPF densitites should be saved
#' @return the file to which the table was saved
#'
#' @details The Ribosome densities are saved in salmon format
#' @export
#' @examples
#' get_exprfile(ribobam, ribofasta, outfile)


get_exprfile <- function(ribobam, ribofasta, outfile) {
  #
  anno <- get_ribofasta_anno(ribofasta)
  #
  rpfs <- get_readgr(ribobam, anno)
  # determine offsets by maximum CDS occupancy
  offsets_df <- get_offsets(rpfs, anno)
  # use our offsets to determine p-site locations
  psites <- get_psite_gr(rpfs, offsets_df, anno)
  # now, use Stan to estimate normalized p-site densities for our data
  ritpms <- get_ritpms(psites, anno)
  #
  n_reads <- n_distinct(names(rpfs))
  lengths <- width(anno$trspacecds[names(ritpms)])
  nucfracs <- (ritpms * lengths)
  nucfracs <- nucfracs / sum(nucfracs)
  counts <- n_reads * nucfracs
  counts <- tibble::enframe(counts, "Name", "NumReads")
  ritpmdf <- tibble::enframe(ritpms, "Name", "ritpm")
  #
  cdslens <- anno$trspacecds %>%
    width() %>%
    setNames(names(anno$trspacecds)) %>%
    tibble::enframe("Name", "Length") %>%
    mutate(EffectiveLength = .data$Length)
  output <- cdslens %>%
    left_join(ritpmdf) %>%
    left_join(counts)
  #
  output <- output %>% 
    select('Name', 'Length', 'EffectiveLength', 'ritpm', 'NumReads')
  output %>% write_tsv(outfile)
}
