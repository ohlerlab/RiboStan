#' @import GenomicFeatures
#' @import GenomicFiles
NULL

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
  gr$phase <- unlist((start(gr) - cdsstarts) %% 3)
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
    Filter(f=function(x)x > 0.025) %>%
    {
      Filter(f=function(x)  1 - x > 0.025,.)
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
  bparam <- ScanBamParam(simpleCigar = TRUE, scanBamFlag(
    isUnmappedQuery = FALSE,
    isMinusStrand = FALSE
  ))
  ribogr <- GenomicAlignments::readGAlignments(
    ribobam, use.names = T, param = bparam)
  #
  readlens <- get_readlens(qwidth(ribogr))
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
  genomicbam = seqnames(cov) %>% head(1) %>% str_detect("chr")
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
    unlist%>%
    intersect(unique(seqnames(cov)))
  cov <- cov %>% keepSeqlevels(sharedseqnames, pruning = "coarse")
  seqlevels(cov) <- seqlevels(trspacecds)
  seqinfo(cov) <- seqinfo(trspacecds)
  cov <- cov %>% subsetByOverlaps(trspacecds)
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

mergeseqlevels <- function(gr1,gr2){
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
  mergeseqinfo <- as(mergeseqinfo, 'Seqinfo')
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

get_psite_gr <- function(rpfs, offsets_df, anno) {
  orfs <- c(anno$trspacecds)
  #
  rpfs <- mergeseqlevels(rpfs, orfs)
  orfs <- mergeseqlevels(orfs, rpfs)
  #
  orfov_ind <- findOverlaps(rpfs,orfs, select='first', ignore.strand=T)
  rpfs$orf <- Rle(names(orfs)[orfov_ind])
  rpfs$phase <- start(rpfs) - start(orfs)[orfov_ind]
  rpfs$phase <- rpfs$phase%%3
  #
  orfov <- countOverlaps(rpfs,orfs, ignore.strand=T)
  is_mov_rpf <- orfov>1
  #get phase for multi orf alignments
  mov_rpfs <- rpfs[is_mov_rpf]
  mov_rpfs_ov <- orfov[is_mov_rpf]
  mov_rpfs_ind <- findOverlaps(mov_rpfs,orfs, ignore.strand=T)
  mov_rpfs <- rep(mov_rpfs, mov_rpfs_ov)
  mov_rpfs$orf <- names(orfs)[subjectHits(mov_rpfs_ind)]
  mov_rpfs$phase <- start(mov_rpfs) - start(orfs)[subjectHits(mov_rpfs_ind)]
  mov_rpfs$phase <- mov_rpfs$phase %%3
  offsetcols <- c('readlen', 'phase', 'p_offset')
  #get offsets for our multi orf alignments
  mov_rpfs$p_offset <- mov_rpfs%>%{
        as.data.frame(mcols(.)[, c("readlen", "phase")])
      } %>%
      left_join(
        offsets_df %>% select(one_of(offsetcols)),
        by = c("readlen", "phase")
      ) %>%
      .$p_offset
  #
  #for the remaining ambigous ones choose randomly
  uniq_mov_rpfs<-mov_rpfs%>%
    resize(1, 'start')%>%
    shift(.$p_offset)%>%
    subsetByOverlaps(orfs)%>%
    sample%>%
    {.[match(unique(names(.)),names(.))]}%>%
    .[order(as.numeric(names(.)))]%>%
    shift(-.$offset)%>%
    sort
  #add in phase to the ambiguous ones
  uniq_mov_rpfs$p_offset<-NULL
  rpfs[names(uniq_mov_rpfs)]<-uniq_mov_rpfs
  #now get offsets for all
  rpfs$p_offset <- rpfs%>%{
        as.data.frame(mcols(.)[, c("readlen", "phase")])
      } %>%
      left_join(
        offsets_df %>% select(one_of(offsetcols)),
        by = c("readlen", "phase")
      ) %>%
      .$p_offset
  psites <- rpfs%>%
    subset(!is.na(p_offset))%>%
    resize(1, 'start')%>%
    shift(.,.$p_offset)
  
  #phaseshift the psites

  rl_phs_ <- mcols(psites)[,c('phase','readlen')]%>%
    as.data.frame%>%
    group_by(phase,readlen)%>%tally
  shiftdf <- rl_phs_%>%group_by(readlen)%>%
    mutate(shft = rank(-n)-1)%>%as.data.frame

  phaseshifts <- merge(
    mcols(psites)[,c('phase','readlen')],
    shiftdf, all.x=TRUE)
  #
  psites <- psites%>%
    GenomicRanges::shift(-.$phase)%>%
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
#' \dontrun{
#'
#' }
#' @export

get_readgr <- function(ribobam, anno, offsets_df=NULL, startstop=FALSE, strip_seqnames = TRUE) {
  bamseqnames <- seqinfo(Rsamtools::BamFile(ribobam))@seqnames
  if (strip_seqnames) bamseqnames <- str_replace(bamseqnames, "\\|.*", "")
  stopifnot(
    mean(unlist(unique(seqnames(anno$trspacecds))) %in% bamseqnames) > .9
  )
  if(startstop){
    seltrs = intersect(names(anno$trspacecds),bamseqnames)
    which <- c(anno$trspacecds[seltrs]%>%resize(1, 'start'),
      anno$trspacecds[seltrs]%>%resize(1, 'end'))
    strand(which)<-'+'
  }else{
    which<-NULL
  }
  cov <- read_ribobam(ribobam, which)
  if(!is.null(offsets_df)){
    cov <- get_psite_gr(cov, offsets_df, anno)
  }else{
    cov <- get_cds_reads(cov, anno)    
  }
  cov
}


#' get_read_spmat
#'
#' get optimal ritpms using stan
#'
#' @param cov GRanges object of
#' @return A matrix of the infile

get_read_spmat <- function(cov, anno) {
  stopifnot(length(cov)>0) 
  orfs <- c(anno$trspacecds)
  #
  cov <- mergeseqlevels(cov, orfs)
  orfs <- mergeseqlevels(orfs, cov)
  orfs <- subsetByOverlaps(orfs, cov)
  #
  spmat <- Matrix::sparseMatrix(
    i = names(cov)%>%id,
    j = cov$orf%>%id,
    x = 1
  )
  colnames(spmat) <- cov$orf%>%unique
  rownames(spmat) <- names(cov)%>%unique
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

optimize_ritpms <- function(spmat, anno, iternum = 500, verbose=FALSE) {
  trlens <- anno$trspacecds %>%
    width() %>%
    unlist() %>%
    setNames(names(anno$trspacecds))
  # now let's try the whole shebang in rstan
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
#' @return A matrix of the


sample_cols_spmat <- function(spmat, return_mat = T) {
  #
  # so I think p gives us element which is the final one for each
  matsum <- summary(spmat)
  cs <- cumsum(spmat@x)
  colsums <- colSums(spmat)
  pts <- spmat@p
  prevcs <- c(0, cs[pts %>% tail(-1)])[matsum$j]
  cs <- cs - prevcs
  xold <- spmat@x
  spmat@x <- cs
  spmat
  spmatnm <- t(t(spmat) / colsums)
  nvals <- spmat@x %>% length()
  passrand <- runif(nvals) < spmatnm@x
  # firstpass = diff(c(0,passrand))==1
  # logmat = spmat>0
  # logmat@x = passrand
  # (spmat*logmat)%>%summary%>%subset(x!=0)%>%nrow
  #
  spmat@x <- xold
  spmat@x <- spmat@x * passrand
  pass_summ <- summary(spmat) %>% as.data.frame()
  # pass_summ$ind = 1:nrow(pass_summ)
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
#' @param cov A GRanges object containing RPFs
#' @param anno An annotation object with a gene-transcript table
#' @return a vector of normalized footprint densities
#' @export


get_ritpms <- function(ritpm_opt) {
  #quantify_orfs
  spmat <- get_read_spmat(cov, anno)
  #
  ritpm_opt <- optimize_ritpms(spmat, anno, iternum=100)
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
#' @param cov A GRanges object containing RPFs
#' @param anno An annotation object with a gene-transcript table
#' @param ritpms a vector of ribosoome densities
#' @return a vector of normalized footprint densities

sample_cov_gr <- function(cov, anno, ritpms) {
  spmat <- get_read_spmat(cov, anno)
  #
  matsample <- sample_cols_spmat(t(spmat) * ritpms, return_mat = F)

  #
  iddf <- tibble(
    rind = 1:length(cov), 
    j = names(cov) %>% id(),
    i = as.numeric(id(cov$orf))
  )
  #
  rinds <- iddf %>%
    inner_join(matsample, by = c("j", "i")) %>%
    .$rind
  #
  sampcov <- cov[rinds]
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
#' @export

gene_level_expr <- function(ripms, anno) {
  trgiddf <- anno$trgiddf
  trgiddf <- trgiddf%>%subset(!uORF)
  #
  gn_expr <- left_join(
    trgiddf,
    tibble::enframe(ripms, "orf_id", "ritpm"),
    by = "orf_id"
  )
  #
  gn_expr <- gn_expr %>% group_by(gene_id) %>%
    summarise(expr = sum(replace_na(ritpm, 0)))
  gn_expr %>% select(gene_id, expr)
}

#' This creates a minimal annotation object from a gencode style fasta
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param ribofasta A gencode style fasta to which RPFs were aligned
#' @return An annotation object with cdsstarts/stops
#' @details The Ribosome densities are saved in salmon format
#'
#' @export

get_ribofasta_anno <- function(ribofasta) {
  Rsamtools::indexFa(ribofasta)
  faheadernames <- seqinfo(Rsamtools::FaFile(ribofasta))
  faheaddf <- seqnames(faheadernames) %>%
    as.vector() %>%
    str_split_fixed("\\|", 10)
  anno <- list()
  anno$trspacecds <- GRanges(
    faheaddf[, 1],
    faheaddf[, 9] %>% str_extract("\\d+\\-\\d+") %>%
      str_split_fixed("-", 2) %>% set_colnames(c("start", "end")) %>%
      apply(2, as.numeric) %>% as.data.frame() %>%
      {
        IRanges(start = .$start, end = .$end)
      }
  ) %>%
    setNames(., as.character(seqnames(.)))
  strand(anno$trspacecds)<-'+'
  anno$trgiddf <- tibble(transcript_id = faheaddf[, 1], gene_id = faheaddf[, 2])
  anno$trgiddf$orf_id <- anno$trgiddf$transcript_id
  anno <- c(
    anno,
    list(
      cdsstarts = anno$trspacecds %>% start() %>%
        setNames(names(anno$trspacecds)),
      cds_prestop_st = anno$trspacecds %>% end() %>% `-`(2) %>%
        setNames(names(anno$trspacecds))
    )
  )
  anno
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

# TODOS remove 'ritpm from the package'

get_exprfile <- function(ribobam, ribofasta, outfile) {
  #
  anno <- get_ribofasta_anno(ribofasta)
  #
  cov <- get_readgr(ribobam, anno)
  #
  spmat <- get_read_spmat(cov, anno)
  #
  ritpm_opt <- optimize_ritpms(spmat, anno, iternum = 100)
  #
  ritpms <- get_ritpms(ritpm_opt)
  #
  n_reads <- nrow(spmat)
  counts <- n_reads * (ritpms / 1e6)
  counts <- tibble::enframe(counts, "Name", "NumReads")
  ritpmdf <- tibble::enframe(ritpms, "Name", "ritpm")
  #
  cdslens <- anno$trspacecds %>%
    width() %>%
    setNames(names(anno$trspacecds)) %>%
    tibble::enframe("Name", "Length") %>%
    mutate(EffectiveLength = Length)
  output <- cdslens %>%
    left_join(ritpmdf) %>%
    left_join(counts)
  #
  output <- output %>% select(Name, Length, EffectiveLength, ritpm, NumReads)
  output %>% write_tsv(outfile)
}
