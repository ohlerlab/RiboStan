################################################################################
########## Range manipulation
################################################################################
#' @importFrom stringr str_detect str_replace str_subset str_interp 
#' @importFrom stringr str_split_fixed str_extract
#' @importFrom readr write_tsv read_tsv
# #' @import magrittr
# #' @import purrr
# # ' @import dplyr
#' @import Biostrings 
#' @import testthat
# # ' @import BiocGenerics intersect setdiff unlist table
#' @importFrom BiocGenerics intersect setdiff unlist table union mean
#' @import ggplot2
#' @importFrom GenomeInfoDb Seqinfo keepSeqlevels seqlevels seqlengths seqinfo
#' @importFrom Biostrings subseq codons oligonucleotideFrequency
#' @importFrom GenomicRanges GRanges split GenomicRangesList strand
#' @importFrom rtracklayer import export
#' @importFrom Rsamtools ScanBamParam
#' @importFrom tidyr replace_na unnest
#' @importFrom dplyr mutate select filter lead summarise tally slice %>% lag  
#' @importFrom dplyr left_join group_by ungroup tibble inner_join bind_rows
#' @importFrom dplyr distinct arrange n count between full_join one_of n_distinct 

NULL



#' Index vector for a GRanges list object with the sub-elements ordered 5' to 3'
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param grl GRangesList; a GRangesList object
#' @return Index vector for a GRanges list object with the sub-elements ordered
#'  5' to 3'

order_grl_st <- function(grl) {
  order(start(grl) * (((strand(grl) != "-") + 1) * 2 - 3))
}

#' Sort a GRanges list object with the sub-elements ordered 5' to 3'
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param grl GRangesList; a GRangesList object
#' @return A GRangeList object with the sub-elements ordered 5' to 3'

sort_grl_st <- function(grl) grl[order_grl_st(grl), ]

#' Resize a GRangesList object holding it's 5' end fixed
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param grl GRangesList; a GRangesList object
#' @param width GRangesList; integer/IntegerList to set as new width.
#' @return A GRangesList object shortend/lengthened, respecting exon boundaries

resize_grl_startfix <- function(grl, width) {
  # what follows is some slightly black magic using S4 vectors
  # Integerlist which showings how much we'd need to trim that exon to get to
  # to the desired transcript length
  trim <- cumsum(width(grl)) - width
  # Where trim is greater than the exon width, we drop it
  drop <- trim >= width(grl)
  grl <- grl[!drop]
  # vector showing location of the new 3' end of each transcript
  newends <- cumsum(elementNROWS(grl))
  # vector with the amount we need to trim each new 3' end by
  endtrims <- trim[IntegerList(as.list(elementNROWS(grl)))]@unlistData
  # finally, use these to trim
  grl@unlistData[newends] <- resize(
    grl@unlistData[newends],
    width(grl@unlistData[newends]) - endtrims
  )
  grl
}


#' Resize a GRangesList object holding it's 3' end fixed
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param grl GRangesList; a GRangesList object
#' @return A GRangesList object shortend/lengthened, respecting exon boundaries

resize_grl_endfix <- function(grl, width) {
  grl <- invertStrand(grl) %>% sort_grl_st()
  #
  grl <- resize_grl_startfix(grl, width)
  invertStrand(grl) %>% sort_grl_st()
}


#' Resize a GRangesList object, respecting exon boundaries when shortening
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param grl GRangesList; a GRangesList object
#' @param width GRangesList; integer/IntegerList to set as new width.
#' @return A GRangesList object shortend/lengthened, respecting exon boundaries

resize_grl <- function(grl, gwidth, fix = "start", check = TRUE) {
  stopifnot(all(gwidth > 0))
  stopifnot(all(all(diff(order_grl_st(grl)) == 1)))
  if (fix == "start") {
    grl <- resize_grl_startfix(grl, gwidth)
  } else if (fix == "end") {
    grl <- resize_grl_endfix(grl, gwidth)
  } else if (fix == "center") {
    grlwidths <- sum(width(grl))
    diffs <- (gwidth - grlwidths)
    #
    grl <- resize_grl_startfix(grl, grlwidths + ceiling(diffs / 2))
    grl <- resize_grl_endfix(grl, grlwidths + diffs)
  }
  if (check) {
    startstoolow <- any(start(grl) <= 0)
    if (any(startstoolow)) {
      errortxt <- str_interp(paste0(
        "${sum(startstoolow)} ranges extended below",
        " 1 .. e.g. ${head(which(startstoolow,1))}"
      ))
      stop(errortxt)
    }
    intlistinds <- IntegerList(as.list(rep(1, length(grl))))
    grlseqs <- as.vector(unlist(use.names = FALSE, seqnames(grl)[intlistinds]))
    endhighvect <- (end(grl) > GenomeInfoDb::seqlengths(grl)[grlseqs])
    endstoohigh <- any(endhighvect)
    if (any(endstoohigh)) {
      errortxt <- str_interp(paste0(
        "${sum(endstoohigh)} ranges extended below ",
        "above seqlength .. e.g. ${head(which(endstoohigh,1))}"
      ))
      stop(errortxt)
    }
  }
  grl
}


#' Pick columns from a GRangestList
#'
#' Given a grangelist of say N genes with X_n exons, this yields a
#' length N vector pulled from the mcols of the first element of each list 
#' element
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param grl String; full path to html report file.
#' @return length n vector pulled from mcols of first list elements


fmcols <- function(grl, ...) {
  startinds <- start(grl@partitioning)
  with(as.data.frame(grl@unlistData@elementMetadata), ...)[startinds]
}

#' Check if GRanges elements are out of bounds
#'
#' Given a grangelist of say N genes with X_n exons, this yields a
#' length N vector pulled from the mcols of the first element of each 
#' list element
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param grl String; full path to html report file.
#' @return a logical vector valued TRUE if GRanges
#' elements are out of the chromosome bounds

is_out_of_bounds <- function(gr, si = seqinfo(gr)) {
  if (is(gr, "GenomicRangesList")) {
    grchrs <- as.character(seqnames(gr@unlistData))
    is_out <- end(gr) > split(
      seqlengths(si)[grchrs],
      gr@partitioning
    )
  } else {
    seqinfo(gr) <- si
    is_out <- end(gr) > seqlengths(gr)[as.character(seqnames(gr))]
  }
  start(gr) < 1 | is_out
}

#' Map From a transcript to the genome, splitting elements by exons
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param trspacegr GRanges; an object in transcript space, to be mapped back
#' to the genome
#' @param exonsgrl exonsgrl; exons making up the space element is to be mapped
#' from.
#' @return a granges object containing 1 or more element for each
#' transcript space range, in genome space, corresponding to pieces
#' of each element split by exon boundaries
#' @import S4Vectors

spl_mapFromTranscripts <- function(trspacegr, exons_grl) {
  exons_tr <- exons_grl %>%
    unlist() %>%
    GenomicFeatures::mapToTranscripts(exons_grl) %>%
    .[names(.) == seqnames(.)]
  ov <- findOverlaps(trspacegr, exons_tr)
  # make sure all our elements have exons
  stopifnot(all(unlist(unique(seqnames(trspacegr))) %in% names(exons_grl)))
  stopifnot((seq_along(trspacegr)) %in% queryHits(ov))
  # multiply our ranges
  trspacegr_spl <- suppressWarnings({
    trspacegr[queryHits(ov)]
  })
  # limit them to overlap one exon
  trspacegr_spl <- suppressWarnings({
    pintersect(trspacegr_spl, exons_tr[subjectHits(ov)])
  })
  # now map to the genome
  genomic_trspacegr <- GenomicFeatures::mapFromTranscripts(
    trspacegr_spl,
    exons_grl
  )
  # note the mapping
  genomic_trspacegr$xHits <- queryHits(ov)[genomic_trspacegr$xHits]
  genomic_trspacegr
}


#' Check if a granges list of CDS have start codons
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param trspacegr GRanges; an object in transcript space, to be mapped back
#' to the genome
#' @param exonsgrl exonsgrl; exons making up the space element is to be mapped
#' from.
#' @return a granges object containing 1 or more element for each
#' transcript space range, in genome space, corresponding to pieces
#' of each element split by exon boundaries

# now only those which have M at the start and '*' at the end
hasMstart <- function(cdsgrl, fafileob) {
  cdsseqstarts <- cdsgrl %>%
    sort_grl_st() %>%
    resize_grl(3, "start") %>%
    GenomicFeatures::extractTranscriptSeqs(x = fafileob, .) %>%
    Biostrings::translate(.)
  cdsseqstarts == "M"
}

#' Check if a granges list of CDS have start codons
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param trspacegr GRanges; an object in transcript space, to be mapped back
#' to the genome
#' @param exonsgrl exonsgrl; exons making up the space element is to be mapped
#' from.
#' @return a granges object containing the coding sequence range for each
#' transcript

get_trspace_cds <- function(cdsgrl, exonsgrl) {
  # now lift cds to exons space
  # nouorf <- cdsgrl%>%names%>%str_detect('_')%>%`!`
  trspacecds <- 
    cdsgrl%>%
    # cdsgrl[nouorf]%>% 
    GenomicFeatures::pmapToTranscripts(
      exonsgrl[fmcols(., transcript_id)]
      # exonsgrl['ENST00000000442.11']
    )
  stopifnot(all(elementNROWS(trspacecds)==1))
  trspacecds <- unlist(trspacecds)
  strand(trspacecds) <- '+'
  trspacecds
}


#' Get a set of filtered cds from an imported GTF GRanges
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param filt_anno GRanges; an unfilted imported GTF
#' @param fafileob FaFile; a genomic Fasta file
#' @details This takes only coding sequences which are a multiple of 3bp and
#' have a start and a stop on either end.it always returns coding sequences
#' without the stops, regardless of their extent in the input.
#' @return a GRangesList split by transcript, which contains the filtered coding
#' ranges for each one

# TODO update this for uORFs
get_cdsgrl <- function(filt_anno, fafileob, ignore_orf_validity) {
  # find which cds are multiples of 3bp
  #
  cdsgrl <- filt_anno %>%
    subset(., type == "CDS") %>%
    split(., .$transcript_id)

  is3bp <- cdsgrl %>%
    width() %>%
    sum() %>%
    `%%`(3) %>%
    `==`(0)

  cdsgrl <- cdsgrl[is3bp]
  message(str_interp(paste0(
    "filtered out ${sum(!is3bp)} ORFs for not being",
    " multiples of 3bp long"
  )))
  # chceck if the cds includes the stop codon
  cdsgrl <- sort_grl_st(cdsgrl)
  cdsseqends <- cdsgrl %>%
    resize_grl(3, "end") %>%
    resize_grl(sum(width(.)) + 3) %>%
    GenomicFeatures::extractTranscriptSeqs(x = fafileob, .)

  # some sequences have spaces (ends of chrs i think)
  filterchars <- cdsseqends %>% str_detect("[^ATCG]")
  cdsseqends[filterchars] <- "AAAAAA"
  cdsseqends <- Biostrings::translate(cdsseqends)
  stopifnot(Biostrings::nchar(cdsseqends)%in%2)
  # now determine if the annotations 'cds' include stop codons
  # if they do, fix that.
  end_stop <- table(subseq(cdsseqends, 1, 1)) %>%
    sort() %>%
    {
      . / sum(.)
    } %>%
    .["*"] %>%
    `>`(0.5)
  if (is.na(end_stop)) end_stop <- FALSE
  end_plusone_stop <- table(subseq(cdsseqends, 2, 2)) %>%
    sort() %>%
    {
      . / sum(.)
    } %>%
    .["*"] %>%
    `>`(0.5)
  stopifnot(end_stop | end_plusone_stop)
  if (end_stop) cdsgrl <- cdsgrl %>% resize_grl(sum(width(.)) - 3, "start")
  #
  endseq <- if (end_plusone_stop) {
    Biostrings::subseq(cdsseqends, 2, 2)
  } else {
    Biostrings::subseq(cdsseqends, 1, 1)
  }
  hasstop <- endseq == "*"
  if (!ignore_orf_validity) {
    cdsgrl <- cdsgrl[hasstop]
    message(str_interp("filtered out ${sum(!hasstop)} ORFs not ending with *"))
  }
  hasM <- hasMstart(cdsgrl, fafileob)
  if (!ignore_orf_validity) {
    cdsgrl <- cdsgrl[hasM]
    message(str_interp("filtered out ${sum(!hasM)} ORFs not starting with M"))
  }
  message(str_interp("${length(cdsgrl)} ORFs left"))
  cdsgrl
}

#' given a granges object, convert it to width1 granges, preserving mcols
#'
#' Given rpfs, and named vectors of starts and stops, this adds the 'phase'
#' to the reads
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param gr GRanges object, sorted
#'
#' @return A granges object with the each position in as it's own range

width1grs <- function(gr) {
  stopifnot(Negate(is.unsorted)(gr))
  isw1 <- width(gr) == 1
  broad <- gr[!isw1]
  # vector of integers - 1,2,3 for range 1-3
  narrowstarts <- unlist(as(broad@ranges, "IntegerList"))
  narrow <- {
    GRanges(
      rep(seqnames(broad), width(broad)),
      IRanges(narrowstarts, w = 1)
    )
  }
  binds <- rep(seq_along(broad), width(broad))
  mcols(narrow) <- mcols(broad)[binds, , drop = FALSE]
  sort(c(gr[isw1], narrow))
}

################################################################################
##########
################################################################################


#' Get a set of filtered cds from an imported GTF GRanges
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param gtf A GTF with gene annotation
#' @param fafileob FaFile; a genomic Fasta file
#' @param add_uorfs Whether to look for and include uORFs
#' @param keep_cols columns to save from the gtf metadata
#' @details This takes only coding sequences which are a multiple of 3bp and
#' have a start and a stop on either end. it always returns coding sequences
#' without the stops, regardless of their extent in the input.
#' @return a list containing annotation objects used by other functions.
#' @export
#' @examples
#'  gtf <- system.file('extdata', 'gcv37.anno.chr22.gtf', package='Ribostan',
#'         mustWork=TRUE)
#'  fafile <- system.file('extdata', 'chr22.fa.gz', package='Ribostan',
#'     mustWork=TRUE)
#'  file.copy(fafile, '.')
#'  system('gunzip -f chr22.fa.gz')
#'  fafile = 'chr22.fa'
#'  anno <- load_annotation(gtf, fafile)
#' 

load_annotation <- function(gtf, fafile, add_uorfs=TRUE,
  ignore_orf_validity = FALSE,
  keep_cols= c('gene_id','transcript_id','gene_name','type')) {
  anno <- rtracklayer::import(gtf)
  stopifnot(all(keep_cols%in%colnames(mcols(anno))))
  anno <- anno[,keep_cols]
  fafileob <- Rsamtools::FaFile(fafile)
  Rsamtools::indexFa(fafile)
  seqinfo(anno) <- seqinfo(fafileob)[as.vector(seqlevels(anno))]
  filt_anno <- anno
  #
  trgiddf <- anno %>%
    mcols() %>%
    .[, c("gene_id", "transcript_id")] %>%
    as.data.frame() %>%
    distinct() %>%
    filter(!is.na(transcript_id))
  # get the cds not including stop codons, possibly filtering for valid orfs
  cdsgrl <- get_cdsgrl(filt_anno, fafileob, ignore_orf_validity)
  #
  if(add_uorfs){
    message('adding uORFs..')
    anno$phase <- NULL
    txdb <- GenomicFeatures::makeTxDbFromGRanges(anno)
    fiveutrs = GenomicFeatures::fiveUTRsByTranscript(txdb, use.names=TRUE)
    alluORFs <- ORFik::findUORFs(fiveutrs, fafile)
    alluORFs <- alluORFs%>%{.@unlistData@ranges@NAMES<-NULL;.}%>%unlist
    # uorftrs <- names(alluORFs_ul)%>%str_replace('_\\d+$', '')
    alluORFs$transcript_id <- names(alluORFs)%>%str_replace('_\\d+$', '')
    alluORFs$type <- 'CDS'
    alluORFs$gene_id <- trgiddf$gene_id[
      match(alluORFs$transcript_id,trgiddf$transcript_id)]
    #remove stop codon from uORFs  
    alluORFs <- alluORFs%>%split(.,names(.))
    stopifnot(is(alluORFs,'GRangesList'))
    seqinfo(alluORFs)<-seqinfo(anno)
    alluORFs<-alluORFs%>%resize_grl(sum(width(.))-3,'start')
    #add uorfs to cdsgrl
    cdsgrl <- c(cdsgrl,alluORFs)
    #now modify metadata
    names(cdsgrl@unlistData)<-NULL
    trgiddf <- unlist(cdsgrl)%>%
      {data.frame(
        orf_id = names(.),
        transcript_id = .$transcript_id,
        gene_id = .$gene_id
      )}%>%
      distinct
    trgiddf$uORF <- trgiddf$transcript_id%in%names(alluORFs)
    is_uORF <-  names(cdsgrl)%in%names(alluORFs)
    names(is_uORF) <- names(cdsgrl)
    message('uORFs found')
  }else{
    is_uORF <- rep(FALSE, length(cdsgrl))
  }
  #
  ribocovtrs <- unique(trgiddf$transcript_id)
  # subset cds and anno with these
  filt_anno <- filt_anno %>%
    subset(type != "CDS") %>%
    subset(transcript_id %in% ribocovtrs)
  #
  filt_anno <- c(filt_anno, unlist(cdsgrl))
  exonsgrl <- filt_anno %>%
    subset(type == "exon") %>%
    split(., .$transcript_id)
  exonsgrl <- exonsgrl[ribocovtrs]
  trspacecds <- get_trspace_cds(cdsgrl, exonsgrl)
  #
  cdsstarts <- trspacecds %>%
    start() %>%
    setNames(names(trspacecds))
  #
  longtrs <- width(exonsgrl)%>%
    sum%>%tibble::enframe('transcript_id','width')%>%
    left_join(trgiddf, 'transcript_id')%>%
    group_by(gene_id)%>%
    arrange(-width)%>%
    dplyr::slice(1)%>%
    .$transcript_id
  #
  outanno <- list(
    trspacecds = trspacecds,
    cdsgrl = cdsgrl,
    exonsgrl = exonsgrl,
    trgiddf = trgiddf,
    fafileob = fafileob,
    longtrs = longtrs,
    uORF = is_uORF
  )
  outanno <- c(
    outanno,
    list(
      cdsstarts = outanno$trspacecds %>% start() %>%
        setNames(names(outanno$trspacecds)),
      cds_prestop_st = outanno$trspacecds %>% end() %>% `-`(2) %>%
        setNames(names(outanno$trspacecds)),
      anno = filt_anno %>% subset(transcript_id %in% ribocovtrs)
    )
  )
  return(outanno)
}


#' Subset an annotation object with a vector of ORF ids
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param orfs A vector of ORF ids
#' @param anno an annotation object
#' @return an annotation object, subsetted
#' 
subset_annotation <- function(orfs, anno){
  newanno <- anno
  newanno$trspacecds <- anno$trspacecds[orfs]
  newanno$cdsgrl <- anno$cdsgrl[orfs]
  orftrs <- unique(unlist(fmcols(anno$cdsgrl[orfs], transcript_id)))
  newanno$exonsgrl <- anno$exonsgrl[orftrs]
  newanno$trgiddf <- anno$trgiddf%>%filter(orf_id %in% orfs)
  newanno$fafileob <- anno$fafileob
  newanno$longtrs <- anno$longtrs%>%intersect(orftrs)
  newanno$uORF <- anno$uORF[orfs]
  newanno
}


#


#' Create a fasta file of coding sequences extended nbp on either end.
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param gtf A GTF file with annotation.
#' @param fasta A genomic Fasta file
#' @param outfasta The name of the fasta file to output
#' @param fpext How many bp to extendt the CDS on the 5' end
#' @param tpext How many bp to extendt the CDS on the 3' end
#' @details This creates a fasta file with gencode-style headers. The sequences
#' are created by extending the coding sequences by the specified amount in
#' transcript space, (or past the) end of the transcript if necessary, so that
#' all sequences have the same 'UTRs'.
#' @import GenomicRanges
#' @return the name of the file to which the sequences were output
#' @examples 
#' gtf <- system.file('extdata', 'gcv37.anno.chr22.gtf', package='Ribostan',
#'   mustWork=TRUE)
#' fafile <- system.file('extdata', 'chr22.fa.gz', package='Ribostan',
#'      mustWork=TRUE)
#' file.copy(fafile, '.')
#' system('gunzip -f chr22.fa.gz')
#' fafile = 'chr22.fa'
#' ext_fasta = make_ext_fasta(gtf, fafile, outfasta='tmp.fa', fpext=50, tpext=50)
#'
#'

make_ext_fasta <- function(gtf, fasta, outfasta, fpext = 50, tpext = 50) {
  stopifnot({
    cat("testing", file = outfasta)
    file.remove(outfasta)
  })

  stopifnot(gtf %>% str_detect("\\.(gtf)$"))
  stopifnot(gtf %>% file.exists())

  stopifnot(outfasta %>% str_detect("\\.(fasta|fa)$"))
  outprefix <- outfasta %>% str_replace("\\.(fasta|fa)$", "")

  # get our filtered annotation
  keepcols = c('transcript_id','type',
    'gene_id',
    'havana_gene',
    'havana_transcript',
    'transcript_name',
    'gene_name')
  anno <- load_annotation(gtf, fasta, add_uorfs=FALSE,
    keep_cols = keepcols)
  cdsgrl <- anno$cdsgrl
  exonsgrl <- anno$exonsgrl[names(cdsgrl)]
  cdsexonsgrl <- anno$exonsgrl[names(cdsgrl)]
  cdsstartpos <- start(anno$trspacecds)
  # get exons for our cds
  cdsexonsgrl <- sort_grl_st(cdsexonsgrl)
  # get an object representing the CDS In transript space
  cdstrspace <- anno$trspacecds[names(cdsgrl)]
  endpos <- sum(width(cdsexonsgrl)) - end(cdstrspace)
  # expand our first exon when needed
  startposexpansion <- pmax(0, fpext - cdsstartpos + 1)
  # expand/trim the 5' end of the exons
  startinds <- start(cdsexonsgrl@partitioning)
  cdsexonsgrl@unlistData[startinds] <- cdsexonsgrl@unlistData[startinds]%>%
  resize(width(.) + startposexpansion, "end")
  # expand or trim the last exon when needed
  endposexpansion <- pmax(0, tpext - endpos)
  endinds <- cdsexonsgrl@partitioning@end
  cdsexonsgrl@unlistData[endinds] <- cdsexonsgrl@unlistData[endinds] %>%
  resize(width(.) + endposexpansion, "start")
  cds_exptrspc <- GenomicFeatures::pmapToTranscripts(cdsgrl,
      cdsexonsgrl[names(cdsgrl)])
  stopifnot(cds_exptrspc %>% elementNROWS() %>% `==`(1))


  expcds_exptrspc <- cds_exptrspc
  stopifnot(!any(expcds_exptrspc %>% elementNROWS() %>% `>`(1)))
  expcds_exptrspc <- unlist(expcds_exptrspc)
  # expand our cds exons
  expcds_exptrspc <- expcds_exptrspc%>%resize(.,width(.) + fpext, "end",
    ignore.strand = TRUE)
  # and expand the 3' ends
  expcds_exptrspc <- expcds_exptrspc %>% resize(., width(.) + tpext, "start",
    ignore.strand = TRUE)
  # now back to genome space
  expcdsgenspace <- spl_mapFromTranscripts(expcds_exptrspc, cdsexonsgrl)
  seqinfo(expcdsgenspace) <- seqinfo(cdsexonsgrl)
  expcdsgenspace <- GenomicRanges::split(expcdsgenspace, names(expcdsgenspace))
  # get the sequences
  isoutofbds <- is_out_of_bounds(expcdsgenspace)
  isoutofbds <- any(is_out_of_bounds(expcdsgenspace))
  message(str_interp(paste0(
    "Excluded ${sum(isoutofbds)} orfs because they",
    " extended beyond chromosomal boundaries"
  )))
  expcdsgenspace <- expcdsgenspace[!isoutofbds]
  cdsexonsgrl <- cdsexonsgrl[names(expcdsgenspace)]
  cds_exptrspc <- cds_exptrspc[names(expcdsgenspace)]
  expcdsgenspaceseq <-
    expcdsgenspace %>%
    sort_grl_st() %>%
    GenomicFeatures::extractTranscriptSeqs(., x = anno$fafileob)

  cdslens <- sum(width(cdsgrl))[names(expcdsgenspace)]
  fastanames <- paste(
    sep = "|",
    fmcols(cdsexonsgrl, transcript_id),
    fmcols(cdsexonsgrl, gene_id),
    fmcols(cdsexonsgrl, havana_gene),
    fmcols(cdsexonsgrl, havana_transcript),
    fmcols(cdsexonsgrl, transcript_name),
    fmcols(cdsexonsgrl, gene_name),
    sum(width(cdsexonsgrl)),
    paste0("UTR5:1-", fpext),
    paste0("CDS:", fpext + 1, "-", fpext + cdslens),
    paste0("UTR3:", 1 + fpext + cdslens, "-", sum(width(expcdsgenspace))),
    "|"
  )

  names(expcdsgenspaceseq) <- fastanames

  trcdscoordsfile <- paste0(outprefix, "_trcds.tsv")
  as.data.frame(unlist(cds_exptrspc)) %>%
    select(seqnames, start, end) %>%
    write_tsv(trcdscoordsfile)
  message(normalizePath(trcdscoordsfile, mustWork = TRUE))

  # also write our cds coordinates to disk in the new trspace
  nms_cds_exptrspc <- as.character(seqnames(cds_exptrspc))
  new_trspc_anno <- c(
    GRanges(nms_cds_exptrspc, IRanges(1, sum(width(cdsexonsgrl)))) %>%
      {
        .$type <- "exon"
        .
      },
    cds_exptrspc %>% unlist() %>%
      {
        .$type <- "CDS"
        .
      }
  )
  new_trspc_anno %>%
    suppressWarnings({
      rtracklayer::export(paste0(outprefix, "_trspaceanno.gtf"))
    })

  # write the expanded cds exon sequences to disk
  Biostrings::writeXStringSet(expcdsgenspaceseq, outfasta)
  message(normalizePath(outfasta, mustWork = TRUE))

  # now make fasta file with shorter transcript names
  shortheaderfasta <- paste0(outprefix, ".shortheader.fa")
  system(str_interp("sed -e 's/|.*$//' ${outfasta} > ${shortheaderfasta}"))
  message(normalizePath(shortheaderfasta, mustWork = TRUE))
  return(outfasta)
}




#' This creates a minimal annotation object from a gencode style fasta
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param ribofasta A gencode style fasta to which RPFs were aligned
#' @return An annotation object with cdsstarts/stops
#' @details The Ribosome densities are saved in salmon format
#' @examples 
#' gtf <- system.file('extdata', 'gcv37.anno.chr22.gtf', package='Ribostan', mustWork=TRUE)
#' fafile <- system.file('extdata', 'chr22.fa.gz', package='Ribostan', mustWork=TRUE)
#' file.copy(fafile, '.')
#' system('gunzip -f chr22.fa.gz')
#' fafile = 'chr22.fa'
#' ext_fasta = make_ext_fasta(gtf, fafile, outfasta='tmp.fa', fpext=50, tpext=50)
#' get_ribofasta_anno(ext_fasta)
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
      str_split_fixed("-", 2) %>% 
      {colnames(.)<-c("start", "end"); .} %>%
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
