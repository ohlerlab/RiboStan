################################################################################
########## Range manipulation
################################################################################
#' @importFrom stringr str_detect str_replace str_subset str_interp
#' @importFrom stringr str_split_fixed str_extract
#' @importFrom readr write_tsv read_tsv
#' @import testthat
#' @import ggplot2
#' @importFrom IRanges IRanges subsetByOverlaps overlapsAny pintersect
#' @importFrom S4Vectors elementNROWS List subjectHits queryHits `%in%`
#' @importFrom GenomeInfoDb Seqinfo keepSeqlevels seqnames seqlevels seqlengths 
#' @importFrom GenomeInfoDb seqinfo seqlevels<- seqinfo<- seqnames<-
#' @importFrom BiocGenerics intersect setdiff unlist table union mean order 
#' @importFrom Biostrings subseq codons oligonucleotideFrequency nchar
#' @importFrom GenomicRanges GRanges split strand mcols width
#' @importFrom GenomicRanges strand<- mcols<- start resize
#' @importFrom GenomicRanges findOverlaps invertStrand seqnames end
#' @importFrom GenomicRanges coverage shift
#' @importFrom parallel detectCores
#' @importFrom rtracklayer import export
#' @importFrom Rsamtools ScanBamParam
#' @importFrom tidyr replace_na unnest
#' @importFrom dplyr mutate select filter lead summarise tally slice %>% lag
#' @importFrom dplyr left_join group_by ungroup tibble inner_join bind_rows
#' @importFrom dplyr distinct arrange n count between full_join one_of 
#' @importFrom dplyr n_distinct mutate_at

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
  stopifnot(length(grl)>0)
  order(GenomicRanges::start(grl) * (((strand(grl) != "-") + 1) * 2 - 3))
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
  iList <- IRanges::IntegerList(as.list(elementNROWS(grl)))
  endtrims <- trim[iList]@unlistData
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
#' @param grl the width to resize the GRL to; integer
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
#' @param gwidth GRangesList; integer/IntegerList to set as new width.
#' @return A GRangesList object shortend/lengthened, respecting exon boundaries

resize_grl <- function(grl, gwidth, fix = "start", check = TRUE) {

  stopifnot(all(gwidth > 0))
  stopifnot(all(all(diff(order_grl_st(grl)) == 1)))
  stopifnot(is.vector(gwidth))
  if(length(gwidth)==1) gwidth <- rep(gwidth, length(grl))
  stopifnot((length(gwidth)==length(grl)))

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
    startstoolow <- any(GenomicRanges::start(grl) <= 0)
    if (any(startstoolow)) {
      errortxt <- str_interp(paste0(
        "${sum(startstoolow)} ranges extended below",
        " 1 .. e.g. ${head(which(startstoolow,1))}"
      ))
      stop(errortxt)
    }
    intlistinds <- IRanges::IntegerList(as.list(rep(1, length(grl))))
    grlseqs <- as.vector(unlist(use.names = FALSE, seqnames(grl)[intlistinds]))
    endhighvect <- (GenomicRanges::end(grl) > GenomeInfoDb::seqlengths(grl)[grlseqs])
    iscirc <- seqinfo(grl)@is_circular[match(grlseqs,seqinfo(grl)@seqnames)]
    endhighvect[iscirc%in%TRUE]<-FALSE
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
  startinds <- GenomicRanges::start(grl@partitioning)
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
    is_out <- GenomicRanges::end(gr) > split(
      seqlengths(si)[grchrs],
      gr@partitioning
    )
  } else {
    seqinfo(gr) <- si
    is_out <- GenomicRanges::end(gr) > seqlengths(gr)[as.character(seqnames(gr))]
  }
  GenomicRanges::start(gr) < 1 | is_out
}

#' Map From a transcript to the genome, splitting elements by exons
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param trspacegr GRanges; an object in transcript space, to be mapped back
#' to the genome
#' @param exons_grl exonsgrl; exons making up the space element is to be mapped
#' from.
#' @return a granges object containing 1 or more element for each
#' transcript space range, in genome space, corresponding to pieces
#' of each element split by exon boundaries

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
    Biostrings::translate(., if.fuzzy.codon = "solve")
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
    cdsgrl %>%
    # cdsgrl[nouorf]%>%
    GenomicFeatures::pmapToTranscripts(
      exonsgrl[fmcols(., transcript_id)]
      # exonsgrl['ENST00000000442.11']
    )
  stopifnot(all(elementNROWS(trspacecds) == 1))
  trspacecds <- unlist(trspacecds)
  strand(trspacecds) <- "+"
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
    resize_grl(sum(width(.)) + 3)

    cdsgrl%>%width%>%sum
    cdsgrl%>%resize_grl(3,'end')

  cdsseqends <- cdsseqends%>% 
    GenomicFeatures::extractTranscriptSeqs(x = fafileob, .)

  # some sequences have spaces (ends of chrs i think)
  filterchars <- cdsseqends %>% as.vector() %>% str_detect("[^ATCG]")
  cdsseqends[filterchars] <- "AAAAAA"
  cdsseqends <- Biostrings::translate(cdsseqends)
  stopifnot(Biostrings::nchar(cdsseqends) %in% 2)
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
      IRanges(narrowstarts, width = 1)
    )
  }
  binds <- rep(seq_along(broad), width(broad))
  mcols(narrow) <- mcols(broad)[binds, , drop = FALSE]
  sort(c(gr[isw1], narrow))
}

################################################################################
##########
################################################################################

#' Convert old-style gtfs (e.g. those output by rtracklayer) to richer format
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param an annotation GRanges read in by rtracklayer
#' @param keep_cols the columns to keep in the final annotation
#' @return a granges object with just the exon and CDS elements 
#' that ribostan uses
convert_gtf <- function(anno, keep_cols){
  genedf <- mcols(anno)%>%as.data.frame%>%
    filter(.data$type=='gene')%>%
    mutate('gene_id'=.data$ID%>%str_replace('GeneID:',''))%>%
    select('type','gene_id','gene_name'='Name')
  trdf <- mcols(anno)%>%as.data.frame%>%filter(.data$type=='mRNA')%>%
    mutate('transcript_id'=.data$ID)%>%
    mutate('gene_id'=.data$Parent%>%str_replace('GeneID:',''))%>%
    mutate('gene_name'=.data$Parent%>%str_replace('GeneID:',''))%>%
    mutate('type'='transcript')%>%
    select('type','gene_name','gene_id','transcript_id',
      'transcript_name'='Name')
  exondf <- mcols(anno)%>%as.data.frame%>%filter(.data$type=='exon')%>%
    mutate('exon_id'='Name')%>%
    mutate('transcript_id'=
      .data$Parent%>%stringr::str_replace_all('TxID:',''))%>%
    select('type','exon_id','transcript_id','exon_name'='Name')
  #transcript_id is multiple
  cds_df <- mcols(anno)%>%as.data.frame%>%
    filter(.data$type=='CDS')%>%
    mutate('transcript_id'=Parent%>%stringr::str_replace_all('TxID:',''))%>%
    select('type','transcript_id')
  allcds <- anno%>%subset(type=='CDS')
  parentsplit <- stringr::str_split(allcds$Parent,',')
  parentnum <- parentsplit%>%as("CharacterList")%>%elementNROWS
  allcds <- rep(allcds,parentnum)
  allcds$Parent <- unlist(parentsplit)
  allexons <- anno%>%subset(type=='exon')
  parentsplit <- stringr::str_split(allexons$Parent,',')
  parentnum <- parentsplit%>%as("CharacterList")%>%elementNROWS
  allexons <- rep(allexons,parentnum)
  allexons$Parent <- unlist(parentsplit)
  #
  allexons$transcript_id <- allexons$Parent
  allcds$transcript_id <- allcds$Parent
  trmatch<-match(allcds$transcript_id,trdf$transcript_id)
  allcds$transcript_id <- trdf$transcript_id[trmatch]
  allcds$gene_id <- trdf$gene_id[trmatch]
  allcds$gene_name <- trdf$gene_name[trmatch]
  extrmatch<-match(allexons$transcript_id,trdf$transcript_id)
  allexons$transcript_id <- trdf$transcript_id[extrmatch]
  allexons$gene_id <- trdf$gene_id[extrmatch]
  allexons$gene_name <- trdf$gene_name[extrmatch]
  anno <- c(allexons,allcds)
  stopifnot(!any(is.na(anno$transcript_id)))
  stopifnot(!any(is.na(anno$gene_name)))
  stopifnot(!any(is.na(anno$gene_id)))
  stopifnot(!any(is.na(anno$type)))
  anno
}

#' Get a set of filtered cds from an imported GTF GRanges
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param gtf A GTF with gene annotation
#' @param fafileob FaFile; a genomic Fasta file
#' @param add_uorfs Whether to look for and include uORFs
#' @param keep_cols columns to save from the gtf metadata
#' @param DEFAULT_CIRC_SEQS default chromsoomes to treat as circular
#' @param findUORFs_args Additional arguments to pass to ORFik::findUORFs
#' @details This takes only coding sequences which are a multiple of 3bp and
#' have a start and a stop on either end. it always returns coding sequences
#' without the stops, regardless of their extent in the input.
#' @return a list containing annotation objects used by other functions.
#' @export
#' @examples
#' gtf <- system.file("extdata", "gcv37.anno.chr22.gtf",
#'   package = "Ribostan",
#'   mustWork = TRUE
#' )
#' fafile <- system.file("extdata", "chr22.fa.gz",
#'   package = "Ribostan",
#'   mustWork = TRUE
#' )
#' file.copy(fafile, ".")
#' system2("gunzip -f chr22.fa.gz")
#' fafile <- "chr22.fa"
#' anno <- load_annotation(gtf, fafile)
load_annotation <- function(
    gtf, fafile, add_uorfs = TRUE,
    ignore_orf_validity = FALSE,
    keep_cols = 
      c("gene_id", "transcript_id", "gene_name", "type"),
    DEFAULT_CIRC_SEQS=NULL,
    findUORFs_args = c(minimumLength=0)
    ) {
  if(is.null(DEFAULT_CIRC_SEQS)){
    DEFAULT_CIRC_SEQS <- unique(
      c("chrM","MT","MtDNA","mit","Mito","mitochondrion",
      "dmel_mitochondrion_genome","Pltd","ChrC","Pt","chloroplast",
      "Chloro","2micron","2-micron","2uM",
      "Mt", "NC_001879.2", "NC_006581.1","ChrM","mitochondrion_genome"))
  }
  anno <- rtracklayer::import(gtf)
  ancols <- colnames(mcols(anno))
  is_compressedgtf <- ('Parent'%in% ancols)&(!'transcript_id'%in%ancols)
  if(is_compressedgtf){
    message('reformatting gtf to include transcript_id etc in mcols')
    anno <- convert_gtf(anno, keep_cols)
  }
  stopifnot(length(anno)>0)
  #for older gencode annotation
  stopifnot(c('exon','CDS')%in%anno$type)
  stopifnot(all(keep_cols %in% colnames(mcols(anno))))

  anno <- anno[, keep_cols]
  stopifnot(file.exists(fafile))
  fafileob <- Rsamtools::FaFile(fafile)
  Rsamtools::indexFa(fafile)
  tokeep <- seqlevels(anno)%>%intersect(seqlevels(seqinfo(fafileob)))
  stopifnot(length(seqlevels(seqinfo(fafileob)))>0)
  stopifnot(seqlevels(seqinfo(fafileob))>0)
  stopifnot(length(tokeep)>0)
  toremove <- seqlevels(anno)%>%setdiff(seqlevels(seqinfo(fafileob)))
  nonempty <- intersect(toremove,seqnames(anno))
  message(str_interp(paste0('removing ${length(nonempty)} non empty seqlevels',
    ' that are absent from the fasta')))
  anno <- anno %>% GenomeInfoDb::dropSeqlevels(nonempty,pruning.mode='coarse')
  empty <- setdiff(toremove,seqnames(anno))
  message(str_interp(paste0('removing ${length(empty)} non empty seqlevels',
    ' that are absent from the fasta')))
  anno <- anno %>% GenomeInfoDb::dropSeqlevels(empty,pruning.mode='coarse')
  seqinfo(anno) <- seqinfo(fafileob)[as.vector(seqlevels(anno))]
  seqinfo(anno)@is_circular <- seqinfo(anno)@seqnames %in% DEFAULT_CIRC_SEQS
  #
  trgiddf <- anno %>%
    mcols() %>%
    .[, c("gene_id", "transcript_id")] %>%
    as.data.frame() %>%
    distinct() %>%
    filter(!is.na(.data$transcript_id))
  trgiddf$orf_id <- trgiddf$transcript_id
  # get the cds not including stop codons, possibly filtering for valid orfs
  cdsgrl <- get_cdsgrl(anno, fafileob, ignore_orf_validity)
  #add uORFs
  if (add_uorfs) {
    message("adding uORFs..")
    anno$phase <- NULL
    txdb <- GenomicFeatures::makeTxDbFromGRanges(anno)
    fiveutrs <- GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE)
    validutrs <- names(fiveutrs)%>%intersect(names(cdsgrl))
    fiveutrs <- fiveutrs[validutrs]
    alluORFs <- do.call(what=ORFik::findUORFs, args = c(findUORFs_args,list(
      fiveUTRs = fiveutrs, 
      fa = fafile ,
      cds = cdsgrl[validutrs]
    )))
    alluORFs <- alluORFs %>%
      {
        .@unlistData@ranges@NAMES <- NULL
        .
      } %>%
      unlist()
    alluORFs$transcript_id <- names(alluORFs) %>% str_replace("_\\d+$", "")
    alluORFs$type <- "CDS"
    alluORFs$gene_id <- trgiddf$gene_id[
      match(alluORFs$transcript_id, trgiddf$transcript_id)
    ]
    # remove stop codon from uORFs
    alluORFs <- alluORFs %>% split(., names(.))
    stopifnot(is(alluORFs, "GRangesList"))
    seqinfo(alluORFs) <- seqinfo(anno)
    alluORFs <- alluORFs %>% resize_grl(sum(width(.)) - 3, "start")
    # add uorfs to cdsgrl
    cdsgrl <- c(cdsgrl, alluORFs)
    # now modify metadata
    names(cdsgrl@unlistData) <- NULL
    trgiddf <- unlist(cdsgrl) %>%
      {
        data.frame(
          orf_id = names(.),
          transcript_id = .$transcript_id,
          gene_id = .$gene_id
        )
      } %>%
      distinct()
    trgiddf$uORF <- trgiddf$transcript_id %in% names(alluORFs)
    is_uORF <- names(cdsgrl) %in% names(alluORFs)
    names(is_uORF) <- names(cdsgrl)
    message("uORFs found")
  } else {
    is_uORF <- rep(FALSE, length(cdsgrl))
  }
  #
  orf_transcripts <- fmcols(cdsgrl, transcript_id)
  # subset cds and anno with these
  anno <- anno %>%
    subset(type != "CDS") %>%
    subset(transcript_id %in% orf_transcripts)
  #
  anno <- c(anno, unlist(cdsgrl))
  exonsgrl <- anno %>%
    subset(type == "exon") %>%
    sort_grl_st%>%
    split(., .$transcript_id)
  exon_tr_names <- names(exonsgrl)
  stopifnot(all(orf_transcripts %in% exon_tr_names))
  setdiff(orf_transcripts,exon_tr_names)
  
  exonsgrl <- exonsgrl[orf_transcripts]
  trspacecds <- get_trspace_cds(cdsgrl, exonsgrl)
  #
  cdsstarts <- trspacecds %>%
    GenomicRanges::start() %>%
    setNames(names(trspacecds))
  #
  longtrs <- width(exonsgrl) %>%
    sum() %>%
    tibble::enframe("transcript_id", "width") %>%
    left_join(trgiddf, "transcript_id") %>%
    group_by(.data$gene_id) %>%
    arrange(-.data$width) %>%
    dplyr::slice(1) %>%
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
      cdsstarts = outanno$trspacecds %>% GenomicRanges::start() %>%
        setNames(names(outanno$trspacecds)),
      cds_prestop_st = outanno$trspacecds %>% GenomicRanges::end() %>% `-`(2) %>%
        setNames(names(outanno$trspacecds))
    )
  )
  return(outanno)
}


#' Subset an annotation object with a vector of ORF ids
#'
#' @keywords Ribostan
#' @author Dermot Harnett, \email{dermot.p.harnett@gmail.com}
#'
#' @param anno an annotation object
#' @param orfs A vector of ORF ids
#' @return an annotation object, subsetted
#'
subset_annotation <- function(anno, orfs) {
  newanno <- anno
  newanno$trspacecds <- anno$trspacecds[orfs]
  newanno$cdsgrl <- anno$cdsgrl[orfs]
  orftrs <- unique(unlist(fmcols(anno$cdsgrl[orfs], transcript_id)))
  newanno$exonsgrl <- anno$exonsgrl[orftrs]
  newanno$trgiddf <- anno$trgiddf %>% filter(.data$orf_id %in% orfs)
  newanno$fafileob <- anno$fafileob
  newanno$longtrs <- anno$longtrs %>% intersect(orftrs)
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
#' @return the name of the file to which the sequences were output
#' @export
#' @examples
#' gtf <- system.file("extdata", "gcv37.anno.chr22.gtf",
#'   package = "Ribostan",
#'   mustWork = TRUE
#' )
#' fafile <- system.file("extdata", "chr22.fa.gz",
#'   package = "Ribostan",
#'   mustWork = TRUE
#' )
#' file.copy(fafile, ".")
#' system2("gunzip -f chr22.fa.gz")
#' fafile <- "chr22.fa"
#' ext_fasta <- make_ext_fasta(gtf, fafile, outfasta = "tmp.fa", fpext = 50, tpext = 50)
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
  keepcols <- c(
    "transcript_id", "type",
    "gene_id",
    # "transcript_name",
    "gene_name"
  )
  anno <- load_annotation(gtf, fasta,
    add_uorfs = FALSE,
    keep_cols = keepcols
  )
  cdsgrl <- anno$cdsgrl
  exonsgrl <- anno$exonsgrl[names(cdsgrl)]
  cdsexonsgrl <- anno$exonsgrl[names(cdsgrl)]
  cdsstartpos <- GenomicRanges::start(anno$trspacecds)
  # get exons for our cds
  cdsexonsgrl <- sort_grl_st(cdsexonsgrl)
  # get an object representing the CDS In transript space
  cdstrspace <- anno$trspacecds[names(cdsgrl)]
  endpos <- sum(width(cdsexonsgrl)) - GenomicRanges::end(cdstrspace)
  # expand our first exon when needed
  startposexpansion <- pmax(0, fpext - cdsstartpos + 1)
  # expand/trim the 5' end of the exons
  startinds <- GenomicRanges::start(cdsexonsgrl@partitioning)
  cdsexonsgrl@unlistData[startinds] <- cdsexonsgrl@unlistData[startinds] %>%
    resize(width(.) + startposexpansion, "end")
  # expand or trim the last exon when needed
  endposexpansion <- pmax(0, tpext - endpos)
  endinds <- cdsexonsgrl@partitioning@end
  cdsexonsgrl@unlistData[endinds] <- cdsexonsgrl@unlistData[endinds] %>%
    resize(width(.) + endposexpansion, "start")
  cds_exptrspc <- GenomicFeatures::pmapToTranscripts(
    cdsgrl,
    cdsexonsgrl[names(cdsgrl)]
  )
  stopifnot(cds_exptrspc %>% elementNROWS() %>% `==`(1))


  expcds_exptrspc <- cds_exptrspc
  stopifnot(!any(expcds_exptrspc %>% elementNROWS() %>% `>`(1)))
  expcds_exptrspc <- unlist(expcds_exptrspc)
  # expand our cds exons
  expcds_exptrspc <- expcds_exptrspc %>% resize(., width(.) + fpext, "end",
    ignore.strand = TRUE
  )
  # and expand the 3' ends
  expcds_exptrspc <- expcds_exptrspc %>% resize(., width(.) + tpext, "start",
    ignore.strand = TRUE
  )
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
    # fmcols(cdsexonsgrl, havana_gene),
    NA,
    # fmcols(cdsexonsgrl, havana_transcript),
    NA,
    # fmcols(cdsexonsgrl, transcript_name),
    NA,
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
    select('seqnames', 'start', 'end') %>%
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
  names(expcdsgenspaceseq) <- str_extract(names(expcdsgenspaceseq), '[^|]+')
  Biostrings::writeXStringSet(expcdsgenspaceseq, shortheaderfasta)
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
#' gtf <- system.file("extdata", "gcv37.anno.chr22.gtf", package = "Ribostan", mustWork = TRUE)
#' fafile <- system.file("extdata", "chr22.fa.gz", package = "Ribostan", mustWork = TRUE)
#' file.copy(fafile, ".")
#' system2("gunzip -f chr22.fa.gz")
#' fafile <- "chr22.fa"
#' ext_fasta <- make_ext_fasta(gtf, fafile, outfasta = "tmp.fa", fpext = 50, tpext = 50)
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
      {
        colnames(.) <- c("start", "end")
        .
      } %>%
      apply(2, as.numeric) %>% as.data.frame() %>%
      {
        IRanges(start = .$start, end = .$end)
      }
  ) %>%
    setNames(., as.character(seqnames(.)))
  strand(anno$trspacecds) <- "+"
  seqinfo(anno$trspacecds) <- Seqinfo( 
    faheaddf[,1],
    faheaddf[,10]%>%str_extract('(?<=-)\\d+')%>%as.numeric)
  anno$trgiddf <- tibble(transcript_id = faheaddf[, 1], gene_id = faheaddf[, 2])
  anno$trgiddf$orf_id <- anno$trgiddf$transcript_id
  anno <- c(
    anno,
    list(
    cdsstarts = anno$trspacecds %>% GenomicRanges::start() %>%
        setNames(names(anno$trspacecds)),
      cds_prestop_st = anno$trspacecds %>% GenomicRanges::end() %>% `-`(2) %>%
        setNames(names(anno$trspacecds))
    )
  )
  anno
}
