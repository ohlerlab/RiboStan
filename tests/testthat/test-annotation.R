test_that("loading annotation works", 
#(
 {
  #first, let's load up some test data
  anno_file <- here::here('test.gc32.gtf')
  #if(!file.exists(anno_file)){
    library(AnnotationHub)
    ah <- AnnotationHub()
    gencode32 <- ah[['AH75191']]
    seqlevels(gencode32)<-'chr22'
    suppressWarnings({rtracklayer::export(gencode32, anno_file, format='GTF')})
  #}
  fafile <- here::here('chr22.fa')
  library(BSgenome.Hsapiens.UCSC.hg38)
  #if(!file.exists(fafile)){
    seq <- Biostrings::DNAStringSet(BSgenome.Hsapiens.UCSC.hg38[['chr22']])
   names(seq)<-'chr22'
   Biostrings::writeXStringSet(
    seq, fafile)
  #}
  # file.copy(fafile, ".")
  # system("gunzip -f chr22.fa.gz")
  # fafile <- "chr22.fa"
  chr22_anno <- load_annotation(anno_file, fafile)
  # These totals reflect the current (conservative, ATG-only) uORF
  # discovery defaults: find_uorfs() only searches ATG starts and
  # .filter_uorfs() drops any candidate that coincides with the
  # annotated CDS boundary. 1493 real CDS transcripts + 1749 uORFs.
  expect_equal(chr22_anno %>% .$uORF %>% sum(), 1749)
  expect_equal(length(chr22_anno$cdsgrl), 3242)
  expect_equal(length(chr22_anno$longtrs), 442)
  expect_equal(length(chr22_anno$exonsgrl), 3242)
  expect_equal(length(names(chr22_anno$trspacecds)), 3242)
  expect_true(all(names(chr22_anno$exonsgrl) %in%
    Ribostan:::fmcols(chr22_anno$cdsgrl, transcript_id)))
  expect_true(all(names(chr22_anno$cdsgrl) %in% names(chr22_anno$trspacecds)))
  testorfs <- chr22_anno$trspacecds %>% .[c(1, 66, 666, length(chr22_anno$trspacecds))]
  testtrs <- Ribostan:::fmcols(chr22_anno$cdsgrl[names(testorfs)], transcript_id)
  testtrs <- unique(testtrs)
  exontestseq <- GenomicFeatures::extractTranscriptSeqs(chr22_anno$exonsgrl[testtrs],
    x = chr22_anno$fafileob)
  seqlevels(testorfs) <- names(exontestseq)
  testorfs <- resize(testorfs, width(testorfs) + 3)
  orftestseq <- exontestseq[testorfs]
  starts <- Biostrings::subseq(orftestseq, 1, 3) %>% as.character()
  names(starts) <- NULL
  # All ATG: with the ATG-only default there are no non-canonical starts.
  expect_equal(starts, c("ATG", "ATG", "ATG", "ATG"))
  stops <- Biostrings::subseq(orftestseq, -3, -1) %>%
    as.character() %>%
    setNames(NULL)
  expect_equal(stops, c("TAA", "TGA", "TGA","TAG"))
  expect_equal(
    chr22_anno$cdsgrl %>% .@unlistData %>% mcols() %>% colnames(),
    c("gene_id", "transcript_id", "gene_name", "type", "start_codon")
  )
  extfasta <- make_ext_fasta(anno_file, fafile, outfasta = "tmp.fa", fpext = 50, tpext = 50)
  extseqs <- rtracklayer::import(extfasta, format = "fasta")
  atgstarts <- extseqs %>%
    subseq(51, 53) %>%
    setNames(NULL) %>%
    as.character() %>%
    `%in%`("ATG")
  expect_true(mean(atgstarts) > 0.99)
  extorfs <- names(extseqs) %>% stringr::str_extract("[^|]+")
  expect_true(all(extorfs %in% (chr22_anno$uORF %>% .[!.] %>% names())))
}
)
