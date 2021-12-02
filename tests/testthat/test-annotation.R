test_that("loading annotation works", {
  gtf <- system.file('extdata', 'gcv37.anno.chr22.gtf', package='Ribostan', mustWork=TRUE)
  fafile <- system.file('extdata', 'chr22.fa.gz', package='Ribostan', mustWork=TRUE)
  file.copy(fafile, '.')
  system('gunzip -f chr22.fa.gz')
  fafile = 'chr22.fa'
  chr22_anno <- load_annotation(gtf, fafile)
  expect_equal(chr22_anno%>%.$uORF%>%sum, 5783)
  expect_equal(length(chr22_anno$cdsgrl), 7346)
  expect_equal(length(chr22_anno$longtrs), 451)
  expect_equal(length(chr22_anno$exonsgrl), 1923)
  expect_equal(length(names(chr22_anno$trspacecds)), 7346)
  expect_true(all(names(chr22_anno$exonsgrl)%in%fmcols(chr22_anno$cdsgrl, transcript_id)))
  expect_true(all(names(chr22_anno$cdsgrl)%in%names(chr22_anno$trspacecds)))
  testorfs <- chr22_anno$trspacecds%>%.[c(1,66,666,6666)]
  testtrs <- fmcols(chr22_anno$cdsgrl[names(testorfs)], transcript_id)
  exontestseq = extractTranscriptSeqs(chr22_anno$exonsgrl[testtrs],x=chr22_anno$fafileob)
  seqlevels(testorfs)<-names(exontestseq)
  testorfs <- resize(testorfs, width(testorfs)+3)
  orftestseq = exontestseq[testorfs]
  starts = Biostrings::subseq(orftestseq,1,3)%>%as.character
  names(starts)<-NULL
  expect_equal(starts, c('ATG','ATG','ATG','CTG'))
  stops = Biostrings::subseq(orftestseq,-3,-1)%>%as.character%>%
    setNames(NULL)
  expect_equal(stops, c('TGA','TAG','TGA','TGA'))
  expect_equal(chr22_anno$cdsgrl%>%.@unlistData%>%mcols%>%colnames,
    c("gene_id", "transcript_id", "gene_name", "type", "names"))


  extfasta = make_ext_fasta(gtf, fafile, outfasta='tmp.fa', fpext=50, tpext=50)
  extseqs = rtracklayer::import(extfasta, format='fasta')
  atgstarts <- extseqs%>%subseq(51,53)%>%setNames(NULL)%>%
    as.character%>%is_in('ATG')
  expect_true(mean(atgstarts)>0.99)
  extorfs <- names(extseqs)%>%str_extract('[^|]+')
  expect_true(all(extorfs%in%(chr22_anno$uORF%>%.[!.]%>%names)))

})
