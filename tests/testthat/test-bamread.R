test_that("reading bam file works", {
  data(chr22_anno)
  testbam <- system.file("extdata", "chr22.bam", package = "Ribostan", mustWork = TRUE)
  rpfs <- get_readgr(testbam, chr22_anno)
  expect_equal(length(rpfs), 65788)
  expect_true(all(seqnames(rpfs) %in% names(anno$exonsgrl)))
  expect_true(all(colnames(mcols(rpfs)) == "readlen"))
  expect_true(all(strand(rpfs) == "+"))
})
